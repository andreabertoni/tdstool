#include <cuda.h>
#include <cublas.h>
#include <cufft.h>

extern "C" void TESTCUDA(void)
{
void *pt;
  cudaMalloc(&pt, 100000);
  cudaFree(pt);
}

typedef float2 Complex;

int c_numx, c_numy;
float c_elch, c_beta;
Complex c_k0, c_k1;
float *c_tmpbuf;
float *g_pot, *g_pot_file_static, *g_pot_filelist, *g_potx, *g_poty, *g_Kx, *g_Ky;
Complex *g_psi, *g_psik;
cufftHandle g_fft_hand;

void double_to_float(float *fpt, double *dpt, int n)
{
int i;

  for (i = 0; i < n; i++)
    fpt[i] = (float)dpt[i];
}

void float_to_double(double *dpt, float *fpt, int n)
{
int i;

  for (i = 0; i < n; i++)
    dpt[i] = (double)fpt[i];
}

int iDivUp(int a, int b)
{
  return (a % b != 0) ? (a / b + 1) : (a / b);
}

__global__ void
add_potential(int numx, int numy, float *pot, float *pot_file_static, float *pot_filelist, float *potx, float *poty, float ELCH)
{
  const int ix = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
  const int iy = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
  int idx = mul24(ix, numy) + iy;

  pot[idx] = (pot[idx] + pot_file_static[idx] + pot_filelist[idx] + potx[ix] + poty[iy])*ELCH;
}

__global__ void
nonlinear_half_step(int numx, int numy, Complex *psi, float *pot, Complex k0, float beta)
{
  const int ix = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
  const int iy = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
  int idx = mul24(ix, numy) + iy;
  Complex cpl, tmp;
  float v1, v2;
  v1 = pot[idx] + beta*(psi[idx].x*psi[idx].x + psi[idx].y*psi[idx].y);
  v2 = exp(v1*k0.x);
  cpl.x = v2 * cos(v1*k0.y);
  cpl.y = v2 * sin(v1*k0.y);
  tmp.x = psi[idx].x*cpl.x - psi[idx].y*cpl.y;
  tmp.y = psi[idx].y*cpl.x + psi[idx].x*cpl.y;
  psi[idx] = tmp;
}

__global__ void
nonlinear_half_step_scaled(int numx, int numy, Complex *psi, float *pot, Complex k0, float beta, float scaling)
{
  const int ix = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
  const int iy = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
  int idx = mul24(ix, numy) + iy;
  Complex cpl, tmp;
  float v1, v2;
  v1 = pot[idx] + beta*(psi[idx].x*psi[idx].x + psi[idx].y*psi[idx].y);
  v2 = scaling*exp(v1*k0.x);
  cpl.x = v2 * cos(v1*k0.y);
  cpl.y = v2 * sin(v1*k0.y);
  tmp.x = psi[idx].x*cpl.x - psi[idx].y*cpl.y;
  tmp.y = psi[idx].y*cpl.x + psi[idx].x*cpl.y;
  psi[idx] = tmp;
}

__global__ void
linear_step(int numx, int numy, Complex *psi, float *kx, float *ky, Complex k1)
{
  const int ix = __mul24(blockDim.x, blockIdx.x) + threadIdx.x;
  const int iy = __mul24(blockDim.y, blockIdx.y) + threadIdx.y;
  int idx = mul24(ix, numy) + iy;
  Complex cpl, tmp;
  float v1, v2;
  v1 = kx[ix] + ky[iy];
  v2 = exp(v1*k1.x);
  cpl.x = v2 * cos(v1*k1.y);
  cpl.y = v2 * sin(v1*k1.y);
  tmp.x = psi[idx].x*cpl.x - psi[idx].y*cpl.y;
  tmp.y = psi[idx].y*cpl.x + psi[idx].x*cpl.y;
  psi[idx] = tmp;
}

extern "C" void GPUSPLIT_INIT(int *numx, int *numy, double *elch, double *k0_real, double *k0_img,
		double *k1_real, double *k1_img, double *beta, double *Kx, double *Ky,
		double *pot_file_static, double *psi)
{
  c_numx = *numx;
  c_numy = *numy;
  c_elch = (float)(*elch);
  c_k0.x = (float)(*k0_real);
  c_k0.y = (float)(*k0_img);
  c_k1.x = (float)(*k1_real);
  c_k1.y = (float)(*k1_img);
  c_beta = (float)(*beta);
  c_tmpbuf = (float *)malloc(4*c_numx*c_numy*sizeof(float));
  cudaMalloc(&g_pot, 4*c_numx*c_numy);
  cudaMalloc(&g_pot_file_static, 4*c_numx*c_numy);
  cudaMalloc(&g_pot_filelist, 4*c_numx*c_numy);
  cudaMalloc(&g_potx, 4*c_numx);
  cudaMalloc(&g_poty, 4*c_numy);
  cudaMalloc(&g_psi, 2*4*c_numx*c_numy);
  cudaMalloc(&g_Kx, 4*c_numx);
  cudaMalloc(&g_Ky, 4*c_numy);
  
  double_to_float(c_tmpbuf, Kx, c_numx);
  cudaMemcpy(g_Kx, c_tmpbuf, c_numx*sizeof(float), cudaMemcpyHostToDevice);
  double_to_float(c_tmpbuf, Ky, c_numy);
  cudaMemcpy(g_Ky, c_tmpbuf, c_numy*sizeof(float), cudaMemcpyHostToDevice);
  double_to_float(c_tmpbuf, pot_file_static, c_numx*c_numy);
  cudaMemcpy(g_pot_file_static, c_tmpbuf, c_numx*c_numy*sizeof(float), cudaMemcpyHostToDevice);
  double_to_float(c_tmpbuf, psi, 2*c_numx*c_numy);
  cudaMemcpy(g_psi, c_tmpbuf, 2*c_numx*c_numy*sizeof(float), cudaMemcpyHostToDevice);
  
  cudaMemset(g_pot_filelist, 0, 4*c_numx*c_numy);
  
  cufftPlan2d(&g_fft_hand, c_numx, c_numy, CUFFT_C2C);
}

extern "C" void GPUSPLIT_DO_STEP(double *pot, double *pot_filelist, double *potx, double *poty, int *filelist_changed)
{
  double_to_float(c_tmpbuf, pot, c_numx*c_numy);
  cudaMemcpy(g_pot, c_tmpbuf, c_numx*c_numy*sizeof(float), cudaMemcpyHostToDevice);
  if (*filelist_changed)
  {
    double_to_float(c_tmpbuf, pot_filelist, c_numx*c_numy);
    cudaMemcpy(g_pot_filelist, c_tmpbuf, c_numx*c_numy*sizeof(float), cudaMemcpyHostToDevice);
  }
  double_to_float(c_tmpbuf, potx, c_numx);
  cudaMemcpy(g_potx, c_tmpbuf, c_numx*sizeof(float), cudaMemcpyHostToDevice);
  double_to_float(c_tmpbuf, poty, c_numy);
  cudaMemcpy(g_poty, c_tmpbuf, c_numy*sizeof(float), cudaMemcpyHostToDevice);

  dim3 threadBlock(16, 16);
  dim3 kernelBlockGrid(iDivUp(c_numx, threadBlock.x), iDivUp(c_numy, threadBlock.y));
  add_potential<<<kernelBlockGrid, threadBlock>>>(c_numx, c_numy, g_pot, g_pot_file_static, g_pot_filelist, g_potx, g_poty, c_elch);

  nonlinear_half_step<<<kernelBlockGrid, threadBlock>>>(c_numx, c_numy, g_psi, g_pot, c_k0, c_beta);
  cufftExecC2C(g_fft_hand, (cufftComplex *)g_psi, (cufftComplex *)g_psi, CUFFT_FORWARD);
  linear_step<<<kernelBlockGrid, threadBlock>>>(c_numx, c_numy, g_psi, g_Kx, g_Ky, c_k1);
  cufftExecC2C(g_fft_hand, (cufftComplex *)g_psi, (cufftComplex *)g_psi, CUFFT_INVERSE);

  nonlinear_half_step_scaled<<<kernelBlockGrid, threadBlock>>>(c_numx, c_numy, g_psi, g_pot, c_k0, c_beta, 1.0F/((float)(c_numx*c_numy)));
}

extern "C" void GPUSPLIT_GET_PSI(double *psi)
{
  cudaMemcpy(c_tmpbuf, g_psi, 2*c_numx*c_numy*sizeof(float), cudaMemcpyDeviceToHost);
  float_to_double(psi, c_tmpbuf, 2*c_numx*c_numy);
}

extern "C" void GPUSPLIT_GET_POT(double *pot)
{
  cudaMemcpy(c_tmpbuf, g_pot, c_numx*c_numy*sizeof(float), cudaMemcpyDeviceToHost);
  float_to_double(pot, c_tmpbuf, c_numx*c_numy);
}

extern "C" void GPUSPLIT_DESTROY(void)
{
  free(c_tmpbuf);
  cudaFree(&g_pot);
  cudaFree(&g_pot_file_static);
  cudaFree(&g_pot_filelist);
  cudaFree(&g_potx);
  cudaFree(&g_poty);
  cudaFree(&g_psi);
  cudaFree(&g_Kx);
  cudaFree(&g_Ky);
  cufftDestroy(g_fft_hand);
}
