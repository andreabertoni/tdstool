MODULE mod_splitstep
  USE mod_staticdata
  USE mod_sparse_fun
  USE mod_indata
  USE mod_filewrite
!  USE mod_fftw3params
  Use MKL_DFTI
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE SPLITSTEP_ALGO

  INTEGER nx, ny, irow, pt, iter, pot_file_changed
  COMPLEX*16 :: alpha, beta, k0, k1
  COMPLEX*16, ALLOCATABLE :: psik(:)
  REAL*8 dkx, dky, fftnorm
  REAL*8, ALLOCATABLE :: kx(:), ky(:)
  INTEGER INFO, file_list_index
  REAL*8 next_write_time
  REAL*8, ALLOCATABLE :: potx(:), poty(:)
  real elapsed, elapsed2, timess(4)
  
!  INTEGER*8 :: planfw, planbk
  type(DFTI_DESCRIPTOR), POINTER :: planfw
  integer   lengths(2)
  
  call INDATA_COMPUTE(INFO, 1)
  if (INFO == 1) then
    return
  end if

    ! Write Grid
  if (write_grid == "txt" .or. write_grid == "both") then
    call write_2D_grid(write_folder, "grid", xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
  end if
  if (write_grid == "bin" .or. write_grid == "both") then
    call write_2D_grid_bin(write_folder, "grid", xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
  end if

  alpha = -HBAR*HBAR/(2*mstar)
  beta = 4*PIG * nonlin_as * HBAR*HBAR / mstar
  k0 = -i * 0.5 * dt / HBAR
  k1 =  i * dt * alpha / HBAR
  dkx = 2*PIG/size_x
  dky = 2*PIG/size_y
  fftnorm = 1.0 / REAL(numx*numy)

  ALLOCATE (potx(1:numx))
  ALLOCATE (poty(1:numy))
  ALLOCATE (kx(numx))
  ALLOCATE (ky(numy))
  ALLOCATE (psik(numx*numy))
  kx(1:numx/2) = (/ ((nx-1)*dkx , nx= 1, numx/2) /)
  kx(numx/2+1:numx) = (/ (-(numx-nx+1)*dkx , nx= numx/2+1, numx) /)
  kx = kx**2;
  ky(1:numy/2) = (/ ((ny-1)*dky , ny= 1, numy/2) /)
  ky(numy/2+1:numy) = (/ (-(numy-ny+1)*dky , ny= numy/2+1, numy) /)
  ky = ky**2;

  call GPUSPLIT_INIT(numx, numy, ELCH, REAL(k0), IMAG(k0), REAL(k1), IMAG(k1), beta, kx, ky, pot_file_static, psi);

  call CPU_TIME(elapsed)

! **** Time steps
  next_write_time = 0.
  file_list_index = 0
  DO iter = 1, MAXIT

    ! Build time dependent potential
    call manage_pot_filelist((iter-1)*dt, file_list_index, pot_file_changed)
    call strtopot1D_time((iter-1) * dt, 1, mstar, strpotentialX, numx, xnodes, potx)
    call strtopot1D_time((iter-1) * dt, 1, mstar, strpotentialY, numy, ynodes, poty)
    call strtopot2D(1, mstar, strpotentialXY, numx, numy, xnodes, ynodes, pot)

    call GPUSPLIT_DO_STEP(pot, pot_filelist, potx, poty, pot_file_changed)
    print *, iter
    if (((iter-1)*dt) >= next_write_time) then
      call GPUSPLIT_GET_PSI(psi)
      call GPUSPLIT_GET_POT(pot)
      if (write_psi == "txt" .or. write_psi == "both") then
        call write_2D_cplx_in_file_gnuplot(write_folder, "psi", iter-1, psi, xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
      end if
      if (write_psi == "bin" .or. write_psi == "both") then
        call write_2D_cplx_in_file_bin(write_folder, "psi", iter-1, psi, numx, numy, write_downsample_x, write_downsample_y)
      end if

      if (write_pot == "txt" .or. write_pot == "both") then
        call write_2D_real_in_file_gnuplot(write_folder, "potential", iter-1, pot, xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
      end if
      if (write_pot == "bin" .or. write_pot == "both") then
        call write_2D_real_in_file_bin(write_folder, "potential", iter-1, pot, numx, numy, write_downsample_x, write_downsample_y)
      end if

      next_write_time = next_write_time + write_timestep
    end if
  END DO

  call CPU_TIME(elapsed2)
  print *, 'Tempo', elapsed2 - elapsed

  DEALLOCATE (poty)
  DEALLOCATE (potx)
  DEALLOCATE (psik)
  DEALLOCATE (ky)
  DEALLOCATE (kx)
  DEALLOCATE (psi)
  DEALLOCATE (pot)
  DEALLOCATE (pot_file_static)
  DEALLOCATE (pot_filelist)
  DEALLOCATE (xnodes)
  DEALLOCATE (ynodes)

END SUBROUTINE

END MODULE

