MODULE mod_mainalgo
  USE mod_staticdata
  USE mod_sparse_fun
  USE mod_indata
  USE mod_filewrite
  USE mod_fftw3params
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE MAIN_ALGO

  INTEGER nx, ny, irow, pt, iter
  COMPLEX*16 :: alpha, beta, k0, k1
  COMPLEX*16, ALLOCATABLE :: psik(:)
  REAL*8 dkx, dky, fftnorm
  REAL*8, ALLOCATABLE :: kx(:), ky(:)
  INTEGER INFO
  REAL*8 next_write_time
  INTEGER*8 :: planfw, planbk

  call INDATA_COMPUTE

    ! Write Potential
  if (write_pot == "txt" .or. write_pot == "both") then
    call write_2D_real_in_file_gnuplot(write_folder, "potential", 0, pot, xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
  end if
  if (write_pot == "bin" .or. write_pot == "both") then
    call write_2D_real_in_file_bin(write_folder, "potential", 0, pot, numx, numy, write_downsample_x, write_downsample_y)
  end if

    ! Write Grid
  if (write_grid) then
    call write_2D_grid(write_folder, "grid", xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
  end if

  if (write_grid == 1) then
  end if

  alpha = -HBAR*HBAR/(2*mstar)
  beta = 4*PIG * nonlin_as * HBAR*HBAR / mstar
  k0 = -i * 0.5 * dt / HBAR
  k1 =  i * dt * alpha / HBAR
  dkx = 2*PIG/size_x
  dky = 2*PIG/size_y
  fftnorm = 1.0 / REAL(numx*numy)

  ALLOCATE (kx(numx))
  ALLOCATE (ky(numy))
  ALLOCATE (psik(numx*numy))
  kx(1:numx/2) = (/ ((nx-1)*dkx , nx= 1, numx/2) /)
  kx(numx/2+1:numx) = (/ (-(numx-nx+1)*dkx , nx= numx/2+1, numx) /)
  kx = kx**2;
  ky(1:numy/2) = (/ ((ny-1)*dky , ny= 1, numy/2) /)
  ky(numy/2+1:numy) = (/ (-(numy-ny+1)*dky , ny= numy/2+1, numy) /)
  ky = ky**2;

  CALL dfftw_plan_dft_2d( planfw, numx, numy, psi, psik,           &
     &                                   FFTW_FORWARD,  FFTW_ESTIMATE )
  CALL dfftw_plan_dft_2d( planbk, numx, numy, psik, psi,           &
     &                                   FFTW_BACKWARD, FFTW_ESTIMATE )

! **** Time steps
  next_write_time = 0.
  DO iter = 1, MAXIT

    ! ****** Nonlinear half step
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        psi(pt) = psi(pt) * exp(k0 * (pot(ny, nx) + beta*(abs(psi(pt))**2)));
        pt = pt + 1
      END DO
    END DO

    ! ****** Linear step
    CALL dfftw_execute( planfw )
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        psik(pt) = psik(pt) * exp(k1 * (kx(nx) + ky(ny)));
        pt = pt + 1
      END DO
    END DO
    CALL dfftw_execute( planbk )
    psi = psi * fftnorm;

    ! ****** Nonlinear half step
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        psi(pt) = psi(pt) * exp(k0 * (pot(ny, nx) + beta*(abs(psi(pt))**2)));
        pt = pt + 1
      END DO
    END DO

    write (*,*) iter
    if (((iter-1)*dt) >= next_write_time) then
      if (write_psi == "txt" .or. write_psi == "both") then
        call write_2D_cplx_in_file_gnuplot(write_folder, "psi", iter-1, psi, xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
      end if
      if (write_psi == "bin" .or. write_psi == "both") then
        call write_2D_cplx_in_file_bin(write_folder, "psi", iter-1, psi, numx, numy, write_downsample_x, write_downsample_y)
      end if
      next_write_time = next_write_time + write_timestep
    end if
  END DO

  CALL dfftw_destroy_plan( planfw )
  CALL dfftw_destroy_plan( planbk )
  DEALLOCATE (psik)
  DEALLOCATE (ky)
  DEALLOCATE (kx)
  DEALLOCATE (psi)
  DEALLOCATE (pot)
  DEALLOCATE (xnodes)
  DEALLOCATE (ynodes)

END SUBROUTINE

END MODULE
