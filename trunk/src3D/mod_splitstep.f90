MODULE mod_mainalgo
  USE mod_staticdata
  USE mod_indata
  USE mod_filewrite
  USE mod_fftw3params
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE MAIN_ALGO

  INTEGER nx, ny, nz, irow, pt, iter
  COMPLEX*16 :: alpha, beta, k0, k1
  COMPLEX*16, ALLOCATABLE :: psik(:)
  REAL*8 dkx, dky, dkz, fftnorm
  REAL*8, ALLOCATABLE :: kx(:), ky(:), kz(:)
  INTEGER INFO
  REAL*8 next_write_time
  INTEGER*8 :: planfw, planbk

  call INDATA_COMPUTE

    ! Write Potential
  if (write_pot == "txt" .or. write_pot == "both") then
    call write_3D_real_in_file_gnuplot(write_folder, "potential", 0, pot, xnodes, ynodes, znodes, numx, numy, numz, write_downsample_x, write_downsample_y, write_downsample_z)
  end if
  if (write_pot == "bin" .or. write_pot == "both") then
    call write_3D_real_in_file_bin(write_folder, "potential", 0, pot, numx, numy, numz, write_downsample_x, write_downsample_y, write_downsample_z)
  end if

    ! Write Grid
  if (write_grid) then
    call write_3D_grid(write_folder, "grid", xnodes, ynodes, znodes, numx, numy, numz, write_downsample_x, write_downsample_y, write_downsample_z)
  end if

  if (write_grid == 1) then
  end if

  alpha = -HBAR*HBAR/(2*mstar)
  beta = 4*PIG * nonlin_as * HBAR*HBAR / mstar
  k0 = -i * 0.5 * dt / HBAR
  k1 =  i * dt * alpha / HBAR
  dkx = 2*PIG/size_x
  dky = 2*PIG/size_y
  dkz = 2*PIG/size_z
  fftnorm = 1.0 / REAL(numx*numy*numz)

  ALLOCATE (kx(numx))
  ALLOCATE (ky(numy))
  ALLOCATE (kz(numz))
  ALLOCATE (psik(numx*numy*numz))
  kx(1:numx/2) = (/ ((nx-1)*dkx , nx= 1, numx/2) /)
  kx(numx/2+1:numx) = (/ (-(numx-nx+1)*dkx , nx= numx/2+1, numx) /)
  kx = kx**2;
  ky(1:numy/2) = (/ ((ny-1)*dky , ny= 1, numy/2) /)
  ky(numy/2+1:numy) = (/ (-(numy-ny+1)*dky , ny= numy/2+1, numy) /)
  ky = ky**2;
  kz(1:numz/2) = (/ ((nz-1)*dkz , nz= 1, numz/2) /)
  kz(numz/2+1:numz) = (/ (-(numz-nz+1)*dkz , nz= numz/2+1, numz) /)
  kz = kz**2;

  CALL dfftw_plan_dft_3d( planfw, numx, numy, numz, psi, psik,           &
     &                                   FFTW_FORWARD,  FFTW_ESTIMATE )
  CALL dfftw_plan_dft_3d( planbk, numx, numy, numz, psik, psi,           &
     &                                   FFTW_BACKWARD, FFTW_ESTIMATE )

! **** Time steps
  next_write_time = 0.
  DO iter = 1, MAXIT

    ! ****** Nonlinear half step
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        DO nz = 1, numz
          psi(pt) = psi(pt) * exp(k0 * (pot(nz, ny, nx) + beta*(abs(psi(pt))**2)));
          pt = pt + 1
        END DO
      END DO
    END DO

    ! ****** Linear step
    CALL dfftw_execute( planfw )
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        DO nz = 1, numz
          psik(pt) = psik(pt) * exp(k1 * (kx(nx) + ky(ny) + kz(nz)));
          pt = pt + 1
        END DO
      END DO
    END DO
    CALL dfftw_execute( planbk )
    psi = psi * fftnorm;

    ! ****** Nonlinear half step
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        DO nz = 1, numz
          psi(pt) = psi(pt) * exp(k0 * (pot(nz, ny, nx) + beta*(abs(psi(pt))**2)));
          pt = pt + 1
        END DO
      END DO
    END DO

    write (*,*) iter
    if (((iter-1)*dt) >= next_write_time) then
      if (write_psi == "txt" .or. write_psi == "both") then
        call write_3D_cplx_in_file_gnuplot(write_folder, "psi", iter-1, psi, xnodes, ynodes, znodes, numx, numy, numz, write_downsample_x, write_downsample_y, write_downsample_z)
      end if
      if (write_psi == "bin" .or. write_psi == "both") then
        call write_3D_cplx_in_file_bin(write_folder, "psi", iter-1, psi, numx, numy, numz, write_downsample_x, write_downsample_y, write_downsample_z)
      end if
      next_write_time = next_write_time + write_timestep
    end if
  END DO

  CALL dfftw_destroy_plan( planfw )
  CALL dfftw_destroy_plan( planbk )
  DEALLOCATE (psik)
  DEALLOCATE (kz)
  DEALLOCATE (ky)
  DEALLOCATE (kx)
  DEALLOCATE (psi)
  DEALLOCATE (pot)
  DEALLOCATE (znodes)
  DEALLOCATE (ynodes)
  DEALLOCATE (xnodes)

END SUBROUTINE

END MODULE

