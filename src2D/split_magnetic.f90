MODULE mod_splitmagnetic
  USE mod_staticdata
  USE mod_sparse_fun
  USE mod_indata
  USE mod_filewrite
!  USE mod_fftw3params
  Use MKL_DFTI
  USE mod_solve1dLAPACK
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE SPLITMAGNETIC_ALGO

  INTEGER pt, iter
  COMPLEX*16 :: alpha, beta, k0, k1
  INTEGER INFO, file_list_index
  REAL*8 next_write_time
  REAL*8, ALLOCATABLE :: potx(:), poty(:)

  REAL*8, ALLOCATABLE :: xgrid(:), ygrid(:)
  REAL*8, ALLOCATABLE :: kxgrid(:), kygrid(:)
  COMPLEX*16, ALLOCATABLE :: psik(:,:)
  REAL*8, ALLOCATABLE :: phiho(:,:, :)
  REAL*8, ALLOCATABLE :: phiho5(:,:)
  REAL*8, ALLOCATABLE :: shiftedxgrid2(:)
  REAL*8, ALLOCATABLE :: potho(:), potho5(:)
  REAL*8, ALLOCATABLE :: energiesho(:,:)
  REAL*8 :: cnormx, cnormy, constho
  REAL*8 :: dx, dy, ky
  INTEGER :: nn, nx, ny, nkx, nky, numkx, numky
  REAL*8 :: scale, omegac, lb2, dkx, dky

  type(DFTI_DESCRIPTOR), POINTER :: planfw
  integer   lengths(1)

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

  omegac= ELCH*magnetic/mstar
  lb2= HBAR/(ELCH*magnetic)
  numkx= numx
  numky= numy
  dx = size_x / numx
  dy = size_y / numy
  dkx= 2*PIG/size_x
  dky= 2*PIG/size_y

  ALLOCATE( psik(0:numx-1, 0:numy-1) )
  ALLOCATE( phiho(0:numx-1, 0:numkx-1, 0:numky-1) )
  ALLOCATE( phiho5(0:1*numx-1, 0:numkx-1) )
  ALLOCATE( shiftedxgrid2(0:1*numx-1) )
  ALLOCATE( potho(0:numx-1) )
  ALLOCATE( potho5(0:1*numx-1) )
  ALLOCATE( energiesho(0:numkx-1, 0:numky-1) )
  ALLOCATE( xgrid(0:numx-1) )
  ALLOCATE( ygrid(0:numy-1) )
  ALLOCATE( kxgrid(0:numkx-1) )
  ALLOCATE( kygrid(0:numky-1) )

  xgrid= (/ (nn*dx+dx/2., nn= 0,numx-1 ) /)
  ygrid= (/ (nn*dy+dx/2., nn= 0,numy-1 ) /)
  kxgrid(0:numkx/2)= (/ (nn*dkx, nn= 0,numkx/2 ) /)
  kygrid(0:numky/2)= (/ (nn*dky, nn= 0,numky/2 ) /)
  kxgrid(numkx/2+1:numkx-1)= (/ ( -(numkx-nn)*dkx, nn= numkx/2+1,numkx-1 ) /)
  kygrid(numky/2+1:numky-1)= (/ ( -(numky-nn)*dky, nn= numky/2+1,numky-1 ) /)

  cnormy= ( 1. / SQRT(REAL(numy,8)) )
  constho= (HBAR/(2.*mstar)) * HBAR/(dx*dx)

  lengths(1) = numy
  INFO = DftiCreateDescriptor( planfw, DFTI_DOUBLE, DFTI_COMPLEX, 1, lengths)
  INFO = DftiCommitDescriptor( planfw )

  DO nky= 0, numky-1
    ky= kygrid(nky)
    DO nx= 0, 1*numx-1
      shiftedxgrid2(nx)= (nx-0*numx)*dx+dx/2 - lb2*ky
      potho5(nx)= mstar/2. * (omegac * shiftedxgrid2(nx))**2
    END DO
    CALL solve1dLAPACK( 1*numx-1, numkx, numkx, constho, potho5, energiesho(:, nky), phiho5 )
    DO nkx = 0, numkx-1
      DO nx= 0, numx-1
        phiho(nx, nkx, nky) = phiho5(nx+0*numx, nkx);
      END DO
    END DO
  END DO


  alpha = -HBAR*HBAR/(2*mstar)
  beta = 4*PIG * nonlin_as * HBAR*HBAR / mstar
  k0 = -i * 0.5 * dt / HBAR
  k1 =  i * dt * alpha / HBAR

  ALLOCATE (potx(1:numx))
  ALLOCATE (poty(1:numy))

! **** Time steps
  next_write_time = 0.
  file_list_index = 0
  DO iter = 1, MAXIT

    ! Build time dependent potential
    call manage_pot_filelist((iter-1)*dt, file_list_index)
    call strtopot1D_time((iter-1) * dt, 1, mstar, strpotentialX, numx, xnodes, potx)
    call strtopot1D_time((iter-1) * dt, 1, mstar, strpotentialY, numy, ynodes, poty)
    call strtopot2D(1, mstar, strpotentialXY, numx, numy, xnodes, ynodes, pot)
    DO nx = 1, numx
      DO ny = 1, numy
        pot(ny, nx) = (pot(ny, nx) + pot_file_static(ny, nx) + pot_filelist(ny, nx) + potx(nx) + poty(ny))*ELCH
      END DO
    END DO

    ! ****** Nonlinear half step
    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        psi(pt) = psi(pt) * exp(k0 * (pot(ny, nx) + beta*(abs(psi(pt))**2)));
        pt = pt + 1
      END DO
    END DO

    ! ****** Linear step
  DO nx = 0, numx-1
    INFO = DftiComputeForward(planfw, psi(nx*numy+1:(nx+1)*numy))
  END DO

  DO nky= 0, numky-1
    DO nkx= 0, numkx-1
      psik(nkx, nky) = 0.
      DO nx= 0, numx-1
        psik(nkx, nky) = psik(nkx, nky) + phiho(nx, nkx, nky)*psi(nx*numy+nky+1) ! psi(nky, nx)
      END DO
    END DO
  END DO

  DO nky= 0, numky-1
    DO nkx= 0, numkx-1
      psik(nkx,nky) = psik(nkx,nky) * EXP(-i*dt * energiesho(nkx, nky) / HBAR)
    END DO
  END DO

  DO nky = 0, numky-1
    DO nx = 0, numx-1
      psi(nx*numy+nky+1) = 0.
      DO nkx = 0, numkx-1
        psi(nx*numy+nky+1) = psi(nx*numy+nky+1) + phiho(nx, nkx, nky)*psik(nkx, nky)
      END DO
    END DO
  END DO

  scale = 1.0 / REAL(lengths(1), 8)
  DO nx = 0, numx-1
    INFO = DftiComputeBackward(planfw, psi(nx*numy+1:(nx+1)*numy))
  END DO
  psi = psi * scale

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

      if (write_pot == "txt" .or. write_pot == "both") then
        call write_2D_real_in_file_gnuplot(write_folder, "potential", iter-1, pot, xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
      end if
      if (write_pot == "bin" .or. write_pot == "both") then
        call write_2D_real_in_file_bin(write_folder, "potential", iter-1, pot, numx, numy, write_downsample_x, write_downsample_y)
      end if

      next_write_time = next_write_time + write_timestep
    end if
  END DO

!  CALL dfftw_destroy_plan( planfw )
!  CALL dfftw_destroy_plan( planbk )
  INFO = DftiFreeDescriptor(planfw)
  
  DEALLOCATE (psik)
  DEALLOCATE (phiho)
  DEALLOCATE (phiho5)
  DEALLOCATE (shiftedxgrid2)
  DEALLOCATE (xgrid)
  DEALLOCATE (ygrid)
  DEALLOCATE (kxgrid)
  DEALLOCATE (kygrid)
  DEALLOCATE (potho)
  DEALLOCATE (potho5)
  DEALLOCATE (energiesho)
  DEALLOCATE (poty)
  DEALLOCATE (potx)
  DEALLOCATE (psi)
  DEALLOCATE (pot)
  DEALLOCATE (pot_file_static)
  DEALLOCATE (pot_filelist)
  DEALLOCATE (xnodes)
  DEALLOCATE (ynodes)

END SUBROUTINE

END MODULE
