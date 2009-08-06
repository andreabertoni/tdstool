MODULE mod_boxintegration
  USE mod_staticdata
  USE mod_sparse_fun
  USE mod_indata
  USE mod_filewrite
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE BOXINTEGRATION_ALGO

  INTEGER nx, ny, irow, pt, iter
  COMPLEX*16 k0, k1;
  COMPLEX*16, ALLOCATABLE :: b(:)
  INTEGER INFO, file_list_index
  REAL*8 next_write_time, zero
  REAL*8, ALLOCATABLE :: potx(:), poty(:)
  COMPLEX*16, ALLOCATABLE :: AA(:)

! PARDISO data
  INTEGER*8  pt_prd(64)
  INTEGER :: iparm_prd(64), msglvl_prd
  INTEGER, ALLOCATABLE :: perm_prd(:)
  
  zero = 0.

  if (nonlin_as /= 0.) then
    call dwgbut("Nonlinear term will be ignored in box integration method");
  end if
  
  call INDATA_COMPUTE(INFO, 0)
  if (INFO == 1) then
    return
  end if

  ALLOCATE (AA(5*numx*numy))
  ALLOCATE (b(numx*numy))
  ALLOCATE (potx(1:numx))
  ALLOCATE (poty(1:numy))

    ! Write Grid
  if (write_grid == "txt" .or. write_grid == "both") then
    call write_2D_grid(write_folder, "grid", xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
  end if
  if (write_grid == "bin" .or. write_grid == "both") then
    call write_2D_grid_bin(write_folder, "grid", xnodes, ynodes, numx, numy, write_downsample_x, write_downsample_y)
  end if

  k0 = -i*HBAR/(2*mstar)
  k1 = -i/HBAR

  if (magnetic /= 0.) then
    call make_box_stiffness_2D(S, A, ia, ja, numx, numy, xnodes, ynodes, ELCH*magnetic/(numx*numy*2*PIG*HBAR))
  else
    call make_box_stiffness_2D_symmetric(S, A, ia, ja, numx, numy, xnodes, ynodes, zero)
  end if
  pt = ia(1)
  irow = 1
  DO nx = 1,numx
    DO ny = 1,numy
      DO WHILE (pt < ia(irow+1))
        if (ja(pt) == irow) then
!          A(pt) = S(irow) - (A(pt)*k0 + pot(ny, nx)*k1*S(irow))*0.5*dt
          A(pt) = S(irow) - (A(pt)*k0)*0.5*dt
        else
          A(pt) = -A(pt) * (0.5*k0*dt)
        end if
        pt = pt+1
      END DO
    irow = irow + 1
    END DO
  END DO

! PARDISO
  ALLOCATE (perm_prd(numx*numy))
  pt_prd(:) = 0
  perm_prd(:) = 0
  iparm_prd(:) = 0
  iparm_prd(1) = 0
  iparm_prd(3) = 4   ! number of threads
  msglvl_prd = 0

  next_write_time = 0.
  file_list_index = 0
! **** Time steps
  DO iter = 1, MAXIT

    ! Build time dependent potential
    call manage_pot_filelist((iter-1)*dt, file_list_index)
    call strtopot1D_time((iter-1) * dt, 1, mstar, strpotentialX, numx, xnodes, potx)
    call strtopot1D_time((iter-1) * dt, 1, mstar, strpotentialY, numy, ynodes, poty)
    call strtopot2D(1, mstar, strpotentialXY, numx, numy, xnodes, ynodes, pot)

    pt = ia(1)
    irow = 1
    DO nx = 1, numx
      DO ny = 1, numy
        pot(ny, nx) = (pot(ny, nx) + pot_file_static(ny, nx) + potx(nx) + poty(ny))*ELCH
        DO WHILE (pt < ia(irow+1))
          if (ja(pt) == irow) then
            AA(pt) = A(pt) - pot(ny, nx)*k1*S(irow)*0.5*dt
          else
            AA(pt) = A(pt)
          end if
          pt = pt+1
        END DO
      irow = irow + 1
      END DO
    END DO

!    psi = (S-A) \ ((S+A)*psi);
    if (magnetic /= 0.) then
      call z_sparse_matvet(numx*numy, AA, ia, ja, psi, b)
    else
      call z_sparse_matvet_symmetric(numx*numy, AA, ia, ja, psi, b)
    end if
    DO nx = 1, numx*numy
      b(nx) = 2.*S(nx)*psi(nx)-b(nx)
    END DO

!  open(30, FILE="linprob.bin", FORM="UNFORMATTED", STATUS='REPLACE')
!  write(30) numx*numy
!  write(30) A
!  write(30) ia
!  write(30) ja
!  write(30) b
!  close(30)
!
    print*, 'ITER: ', iter
    if (magnetic /= 0.) then
      call pardiso(pt_prd, 1, 1, 13, 13, numx*numy, AA, ia, ja, perm_prd, 1, iparm_prd, 0, b, psi, INFO)
    else
      call pardiso(pt_prd, 1, 1, 6, 13, numx*numy, AA, ia, ja, perm_prd, 1, iparm_prd, 0, b, psi, INFO)
    end if
    IF (INFO /= 0) THEN
      PRINT*, "pardiso fallita"
    END IF

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

!open(10, FILE='schro_out2.dat', STATUS='REPLACE')
!DO ii = 1, nx, 100
!  WRITE(10,*) ii*dx, abs(y(ii))
!END DO
!close(10)

  if (nonlin_as /= 0.) then
    call pardiso(pt_prd, 1, 1, -4, -1, numx*numy, A, ia, ja, perm_prd, 1, iparm_prd, 0, b, b, INFO)
  else
    call pardiso(pt_prd, 1, 1, 6, -1, numx*numy, A, ia, ja, perm_prd, 1, iparm_prd, 0, b, b, INFO)
  end if
  
  DEALLOCATE(perm_prd)
  DEALLOCATE (poty)
  DEALLOCATE (potx)
  DEALLOCATE (b)
  DEALLOCATE (ia)
  DEALLOCATE (ja)
  DEALLOCATE (AA)
  DEALLOCATE (A)
  DEALLOCATE (S)
  DEALLOCATE (psi)
  DEALLOCATE (pot)
  DEALLOCATE (pot_file_static)
  DEALLOCATE (pot_filelist)
  DEALLOCATE (xnodes)
  DEALLOCATE (ynodes)

END SUBROUTINE

END MODULE

