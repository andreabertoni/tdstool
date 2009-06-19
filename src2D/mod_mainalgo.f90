MODULE mod_mainalgo
  USE mod_staticdata
  USE mod_sparse_fun
  USE mod_indata
  USE mod_filewrite
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE MAIN_ALGO

  INTEGER nx, ny, irow, pt, iter
  COMPLEX*16 k0, k1;
  COMPLEX*16, ALLOCATABLE :: b(:)
  INTEGER INFO
  REAL*8 next_write_time

! PARDISO data
  INTEGER*8  pt_prd(64)
  INTEGER :: iparm_prd(64), msglvl_prd
  INTEGER, ALLOCATABLE :: perm_prd(:)

  call INDATA_COMPUTE(INFO)
  if (INFO == 1) then
    return
  end if

print *, "Kx = ", sqrt(2*mstar*xenergy*ELCH)/HBAR
print *, "Ky = ", sqrt(2*mstar*yenergy*ELCH)/HBAR
print *, "x_period = ", 2*PIG/(sqrt(2*mstar*xenergy*ELCH)/HBAR)
print *, "y_period = ", 2*PIG/(sqrt(2*mstar*yenergy*ELCH)/HBAR)
print *, "dx = ", xnodes(10)-xnodes(9)
print *, "dy = ", ynodes(10)-ynodes(9)

  ALLOCATE (b(numx*numy))

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

  k0 = -i*HBAR/(2*mstar)
  k1 = -i/HBAR

  call make_box_stiffness_2D(S, A, ia, ja, numx, numy, xnodes, ynodes)
  pt = ia(1)
  irow = 1
  DO nx = 1,numx
    DO ny = 1,numy
      DO WHILE (pt < ia(irow+1))
        if (ja(pt) == irow) then
          A(pt) = S(irow) - (A(pt)*k0 + pot(ny, nx)*k1*S(irow))*0.5*dt
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

  ! 13 = complex and unsymmetric, 12 = phases 1 & 2
  call pardiso(pt_prd, 1, 1, 6, 12, numx*numy, A, ia, ja, perm_prd, 1, iparm_prd, 1, b, b, INFO)

  next_write_time = 0.
! **** Time steps
  DO iter = 1, MAXIT
!    psi = (S-A) \ ((S+A)*psi);
    call z_sparse_matvet(numx*numy, A, ia, ja, psi, b)
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
    call pardiso(pt_prd, 1, 1, 6, 33, numx*numy, A, ia, ja, perm_prd, 1, iparm_prd, 0, b, psi, INFO)
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
      next_write_time = next_write_time + write_timestep
    end if
  END DO

!open(10, FILE='schro_out2.dat', STATUS='REPLACE')
!DO ii = 1, nx, 100
!  WRITE(10,*) ii*dx, abs(y(ii))
!END DO
!close(10)

  call pardiso(pt_prd, 1, 1, 6, -1, numx*numy, A, ia, ja, perm_prd, 1, iparm_prd, 0, b, b, INFO)

  DEALLOCATE(perm_prd)
  DEALLOCATE (b)
  DEALLOCATE (ia)
  DEALLOCATE (ja)
  DEALLOCATE (A)
  DEALLOCATE (S)
  DEALLOCATE (psi)
  DEALLOCATE (pot)
  DEALLOCATE (xnodes)
  DEALLOCATE (ynodes)

END SUBROUTINE

END MODULE

