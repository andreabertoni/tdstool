MODULE mod_sparse_fun
  USE mod_staticdata
  IMPLICIT NONE
  CONTAINS

! y = A*x
SUBROUTINE z_sparse_matvet(n, A, ia, ja, x, y)
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, INTENT(IN) :: A(:), x(:)
  COMPLEX*16, INTENT(OUT) :: y(:)
  INTEGER, INTENT(IN) :: ia(:), ja(:)

  INTEGER irow, pt, ii

  y(:) = 0
  pt = ia(1)
  DO irow = 1,n
    DO WHILE (pt < ia(irow+1))
      y(irow) = y(irow) + A(pt)*x(ja(pt))
      pt = pt+1
    END DO
  END DO

END SUBROUTINE

SUBROUTINE z_sparse_matvet_symmetric(n, A, ia, ja, x, y)
  INTEGER, INTENT(IN) :: n
  COMPLEX*16, INTENT(IN) :: A(:), x(:)
  COMPLEX*16, INTENT(OUT) :: y(:)
  INTEGER, INTENT(IN) :: ia(:), ja(:)

  INTEGER irow, pt, ii

  y(:) = 0
  pt = ia(1)
  DO irow = 1,n
    DO WHILE (pt < ia(irow+1))
      y(irow) = y(irow) + A(pt)*x(ja(pt))
      if (irow /= ja(pt)) then
        y(ja(pt)) = y(ja(pt)) + A(pt)*x(irow)
      end if
      pt = pt+1
    END DO
  END DO

END SUBROUTINE

! Creates the stiffness matrix, in column-sparse format, using the box integration method
! with a non uniform orthogonal grid of size (numx+2) x (numy+2)
! whose coordinates are in xnodes e ynodes.
! It uses null Dirichlet border conditions.

SUBROUTINE make_box_stiffness_2D(S, A, ia, ja, numx, numy, xnodes, ynodes, alpha_magnetic)
  COMPLEX*16, ALLOCATABLE, INTENT(OUT) :: A(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: S(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: ia(:), ja(:)
  INTEGER, INTENT(IN) :: numx, numy
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:)
  REAL*8, INTENT(IN) :: alpha_magnetic

  REAL*8, ALLOCATABLE :: anodes(:), bnodes(:)
  INTEGER nx, ny
  INTEGER row  ! row index into the stiffness matrix
  INTEGER pt   ! index into the A values array
  COMPLEX*16 k1, k2, k3, k4, expc, expcm

  ALLOCATE (S(numx*numy))
  ALLOCATE (A(5*numx*numy))
  ALLOCATE (ja(5*numx*numy))
  ALLOCATE (ia(numx*numy+1))
  ALLOCATE (anodes(numx+1))
  ALLOCATE (bnodes(numy+1))

  anodes = xnodes(2:numx+2) - xnodes(1:numx+1)
  bnodes = ynodes(2:numy+2) - ynodes(1:numy+1)

  pt = 1;
  row = 1;
  DO nx = 1,numx
    DO ny = 1,numy
      !row = (nx-1)*numy+ny
      S(row) = (anodes(nx)+anodes(nx+1)) * (bnodes(ny)+bnodes(ny+1)) / 4
      k1 = (bnodes(ny)+bnodes(ny+1)) / (2*anodes(nx))
      k2 = (bnodes(ny)+bnodes(ny+1)) / (2*anodes(nx+1))
      k3 = (anodes(nx)+anodes(nx+1)) / (2*bnodes(ny))
      k4 = (anodes(nx)+anodes(nx+1)) / (2*bnodes(ny+1))
      expc = exp(i*(ny-1)*alpha_magnetic);
      expcm = exp(-i*(ny-1)*alpha_magnetic);
!      write (99, *) ny, nx, expc

      ia(row) = pt
      if (nx == 1) then
        if (ny == 1) then
          ! angolo alto a sinistra
          A(pt) = k1+k2+k3+k4
          ja(pt) = row
          A(pt+1) = -k4
          ja(pt+1) = row+1;
          A(pt+2) = -k2 * expc
          ja(pt+2) = row+numy;
          pt = pt + 3;
	    else if (ny == numy) then
          ! angolo basso a sinistra
          A(pt) = -k3
          ja(pt) = row-1;
          A(pt+1) = k1+k2+k3+k4
          ja(pt+1) = row;
          A(pt+2) = -k2 * expc
          ja(pt+2) = row+numy;
          pt = pt + 3;
	else
          ! bordo sinistro
          A(pt) = -k3
          ja(pt) = row-1
          A(pt+1) = k1+k2+k3+k4
          ja(pt+1) = row;
          A(pt+2) = -k4
          ja(pt+2) = row+1
          A(pt+3) = -k2 * expc
          ja(pt+3) = row+numy;
          pt = pt + 4;
	end if
      else if (nx == numx) then
        if (ny == 1) then
	  ! angolo alto a destra
          A(pt) = -k1 * expcm
          ja(pt) = row-numy;
          A(pt+1) = k1+k2+k3+k4
          ja(pt+1) = row;
          A(pt+2) = -k4
          ja(pt+2) = row+1;
          pt = pt + 3;
        else if (ny == numy) then
	  ! angolo basso a destra
          A(pt) = -k1 * expcm
          ja(pt) = row-numy;
          A(pt+1) = -k3
          ja(pt+1) = row-1;
          A(pt+2) = k1+k2+k3+k4
          ja(pt+2) = row;
          pt = pt + 3;
        else
	  ! bordo destro
          A(pt) = -k1 * expcm
          ja(pt) = row-numy;
          A(pt+1) = -k3
          ja(pt+1) = row-1;
          A(pt+2) = k1+k2+k3+k4
          ja(pt+2) = row;
          A(pt+3) = -k4
          ja(pt+3) = row+1;
          pt = pt + 4;
        end if
      else
        if (ny == 1) then
          ! bordo alto
          A(pt) = -k1 * expcm
          ja(pt) = row-numy;
          A(pt+1) = k1+k2+k3+k4
          ja(pt+1) = row;
          A(pt+2) = -k4
          ja(pt+2) = row+1;
          A(pt+3) = -k2 * expc
          ja(pt+3) = row+numy;
          pt = pt + 4;
        else if (ny == numy) then
          ! bordo basso
          A(pt) = -k1 * expcm
          ja(pt) = row-numy;
          A(pt+1) = -k3
          ja(pt+1) = row-1;
          A(pt+2) = k1+k2+k3+k4
          ja(pt+2) = row;
          A(pt+3) = -k2 * expc
          ja(pt+3) = row+numy;
          pt = pt + 4;
        else
          ! interno
          A(pt) = -k1 * expcm
          ja(pt) = row-numy;
          A(pt+1) = -k3
          ja(pt+1) = row-1;
          A(pt+2) = k1+k2+k3+k4
          ja(pt+2) = row;
          A(pt+3) = -k4
          ja(pt+3) = row+1;
          A(pt+4) = -k2 * expc
          ja(pt+4) = row+numy;
          pt = pt + 5;
        end if
      end if
      row = row + 1
    end do
  end do

  ia(row) = pt
  DEALLOCATE (bnodes)
  DEALLOCATE (anodes)

END SUBROUTINE



SUBROUTINE make_box_stiffness_2D_symmetric(S, A, ia, ja, numx, numy, xnodes, ynodes, alpha_magnetic)
  COMPLEX*16, ALLOCATABLE, INTENT(OUT) :: A(:)
  REAL*8, ALLOCATABLE, INTENT(OUT) :: S(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: ia(:), ja(:)
  INTEGER, INTENT(IN) :: numx, numy
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:)
  REAL*8, INTENT(IN) :: alpha_magnetic

  REAL*8, ALLOCATABLE :: anodes(:), bnodes(:)
  INTEGER nx, ny
  INTEGER row  ! row index into the stiffness matrix
  INTEGER pt   ! index into the A values array
  COMPLEX*16 k1, k2, k3, k4, expc

  ALLOCATE (S(numx*numy))
  ALLOCATE (A(3*numx*numy))
  ALLOCATE (ja(3*numx*numy))
  ALLOCATE (ia(numx*numy+1))
  ALLOCATE (anodes(numx+1))
  ALLOCATE (bnodes(numy+1))

  anodes = xnodes(2:numx+2) - xnodes(1:numx+1)
  bnodes = ynodes(2:numy+2) - ynodes(1:numy+1)

  pt = 1;
  row = 1;
  DO nx = 1,numx
    DO ny = 1,numy
      !row = (nx-1)*numy+ny
      S(row) = (anodes(nx)+anodes(nx+1)) * (bnodes(ny)+bnodes(ny+1)) / 4
      k1 = (bnodes(ny)+bnodes(ny+1)) / (2*anodes(nx))
      k2 = (bnodes(ny)+bnodes(ny+1)) / (2*anodes(nx+1))
      k3 = (anodes(nx)+anodes(nx+1)) / (2*bnodes(ny))
      k4 = (anodes(nx)+anodes(nx+1)) / (2*bnodes(ny+1))
      expc = exp(i*(ny-1)*alpha_magnetic);

      ia(row) = pt
      if (nx == 1) then
        if (ny == 1) then
          ! angolo alto a sinistra
          A(pt) = k1+k2+k3+k4
          ja(pt) = row
          A(pt+1) = -k4
          ja(pt+1) = row+1;
          A(pt+2) = -k2 * expc
          ja(pt+2) = row+numy;
          pt = pt + 3;
	    else if (ny == numy) then
          ! angolo basso a sinistra
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k2 * expc
          ja(pt+1) = row+numy;
          pt = pt + 2;
	else
          ! bordo sinistro
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k4
          ja(pt+1) = row+1
          A(pt+2) = -k2 * expc
          ja(pt+2) = row+numy;
          pt = pt + 3;
	end if
      else if (nx == numx) then
        if (ny == 1) then
	  ! angolo alto a destra
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k4
          ja(pt+1) = row+1;
          pt = pt + 2;
        else if (ny == numy) then
	  ! angolo basso a destra
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          pt = pt + 1;
        else
	  ! bordo destro
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k4
          ja(pt+1) = row+1;
          pt = pt + 2;
        end if
      else
        if (ny == 1) then
          ! bordo alto
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k4
          ja(pt+1) = row+1;
          A(pt+2) = -k2 * expc
          ja(pt+2) = row+numy;
          pt = pt + 3;
        else if (ny == numy) then
          ! bordo basso
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k2 * expc
          ja(pt+1) = row+numy;
          pt = pt + 2;
        else
          ! interno
          A(pt) = k1+k2+k3+k4
          ja(pt) = row;
          A(pt+1) = -k4
          ja(pt+1) = row+1;
          A(pt+2) = -k2 * expc
          ja(pt+2) = row+numy;
          pt = pt + 3;
        end if
      end if
      row = row + 1
    end do
  end do

  ia(row) = pt
  DEALLOCATE (bnodes)
  DEALLOCATE (anodes)

END SUBROUTINE




! Creates the stiffness matrix, in column-sparse format, using the box integration method
! with a non uniform orthogonal grid of size (numx+2) x (numy+2) x (numz+2)
! whose coordinates are in xnodes, ynodes and znodes.
! It uses null Dirichlet border conditions.
SUBROUTINE make_box_stiffness_3D(A, ia, ja, numx, numy, numz, xnodes, ynodes, znodes)
  COMPLEX*16, ALLOCATABLE, INTENT(OUT) :: A(:)
  INTEGER, ALLOCATABLE, INTENT(OUT) :: ia(:), ja(:)
  INTEGER, INTENT(IN) :: numx, numy, numz
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:), znodes(:)

  REAL*8, ALLOCATABLE :: anodes(:), bnodes(:), cnodes(:)
  INTEGER nx, ny, nz
  INTEGER row  ! row index into the stiffness matrix
  INTEGER pt   ! index into the A values array
  COMPLEX*16 k1, k2, k3, k4, k5, k6

  ALLOCATE (A(7*numx*numy*numz))
  ALLOCATE (ja(7*numx*numy*numz))
  ALLOCATE (ia(numx*numy*numz+1))
  ALLOCATE (anodes(numx+1))
  ALLOCATE (bnodes(numy+1))
  ALLOCATE (cnodes(numz+1))

  anodes = xnodes(2:numx+2) - xnodes(1:numx+1)
  bnodes = ynodes(2:numy+2) - ynodes(1:numy+1)
  cnodes = znodes(2:numz+2) - znodes(1:numz+1)

  pt = 1;
  row = 1;
  DO nx = 1,numx
    DO ny = 1,numy
      DO nz = 1,numz
        k1 = 2 /(anodes(nx)*(anodes(nx)+anodes(nx+1)))   ! west
        k2 = 2 /(anodes(nx+1)*(anodes(nx)+anodes(nx+1))) ! east
        k3 = 2 /(bnodes(ny)*(bnodes(ny)+bnodes(ny+1)))   ! north
        k4 = 2 /(bnodes(ny+1)*(bnodes(ny)+bnodes(ny+1))) ! south
        k5 = 2 /(bnodes(nz)*(bnodes(nz)+bnodes(nz+1)))   ! top
        k6 = 2 /(bnodes(nz+1)*(bnodes(nz)+bnodes(nz+1))) ! bottom

        ia(row) = pt
        if (nx > 1) then
          A(pt) = -k1
          ja(pt) = row - numy*numz
          pt = pt + 1
        end if
        if (ny > 1) then
          A(pt) = -k3
          ja(pt) = row - numz
          pt = pt + 1
        end if
        if (nz > 1) then
          A(pt) = -k5
          ja(pt) = row - 1;
        end if

        A(pt) = k1+k2+k3+k4+k5+k6;
        ja(pt) = row;
        pt = pt + 1;

        if (nz < numz) then
          A(pt) = -k6
          ja(pt) = row + 1;
        end if
        if (ny < numy) then
          A(pt) = -k4
          ja(pt) = row + numz
          pt = pt + 1
        end if
        if (nx < numx) then
          A(pt) = -k2
          ja(pt) = row + numy*numz
          pt = pt + 1
        end if

        row = row + 1
      end do
    end do
  end do

  ia(row+1) = pt
  DEALLOCATE (cnodes)
  DEALLOCATE (bnodes)
  DEALLOCATE (anodes)

END SUBROUTINE


END MODULE

