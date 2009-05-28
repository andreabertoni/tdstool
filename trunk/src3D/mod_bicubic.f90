MODULE mod_bicubic
  IMPLICIT NONE
  SAVE

CONTAINS

SUBROUTINE bcucof(y, y1, y2, y12, d1, d2, c)
  REAL*8, INTENT(IN) :: y(4), y1(4), y2(4), y12(4), d1, d2
  REAL*8, INTENT(OUT) :: c(4,4)

  INTEGER i, j, k, l
  REAL*8 d1d2, xx, cl(16), wt(16,16), x(16)
  SAVE wt
  DATA wt/1, 0, -3, 2, 4*0, -3, 0, 9, -6, 2, 0, -6, 4, 8*0, 3, 0, -9, 6, -2, 0, 6, -4, &
   10*0, 9, -6, 2*0, -6, 4, 2*0, 3, -2, 6*0, -9, 6, 2*0, 6, -4, &
   4*0, 1, 0, -3, 2, -2, 0, 6, -4, 1, 0, -3, 2, 8*0, -1, 0, 3, -2, 1, 0, -3, 2, &
   10*0, -3, 2, 2*0, 3, -2, 6*0, 3, -2, 2*0, -6, 4, 2*0, 3, -2, &
   0, 1, -2, 1, 5*0, -3, 6, -3, 0, 2, -4, 2, 9*0, 3, -6, 3, 0, -2, 4, -2, &
   10*0, -3, 3, 2*0, 2, -2, 2*0, -1, 1, 6*0, 3, -3, 2*0, -2, 2, &
    5*0, 1, -2, 1, 0, -2, 4, -2, 0, 1, -2, 1, 9*0, -1, 2, -1, 0, 1, -2, 1, &
   10*0, 1, -1, 2*0, -1, 1, 6*0, -1, 1, 2*0, 2, -2, 2*0, -1, 1 /

  d1d2 = d1*d2
  do i = 1,4
    x(i)=y(i)
    x(i+4) = y1(i)*d1
    x(i+8) = y2(i)*d2
    x(i+12) = y12(i)*d1d2
  end do

  do i = 1,16
    xx = 0.
    do k = 1,16
      xx = xx + wt(i, k)*x(k)
    end do
    cl(i) = xx
  end do

  l = 0
  do i = 1,4
    do j = 1,4
      l = l + 1
      c(i,j)=cl(l)
    end do
  end do
  return
END SUBROUTINE

SUBROUTINE bicubic(y, x1nod, x2nod, x1, x2, ansy)
  REAL*8, INTENT(IN) :: y(4, 4), x1nod(4), x2nod(4), x1, x2
  REAL*8, INTENT(OUT) :: ansy

  INTEGER i, j, k
  REAL*8 t, u, c(4,4), ya(4), y1a(4), y2a(4), y12a(4)

  ! Estimate derivatives
  j = 2
  k = 2
  ya(1) = y(j, k)
  y1a(1) = (y(j,k+1)-y(j,k-1)) / (x1nod(k+1)-x1nod(k-1))
  y2a(1) = (y(j+1,k)-y(j-1,k)) / (x2nod(j+1)-x2nod(j-1))
  y12a(1) = (y(j+1,k+1)-y(j+1,k-1)-y(j-1,k+1)+y(j-1,k-1)) / ((x1nod(k+1)-x1nod(k-1))*(x2nod(j+1)-x2nod(j-1)))
  j = 2
  k = 3
  ya(2) = y(j, k)
  y1a(2) = (y(j,k+1)-y(j,k-1)) / (x1nod(k+1)-x1nod(k-1))
  y2a(2) = (y(j+1,k)-y(j-1,k)) / (x2nod(j+1)-x2nod(j-1))
  y12a(2) = (y(j+1,k+1)-y(j+1,k-1)-y(j-1,k+1)+y(j-1,k-1)) / ((x1nod(k+1)-x1nod(k-1))*(x2nod(j+1)-x2nod(j-1)))
  j = 3
  k = 3
  ya(3) = y(j, k)
  y1a(3) = (y(j,k+1)-y(j,k-1)) / (x1nod(k+1)-x1nod(k-1))
  y2a(3) = (y(j+1,k)-y(j-1,k)) / (x2nod(j+1)-x2nod(j-1))
  y12a(3) = (y(j+1,k+1)-y(j+1,k-1)-y(j-1,k+1)+y(j-1,k-1)) / ((x1nod(k+1)-x1nod(k-1))*(x2nod(j+1)-x2nod(j-1)))
  j = 3
  k = 2
  ya(4) = y(j, k)
  y1a(4) = (y(j,k+1)-y(j,k-1)) / (x1nod(k+1)-x1nod(k-1))
  y2a(4) = (y(j+1,k)-y(j-1,k)) / (x2nod(j+1)-x2nod(j-1))
  y12a(4) = (y(j+1,k+1)-y(j+1,k-1)-y(j-1,k+1)+y(j-1,k-1)) / ((x1nod(k+1)-x1nod(k-1))*(x2nod(j+1)-x2nod(j-1)))

  CALL bcucof(ya, y1a, y2a, y12a, x1nod(3)-x1nod(2), x2nod(3)-x2nod(2), c)
  u = (x1-x1nod(2)) / (x1nod(3)-x1nod(2))
  t = (x2-x2nod(2)) / (x2nod(3)-x2nod(2))
  ansy = 0
  do i = 4,1,-1
    ansy = t*ansy + ((c(i,4)*u + c(i,3))*u + c(i,2))*u + c(i, 1)
  end do

END SUBROUTINE


! all the nodes should be ordered ascendingly
SUBROUTINE resample_matrix_2D(A, xnodes, ynodes, numx, numy, B, new_xnodes, new_ynodes, new_numx, new_numy)
  REAL*8, INTENT(IN) :: A(:,:), xnodes(:), ynodes(:), new_xnodes(:), new_ynodes(:)
  REAL*8, INTENT(OUT) :: B(:,:)
  INTEGER, INTENT(IN) :: numx, numy, new_numx, new_numy

  INTEGER nx, ny, nxx, nyy, idxx, idxy
  REAL*8 y(4, 4), ansy

  idxx = 1
  DO nx = 1, new_numx
    DO WHILE ((idxx <= numx) .and. (xnodes(idxx) <= new_xnodes(nx)))
      idxx = idxx + 1
    END DO
    idxx = idxx - 1
    if (idxx > 1 .and. idxx < (numx-1)) then
      idxy = 1
      DO ny = 1, new_numy
        DO WHILE ((idxy <= numy) .and. (ynodes(idxy) <= new_ynodes(ny)))
          idxy = idxy + 1
        END DO
        idxy = idxy - 1
        if (idxy > 1 .and. idxy < (numy-1)) then
          DO nxx = idxx-1, idxx+2
            DO nyy = idxy-1, idxy+2
              y(nyy-idxy+2, nxx-idxx+2) = A(nyy, nxx)
            END DO
          END DO
          CALL bicubic(y, xnodes(idxx-1:idxx+2), ynodes(idxy-1:idxy+2), new_xnodes(nx), new_ynodes(ny), ansy)
          B(ny, nx) = ansy
        else
          B(ny, nx) = 0
        end if
      END DO
    else
      DO ny = 1, new_numy
        B(ny, nx) = 0
      END DO
    end if
  END DO

END SUBROUTINE



SUBROUTINE resample_matrix_1D(A, xnodes, numx, B, new_xnodes, new_numx)
  REAL*8, INTENT(IN) :: A(:), xnodes(:), new_xnodes(:)
  REAL*8, INTENT(OUT) :: B(:)
  INTEGER, INTENT(IN) :: numx, new_numx

  INTEGER :: idxx, nx, i, j
  REAL*8 :: v

  idxx = 1
  DO nx = 1, new_numx
    DO WHILE ((idxx <= numx) .and. (xnodes(idxx) <= new_xnodes(nx)))
      idxx = idxx + 1
    END DO
    idxx = idxx - 1
    if (idxx == 1 .or. idxx == (numx-1)) then
      ! Use Linear interpolation
      B(nx) = ((new_xnodes(nx)-xnodes(idxx+1))*A(idxx) - (new_xnodes(nx)-xnodes(idxx))*A(idxx+1)) / (xnodes(idxx)-xnodes(idxx+1))
    elseif (idxx > 1 .and. idxx < (numx-1)) then
      ! Use cubic interpolation
      B(nx) = 0
      do i = -1, 2
        v = 1
        do j = -1, 2
          if (j /= i) then
            v = v * (new_xnodes(nx)-xnodes(idxx+j)) / (xnodes(idxx+i)-xnodes(idxx+j))
          end if
        end do
        B(nx) = v * A(idxx+i)
      end do
    else
      B(nx) = 0
    end if
  END DO

END SUBROUTINE



SUBROUTINE Integrate_2D(integ, psi_in, NUMX, NUMY, xnodes, ynodes, intx0, inty0, intx1, inty1)
  REAL*8, INTENT(OUT) :: integ
  REAL*8, INTENT(IN) :: psi_in(:,:)
  INTEGER, INTENT(IN) :: NUMX, NUMY
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:)
  REAL*8, INTENT(IN) :: intx0, inty0, intx1, inty1

  INTEGER :: idxx0, idxy0, idxx1, idxy1, nx, ny
  REAL*8 :: vsum, w, h

  idxx0 = 1
  DO WHILE ((idxx0 < numx) .and. (xnodes(idxx0) < intx0))
    idxx0 = idxx0 + 1
  END DO
  idxy0 = 1
  DO WHILE ((idxy0 < numy) .and. (ynodes(idxy0) < inty0))
    idxy0 = idxy0 + 1
  END DO
  idxx1 = 1
  DO WHILE ((idxx1 < numx) .and. (xnodes(idxx1) < intx1))
    idxx1 = idxx1 + 1
  END DO
  idxy1 = 1
  DO WHILE ((idxy1 < numy) .and. (ynodes(idxy1) < inty1))
    idxy1 = idxy1 + 1
  END DO

  integ = 0
  do nx = idxx0, idxx1
    if (nx < numx) then
      w = xnodes(nx+1)-xnodes(nx)
    else
      w = xnodes(nx)-xnodes(nx-1)
    endif
    vsum = 0
    do ny = idxy0, idxy1
      if (ny < numy) then
        h = ynodes(ny+1)-ynodes(ny)
      else
        h = ynodes(ny)-ynodes(ny-1)
      endif
      vsum = vsum + psi_in(ny, nx) * w * h
    end do
    integ = integ + vsum
  end do

END SUBROUTINE


END MODULE mod_bicubic

