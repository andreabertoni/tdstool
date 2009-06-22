MODULE mod_filewrite
  USE mod_staticdata
  IMPLICIT NONE
  CONTAINS

FUNCTION STRING(inn)
IMPLICIT NONE
!  converts INTEGER "inn" in an ascii string with POS characters
INTEGER, PARAMETER :: POS= 4
INTEGER, INTENT(IN) :: inn
CHARACTER(LEN=POS) :: STRING
!.........................................................Tipi locali
INTEGER :: cifra, np, mm, num
IF (inn > (10**POS)-1) STOP "ERROR: (inn > (10**4)-1)  in STRING"
num= inn
DO np= 1, POS
  mm= pos-np
  cifra= num/(10**mm)            !  integer division
  STRING(np:np)= ACHAR(48+cifra)
  num= num - cifra*(10**mm)
END DO
END FUNCTION STRING


SUBROUTINE write_2D_in_file_matlab(path, fname, iter, psi, numx, numy)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: iter, numx, numy
  COMPLEX*16, INTENT(IN) :: psi(:)

  CHARACTER(4) :: filenum
  INTEGER nx, ny

  filenum = STRING(iter)
  open(10, FILE=TRIM(path)//"/"//fname//filenum//".dat", STATUS='REPLACE')
  DO nx = 1, numx
!    WRITE(10,*) (abs(psi(ny)), ny=(nx-1)*numy+1,nx*numy)
    DO ny = 1, numy
      IF (ny < numy) THEN
        WRITE (10, "(E12.6,A)", ADVANCE="NO") abs(psi((nx-1)*numy+ny)), " "
      ELSE
        WRITE (10, "(E12.6)") abs(psi((nx-1)*numy+ny))
      END IF
    END DO
  END DO
  close(10)

END SUBROUTINE


SUBROUTINE write_2D_cplx_in_file_gnuplot(path, fname, iter, psi, xnodes, ynodes, numx, numy, downsample_x, downsample_y)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: iter, numx, numy, downsample_x, downsample_y
  COMPLEX*16, INTENT(IN) :: psi(:)
  REAL*8 xnodes(:), ynodes(:)

  CHARACTER(4) :: filenum
  INTEGER nx, ny

  filenum = STRING(iter)
  open(10, FILE=TRIM(path)//"/"//fname//filenum//".dat", STATUS='REPLACE')
  DO nx = 1, numx, downsample_x
    DO ny = 1, numy, downsample_y
        WRITE (10, *) xnodes(nx+1), ynodes(ny+1), abs(psi((nx-1)*numy+ny))**2
    END DO
  END DO
  close(10)

END SUBROUTINE

SUBROUTINE write_2D_real_in_file_gnuplot(path, fname, iter, psi, xnodes, ynodes, numx, numy, downsample_x, downsample_y)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: iter, numx, numy, downsample_x, downsample_y
  REAL*8, INTENT(IN) :: psi(:,:)
  REAL*8 xnodes(:), ynodes(:)

  CHARACTER(4) :: filenum
  INTEGER nx, ny

  filenum = STRING(iter)
  open(10, FILE=TRIM(path)//"/"//fname//filenum//".dat", STATUS='REPLACE')
  DO nx = 1, numx, downsample_x
    DO ny = 1, numy, downsample_y
        WRITE (10, *) xnodes(nx+1), ynodes(ny+1), psi(ny, nx)
    END DO
  END DO
  close(10)

END SUBROUTINE


SUBROUTINE write_2D_cplx_in_file_bin(path, fname, iter, psi, numx, numy, downsample_x, downsample_y)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: iter, numx, numy, downsample_x, downsample_y
  COMPLEX*16, INTENT(IN) :: psi(:)
  COMPLEX*16, ALLOCATABLE :: psi_down(:)

  CHARACTER(4) :: filenum
  INTEGER nx, ny, counter

  filenum = STRING(iter)
  open(30, FILE=TRIM(path)//"/"//fname//filenum//".bin", FORM="UNFORMATTED", STATUS='REPLACE')

  if (downsample_x == 1 .and. downsample_y == 1) then
    WRITE (30) psi
  else
    ALLOCATE (psi_down((numx/downsample_x)*(numy/downsample_y)))
    counter = 1;
    do nx = 1,numx,downsample_x
      do ny = 1,numy,downsample_y
        psi_down(counter) = psi((nx-1)*numy + ny)
        counter = counter+1
      end do
    end do
    WRITE (30) psi_down;
    DEALLOCATE(psi_down)
  end if
  close(30)

END SUBROUTINE

SUBROUTINE write_2D_real_in_file_bin(path, fname, iter, psi, numx, numy, downsample_x, downsample_y)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: iter, numx, numy, downsample_x, downsample_y
  REAL*8, INTENT(IN) :: psi(:,:)
  REAL*8, ALLOCATABLE :: psi_down(:)

  CHARACTER(4) :: filenum
  INTEGER nx, ny, counter

  filenum = STRING(iter)
  open(30, FILE=TRIM(path)//"/"//fname//filenum//".bin", FORM="UNFORMATTED", STATUS='REPLACE')

  ALLOCATE (psi_down((numx/downsample_x)*(numy/downsample_y)))
  counter = 1;
  do nx = 1,numx,downsample_x
    do ny = 1,numy,downsample_y
      psi_down(counter) = psi(ny, nx)
      counter = counter+1
    end do
  end do
  WRITE (30) psi_down;
  DEALLOCATE(psi_down)
  close(30)

END SUBROUTINE



SUBROUTINE write_2D_grid(path, fname, xnodes, ynodes, numx, numy, downsample_x, downsample_y)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: numx, numy, downsample_x, downsample_y
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:)

!  CHARACTER(4) :: filenum
  INTEGER nx, ny

!  filenum = STRING(iter)
  open(10, FILE=TRIM(path)//"/"//fname//".dat", STATUS='REPLACE')
  write (10, *) numx/downsample_x, numy/downsample_y
  do nx = 2,numx+1,downsample_x
    write (10, *) xnodes(nx)
  end do
  do ny = 2,numy+1,downsample_y
    write (10, *) ynodes(ny)
  end do
  close(10)

END SUBROUTINE

SUBROUTINE write_2D_grid_bin(path, fname, xnodes, ynodes, numx, numy, downsample_x, downsample_y)
  CHARACTER(*), INTENT(IN) :: path, fname
  INTEGER, INTENT(IN) :: numx, numy, downsample_x, downsample_y
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:)
  REAL*8, ALLOCATABLE :: nodes_down(:)

  INTEGER nx, ny, pt

!  filenum = STRING(iter)
  open(10, FILE=TRIM(path)//"/"//fname//".bin", FORM="UNFORMATTED", STATUS='REPLACE')
  nx = numx/downsample_x
  ny = numy/downsample_y
  write (10) nx, ny
  ALLOCATE(nodes_down(max(nx, ny)))
  
  pt = 1
  do nx = 2,numx+1,downsample_x
    nodes_down(pt) = xnodes(nx)
    pt = pt + 1
  end do
  write(10) nodes_down
  
  pt = 1
  do ny = 2,numy+1,downsample_y
    nodes_down(pt) = ynodes(ny)
    pt = pt + 1
  end do
  write(10) nodes_down

  DEALLOCATE(nodes_down)
  close(10)

END SUBROUTINE


SUBROUTINE th_copy_file(src, dest)
  CHARACTER(*), INTENT(IN) :: src, dest
  CHARACTER(512) :: src_c, dest_c, full
  INTEGER i;
  
  if (OSV == 1) then
    src_c = REPEAT(' ', 512)
    dest_c = REPEAT(' ', 512)
    DO i = 1, LEN_TRIM(src)
      if (src(i:i) == '/') then
        src_c(i:i) = '\'
      else
        src_c(i:i) = src(i:i)
      endif
    end do
    DO i = 1, LEN_TRIM(dest)
      if (dest(i:i) == '/') then
        dest_c(i:i) = '\'
      else
        dest_c(i:i) = dest(i:i)
      endif
    end do
    
    full = TRIM("copy " // src_c) // TRIM(" " // dest_c)
    call systemqq(full)
  else
    full = TRIM("cp " // src) // TRIM(" " // dest)
    call systemqq(full)
  end if
  
END SUBROUTINE

END MODULE

