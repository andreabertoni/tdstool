MODULE mod_indata
  USE mod_staticdata
  USE mod_strtopot
  USE mod_bicubic
  USE IFPOSIX
  IMPLICIT NONE
  SAVE

! Simulation parameters
  REAL*8 :: dt
  REAL*8, ALLOCATABLE :: xnodes(:), ynodes(:)
  INTEGER :: numx, numy
  INTEGER :: MAXIT
  REAL*8 :: x0, y0, sigmax, sigmay, xenergy, yenergy
  REAL*8 :: electronmass, mstar

  REAL*8 :: size_x, size_y
  CHARACTER(32) :: psi_mode  ! gauss, file
  CHARACTER(32) :: grid_mode ! uniform, file, estimate, pot, adaptive
  CHARACTER(32) :: pot_mode  ! zero, string, file

  REAL*8, ALLOCATABLE :: pot(:,:)
  COMPLEX*16, ALLOCATABLE :: psi(:)
  COMPLEX*16, ALLOCATABLE :: A(:), S(:)
  INTEGER, ALLOCATABLE :: ia(:), ja(:)

  CHARACTER(80) :: psi_file_in   ! name of the psi (x,y) file
  CHARACTER(80) :: grid_file_in  ! name of the xnodes and ynodes file
  CHARACTER(80) :: pot_file_in   ! name of the potential file
  CHARACTER(999) :: strpotentialX = ""
  CHARACTER(999) :: strpotentialY = ""
  CHARACTER(999) :: strpotentialXY = ""

    ! ATTENTION: in the potential from file, the sizes pot_numx_file and pot_numy_file
    ! are the total number of nodes and not only the number of internal nodes, like in the numx and numy case
  REAL*8, ALLOCATABLE :: pot_file(:,:)
  REAL*8, ALLOCATABLE :: pot_x_file(:), pot_y_file(:)
  INTEGER pot_numx_file, pot_numy_file
  INTEGER :: allow_pot_interpolation

    ! ATTENTION: in the psi from file, the sizes psi_numx_file and psi_numy_file
    ! are the total number of nodes and not only the number of internal nodes, like in the numx and numy case
  COMPLEX*16, ALLOCATABLE :: psi_file(:,:)
  REAL*8, ALLOCATABLE :: psi_x_file(:), psi_y_file(:)
  INTEGER psi_numx_file, psi_numy_file

    ! Output Parameters
  CHARACTER(256) :: write_folder
  INTEGER :: write_grid   ! 0 = do not write grid, 1 = write grid
  CHARACTER(30) :: write_pot, write_psi   ! none, txt, bin, both
  REAL*8 :: write_timestep
  INTEGER :: write_downsample_x, write_downsample_y


  NAMELIST /device/ &
    & electronmass

  NAMELIST /initial_psi/ &
    & psi_mode, &
    & x0, &
    & y0, &
    & sigmax, &
    & sigmay, &
    & xenergy, &
    & yenergy, &
    & psi_file_in

  NAMELIST /time/ &
    & dt, &
    & MAXIT

  NAMELIST /grid/ &
    & grid_mode, &
    & numx, &
    & numy, &
    & size_x, &
    & size_y, &
    & grid_file_in

  NAMELIST /potential/ &
    & pot_mode, &
    & allow_pot_interpolation, &
    & strpotentialX, &
    & strpotentialY, &
    & strpotentialXY, &
    & pot_file_in

  NAMELIST /write_out/ &
    & write_folder, &
    & write_grid, &
    & write_pot, &
    & write_psi, &
    & write_timestep, &
    & write_downsample_x, &
    & write_downsample_y


CONTAINS

!===================================================================
SUBROUTINE INDATA_GET(nmlfile, INFO)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: nmlfile
  INTEGER, INTENT(OUT) :: INFO

  LOGICAL :: file_exists
  INQUIRE(FILE=nmlfile, EXIST=file_exists)

  if (file_exists) then
    OPEN(33, FILE=TRIM(nmlfile), FORM="FORMATTED", ACTION="READ",  &
         &   STATUS="OLD")
    READ(33,NML=device)
    READ(33,NML=initial_psi)
    READ(33,NML=time)
    READ(33,NML=grid)
    READ(33,NML=potential)
    READ(33,NML=write_out)
    CLOSE(33)
    mstar = electronmass * ELMASS0
    INFO = 0
  else
    INFO = 1
  end if

END SUBROUTINE INDATA_GET

SUBROUTINE INDATA_FILL_WITH_DEFAULT
  electronmass = 0.067
  psi_mode = 'gauss'
  x0 = 1e-7
  y0 = 1e-7
  sigmax = 3e-8
  sigmay = 3e-8
  xenergy = 8e-3
  yenergy = 8e-3
  psi_file_in = ""
  dt = 2e-13
  MAXIT = 100
  grid_mode = 'pot'
  numx = 100
  numy = 100
  size_x = 1e-6
  size_y = 1e-6
  grid_file_in = ""
  pot_mode = "file"
  allow_pot_interpolation = 1
  strpotentialX = "const  -10.e-9   100.e-9   252.8e-3"
  strpotentialY = "const  -10.e-9   100.e-9   252.8e-3"
  pot_file_in = "pot.dat"
  write_folder = "out"
  write_grid = 1
  write_pot = "both"
  write_psi = "both"
  write_timestep = 2e-12
  write_downsample_x = 1
  write_downsample_y = 1
END SUBROUTINE

SUBROUTINE INDATA_SAVE(nmlfile)
  IMPLICIT NONE
  CHARACTER(*), INTENT(IN) :: nmlfile

  OPEN(33, FILE=TRIM(nmlfile), FORM="FORMATTED", DELIM='quote', STATUS='REPLACE')
  WRITE(33,NML=device)
  WRITE(33,NML=initial_psi)
  WRITE(33,NML=time)
  WRITE(33,NML=grid)
  WRITE(33,NML=potential)
  WRITE(33,NML=write_out)
  CLOSE(33)

END SUBROUTINE INDATA_SAVE


!===================================================================
SUBROUTINE INDATA_COMPUTE
  INTEGER :: nx, ny, grid_numx, grid_numy
  REAL*8, ALLOCATABLE :: potx(:), poty(:)
  INTEGER :: lenstr, INFO, pt

  CALL PXFMKDIR(TRIM(write_folder), LEN_TRIM(write_folder), 8*8*8-1, INFO)
print *, INFO

!*******************************************************
!    Preload potential file
!*******************************************************
  if (pot_mode == "file") then
    call READ_POT_FILE(pot_file_in, INFO);
    if (INFO == 0) then
      print *, "Error reading potential file."
      STOP
    end if
  end if

!*******************************************************
!    Grid init
!*******************************************************
  if (grid_mode == "file") then
    lenstr = LEN_TRIM(grid_file_in)
    OPEN(22, FILE=grid_file_in, FORM="FORMATTED", STATUS="OLD")
    READ(22) grid_numx, grid_numy
    if (numx /= grid_numx .or. numy /= grid_numy) then
      print *, "Warning: grid size in namelist differs from grid size in file "//grid_file_in
      print *, "Taking sizes in the file."
      numx = grid_numx;
      numy = grid_numy;
    end if
    ALLOCATE (pot(numy,numx))
    ALLOCATE (psi(numx*numy))
    ALLOCATE (xnodes(numx+2))
    ALLOCATE (ynodes(numy+2))
    DO nx = 1,numx+2
      READ(22,*) xnodes(nx)
    END DO
    DO ny = 1,numy+2
      READ(22,*) ynodes(ny)
    END DO
    CLOSE(22)

  else if (grid_mode == "pot") then
    numx = pot_numx_file - 2
    numy = pot_numy_file - 2
    ALLOCATE (pot(numy,numx))
    ALLOCATE (psi(numx*numy))
    ALLOCATE (xnodes(numx+2))
    ALLOCATE (ynodes(numy+2))
    xnodes = pot_x_file
    ynodes = pot_y_file

  else  ! uniform
    ALLOCATE (pot(numy,numx))
    ALLOCATE (psi(numx*numy))
    ALLOCATE (xnodes(numx+2))
    ALLOCATE (ynodes(numy+2))
    DO nx = 1,numx+2
      xnodes(nx) = size_x * (nx-1) / (numx+1)
    END DO
    DO ny = 1,numy+2
      ynodes(ny) = size_y * (ny-1) / (numy+1)
    END DO
  end if


!*******************************************************
!    Potential init
!*******************************************************
  pot(:,:) = 0;
  if (pot_mode == "string") then
    ALLOCATE (potx(1:numx))
    ALLOCATE (poty(1:numy))
    call strtopot1D(1, mstar, strpotentialX, numx, xnodes, potx)
    call strtopot1D(1, mstar, strpotentialY, numy, ynodes, poty)
    call strtopot2D(1, mstar, strpotentialXY, numx, numy, xnodes, ynodes, pot)
    DO nx = 1, numx
      DO ny = 1, numy
        pot(ny, nx) = (pot(ny, nx) + potx(nx) + poty(ny))*ELCH
      END DO
    END DO
    DEALLOCATE (poty)
    DEALLOCATE (potx)
  else if (pot_mode == "file") then
    if (grid_mode == "pot") then
      DO nx = 1, numx
        DO ny = 1, numy
          pot(ny, nx)= pot_file(ny+1, nx+1)
        END DO
      END DO
    else
      ! check if resampling i really needed
      if ((pot_numx_file-2) == numx .and. (pot_numy_file-2) == numy) then
        DO nx = 1, numx
          DO ny = 1, numy
            pot(ny, nx)= pot_file(ny+1, nx+1)
          END DO
        END DO
      else if (allow_pot_interpolation > 0) then
        CALL resample_matrix_2D(pot_file, pot_x_file, pot_y_file, pot_numx_file, pot_numy_file, pot, xnodes(2:numx), ynodes(2:numy), numx, numy)
      else
        print *, "Potential interpolation needed, but not allowed by user"
        STOP
      end if
    end if

    DEALLOCATE (pot_file)
    DEALLOCATE (pot_x_file)
    DEALLOCATE (pot_y_file)
  end if


!*******************************************************
!    Wave function init
!*******************************************************
  if (psi_mode == "gauss") then

    call create_gaussian_packet(psi, numx, numy, xnodes, ynodes, x0, y0, sigmax, sigmay, xenergy, yenergy)

  else if (psi_mode == "file") then

    call READ_PSI_FILE(psi_file_in, INFO);
    if (INFO == 0) then
      print *, "Error reading psi file."
      STOP
    end if

    if ((psi_numx_file - 2) /= numx .or. (psi_numy_file - 2) /= numy) then
      print *, "psi file grid not consistent with problem grid"
      STOP
    end if

    pt = 1;
    DO nx = 1, numx
      DO ny = 1, numy
        psi(pt)= psi_file(ny+1, nx+1)
        pt = pt + 1
      END DO
    END DO
    DEALLOCATE (psi_file)
    DEALLOCATE (psi_x_file)
    DEALLOCATE (psi_y_file)
  end if

END SUBROUTINE INDATA_COMPUTE

SUBROUTINE create_gaussian_packet(psi, numx, numy, xnodes, ynodes, x0, y0, sigmax, sigmay, xenergy, yenergy)
  COMPLEX*16, INTENT(OUT) :: psi(:)
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:)
  REAL*8, INTENT(IN) :: x0, y0, sigmax, sigmay, xenergy, yenergy
  INTEGER, INTENT(IN) :: numx, numy

  REAL*8 k0x, k0y, norma
  INTEGER nx, ny, pt

  k0x = sqrt(2*mstar*ELCH*xenergy) / HBAR
  k0y = sqrt(2*mstar*ELCH*yenergy) / HBAR

  norma = 1. / (sqrt(sigmax * sqrt(2.*PIG)) * sqrt(sigmay * sqrt(2.*PIG)))

  pt = 1
  DO nx = 2, numx+1
    DO ny = 2, numy+1
      psi(pt) = norma * exp(-(((xnodes(nx)-x0)/(2.*sigmax))**2) + i*k0x*(xnodes(nx)-x0)) * exp(-(((ynodes(ny)-y0)/(2.*sigmay))**2) + i*k0y*(ynodes(ny)-y0))
    pt = pt + 1
    end do
  end do
END SUBROUTINE

SUBROUTINE READ_POT_FILE(fname, INFO)
  CHARACTER(80), INTENT(IN) :: fname
  INTEGER, INTENT(OUT) :: INFO
  INTEGER nx, ny
  REAL*8 x, y, val
  CHARACTER(512) :: line, str1

  INFO = 0
  OPEN(22, FILE=TRIM(fname), FORM="FORMATTED", STATUS="OLD")
  do
    read (22, '(A)', end=999) line
    if (line(1:5) == "#TDS " .or. line(1:5) == "#tds ") then
      READ(line, *) str1, pot_numx_file, pot_numy_file
      EXIT
    end if
  end do

  ALLOCATE (pot_file(pot_numy_file, pot_numx_file))
  ALLOCATE (pot_x_file(pot_numx_file))
  ALLOCATE (pot_y_file(pot_numy_file))

  DO nx = 1, pot_numx_file
  DO ny = 1, pot_numy_file
    do
      read (22, '(A)', end=999) line
      if (line(1:1) /= "#") then
        EXIT
      end if
    end do
    READ (line, *, end=999), x, y, val
    if (nx == 1) then
      pot_y_file(ny) = y
    else
      if (pot_y_file(ny) /= y) then
        PRINT *, "Not consistent potential grid"
        STOP
      end if
    end if

    if (ny == 1) then
      pot_x_file(nx) = x
    else
      if (pot_x_file(nx) /= x) then
        PRINT *, "Not consistent potential grid"
        STOP
      end if
    end if

    pot_file(ny, nx) = val
  END DO
  END DO
  INFO = 1

999 CLOSE(22)
END SUBROUTINE


SUBROUTINE READ_PSI_FILE(fname, INFO)
  CHARACTER(80), INTENT(IN) :: fname
  INTEGER, INTENT(OUT) :: INFO
  INTEGER nx, ny
  REAL*8 x, y, val
  CHARACTER(512) :: line, str1

  INFO = 0
  OPEN(22, FILE=TRIM(fname), FORM="FORMATTED", STATUS="OLD")
  do
    read (22, '(A)', end=999) line
    if (line(1:5) == "#TDS " .or. line(1:5) == "#tds ") then
      READ(line, *) str1, psi_numx_file, psi_numy_file
      EXIT
    end if
  end do

  ALLOCATE (psi_file(psi_numy_file, psi_numx_file))
  ALLOCATE (psi_x_file(psi_numx_file))
  ALLOCATE (psi_y_file(psi_numy_file))

  DO nx = 1, psi_numx_file
  DO ny = 1, psi_numy_file
    do
      read (22, '(A)', end=999) line
      if (line(1:1) /= "#") then
        EXIT
      end if
    end do
    READ (line, *, end=999), x, y, val
    if (nx == 1) then
      psi_y_file(ny) = y
    else
      if (psi_y_file(ny) /= y) then
        PRINT *, "Not consistent psi grid"
        STOP
      end if
    end if

    if (ny == 1) then
      psi_x_file(nx) = x
    else
      if (psi_x_file(nx) /= x) then
        PRINT *, "Not consistent psi grid"
        STOP
      end if
    end if

    psi_file(ny, nx) = val
  END DO
  END DO
  INFO = 1

999 CLOSE(22)
END SUBROUTINE

END MODULE mod_indata

