MODULE mod_indata
  USE mod_staticdata
  USE mod_strtopot
  USE IFPOSIX
  IMPLICIT NONE
  SAVE

! Simulation parameters
  REAL*8 :: dt
  REAL*8, ALLOCATABLE :: xnodes(:), ynodes(:), znodes(:)
  INTEGER :: numx, numy, numz
  INTEGER :: MAXIT
  REAL*8 :: x0, y0, z0, sigmax, sigmay, sigmaz, xenergy, yenergy, zenergy
  REAL*8 :: electronmass, mstar
  REAL*8 :: nonlin_as, magnetic

  REAL*8 :: size_x, size_y, size_z
  CHARACTER(32) :: psi_mode  ! gauss, file
  CHARACTER(32) :: pot_mode  ! zero, string, file

  REAL*8, ALLOCATABLE :: pot(:,:,:)
  COMPLEX*16, ALLOCATABLE :: psi(:)

  CHARACTER(80) :: psi_file_in   ! name of the psi (x,y) file
  CHARACTER(80) :: pot_file_in   ! name of the potential file
  CHARACTER(999) :: strpotentialX = ""
  CHARACTER(999) :: strpotentialY = ""
  CHARACTER(999) :: strpotentialZ = ""
  CHARACTER(999) :: strpotentialXYZ = ""

    ! ATTENTION: in the potential from file, the sizes pot_numx_file and pot_numy_file
    ! are the total number of nodes and not only the number of internal nodes, like in the numx and numy case
  REAL*8, ALLOCATABLE :: pot_file(:,:,:)
  REAL*8, ALLOCATABLE :: pot_x_file(:), pot_y_file(:), pot_z_file(:)
  INTEGER pot_numx_file, pot_numy_file, pot_numz_file
  INTEGER :: allow_pot_interpolation

    ! ATTENTION: in the psi from file, the sizes psi_numx_file and psi_numy_file
    ! are the total number of nodes and not only the number of internal nodes, like in the numx and numy case
  COMPLEX*16, ALLOCATABLE :: psi_file(:,:,:)
  REAL*8, ALLOCATABLE :: psi_x_file(:), psi_y_file(:), psi_z_file(:)
  INTEGER psi_numx_file, psi_numy_file, psi_numz_file

    ! Output Parameters
  CHARACTER(256) :: write_folder
  INTEGER :: write_grid   ! 0 = do not write grid, 1 = write grid
  CHARACTER(30) :: write_pot, write_psi   ! none, txt, bin, both
  REAL*8 :: write_timestep
  INTEGER :: write_downsample_x, write_downsample_y, write_downsample_z


  NAMELIST /device/ &
    & electronmass, &
    & nonlin_as,    &
    & magnetic

  NAMELIST /initial_psi/ &
    & psi_mode, &
    & x0, &
    & y0, &
    & z0, &
    & sigmax, &
    & sigmay, &
    & sigmaz, &
    & xenergy, &
    & yenergy, &
    & zenergy, &
    & psi_file_in

  NAMELIST /time/ &
    & dt, &
    & MAXIT

  NAMELIST /grid/ &
    & numx, &
    & numy, &
    & numz, &
    & size_x, &
    & size_y, &
    & size_z

  NAMELIST /potential/ &
    & pot_mode, &
    & allow_pot_interpolation, &
    & strpotentialX, &
    & strpotentialY, &
    & strpotentialZ, &
    & strpotentialXYZ, &
    & pot_file_in

  NAMELIST /write_out/ &
    & write_folder, &
    & write_grid, &
    & write_pot, &
    & write_psi, &
    & write_timestep, &
    & write_downsample_x, &
    & write_downsample_y, &
    & write_downsample_z


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
  magnetic = 0
  nonlin_as = 1e-3
  psi_mode = 'gauss'
  x0 = 1e-7
  y0 = 1e-7
  z0 = 1e-7
  sigmax = 3e-8
  sigmay = 3e-8
  sigmaz = 3e-8
  xenergy = 8e-3
  yenergy = 8e-3
  zenergy = 8e-3
  psi_file_in = ""
  dt = 2e-13
  MAXIT = 100
  numx = 128
  numy = 128
  numz = 32
  size_x = 1e-6
  size_y = 1e-6
  size_z = 0.5e-6
  allow_pot_interpolation = 1
  strpotentialX = "const  -10.e-9   100.e-9   252.8e-3"
  strpotentialY = "const  -10.e-9   100.e-9   252.8e-3"
  strpotentialZ = "const  -10.e-9   100.e-9   252.8e-3"
  strpotentialXYZ = ""
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
  INTEGER :: nx, ny, nz, grid_numx, grid_numy, grid_numz
  REAL*8, ALLOCATABLE :: potx(:), poty(:), potz(:)
  INTEGER :: lenstr, INFO, pt

  CALL PXFMKDIR(TRIM(write_folder), LEN_TRIM(write_folder), 8*8*8-1, INFO)

!*******************************************************
!    Preload potential file
!*******************************************************
  if (pot_mode == "file") then
    call READ_POT_FILE_3D(pot_file_in, INFO);
    if (INFO == 0) then
      print *, "Error reading potential file."
      STOP
    end if
  end if

!*******************************************************
!    Grid init
!*******************************************************
  ALLOCATE (pot(numz, numy,numx))
  ALLOCATE (psi(numx*numy*numz))
  ALLOCATE (xnodes(numx+2))
  ALLOCATE (ynodes(numy+2))
  ALLOCATE (znodes(numz+2))
  DO nx = 1,numx
    xnodes(nx+1) = size_x * (nx-1) / (numx-1)
  END DO
  xnodes(1) = 2*xnodes(2) - xnodes(3)
  xnodes(numx+2) = 2*xnodes(numx+1) - xnodes(numx)
  DO ny = 1,numy
    ynodes(ny+1) = size_y * (ny-1) / (numy-1)
  END DO
  ynodes(1) = 2*ynodes(2) - ynodes(3)
  ynodes(numy+2) = 2*ynodes(numy+1) - ynodes(numy)
  DO nz = 1,numz
    znodes(nz+1) = size_z * (nz-1) / (numz-1)
  END DO
  znodes(1) = 2*znodes(2) - znodes(3)
  znodes(numz+2) = 2*znodes(numz+1) - znodes(numz)

!*******************************************************
!    Potential init
!*******************************************************
  pot(:,:,:) = 0;
  if (pot_mode == "string") then
    ALLOCATE (potx(1:numx))
    ALLOCATE (poty(1:numy))
    ALLOCATE (potz(1:numz))
    call strtopot1D(1, mstar, strpotentialX, numx, xnodes, potx)
    call strtopot1D(1, mstar, strpotentialY, numy, ynodes, poty)
    call strtopot1D(1, mstar, strpotentialZ, numz, znodes, potz)
    call strtopot3D(1, mstar, strpotentialXYZ, numx, numy, numz, xnodes, ynodes, znodes, pot)
    DO nx = 1, numx
      DO ny = 1, numy
        DO nz = 1, numz
          pot(nz, ny, nx) = (pot(nz, ny, nx) + potx(nx) + poty(ny) + potz(nz))*ELCH
        END DO
      END DO
    END DO
    DEALLOCATE (potz)
    DEALLOCATE (poty)
    DEALLOCATE (potx)
  else if (pot_mode == "file") then
      ! check if resampling i really needed
    if ((pot_numx_file-2) == numx .and. (pot_numy_file-2) == numy) then
      DO nx = 1, numx
        DO ny = 1, numy
          DO nz = 1, numz
            pot(nz, ny, nx)= pot_file(nz+1, ny+1, nx+1)
          END DO
        END DO
      END DO
    else
      print *, "Potential grid is different than TDS tool grid."
      STOP
    end if

    DEALLOCATE (pot_file)
    DEALLOCATE (pot_x_file)
    DEALLOCATE (pot_y_file)
    DEALLOCATE (pot_z_file)
  end if


!*******************************************************
!    Wave function init
!*******************************************************
  if (psi_mode == "gauss") then

    call create_gaussian_packet_3D(psi, numx, numy, numz, xnodes, ynodes, znodes, x0, y0, z0, sigmax, sigmay, sigmaz, xenergy, yenergy, zenergy)

  else if (psi_mode == "file") then

    call READ_PSI_FILE_3D(psi_file_in, INFO);
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
        DO nz = 1, numz
          psi(pt)= psi_file(nz+1, ny+1, nx+1)
          pt = pt + 1
        END DO
      END DO
    END DO
    DEALLOCATE (psi_file)
    DEALLOCATE (psi_x_file)
    DEALLOCATE (psi_y_file)
    DEALLOCATE (psi_z_file)
  end if

END SUBROUTINE INDATA_COMPUTE

SUBROUTINE create_gaussian_packet_3D(psi, numx, numy, numz, xnodes, ynodes, znodes, x0, y0, z0, sigmax, sigmay, sigmaz, xenergy, yenergy, zenergy)
  COMPLEX*16, INTENT(OUT) :: psi(:)
  REAL*8, INTENT(IN) :: xnodes(:), ynodes(:), znodes(:)
  REAL*8, INTENT(IN) :: x0, y0, z0, sigmax, sigmay, sigmaz, xenergy, yenergy, zenergy
  INTEGER, INTENT(IN) :: numx, numy, numz

  REAL*8 k0x, k0y, k0z, norma
  INTEGER nx, ny, nz, pt

  k0x = sqrt(2*mstar*ELCH*xenergy) / HBAR
  k0y = sqrt(2*mstar*ELCH*yenergy) / HBAR
  k0z = sqrt(2*mstar*ELCH*zenergy) / HBAR

  norma = 1. / (sqrt(sigmax * sqrt(2.*PIG)) * sqrt(sigmay * sqrt(2.*PIG)) * sqrt(sigmaz * sqrt(2.*PIG)))

  pt = 1
  DO nx = 2, numx+1
    DO ny = 2, numy+1
      DO nz = 2, numz+1
        psi(pt) = norma * exp(-(((xnodes(nx)-x0)/(2.*sigmax))**2) + i*k0x*(xnodes(nx)-x0)) * &
                      &   exp(-(((ynodes(ny)-y0)/(2.*sigmay))**2) + i*k0y*(ynodes(ny)-y0)) * &
                      &   exp(-(((znodes(nz)-z0)/(2.*sigmaz))**2) + i*k0z*(znodes(nz)-z0))
        pt = pt + 1
      end do
    end do
  end do
END SUBROUTINE

SUBROUTINE READ_POT_FILE_3D(fname, INFO)
  CHARACTER(80), INTENT(IN) :: fname
  INTEGER, INTENT(OUT) :: INFO
  INTEGER nx, ny, nz
  REAL*8 x, y, z, val
  CHARACTER(512) :: line, str1

  INFO = 0
  OPEN(22, FILE=TRIM(fname), FORM="FORMATTED", STATUS="OLD")
  do
    read (22, '(A)', end=999) line
    if (line(1:5) == "#TDS " .or. line(1:5) == "#tds ") then
      READ(line, *) str1, pot_numx_file, pot_numy_file, pot_numz_file
      EXIT
    end if
  end do

  ALLOCATE (pot_file(pot_numz_file, pot_numy_file, pot_numx_file))
  ALLOCATE (pot_x_file(pot_numx_file))
  ALLOCATE (pot_y_file(pot_numy_file))
  ALLOCATE (pot_z_file(pot_numz_file))

  DO nx = 1, pot_numx_file
  DO ny = 1, pot_numy_file
  DO nz = 1, pot_numz_file
    do
      read (22, '(A)', end=999) line
      if (line(1:1) /= "#") then
        EXIT
      end if
    end do
    READ (line, *, end=999), x, y, z, val
    if (nx == 1 .and. ny == 1) then
      pot_z_file(nz) = z
    else
      if (pot_z_file(nz) /= z) then
        PRINT *, "Not consistent potential grid"
        STOP
      end if
    end if

    if (nx == 1 .and. nz == 1) then
      pot_y_file(ny) = y
    else
      if (pot_y_file(ny) /= y) then
        PRINT *, "Not consistent potential grid"
        STOP
      end if
    end if

    if (ny == 1 .and. nz == 1) then
      pot_x_file(nx) = x
    else
      if (pot_x_file(nx) /= x) then
        PRINT *, "Not consistent potential grid"
        STOP
      end if
    end if

    pot_file(nz, ny, nx) = val
  END DO
  END DO
  END DO
  INFO = 1

999 CLOSE(22)
END SUBROUTINE


SUBROUTINE READ_PSI_FILE_3D(fname, INFO)
  CHARACTER(80), INTENT(IN) :: fname
  INTEGER, INTENT(OUT) :: INFO
  INTEGER nx, ny, nz
  REAL*8 x, y, z, val
  CHARACTER(512) :: line, str1

  INFO = 0
  OPEN(22, FILE=TRIM(fname), FORM="FORMATTED", STATUS="OLD")
  do
    read (22, '(A)', end=999) line
    if (line(1:5) == "#TDS " .or. line(1:5) == "#tds ") then
      READ(line, *) str1, psi_numx_file, psi_numy_file, psi_numz_file
      EXIT
    end if
  end do

  ALLOCATE (psi_file(psi_numz_file, psi_numy_file, psi_numx_file))
  ALLOCATE (psi_x_file(psi_numx_file))
  ALLOCATE (psi_y_file(psi_numy_file))
  ALLOCATE (psi_z_file(psi_numz_file))

  DO nx = 1, psi_numx_file
  DO ny = 1, psi_numy_file
  DO nz = 1, psi_numz_file
    do
      read (22, '(A)', end=999) line
      if (line(1:1) /= "#") then
        EXIT
      end if
    end do
    READ (line, *, end=999), x, y, z, val

    if (nx == 1 .and. ny == 1) then
      psi_z_file(nz) = z
    else
      if (psi_z_file(nz) /= z) then
        PRINT *, "Not consistent psi grid"
        STOP
      end if
    end if

    if (nx == 1 .and. nz == 1) then
      psi_y_file(ny) = y
    else
      if (psi_y_file(ny) /= y) then
        PRINT *, "Not consistent psi grid"
        STOP
      end if
    end if

    if (ny == 1 .and. nz == 1) then
      psi_x_file(nx) = x
    else
      if (psi_x_file(nx) /= x) then
        PRINT *, "Not consistent psi grid"
        STOP
      end if
    end if

    psi_file(nz, ny, nx) = val
  END DO
  END DO
  END DO
  INFO = 1

999 CLOSE(22)
END SUBROUTINE

END MODULE mod_indata

