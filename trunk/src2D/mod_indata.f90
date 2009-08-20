MODULE mod_indata
  USE mod_staticdata
  USE mod_strtopot
  USE mod_bicubic
  USE mod_filewrite
  USE IFPOSIX
  IMPLICIT NONE
  SAVE

  TYPE filelist_struct
    REAL*8 :: time
    CHARACTER(256) :: file
  END TYPE filelist_struct
  
! Simulation parameters
  REAL*8 :: dt
  REAL*8, ALLOCATABLE :: xnodes(:), ynodes(:)
  INTEGER :: numx, numy
  INTEGER :: MAXIT
  REAL*8 :: x0, y0, sigmax, sigmay, xenergy, yenergy
  REAL*8 :: electronmass, mstar
  REAL*8 :: nonlin_as, magnetic

  REAL*8 :: size_x, size_y
  CHARACTER(32) :: psi_mode  ! gauss, file
  CHARACTER(32) :: grid_mode ! uniform, file, estimate, pot, adaptive

  REAL*8, ALLOCATABLE :: pot(:,:)
  COMPLEX*16, ALLOCATABLE :: psi(:)
  COMPLEX*16, ALLOCATABLE :: A(:)
  REAL*8, ALLOCATABLE :: S(:)
  INTEGER, ALLOCATABLE :: ia(:), ja(:)

  CHARACTER(80) :: psi_file_in   ! name of the psi (x,y) file
  CHARACTER(80) :: grid_file_in  ! name of the xnodes and ynodes file
  CHARACTER(80) :: pot_filelist_name  ! name of the file containing the potential file list
  CHARACTER(80) :: pot_file_in   ! name of the potential file
  CHARACTER(999) :: strpotentialX = ""
  CHARACTER(999) :: strpotentialY = ""
  CHARACTER(999) :: strpotentialXY = ""

    ! ATTENTION: in the potential from file, the sizes pot_numx_file and pot_numy_file
    ! are the total number of nodes and not only the number of internal nodes, like in the numx and numy case
  REAL*8, ALLOCATABLE :: pot_file_raw(:,:)
  REAL*8, ALLOCATABLE :: pot_x_file(:), pot_y_file(:)
  INTEGER pot_numx_file, pot_numy_file
  INTEGER :: allow_pot_interpolation

  REAL*8, ALLOCATABLE :: pot_file_static(:,:), pot_filelist(:,:)
  
  TYPE(filelist_struct), ALLOCATABLE :: filelist_data(:)
  INTEGER :: filelist_count;

    ! ATTENTION: in the psi from file, the sizes psi_numx_file and psi_numy_file
    ! are the total number of nodes and not only the number of internal nodes, like in the numx and numy case
  COMPLEX*16, ALLOCATABLE :: psi_file(:,:)
  REAL*8, ALLOCATABLE :: psi_x_file(:), psi_y_file(:)
  INTEGER psi_numx_file, psi_numy_file

    ! Output Parameters
  CHARACTER(512) :: write_folder, run_name, nmlfile_name
  CHARACTER(30) :: write_grid, write_pot, write_psi   ! none, txt, bin, both
  REAL*8 :: write_timestep
  INTEGER :: write_downsample_x, write_downsample_y

  LOGICAL :: output_nested

  NAMELIST /device/ &
    & electronmass, &
    & nonlin_as,    &
    & magnetic

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
    & allow_pot_interpolation, &
    & strpotentialX, &
    & strpotentialY, &
    & strpotentialXY, &
    & pot_file_in, &
    & pot_filelist_name

  NAMELIST /write_out/ &
    & run_name, &
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

  CHARACTER(256) :: nml_folder
  INTEGER :: pt, pt2, sl_count
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
    nmlfile_name = nmlfile
    INFO = 0

    output_nested = 0
    write_folder = TRIM(run_name)//"_out"
  ! Now check if the namelist folder is equal to the write folder. In this case,
  ! do not create the write folder
    sl_count = 0
    DO pt = LEN_TRIM(nmlfile_name), 1, -1
      if (nmlfile_name(pt:pt) == '/' .or. nmlfile_name(pt:pt) == '\') then
        sl_count = sl_count + 1
      end if
      if (sl_count == 2) then
        exit
      end if
    end do
    if (pt > 1) then
      pt = pt + 1
      do pt2 = pt, LEN_TRIM(nmlfile_name)
        if (nmlfile_name(pt2:pt2) == '/' .or. nmlfile_name(pt2:pt2) == '\') then
          exit
        end if
      end do
      pt2 = pt2 - 1
      nml_folder = nmlfile_name(pt:pt2)
      if (TRIM(nml_folder) == TRIM(write_folder)) then
        output_nested = 1
      end if
    end if
  
  ! concat the write folder name to the namelist name
    DO pt = LEN_TRIM(nmlfile_name), 1, -1
      if (nmlfile_name(pt:pt) == '/' .or. nmlfile_name(pt:pt) == '\') then
        exit
      end if
    end do
    if (pt > 1) then
      if (output_nested == 0) then
        write_folder = TRIM(nmlfile_name(1:pt))//TRIM(run_name)//"_out"
      else
        write_folder = TRIM(nmlfile_name(1:pt-1))
      end if
    end if

  else
    INFO = 1
  end if

END SUBROUTINE INDATA_GET

SUBROUTINE INDATA_FILL_WITH_DEFAULT
  electronmass = 0.067
  mstar = electronmass * ELMASS0
  magnetic = 0
  nonlin_as = 1e-6
  psi_mode = 'gauss'
  x0 = 2e-7
  y0 = 5e-7
  sigmax = 5e-8
  sigmay = 5e-8
  xenergy = 3.7e-2
  yenergy = 0.0
  psi_file_in = ""
  dt = 5.0e-15
  MAXIT = 300
  grid_mode = 'uniform'
  numx = 512
  numy = 512
  size_x = 1e-6
  size_y = 1e-6
  grid_file_in = ""
  allow_pot_interpolation = 1
  strpotentialX = ""
  strpotentialY = ""
  strpotentialxy = "box 5.0e-7 5.5e-7 5.2e-7 1.0e-6 30.0e-3;box 5.0e-7 0.0e-7 5.2e-7 4.5e-7 30.0e-3;"
  pot_file_in = ""
  pot_filelist_name = ""
  run_name = "default"
  write_folder = "default_out"
  output_nested = 0
  write_grid = "bin"
  write_pot = "bin"
  write_psi = "bin"
  write_timestep = 1e-13
  write_downsample_x = 2
  write_downsample_y = 2
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
SUBROUTINE INDATA_COMPUTE(OKCANCEL, SOLVE_METHOD)
  INTEGER, INTENT(OUT) :: OKCANCEL
  INTEGER, INTENT(IN) :: SOLVE_METHOD
  INTEGER :: nx, ny, grid_numx, grid_numy
  INTEGER :: lenstr, INFO, pt, pt2
  LOGICAL :: file_exists
  
  OKCANCEL = 0
  CALL PXFMKDIR(TRIM(write_folder), LEN_TRIM(write_folder), 8*8*8-1, INFO)

  INQUIRE(FILE=TRIM(write_folder)//'/grid.dat', EXIST=file_exists)
  if (file_exists) then
    call dwgbut("Output folder contains the output of a previous run. Do you want to overwrite it?", INFO)
    if (INFO == 1) then
    else
      OKCANCEL = 1;
      return;
    end if
  end if

  if (.not. output_nested) then
    DO pt = LEN_TRIM(nmlfile_name), 1, -1
      if (nmlfile_name(pt:pt) == '/' .or. nmlfile_name(pt:pt) == '\') then
        exit
      end if
    end do

    call th_copy_file(TRIM(nmlfile_name), TRIM(write_folder)//TRIM(nmlfile_name(pt:LEN_TRIM(nmlfile_name))))

    if (TRIM(pot_file_in) /= "") then
      DO pt2 = LEN_TRIM(pot_file_in), 1, -1
        if (pot_file_in(pt2:pt2) == '/' .or. pot_file_in(pt2:pt2) == '\') then
          exit
        end if
      end do
      pt2 = pt2 + 1
      call th_copy_file(pot_file_in, TRIM(write_folder)//TRIM('/'//pot_file_in(pt2:LEN_TRIM(pot_file_in))))
    end if

    if (TRIM(psi_file_in) /= "") then
      DO pt2 = LEN_TRIM(psi_file_in), 1, -1
        if (psi_file_in(pt2:pt2) == '/' .or. psi_file_in(pt2:pt2) == '\') then
          exit
        end if
      end do
      pt2 = pt2 + 1
      call th_copy_file(psi_file_in, TRIM(write_folder)//TRIM('/'//psi_file_in(pt2:LEN_TRIM(psi_file_in))))
    end if

    if (TRIM(grid_file_in) /= "") then
      DO pt2 = LEN_TRIM(grid_file_in), 1, -1
        if (grid_file_in(pt2:pt2) == '/' .or. grid_file_in(pt2:pt2) == '\') then
          exit
        end if
      end do
      pt2 = pt2 + 1
      call th_copy_file(grid_file_in, TRIM(write_folder)//TRIM('/'//grid_file_in(pt2:LEN_TRIM(grid_file_in))))
    end if

  endif
  
!*******************************************************
!    Preload potential file
!*******************************************************
  if (TRIM(pot_file_in) /= "") then
    call READ_POT_FILE(pot_file_in, INFO)
    if (INFO == 0) then
      print *, "Error reading potential file."
      STOP
    end if
  end if
  
  if (pot_filelist_name /= "") then
    call READ_POT_FILELIST(INFO)
    if (INFO == 1) then
      print *, "Error reading potential file-list."
      STOP
    end if
  end if

!*******************************************************
!    Grid init
!*******************************************************
  if (SOLVE_METHOD == 0) then
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
      ALLOCATE (pot_file_static(numy,numx))
      ALLOCATE (pot_filelist(numy,numx))
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
      ALLOCATE (pot_file_static(numy,numx))
      ALLOCATE (pot_filelist(numy,numx))
      ALLOCATE (psi(numx*numy))
      ALLOCATE (xnodes(numx+2))
      ALLOCATE (ynodes(numy+2))
      xnodes = pot_x_file
      ynodes = pot_y_file

    else  ! uniform
      ALLOCATE (pot(numy,numx))
      ALLOCATE (pot_file_static(numy,numx))
      ALLOCATE (pot_filelist(numy,numx))
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

  else if (SOLVE_METHOD == 1) then
    ! Split Step
    ALLOCATE (pot(numy,numx))
    ALLOCATE (pot_file_static(numy,numx))
    ALLOCATE (pot_filelist(numy,numx))
    ALLOCATE (psi(numx*numy))
    ALLOCATE (xnodes(numx+2))
    ALLOCATE (ynodes(numy+2))
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
  end if

!*******************************************************
!    Potential init
!*******************************************************

  pot_filelist = 0
  if (TRIM(pot_file_in) /= "") then
    if (grid_mode == "pot") then
      DO nx = 1, numx
        DO ny = 1, numy
          pot_file_static(ny, nx)= pot_file_raw(ny+1, nx+1)
        END DO
      END DO
    else
      ! check if resampling is really needed
      if ((pot_numx_file-2) == numx .and. (pot_numy_file-2) == numy) then
        DO nx = 1, numx
          DO ny = 1, numy
            pot_file_static(ny, nx)= pot_file_raw(ny+1, nx+1)
          END DO
        END DO
      else if (allow_pot_interpolation > 0) then
        CALL resample_matrix_2D(pot_file_raw, pot_x_file, pot_y_file, pot_numx_file, pot_numy_file, pot_file_static, xnodes(2:numx), ynodes(2:numy), numx, numy)
      else
        print *, "Potential interpolation needed, but not allowed by user"
        STOP
      end if
    end if

    DEALLOCATE (pot_file_raw)
    DEALLOCATE (pot_x_file)
    DEALLOCATE (pot_y_file)
  else
    pot_file_static = 0.
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
  LOGICAL :: file_exists

  INFO = 0
  INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
  if (.not. file_exists) then
    return
  end if
  OPEN(22, FILE=TRIM(fname), FORM="FORMATTED", STATUS="OLD")
  do
    read (22, '(A)', end=999) line
    if (line(1:5) == "#TDS " .or. line(1:5) == "#tds ") then
      READ(line, *) str1, pot_numx_file, pot_numy_file
      EXIT
    end if
  end do

  ALLOCATE (pot_file_raw(pot_numy_file, pot_numx_file))
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

    pot_file_raw(ny, nx) = val
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
  LOGICAL :: file_exists

  INFO = 0
  INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
  if (.not. file_exists) then
    return
  end if
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

SUBROUTINE READ_POT_FILELIST(INFO)
  INTEGER, INTENT(OUT) :: INFO
  CHARACTER(512) :: line, str1
  LOGICAL :: file_exists
  INTEGER :: i

  INFO = 0
  INQUIRE(FILE=TRIM(pot_filelist_name), EXIST=file_exists)
  if (.not. file_exists) then
    return
  end if

  OPEN(22, FILE=TRIM(pot_filelist_name), FORM="FORMATTED", STATUS="OLD")
  filelist_count = 0
  do
    read (22, '(A)', end=998) line
    filelist_count = filelist_count + 1
  end do

998 CLOSE(22)
  ALLOCATE(filelist_data(filelist_count))

  OPEN(22, FILE=TRIM(pot_filelist_name), FORM="FORMATTED", STATUS="OLD")
  do i = 1, filelist_count
    read(22, *) filelist_data(i)%time, filelist_data(i)%file
  end do
  CLOSE(22)

END SUBROUTINE


SUBROUTINE manage_pot_filelist(time, index)
  INTEGER, INTENT(INOUT) :: index
  REAL*8, INTENT(IN) :: time
  INTEGER :: INFO, nx, ny
  
  if (filelist_count == 0 .or. index == filelist_count) then
    return
  end if
  
  if (time < filelist_data(index+1)%time) then
    return
  end if

  index = index + 1
  call READ_POT_FILE(filelist_data(index)%file, INFO)
  if (INFO == 0) then
    return
  end if

  if (grid_mode == "pot") then
    DO nx = 1, numx
      DO ny = 1, numy
        pot_filelist(ny, nx)= pot_file_raw(ny+1, nx+1)
      END DO
    END DO
  else
    ! check if resampling is really needed
    if ((pot_numx_file-2) == numx .and. (pot_numy_file-2) == numy) then
      DO nx = 1, numx
        DO ny = 1, numy
          pot_filelist(ny, nx)= pot_file_raw(ny+1, nx+1)
        END DO
      END DO
    else if (allow_pot_interpolation > 0) then
      CALL resample_matrix_2D(pot_file_raw, pot_x_file, pot_y_file, pot_numx_file, pot_numy_file, pot_filelist, xnodes(2:numx), ynodes(2:numy), numx, numy)
    else
      print *, "Potential interpolation needed, but not allowed by user"
      STOP
    end if
  end if

  DEALLOCATE (pot_file_raw)
  DEALLOCATE (pot_x_file)
  DEALLOCATE (pot_y_file)
END SUBROUTINE

END MODULE mod_indata
