MODULE mod_viewgraph
  USE dislin
  USE mod_staticdata
  USE mod_indata
  USE mod_filewrite
  USE mod_bicubic
  IMPLICIT NONE
  SAVE

  INTEGER :: hWaveGraph, hCutoff, hViewList, hShowModulus
  INTEGER :: hIntegX0, hIntegY0, hIntegX1, hIntegY1, hIntegrate, hIntegral1, hIntegral2, hSlice, hSliceCoord
  INTEGER :: IsViewAllocated

CONTAINS



SUBROUTINE GetSelectedTimestep(str_num)
  CHARACTER(4), INTENT(OUT) :: str_num
  REAL*8 :: next_write_time
  INTEGER nf, cnt, iter

  str_num = "";
  CALL GWGLIS(hViewList, nf)
  if (nf <= 0 .or. nf > 9999) then
    return;
  end if
  cnt = 1;
  next_write_time = 0.
  DO iter = 1, MAXIT
    if (((iter-1)*dt) >= next_write_time) then
      if (cnt == nf) then
        exit
      end if
      next_write_time = next_write_time + write_timestep
      cnt = cnt + 1
    end if
  END DO
  str_num = STRING(iter-1)
END SUBROUTINE





SUBROUTINE ViewListCbk(ID)
  INTEGER, INTENT(IN) :: ID
  INTEGER :: read_NUMX, read_NUMY, read_NUMZ, numx, numy, idx
  REAL :: psimax, vmax, cutoff, slice_coord
  REAL*8, ALLOCATABLE :: read_xnodes(:), read_ynodes(:), read_znodes(:), new_xnodes(:), new_ynodes(:)

  REAL*8, ALLOCATABLE :: vstr_in(:,:,:)
  COMPLEX*16, ALLOCATABLE :: psic(:,:,:)
  REAL, ALLOCATABLE :: vstr_t(:,:), psi_t(:,:)

  REAL*8 :: intx0, inty0, intx1, inty1, int_box, int_all
  INTEGER :: nx, ny, nz, iter, show_mod, slice
  LOGICAL :: file_exists
  CHARACTER(4) :: str_num
  CHARACTER(80) :: valstr

  CALL GetSelectedTimestep(str_num)

  CALL GWGTXT(hCutoff, valstr)
  read (valstr, *) cutoff
  if (cutoff > 100) then
    cutoff = 100
  end if
  cutoff = cutoff / 100

  INQUIRE(FILE=TRIM(write_folder)//'/grid.dat', EXIST=file_exists)
  if (file_exists) then
    OPEN(33, FILE=TRIM(write_folder)//'/grid.dat', FORM="FORMATTED", ACTION="READ",  &
         &   STATUS="OLD")
    READ(33,*) read_NUMX, read_NUMY, read_NUMZ
    ALLOCATE(read_xnodes(read_NUMX))
    ALLOCATE(read_ynodes(read_NUMY))
    ALLOCATE(read_znodes(read_NUMZ))
    read (33, *) read_xnodes(1:read_NUMX)
    read (33, *) read_ynodes(1:read_NUMY)
    read (33, *) read_znodes(1:read_NUMZ)
    CLOSE(33)
  else
    return;
  end if

  ALLOCATE(vstr_in(read_NUMZ, read_NUMY, read_NUMX))
  ALLOCATE(psic(read_NUMZ, read_NUMY, read_NUMX))

  INQUIRE(FILE=TRIM(write_folder)//'/potential'//str_num//'.bin', EXIST=file_exists)
  if (file_exists) then
    OPEN(33, FILE=TRIM(write_folder)//'/potential'//str_num//'.bin', FORM="UNFORMATTED", ACTION="READ",  &
         &   STATUS="OLD")
    read (33) vstr_in
    CLOSE(33)
  else
    return
  end if

  INQUIRE(FILE=TRIM(write_folder)//'/psi'//str_num//'.bin', EXIST=file_exists)
  if (file_exists) then
    OPEN(33, FILE=TRIM(write_folder)//'/psi'//str_num//'.bin', FORM="UNFORMATTED", ACTION="READ",  &
         &   STATUS="OLD")
    read (33) psic
    CLOSE(33)
  else
    return
  end if

  CALL GWGBOX(hShowModulus, show_mod)
  CALL GWGLIS(hSlice, slice)
  CALL GWGTXT(hSliceCoord, valstr)
  read (valstr, *) slice_coord

  if (slice == 1) then
    ! section along X axis
    numx = read_numy
    numy = read_numz
    ALLOCATE(new_xnodes(numx))
    ALLOCATE(new_ynodes(numy))
    new_xnodes = read_ynodes;
    new_ynodes = read_znodes;

    idx = 1
    DO WHILE ((idx < read_numx) .and. (read_xnodes(idx) < slice_coord))
      idx = idx + 1
    END DO

    ALLOCATE(psi_t(NUMX, NUMY))
    ALLOCATE(vstr_t(NUMX, NUMY))
    DO NX = 1, NUMX
      DO NY = 1, NUMY
        if (show_mod == 2) then
          psi_t(NX, NY) = REAL(psic(NY, NX, idx))
        else if (show_mod == 3) then
          psi_t(NX, NY) = IMAG(psic(NY, NX, idx))
        else
          psi_t(NX, NY) = REAL(abs(psic(NY, NX, idx))**2)
        end if
        vstr_t(NX, NY) = REAL(vstr_in(NY, NX, idx)) / ELCH
      END DO
    END DO

  else if (slice == 2) then
    ! section along Y axis
    numx = read_numx
    numy = read_numz
    ALLOCATE(new_xnodes(numx))
    ALLOCATE(new_ynodes(numy))
    new_xnodes = read_xnodes;
    new_ynodes = read_znodes;

    idx = 1
    DO WHILE ((idx < read_numy) .and. (read_ynodes(idx) < slice_coord))
      idx = idx + 1
    END DO

    ALLOCATE(psi_t(NUMX, NUMY))
    ALLOCATE(vstr_t(NUMX, NUMY))
    DO NX = 1, NUMX
      DO NY = 1, NUMY
        if (show_mod == 2) then
          psi_t(NX, NY) = REAL(psic(NY, idx, NX))
        else if (show_mod == 3) then
          psi_t(NX, NY) = IMAG(psic(NY, idx, NX))
        else
          psi_t(NX, NY) = REAL(abs(psic(NY, idx, NX))**2)
        end if
        vstr_t(NX, NY) = REAL(vstr_in(NY, idx, NX)) / ELCH
      END DO
    END DO

  else if (slice == 3) then
    ! section along Z axis
    numx = read_numx
    numy = read_numy
    ALLOCATE(new_xnodes(numx))
    ALLOCATE(new_ynodes(numy))
    new_xnodes = read_xnodes;
    new_ynodes = read_ynodes;

    idx = 1
    DO WHILE ((idx < read_numz) .and. (read_znodes(idx) < slice_coord))
      idx = idx + 1
    END DO

    ALLOCATE(psi_t(NUMX, NUMY))
    ALLOCATE(vstr_t(NUMX, NUMY))
    DO NX = 1, NUMX
      DO NY = 1, NUMY
        if (show_mod == 2) then
          psi_t(NX, NY) = REAL(psic(idx, NY, NX))
        else if (show_mod == 3) then
          psi_t(NX, NY) = IMAG(psic(idx, NY, NX))
        else
          psi_t(NX, NY) = REAL(abs(psic(idx, NY, NX))**2)
        end if
        vstr_t(NX, NY) = REAL(vstr_in(idx, NY, NX)) / ELCH
      END DO
    END DO

  endif

  DEALLOCATE(vstr_in)
  DEALLOCATE(psic)
  DEALLOCATE(read_xnodes)
  DEALLOCATE(read_ynodes)
  DEALLOCATE(read_znodes)

  psimax= MAXVAL(psi_t(:,:))
  vmax= MAXVAL(vstr_t(:,:))

!  CALL METAFL('PNG')
!  CALL SETPAG('DA4P')
!  CALL FILMOD('DELETE')
!  CALL SETFIL(TRIM(write_folder)//'/img'//str_num//'.png')

  CALL METAFL('CONS')

  CALL SETXID(hWaveGraph, 'WIDGET')
  CALL page(3500, 2600)
  CALL PAGMOD('LAND')
  CALL X11MOD('AUTO')
  CALL DISINI()
  CALL errmod('ALL', 'OFF')    
  CALL ERASE()
  CALL PAGERA()
  CALL complx

!...............................................crea gli assi X e Y (e Z per il pot)
  CALL WINMOD('NOERASE')
  CALL AUTRES(NUMX, NUMY)
  CALL AXSPOS(280, 2400)
  CALL AX3LEN(2200, 2200, 1800)
  CALL DIGITS(-1, 'XY')
  CALL DIGITS(2, 'Z')

  if (slice == 1) then
    CALL NAME('y (nm)', 'X')
    CALL NAME('z (nm)', 'Y')
  else if (slice == 2) then
    CALL NAME('x (nm)', 'X')
    CALL NAME('z (nm)', 'Y')
  else if (slice == 3) then
    CALL NAME('x (nm)', 'X')
    CALL NAME('y (nm)', 'Y')
  end if

  CALL LABELS ('FEXP', 'XYZ')
  CALL NAME('potential (eV)', 'Z')
  CALL SETVLT('GREY')
  CALL COLRAN(20, 120)
  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
           & 0., vmax, 0., vmax / 5.)
  CALL ENDGRF

!..........................................................crea l'asse Z per la psi
  CALL NOGRAF()
  CALL MYVLT( (/ (nx/250., nx= 1, 250) /), (/ (nx*1e-3, nx= 1, 250) /),      &
       &      (/ (nx*1e-3, nx= 1, 250) /), 250)
  CALL COLRAN(50, 250)
  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
           & psimax*cutoff, psimax, psimax*cutoff, psimax / 50.)
  CALL SMXALF('GREEK', '@', '@', 1)
  CALL MIXALF
  CALL LABELS ('FEXP', 'Z')
  CALL ZAXIS(0., psimax, 0., psimax / 5., 1800,'|@u@|[2$ probability density',0,0,3020,2400 )
  CALL ENDGRF

  !..................................grafica il potenziale
  CALL SETVLT('GREY')
  CALL COLRAN(20, 120)
!  CALL SETIND(255, 0.4, 0.3, 0.6)  ! crea il colore del pot infinito
!  CALL SETCLR(254)                 ! colore da usare per le scritte
  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
           & 0., vmax, 0., vmax / 5.)
  CALL CRVMAT(vstr_t, NUMX, NUMY, 1, 1)
  CALL ENDGRF

  !.............................................grafica le f. d'onda
  CALL MYVLT( (/ (nx/250., nx= 1, 250) /), (/ (nx*1e-3, nx= 1, 250) /),      &
       &      (/ (nx*1e-3, nx= 1, 250) /), 250)
  CALL COLRAN(50, 240)
  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
           & psimax*cutoff, psimax, psimax*0.02, psimax / 50.)

  CALL NOBGD
  CALL CRVMAT(psi_t, NUMX, NUMY, 2, 2)
  CALL RESET('NOBGD')
  CALL ENDGRF

  CALL SETVLT("RAIN")
  CALL COLOR("GREEN")
  CALL RECTAN (INT(280. + 2200.*intx0/size_x), INT(200. + 2200.*inty0/size_y),  &
    INT(2200.*(intx1-intx0)/size_x), INT(2200.*(inty1-inty0)/size_y))

  CALL DISFIN()

  DEALLOCATE(vstr_t)
  DEALLOCATE(psi_t)
  DEALLOCATE(new_xnodes)
  DEALLOCATE(new_ynodes)
END SUBROUTINE



SUBROUTINE ViewListCbk2(ID)
  INTEGER, INTENT(IN) :: ID
  INTEGER :: N
  REAL, DIMENSION (100, 100) :: ZMAT
  REAL    :: FPI,STEP,X,Y
  INTEGER :: I,J

  N = 100
  FPI=3.1415927/180.
      STEP=360./(N-1)
      DO I=1,N
        X=(I-1.)*STEP
        DO J=1,N
          Y=(J-1.)*STEP
          ZMAT(I,J)=2*SIN(X*FPI)*SIN(Y*FPI)
        END DO
      END DO
  CALL METAFL('CONS')
  CALL SETXID(hWaveGraph, 'WIDGET')
  CALL page(2500, 2500)
  CALL metafl('CONS')
  CALL DISINI()
  CALL errmod('ALL', 'OFF')    
  CALL ERASE()
  CALL PAGERA()
  CALL complx
!  CALL HWFONT()

!      CALL TITLIN('3-D Colour Plot of the Function',2)
!      CALL TITLIN('F(X,Y) = 2 * SIN(X) * SIN(Y)',4)

!      CALL NAME('X-axis','X')
!      CALL NAME('Y-axis','Y')
!      CALL NAME('Z-axis','Z')

      CALL INTAX()
      CALL AUTRES(N,N)
      CALL AXSPOS(300,1850)
      CALL AX3LEN(2200,1400,1400)

      CALL GRAF3(0.,360.,0.,90.,0.,360.,0.,90.,-2.,2.,-2.,1.)
      CALL CRVMAT(ZMAT,N,N,1,1)

      CALL HEIGHT(50)
      CALL TITLE()
      CALL MPAEPL(3)
      CALL DISFIN()

END SUBROUTINE

!
!END MODULE
!
!
!MODULE mod_viewgraph
!  USE mod_staticdata
!  USE mod_indata
!  USE mod_viewgraphcbk
!  IMPLICIT NONE
!  SAVE
!
!CONTAINS

SUBROUTINE SHOW_GRAPH

  INTEGER :: hViewForm, hNull, hOk, ctrlh
  REAL*8 :: next_write_time
  INTEGER :: iter
  CHARACTER(40) :: vstr
  CHARACTER(2000) :: list_str

  list_str = ""
  next_write_time = 0.
  DO iter = 1, MAXIT
    if (((iter-1)*dt) >= next_write_time) then
      WRITE (vstr, '(ES10.4)') (iter-1)*dt
      if (iter == 1) then
        list_str = TRIM(vstr)
      else
        list_str = TRIM(list_str)//'|'//TRIM(vstr)
      end if
      next_write_time = next_write_time + write_timestep
    end if
  END DO

    ! Main navigation
  CALL SWGTIT('Display Wave Function')
  CALL WGINI('FORM', hViewForm)
  if (OSV == 0) then
    CALL SWGFNT ("fixed", 10)  ! lunux
    ctrlh = 30
  else
    CALL SWGFNT ("fixed", 14)  ! windows
    ctrlh = 22
  end if

  CALL SWGWIN(5, 5, 100, 35)
  CALL WGDLIS(hViewForm, list_str, 1, hViewList)
  CALL SWGCBK(hViewList, ViewListCbk)
  CALL SWGWIN(110, 5, 50, 35)
  CALL WGLAB(hViewForm, 'Cut off:', hNull)
  CALL SWGWIN(160, 5, 30, ctrlh)
  CALL WGTXT(hViewForm, '10', hCutoff)
  CALL SWGWIN(190, 5, 30, 35)
  CALL WGLAB(hViewForm, '%', hNull)
  CALL SWGWIN(230, 8, 50, 28)
  CALL WGOK(hViewForm, hOk)

  if (OSV == 0) then
    CALL SWGWIN(310, 2, 100, 125)
  else
    CALL SWGWIN(310, 2, 100, 25)
  end if
  CALL WGBOX(hViewForm, 'Modulus|Real|Imag', 1, hShowModulus)
  CALL SWGCBK(hShowModulus, ViewListCbk)

    ! Integral bar
!  CALL SWGWIN(5, 45, 25, 30)
!  CALL WGLAB(hViewForm, 'X0:', hNull)
!  CALL SWGWIN(30, 45, 70, ctrlh)
!  CALL WGTXT(hViewForm, '0', hIntegX0)
!  CALL SWGWIN(80, 45, 25, 30)
!  CALL WGLAB(hViewForm, 'Y0:', hNull)
!  CALL SWGWIN(105, 45, 70, ctrlh)
!  CALL WGTXT(hViewForm, '0', hIntegY0)
!  CALL SWGWIN(180, 45, 25, 30)
!  CALL WGLAB(hViewForm, 'X1:', hNull)
!  CALL SWGWIN(205, 45, 70, ctrlh)
!  write (vstr, '(ES10.4)') size_x
!  CALL WGTXT(hViewForm, TRIM(vstr), hIntegX1)
!  CALL SWGWIN(280, 45, 25, 30)
!  CALL WGLAB(hViewForm, 'X1:', hNull)
!  CALL SWGWIN(305, 45, 70, ctrlh)
!  write (vstr, '(ES10.4)') size_y
!  CALL WGTXT(hViewForm, TRIM(vstr), hIntegY1)

!  CALL SWGWIN(380, 46, 70, 28)
!  CALL WGPBUT(hViewForm, "Integrate", hIntegrate)
!  CALL SWGCBK(hIntegrate, ViewListCbk)
!  CALL SWGWIN(460, 45, 100, 30)
!  CALL WGLAB(hViewForm, '', hIntegral1)
!  CALL SWGWIN(570, 45, 100, 30)
!  CALL WGLAB(hViewForm, '', hIntegral2)

  CALL SWGWIN(5, 45, 35, 30)
  CALL WGLAB(hViewForm, 'Slice:', hNull)
  CALL SWGWIN(45, 45, 60, 35)
  CALL WGDLIS(hViewForm, "X|Y|Z", 3, hSlice)
  CALL SWGCBK(hSlice, ViewListCbk)
  write (vstr, '(ES11.4)') size_z/2.
  CALL SWGWIN(120, 45, 100, 30)
  CALL WGTXT(hViewForm, TRIM(vstr), hSliceCoord)

  CALL SWGWIN(5, 100, 800, 600)
  CALL WGDRAW(hViewForm, hWaveGraph)
  CALL WGFIN

END SUBROUTINE

END MODULE

