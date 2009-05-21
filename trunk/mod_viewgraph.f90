MODULE mod_viewgraph
  USE dislin
  USE mod_staticdata
  USE mod_indata
  USE mod_filewrite
  USE mod_bicubic
  IMPLICIT NONE
  SAVE

  INTEGER :: hWaveGraph, hCutoff, hViewList, hShowModulus
  INTEGER :: hIntegX0, hIntegY0, hIntegX1, hIntegY1, hIntegrate, hIntegral1, hIntegral2
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
  INTEGER :: NUMX, NUMY
  REAL :: psimax, vmax, cutoff
  REAL*8, ALLOCATABLE :: xnodes(:), ynodes(:), new_xnodes(:), new_ynodes(:)

  REAL*8, ALLOCATABLE :: vstr_in(:,:), vstr(:,:), psi_in(:,:), psi(:,:)
  COMPLEX*16, ALLOCATABLE :: psic(:,:)
  REAL, ALLOCATABLE :: vstr_t(:,:), psi_t(:,:)

  REAL*8 :: intx0, inty0, intx1, inty1, int_box, int_all
  INTEGER :: nx, ny, is_uniform, iter, show_mod
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
    READ(33,*) NUMX, NUMY
    ALLOCATE(xnodes(NUMX+2))
    ALLOCATE(ynodes(NUMY+2))
    read (33, *) xnodes(2:NUMX+1)
    read (33, *) ynodes(2:NUMY+1)
    xnodes(1) = 2.*xnodes(2) - xnodes(3)
    xnodes(NUMX+2) = 2.*xnodes(NUMX+1) - xnodes(NUMX)
    ynodes(1) = 2.*ynodes(2) - ynodes(3)
    ynodes(NUMY+2) = 2.*ynodes(NUMY+1) - ynodes(NUMY)
    CLOSE(33)
  end if

    ! Check if grid is uniform (if not interpolation will be required)
  is_uniform = 1
  DO NX = 2, NUMX
    if (abs((xnodes(NX+1)-xnodes(NX)) - (xnodes(NX+2)-xnodes(NX+1))) > 1e-5*xnodes(NX+1)) then
      is_uniform = 0
      EXIT
    end if
  end do
  DO NY = 2, NUMY
    if (abs((xnodes(NY+1)-xnodes(NY)) - (xnodes(NY+2)-xnodes(NY+1))) > 1e-5*xnodes(NY+1)) then
      is_uniform = 0
      EXIT
    end if
  end do

  ALLOCATE(vstr_in(NUMY, NUMX))
  ALLOCATE(psi_in(NUMY, NUMX))
  ALLOCATE(psic(NUMY, NUMX))

  INQUIRE(FILE=TRIM(write_folder)//'/potential0000.bin', EXIST=file_exists)
  if (file_exists) then
    OPEN(33, FILE=TRIM(write_folder)//'/potential0000.bin', FORM="UNFORMATTED", ACTION="READ",  &
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

  CALL GWGBUT(hShowModulus, show_mod)
  if (show_mod == 0) then
    DO NX = 1, NUMX
      DO NY = 1, NUMY
        psi_in(NY, NX) = REAL(psic(NY, NX))
      END DO
    END DO
  else
    DO NX = 1, NUMX
      DO NY = 1, NUMY
        psi_in(NY, NX) = abs(psic(NY, NX))**2
      END DO
    END DO
  end if

    ! compute integrals
  CALL GWGTXT(hIntegX0, valstr)
  read (valstr, *) intx0
  intx0 = min(size_x, max(0., intx0))
  CALL GWGTXT(hIntegY0, valstr)
  read (valstr, *) inty0
  inty0 = min(size_y, max(0., inty0))
  CALL GWGTXT(hIntegX1, valstr)
  read (valstr, *) intx1
  intx1 = min(size_x, max(0., intx1))
  CALL GWGTXT(hIntegY1, valstr)
  read (valstr, *) inty1
  inty1 = min(size_y, max(0., inty1))

  CALL Integrate_2D(int_box, psi_in, NUMX, NUMY, xnodes(2:NUMX+1), ynodes(2:NUMY+1), intx0, inty0, intx1, inty1)
  CALL Integrate_2D(int_all, psi_in, NUMX, NUMY, xnodes(2:NUMX+1), ynodes(2:NUMY+1), xnodes(1), ynodes(1), xnodes(NUMX+2), ynodes(NUMY+2))
  if (abs(int_all) > 1.) then
    CALL Integrate_2D(int_box, psi_in, NUMX, NUMY, xnodes(2:NUMX+1), ynodes(2:NUMY+1), intx0, inty0, intx1, inty1)
  endif

  write (valstr, '(ES11.4)') int_box
  CALL SWGTXT(hIntegral1, TRIM(valstr))
  write (valstr, '(ES11.4)') int_all
  CALL SWGTXT(hIntegral2, TRIM(valstr))

    ! Draws graph
  if (is_uniform == 0) then
!    print *, 'Non uniform grid: performing interpolation'
    ALLOCATE(vstr(400, 400))
    ALLOCATE(psi(400, 400))
    ALLOCATE(new_xnodes(400))
    ALLOCATE(new_ynodes(400))
    DO NX = 1,400
      new_xnodes(NX) = xnodes(2) + REAL(NX-1) * (xnodes(numx+1)-xnodes(2)) / 399.
      new_ynodes(NX) = ynodes(2) + REAL(NX-1) * (ynodes(numy+1)-ynodes(2)) / 399.
    END DO
    call resample_matrix_2D(vstr_in, xnodes, ynodes, numx, numy, vstr, new_xnodes, new_ynodes, 400, 400)
    call resample_matrix_2D(psi_in, xnodes, ynodes, numx, numy, psi, new_xnodes, new_ynodes, 400, 400)
    numx = 400
    numy = 400
  else
!    print *, 'Uniform grid: no interpolation needed'
    ALLOCATE(vstr(NUMY, NUMX))
    ALLOCATE(psi(NUMY, NUMX))
    ALLOCATE(new_xnodes(NUMX))
    ALLOCATE(new_ynodes(NUMY))
    vstr = vstr_in
    psi = psi_in
    new_xnodes = xnodes(2:NUMX+1)
    new_ynodes = ynodes(2:NUMY+1)
  end if
  DEALLOCATE(vstr_in)
  DEALLOCATE(psi_in)
  DEALLOCATE(xnodes)
  DEALLOCATE(ynodes)

  ALLOCATE(vstr_t(numx,numy))
  ALLOCATE(psi_t(numx, numy))
  do nx = 1, numx
    do ny = 1, numy
      vstr_t(nx, ny) = REAL(vstr(ny, nx)) / ELCH
      psi_t(nx, ny) = REAL(psi(ny, nx))
    end do
  end do
  DEALLOCATE(vstr)
  DEALLOCATE(psi)

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
  CALL NAME('x (nm)', 'X')
  CALL NAME('y (nm)', 'Y')

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
  CALL SWGWIN(290, 5, 100, 28)
  CALL WGBUT(hViewForm, "Show Modulus", 1, hShowModulus)
  CALL SWGCBK(hShowModulus, ViewListCbk)

    ! Integral bar
  CALL SWGWIN(5, 45, 25, 30)
  CALL WGLAB(hViewForm, 'X0:', hNull)
  CALL SWGWIN(30, 45, 70, ctrlh)
  CALL WGTXT(hViewForm, '0', hIntegX0)
  CALL SWGWIN(80, 45, 25, 30)
  CALL WGLAB(hViewForm, 'Y0:', hNull)
  CALL SWGWIN(105, 45, 70, ctrlh)
  CALL WGTXT(hViewForm, '0', hIntegY0)
  CALL SWGWIN(180, 45, 25, 30)
  CALL WGLAB(hViewForm, 'X1:', hNull)
  CALL SWGWIN(205, 45, 70, ctrlh)
  write (vstr, '(ES10.4)') size_x
  CALL WGTXT(hViewForm, TRIM(vstr), hIntegX1)
  CALL SWGWIN(280, 45, 25, 30)
  CALL WGLAB(hViewForm, 'X1:', hNull)
  CALL SWGWIN(305, 45, 70, ctrlh)
  write (vstr, '(ES10.4)') size_y
  CALL WGTXT(hViewForm, TRIM(vstr), hIntegY1)

  CALL SWGWIN(380, 46, 70, 28)
  CALL WGPBUT(hViewForm, "Integrate", hIntegrate)
  CALL SWGCBK(hIntegrate, ViewListCbk)
  CALL SWGWIN(460, 45, 100, 30)
  CALL WGLAB(hViewForm, '', hIntegral1)
  CALL SWGWIN(570, 45, 100, 30)
  CALL WGLAB(hViewForm, '', hIntegral2)

  CALL SWGWIN(5, 80, 800, 600)
  CALL WGDRAW(hViewForm, hWaveGraph)
  CALL WGFIN

!  CALL SWGTIT('Display Wave Function')
!  CALL WGINI('VERT', hViewForm)
!  CALL WGDLIS(hViewForm, list_str, 1, hViewList)
!  CALL SWGCBK(hViewList, ViewListCbk)
!  CALL WGOK(hViewForm, hOk)
!  CALL WGDRAW(hViewForm, hWaveGraph)
!  CALL WGFIN

END SUBROUTINE

END MODULE

