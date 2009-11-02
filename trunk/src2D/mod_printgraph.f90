MODULE mod_printgraph
  USE dislin
  USE mod_staticdata
  USE mod_indata
  USE mod_filewrite
  USE mod_bicubic
  IMPLICIT NONE
  SAVE

CONTAINS


SUBROUTINE PRINT_GRAPH

  INTEGER :: hWaveGraph
  INTEGER :: hViewForm, hNull, hOk
  REAL*8 :: next_write_time
  CHARACTER(40) :: valstr
  REAL :: psimax_l, vmax_l, psimax, vmax, l_cutoff, u_cutoff
  REAL*8, ALLOCATABLE :: xnodes(:), ynodes(:), new_xnodes(:), new_ynodes(:)

  REAL*8, ALLOCATABLE :: vstr_in(:,:), vstr(:,:), psi_in(:,:), psi(:,:)
  COMPLEX*16, ALLOCATABLE :: psic(:,:)
  REAL, ALLOCATABLE :: vstr_t(:,:), psi_t(:,:)

  INTEGER :: nx, ny, is_uniform, iter, show_mod
  LOGICAL :: file_exists
  REAL :: red(240), green(240), blue(240)
  
  do nx = 1, 120
    red(nx) = 0.7*REAL(nx)/120
    green(nx) = 0.7*REAL(nx)/120
    blue(nx) = 0.7*REAL(nx)/120
  end do
  do nx = 121, 240
    red(nx) = REAL(nx-120)/120
    green(nx) = REAL(nx-120)*3e-3
    blue(nx) = REAL(nx-120)*3e-3
  end do

  l_cutoff = 0.02
  u_cutoff = 0.25
  
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
  else
    return;
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
    if (abs((ynodes(NY+1)-ynodes(NY)) - (ynodes(NY+2)-ynodes(NY+1))) > 1e-5*ynodes(NY+1)) then
      is_uniform = 0
      EXIT
    end if
  end do

  ALLOCATE(vstr_in(NUMY, NUMX))
  ALLOCATE(psi_in(NUMY, NUMX))
  ALLOCATE(psic(NUMY, NUMX))



! Find the maximum of the potential and of the wave function along the frames
  psimax = 0.
  vmax = 0.
  next_write_time = 0.
  DO iter = 1, MAXIT
    if (((iter-1)*dt) >= next_write_time) then
      next_write_time = next_write_time + write_timestep

      INQUIRE(FILE=TRIM(write_folder)//'/potential'//STRING(iter-1)//'.bin', EXIST=file_exists)
      if (file_exists) then
        OPEN(33, FILE=TRIM(write_folder)//'/potential'//STRING(iter-1)//'.bin', FORM="UNFORMATTED", ACTION="READ",  &
             &   STATUS="OLD")
        read (33) vstr_in
        CLOSE(33)
      else
        return
      end if

      INQUIRE(FILE=TRIM(write_folder)//'/psi'//STRING(iter-1)//'.bin', EXIST=file_exists)
      if (file_exists) then
        OPEN(33, FILE=TRIM(write_folder)//'/psi'//STRING(iter-1)//'.bin', FORM="UNFORMATTED", ACTION="READ",  &
             &   STATUS="OLD")
        read (33) psic
        CLOSE(33)
      else
        return
      end if
      psimax_l= REAL(MAXVAL(ABS(psic(:,:)) ** 2))
      vmax_l= REAL(MAXVAL(vstr_in(:,:))) / ELCH
      if (psimax_l > psimax) psimax = psimax_l
      if (vmax_l > vmax) vmax = vmax_l
    end if
  end do



! Draw loop
  next_write_time = 0.
  DO iter = 1, MAXIT
    if (((iter-1)*dt) >= next_write_time) then
      WRITE (valstr, '(ES10.4)') (iter-1)*dt
      next_write_time = next_write_time + write_timestep

  vstr_in = 0.
  psi_in = 0.
  psic = 0.

  INQUIRE(FILE=TRIM(write_folder)//'/potential'//STRING(iter-1)//'.bin', EXIST=file_exists)
  if (file_exists) then
    OPEN(33, FILE=TRIM(write_folder)//'/potential'//STRING(iter-1)//'.bin', FORM="UNFORMATTED", ACTION="READ",  &
         &   STATUS="OLD")
    read (33) vstr_in
    CLOSE(33)
  else
    return
  end if

  INQUIRE(FILE=TRIM(write_folder)//'/psi'//STRING(iter-1)//'.bin', EXIST=file_exists)
  if (file_exists) then
    OPEN(33, FILE=TRIM(write_folder)//'/psi'//STRING(iter-1)//'.bin', FORM="UNFORMATTED", ACTION="READ",  &
         &   STATUS="OLD")
    read (33) psic
    CLOSE(33)
  else
    return
  end if

  DO NX = 1, NUMX
    DO NY = 1, NUMY
      psi_in(NY, NX) = abs(psic(NY, NX))**2
    END DO
  END DO

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

  do nx = 1, numx
    do ny = 1, numy
      if (psi_t(nx, ny) < l_cutoff*psimax) psi_t(nx, ny) = 0.
      if (psi_t(nx, ny) > u_cutoff*psimax) psi_t(nx, ny) = u_cutoff*psimax
    end do
  end do
  
  CALL METAFL('PNG')
  CALL FILMOD('DELETE')

!  CALL SETXID(hWaveGraph, 'WIDGET')
  CALL page(3500, 2600)
  CALL SETFIL(TRIM(write_folder)//'/psi'//STRING(iter-1)//'.png')
  CALL PAGMOD('LAND')
  CALL X11MOD('AUTO')
  CALL DISINI()
  CALL errmod('ALL', 'OFF')    
  CALL ERASE()
  CALL PAGERA()
  CALL complx

!...............................................crea gli assi X e Y (e Z per il pot)
!  CALL WINMOD('NOERASE')
  CALL AUTRES(NUMX, NUMY)
  CALL AXSPOS(280, 2400)
  CALL AX3LEN(2200, 2200, 1800)
  CALL DIGITS(-1, 'XY')
  CALL DIGITS(2, 'Z')
  CALL NAME('x (nm)', 'X')
  CALL NAME('y (nm)', 'Y')

  CALL LABELS ('FEXP', 'XYZ')
  CALL NAME('potential (eV)', 'Z')
  CALL MYVLT(red, green, blue, 240)
  CALL COLRAN(1, 119)
  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
           & 0., vmax, 0., vmax / 5.)
  CALL CRVMAT(vstr_t, NUMX, NUMY, 2, 2)
  CALL ENDGRF

!..........................................................crea l'asse Z per la psi
  CALL NOGRAF()
  CALL COLRAN(121, 239)
  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
           & psimax*l_cutoff, u_cutoff*psimax, psimax*l_cutoff, u_cutoff*psimax / 5.)
  CALL SMXALF('GREEK', '@', '@', 1)
  CALL MIXALF
  CALL LABELS ('FEXP', 'Z')
  CALL ZAXIS(0., u_cutoff*psimax, 0., u_cutoff*psimax / 5., 1800,'|@u@|[2$ probability density',0,0,3020,2400 )
  CALL NOBGD
  CALL CRVMAT(psi_t, NUMX, NUMY, 2, 2)
  CALL RESET('NOBGD')
  CALL ENDGRF

  !..................................grafica il potenziale
!  CALL MYVLT(red, green, blue, 240)
!  CALL COLRAN(1, 120)
!  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
!           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
!           & 0., vmax, 0., vmax / 5.)
!  CALL CRVMAT(vstr_t, NUMX, NUMY, 1, 1)
!  CALL ENDGRF

  !.............................................grafica le f. d'onda
!  CALL MYVLT(red, green, blue, 240)
!  CALL COLRAN(121, 240)
!  CALL GRAF3(REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9, REAL(new_xnodes(1))*1e9, REAL(new_xnodes(NUMX))*1e9/5., &
!           & REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9, REAL(new_ynodes(1))*1e9, REAL(new_ynodes(NUMY))*1e9/5., &
!           & psimax*cutoff, psimax, psimax*0.02, psimax / 5.)

!  CALL NOBGD
!  CALL CRVMAT(psi_t, NUMX, NUMY, 2, 2)
!  CALL RESET('NOBGD')
!  CALL ENDGRF

  CALL DISFIN()

  DEALLOCATE(vstr_t)
  DEALLOCATE(psi_t)
  DEALLOCATE(new_xnodes)
  DEALLOCATE(new_ynodes)

    end if
  END DO

  DEALLOCATE(vstr_in)
  DEALLOCATE(psi_in)
  DEALLOCATE(xnodes)
  DEALLOCATE(ynodes)
  DEALLOCATE(psic)

END SUBROUTINE

END MODULE

