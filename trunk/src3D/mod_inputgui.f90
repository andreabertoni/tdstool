MODULE mod_inputgui
  USE mod_staticdata
  USE dislin
  USE mod_indata
  USE mod_strtopot
  USE mod_mainalgo
  USE mod_viewgraph
  IMPLICIT NONE
  SAVE

  INTEGER :: hMainForm, hMainCmd, hMainOk

  INTEGER :: hForm, hOk, hNull, hQuit
  INTEGER :: hPsi, hPsiLabel, hPsiMode, hPsiX0, hPsiY0, hPsiZ0, hPsiSigX, hPsiSigY, hPsiSigZ, hPsiXenergy, hPsiYenergy, hPsiZenergy, hPsiFile
  INTEGER :: hTime, hTimeDelta, hTimeSteps
  INTEGER :: hGrid, hGridNumX, hGridNumY, hGridNumZ, hGridSizeX, hGridSizeY, hGridSizeZ
  INTEGER :: hPot, hPotMode, hPotFile, hPotFileList, hPotStrX, hPotStrY, hPotStrZ, hPotStrXYZ, hPotStrXEdit, hPotStrYEdit, hPotStrZEdit, hPotStrXYZEdit
  INTEGER :: hOut, hOutDir, hOutWGridTxt, hOutWGridBin, hOutWPotTxt, hOutWPotBin, hOutWPsiTxt, hOutWPsiBin, hOutTime, hOutDownX, hOutDownY, hOutDownZ
  INTEGER :: hData, hData2, hElMass, hNonlin

  CHARACTER(260) nml_file_name

CONTAINS

SUBROUTINE MAIN_GUI

  INTEGER :: MainMenuCmd, INFO

  CALL DWGFIL('Select the data file to be opened/created', nml_file_name, '*.nml')
  if (nml_file_name == "") then
    return
  end if

  CALL INDATA_GET(nml_file_name, INFO)
  if (INFO /= 0) then
    call INDATA_FILL_WITH_DEFAULT
  end if

  MainMenuCmd = 1
  DO
    CALL SWGTIT('TDS Tool')
    CALL WGINI('FORM', hMainForm)
    if (OSV == 0) then
      CALL SWGFNT ("fixed", 10)  ! linux
    else
      CALL SWGFNT ("fixed", 14)  ! windows
    end if

    CALL SWGWIN(0, 0, 400, 28)
    CALL WGLAB(hMainForm, nml_file_name, hNull)

    if (OSV == 0) then
      CALL SWGWIN(5, 35, 70, 125)
    else
      CALL SWGWIN(5, 35, 70, 25)
    end if
    CALL WGBOX(hMainForm, 'Edit|Run|View|Quit', MainMenuCmd, hMainCmd)
    CALL SWGWIN(100, 70, 50, 28)
    CALL WGOK(hMainForm, hMainOk)

    CALL WGFIN

    CALL GWGBOX(hMainCmd, MainMenuCmd)

    if (MainMenuCmd == 1) then
      CALL SHOW_IN_GUI
    else if (MainMenuCmd == 2) then
      if (write_pot == "none" .and. write_psi == "none") then
        CALL DWGBUT("No output is selected for this simulation.|Do you want to continue?", INFO)
        if (INFO == 1) then
          CALL MAIN_ALGO
        end if
      else
        CALL MAIN_ALGO
      endif
    else if (MainMenuCmd == 3) then
      CALL SHOW_GRAPH
    else if (MainMenuCmd == 4) then
      EXIT
    end if
    if (MainMenuCmd < 4) then
      MainMenuCmd = MainMenuCmd + 1
    end if
  END DO

END SUBROUTINE

SUBROUTINE SHOW_IN_GUI

  CHARACTER(80) :: vstr
  INTEGER :: mode, mode2, ctrlh
  REAL tv

  CALL SWGTIT('TDS Tools')
  CALL WGINI('FORM', hForm)
  if (OSV == 0) then
    CALL SWGFNT ("fixed", 10)  ! lunux
    ctrlh = 30
  else
    CALL SWGFNT ("fixed", 14)  ! windows
    ctrlh = 24
  end if

    ! Setup frames
  CALL SWGWIN(10, 10, 390, 40)
  CALL WGBAS(hForm, 'FORM', hData)
  CALL SWGWIN(10, 60, 390, 200)
  CALL WGBAS(hForm, 'FORM', hPsi)
  CALL SWGWIN(10, 270, 390, 70)
  CALL WGBAS(hForm, 'FORM', hTime)
  CALL SWGWIN(10, 355, 390, 150)
  CALL WGBAS(hForm, 'FORM', hGrid)

  CALL SWGWIN(440, 10, 250, 40)
  CALL WGBAS(hForm, 'FORM', hData2)
  CALL SWGWIN(440, 50, 250, 240)
  CALL WGBAS(hForm, 'FORM', hPot)
  CALL SWGWIN(440, 300, 250, 220)
  CALL WGBAS(hForm, 'FORM', hOut)

    ! Graphics around frames
  CALL SWGWIN(10, 40, 390, 10)
  CALL WGLAB(hForm, "---------------------------------------", hNull)
  CALL SWGWIN(440, 40, 250, 10)
  CALL WGLAB(hForm, "---------------------------------------", hNull)
  CALL SWGWIN(10, 250, 390, 10)
  CALL WGLAB(hForm, "---------------------------------------", hNull)
  CALL SWGWIN(10, 342, 390, 10)
  CALL WGLAB(hForm, "---------------------------------------", hNull)
  CALL SWGWIN(440, 292, 250, 10)
  CALL WGLAB(hForm, "---------------------------------------", hNull)

  CALL SWGWIN(20, 520, 50, 28)
  CALL WGOk(hForm, hOk)
  CALL SWGWIN(80, 520, 50, 28)
  CALL WGQUIT(hForm, hQuit)

    ! El Mass
  CALL SWGWIN(0, 0, 95, 30)
  CALL WGLAB(hData, 'Electron Mass:', hNull)
  CALL SWGWIN(100, 0, 90, ctrlh)
  write (vstr, '(ES11.5)') electronmass
  CALL WGTXT(hData, vstr, hElMass)

  CALL SWGWIN(0, 0, 145, 30)
  CALL WGLAB(hData2, '3D scattering length:', hNull)
  CALL SWGWIN(150, 0, 90, ctrlh)
  write (vstr, '(ES11.5)') nonlin_as
  CALL WGTXT(hData2, vstr, hNonlin)

    ! Wave Function Frame
  CALL SWGWIN(0, 0, 80, 35)
  CALL WGLAB(hPsi, "WAVE FUNCTION", hNull)
  CALL SWGWIN(110, 0, 40, 35)
  CALL WGLAB(hPsi, "Mode:", hNull)
  if (psi_mode == "file") then
    mode = 2
  else
    mode = 1
  endif
  CALL SWGWIN(150, 0, 100, 33)
  CALL WGDLIS(hPsi, "Gaussian|File", mode, hPsiMode)
  write (vstr, '(ES11.5)') x0
  CALL SWGWIN(0, 40, 120, ctrlh)
  CALL WGLTXT(hPsi, 'x0', TRIM(vstr), 70, hPsiX0)
  write (vstr, '(ES11.5)') y0
  CALL SWGWIN(130, 40, 120, ctrlh)
  CALL WGLTXT(hPsi, 'y0', TRIM(vstr), 70, hPsiY0)
  write (vstr, '(ES11.5)') z0
  CALL SWGWIN(250, 40, 120, ctrlh)
  CALL WGLTXT(hPsi, 'z0', TRIM(vstr), 70, hPsiZ0)
  write (vstr, '(ES11.5)') sigmax
  CALL SWGWIN(0, 75, 120, ctrlh)
  CALL WGLTXT(hPsi, 'SigX', TRIM(vstr), 70, hPsiSigX)
  write (vstr, '(ES11.5)') sigmay
  CALL SWGWIN(130, 75, 120, ctrlh)
  CALL WGLTXT(hPsi, 'SigY', TRIM(vstr), 70, hPsiSigY)
  write (vstr, '(ES11.5)') sigmaz
  CALL SWGWIN(250, 75, 120, ctrlh)
  CALL WGLTXT(hPsi, 'SigZ', TRIM(vstr), 70, hPsiSigZ)
  write (vstr, '(ES11.5)') xenergy
  CALL SWGWIN(0, 110, 120, ctrlh)
  CALL WGLTXT(hPsi, 'Xnrg', TRIM(vstr), 70, hPsiXenergy)
  write (vstr, '(ES11.5)') yenergy
  CALL SWGWIN(130, 110, 120, ctrlh)
  CALL WGLTXT(hPsi, 'Ynrg', TRIM(vstr), 70, hPsiYenergy)
  write (vstr, '(ES11.5)') zenergy
  CALL SWGWIN(250, 110, 120, ctrlh)
  CALL WGLTXT(hPsi, 'Znrg', TRIM(vstr), 70, hPsiZenergy)
  CALL SWGWIN(0, 145, 37, 33)
  CALL WGLAB(hPsi, "File:", hNull)
  CALL SWGWIN(37, 145, 213, ctrlh)
  CALL WGFIL (hPsi, 'Select Initial Wave File', psi_file_in, '*.dat', hPsiFile)

    ! Time Frame
  write (vstr, '(ES11.5)') dt
  CALL SWGWIN(0, 0, 160, ctrlh)
  CALL WGLTXT(hTime, 'Time Step:', TRIM(vstr), 55, hTimeDelta)
  write (vstr, '(I10)') MAXIT
  CALL SWGWIN(0, 35, 160, ctrlh)
  CALL WGLTXT(hTime, 'Num Steps:', TRIM(vstr), 55, hTimeSteps)

    ! Grid Frame
  CALL SWGWIN(0, 0, 80, 33)
  CALL WGLAB(hGrid, "GRID", hNull)
  write (vstr, '(I5)') numx
  CALL SWGWIN(0, 40, 120, ctrlh)
  CALL WGLTXT(hGrid, 'NumX', TRIM(vstr), 70, hGridNumX)
  write (vstr, '(I5)') numy
  CALL SWGWIN(130, 40, 120, ctrlh)
  CALL WGLTXT(hGrid, 'NumY', TRIM(vstr), 70, hGridNumY)
  write (vstr, '(I5)') numz
  CALL SWGWIN(250, 40, 120, ctrlh)
  CALL WGLTXT(hGrid, 'NumZ', TRIM(vstr), 70, hGridNumZ)
  write (vstr, '(ES11.5)') size_x
  CALL SWGWIN(0, 75, 120, ctrlh)
  CALL WGLTXT(hGrid, 'SizeX', TRIM(vstr), 70, hGridSizeX)
  write (vstr, '(ES11.5)') size_y
  CALL SWGWIN(130, 75, 120, ctrlh)
  CALL WGLTXT(hGrid, 'SizeY', TRIM(vstr), 70, hGridSizeY)
  write (vstr, '(ES11.5)') size_z
  CALL SWGWIN(250, 75, 120, ctrlh)
  CALL WGLTXT(hGrid, 'SizeZ', TRIM(vstr), 70, hGridSizeZ)

    ! Potential Frame
  CALL SWGWIN(0, 0, 80, 33)
  CALL WGLAB(hPot, "POTENTIAL", hNull)
  CALL SWGWIN(0, 10, 40, 33)
  CALL WGLAB(hPot, "File:", hNull)
  CALL SWGWIN(40, 10, 200, ctrlh)
  CALL WGFIL(hPot, 'Select Potential File', pot_file_in, '*.dat', hPotFile)
  CALL SWGWIN(0, 45, 40, 33)
  CALL WGLAB(hPot, "Filelist:", hNull)
  CALL SWGWIN(40, 45, 200, ctrlh)
  CALL WGFIL(hPot, 'Select Pot Filelist', pot_filelist_name, '*.dat', hPotFilelist)

  CALL strtopot_countcmds(strpotentialX, mode)
  write (vstr, '(I2)') mode
  CALL SWGWIN(0, 80, 150, ctrlh)
  CALL WGLTXT(hPot, 'Cmd X num of lines:', vstr, 23, hPotStrX)
  CALL SWGWIN(155, 80, 40, ctrlh)
  CALL WGPBUT(hPot, 'Edit', hPotStrXEdit)
  CALL SWGCBK(hPotStrXEdit, PotStrCbk)

  CALL strtopot_countcmds(strpotentialy, mode)
  write (vstr, '(I2)') mode
  CALL SWGWIN(0, 115, 150, ctrlh)
  CALL WGLTXT(hPot, 'Cmd Y num of lines:', vstr, 23, hPotStrY)
  CALL SWGWIN(155, 115, 40, ctrlh)
  CALL WGPBUT(hPot, 'Edit', hPotStrYEdit)
  CALL SWGCBK(hPotStrYEdit, PotStrCbk)

  CALL strtopot_countcmds(strpotentialz, mode)
  write (vstr, '(I2)') mode
  CALL SWGWIN(0, 150, 150, ctrlh)
  CALL WGLTXT(hPot, 'Cmd Z num of lines:', vstr, 23, hPotStrZ)
  CALL SWGWIN(155, 150, 40, ctrlh)
  CALL WGPBUT(hPot, 'Edit', hPotStrZEdit)
  CALL SWGCBK(hPotStrZEdit, PotStrCbk)

  CALL strtopot_countcmds(strpotentialXYZ, mode)
  write (vstr, '(I2)') mode
  CALL SWGWIN(0, 185, 150, ctrlh)
  CALL WGLTXT(hPot, 'Cmd XYZ num of lines:', vstr, 23, hPotStrXYZ)
  CALL SWGWIN(155, 185, 40, ctrlh)
  CALL WGPBUT(hPot, 'Edit', hPotStrXYZEdit)
  CALL SWGCBK(hPotStrXYZEdit, PotStrCbk)

    ! Output Frame
  CALL SWGWIN(0, 0, 80, 33)
  CALL WGLAB(hOut, "OUTPUT", hNull)
  CALL SWGWIN(0, 30, 60, 33)
  CALL WGLAB(hOut, "Run name:", hNull)
  CALL SWGWIN(60, 30, 180, ctrlh)
  CALL WGTXT(hOut, TRIM(ADJUSTL(run_name)), hOutDir)
  CALL SWGWIN(0, 65, 100, 20)

  CALL SWGWIN(0, 65, 100, 20)
  CALL WGLAB(hOut, "Write Grid:", hNull)
  if (write_grid == 'txt' .or. write_grid == 'both') then
    mode = 1
  else
    mode = 0
  endif
  CALL SWGWIN(100, 65, 55, 20)
  CALL WGBUT(hOut, 'Txt', mode, hOutWGridTxt)
  if (write_grid == 'bin' .or. write_grid == 'both') then
    mode = 1
  else
    mode = 0
  endif
  CALL SWGWIN(160, 65, 55, 20)
  CALL WGBUT(hOut, 'Bin', mode, hOutWGridBin)

  CALL WGLAB(hOut, "Write Pot", hNull)
  if (write_pot == 'txt' .or. write_pot == 'both') then
    mode = 1
  else
    mode = 0
  endif
  CALL SWGWIN(100, 90, 55, 20)
  CALL WGBUT(hOut, 'Txt', mode, hOutWPotTxt)
  if (write_pot == 'bin' .or. write_pot == 'both') then
    mode = 1
  else
    mode = 0
  endif
  CALL SWGWIN(160, 90, 55, 20)
  CALL WGBUT(hOut, 'Bin', mode, hOutWPotBin)
  CALL SWGWIN(0, 115, 90, 20)
  CALL WGLAB(hOut, "Write Psi", hNull)
  if (write_psi == 'txt' .or. write_psi == 'both') then
    mode = 1
  else
    mode = 0
  endif
  CALL SWGWIN(100, 115, 55, 20)
  CALL WGBUT(hOut, 'Txt', mode, hOutWPsiTxt)
  if (write_psi == 'bin' .or. write_psi == 'both') then
    mode = 1
  else
    mode = 0
  endif
  CALL SWGWIN(160, 115, 55, 20)
  CALL WGBUT(hOut, 'Bin', mode, hOutWPsiBin)
  CALL SWGWIN(0, 140, 100, ctrlh)
  CALL WGLAB(hOut, 'Output Time Step:', hNull)
  write (vstr, '(ES11.5)') write_timestep
  CALL SWGWIN(100, 140, 90, 33)
  CALL WGTXT(hOut, TRIM(vstr), hOutTime)
  CALL SWGWIN(0, 175, 100, 30)
  CALL WGLAB(hOut, "Downsample:    X", hNull)
  write (vstr, '(I2)') write_downsample_x
  CALL SWGWIN(100, 175, 30, ctrlh)
  CALL WGTXT(hOut, TRIM(vstr), hOutDownX)
  CALL SWGWIN(140, 175, 10, 30)
  CALL WGLAB(hOut, "Y", hNull)
  write (vstr, '(I2)') write_downsample_y
  CALL SWGWIN(150, 175, 30, ctrlh)
  CALL WGTXT(hOut, TRIM(vstr), hOutDownY)
  CALL SWGWIN(180, 175, 10, 30)
  CALL WGLAB(hOut, "Z", hNull)
  write (vstr, '(I2)') write_downsample_z
  CALL SWGWIN(190, 175, 30, ctrlh)
  CALL WGTXT(hOut, TRIM(vstr), hOutDownZ)

    ! Displays the dialog
  CALL WGFIN

  CALL GWGTXT(hElMass, vstr)
  read (vstr, *) electronmass
  CALL GWGTXT(hNonlin, vstr)
  read (vstr, *) nonlin_as

    ! Get Wave Function data
  CALL GWGLIS(hPsiMode, mode)
  if (mode == 1) then
    psi_mode = 'gauss'
  else
    psi_mode = 'file'
  endif

  CALL GWGTXT(hPsiX0, vstr)
  read (vstr, *) x0
  CALL GWGTXT(hPsiY0, vstr)
  read (vstr, *) y0
  CALL GWGTXT(hPsiZ0, vstr)
  read (vstr, *) z0
  CALL GWGTXT(hPsiSigX, vstr)
  read (vstr, *) sigmax
  CALL GWGTXT(hPsiSigY, vstr)
  read (vstr, *) sigmay
  CALL GWGTXT(hPsiSigZ, vstr)
  read (vstr, *) sigmaz
  CALL GWGTXT(hPsiXenergy, vstr)
  read (vstr, *) xenergy
  CALL GWGTXT(hPsiYenergy, vstr)
  read (vstr, *) yenergy
  CALL GWGTXT(hPsiZenergy, vstr)
  read (vstr, *) zenergy
  CALL GWGFIL(hPsiFile, psi_file_in)

    ! Get Time Frame
  CALL GWGTXT(hTimeDelta, vstr)
  read (vstr, *) dt
  CALL GWGINT(hTimeSteps, MAXIT)

    ! Get Grid Data
  CALL GWGINT(hGridNumX, numx)
  CALL GWGINT(hGridNumY, numy)
  CALL GWGINT(hGridNumZ, numz)
  CALL GWGTXT(hGridSizeX, vstr)
  read (vstr, *) size_x
  CALL GWGTXT(hGridSizeY, vstr)
  read (vstr, *) size_y
  CALL GWGTXT(hGridSizeZ, vstr)
  read (vstr, *) size_z

    ! Potential Frame
  CALL GWGFIL(hPotFile, pot_file_in)
  CALL GWGFIL(hPotFileList, pot_filelist_name)
  
    ! Output Data
  CALL GWGTXT(hOutDir, write_folder)

  CALL GWGBUT(hOutWGridTxt, mode)
  CALL GWGBUT(hOutWGridBin, mode2)
  SELECT CASE (mode*2 + mode2)
  CASE (0)
    write_grid = 'none'
  CASE (1)
    write_grid = 'bin'
  CASE (2)
    write_grid = 'txt'
  CASE (3)
    write_grid = 'both'
  END SELECT

  CALL GWGBUT(hOutWPotTxt, mode)
  CALL GWGBUT(hOutWPotBin, mode2)
  SELECT CASE (mode*2 + mode2)
  CASE (0)
    write_pot = 'none'
  CASE (1)
    write_pot = 'bin'
  CASE (2)
    write_pot = 'txt'
  CASE (3)
    write_pot = 'both'
  END SELECT
  CALL GWGBUT(hOutWPsiTxt, mode)
  CALL GWGBUT(hOutWPsiBin, mode2)
  SELECT CASE (mode*2 + mode2)
  CASE (0)
    write_psi = 'none'
  CASE (1)
    write_psi = 'bin'
  CASE (2)
    write_psi = 'txt'
  CASE (3)
    write_psi = 'both'
  END SELECT
  CALL GWGTXT(hOutTime, vstr)
  read (vstr, *) write_timestep
  CALL GWGINT(hOutDownX, write_downsample_x)
  CALL GWGINT(hOutDownY, write_downsample_y)
  CALL GWGINT(hOutDownZ, write_downsample_z)

  call INDATA_SAVE(nml_file_name)

END SUBROUTINE


SUBROUTINE PotStrCbk(ID)
  INTEGER, INTENT(IN) :: ID
  INTEGER hFormPot, ii, num, hCmd(99), ctrlh
  CHARACTER(500) strcmd
  CHARACTER(999) fullcmd

  if (ID == hPotStrXEdit) then
    CALL GWGINT(hPotStrX, num)
    CALL SWGTIT('Potential X Commands')
  else if (ID == hPotStrYEdit) then
    CALL GWGINT(hPotStrY, num)
    CALL SWGTIT('Potential Y Commands')
  else if (ID == hPotStrZEdit) then
    CALL GWGINT(hPotStrZ, num)
    CALL SWGTIT('Potential Z Commands')
  else
    CALL GWGINT(hPotStrXYZ, num)
    CALL SWGTIT('Potential XYZ Commands')
  end if
  CALL WGINI('FORM', hFormPot)
  if (OSV == 0) then
    CALL SWGFNT ("fixed", 10)  ! lunux
    ctrlh = 30
  else
    CALL SWGFNT ("fixed", 14)  ! windows
    ctrlh = 24
  end if

  DO ii = 1,num
    if (ID == hPotStrXEdit) then
      CALL strtopot_extractcmd(strcmd, strpotentialX, ii)
    else if (ID == hPotStrYEdit) then
      CALL strtopot_extractcmd(strcmd, strpotentialY, ii)
    else if (ID == hPotStrZEdit) then
      CALL strtopot_extractcmd(strcmd, strpotentialZ, ii)
    else
      CALL strtopot_extractcmd(strcmd, strpotentialXYZ, ii)
    end if
    CALL SWGWIN(0, 32*(ii-1), 400, ctrlh)
    CALL WGTXT(hFormPot, TRIM(ADJUSTL(strcmd)), hCmd(ii))
  END DO
  CALL SWGWIN(0, 4 + 32*num, 50, 28)
  CALL WGOK(hFormPot, hOk)

  CALL WGFIN

  fullcmd = ""
  DO ii = 1,num
    CALL GWGTXT(hCmd(ii), strcmd)
    if (ii < num) then
      fullcmd = TRIM(fullcmd) // TRIM(strcmd) // ';'
    else
      fullcmd = TRIM(fullcmd) // TRIM(strcmd)
    end if
  END DO

  if (ID == hPotStrXEdit) then
    strpotentialX = fullcmd
  else if (ID == hPotStrYEdit) then
    strpotentialY = fullcmd
  else if (ID == hPotStrZEdit) then
    strpotentialZ = fullcmd
  else
    strpotentialXYZ = fullcmd
  end if

END SUBROUTINE

END MODULE
