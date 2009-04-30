!$$$$$$$$$$$$$$$$$$$$$$$$$$  stic2D  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_strtopot
  USE mod_staticdata
  IMPLICIT  NONE
  SAVE
  
CONTAINS

SUBROUTINE strtopot1d( npotloop, mstar, strpotential, nump, p, pot )
  IMPLICIT NONE
! Creates the 1d structure potential pot (in eV), reading values in
! eV from the definition string strpotential
!
! keywords:
! const    from  to  val
! linear   from  to  val_from  val_to
! poly3    from  to  val_from  val_to
! poly5    from  to  val_from  val_to
! harmonic from  to  hbaromega potoffset (eV)
!
! each argument "loop" is substituted with 
!       potloopvar_from + (potloopvar_step * npotloop)
!


INTEGER, INTENT(IN) :: npotloop     ! n. of cycle loop on potential
CHARACTER(999), INTENT(INOUT) :: strpotential
INTEGER, INTENT(IN) :: nump
REAL*8, INTENT(IN) :: mstar
REAL*8, INTENT(in) :: p(1:nump)
REAL*8, INTENT(OUT) :: pot(1:nump) ! the potential

CHARACTER(200) :: strcommand, strkeyword, strvalues
INTEGER :: lenstr
REAL*8 :: potdefvals(10)
REAL*8 :: paux
INTEGER :: np, nn, nnp


pot(:)= 0.
lenstr= LEN_TRIM(strpotential)
IF ( lenstr > 0 ) THEN
  IF ( strpotential(lenstr:lenstr) /= ";" ) THEN
    lenstr= lenstr+1
    strpotential(lenstr:lenstr)= ";"
  END IF
END IF
!PRINT*, '---------------------------'

nnp=1
DO nn= 1, lenstr
  IF ( strpotential(nn:nn)==";" ) THEN
!    PRINT*, " --- ", nn
    strcommand= REPEAT(" ",LEN(strcommand))
    strcommand= ADJUSTL( strpotential(nnp:nn-1) )
!    PRINT*, strcommand
    strkeyword= strcommand( : SCAN(strcommand," ") )
!    PRINT*, strcommand
    strvalues= strcommand( SCAN(strcommand," ")+1 : )
    nnp=nn+1
    SELECT CASE (strkeyword)
    CASE ("const", "CONST", "constant", "CONSTANT")
      CALL strtopot_valsread( npotloop, strvalues, 3, potdefvals(:) )
      ! PRINT*, "ccc", potdefvals(:3)
      DO np= 1, nump
        IF (p(np)>=potdefvals(1) .AND. p(np)<=potdefvals(2)) THEN
          pot(np)= potdefvals(3)
        END IF
      END DO
    CASE ("linear", "LINEAR")
      CALL strtopot_valsread( npotloop, strvalues, 4, potdefvals(:) )
      ! PRINT*, "lll", potdefvals(:4)
      DO np= 1, nump
        IF (p(np)>=potdefvals(1) .AND. p(np)<=potdefvals(2)) THEN
          pot(np)= potdefvals(3) + (potdefvals(4)-potdefvals(3))*     &
               &  (p(np)-potdefvals(1))/(potdefvals(2)-potdefvals(1))
        END IF
      END DO
    CASE ("poly3", "POLY3")
      CALL strtopot_valsread( npotloop, strvalues, 4, potdefvals(:) )
      ! PRINT*, "ppp", potdefvals(:4)
      IF (potdefvals(1) >= potdefvals(2)) STOP "mk1dpot: values error"
      DO np= 1, nump
        IF (p(np)>=potdefvals(1) .AND. p(np)<=potdefvals(2)) THEN
          paux= (p(np)-potdefvals(1)) / (potdefvals(2)-potdefvals(1))
          pot(np)= potdefvals(3) + (potdefvals(4)-potdefvals(3))*     &
               &  ( 3.*(paux**2) - 2.*(paux**3) )
        END IF
      END DO
    CASE ("poly5", "POLY5")
      CALL strtopot_valsread( npotloop, strvalues, 4, potdefvals(:) )
      ! PRINT*, "ppp", potdefvals(:4)
      IF (potdefvals(1) >= potdefvals(2)) STOP "mk1dpot: values error"
      DO np= 1, nump
        IF (p(np)>=potdefvals(1) .AND. p(np)<=potdefvals(2)) THEN
          paux= (p(np)-potdefvals(1)) / (potdefvals(2)-potdefvals(1))
          pot(np)= potdefvals(3) + (potdefvals(4)-potdefvals(3))*     &
               &  ( 10.*(paux**3) - 15.*(paux**4) + 6.*(paux**5) )
        END IF
      END DO

    CASE ("harmonic", "HARMONIC")
      CALL strtopot_valsread( npotloop, strvalues, 4, potdefvals(:) )
      ! PRINT*, "ppp", potdefvals(:4)
      IF (potdefvals(1) >= potdefvals(2)) STOP "mk1dpot: values error"
      DO np= 1, nump
        IF (p(np)>=potdefvals(1) .AND. p(np)<=potdefvals(2)) THEN
          paux= p(np) - potdefvals(1) - (potdefvals(2)-potdefvals(1))/2.
! TODO: controllare
          pot(np)= (mstar/2.)*((potdefvals(3)/HBAR_eV)**2) / ELCH
          pot(np)= pot(np)*(paux**2) + potdefvals(4)
        END IF
      END DO
    CASE DEFAULT 
      IF ( strkeyword(1:1)/="#" ) THEN
        print *, strpotential
        STOP "mk1dpot: unrecognised keyword"
      END IF
    END SELECT
  END IF
END DO

END SUBROUTINE strtopot1d

!************************************************************************
SUBROUTINE strtopot2d( npotloop, mstar, strpotential, numpx, numpy, px, py, pot )
  IMPLICIT NONE
! Creates the 2d structure potential pot (in eV), reading values in
! eV from the definition string strpotential
!
! keywords: 
! box      xLowerLeft  yLowerLeft  xUpperRight yUpperRight  val(eV)
!
! each argument "loop" is substituted with 
!       potloopvar_from + (potloopvar_step * npotloop)
!


INTEGER, INTENT(IN) :: npotloop     ! n. of cycle loop on potential
CHARACTER(999), INTENT(INOUT) :: strpotential
INTEGER, INTENT(IN) :: numpx, numpy
REAL*8, INTENT(IN) :: mstar
REAL*8, INTENT(IN) :: px(1:numpx), py(1:numpy)
REAL*8, INTENT(OUT) :: pot(1:numpy, 1:numpx)  ! the potential

CHARACTER(200) :: strcommand, strkeyword, strvalues
INTEGER :: lenstr
REAL*8 :: potdefvals(10)
!REAL*8 :: paux
INTEGER :: npx, npy, nn, nnp


pot(:,:)= 0.

lenstr= LEN_TRIM(strpotential)
IF ( lenstr > 0 ) THEN
  !PRINT*, strpotential
  IF ( strpotential(lenstr:lenstr) /= ";" ) THEN
    lenstr= lenstr+1
    strpotential(lenstr:lenstr)= ";"
  END IF
END IF
!PRINT*, '---------------------------', lenstr

nnp=1
DO nn= 1, lenstr
  IF ( strpotential(nn:nn)==";" ) THEN
    !PRINT*, " --- ", nn
    strcommand= REPEAT(" ",LEN(strcommand))
    strcommand= ADJUSTL( strpotential(nnp:nn-1) )
    !PRINT*, strcommand
    strkeyword= strcommand( : SCAN(strcommand," ") )
    !PRINT*, strcommand
    strvalues= strcommand( SCAN(strcommand," ")+1 : )
    nnp=nn+1
    SELECT CASE (strkeyword)
    CASE ("box", "BOX")
      CALL strtopot_valsread( npotloop, strvalues, 5, potdefvals(:) )
      DO npx= 1, numpx
        IF (px(npx)>=potdefvals(1) .AND. px(npx)<=potdefvals(3)) THEN
          DO npy= 1, numpy
            IF (py(npy)>=potdefvals(2) .AND. py(npy)<=potdefvals(4)) THEN
              pot(npy,npx)= potdefvals(5)
            END IF
          END DO
        END IF
      END DO
    CASE DEFAULT 
      IF ( strkeyword(1:1)/="#" ) THEN
        print *, strpotential
        STOP "strtopot2d: unrecognized keyword"
      END IF
    END SELECT
  END IF
END DO

END SUBROUTINE strtopot2d

!************************************************************************
SUBROUTINE strtopot_valsread( npotloop, strvalues, numvalues, potdefvals )
  IMPLICIT NONE
! 

INTEGER, INTENT(IN) :: npotloop      ! n. of cycle loop on potential
CHARACTER(*), INTENT(IN) :: strvalues
INTEGER, INTENT(IN) :: numvalues
REAL*8, INTENT(OUT) :: potdefvals(1:)

CHARACTER(200) :: straux, strval
INTEGER :: readerror, npos1, nv


IF ( numvalues > SIZE(potdefvals(:)) ) stop "strtopot_valsread: numvalues"

straux= strvalues
DO nv= 1, numvalues
  straux= ADJUSTL(straux)
  npos1= SCAN(straux," ")
  IF (npos1==0) STOP "strtopot_valsread: not enough values"
  strval= straux( 1:npos1 )
  straux( 1:npos1 )= REPEAT(" ",npos1)
!  IF (strval=="loop" .OR. strval=="LOOP") THEN
!    potdefvals(nv)= potloopvar_from + (potloopvar_step * npotloop)
    ! potloopvar= potdefvals(nv)
!  ELSE
    READ(strval,*,IOSTAT=readerror) potdefvals(nv)
    IF (readerror/=0) STOP "strtopot_valsread: Err in parsing pot values"
!  END IF
END DO
IF (straux/=" ") STOP "strtopot_valsread: too many values"

END SUBROUTINE strtopot_valsread


SUBROUTINE strtopot_countcmds(strpot, num_cmd)
  CHARACTER(*), INTENT(IN) :: strpot
  INTEGER, INTENT(OUT) :: num_cmd

  INTEGER lenstr, nn

  num_cmd = 0
  lenstr= LEN_TRIM(strpot)
  IF ( strpot(lenstr:lenstr) /= ";" ) THEN
    num_cmd = 1
  END IF

  DO nn= 1, lenstr
    IF ( strpot(nn:nn) == ';' ) THEN
      num_cmd = num_cmd + 1
    END IF
  END DO

END SUBROUTINE strtopot_countcmds


! Extracts the num_cmd-th  command from strpot and put it into str_single_cmd
SUBROUTINE strtopot_extractcmd(str_single_cmd, strpot, num_cmd)
  CHARACTER(*), INTENT(INOUT) :: strpot
  CHARACTER(*), INTENT(OUT) :: str_single_cmd
  INTEGER, INTENT(IN) :: num_cmd

  INTEGER :: lenstr, nn, nnp, sent

  lenstr= LEN_TRIM(strpot)
  IF ( strpot(lenstr:lenstr) /= ";" ) THEN
    lenstr= lenstr+1
    strpot(lenstr:lenstr)= ";"
  END IF

  str_single_cmd = ""
  nnp = 1
  sent = 1
  DO nn = 1, lenstr
    IF ( strpot(nn:nn)==";" ) THEN
      IF (sent == num_cmd) THEN
        str_single_cmd = TRIM(strpot(nnp:nn-1))
        RETURN
      END IF
      sent = sent + 1
      nnp = nn+1
    END IF
  END DO

END SUBROUTINE

END MODULE mod_strtopot
