!$$$$$$$$$$$$$$$$$$$$$$$$$$$$  stio3D_qtbm  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_staticdata

  IMPLICIT  NONE
  SAVE
  COMPLEX*16, PARAMETER :: i=(0.,1.)

! physics constants   (in MKS)
  REAL*8, PARAMETER  :: ELMASS0= 9.1095e-31,                       &
       &                ELCH= 1.602176e-19,                        &
       &                HBAR= 1.05459e-34,                         &
       &                HBAR_eV= HBAR / ELCH,                      &
       &                BOLTZ= 1.3807e-23,                         &
       &                PIG=3.14159265,                            &
       &                EPSILON0= 8.854e-12
  
  INTEGER, PARAMETER :: OSV = 1   ! OS Version: 0 = linux, 1 = Windows

END MODULE mod_staticdata
