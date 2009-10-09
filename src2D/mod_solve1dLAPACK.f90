!$$$$$$$$$$$$$$$$$$$$$$$$$$$$  stic2P1D  $$$$$$$$$$$$$$$$$$$$$$$$$$

MODULE mod_solve1dLAPACK
  IMPLICIT  NONE
  SAVE
  
CONTAINS

SUBROUTINE solve1dLAPACK( numx,numevals,numevecs,const,potx,energies,psi )
  IMPLICIT NONE
! Solves the 1-electron 1-dimension Schr eq 
! with closed boundary with potential pot
! the eigenvectors evec are orthonormalised (norm=1)
! The |constant| that multiply the numerical II derivative is given in const
! usually: const= (hbar^2 / (2 emme dx^2)).
! Uses LAPACK double prec lib

INTEGER, INTENT(IN) :: numx        ! num of points (0:numx)
INTEGER, INTENT(IN) :: numevals    ! num of eigenenergies to be computed
INTEGER, INTENT(IN) :: numevecs    ! num of eigenfunctions to be computed
REAL*8, INTENT(IN) :: const        ! const in front of the II derivative
REAL*8, INTENT(IN) :: potx(0:numx) ! the 1d potential
REAL*8, INTENT(OUT) :: energies(numevals)       ! computed eigenvalues
REAL*8, INTENT(OUT) :: psi(0:numx,numevecs)     ! computed eigenvectors

REAL*8 :: potxconst(0:numx)
REAL*8 :: cnorm
INTEGER :: nx, nn

REAL*8, ALLOCATABLE :: d_m(:)         ! matrix diagonal
REAL*8, ALLOCATABLE :: e_m(:)         ! matrix super/sub diagonal
REAL*8, ALLOCATABLE :: w_m(:)           ! eigenvalues
REAL*8, ALLOCATABLE :: z_m(:,:)         ! eigenvectors
REAL*8, ALLOCATABLE :: work_m(:)
INTEGER, ALLOCATABLE :: iwork_m(:)
INTEGER, ALLOCATABLE :: ifail_m(:)
INTEGER :: n_m, numevalsfound
INTEGER :: info_m
REAL*8 :: abstol

REAL*8, EXTERNAL :: DLAMCH            ! DOUBLE PRECISION machine parameters

!............................................................initializations
potxconst(:)= potx(:)/const

n_m= numx+1
abstol= 2.*DLAMCH('S')
!lwork_m= 16*n_m

ALLOCATE( d_m(n_m) )
ALLOCATE( e_m(n_m) )
ALLOCATE( w_m(n_m) )
ALLOCATE( z_m(n_m,numevals) )
ALLOCATE( work_m(5*n_m) )
ALLOCATE( iwork_m(5*n_m) )
ALLOCATE( ifail_m(n_m) )

!......................................................creates the matrix
DO nx= 0, numx
  d_m( nx+1 )= 2. + potxconst(nx)
END DO
DO nx= 0, numx
  e_m( nx+1 )= -1.
END DO


CALL DSTEVX( 'V', 'I', n_m, d_m, e_m, 0d0, 0d0, 1, numevals,         &
     &       abstol, numevalsfound, w_m, z_m, n_m,                   &
     &       work_m, iwork_m, ifail_m, info_m  )

IF (info_m /= 0) THEN
  PRINT*, "mod_solve1dLAPACK: in dseupd info=", info_m
  STOP
END IF
IF (numevalsfound < MIN(numevals,numevecs)) THEN
  PRINT*, "mod_solve1dLAPACK: numevalsfound=", numevalsfound
  STOP
END IF

!..................................................fills energies and psi
energies(1:numevals)= w_m(1:numevals) * const

psi(:,:)= 0.
DO nn= 1, numevecs
  DO nx= 0, numx
    psi(nx,nn)= z_m(nx+1,nn) 
  END DO
cnorm= SUM(psi(:,nn)**2)
psi(:,nn)= psi(:,nn) / SQRT(cnorm)
END DO

!............................................................finalizations
DEALLOCATE( d_m )
DEALLOCATE( e_m )
DEALLOCATE( w_m )
DEALLOCATE( z_m )
DEALLOCATE( work_m )
DEALLOCATE( iwork_m )
DEALLOCATE( ifail_m )

END SUBROUTINE solve1dLAPACK


END MODULE mod_solve1dLAPACK
