!------------------------------------------------------------------------------
! MODULE: Modulename
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains some general parameters critical to controlling
!!the code, grid variables and some mathematical and physical constants 
!!used everywhere.
!------------------------------------------------------------------------------

MODULE paramr

IMPLICIT NONE

LOGICAL,PARAMETER :: tpairon=.true.      !<Switch for using DDDI contrib. to mean field (strongly recommended).
LOGICAL,PARAMETER :: tadsic=.false.      !<Remove self-interaction in Coulomb term
REAL(8),PARAMETER :: einitbubble=0.0D0   !<Tries to drive bubble by adding a -m**2 potential for a few initial iterations.
INTEGER,PARAMETER :: nzd= 200            !<Max. no. of mesh points in z direction.
INTEGER,PARAMETER :: nrd= 100            !<Max. no. of mesh points in r direction.
INTEGER,PARAMETER :: nsd=2               !<No. of spinor components.
INTEGER,PARAMETER :: nhwd=300            !<Maximum no. of oscillator quanta.
!INTEGER,PARAMETER :: nhwd=120            !<Maximum no. of oscillator quanta.
INTEGER,PARAMETER :: kqm=  2             !<No. of isospins (standard: 2).
!INTEGER,PARAMETER :: kmpsi=400            !<Maximal number of mpsi for damping in r-direction.
INTEGER,PARAMETER :: kmpsi=50            !<Maximal number of mpsi for damping in r-direction.
INTEGER,PARAMETER :: klocal=0            !<Switch to print localization criterion.
INTEGER,PARAMETER :: ifoccs=0            !<Switch to print occ.numbers.

!Derived parameters
INTEGER,PARAMETER :: mxd=nhwd+1          !<Maximum no. of allowed m values (=nhwd+1).
INTEGER,PARAMETER :: nnzxd=nhwd          !<Maximum principle quantum number for the cylindrical oscillator functions (z-direction) (nhwd)
INTEGER,PARAMETER :: nnrxd=(nhwd+1)/2    !<Maximum principle quantum number for the cylindrical oscillator functions (r-direction) ((nhwd+1)/2).

!**********************************************************************
!                                                                     *
!        physical and mathematical constants                          *
!                                                                     *
!**********************************************************************
REAL(8),PARAMETER :: zero=0D0            !<0.0
REAL(8),PARAMETER :: one=1D0             !<1.0
REAL(8),PARAMETER :: two=2D0             !<2.0
REAL(8),PARAMETER :: xxxl=HUGE(1D0)      !<Largest REAL(8).
REAL(8),PARAMETER :: xsmall=1D-30        ! some small value
REAL(8),PARAMETER :: pi=4.0d0*atan(1.0d0)!<The constant \f$ \pi \f$.
REAL(8),PARAMETER :: edr=1.0d0/3.0d0     !<1/3
REAL(8),PARAMETER :: oh=0.5D0            !<1/2
REAL(8),PARAMETER :: eighty=80D0         !<80.0
REAL(8),PARAMETER :: smal=1D-30          !<Small real number.
INTEGER,PARAMETER :: ixxxl=999999        !<Huge integer to signal problems.

REAL(8),PARAMETER :: hc=197.3285879D0    !<The constant \f$ \hbar c \f$ in units of MeV*fm.
REAL(8),PARAMETER :: e2=1.43997577D0     !<Electron charge squared.
REAL(8),PARAMETER :: r0=1.4D00           !<Pion range.

!**********************************************************************
!                                                                     *
!        global dynamical variables                                   *
!                                                                     *
!**********************************************************************
INTEGER,PARAMETER :: itzpe=50               !<Controls frequency the c.m. correction is evaluated. 
INTEGER :: itest=0                       !<Controlls the test output frequency.


CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: sqrtb
!> @brief
!!Square root with minimal value \c smal.
!>
!> @param[in] x
!> REAL(8), takes the input value.
!---------------------------------------------------------------------------
REAL(8) PURE FUNCTION sqrtb(x)

REAL(8), INTENT(IN) :: x

sqrtb = SQRT(MAX(x,smal))
RETURN

END FUNCTION sqrtb

END MODULE paramr
