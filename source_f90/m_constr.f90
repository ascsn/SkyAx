!------------------------------------------------------------------------------
! MODULE: constr
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routinesand variables for constraint iterations.
!------------------------------------------------------------------------------
MODULE constr

USE grid, ONLY : nz,nr,dz,dr,tsym

IMPLICIT NONE

INTEGER, PUBLIC :: itypco=-1      !< switches scheme for constrained iteration
!! and derives the 'QFAC' correction as well as 'qlamb' from it.
!! For details, see SUBROUTINE 'update_lambda()'.
REAL(8) :: qerr=1E-4      !<Accuracy of quadrupole deformation.
REAL(8) :: c0=0.8         !<The iteration parameter \f$ c_0 \f$ for constraint calculations.
REAL(8) :: d0=0D0         !<The stabilizing iteration parameter \f$ c_0 \f$ for constraint calculations.
REAL(8) :: qeps=-1D3      !<step size for update \f$ \lambda \f$ in constraint calculation.
REAL(8) :: dalp=3.5D0     !<Radius of Wood-Saxon profile to damp quadrupole operator at large distances.
REAL(8) :: dgam=0.5D0     !<Width of Wood-Saxon profile to damp quadrupole operator at large distances.
REAL(8) :: swifac=0.05D0  !<Defines the switching point from quadratic to linear shape. It is at 
!!\f$ \texttt{swifac} \cdot 1.2 \cdot A^{1/3}\f$ (only on z-direction).
LOGICAL :: tswicon=.FALSE. !< Switch for any constraint
LOGICAL :: tswicon0=.FALSE. !< Switches constraint calculations L=0
LOGICAL :: tfixcon0=.FALSE. !< Switches constraint calculations with fixed force L=0
REAL(8) :: swicon0=0D0    !<Switches constraint calculations L=0: 0=off, +1=inhibitd for L=0, -1=fix constraining force
REAL(8) :: betmin0=0D0    !<Minimal deformation parameter L=0.
REAL(8) :: betmax0=0D0    !<Maximal deformation parameter L=0.
REAL(8) :: betaw0=0.1D0   !<Deformation parameter step width L=0.
LOGICAL :: tswicon2=.FALSE. !< Switches constraint calculations L=2
LOGICAL :: tfixcon2=.FALSE. !< Switches constraint calculations with fixed force L=2
REAL(8) :: swicon2=0D0    !<Switches constraint calculations L=2: 0=off, +1=fix deformation, -1=fix constraining force
REAL(8) :: betmin2=0D0    !<Minimal deformation parameter L=2.
REAL(8) :: betmax2=0D0    !<Maximal deformation parameter L=2.
REAL(8) :: betaw2=0.1D0   !<Deformation parameter step width L=2.
LOGICAL :: tswicon3=.FALSE. !< Switches constraint calculations L=3
LOGICAL :: tfixcon3=.FALSE. !< Switches constraint calculations with fixed force L=3
REAL(8) :: swicon3=0D0    !<Switches constraint calculations L=3: 0=off, +1=fix deformation, -1=fix constraining force
REAL(8) :: betmin3=0D0    !<Minimal deformation parameter L=3.
REAL(8) :: betmax3=0D0    !<Maximal deformation parameter L=3.
REAL(8) :: betaw3=0.1D0   !<Deformation parameter step width L=3.
LOGICAL :: tswicon4=.FALSE. !< Switches constraint calculations L=4 
LOGICAL :: tfixcon4=.FALSE. !< Switches constraint calculations with fixed force L=4
REAL(8) :: swicon4=0D0    !<Switches constraint calculations L=4: 0=off, +1=fix deformation, -1=fix constraining force
REAL(8) :: betmin4=0D0    !<Minimal deformation parameter L=4.
REAL(8) :: betmax4=0D0    !<Maximal deformation parameter L=4.
REAL(8) :: betaw4=0.1D0   !<Deformation parameter step width L=4.
REAL(8) :: octfac         ! fore-factor on octupole moment
REAL(8) :: dipin=0D0      !<Strength of external dipole field along z-axis

INTEGER :: itcons0=99999  !<Maximum nr. constraint steps within one iteration (L=0)
INTEGER :: itcons2=99999  !<Maximum nr. constraint steps within one iteration (L=2)
INTEGER :: itcons3=99999  !<Maximum nr. constraint steps within one iteration (L=3)
INTEGER :: itcons4=99999  !<Maximum nr. constraint steps within one iteration (L=4)
INTEGER :: itcons=0       !<indicates whether ther is any dynamic constraint at all
LOGICAL :: tlin=.FALSE.   !<when true switches to linear asymptotic L=2 constraint.
LOGICAL :: tfresh=.FALSE. !<Switch for using new initialized wave function at each deformation. 

LOGICAL, ALLOCATABLE :: lrho(:,:)    !<protocol of changes in density
REAL(8), ALLOCATABLE :: distsf(:,:)  !<Field of distances to the chosen surface line.

REAL(8), ALLOCATABLE :: qfac0(:,:)   !<Damped constraing force L=0
REAL(8), ALLOCATABLE :: qfac2(:,:)   !<Damped constraing force L=2
REAL(8), ALLOCATABLE :: qfac3(:,:)   !<Damped constraing force L=3
REAL(8), ALLOCATABLE :: qfac4(:,:)   !<Damped constraing force L=4
REAL(8), ALLOCATABLE :: qcfac0(:,:)  !<Raw constraing force L=0
REAL(8), ALLOCATABLE :: qcfac2(:,:)  !<Raw constraing force L=2
REAL(8), ALLOCATABLE :: qcfac3(:,:)  !<Raw constraing force L=3
REAL(8), ALLOCATABLE :: qcfac4(:,:)  !<Raw constraing force L=4
REAL(8), ALLOCATABLE :: confield(:,:)!< total constraint contribution
REAL(8) :: qcst                      !< stores previous value of deformation
REAL(8) :: qsoll0=0.0d0              !< wanted values of deformation L=0
REAL(8),PARAMETER :: qsoll1=0.0d0     !< wanted values of deformation L=1 (c.m.)
REAL(8) :: qsoll2=0.0d0              !< wanted values of deformation L=2
REAL(8) :: qsoll3=0.0d0              !< wanted values of deformation L=3
REAL(8) :: qsoll4=0.0d0              !< wanted values of deformation L=4
REAL(8) :: qlamb0=0D0                !< actual Lagrange parameters L=0
REAL(8) :: qlamb1=0D0                !< actual Lagrange parameters L=1
REAL(8) :: qlamb2=0D0                !< actual Lagrange parameters L=2
REAL(8) :: qlamb3=0D0                !< actual Lagrange parameters L=3
REAL(8) :: qlamb4=0D0                !< actual Lagrange parameters L=4
REAL(8) :: qlambi0=0D0               !< intrinsic Lagrange parameters L=0 (for augmented Lagrangian method)
REAL(8) :: qlambi1=0D0               !< intrinsic Lagrange parameters L=1 (for augmented Lagrangian method)
REAL(8) :: qlambi2=0D0               !< intrinsic Lagrange parameters L=2 (for augmented Lagrangian method)
REAL(8) :: qlambi3=0D0               !< intrinsic Lagrange parameters L=3 (for augmented Lagrangian method)
REAL(8) :: qlambi4=0D0               !< intrinsic Lagrange parameters L=4 (for augmented Lagrangian method)
INTEGER :: nqmax0=0                  !< nr. of wanted deformations L=0
INTEGER :: nqmax2=0                  !< nr. of wanted deformations L=2
INTEGER :: nqmax3=0                  !< nr. of wanted deformations L=3
INTEGER :: nqmax4=0                  !< nr. of wanted deformations L=4
REAL(8) :: qmin0=0D0                 !< lower wanted deformation L=0
REAL(8) :: qmin2=0D0                 !< lower wanted deformation L=2
REAL(8) :: qmin3=0D0                 !< lower wanted deformation L=3
REAL(8) :: qmin4=0D0                 !< lower wanted deformation L=4
REAL(8) :: qmax0=0D0                 !< upper wanted deformation L=0
REAL(8) :: qmax2=0D0                 !< upper wanted deformation L=2
REAL(8) :: qmax3=0D0                 !< upper wanted deformation L=3
REAL(8) :: qmax4=0D0                 !< upper wanted deformation L=4
REAL(8) :: qw0=0D0                   !< step in wanted deformation L=0
REAL(8) :: qw2=0D0                   !< step in wanted deformation L=2
REAL(8) :: qw3=0D0                   !< step in wanted deformation L=3
REAL(8) :: qw4=0D0                   !< step in wanted deformation L=4
REAL(8) :: cfact0=0D0                !< conversion to dimensionless deformation L=0
REAL(8) :: cfact2=0D0                !< conversion to dimensionless deformation L=2
REAL(8) :: cfact3=0D0                !< conversion to dimensionless deformation L=3
REAL(8) :: cfact4=0D0                !< conversion to dimensionless deformation L=4

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: init_constraint
!> @brief
!!Initializes the constraint with Wood-Saxon damping at the outer surface. 
!> @details
!!The damping factor as such is stored in \c fda.
!!  - Option \c tlin=.FALSE.: quadrupole as such is damped
!!  - Option \c tlin=.TRUE.: linear asymptotic and softened damping.
!!
!!The linear asymtotic constraint is:
!!\f[ Q=c\cdot \sqrt{1+\frac{4z^2}{c}}-c-\frac{r^2}{\sqrt{1+4z^2/c}} \f]
!!where \f$ c=16\cdot z_{\rm switch}^2 \f$
!!and \c z_switch indicates the z-value where the quadratic behaviour switches 
!!to linear behaviour. \c z_switch is determined as 
!!\f$ \texttt{swifac} \cdot 1.2 \cdot A^{1/3}\f$ and \c swifac comes from 
!!input.
!--------------------------------------------------------------------------- 
SUBROUTINE init_constraint

! Note: adjustment of damping presently optimized to quadrupole only

USE paramr
USE levels, ONLY : npart
USE grid, ONLY :nz,nr,test_print,r,z

IMPLICIT NONE

REAL(8),SAVE :: qr=0D0                      ! to make it static
REAL(8),PARAMETER :: y40fac=0.846284375D0   ! 24/(16*SQRT(PI))

INTEGER :: iz,ir

REAL(8) :: a,besoll,cswitc,dterm,radd,radd1,radd2,rcd,rmin,&
           rmin2,rr,sqrtrr,sqswit,zswitc,dampfac
           
!initialisation of the constraint



IF(.NOT.ALLOCATED(confield)) THEN
  ALLOCATE(confield(0:nz,0:nr))
  confield=0D0
END IF

tswicon = tswicon0 .OR. tswicon2 .OR. tswicon3 .OR. tswicon4

itcons=MAX(itcons0,itcons2,itcons3,itcons4)

IF(.NOT.ALLOCATED(lrho)) THEN
  ALLOCATE(lrho(0:nz,0:nr))
  lrho=.FALSE.
END IF
IF(.NOT.ALLOCATED(distsf)) THEN
  ALLOCATE(distsf(0:nz,0:nr))
  distsf=0D0
END IF

IF(tswicon0 .AND. .NOT.ALLOCATED(qfac0)) THEN
  ALLOCATE(qfac0(0:nz,0:nr),qcfac0(0:nz,0:nr))
  qfac0=0D0
  qcfac0=0D0
END IF
IF(tswicon2 .AND. .NOT.ALLOCATED(qfac2)) THEN
  ALLOCATE(qfac2(0:nz,0:nr),qcfac2(0:nz,0:nr))
  qfac2=0D0
  qcfac2=0D0
END IF
IF(tswicon3 .AND. .NOT.ALLOCATED(qfac3)) THEN
  ALLOCATE(qfac3(0:nz,0:nr),qcfac3(0:nz,0:nr))
  qfac3=0D0
  qcfac3=0D0
END IF
IF(tswicon4 .AND. .NOT.ALLOCATED(qfac4)) THEN
  ALLOCATE(qfac4(0:nz,0:nr),qcfac4(0:nz,0:nr))
  qfac4=0D0
  qcfac4=0D0
END IF

! some factors for defining and converting multipoles
octfac = 0.5D0 * SQRT(7D0/pi)
a      = 1D0*(npart(1)+npart(2))
rcd    = 1.2D0*a**(edr)
besoll = SQRT(5D0*pi)*qsoll2/(3.0*rcd*rcd*a)
rmin   = rcd+dalp
rmin2  = rmin*rmin
radd1  = besoll*rcd*SQRT(5D0/pi)/(2.0)
zswitc = swifac*rcd
cswitc = 16D0*zswitc*zswitc
IF (qsoll2 > 0D0) THEN
  radd2  = SQRT((qsoll2)/(0.5D0*a) )*0.5D0
  IF (radd1 <= radd2) THEN
    radd = radd1
  ELSE
    radd = radd2
  END IF
ELSE IF (qsoll2 == zero) THEN
  radd  = 0D0
ELSE
  radd  =-0.5D0*radd1
END IF
dterm = rmin2-radd*radd
DO iz=0,nz
  IF (qsoll2 >= zero) THEN
    IF (ABS(z(iz)) >= radd) THEN
      qr         = SQRT(dterm+2.0*radd*ABS(z(iz)))
    ELSE
      qr         = SQRT(z(iz)**2+rmin2)
    END IF
  END IF
  DO ir=0,nr

    rr     = z(iz)**2+r(ir)**2
    sqrtrr = SQRT(rr)
    
!the damping functional
    IF (qsoll2 < zero) THEN
      IF (r(ir) >= radd) THEN
        qr = SQRT(dterm+2.0*radd*r(ir))
      ELSE
        qr = SQRT(rmin2+r(ir)**2)
      END IF
    END IF
    dampfac  = one/(one+EXP((sqrtrr-qr)/dgam)) 

! the bare constraining force on save field 'qcfac'
    
    IF(tswicon0) THEN
      qcfac0(iz,ir) = rr               ! monopole
      qfac0(iz,ir) = dampfac*qcfac0(iz,ir)
    END IF
    IF(tswicon2) THEN
      IF(tlin) THEN
        sqswit = SQRT(one+4*z(iz)**2/cswitc)
        qcfac2(iz,ir) = (cswitc*sqswit-cswitc-r(ir)**2/sqswit**3) ! asymp. linear
      ELSE
        qcfac2(iz,ir) = (2*z(iz)**2 - r(ir)**2)       ! cart. quadrupole
      END IF
      qfac2(iz,ir) = dampfac*qcfac2(iz,ir)
    END IF
    IF(tswicon3) THEN
      qcfac3(iz,ir) = octfac*z(iz)*(z(iz)**2 - 1.5*r(ir)**2)   ! spher. octupole
      qfac3(iz,ir) = dampfac*qcfac3(iz,ir)
    END IF
    IF(tswicon4) THEN
      qcfac4(iz,ir) = y40fac*(z(iz)**4-3.0*z(iz)**2*r(ir)**2+0.375*r(ir)**4)
      qfac4(iz,ir) = dampfac*qcfac4(iz,ir)
    END IF
  END DO
END DO
IF(itest > 0) THEN
  IF(tswicon0) CALL test_print(qfac0(0:nz,0:nr),0,nz,0,9,'qfac0')
  IF(tswicon2) CALL test_print(qfac2(0:nz,0:nr),0,nz,0,9,'qfac2')
  IF(tswicon3) CALL test_print(qfac3(0:nz,0:nr),0,nz,0,9,'qfac3')
  IF(tswicon4) CALL test_print(qfac4(0:nz,0:nr),0,nz,0,9,'qfac4')
END IF
RETURN
END SUBROUTINE init_constraint

!---------------------------------------------------------------------------  
! DESCRIPTION: prep_qsoll
!> @brief
!! Prepares pool of wanted deformations by translating dimensionless 
!! deformation to multipole moments. Determines also values for \c nz and 
!! \c nr.
!> @details
!! Convert deformation, using the relations:
!!   - \f$ Q_{\rm cartesian} = 1.09  \cdot A^{5/3}  \cdot \beta_2 \f$
!!   - \f$ Q_{40}            = 0.495 \cdot A^{7/3}  \cdot \beta_4 \f$
!!   - \f$ Q_{30}            = 0.41  \cdot A^{2}    \cdot \beta_3 \f$
!!   - \f$ Q_{20}            = 0.34  \cdot A^{5/3}  \cdot \beta_2 \f$
!> @param[in] amass 
!> REAL(8), takes the mass number.
!--------------------------------------------------------------------------- 
SUBROUTINE prep_qsoll(amass)

USE paramr
USE levels, ONLY: octini,betain

IMPLICIT NONE

REAL(8),INTENT(IN) :: amass

REAL(8) :: reff,reffr,reffz,b2r,b2z,b3max

! the scaling factors
  cfact0  = 1D0                      ! scaling yet to be defined  
  cfact2  = 1.09D0 * amass**1.666666667
  cfact3  = 0.41D0 * amass**2
  cfact4  = 0.495D0 * amass**2.3333333333

! convert constraint switch to logical
! note that constraint is deactivated with 'betaw..==0'
  tswicon0 = ABS(swicon0)>0.5D0 .AND. betaw0.NE.0
  tswicon2 = ABS(swicon2)>0.5D0 .AND. betaw2.NE.0
  tswicon3 = ABS(swicon3)>0.5D0 .AND. betaw3.NE.0
  tswicon4 = ABS(swicon4)>0.5D0 .AND. betaw4.NE.0
  IF(tswicon0) tfixcon0=swicon0<-0.5D0
  IF(tswicon2) tfixcon2=swicon2<-0.5D0
  IF(tswicon3) tfixcon3=swicon3<-0.5D0
  IF(tswicon4) tfixcon4=swicon4<-0.5D0
  IF(tswicon0 .AND. .NOT.tfixcon0) STOP " L=0 allows only fixed constraining force"

! define limits in terms of multipole moments
  IF(tswicon0) THEN
    qmin0   = betmin0*cfact0
    qmax0   = betmax0*cfact0
    qw0     = betaw0 *cfact0
    nqmax0  = NINT((qmax0-qmin0)/qw0)
  END IF
  IF(tswicon2) THEN
    qmin2   = betmin2*cfact2
    qmax2   = betmax2*cfact2
    qw2     = betaw2 *cfact2
    nqmax2  = NINT((qmax2-qmin2)/qw2)
  END IF
  IF(tswicon3) THEN
    qmin3   = betmin3*cfact3
    qmax3   = betmax3*cfact3
    qw3     = betaw3 *cfact3
    nqmax3  = NINT((qmax3-qmin3)/qw3)
  END IF
  IF(tswicon4) THEN
    qmin4   = betmin4*cfact4
    qmax4   = betmax4*cfact4
    qw4     = betaw4 *cfact4
    nqmax4  = NINT((qmax4-qmin4)/qw4)
  END IF
  
! automatic setting of grid size to cover all conceivable deformations
! also for the case of no constraint
  IF(nz == 0) THEN
    IF(nr == 0) THEN
      reff  = 1.2D0*amass**0.3333333D0
      b2z   = MAX(betain,0D0)
      b2r   = MAX(-betain,0D0)
      b3max = ABS(octini)
      IF(tswicon2) THEN
        b2z = MAX(betmin2,betmax2,b2z)
        b2r = MIN(betmin2,betmax2,b2r)
      END IF
      IF(tswicon3) b3max = MAX(ABS(betmin3),ABS(betmax3),b3max)
    ELSE
      STOP ' NZ, NR must be both 0 or both nonzero'
    END IF
    reffz  = reff*(1.0+b2z+0.75*b3max)
    reffr  = reff*(1.0-b2r/3.0+0.65*b3max)
    nz     = INT((2*reffz+6)/dz)
    IF(.NOT.tsym) nz = 2*nz+1
    nr     = NINT((2D0*reffr+6D0)/dr)
    IF(itest>0) WRITE(6,'(a/6f8.2,2i5)') ' b2z,b2r,b3max,reff,reffz,reffr,nz,nr=',  &
        b2z,b2r,b3max,reff,reffz,reffr,nz,nr
  ELSE
    IF(nr == 0) STOP ' NZ, NR must be both 0 or both nonzero'
  END IF

END SUBROUTINE prep_qsoll

END MODULE constr
