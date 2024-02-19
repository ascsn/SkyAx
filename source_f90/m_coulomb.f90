!------------------------------------------------------------------------------
! MODULE: coulomb
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routines and variables to solve the Coulomb problem. 
!> @details
!!The Coulomb field is determined from the charge density by the Poisson 
!!equation
!!\f[\Delta\Phi(\mathbf{r}) = 4\pi e^2\rho(\mathbf{r}) \quad.\f]
!!Although the given density \f$ \rho(\mathbf{r}) \f$ is well localized, the
!!resulting Coulomb field \f$\Phi(\mathbf{r})\f$ is long ranged. The behavior at 
!!large distance is determined by the multipole moments
!!\f[ q_{lm}=\int d^3r\,r^lY_{lm}\rho(\mathbf{r})\quad.\f]
!!We deal with axially symmetric systems which have only non-zero moments for 
!!\f$m=0\f$ and ignore \f$m\neq 0\f$ in the following.
!! 
!!To cope with the long range behavior, we separate the problem into a
!!longe-range part and a short-range part.
!!\f[\rho=\rho^\mathrm{(long)}+\rho^\mathrm{(short)}\quad,\f]
!!\f[ \rho^\mathrm{(long)}(\mathbf{r})=\sum_{l=0}^4\sum_{m=-l}^lq_{lm}
!!\rho_l(r)Y_{lm}(\Omega_\mathbf{r})\quad,\f]
!!where the long-range reference densities are defined as
!!\f[\rho_l(r)=\mathcal{N}_l\frac{r^{l+2}}{R_\mathrm{max}+r^l}
!!e^{-r^2/\mu_l^2}\quad,\f]
!!\f[\mathcal{N}_l: 1=4\pi\int_0^\infty dr\,r^{l+2}\,\rho_l(r)\quad.\f]
!!The choice of the width parameter is \f$\mu_l=R_\mathrm{min}/(6\log2)\f$  
!!The reference radii are \f$R_\mathrm{min}\f$ the smaller of the extensions 
!!in \f$z\f$ and \f$r\f$ direction and \f$R_\mathrm{max}\f$ the larger one.
!!
!!The Coulomb problem for the reference densities becomes a simple radial 
!!Poisson equation which can be solved in inexpensive manner by a couple 
!!of radial integrals. The solution for the long-range part is thus composed as
!!\f[\Phi=\Phi^\mathrm{(long)}+\Phi^\mathrm{(short)}\quad,\f]
!!\f[\Phi^\mathrm{(long)}(\mathbf{r})=\sum_{l=0}^4q_{l0}
!!\Phi_l(r)Y_{l0}(\Omega_\mathbf{r})\quad,\f]
!!\f[\Phi_l(r)=\frac{e^2}{r^{l+1}}\int_0^r dr'\,{r'}^l\,\rho_l(r')+
!!e^2{r}^l\int_r^\infty dr'\,\frac{1}{{r'}^{l-1}}\,\rho_l(r')\quad.\f]
!!
!!There remains the Poisson equation for the short-range part which is 
!!simplified by the boundary conditions \f$\Phi^\mathrm{(short)}
!!\big|_\mathrm{bounds}=0\f$, i.e. that the field vanishes everywhere at the 
!!bounds of the box. We solve it by damped gradient iteration similarly as for 
!!the s.p. wavefunctions. 
!!
!------------------------------------------------------------------------------
MODULE coulomb
USE paramr

IMPLICIT NONE

INTEGER,PUBLIC :: nqrun=1          !<Switch for Coulomb quadrupole field.
INTEGER,PUBLIC :: nhrun=1          !<Switch for Coulomb hexadecupole field.
INTEGER,PUBLIC :: ndrun=1          !<Switch for Coulomb dipole field.
INTEGER,PUBLIC :: norun=1          !<Switch for Coulomb octupole field.
INTEGER,PUBLIC :: irepc=50         !<Max number of short-range Coulomb iterations.

REAL(8) :: colerr=1D-7             !<Allowed Coulomb error. 
REAL(8) :: cmssqr=0.6D0            !<Step width for Coulomb iterations.
REAL(8) :: dmsfac=1D0              !<Factor for meson damping.

REAL(8) :: cmono                   !<Coulomb monopole moment.
REAL(8) :: cquad                   !<Coulomb quadrupole moment.
REAL(8) :: chex                    !<Coulomb hexadecupole moment.
REAL(8) :: cdip                    !<Coulomb dipole moment.
REAL(8) :: coct                    !<Coulomb octupole moment.

REAL(8),ALLOCATABLE,TARGET,PRIVATE :: colref(:,:,:)  !<Super-array 
!!which contains all reference fields:
!!l=0 density: lpoint = 0; l=0 potential: lpoint = 1; l=1 density: lpoint = 2...
!!l=4 potential: lpoint = 9.
REAL(8),POINTER,PRIVATE :: crhom(:,:)    !<Coulomb monopole density.
REAL(8),POINTER,PRIVATE :: coulm(:,:)    !<Coulomb monopole potential.
REAL(8),POINTER,PRIVATE :: crhoq(:,:)    !<Coulomb quadrupole density.
REAL(8),POINTER,PRIVATE :: coulq(:,:)    !<Coulomb quadrupole potential.
REAL(8),POINTER,PRIVATE :: crhoh(:,:)    !<Coulomb hexadecupole density.
REAL(8),POINTER,PRIVATE :: coulh(:,:)    !<Coulomb hexadecupole potential.
REAL(8),POINTER,PRIVATE :: crhod(:,:)    !<Coulomb dipole density.
REAL(8),POINTER,PRIVATE :: could(:,:)    !<Coulomb dipole potential.
REAL(8),POINTER,PRIVATE :: crhoo(:,:)    !<Coulomb octupole density.
REAL(8),POINTER,PRIVATE :: coulo(:,:)    !<Coulomb octupole potential.

REAL(8),ALLOCATABLE,PUBLIC  :: wcolt(:,:)  !<Coulomb potential.
REAL(8),ALLOCATABLE,PRIVATE :: sqrm(:,:)   !<Damping matrix for Coulomb

! line 1 for all symmetries (sym=0 and 1)
REAL(8),PRIVATE :: am2            !<factor on long-range monopole solution
REAL(8),PRIVATE :: aq2            !<factor on long-range quadrupole solution
REAL(8),PRIVATE :: ah2            !<factor on long-range hexadecapole solution
REAL(8),PRIVATE :: rad0           !<reference radius in long-range Coulomb
REAL(8),PRIVATE :: rad02          !< rad0**2
REAL(8),PRIVATE :: rad03          !< rad0**3
REAL(8),PRIVATE :: rad04          !< rad0**4
REAL(8),PRIVATE :: radpow         !<Power of prefactor for integration, transfered to function call.
! 2 for sym=0
REAL(8),PRIVATE :: ad2            !<factor on long-range dipole solution
REAL(8),PRIVATE :: ao2            !<factor on long-range octupole solution

REAL(8),PRIVATE :: aexp(0:4)      !<Width parameter for the Gaussian.
REAL(8),PRIVATE :: bdiv(0:4)      !<Weight parameter for the denominator.

INTEGER,PRIVATE :: lpow           !<Power of actual 'l', transfered to function call.

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: coulmb
!> @brief
!!Main routine for solving the Coulomb routine.
!> @details
!!Calls initialization routines if \c iter = 0 and calls routine to solve 
!!the Coulomb problem.
!>
!> @param[in] iter
!> INTEGER, takes the current iteration number.
!--------------------------------------------------------------------------- 
SUBROUTINE coulmb(iter)

USE fields
USE wavefr
USE levels
USE grid
USE observ, ONLY : iprint
IMPLICIT NONE

INTEGER, INTENT(IN) :: iter

IF (iter == 0) THEN
! define pointers
  crhom(0:,0:) => colref(0:nz,0:nr,0)
  coulm(0:,0:) => colref(0:nz,0:nr,1)
  crhoq(0:,0:) => colref(0:nz,0:nr,2)
  coulq(0:,0:) => colref(0:nz,0:nr,3)
  crhoh(0:,0:) => colref(0:nz,0:nr,4)
  coulh(0:,0:) => colref(0:nz,0:nr,5)
  crhod(0:,0:) => colref(0:nz,0:nr,6)
  could(0:,0:) => colref(0:nz,0:nr,7)
  crhoo(0:,0:) => colref(0:nz,0:nr,8)
  coulo(0:,0:) => colref(0:nz,0:nr,9)

  CALL colin2(coulpt(0:nz,0:nr))
END IF
CALL colrun(crho(0:nz,0:nr),coulpt(0:nz,0:nr),iter)
IF(itest > 1 .AND. MOD(iter,iprint) == 0)  &
    CALL test_print(coulpt(0:nz,0:nr),0,15,0,9,' coulpt')

RETURN
END SUBROUTINE coulmb

!---------------------------------------------------------------------------  
! DESCRIPTION: inicms
!> @brief
!!Initialises the damping mass for the Coulomb potential.
!> @details
!!This variant is oriented towards the boundaries of the box and carries the 
!!damping \c d2 = 2 fm inside.
!--------------------------------------------------------------------------- 
SUBROUTINE inicms

USE paramr
USE constr
USE dynpar
USE fields
USE levels
USE grid
IMPLICIT NONE

INTEGER :: ir, iz
REAL(8) :: drmin,drmin2,drr,dzadd,dradd,dtermz,dtermr,rmax,zmax

REAL(8),PARAMETER :: ddw=0.4D0        ! width of damping
REAL(8),PARAMETER :: ddmass=1D0       ! factor on damping mass
REAL(8),PARAMETER :: d2=2D0           ! thickness of surface zone
REAL(8)           :: dqr=0D0          ! to make this variable static

IF(.NOT.ALLOCATED(colref)) THEN
  ALLOCATE(colref(0:nz,0:nr,0:9),wcolt(0:nz,0:nr),sqrm(0:nz,0:nr))
  colref=0D0
  wcolt=0D0
  sqrm=0D0
END IF

zmax    = ((nz-icm)+sz)*dz
rmax    = ((nr)+sr)*dr
drmin = MIN(zmax,rmax) - d2
drmin2 = drmin * drmin
dzadd = zmax - drmin - d2
dradd = rmax - drmin - d2
dtermr = drmin2 - dradd*dradd
dtermz = drmin2 - dzadd*dzadd

DO iz = 0, nz
  IF( dzadd >= dradd ) THEN
    IF (z(iz) >= dzadd) THEN
      dqr = SQRT(dtermz + d2*dzadd*z(iz))
    ELSE
      dqr = SQRT(z(iz)**2 + drmin2)
    END IF
  END IF
  DO ir = 0, nr
    drr = DSQRT(z(iz)**2 + r(ir)**2)
    IF( dradd > dzadd ) THEN
      IF (r(ir) >= dradd) THEN
        dqr = SQRT(dtermr + d2*dradd*r(ir))
      ELSE
        dqr = SQRT(r(ir)**2 + drmin2)
      END IF
    END IF
    sqrm(iz,ir) = ddmass / (one + EXP((dqr-drr)/ddw))
  END DO
END DO

IF(itest > 0) CALL test_print(sqrm(0:nz,0:nr),0,nz,0,nr,'sqrm')

RETURN
END SUBROUTINE inicms

!---------------------------------------------------------------------------  
! DESCRIPTION: colrun
!> @brief
!!Solves the Poisson equation for charged particles
!> @details
!!Uses the FALR technique, see \cite Lauritsch94 .
!!An exactly solvable long-range part of the Coulomb problem is
!!separated and the remaining short-range part, which is then free
!!of singularities in k-space, is solved by FFT.
!>
!> @param[in] rhoin
!> REAL(8), array, takes the proton density.
!> @param[out] wcol
!> REAL(8), array, returns the Coulomb potential.
!> @param[in] iter
!> INTEGER, takes the current iteration number.
!--------------------------------------------------------------------------- 
SUBROUTINE colrun(rhoin,wcol,iter)

USE paramr
USE fields
USE wavefr
USE levels
USE grid
USE dynpar, ONLY: chg
USE observ, ONLY: obs
IMPLICIT NONE

REAL(8), INTENT(IN)                  :: rhoin(0:nz,0:nr)
REAL(8), INTENT(OUT)                 :: wcol(0:nz,0:nr)
INTEGER, INTENT(IN)                  :: iter

REAL(8),ALLOCATABLE ::  rhoa(:,:)
REAL(8),SAVE        :: ephimx=0D0     ! to make it static

INTEGER :: irepmx,irep

REAL(8) :: alpha,atnump,dfac,dp,facd,fach,facm,faco,facq,hp,op,qp,qsqp,rmsp

ALLOCATE(rhoa(0:nz,0:nr))
rhoa=0D0

IF (iter <= 10) THEN
  irepmx = 10*irepc
ELSE
  irepmx = irepc
END IF

alpha   = 4.0*pi*chg
dfac    = one/cmssqr

!--------------------------------------------------------------------------- 
!Step 1: Calculates the moments $l\leq 4$ of proton density and determines 
!        parameters to project out those moments with the densities crhom, 
!        crhod, crhoq and crhoh
!--------------------------------------------------------------------------- 
CALL obs(rhoin(0:nz,0:nr),atnump,rmsp,dp,qp,op,hp,qsqp,zero)

IF(itest >= 1) WRITE(6,'(a,/4g16.6)') ' colrun with dp,qp,op,hp=',dp,qp,op,hp
IF(nhrun == 1) fach   = hp/chex
IF(nqrun == 1) facq   = qp/cquad

IF(.NOT.tsym .AND.ndrun == 1) THEN
  facd   = dp/cdip
ELSE
  facd   = zero
END IF

IF(.NOT.tsym .AND.nhrun == 1) THEN
  faco   = op/coct
ELSE
  faco   = zero
END IF

facm   = atnump/cmono

!--------------------------------------------------------------------------- 
!Step 2: Project out low l long-range part. 
!--------------------------------------------------------------------------- 

IF (iter == 0) THEN
  rhoa(0:nz,0:nr)  = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr))
ELSE IF (tsym) THEN
  IF (nqrun == 0 .AND. nhrun == 0) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE IF (nqrun == 1 .AND. nhrun == 0) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr)  &
                    -facq*crhoq(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE IF (nqrun == 1 .AND. nhrun == 1) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr)  &
                    -facq*crhoq(0:nz,0:nr) -fach*crhoh(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE
    STOP ' no such case for "nqrun" or "nhrun"'
  END IF
ELSE IF (.NOT.tsym) THEN
  IF (ndrun == 0) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE IF (ndrun == 1 .AND. nqrun == 0) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr)  &
                    -facd*crhod(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE IF (nqrun == 1 .AND. norun == 0) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr)  &
                    -facd*crhod(0:nz,0:nr) -facq*crhoq(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE IF (norun == 1 .AND. nhrun == 0) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) -facm*crhom(0:nz,0:nr)  &
                    - facd*crhod(0:nz,0:nr) - facq*crhoq(0:nz,0:nr)  &
                    - faco*crhoo(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  ELSE IF (nhrun == 1) THEN
    rhoa(0:nz,0:nr) = alpha*(rhoin(0:nz,0:nr) - facm*crhom(0:nz,0:nr)  &
                    - facd*crhod(0:nz,0:nr) - facq*crhoq(0:nz,0:nr)  &
                    - faco*crhoo(0:nz,0:nr) - fach*crhoh(0:nz,0:nr))
    wcol(0:nz,0:nr) = wcolt(0:nz,0:nr)
  END IF
END IF

!--------------------------------------------------------------------------- 
!Step 3: Solve remaining short-range potential.
!--------------------------------------------------------------------------- 

DO  irep=1,irepmx
  
  CALL laplace_scalar(wcol(0:nz,0:nr),wcolt(0:nz,0:nr),+1,0,1,-one)
  wcolt(0:nz,0:nr) = wcolt(0:nz,0:nr) + sqrm(0:nz,0:nr)*wcol(0:nz,0:nr) - rhoa(0:nz,0:nr)

! 
! damped gradient step for iterating short-range Coulomb field
!
  CALL damp_meson(wcolt(0:nz,0:nr),wcolt(0:nz,0:nr),one)
  
  ephimx = MAXVAL(ABS(wcolt(:,:)))
  wcol(:,:) = wcol(:,:) - dfac*wcolt(:,:)
  IF (ABS(ephimx) <= colerr) EXIT
END DO


!--------------------------------------------------------------------------- 
!Step 4: Add long-range part of potential with solutions coulm, could, coulq, 
!        coulh, already determined at initialization.
!--------------------------------------------------------------------------- 
IF (tsym) THEN
  IF (nqrun == 0.OR.iter == 0) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol (0:nz,0:nr)= wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)
  ELSE IF (nqrun == 1.AND.nhrun == 0) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol(0:nz,0:nr) = wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)  &
                + facq*coulq(0:nz,0:nr)
  ELSE IF (nqrun == 1.AND.nhrun == 1) THEN
        wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
        wcol(0:nz,0:nr) = wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)  &
            + facq*coulq(0:nz,0:nr) + fach*coulh(0:nz,0:nr)
  END IF
ELSE IF (.NOT.tsym) THEN
  IF (ndrun == 0.OR.iter == 0) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol (0:nz,0:nr)= wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)
  ELSE IF (ndrun == 1.AND.nqrun == 0) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol (0:nz,0:nr)= wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)  &
            + facd*could(0:nz,0:nr)
  ELSE IF (nqrun == 1.AND.norun == 0) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol (0:nz,0:nr)= wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)  &
                    + facd*could(0:nz,0:nr) + facq*coulq(0:nz,0:nr)
  ELSE IF (norun == 1.AND.nhrun == 0) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol (0:nz,0:nr)= wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)  &
                    + facd*could(0:nz,0:nr) + facq*coulq(0:nz,0:nr)&
                    + faco*coulo(0:nz,0:nr)
  ELSE IF (nhrun == 1) THEN
    wcolt(0:nz,0:nr)= wcol(0:nz,0:nr)
    wcol (0:nz,0:nr)= wcol(0:nz,0:nr) + facm*coulm(0:nz,0:nr)  &
                    + facd*could(0:nz,0:nr) + facq*coulq(0:nz,0:nr)&
                    + faco*coulo(0:nz,0:nr) + fach*coulh(0:nz,0:nr)
  END IF
  IF (itest >= 1) THEN
    WRITE(6,'(a,5g15.5)') 'fm,fd,fq,fo,fh =',facm,facd,facq,faco, fach
  END IF
END IF

IF (ABS(ephimx) >= colerr) THEN
  WRITE(6,'(/2(a,g15.5),i3,a/)') ' colrun: max error=',  &
      colerr,' actual error=',ABS(ephimx), irepmx,' iterations'
END IF
IF(itest >= 1) CALL coldif(wcol(0:nz,0:nr),rhoin(0:nz,0:nr),wcolt(0:nz,0:nr),' full solution')

DEALLOCATE(rhoa)


RETURN
END SUBROUTINE colrun

!---------------------------------------------------------------------------  
! DESCRIPTION: coldif
!> @brief
!!Tests error in Poisson equation.
!>
!> @param[in] rhoin
!> REAL(8), array, takes the proton density.
!> @param[in] wcol
!> REAL(8), array, takes the Coulomb potential.
!> @param[out] work
!> REAL(8), array, returns the array with errors.
!> @param[in] text
!> CHARACTER, takes a text for print out.
!--------------------------------------------------------------------------- 
SUBROUTINE coldif(colin,rhoin,work,text)

USE paramr
USE fields
USE levels
USE grid
USE dynpar, ONLY: chg
IMPLICIT NONE

REAL(8),           INTENT(IN)        :: colin(0:nz,0:nr)
REAL(8),           INTENT(IN)        :: rhoin(0:nz,0:nr)
REAL(8),           INTENT(OUT)       :: work(0:nz,0:nr)
CHARACTER (LEN=*), INTENT(IN)        :: text

INTEGER :: izmax,iz,irmax,ir

REAL(8) :: alpha, differ, errav, errmax, errwgt, rhowgt

WRITE(6,'(2a)') ' test poisson equation for field:',text

CALL laplace_scalar(colin(0:nz,0:nr),work(0:nz,0:nr),+1,0,1,-one)
alpha   = 4.0*pi*chg
errav   = zero
errwgt  = zero
errmax  = zero
irmax   = 0
izmax   = 0
DO ir=0,nr
  DO iz=0,nz
    rhowgt = rhoin(iz,ir)
    differ  = ABS(work(iz,ir)-alpha*rhowgt)
    work(iz,ir) = differ
    IF(differ > errmax) THEN
      irmax  = ir
      izmax  = iz
      errmax = differ
    END IF
    errwgt = rhowgt + errwgt
    errav  = differ*differ*rhowgt+errav
  END DO
END DO
errav  = errav/errwgt

WRITE(6,'(a,g12.4,a,g12.4,a,2i5)') ' average error=',SQRT(ABS(errav)),  &
    ' max. error=',errmax,' at ir,iz=',irmax,izmax
IF(itest >= 2) THEN
  WRITE(6,'(a,20(/1x,10(1pg11.3)))')  &
      ' error along r-axis:',(work(0,ir),ir=0,nr)
  WRITE(6,'(a,20(/1x,10(1pg11.3)))')  &
      ' error along z-axis:',(work(iz,0),iz=0,nz)
END IF

RETURN
END SUBROUTINE coldif

!---------------------------------------------------------------------------  
! DESCRIPTION: colin2
!> @brief
!!Initialises the monpole and quadrupole background fields for the 
!!subroutine colrun, here using solution by radial integration of the 
!!Green's function.
!> @details
!!Prepares the exactly solvable part (densities and corresponding
!!Coulomb potentials) for the FALR technique, see \cite Lauritsch94 .
!!
!> @param[in] rhoin
!> REAL(8), array, takes the proton density.
!> @param[out] wcol
!> REAL(8), array, returns the Coulomb potential and sets it to zero.
!--------------------------------------------------------------------------- 
SUBROUTINE colin2(wcol)

USE paramr
USE fields
USE levels
USE grid, ONLY: nr,nz,dr,dz,icm,tsym,sr,sz
USE observ, ONLY: obs
USE dynpar, ONLY: chg
IMPLICIT NONE

REAL(8), INTENT(OUT)      :: wcol(0:nz,0:nr)

INTEGER :: ir,iz

REAL(8) :: ad,ah,am,am3,ao,aq,bh,bin1,binup,bm,bq,cdhex,cdmono,&
           cdoct,cdquad,chdip,chmono,choct,chquad,cmdip,cmhex,cmoct,cmquad,&
           codip,cohex,comono,coquad,cqdip,cqhex,cqmono,cqoct,dh,dh0,&
           dinteg,dm0,do,dp,dq0,fd0,fh0,fm0,fo0,fq0,fzrm,hint,hp,oint,&
           op,pfy,pifac,pro,qinteg,qp,qqt2,qsqcd,qsqch,qsqcm,qsqco,&
           qsqcq,qt1,qt2,r,rad,rad2,rad3,rad4,rad5,rad7,rad9,rmax,rmsd,&
           rmsh,rmsm,rmso,rmsq,y10p,y20p,y30p,y40p,z,zmax,zrmax,zrmin

wcol(0:nz,0:nr) = zero

!     general parameters for background fields
zmax    = ((nz-icm)+sz)*dz
rmax    = ((nr)+sr)*dr
zrmin   = MIN(zmax,rmax)
zrmax   = MAX(zmax,rmax)
binup   = zrmax+5.0D0
fzrm    = zrmin/5.5D0
rad0    = zrmax
rad02   = rad0 *rad0

!     monopole field and density
pro     = one
am      = 6.0*SQRT(-LOG(oh))/zrmin
am2     = am*am
am3     = am2*am
bm      = two*am3*pro/pi
pifac   = two*SQRT(pi)
dm0     = bm/pifac
fm0     = bm*pifac*chg
DO ir=0,nr
  r = (ir+sr)*dr
  DO iz=0,nz
    z     = (iz+sz-icm)*dz
    rad2  = z*z+r*r
    rad   = SQRT(rad2)
    radpow = zero                              ! setting for 'fm1'
    crhom(iz,ir) = dm0*fm1(rad)
    
    bin1 = am * rad
    IF (bin1 < 1.d-2) THEN
      qt1 = zero
    ELSE
      radpow = two                               ! setting for 'fm1'
      CALL qgaus (fm1, zero, rad, qt1)
      qt1 = qt1/(rad)
    END IF
    radpow = one
    CALL qgaus(fm1, rad , binup, qt2)
    coulm(iz,ir) = fm0*(qt1+qt2)
  END DO
END DO
CALL obs(crhom(0:nz,0:nr),cmono,rmsm,cmdip,cmquad,cmoct,cmhex,qsqcm,zero)
IF(itest>0)WRITE(6,'(a/3x,7(1pg11.3))')  &
    ' from "crhom": cmono,rmsm,cmdip,cmquad,cmoct,cmhex,qsqcm=',  &
    cmono,rmsm,cmdip,cmquad,cmoct,cmhex,qsqcm
IF(itest >= 1) THEN
  CALL coldif(coulm(0:nz,0:nr),crhom(0:nz,0:nr),wcolt(0:nz,0:nr),' monopole background')
END IF

IF(nqrun == 1) THEN
  qp      = one
  aq      = SQRT( -LOG( oh ) )/fzrm
  aq2     = aq*aq
  radpow = 6.0
  CALL qgaus(fq1, zero, binup, qinteg)
  pifac   = 0.25 * SQRT(5.0/pi)
  bq      = qp / qinteg * pifac
  dq0     = bq
  fq0     = bq*0.8D0*pi*chg
  pfy     = pifac
  DO ir=0,nr
    r = (ir+sr)*dr
    DO iz=0,nz
      z     = (iz+sz-icm)*dz
      rad2  = z*z+r*r
      rad   = SQRT(rad2)
      rad5  = rad2*rad2*rad
      y20p  = pfy*(two*z*z-r*r)
      radpow = zero
      crhoq(iz,ir) = dq0*fq1(rad)*y20p
      
      IF (rad < 1.d-2) THEN
        qt1 = zero
      ELSE
        radpow = 6.0
        CALL qgaus(fq1, zero, rad, qt1)
        qt1 = qt1/rad5
      END IF
      radpow = one
      CALL qgaus(fq1, rad , binup, qqt2)
      coulq(iz,ir) = fq0*(qt1+qqt2)*y20p
    END DO
  END DO
  CALL obs(crhoq(0:nz,0:nr),cqmono,rmsq,cqdip,cquad,cqoct,cqhex,qsqcq, zero)
  IF(itest>0) WRITE(6,'(a/3x,7(1pg11.3))')  &
      ' from "crhoq": cqmono,rmsq,cqdip,cquad,cqoct,cqhex,qsqcq=',  &
      cqmono,rmsq,cqdip,cquad,cqoct,cqhex,qsqcq
  IF(itest >= 1) CALL coldif(coulq(0:nz,0:nr),crhoq(0:nz,0:nr),wcolt(0:nz,0:nr),' quadrupole background')
ELSE
  coulq(0:nz,0:nr) = 0.0d0
  crhoq(0:nz,0:nr) = 0.0d0
END IF

!     optionally initialize the hexadecupole fields
IF(nhrun == 1) THEN
  hp      = one
  rad04   = rad02*rad02
  ah      = SQRT( -LOG( oh ) )/fzrm
  ah2     = ah*ah
  radpow = 10.0
  CALL qgaus(fh1, zero, binup, hint)
  bh      = hp / hint
  dh0     = bh
  fh0     = bh*4.0D0/9.0D0*pi*chg
  pfy     = 3.0 / 16.0 / SQRT(pi)
  DO ir=0,nr
    r = (ir+sr)*dr
    DO iz=0,nz
      z     = (iz+sz-icm)*dz
      rad2  = z*z+r*r
      rad4  = rad2*rad2
      rad   = SQRT(rad2)
      rad5  = rad4*rad
      rad9  = rad4*rad5
      y40p  = pfy*(8.0D0*z*z*z*z + 3.0D0 *r*r*r*r - 24.0D0 *z*z*r*r)
      radpow = zero
      crhoh(iz,ir) = dh0*fh1(rad)*y40p
      
      IF (rad < 1.d-2) THEN
        qt1 = zero
      ELSE
        radpow = 10.0
        CALL qgaus(fh1, zero, rad, qt1)
        qt1 = qt1/rad9
      END IF
      radpow = one
      CALL qgaus(fh1, rad , binup, qqt2)
      coulh(iz,ir) = fh0*(qt1+qqt2)*y40p
    END DO
  END DO
  CALL obs(crhoh(0:nz,0:nr),chmono,rmsh,chdip,chquad,choct,chex,qsqch, zero)
  IF(itest>0)WRITE(6,'(a/3x,7(1pg11.3))')  &
      ' from "crhoh": chmono,rmsh,chdip,chquad,choct,chex,qsqch=',  &
      chmono,rmsh,chdip,chquad,choct,chex,qsqch
  IF(itest >= 1) CALL coldif(coulh(0:nz,0:nr),crhoh(0:nz,0:nr),wcolt(0:nz,0:nr),' hexadecupole background')
ELSE
  coulh(0:nz,0:nr) = 0.0d0
  crhoh(0:nz,0:nr) = 0.0d0
END IF

! initialize dipole and octupole fields in case of tsym being false
IF (.NOT.tsym) THEN

! optionally initialize the dipole fields
  IF (ndrun == 1) THEN
    dp      = one
    ad      = 1.35*SQRT( -LOG( oh ) )/fzrm
    ad2     = ad * ad
    radpow  = 4.0
    CALL qgaus(fd1, zero, binup, dinteg)
    dh      = dp / dinteg
    fd0     = dh*4.0D0/3.0D0*pi*chg
    pfy     = oh * SQRT(3./pi)
    DO ir=0,nr
      r = (ir+sr)*dr
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
        rad2  = z*z+r*r
        rad   = SQRT(rad2)
        rad3  = rad2 * rad
        y10p  = pfy * z
        radpow = zero
        crhod(iz,ir) = dh*fd1(rad)*y10p
!   mind periodic boundaries:
        IF (iz == nz) crhod(iz,ir) = zero
        IF (iz == 0.AND..NOT.tsym) crhod(iz,ir) = zero
       
        IF (rad < 1.d-2) THEN
          qt1 = zero
        ELSE
          radpow = 4.0
          CALL qgaus(fd1, zero, rad, qt1)
          qt1 = qt1/rad3
        END IF
        radpow = one
        CALL qgaus(fd1, rad , binup, qt2)
        could(iz,ir) = fd0*(qt1+qt2)*y10p
!   mind periodic boundaries:
        IF (iz == nz) could(iz,ir) = zero
        IF (iz == 0.AND..NOT.tsym) could(iz,ir) = zero
      END DO
    END DO
    CALL obs(crhod(0:nz,0:nr),cdmono,rmsd,cdip,cdquad,cdoct,cdhex,qsqcd, zero)
    IF(itest>0) WRITE(6,'(a/3x,7(1pg11.3))')  &
        ' from "crhod": cdmono,rmsd,cdip,cdquad,cdoct,cdhex,qsqcd=',  &
        cdmono,rmsd,cdip,cdquad,cdoct,cdhex,qsqcd
    IF(itest >= 1) THEN
      CALL coldif(could(0:nz,0:nr),crhod(0:nz,0:nr),wcolt(0:nz,0:nr),' dipole background')
    END IF
  ELSE
    could(0:nz,0:nr) = 0.0d0
    crhod(0:nz,0:nr) = 0.0d0
  END IF

! optionally initialize the octupole fields
  IF (norun == 1) THEN
    rad03   = rad0 * rad02
    op      = one
    ao      = SQRT( -LOG( oh ) )/fzrm
    ao2     = ao * ao
    radpow  = 8.0
    CALL qgaus(fo1, zero, binup, oint)
    DO      = op / oint
      fo0     = DO * 4.0D0/7.0D0 *pi*chg
        pfy     = 0.25D0 * SQRT(7.0D0/pi)
        DO ir=0,nr
          r = (ir+sr)*dr
          DO iz=0,nz
            z     = (iz+sz-icm)*dz
            rad2  = z*z+r*r
            rad   = SQRT(rad2)
            rad3  = rad2 * rad
            rad7  = rad3 * rad3 * rad
            y30p  = pfy * (two*z*z*z - 3.*z*r*r)
            radpow = zero
            crhoo(iz,ir) = DO*fo1(rad)*y30p
!   mind periodic boundaries:
            IF (iz == nz) crhoo(iz,ir) = zero
            IF (iz == 0.AND..NOT.tsym) crhoo(iz,ir) = zero
            
            IF (rad < 1.d-2) THEN
              qt1 = zero
            ELSE
              radpow = 8.0
              CALL qgaus(fo1, zero, rad, qt1)
              qt1 = qt1/rad7
            END IF
            radpow = one
            CALL qgaus(fo1, rad , binup, qt2)
            coulo(iz,ir) = fo0*(qt1+qt2)*y30p
!   mind periodic boundaries:
            IF (iz == nz) coulo(iz,ir) = zero
            IF (iz == 0.AND..NOT.tsym) coulo(iz,ir) = zero
          END DO
        END DO
        CALL obs(crhoo(0:nz,0:nr),comono,rmso,codip,coquad,coct,cohex,qsqco ,zero)
        IF(itest>0) WRITE(6,'(a/3x,7(1pg11.3))')  &
            ' from "crhoo": comono,rmso,codip,coquad,coct,cohex,qsqco=',  &
            comono,rmso,codip,coquad,coct,cohex,qsqco
        IF(itest >= 1) THEN
          CALL coldif(coulo(0:nz,0:nr),crhoo(0:nz,0:nr),wcolt(0:nz,0:nr),' octupole background')
        END IF
      ELSE
        coulo(0:nz,0:nr) = 0.0d0
        crhoo(0:nz,0:nr) = 0.0d0
      END IF
    END IF
    
    RETURN
END SUBROUTINE colin2

!---------------------------------------------------------------------------  
! DESCRIPTION: fm1
!> @brief
!!Helper function for Coulomb monopole.
!>
!> @param[in] x
!> REAL(8), takes the function argument.
!--------------------------------------------------------------------------- 
REAL(8) PURE FUNCTION fm1(x)
USE paramr

IMPLICIT NONE
  
REAL(8), INTENT(IN)   :: x

REAL(8) :: expa

  expa  = x*x*am2
  IF (expa <= eighty) THEN
    expa = EXP(-expa)
  ELSE
    expa = zero
  END IF
  IF(x /= zero) THEN
    fm1 =  x**radpow * expa
  ELSE IF(radpow == zero) THEN
    fm1 = expa
  ELSE
    fm1 = zero
  END IF
  
  RETURN
END FUNCTION fm1

!---------------------------------------------------------------------------  
! DESCRIPTION: fq1
!> @brief
!!Helper function for Coulomb quadrupole.
!>
!> @param[in] x
!> REAL(8), takes the function argument.
!--------------------------------------------------------------------------- 
REAL(8) PURE FUNCTION fq1(x)
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)     :: x

REAL(8) :: expa

expa=x*x*aq2
IF (expa <= eighty) THEN
  expa = EXP(-expa)
ELSE
  expa = zero
END IF
IF(radpow == zero) THEN
  fq1 =  expa/(rad02+x**2)
ELSE IF(x /= zero) THEN
  fq1 =  x**radpow*expa/(rad02+x**2)
ELSE
  fq1 = zero
END IF

RETURN
END FUNCTION fq1

!---------------------------------------------------------------------------  
! DESCRIPTION: fh1
!> @brief
!!Helper function for Coulomb hexadecupole.
!>
!> @param[in] x
!> REAL(8), takes the function argument.
!---------------------------------------------------------------------------
REAL(8) PURE FUNCTION fh1(x)
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)       :: x

REAL(8) :: expa

expa=x*x*ah2
IF (expa <= eighty) THEN
  expa = EXP(-expa)
ELSE
  expa = zero
END IF
IF(radpow == zero) THEN
  fh1   =    expa/(rad04+x**4)
ELSE IF(x /= zero) THEN
  fh1   =    x**radpow*expa/(rad04+x**4)
ELSE
  fh1   = zero
END IF

RETURN
END FUNCTION fh1

!---------------------------------------------------------------------------  
! DESCRIPTION: fd1
!> @brief
!!Helper function for Coulomb dipole.
!>
!> @param[in] x
!> REAL(8), takes the function argument.
!---------------------------------------------------------------------------
REAL(8) PURE FUNCTION fd1(x)
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)       :: x

REAL(8) :: expa

expa  = (x-5.0d0)**2*ad2
IF (expa <= eighty) THEN
  expa = EXP(-expa)
ELSE
  expa = zero
END IF
IF(radpow == zero) THEN
  fd1   =    expa/(rad0+x)
ELSE IF(x /= zero) THEN
  fd1   =    x**radpow*expa/(rad0+x)
ELSE
  fd1   = zero
END IF

RETURN
END FUNCTION fd1

!---------------------------------------------------------------------------  
! DESCRIPTION: fo1
!> @brief
!!Helper function for Coulomb octupole.
!>
!> @param[in] x
!> REAL(8), takes the function argument.
!---------------------------------------------------------------------------
REAL(8) PURE FUNCTION fo1(x)
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)        :: x

REAL(8) :: expa

expa  = x**2*ao2
IF (expa <= eighty) THEN
  expa = EXP(-expa)
ELSE
  expa = zero
END IF
IF(radpow == zero) THEN
  fo1   =    expa/(rad03 + x**3)
ELSE IF(x /= zero) THEN
  fo1   =    x**radpow*expa/(rad03 + x**3)
ELSE
  fo1   = zero
END IF

RETURN
END FUNCTION fo1

!---------------------------------------------------------------------------  
! DESCRIPTION: qgaus
!> @brief
!!Performs a 20 point Gauss-Legendre integration.
!>
!> @param[in] func
!> REAL(8), takes the function.
!> @param[in] a
!> REAL(8), takes the lower integration limit.
!> @param[in] b
!> REAL(8), takes the upper integration limit.
!> @param[out] func
!> REAL(8), returns the result.
!---------------------------------------------------------------------------
SUBROUTINE qgaus(func,a,b,ss)
USE paramr

IMPLICIT NONE

REAL(8), EXTERNAL          :: func
REAL(8), INTENT(IN)        :: a
REAL(8), INTENT(IN)        :: b
REAL(8), INTENT(OUT)       :: ss
!
INTEGER :: j

REAL(8),PARAMETER :: x(10)=&
  (/.076526521133497, .227785851141645, .373706088715419,  &
    .510867001950827, .636053680726515, .746331906460150,  &
    .839116971822218, .912234428251325, .963971927277913, .993128599185094/)
REAL(8),PARAMETER :: w(10)=&
  (/.152753387130725, .149172986472603, .142096109318382,  &
    .131688638449176, .118194531961518, .101930119817240,  &
    .083276741576704, .062672048334109, .040601429800386, .017614007139152/)
    
REAL(8) :: dx,sum,xm,xr

xm = 0.5D0*(b+a)
xr = 0.5D0*(b-a)
sum= 0.0D0
DO j=1,10
  dx = xr*x(j)
  sum= sum+w(j) * (func(xm+dx) + func(xm-dx))
END DO
ss = xr*sum
RETURN
END SUBROUTINE qgaus



!---------------------------------------------------------------------------  
! DESCRIPTION: refrho
!> @brief
!!Reference density at x for given angular momentum \c lpow with 
!!prefactor \c radpow.
!>
!> @param[in] x
!> REAL(8), takes function argument.
!---------------------------------------------------------------------------
REAL(8) PURE FUNCTION refrho(x)

USE paramr
IMPLICIT NONE

REAL(8), INTENT(IN)     :: x

REAL(8) :: arg,expa,x2

!x2    = x-rshif(lpow)                         obsolete?
!x2    = x2*x2
x2    = x*x
arg   = x2*aexp(lpow)
IF (arg <= eighty) THEN
  expa = EXP(-arg)
ELSE
  expa = zero
END IF
IF(lpow > 0) THEN
  IF(radpow == zero) THEN
    refrho =  expa/(bdiv(lpow)+x**lpow)
  ELSE IF(x /= zero) THEN
    refrho =  x**radpow*expa/(bdiv(lpow)+x**lpow)
  ELSE
    refrho = zero
  END IF
ELSE
  expa   = expa*(one+arg+oh*arg*arg*(one+edr*arg))
  IF(x /= zero) THEN
    refrho =  x**radpow * expa
  ELSE IF(radpow == zero) THEN
    refrho = expa
  ELSE
    refrho = zero
  END IF
END IF

RETURN
END FUNCTION refrho

!---------------------------------------------------------------------------  
! DESCRIPTION: ylm
!> @brief
!!Spheriacal Harmonics. \f$ R^l\cdot Y_{lm}(r,z),\quad R=\sqrt{r^2+z^2} \f$
!>
!> @param[in] z
!> REAL(8), takes z-coordinate.
!> @param[in] x
!> REAL(8), takes r-coordinate.
!> @param[in] l
!> INTEGER, takes l quantum number.
!> @param[in] min
!> INTEGER, takes m quantum number.
!---------------------------------------------------------------------------
REAL(8) FUNCTION ylm(z,r,l,min)

USE paramr
IMPLICIT NONE

REAL(8), INTENT(IN)         :: z
REAL(8), INTENT(IN)         :: r
INTEGER, INTENT(IN)         :: l
INTEGER, INTENT(IN)         :: min

REAL(8),PARAMETER :: y00fac=0.282094791D0    !     1/(2*SQRT(PI))
REAL(8),PARAMETER :: y10fac=0.488602512D0    !   -SQRT(3)/(2*SQRT(PI))
REAL(8),PARAMETER :: y11fac=0.345494149D0    !    SQRT(6)/(4*SQRT(PI))
REAL(8),PARAMETER :: y20fac=0.315391565D0    !    SQRT(5)/(4*SQRT(PI))
REAL(8),PARAMETER :: y21fac=0.772548404D0    !   -SQRT(30)/(4*SQRT(PI))
REAL(8),PARAMETER :: y22fac=0.386274200D0    !
REAL(8),PARAMETER :: y30fac=0.746352665D0    !   -2*SQRT(7)/(4*SQRT(PI))
REAL(8),PARAMETER :: y31fac=0.323180184D0    !    SQRT(21)/(8*SQRT(PI))
REAL(8),PARAMETER :: y32fac=1D0    !
REAL(8),PARAMETER :: y33fac=1D0    !
REAL(8),PARAMETER :: y40fac=0.846284375D0    !    24/(16*SQRT(PI))
REAL(8),PARAMETER :: y41fac=1.892349392D0    !   -12*SQRT(5)/(8*SQRT(PI))
REAL(8),PARAMETER :: y42fac=1D0    !
REAL(8),PARAMETER :: y43fac=1D0    !
REAL(8),PARAMETER :: y44fac=1D0    !

INTEGER :: m

REAL(8) phasefac

!----------------------------------------------------------------------

m=min

! Check forbidden case=
IF (l < 0) THEN
  STOP ' l<0 is nonsense in FUNCTION YLM '
END IF

IF (m > l) THEN
  STOP ' m>l is nonsense in FUNCTION YLM '
END IF

IF (l > 4) THEN
  STOP ' l>4 not yet provided in FUNCTION YLM '
END IF

! Determin phasefactor for negative m

IF (m < 0) THEN
  phasefac = (-1)**m
  m        = -m
ELSE
  phasefac = +1
END IF

! Spherical harmonics in separate cases

IF (l == 0) THEN
  ylm = y00fac
  
ELSE IF (l == 1) THEN
  IF (m == 0) THEN
    ylm = y10fac*z
  ELSE
    ylm = y11fac*r
  END IF
  
ELSE IF (l == 2) THEN
  IF (m == 0) THEN
    ylm = y20fac*(2.0D0*(z*z)-(r*r))
  ELSE IF (m == 1) THEN
    ylm = y21fac*z*r
  ELSE
    ylm = y22fac*r*r
  END IF
  
ELSE IF(l == 3) THEN
  IF (m == 0) THEN
    ylm = y30fac*z*((z*z)-1.5D0*(r*r))
  ELSE IF (m == 1) THEN
    ylm = y31fac*r*(4.0D0*(z*z)-(r*r))
  ELSE IF (m == 2) THEN
    ylm = y32fac*(r*r)*z
  ELSE
    ylm = y33fac*(r*r)*r
  END IF
  
ELSE
  IF (m == 0) THEN
    ylm = y40fac*((z*z)*(z*z) -3.0D0*(z*z)*(r*r)+0.375D0*(r*r)*(r*r))
  ELSE IF (m == 1) THEN
    ylm = y41fac*z*r*((z*z)-0.75D0*(r*r))
  ELSE IF (m == 2) THEN
    ylm = y42fac*1.0D0 !?
  ELSE IF (m == 3) THEN
    ylm = y43fac*1.0D0 !?
  ELSE
    ylm = y44fac*1.0D0 !?
  END IF
  
END IF

ylm = ylm * phasefac

RETURN
END FUNCTION ylm


END MODULE coulomb
