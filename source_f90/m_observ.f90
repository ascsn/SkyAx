!------------------------------------------------------------------------------
! MODULE: observ
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the variables for most observables and routines to 
!!calculate them. \c sinfo is the main routine for calculating and printing 
!!observables.
!------------------------------------------------------------------------------
MODULE observ

IMPLICIT NONE

REAL(8) :: tempmn=0D0    !<Minimal temperature.
REAL(8) :: tempmx=0D0    !<Maximal temoerature.
REAL(8) :: tempw=1D0     !<Step width for temperature loop.
REAL(8) :: tempac=1D0    !<Actual temperature for current step.
INTEGER :: itmax=1000    !<Maximum number of iterations for one deformation.
INTEGER :: iprint=1000   !<Modulus for printing output.
!NEW ODD
INTEGER, PARAMETER :: mod_upd_conjstates=1000  !< Modulus for updating table of conjugated states. 
!END NEW ODD

REAL(8) :: atnum    !<Mass number for nucleus.
REAL(8) :: apnum    !<Proton number.
REAL(8) :: annum    !<Neutron number.
REAL(8) :: rmsc     !<Charge rms radius.
REAL(8) :: rmsc0    !<Charge rms radius. no l*s folding
REAL(8) :: rms      !<Total rms radius.
REAL(8) :: rmsp     !<Proton rms radius.
REAL(8) :: rmsn     !<Neutron rms radius.
REAL(8) :: r4c      !<charge radius from 4. radial moment
REAL(8) :: r4c0     !<charge radius from 4. radial moment, no l*s folding
REAL(8) :: r4       !<total radius from 4. radial moment
REAL(8) :: r4p      !<proton radius from 4. radial moment
REAL(8) :: r4n      !<neutron radius from 4. radial moment
REAL(8) :: d        !<Dipole moment.
REAL(8) :: dn       !<Neutron dipole moment
REAL(8) :: dp       !<Proton dipole moment
REAL(8) :: q        !<Quadrupole moment.
REAL(8) :: qp       !<Proton quadrupole moment.
REAL(8) :: qn       !<Neutron quadrupole moment.
REAL(8) :: o        !<Octupole moment.
REAL(8) :: on       !<Neutron octupole moment.
REAL(8) :: op       !<Proton octupole moment.
REAL(8) :: h        !<Hexadecupole moment.
REAL(8) :: hp       !<Proton hexadecupole moment.
REAL(8) :: hn       !<Neutron hexadecupole moment.
REAL(8) :: beta1    !<Dipole deformation paramter.
REAL(8) :: beta1p   !<Proton dipole deforamtion parameter.
REAL(8) :: beta1n   !<Neutron dipole deformation parameter.
REAL(8) :: beta2    !<Quadrupole deformation parameter.
REAL(8) :: beta2p   !<Proton quadrupole deformation parameter.
REAL(8) :: beta2n   !<Neutron quadrupole deformation parameter.
REAL(8) :: beta3=0D0    !<Octupole deformation parameter.
REAL(8) :: beta3p   !<Proton octupole deformation parameter.
REAL(8) :: beta3n   !<Neutron octupole deformation parameter.
REAL(8) :: beta4    !<Hexadecupole deformation parameter.
REAL(8) :: beta4p   !<Proton hexadecupole deformation parameter.
REAL(8) :: beta4n   !<Neutron hexadecupole deformation parameter.
REAL(8) :: betafull !<Quadrupole including variance \f$\beta_\mathrm{f}=\sqrt{\beta_{g.s.}^2+\Delta^2\beta}\f$.
REAL(8) :: rsepfrag !<distance of two separated fragments.
REAL(8) :: zcm      !<Center of mass on z-axis.
REAL(8) :: rdmsc    !<Charge diffraction radius.
REAL(8) :: rdmsp    !<Charge diffraction radius.
REAL(8) :: rdmsc0   !<Charge diffraction radius.
REAL(8) :: surfc    !<Surface thickness of charge distribution.
REAL(8) :: surfp    !<Surface thickness of charge distribution.
REAL(8) :: surfc0   !<Surface thickness of charge distribution.
REAL(8) :: xinert   !<Rigid-body momentum of inertia orthognal to z-axis.
REAL(8) :: bmassc   !<Cartesian collective quadrupole mass in scaling approach.
REAL(8) :: bmassr   !<Radial collective quadrupole mass in scaling approach.
REAL(8) :: zpetri   !<Estimate of zero-point energy in triaxial direction.
REAL(8) :: omgtri   !<Estimate of collective energy in triaxial direction.
REAL(8) :: widdxy   !<Width of center-of-mass momentum along r direction.
REAL(8) :: widdz    !<Width of center-of-mass momentum along z direction.
!REAL(8) :: zpemod   !<Difference of self-consistent c.m. correction and simple estimate (from 'widdz' and widdxy').
REAL(8) :: ezp      !<full ZPE from <P**2>/2mA.
REAL(8) :: cmwid    !<width for c.m. correction
REAL(8) :: zpecmx   !<Alternative center-of-mass energy along x direction.
REAL(8) :: widcmx   !<Alternative width of center-of-mass momentum along x direction.
REAL(8) :: zpecmz   !<Alternative center-of-mass energy along z direction.
REAL(8) :: widcmz   !<Alternative width of center-of-mass momentum along x direction.
REAL(8) :: bubdipp  !<Proton bubble parameter.
REAL(8) :: bubdipn  !<Neutron bubble parameter
REAL(8) :: bubdipt  !<Bubble parameter.
REAL(8) :: boundnum(2)!<Number of bound nucleons.
REAL(8), PUBLIC :: qcorrs(0:4)    ! collection of Q corrections in CHF
REAL(8), PUBLIC :: qlcstps(0:4)   ! collection of lambda_Q steps in CHF
REAL(8), PUBLIC :: varqscs(0:4)   ! collection of Q variance in CHF
REAL(8), PUBLIC :: qcons(0:4)     ! collection of Q deviation in CHF

!     the observables
 
REAL(8) :: econ     !<Energy convergence parameter.
REAL(8) :: eold     !<Energy from previous iteration.
REAL(8) :: energy   !<Total energy directly from densities.
REAL(8) :: energyFy !<Total energy directly from densities, Fayans functional
REAL(8) :: energr   !<Total energy from shell and rearrangement energy.
REAL(8) :: ekin     !<Kinetic energy.
REAL(8) :: epot     !<Potential energy.
REAL(8) :: ecoul    !<Coulomb energy.
REAL(8) :: etpair   !<Total pairing energy (or -T*S in case of thermal).
REAL(8) :: eshell   !<Shell energy.
REAL(8) :: erearr   !<Rearrangement energy.
REAL(8) :: equad    !<Energy contained in constraint.
REAL(8) :: qold     !<Quadrupole from previous iteration.
REAL(8) :: qcon     !<Change of quadrupole momentum in one iteration.
REAL(8) :: ecolex   !<Energy in Coulomb exchange.
REAL(8) :: devmf2   !<Average variance of s.p. energies.
REAL(8) :: devh2    !<Variance of mean-field operator.
REAL(8) :: devoc2   !<Content of h*psi in given wavefunctions.
REAL(8) :: devgap   !<Average deviation from gap equation.
REAL(8) :: varqsq   !<Variance of constraint momentum (for constraint iteration).
REAL(8) :: qq0      !<Variance of monopole moment
REAL(8) :: qq1      !<Variance of dipole moment
REAL(8) :: qq2      !<Variance of quarupole moment
REAL(8) :: qq3      !<Variance of octupole moment
REAL(8) :: qq4      !<Variance of hexadecupole moment
REAL(8) :: qdh      !<Element <Q D h> for L=2
REAL(8) :: symcondiso(2) !< HFB symmetry conditions for protons and neutrons
REAL(8) :: symcondi(2)   !< initial HFB symmetry condition in 'hfb_natural'
REAL(8) :: symcondf(2)   !< final HFB symmetry condition in 'hfb_natural'

REAL(8),ALLOCATABLE,SAVE :: rmsfield(:,:) ! field containing r**2


INTEGER,PARAMETER,PRIVATE :: isign_rso=1  !< +-1 to switch new/old magnetic contribution to radii
INTEGER,PARAMETER,PRIVATE :: kfrm=512     !<Maximum dimension of formfactor fields.

! parameters for nucleon contribution to charge rms radius
REAL(8),PARAMETER,PRIVATE :: xmgne=-0.04197d0,xmgpr=0.06131d0,d4m=0.010989d0
REAL(8),PARAMETER,PRIVATE :: r2intrp=0.7297D0,r2intrn=-0.1174D0

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: sinfo
!> @brief
!!Prints and calculates observables such as total energy, pairing energy, 
!!kinetic energy, deformation parameters, radii etc. If \c iter matches print 
!!modulus prints out detailed information and s.p. states. 
!> @details
!!The pairing energy is evaluated as follows:
!!
!!\f$ E_\mathrm{pair} \f$: the pairing energy. It consists of a
!!contact pairing interaction involving the pairing densities \f$ \xi_q \f$
!!augmented by an optional density dependence. The formula is
!!\f[ E_{\rm pair} = \frac{1}{4} \sum_{q\in\{p,n\}}V_\mathrm{pair,q}
!!\int \D^3r |\xi_q|^2 \left[1 -\frac{\rho}{\rho_{0,\mathrm{pair}}}\right]\;. \f]     
!!It contains a continuous switch, the parameter
!!\f$ \rho_{0,\mathrm{pair}} \f$. A pure \f$ \delta \f$-interaction (DI), also
!!called volume pairing, is recovered for
!!\f$ \rho_{0,\mathrm{pair}}\longrightarrow\infty \f$. The general case is
!!the density dependent \f$ \delta \f$-interaction (DDDI). A typical
!!value near matter equilibrium density
!!\f$ \rho_{0,\mathrm{pair}}=0.16 \f$ fm\f$ ^{-3} \f$ concentrates pairing to
!!the surface. The most flexible choice is to consider
!!\f$ \rho_{0,\mathrm{pair}} \f$ as an additional free parameter. Actual
!!adjustments with this option deliver a form of the pairing
!!functional which stays in between the extremes of volume and
!!surface pairing.
!>
!> @param[in] iter
!> INTEGER, takes the current iteration number; iter=itmax signals final evaluation.
!--------------------------------------------------------------------------- 

SUBROUTINE sinfo(iter)

USE paramr
USE dynpar
USE grid, ONLY: nr,nz,icm,tsym,sr,sz,f1int
USE constr
USE wavefr, ONLY: psi,psifold
USE fields
USE levels
!USE steputil, ONLY: cpairpot


IMPLICIT NONE

INTEGER, INTENT(IN)                  :: iter

INTEGER, PARAMETER :: iq2old=0    ! switch r_rms <--> R_dms in beta_i
CHARACTER(LEN=4),PARAMETER :: itype(2,2)=RESHAPE((/'prot','on  ','neut','ron '/),(/2,2/)) 
CHARACTER(LEN=40) :: format_conver

REAL(8) :: dipout,z,v2,uv,epairtestn,epairtestp,etot1,fac,occ,occd,qs,&
           qsn,qsp,qsqn,qsqp,r2h,r2t,rmsb,rmsbn,rmsbp,r,qsq
INTEGER :: ind,ir,iq,ist,iz,jst,nstm,ine,lst
LOGICAL :: topened

INTEGER,SAVE :: lconp(0:4)     ! pointer to active constraints
INTEGER,SAVE :: maxcon       ! max. nr. of active constraints
REAL(8) :: qval(4)           ! storage for actual multipole mooments
REAL(8),ALLOCATABLE ::  work(:,:)
REAL(8), ALLOCATABLE :: paircrtest(:,:,:),pairpttest(:,:,:)
REAL(8) :: densint0,densintfold,densintdiff, epairtestnfold
REAL(8) :: DelQ2(-2:2)
LOGICAL, PARAMETER :: testfolding=.FALSE.

!**********************************************************************
!         block to compute the shell energy, the field energy         *
!         and the rearrangement energy                                *
!**********************************************************************

! repair possibly wrong signs between pairing conjugated states
! here as an occasional safety measure
IF(toddterms .AND. iter==iterbcs+1 .AND. iter==iterhfb+1)&
  CALL testconj_all('in SINFO')


ALLOCATE(work(0:nz,0:nr))
IF(testfolding) ALLOCATE(paircrtest(0:nz,0:nr,kqm),pairpttest(0:nz,0:nr,kqm))
!  constraint energy

work=0D0
IF(tswicon0) work = work + qlamb0*qfac0(0:nz,0:nr)
IF(tswicon2) work = work + qlamb2*qfac2(0:nz,0:nr)
IF(tswicon3) work = work + qlamb3*qfac3(0:nz,0:nr)
IF(tswicon4) work = work + qlamb4*qfac4(0:nz,0:nr)

work = (rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2))*work

equad  = f1int(work(0:nz,0:nr))

!     the zero point energy in a rough harmonic oscillator approximation

IF(ifzpe==1) THEN
  IF(iter<itzpe) THEN
    ezp = 17.3/(npart(1)+npart(2))**0.2 ! simple estimate for inital iterations
!    WRITE(7,*) ' approximate ZPE=',iter,ezp
  ELSE IF(iter==itmax .OR. MOD(iter,itzpe)==0) THEN
    CALL zpe_cm()                       ! true c.m. correction 
!    WRITE(7,*) ' full ZPE=',iter,ezp
  END IF
ELSE
  ezp = 0.0D0
END IF

!**********************************************************************
!                                                                     *
!     Calculate entropy for ipair=4 (thermal CHF)                     *
!                                                                     *
!                                                                     *
!**********************************************************************

IF(ipair == 4) THEN                         ! compute entropy
  epair(1:2) = zero
  DO ist=1,nst(1)+nst(2)
    iq   = icharg(ist)
    occ  = MAX(vsqf(ist),smal)
    occd = MAX(one-occ,smal)
    epair(iq) = wate(ist)*(occ*LOG(occ)+occd*LOG(occd)) +epair(iq)
  END DO
  epair(1) = gap(1)*epair(1)                 ! -T*S
  epair(2) = gap(2)*epair(2)
END IF

paircr(0:nz,0:nr,1:2)  = 0D0
IF(testfolding) paircrtest(0:nz,0:nr,1:2)  = 0D0
!     accumulate new pair-density
!NEW ODD
IF(thfbfold) THEN
  DO jst=1,nst(1)+nst(2)
    iq    = icharg(jst)
    IF(jst==istblo(iq) .OR. jst==istcblo(iq)) CYCLE    !NEW ODD
    v2    = vsqf(jst)
    uv    = wtpair_uv(jst)*wtpair_v2(jst)*SQRT(MAX(v2-v2*v2,1D-11))*0.5
!    CALL Gaussfolding_spin(psi(0:nz,0:nr,1:2,jst),psifold(0:nz,0:nr,1:2,jst),nparz(1,jst),mpsi(jst),1)
    IF(testfolding) CALL sdenti(psi(0:nz,0:nr,1:2,jst),paircrtest(0:nz,0:nr,iq),uv,jst)
    IF(toddterms) THEN
      IF(mpsi(jst).lt.0) uv = -uv
      CALL cpairpot(psifold(0:nz,0:nr,1:2,jst),psifold(0:nz,0:nr,1:2,itpartn(jst)),paircr(0,0,iq),uv,jst)
    ELSE
      CALL sdenti(psifold(0:nz,0:nr,1:2,jst),paircr(0:nz,0:nr,iq),uv,jst)
    END IF
  END DO
ELSE
  DO jst=1,nst(1)+nst(2)
    iq    = icharg(jst)
    IF(jst==istblo(iq) .OR. jst==istcblo(iq)) CYCLE    !NEW ODD
    v2    = vsqf(jst)
    uv    = wtpair_uv(jst)*wtpair_v2(jst)*SQRT(MAX(v2-v2*v2,1D-11))*0.5
    IF(toddterms) THEN
      IF(mpsi(jst).lt.0) uv = -uv
      CALL cpairpot(psi(0:nz,0:nr,1:2,jst),psi(0:nz,0:nr,1:2,itpartn(jst)),paircr(0,0,iq),uv,jst)
    ELSE
      CALL sdenti(psi(0:nz,0:nr,1:2,jst),paircr(0:nz,0:nr,iq),uv,jst)
    END IF
  END DO
END IF
!END NEW ODD

!c     now multiply with strength to obtain local pair-potential
IF(.NOT.tFayans) THEN
  IF(ipair <= 6 .OR. ipair == 9) THEN
    pairpt(0:nz,0:nr,1)  = -pairfc(1) * paircr(0:nz,0:nr,1)
    pairpt(0:nz,0:nr,2)  = -pairfc(2) * paircr(0:nz,0:nr,2)
  ELSE
    work(0:nz,0:nr)=one-(rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2))/rho0pr
    pairpt(0:nz,0:nr,1)  = -pairfc(1) * work * paircr(0:nz,0:nr,1)
    pairpt(0:nz,0:nr,2)  = -pairfc(2) * work * paircr(0:nz,0:nr,2)
    IF(testfolding) THEN
      pairpttest(:,:,1)  = -pairfc(1) * work * paircrtest(:,:,1)
      pairpttest(:,:,2)  = -pairfc(2) * work * paircrtest(:,:,2)
    END IF
  END IF

  IF(ipair <= 6 .OR. ipair == 9) THEN
    work=-pairfc(1)* paircr(0:nz,0:nr,1)* paircr(0:nz,0:nr,1)
    epairtestp=f1int(work)  !0.25 oder 0.5,was ist v0 ?
    work=-pairfc(2)* paircr(0:nz,0:nr,2)* paircr(0:nz,0:nr,2)
    epairtestn = f1int(work)  !0.25 oder 0.5,was ist v0 ? 
  ELSE
    work=one-(rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2))/rho0pr
    work=-pairfc(1)* paircr(0:nz,0:nr,1)* paircr(0:nz,0:nr,1)* work
    epairtestp=f1int(work)  !0.25 oder 0.5,was ist v0 ?
    work=one-(rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2))/rho0pr
    work=-pairfc(2)* paircr(0:nz,0:nr,2)* paircr(0:nz,0:nr,2)* work
    epairtestn = f1int(work)  !0.25 oder 0.5,was ist v0 ?
  END IF
  etpair= epairtestp+ epairtestn
  epair(1) = epairtestp
  epair(2) = epairtestn
  IF(ipair <= 4)  etpair=epair(1)+epair(2)
END IF

IF(testfolding) THEN
  densint0 = f1int(paircrtest(:,:,2))
  densintfold = f1int(paircr(:,:,2))
  work(:,:)=paircrtest(:,:,2)*pairpttest(:,:,2)
  epairtestnfold = f1int(work)  
  pairpttest = (pairpt-pairpttest)**2
  paircrtest = (paircr-paircrtest)**2
  densintdiff = SQRT(f1int(paircrtest(:,:,2)))
  WRITE(6,'(a,8(1pg15.7))') 'COMPOBS densities:', &
     densint0,densintfold,densintdiff,densint0-densintfold, &
     epairtestn,epairtestnfold,epairtestn-epairtestnfold, &
     f1int(pairpttest(:,:,2))
  OPEN(317,file='pairpttest')
    WRITE(317,'(a)') '# paipttest observ:'
    DO ir=0,nr
      DO iz=0,nz
        WRITE(317,'(2i3,2(1pg13.5))') iz,ir,pairpttest(iz,ir,2),paircrtest(iz,ir,2)
      END DO
      WRITE(317,'(1x)')
    END DO
  CLOSE(317)
  DEALLOCATE(paircrtest,pairpttest)
END IF


!**********************************************************************
!     finally the total energy                                        *
!**********************************************************************
IF(itest >= 1) THEN
  CALL energ2()         ! energy using rearrangement energy
  energr = energr - ezp + etpair
END IF
IF(tFayans) THEN
  etpair=SUM(epair)
  ekin   = h2meff(1)*f1int(tau(:,:,1))+h2meff(2)*f1int(tau(:,:,2))
  etot1  = ekin+energyFy+ecoul/2D0
  energy = ekin+energyFy+ecoul/2D0-ezp+etpair
!  IF(itest>0) WRITE(7,'(a,i5,6(1pg13.5))') &
!  WRITE(7,'(a,i5,6(1pg13.5))') &
!      'energyFy=',iter,etot1,ekin,energyFy,epair,ezp,f1int(rho(:,:,1))
ELSE
  CALL energ()           ! energy directly from density
  etot1  = energy
  energy = energy-ezp+etpair
END IF
IF(itest >= 1 .OR. iter==itmax) THEN
  WRITE(6,'(a,i5,a,2f14.4)') ' iter=',iter,  &
      ', energies: from rearr. vs. from densities=',energr,energy
  WRITE(7,'(a,2f12.4/a,4(1pg12.5)/a,5(1pg12.5)/a,i5,1pg13.5)')  &
      '  energies: from rearr. vs. from densities=',energr,energy,  &
      '  ekin,epot,etot=',ekin,epot,etot1,energy,  &
      '  epair-direct  =',epairtestp,epairtestn,etpair,epair, &
      '  sum_s.p.energ.=',iter,SUM(spenrg(1:nst(1))*wate(1:nst(1))*vsqf(1:nst(1)))
END IF
!Copy previous s.p. energies to speold
lst = nst(1)+nst(2)
speold(1:lst) = spenrg(1:lst)

!**********************************************************************
!        Calculate deformations                                       *
!**********************************************************************
zcm = zcmobs(rho(0:nz,0:nr,:))
CALL obs(rho(0:nz,0:nr,1),apnum,rmsp,dp,qp,op,hp,qsqp,zcm)
CALL obs(rho(0:nz,0:nr,2),annum,rmsn,dn,qn,on,hn,qsqn,zcm)
!Determine deformations for total nucleus
atnum = apnum+annum
rms   = SQRT((rmsp*rmsp*apnum+rmsn*rmsn*annum)/atnum)
d     = dp+dn
q     = qp+qn
o     = op+on
h     = hp+hn
qsq   = qsqp+qsqn
qval(1) = rms
qval(2) = q
qval(3) = o
qval(4) = h
!WRITE(*,*) 'd,q,o,h=',d,q,o,h
!WRITE(7,'(a,3f6.1,3f7.4)') 'rms:',apnum,annum,atnum,rmsp,rmsn,rms
IF(itest > 0) THEN
   WRITE(6,'(a,7(1pg12.4))')  &
    'pr: Z,rms,d,q,o,h=',apnum,rmsp,dp,qp,op,hp
   WRITE(6,'(a,7(1pg12.4))')  &
    'ne: N,rms,d,q,o,h=',annum,rmsn,dn,qn,on,hn
END IF
!Calculate deformation parameters
fac     = SQRT(0.2D0*pi)
IF(iq2old == 1) THEN
  rmsb   = rms
  rmsbp  = rmsp
  rmsbn  = rmsn
ELSE
  rmsb   = 0.929516003D0*atnum**0.333333333333D0
  rmsbp  = rmsb
  rmsbn  = rmsb
END IF
beta2  = fac*q /MAX(rmsb*rmsb*atnum,smal)
beta2p = fac*qp/MAX(rmsbp*rmsbp*apnum,smal)
beta2n = fac*qn/MAX(rmsbn*rmsbn*annum,smal)
beta4  = h * 1.507964474/(MAX(atnum*rmsb**4,smal))
beta4p = hp* 1.507964474/(MAX(apnum*rmsbp**4,smal))
beta4n = hn* 1.507964474/(MAX(annum*rmsbn**4,smal))
IF(.NOT.tsym) THEN
  beta1  = d * 3.244622941/MAX(atnum*rmsb,smal)
  beta1p = dp* 3.244622941/MAX(apnum*rmsbp,smal)
  beta1n = dn* 3.244622941/MAX(annum*rmsbn,smal)
  beta3  = o * 1.946773765/MAX(atnum*rmsb**3,smal)
  beta3p = op* 1.946773765/MAX(apnum*rmsbp**3,smal)
  beta3n = on* 1.946773765/MAX(annum*rmsbn**3,smal)
ELSE
  beta1  = 0D0
  beta1p = 0D0
  beta1n = 0D0
  beta3  = 0D0
  beta3p = 0D0
  beta3n = 0D0
END IF
! fragment separation
IF(iter==itmax .OR. (iter>0 .AND. MOD(iter,iprint)==0)) &
               rsepfrag = fragdist(rho(0:nz,0:nr,:))
! cartesian isovector dipole
dipout = SQRT(4D0*pi/3D0)/anumb*(npart(2)*dp-npart(1)*dn)
! bubbleness
CALL bubble()
! rigid-body momentum of inertia for rotation about the x-axis
DO ir=0,nr
  r     = dr*(ir+sr)
  r2h   = r*r*0.5
  DO iz=0,nz
    z     = ((iz)+sz-icm)*dz
    r2t   = z*z+r2h
    work(iz,ir) = (rho(iz,ir,1)+rho(iz,ir,2))*r2t
  END DO
END DO
xinert  = f1int(work(0:nz,0:nr))/h2mav
!     collective quadrupole mass in scaling approach
bmassc = 8.0*h2mav*(2.0*rms*rms*atnum+q)              ! cartesian
bmassr = bmassc*(fac/(rms*rms*atnum))**2              ! r-quadr.

DEALLOCATE(work)
!**********************************************************************
!        Formfactor analysis and radii wih c.m. correction            *
!**********************************************************************
IF(iter==itmax) CALL radius()
!**********************************************************************
!        begin block to write out the observables                     *
!**********************************************************************

IF(iter == 0 .AND. tpri_conver) THEN
    WRITE(71,1010) npart(1)+npart(2),npart(1), dz,dr,nz,nr,  &
        addrho,addpot,addocc,c0,qeps, x0inv,e0inv
    1010   FORMAT(' Nmass=',i5,', Ncharge=',i5,/,  &
        '  dz,dr: ',2F6.3,'   nz,nr=',2I5,/,  &
        ' addrho,addpot,addocc: ',3F6.3,'   c0,qeps=',2F6.3,/,  &
        ' x0inv(initial),e0inv(mev): ',2F7.2)
  IF(.NOT.tswicon) THEN
    WRITE(71,'(a,7x,a)')  &
        ' iter   energy    variance  difference    qtotal',  &
        '   devoc2    devmf2      rmsp     x0inv   epair(1:2)'
!        '   devoc2    devmf2     devgap    x0inv   epair(1:2)'
  ELSE
    ! prepare pointers for printing
    lconp = 0
    maxcon = 0
    IF(tswicon0) THEN
      maxcon=1+maxcon
      lconp(maxcon)=0
    ELSE IF(tswicon2) THEN
      maxcon=1+maxcon
      lconp(maxcon)=2
    ELSE IF(tswicon3) THEN
      maxcon=1+maxcon
      lconp(maxcon)=3
    ELSE IF(tswicon4) THEN
      maxcon=1+maxcon
      lconp(maxcon)=4
    END IF
    IF(maxcon==1) format_conver='(a,4(4(a,i1)),a)'
    IF(maxcon==2) format_conver='(a,8(4(a,i1)),a)'
    IF(maxcon==3) format_conver='(a,12(4(a,i1)),a)'
    IF(maxcon==4) format_conver='(a,16(4(a,i1)),a)'
    WRITE(71,format_conver)  &
        ' iter  energy    variance  difference   x0inv   ',&
        (' qtotal(',lconp(ir),')  qcorr(',lconp(ir),&
         ')   del-qlamb(',lconp(ir),') qcon(',lconp(ir),')',&
        ir=1,maxcon)
  END IF
ELSE IF(tpri_conver) THEN
  IF(.NOT.tswicon) THEN
    IF(MOD(iter,itzpe).NE.1) WRITE(71,'(i5,f11.4,9(1pg13.5))')  &
        iter,energy,devh2,econ,q,devoc2,devmf2,rmsp,x0inv,epair
!        iter,energy,devh2,econ,q,devoc2,devmf2,SQRT(ABS(devgap)),x0inv,epair
  ELSE
    IF(maxcon==1) format_conver='(i5,f11.4,3(1pg10.3),4(X,4(1pg12.3)))'
    IF(maxcon==2) format_conver='(i5,f11.4,3(1pg10.3),8(X,8(1pg12.3)))'
    IF(maxcon==3) format_conver='(i5,f11.4,3(1pg10.3),12(X,12(1pg12.3)))'
    IF(maxcon==4) format_conver='(i5,f11.4,3(1pg10.3),16(X,16(1pg12.3)))'
    WRITE(71,format_conver)  &
        iter,energy,devh2,econ,x0inv,&
        (qval(lconp(ir)),qcorrs(lconp(ir)),qlcstps(lconp(ir)),&
         qcons(lconp(ir)),ir=1,maxcon)
  END IF
END IF
!IF(energy /= zero) econ   = ABS((energy-eold)/energy)
IF(energy /= zero .AND. iter<itmax) econ   = ABS((energy-eold)/energy)

!copy total energy to storage for previous iteration
eold   = energy


!**********************************************************************
!        only write full information every iprint iterations          *
!**********************************************************************
IF((iter == 0 .OR. iprint == 0 .OR. MOD(iter,iprint) /= 0) &
   .AND. iter.NE.itmax) RETURN
WRITE(7,25010) iter,energy,energy/MAX(atnum,smal),eshell,ecoul/2D0,ecolex,erearr,  &
    ezp,epot,ekin,etpair,equad,econ,symcondiso
25010 FORMAT('1'/,' iteration no. = ',i10, /,' binding energy (MeV) = ',  &
    f10.3,' binding energy per nucleon (MeV) = ',f10.3,/, '   eshell = ',f10.3,  &
    '   ecoul  = ',f10.3, '   ecolex = ',f10.3,/,  &
    '   erearr = ',f10.3, '   ezp    = ',f10.3,  &
    '   epot   = ',f10.3,/, '   ekin   = ',f10.3,  &
    '   epair  = ',f10.3, '   equad  = ',f10.3,/,  &
    '   de/e   = ',1PG12.4,' symcond =',2(1pg12.4))
IF(ipair >= 0) WRITE(7,'(3(a,2(1pg10.3)))')  &
    '   gap(1/2)=',gap,'   fer(1/2)=',fer
IF (tswicon0)  WRITE(7,'(a,1pg14.6)')' L=0 Lagrange parameter:',qlamb0
IF (tswicon2)  WRITE(7,'(a,1pg14.6)')' L=2 Lagrange parameter:',qlamb2
IF (tswicon3)  WRITE(7,'(a,1pg14.6)')' L=3 Lagrange parameter:',qlamb3
IF (tswicon4)  WRITE(7,'(a,1pg14.6)')' L=4 Lagrange parameter:',qlamb4

!        cartesian quadrupole moments are transformed
!        in spherical ones (qs,qsp,qsn).
qs  = q*SQRT(5.0/pi)/4.0
qsp = qp*SQRT(5.0/pi)/4.0
qsn = qn*SQRT(5.0/pi)/4.0
WRITE(7,'(/a,a)')  &
    '         part.no.  radius   dipol       quad        octu  ',  &
    '     hexadeca     quad_cart'
WRITE(7,'(a,2f8.3,5(1pg12.4))') '  proton:',apnum,rmsp,dp,qsp,op,hp,qp
WRITE(7,'(a,2f8.3,5(1pg12.4))') ' neutron:',annum,rmsn,dn,qsn,on,hn,qn
WRITE(7,'(a,2f8.3,5(1pg12.4))') '   total:',atnum,rms,d,qs,o,h,q
WRITE(7,'(a,8x,(1pg12.4))') '   diff:',rmsn-rmsp
IF(ABS(dipin)>0D0) &
  WRITE(7,'(/a,2i4,3(1pg20.12))') ' dipoles:  Z,N,cm,isovector,polariz=',  &
    npart(1),npart(2),zcm,dipout,-e2*dipout/dipin
WRITE(7,'(/a)') '            beta_1      beta_2      beta_3      beta_4  '
WRITE(7,'(a,4(1pg12.4))') '  proton:',beta1p,beta2p,beta3p,beta4p
WRITE(7,'(a,4(1pg12.4))') ' neutron:',beta1n,beta2n,beta3n,beta4n
WRITE(7,'(a,4(1pg12.4))') '   total:',beta1,beta2,beta3,beta4
IF(ipair /= -1.AND.iter > itpini) THEN
  WRITE(7,'(/11x,2a)') 'e-pair      e-fermi  ',  &
      '   gap  '
  WRITE(7,'(a,4(1pg12.4))') '  proton:',epair(1),fer(1),gap(1)
  WRITE(7,'(a,4(1pg12.4))') ' neutron:',epair(2),fer(2),gap(2)
  WRITE(7,'(a,5(1pg12.4))') '   total:',etpair
END IF
!     optionally write polarizability

IF(dipin.NE.0D0) THEN
  INQUIRE(49,OPENED=topened)
  IF(.NOT.topened) THEN
    OPEN(49,STATUS='unknown',FORM='formatted',FILE=TRIM(fpath('polariz'))//'polariz.'//outnam)
    WRITE(49,'(a)') 'nprot nneut beta_2 dipin   dipout    polarizability'
  END IF
  WRITE(49,'(2i5,f8.3,f9.5,2(1pg20.12))')  &
    npart(1),npart(2),beta2,dipin,dipout,-e2*dipout/dipin
  CALL flush(49)
END IF

!**********************************************************************
!        begin block to write out the single particle observables     *
!        note that the ordering in energy is not the natural ordering *
!        of the levels                                                *
!**********************************************************************
DO  iq=1,2
  nstm=nst(iq)
  WRITE(7,25030) itype(1,iq),itype(2,iq),nst(iq)
  25030 FORMAT(//2A4,' levels:',5X,i4,'  maximum')
  
!     print levels
  IF(ipair < 6) THEN
    WRITE(7,'(/18x,a,9x,a/a,a)') 'parities','components',  &
        '   Nr Iso   j_z total   r    z  wtpair  v^2  1.   2.', '   energy'
    boundnum(iq) = 0D0
    DO ist=1,nstm
      ine    = iokk(iq)+ist-1
      ind    = iorder(ine)+iokk(iq)-1
      IF(tcutoff_uvspace) THEN
        WRITE(7,'(i5,i4,f6.1,3i5,2f8.5,2f5.2,f9.3,i5)')  &
          ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),  &
          nparz(1,ind),wtpair_uv(ind), vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind),itpartn(ind)
      ELSE
        WRITE(7,'(i5,i4,f6.1,3i5,2f8.5,2f5.2,f9.3,i5)')  &
          ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),  &
          nparz(1,ind),wtpair_v2(ind), vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind),itpartn(ind)
      END IF
      IF(spenrg(ind)<0D0) boundnum(iq)=boundnum(iq)+vsqf(ind)*wate(ind)
    END DO
    WRITE(7,'(a,1pg13.5)') '# nr. of bound nucleons=',boundnum(iq)
    WRITE(7,'(a)') ' '
  ELSE
    WRITE(7,'(/10x,a,16x,a/a,a)') 'parities','components',  &
        '   Nr Iso j_z tot r z  wtpair  v^2      uv  delta',&
        '   energy  E-variance  ang.mom.  s.p.radius conj.state'
    boundnum(iq) = 0D0
    DO ist=1,nstm
      ine    = iokk(iq)+ist-1
      ind    = iorder(ine)+iokk(iq)-1
!      WRITE(6,'(a,3i5,1pg13.5)') 'testconj:',ind,itpartn(ind),mpsi(ind), &
!       testconj(psi(:,:,:,ind),psi(:,:,:,itpartn(ind)))*SIGN(1D0,mpsi(ind)+0.5D0)
      IF(tcutoff_uvspace) THEN
        WRITE(7,'(i5,i2,f5.1,3i3,2f8.5,2f6.2,f9.3,1pg12.4,0pf6.3,f8.3,i5)')  &
          ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),  &
          nparz(1,ind),wtpair_uv(ind),  &
          vsqf(ind),SQRT(MAX(ABS(vsqf(ind)-vsqf(ind)**2),smal)),  &
          deltaf(ind)*wtpair_uv(ind),spenrg(ind),spevar(ind), &
          angularmom(ind), spradius(psi(:,:,:,ind)),itpartn(ind) !spevhh(ind)
      ELSE
        WRITE(7,'(i5,i2,f5.1,3i3,2f8.5,2f6.2,f9.3,1pg12.4,0pf6.3,f8.3,i5)')  &
          ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),  &
          nparz(1,ind),wtpair_v2(ind),  &
          vsqf(ind),SQRT(MAX(ABS(vsqf(ind)-vsqf(ind)**2),smal)),  &
          deltaf(ind),spenrg(ind),spevar(ind),  &
          angularmom(ind),spradius(psi(:,:,:,ind)),itpartn(ind) !spevhh(ind)
      END IF
      IF(spenrg(ind)<0D0) boundnum(iq)=boundnum(iq)+vsqf(ind)*wate(ind)
    END DO
    WRITE(7,'(a,1pg13.5)') '# nr. of bound nucleons=',boundnum(iq)
    WRITE(7,'(a)') ' '
  END IF
  
END DO

! compute and print Q variance
CALL Qvariance(DelQ2)
rmsb = 0.929516003D0*atnum**0.333333333333D0
WRITE(6,'(a,5(1pg13.5))') 'Delta Q^2:',DelQ2
WRITE(7,'(a,5(1pg13.5))') 'Delta Q^2:',DelQ2
!WRITE(7,'(a,1pg13.5)') 'atnum=',atnum
DelQ2 = DelQ2*(4D0*pi/MAX(5D0*rmsb*rmsb*atnum,smal))**2  ! convert to dimensionless
DO ind=6,7
  WRITE(ind,'(a,f7.3,5(1pg13.5),5x,1pg13.5)') 'Delta_beta^2:  ',beta2,DelQ2,SUM(DelQ2)
  qsq = SQRT(SUM(DelQ2)+beta2**2)
  WRITE(ind,'(a,f7.3,5(1pg13.5),5x,2(1pg13.5))') 'SQRT(D_beta^2):',&
     beta2,SQRT(DelQ2),SQRT(SUM(DelQ2)),qsq
!  WRITE(ind,'(a,f7.3,6(1pg13.5))') 'Delta_SQRT(sum_beta^2):',beta2,SQRT(SUM(DelQ2)),qsq
!  WRITE(ind,'(a,f7.3,5(1pg13.5))') 'Delta_SQRT(tot_beta^2):',beta2,SQRT(SUM(DelQ2))
END DO



!     optionally save occupation numbers
IF(ifoccs == 1 .AND. ifocsv == 0) THEN
  DO iq=1,2
    nstm=nst(iq)
    DO ist=1,nstm
      ine     = iokk(iq)+ist-1
      oc(ine) = vsqf(ine)
    END DO
  END DO
  ifocsv = 1
  WRITE(6,'(a)') ' set flag IFOCSV=1'
END IF

RETURN
END SUBROUTINE sinfo

!---------------------------------------------------------------------------  
! DESCRIPTION: write_ensurf()
!> @brief
!! Writes to protocol file for basic PES and other gross properties
!> 
!> @param[in] iter
!> INTEGER, takes the current iteration number.
!---------------------------------------------------------------------------  
SUBROUTINE write_ensurf(iter)

USE paramr
USE grid
USE dynpar
!USE constr, ONLY: iflmb0
USE fields, ONLY: rho
USE levels, ONLY: npart,kblock,angblock

INTEGER, INTENT(IN) :: iter
!LOGICAL :: topen173

REAL(8) :: efree

  OPEN(90,STATUS='unknown',FORM='formatted',FILE=TRIM(fpath('ensurf'))//'ensurf.'//outnam,  &
          POSITION='APPEND')
  IF(ipair /= 4) THEN
    IF(.NOT.tsym) THEN
      WRITE(90,'(12f8.3,f18.11,2f7.1,2f7.2,2f9.3,i6,2f5.1,2f6.2,f16.12,3f8.4,2f11.6,2(1pg12.4))')  &
          beta1,beta1p,beta1n,beta2,beta2p,beta2n,  &
          beta3,beta3p,beta3n,beta4,beta4p,beta4n,  &
          energy,epair,fer,gap,iter,kblock*0.5D0,angblock,rmsc,rms,rmsp,rmsn,boundnum,devh2,econ
    ELSE
      WRITE(90,'(6f7.3,f18.11,2f7.1,2f7.2,2f9.3,i6,2f5.1,2f6.2,f16.12,3f8.4,2f11.6,2(1pg12.4))') &
          beta2,beta2p,beta2n,beta4,beta4p,beta4n,  &
          energy,epair,fer,gap,iter,kblock*0.5D0,angblock,rmsc,rms,rmsp,rmsn,boundnum,devh2,econ
    END IF
    IF(ABS(apnum-npart(1)) > 1.0E-5  &
        .OR. ABS(annum-npart(2)) > 1.0E-5 ) WRITE(90,'(35x,a,2(1pg14.6))')  &
        'part.nr.=',apnum,annum
  ELSE
    efree  = energy
    energy = efree-etpair
    WRITE(90,'(f6.2,4f8.3,3f12.3,f7.2)')  &
        tempac,beta2,beta2p,beta4,beta4p,energy,efree,-etpair,  &
        -3.0*etpair/(tempac*pi)**2
  END IF
!   IF(iflmb0 /= 1) CLOSE(90)
!  CALL flush(90)
  CLOSE(90) 
!  INQUIRE(UNIT=173,OPENED=topen173)
  IF(tpri_bubble) THEN
    CALL bubble()
    WRITE(173,'(f7.3,f10.3,5(1pg13.5))')  &
        beta2,energy,rho(icm,0,1),rho(icm,0,2),bubdipt,bubdipp,bubdipn
    CALL flush(173)
  END IF
  CLOSE(90)      


  IF(tpri_radii)  THEN
    OPEN(76,STATUS='unknown',FORM='formatted',FILE=TRIM(fpath('radii'))//'radii.'//outnam,  &
          POSITION='APPEND')
    WRITE(76,'(f8.4,f12.3,4f8.4)') beta2,energy,rmsc,rdmsc,surfc,r4c
    CLOSE(76)
  END IF


  RETURN

END SUBROUTINE write_ensurf

!---------------------------------------------------------------------------  
! DESCRIPTION: energ
!> @brief
!!Computes the total energy directly from densities.
!> @details
!!The different contribution to the total energy are
!!  - \f$ T \f$: the total kinetic energy calculated as
!!    \f[ T=\sum_q \frac{\hbar^2}{2m_q}\int\D^3 r\, \tau_q \f]
!!    with \f$ \tau_q \f$ the isospin-specific kinetic density.
!!  - \f$ E_0 \f$: The \f$ b_0 \f$ and \f$ b'_0 \f$-dependent part is
!!    \f[ E_0=\int \D^3\!r\,\left(\frac{b_0}{2}\rho^2-\frac{b_0'}{2}\sum_q\rho_q^2\right). \f]
!!  - \f$ E_{1}\f$: kinetic terms containing the coefficients \f$ b_1 \f$ and \f$ b_1' \f$:
!!    \f[ E_1=\int \D^3\!r\,\left(b_1[\rho\tau-{\vec\jmath\,}{}^2]
!!        -b_1'\sum_q[\rho_q\tau_q-{\vec\jmath_q}{}^2]\right).\f]
!!  - \f$ E_{2} \f$: terms containing the coefficients \f$ b_2 \f$ and \f$ b_2' \f$.
!!    They involve the Laplacians of the densities.
!!    \f[ E_2=\int \D^3\!r\,\left(-\frac{b_2}{2}\rho\Delta\rho
!!        +\frac{b_2'}{2}\sum_q\rho_q\Delta\rho_q\right). \f]
!!  - \f$ E_{3} \f$: The many-body contribution is given by
!!    \f[ E_{3}=\int \D^3\!r\,\left(\frac{b_3}{3}\rho^{\alpha+2}
!!        -\frac{b_3'}{3}\rho^\alpha\sum_q\rho_q^2\right). \f]
!!  - \f$ E_\mathrm{ls} \f$: the spin-orbit energy
!!    \f[ E_\mathrm{ls}=\int \D^3\!r\,\biggl(
!!        -b_4[\rho\nabla\cdot\vec J+
!!        \vec s\cdot(\nabla\times\vec\jmath)]-b_4'\sum_q[\rho_q\nabla\cdot\vec J_q+
!!        \vec s_q\cdot(\nabla\times\vec\jmath_q)]\biggr) \f]
!!  - \f$ E_\mathrm{C} \f$: the Coulomb energy. It consists of the
!!    standard expression for a charge distribution in its own field
!!    (Hartree term) plus the exchange term in the Slater approximation. The formula is
!!    \f[ E_\mathrm{C} = \frac{e^2}{2}\int \D^3\!r \D^3\!r'
!!      \frac{\rho_p(\vec{r})\rho_p(\vec{r'})}{|\vec{r}-\vec{r'}|}
!!      -\int \D^3\!r\,\frac{3e^2}{4}\left(\frac3{\pi}\right)^{\tfrac1{3}}\rho_p^{4/3} \f]
!!    where \f$ e \f$ is the elementary charge with \f$ e^2=1.43989 \f$ MeV\f$ \cdot \f$fm.
!--------------------------------------------------------------------------- 

SUBROUTINE energ()

USE paramr
USE grid
USE dynpar
USE fields

IMPLICIT NONE


LOGICAL :: talpinv,tdens2

REAL(8) :: epotcr   !> Correction of potential energy for partial integration on the grid.
REAL(8) :: ekincr   !> Correction of kinetic energy for partial integration on the grid.
REAL(8) :: ts,tp,tn,alpeff,d2n,d2s,d2p,potn,potp,&
           rlsrn,rlsrp,rlszn,rlszp,rn,rsfac,rsinv,rspow,rspow7,ssn,ssp,&
           rs,rp,sss
           
INTEGER :: ir,iz

REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:),wrk3(:,:),wrk4(:,:)

! override by a RETURN if external field is wanted --> via user interfaces

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr),wrk3(0:nz,0:nr),wrk4(0:nz,0:nr))
wrk1=0D0;wrk2=0D0;wrk3=0D0;wrk4=0D0

IF(ABS(b7) > 1D-10) THEN                         !new
  tdens2 = .true.
ELSE
  tdens2 = .false.
END IF
IF(ABS(alpinv) > 1D-10) THEN
  talpinv = .true.
ELSE
  talpinv = .false.
END IF                                             !new

DO ir=0,nr
  DO iz=0,nz
    rp      = rho(iz,ir,1)
    rn      = rho(iz,ir,2)
    rs      = MAX(rn+rp,smal)
    tp      = tau(iz,ir,1)
    tn      = tau(iz,ir,2)
    ts      = tp+tn
    d2p     = d2rho(iz,ir,1)
    d2n     = d2rho(iz,ir,2)
    d2s     = d2p+d2n
    ssp     = drhols(iz,ir,1)
    ssn     = drhols(iz,ir,2)
    sss     = ssp+ssn
    rlsrp   = rholsr(iz,ir,1)
    rlszp   = rholsz(iz,ir,1)
    rlsrn   = rholsr(iz,ir,2)
    rlszn   = rholsz(iz,ir,2)
    
    wrk2(iz,ir) = h2meff(1)*tp+h2meff(2)*tn
    wrk1(iz,ir) = + 0.5*(b0*rs*rs-b0p*(rp*rp+rn*rn))  &
        + b1*rs*ts-b1p*(rp*tp+rn*tn) - 0.5*(b2*rs*d2s-b2p*(rp*d2p+rn*d2n))  &
        - b4*rs*sss-b4p*(rp*ssp+rn*ssn)  &
        - 0.5*t1p2x*((rlsrp+rlsrn)*(rlsrp+rlsrn) +(rlszp+rlszn)*(rlszp+rlszn))  &
        + 0.5*t1m20*(rlsrp*rlsrp+rlsrn*rlsrn +rlszp*rlszp+rlszn*rlszn)
    IF(talpinv) THEN                                      !new
      alpeff = alpinv+betinv
      rsinv = 1D0/(1D0+alpeff*rs)
      rsfac = (rs+betinv*(rs*rs))*rsinv*alpeff
      wrk1(iz,ir) =  wrk1(iz,ir) +0.5D0*rsfac*(t3a*(rs*rs)-t3b*(rp*rp+rn*rn))
    ELSE
      rspow   = rs**power
      wrk1(iz,ir) =  wrk1(iz,ir) + rspow/3.0*(b3*rs*rs-b3p*(rp*rp+rn*rn))
    END IF
    IF(tdens2) THEN
      rspow7 = rs**pow7
      potp =  rspow7*(b7a*rs - b7pa*rp -b7pb*(rp*rp+rn*rn)/MAX(rs,xsmall) )  &
          + potp
      potn =  rspow7*(b7a*rs - b7pa*rn -b7pb*(rp*rp+rn*rn)/MAX(rs,xsmall) )  &
          + potn
      wrk1(iz,ir) =  wrk1(iz,ir)  &
          + rspow7*0.5D0*(t7a*(rs*rs)-t7b*(rp*rp+rn*rn))
    END IF                                                  !new
    IF (chg /= 0.0) wrk1(iz,ir) =  0.5*rp*coulpt(iz,ir) + wrk1(iz,ir)
    IF(chgx /= zero) THEN
      wrk1(iz,ir) =  0.75*excfac*MAX(rp,smal)**(4.0d0/3.0d0) + wrk1(iz,ir)
    END IF
    
  END DO
END DO
ekin   = f1int(wrk2)
epot   = f1int(wrk1)

!     2. radial derivative for compensation of the error in
!     "partial integration on the grid".

CALL laplacer_scalar(rho(0:nz,0:nr,1),wrk1, 0, 1, one)
CALL laplacer_scalar(rho(0:nz,0:nr,2),wrk2, 0, 1, one)
ekincr = -0.5*(h2meff(1)*f1int(wrk1)+h2meff(2)*f1int(wrk2))

wrk1 = (rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2))**2
wrk2 = rho(0:nz,0:nr,1)**2+rho(0:nz,0:nr,2)**2

CALL laplacer_scalar(wrk1,wrk3, 0, 1, one)
CALL laplacer_scalar(wrk2,wrk4, 0, 1, one)

epotcr = -0.125*(b1*f1int(wrk3)-b1p*f1int(wrk4))       ! sign ???

!     final composition

epot    = epot+epotcr
ekin    = ekin+ekincr
energy = epot+ekin

IF(itest >= 1) THEN
  WRITE(7,'(a,2(1pg12.5))') '  ekincr,epotcr=',ekincr,epotcr
END IF
DEALLOCATE(wrk1,wrk2,wrk3,wrk4)
RETURN
END SUBROUTINE energ

!---------------------------------------------------------------------------  
! DESCRIPTION: energ2
!> @brief
!!Computes the total energy from sum over s.p. energies and rearrangement 
!!energy.
!--------------------------------------------------------------------------- 

SUBROUTINE energ2()
USE paramr
USE grid
USE dynpar
USE fields
USE levels

IMPLICIT NONE

REAL(8) :: adparp,ekincr,epotr,parfac,exrfac,rm,rp,rs,rspow,sumspe

INTEGER :: ir,ist,iz

REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
wrk1=0D0;wrk2=0D0

IF(tpair) parfac = 1D0/rho0pr         ! factor on pairing potential

!**********************************************************************
!        sum of s.p. energies                                         *
!**********************************************************************  
sumspe = 0.0
DO ist=1,nst(1)+nst(2)
  IF (vsqf(ist) > 1.0D-9) sumspe = spenrg(ist)*wate(ist)*vsqf(ist) + sumspe
END DO

!**********************************************************************
!        two-body potential energy                                    *
!**********************************************************************                                          
wrk1(:,:) = rho(:,:,1)*pot(:,:,1) +rho(:,:,2)*pot(:,:,2)  &
    +tau(:,:,1)*(bmass(:,:,1)-h2meff(1))  &
    +tau(:,:,2)*(bmass(:,:,2)-h2meff(2))

epotr  = oh*f1int(wrk1)         ! ?? note: misses half spin-orbit part

ekin   = h2meff(1)*f1int(tau(0:nz,0:nr,1))+h2meff(2)*f1int(tau(0:nz,0:nr,2))

!**********************************************************************
!     2. radial derivative for compensation of the error in
!     "partial integration on the grid".
!**********************************************************************
CALL laplacer_scalar(rho(0:nz,0:nr,1),wrk1, 0, 1, one)
CALL laplacer_scalar(rho(0:nz,0:nr,2),wrk2, 0, 1, one)
ekincr = -0.5*(h2meff(1)*f1int(wrk1)+h2meff(2)*f1int(wrk2))

ekin   = ekincr+ekin
IF(itest >= 1) WRITE(7,'(a,3(1pg12.5))') ' ekin,ekincor,sumspe=',ekin,ekincr,sumspe
eshell = oh*(sumspe+ekin+equad)

!**********************************************************************
!     rearrangement energy
!**********************************************************************

exrfac = -excfac*oh
DO ir=0,nr
  DO iz=0,nz
    rp      = rho(iz,ir,1)
    rm      = rho(iz,ir,2)
    rs      = rm+rp
    rs      = MAX(rm+rp,smal)
    rspow   = rs**power
    wrk1(iz,ir) = rspow*(t3as*rs*rs-t3bs*(rp*rp+rm*rm))
    IF(chgx /= zero) wrk1(iz,ir) = exrfac*MAX(rp,smal)**(4.0d0/3.0d0) + wrk1(iz,ir)
    IF(tpair) THEN
      adparp = (pairfc(1)*(paircr(iz,ir,1)*paircr(iz,ir,1))  &
          +pairfc(2)*(paircr(iz,ir,2)*paircr(iz,ir,2)))*parfac
      wrk1(iz,ir) = (rs)*adparp  + wrk1(iz,ir)
    END IF
  END DO
END DO
!                                                  finally integrate
erearr = -f1int(wrk1)*oh
epotr  = epotr+erearr



energr = (eshell+erearr)
DEALLOCATE(wrk1,wrk2)
RETURN
END SUBROUTINE energ2

!---------------------------------------------------------------------------  
! DESCRIPTION: obs
!> @brief
!!This routine computes observables from an input density such as particle 
!!number, radius and moments. 
!> @details
!!The properties are evaluated as integrals.
!>
!> @param[in] rhoini
!> REAL(8), array, takes the input density.
!> @param[out] atnum
!> REAL(8), returns the particle number.
!> @param[out] rms
!> REAL(8), returns the rms radius.
!> @param[out] d
!> REAL(8), returns the spherical dipole moment.
!> @param[out] q
!> REAL(8), returns the cartesian quadrupole moment.
!> @param[out] o
!> REAL(8), returns the spherical octupole moment.
!> @param[out] h
!> REAL(8), returns the spherical hexadecupole moment.
!> @param[out] qsq
!> REAL(8), returns the squared quadrupole moment as an estimate of variance.
!> @param[in] atnum
!> REAL(8), takes the center-of-mass relative to which the above moments are
!!computed.
!--------------------------------------------------------------------------- 
SUBROUTINE obs(rhoini,atnum,rms,d,q,o,h,qsq,zcm)

USE paramr
USE grid

IMPLICIT NONE

REAL(8), INTENT(IN)                         :: rhoini(0:nz,0:nr)
REAL(8), INTENT(OUT)                        :: atnum
REAL(8), INTENT(OUT)                        :: rms
REAL(8), INTENT(OUT)                        :: d
REAL(8), INTENT(OUT)                        :: q
REAL(8), INTENT(OUT)                        :: o
REAL(8), INTENT(OUT)                        :: h
REAL(8), INTENT(OUT)                        :: qsq
REAL(8), INTENT(IN)                         :: zcm

REAL(8),PARAMETER :: three=3D0,seven=7D0,eight=8D0,twfour=24D0,drezw=one/12.0D0

REAL(8),ALLOCATABLE,SAVE :: dfield(:,:),qfield(:,:),ofield(:,:),&
                            hfield(:,:),qsqfield(:,:)
INTEGER :: ir,iz

!REAL(8) ::asum,r2sum,zact,ract,rvol,r2

!Prepare helper fields if called for the first time
IF(.NOT.ALLOCATED(rmsfield)) THEN
  ALLOCATE(rmsfield(0:nz,0:nr),dfield(0:nz,0:nr),qfield(0:nz,0:nr),&
           ofield(0:nz,0:nr),hfield(0:nz,0:nr),qsqfield(0:nz,0:nr))
END IF
  rmsfield=0D0
  dfield=0D0
  qfield=0D0
  ofield=0D0
  hfield=0D0
  qsqfield=0D0
  DO ir=0,nr; DO iz=0,nz
    rmsfield(iz,ir)=r(ir)**2+(z(iz)-zcm)**2
    qfield(iz,ir)=2*(z(iz)-zcm)**2-r(ir)**2
    hfield(iz,ir)=eight*(z(iz)-zcm)**4-twfour*(z(iz)-zcm)**2*r(ir)**2+three*r(ir)**4
    qsqfield(iz,ir)=(2*(z(iz)-zcm)**2-r(ir)**2)**2
  END DO; END DO
  DO ir=0,nr; DO iz=0,nz-1
    dfield(iz,ir)=(z(iz)-zcm)
    ofield(iz,ir)=(z(iz)-zcm)*(two*(z(iz)-zcm)**2 - three*r(ir)**2)
  END DO; END DO

!Execute integrals
atnum= f1int(rhoini(0:nz,0:nr))
rms  = f1int(rhoini(0:nz,0:nr)*rmsfield)
rms  = SQRT(ABS(rms/MAX(atnum,smal)))
q    = f1int(rhoini(0:nz,0:nr)*qfield)
h    = f1int(rhoini(0:nz,0:nr)*hfield)* three/(two*eight*SQRT(pi))
qsq  = f1int(rhoini(0:nz,0:nr)*qsqfield)

IF (tsym) THEN
  d = zero
  o = zero
ELSE
  d    = f1int(rhoini(0:nz,0:nr)*dfield)*oh*SQRT(three/pi)
  o    = f1int(rhoini(0:nz,0:nr)*ofield)* oh * oh * SQRT(seven/pi)
END IF

!! check
!asum=0D0
!r2sum=0D0
!DO ir=0,nr
!  ract = dr*ir
!  IF(ir>0) THEN
!    rvol=dr*ir
!  ELSE
!    rvol = dr*drezw
!  END IF
!  DO iz=0,nz
!   zact = dz*iz
!   r2 = ract**2+zact**2
!   IF(iz>0) THEN
!     asum = asum + rvol*rhoini(iz,ir)
!     r2sum = r2sum + r2*rvol*rhoini(iz,ir)
!   ELSE
!     asum = asum + rvol*rhoini(iz,ir)/2D0
!     r2sum = r2sum + r2*rvol*rhoini(iz,ir)/2D0
!   END IF
!  END DO
!END DO
!asum=2D0*3.1415926D0*asum
!r2sum=2D0*3.1415926D0*r2sum
!
!
!WRITE(*,'(a,5(1pg14.6))') 'OBS: nr.,rms(2)=',asum/2D0,SQRT(r2sum/asum),rms
!WRITE(*,'(a,5(1pg14.6))') 'OBS: zcm,d=',zcm,d/atnum

RETURN
END SUBROUTINE obs

!---------------------------------------------------------------------------  
! DESCRIPTION: obsti
!> @brief
!!This routine computes observables from an input density such as particle 
!!number, radius and moments. Constraint fields.
!> @details
!!The properties are evaluated as integrals.
!>
!> @param[in] rhoin
!> REAL(8), array, takes the input density.
!> @param[in] qfac
!> REAL(8), array, takes the actual damped constraint.
!> @param[in] qcfac
!> REAL(8), array, takes the actual raw constraint.
!> @param[out] qa
!> REAL(8), returns the constraining momentum from \c qcfac.
!> @param[out] qc
!> REAL(8), returns the momentum of \c qfac.
!> @param[out] qa2
!> REAL(8), returns the constraining momentum from \c qcfac squared.
!> @param[out] qc2
!> REAL(8), returns the momentum of \c qfac squared.
!--------------------------------------------------------------------------- 
SUBROUTINE obsti(rhoin,qfac,qcfac,qa,qc,qa2,qc2)

USE paramr
USE grid

IMPLICIT NONE

REAL(8), INTENT(IN)   :: rhoin(0:nz,0:nr,kqm)
REAL(8), INTENT(IN)   :: qfac(0:nz,0:nr)
REAL(8), INTENT(IN)   :: qcfac(0:nz,0:nr)
REAL(8), INTENT(OUT)  :: qa
REAL(8), INTENT(OUT)  :: qc
REAL(8), INTENT(OUT)  :: qa2
REAL(8), INTENT(OUT)  :: qc2

qa =f1int((rhoin(:,:,1)+rhoin(:,:,2))*qcfac(:,:))
qc =f1int((rhoin(:,:,1)+rhoin(:,:,2))*qfac(:,:))
qa2=f1int((rhoin(:,:,1)+rhoin(:,:,2))*qcfac(:,:)**2)
qc2=f1int((rhoin(:,:,1)+rhoin(:,:,2))*qfac(:,:)**2)

RETURN
END SUBROUTINE obsti

!---------------------------------------------------------------------------  
! DESCRIPTION: bubble
!> @brief
!!This routine computes the central depth of neutron, proton and total 
!!densities and stores it in \c bubdipn, \c bubdipp and \c bubdipt.
!--------------------------------------------------------------------------- 

SUBROUTINE bubble()
USE grid
USE paramr
USE fields

IMPLICIT NONE

REAL(8) :: rmaxn,rmaxp,rmaxt

rmaxp = MAXVAL(rho(0:nz,0:nr,1))
rmaxn = MAXVAL(rho(0:nz,0:nr,2))
rmaxt = MAXVAL(rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2))

bubdipp = rho(1,1,1)/rmaxp
bubdipn = rho(1,1,2)/rmaxn
bubdipt = (rho(1,1,1)+rho(1,1,2))/rmaxt

RETURN
END SUBROUTINE bubble

!---------------------------------------------------------------------------  
! DESCRIPTION: fragdist
!> @brief
!!Computes the distance of the centers of mass in the right and
!!left half grid for the input density \c rhoin.
!> @details
!!The properties are evaluated as integrals.
!>
!> @param[in] rhoin
!> REAL(8), array, takes the input density.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION fragdist(rhoin)

USE grid
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN) :: rhoin(0:nz,0:nr,kqm)

REAL(8), PARAMETER :: drezw  = one/12.0
REAL(8), PARAMETER :: cutrho = 0.04D0  ! minimal density for search

REAL(8) :: acc0,acc1,dir,rhotot,rhos(3),rhominimum,rhomaximum
REAL(8) :: acc0p,acc0n
REAL(8) :: partnm1,partnm2,partp1,partp2,partn1,partn2,zcm1,zcm2

INTEGER :: ir,iz,izmin,izmax,izminimum

IF(tsym) THEN

! reflection symmetric case
  rhomaximum=-1D3
  acc0  = 0D0
  acc0p  = 0D0
  acc0n  = 0D0
  acc1  = 0D0
  DO iz=1,nz
    rhotot= (rhoin(iz,0,1)+rhoin(iz,0,2))*drezw
    rhomaximum = MAX(rhomaximum,rhotot)
    acc0 = rhotot + acc0
    acc0p = rhoin(iz,0,1)*drezw + acc0p
    acc0n = rhoin(iz,0,2)*drezw + acc0n
    acc1 = rhotot*ABS(z(iz)) + acc1
  END DO

  DO ir=1,nr
    dir   = (ir)+sr
    rhotot= (rhoin(0,ir,1)+rhoin(0,ir,2))*dir*0.5D0
    rhomaximum = MAX(rhomaximum,rhotot)
    acc0p = rhoin(0,ir,1)*dir*0.5D0 + acc0p
    acc0n = rhoin(0,ir,2)*dir*0.5D0 + acc0n
    acc0 = rhotot + acc0
  END DO


  DO ir=1,nr
    dir   = (ir)+sr
    DO iz=1,nz
      rhotot= (rhoin(iz,ir,1)+rhoin(iz,ir,2))*dir
      rhomaximum = MAX(rhomaximum,rhotot)
      acc0p = rhoin(iz,ir,1)*dir + acc0p
      acc0n = rhoin(iz,ir,2)*dir + acc0n
      acc0 = rhotot + acc0
      acc1 = rhotot*ABS(z(iz)) + acc1
    END DO
  END DO

  acc0 = acc0 *(sym+one)*dvol/2D0
  acc0p = acc0p *(sym+one)*dvol/2D0
  acc0n = acc0n *(sym+one)*dvol/2d0
  acc1 = acc1 *(sym+one)*dvol/2D0
  zcm1 = acc1/acc0
  fragdist = 2D0*zcm1

! print and define return value
  WRITE(7,'(/(a,3(1pg13.5)/),2(a,4(1pg13.5)/))') &
      ' fragments:  neck,ratio,distance=',z(0),&
                 (rhoin(0,0,1)+rhoin(0,0,2))/rhomaximum,2D0*zcm1,&
      '            partnum(tot,p,n),zcm=',acc0,acc0p,acc0n,zcm1

ELSE

! no reflection symmetry imposed

! find upper an lower z value for search of minima
  izmin = 0
  DO iz=1,nz
    IF(rhoin(iz,0,1)+rhoin(iz,0,2)>cutrho) EXIT
    izmin = iz    
  END DO
  izmax = nz
  DO iz=nz,1,-1
    IF(rhoin(iz,0,1)+rhoin(iz,0,2)>cutrho) EXIT
    izmax = iz    
  END DO

! search deepest minimum in given interval
  rhominimum=1D3
  rhomaximum=-1D3
  rhos(1)=-2D0
  rhos(2)=-1D0
  izminimum=0
  DO iz=izmin,izmax
    rhos(3)=rhoin(iz,0,1)+rhoin(iz,0,2)
    rhomaximum = MAX(rhomaximum,rhos(3))
    IF(rhos(3)>rhos(2).AND.rhos(1)>rhos(2)) THEN
      IF(rhos(2)<rhominimum) THEN
        izminimum = iz-1
        rhominimum = rhos(2)
      END IF
    END IF
    rhos(1)=rhos(2)
    rhos(2)=rhos(3)
  END DO
!  WRITE(7,*) 'izminimum=',izminimum

! integrate density and c.m. in lower space

  acc0  = 0D0
  acc0p  = 0D0
  acc0n  = 0D0
  acc1  = 0D0
  DO ir=0,nr; DO iz=0,izminimum
    rhotot = (rhoin(iz,ir,1)+rhoin(iz,ir,2))*int_weight(iz,ir)
    IF(iz==izminimum) THEN
      rhotot=rhotot*0.5D0
      acc0p = rhoin(iz,ir,1)*int_weight(iz,ir)*0.5D0 + acc0p
      acc0n = rhoin(iz,ir,2)*int_weight(iz,ir)*0.5D0 + acc0n
    ELSE
      acc0p = rhoin(iz,ir,1)*int_weight(iz,ir) + acc0p
      acc0n = rhoin(iz,ir,2)*int_weight(iz,ir) + acc0n
    END IF
    acc0 = rhotot + acc0
    acc1 = rhotot*z(iz) + acc1
  END DO; END DO

  partnm1 = acc0
  partp1 = acc0p
  partn1 = acc0n
  zcm1 = acc1/acc0

! integrate density and c.m. in upper space

  acc0  = 0D0
  acc0p  = 0D0
  acc0n  = 0D0
  acc1  = 0D0
  DO ir=0,nr; DO iz=izminimum,nz
    rhotot = (rhoin(iz,ir,1)+rhoin(iz,ir,2))*int_weight(iz,ir)
    IF(iz==izminimum) THEN
      rhotot=rhotot*0.5D0
      acc0p = rhoin(iz,ir,1)*int_weight(iz,ir)*0.5D0 + acc0p
      acc0n = rhoin(iz,ir,2)*int_weight(iz,ir)*0.5D0 + acc0n
    ELSE
      acc0p = rhoin(iz,ir,1)*int_weight(iz,ir) + acc0p
      acc0n = rhoin(iz,ir,2)*int_weight(iz,ir) + acc0n
    END IF
    acc0 = rhotot + acc0
    acc1 = rhotot*z(iz) + acc1
  END DO; END DO

  partnm2 = acc0
  partp2 = acc0p
  partn2 = acc0n
  zcm2 = acc1/acc0

! print and define return value
  WRITE(7,'(/(a,3(1pg13.5)/),2(a,4(1pg13.5)/))') &
      ' fragments:  neck,ratio,distance=',z(izminimum),&
                      rhominimum/rhomaximum,zcm2-zcm1,&
      '      left: partnum(tot,p,n),zcm=',partnm1,partp1,partn1,zcm1,&
      '     right: partnum(tot,p,n),zcm=',partnm2,partp2,partn2,zcm2

  fragdist = ABS(zcm2-zcm1)

END IF


RETURN

END FUNCTION fragdist


!---------------------------------------------------------------------------  
! DESCRIPTION: occ_lz
!> @brief
!!This routine computes the total occupation of each \f$ l_z \f$ band and 
!!prints it.
!> @details
!!The properties are evaluated as integrals.
!--------------------------------------------------------------------------- 
SUBROUTINE occ_lz()

USE grid
USE paramr
USE levels
USE dynpar, ONLY: tpri_occup,outnam,fpath
IMPLICIT NONE

REAL(8) :: occuplz(0:kmpsi,2)

INTEGER, SAVE :: lzmax(2)
LOGICAL, SAVE :: tfirst=.true.

INTEGER :: im,ist

IF(tfirst) THEN
  lzmax=0
  DO ist=1,nst(1)+nst(2)
    IF(mpsi(ist)>lzmax(icharg(ist))) lzmax(icharg(ist))=mpsi(ist)
  END DO
  IF(tpri_occup) THEN
    OPEN(171,FILE=TRIM(fpath('occup_lz_p'))//'occup_lz_p.'//outnam)
    WRITE(171,'(a)') '# occupation of l_z bands protons'
    WRITE(171,'(a,900i3)') '# beta2, m=',(im,im=0,lzmax(1))
    OPEN(172,FILE=TRIM(fpath('occup_lz_n'))//'occup_lz_n.'//outnam)
    WRITE(172,'(a)') '# occupation of l_z bands neutrons'
    WRITE(172,'(a,900i3)') '# beta2, m=',(im,im=0,lzmax(2))
  END IF
  tfirst = .false.
END IF

IF(tpri_occup) THEN
  occuplz=0D0
  DO ist=1,nst(1)
    occuplz(mpsi(ist),1) = wate(ist)*vsqf(ist) + occuplz(mpsi(ist),1)
  END DO
  DO ist=nst(1)+1,nst(1)+nst(2)
    occuplz(mpsi(ist),2) = wate(ist)*vsqf(ist) + occuplz(mpsi(ist),2)
  END DO
  WRITE(171,'(f9.3,500f5.1)') beta2,(occuplz(im,1),im=0,lzmax(1))
  WRITE(172,'(f9.3,500f5.1)') beta2,(occuplz(im,2),im=0,lzmax(2))
END IF


RETURN
END SUBROUTINE occ_lz

!---------------------------------------------------------------------------  
! DESCRIPTION: radius
!> @brief
!!Compute formfactor and deduced radius-observables.
!> @details
!!First the grid of the form factors are determined, then the width 
!!parameter for the c.m. corrections are calculated. Next all formfactors 
!!are computed. From those the (charge?) surface width and the rms radii 
!!are determined.
!--------------------------------------------------------------------------- 
SUBROUTINE radius()

USE dynpar
USE levels
USE grid
USE paramr
USE fields, ONLY: drhols

IMPLICIT NONE

REAL(8), PARAMETER :: delrad=0.3D0    ! assumed spherical grid spacing
!REAL(8), PARAMETER :: cutk=4D0,cutdel=0.4D0    ! soft cutoff in k space
REAL(8), PARAMETER :: cutk=3.5D0,cutdel=0.3D0    ! soft cutoff in k space

REAL(8) :: foc(kfrm),fop(kfrm),fon(kfrm),fo(kfrm),foc0(kfrm)
REAL(8) :: rhoback(kfrm),fob(kfrm)   ! auxiliary arrays for back transformation
!REAL(8) :: r4back                    ! fourth radial moment from Fourier back transform
REAL(8) :: cmwid,delk,delrac,q1,qr,rleng,rmax,zmax,zpeact
!REAL(8) :: r2lsp,r2lsn,r2mcor,r2intrnact,rlsmomp2,rlsmomn2,darcor,r2cmcor
REAL(8) :: r2lsp,r2lsn,r2intrnact,rlsmomp2,rlsmomn2,darcor,r2cmcor
REAL(8) :: rmspback,rmsnback,rmsback
INTEGER :: iq,ikmax

!     determine grid in spherical coordinate space from backward F.B. trafo
zmax   = ((nz)+sz-icm)*dz-zcm
rmax   = (nr+sr)*dr
rleng  = MIN(zmax,rmax)
!ikmax  = 2*INT(rleng/delrad)+1
ikmax  = INT(rleng/delrad)+1
delrac = rleng/ikmax
delk   =  pi/rleng  !/2D0
IF(ikmax > kfrm) STOP ' in RADIUS: KFRM too small'

!     width parameter for c.m. correction
IF(ifzpe == 0) THEN
  cmwid  = 1.58D0/(6D0*atnum**0.66666667D0)
ELSE IF(ifzpe == 1) THEN
  zpeact = (widdz+widdxy)*0.5D0*(h2m(1)+h2m(2))/(npart(1)+npart(2))
  cmwid  = 0.375D0*h2mav/(atnum*zpeact)
ELSE
  cmwid  = 0D0
END IF

!     compute all formfactors
!WRITE(7,'(a)') ' iq   q    formfactors '
DO iq=1,ikmax
  q1     = iq*delk-delk
  CALL fourf(q1,cmwid,foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq))
END DO

rdmsc  = 4.493/qzero(foc,1,delk,rleng,ikmax)
rdmsp  = 4.493/qzero(fop,1,delk,rleng,ikmax)
rdmsc0 = 4.493/qzero(foc0,1,delk,rleng,ikmax)

!      evaluating surface width from height of first maximum.
!      the 'q1'="momentum_of_first_maximum" is estimated from
!      'rdmsc' as in a square-well density-distribution.
q1     = 5.600/rdmsc
qr     = q1*rdmsc
surfc  = SQRT(ABS(2./(q1*q1)*LOG(ABS((SIN(qr)/qr-COS(qr))  &
        *3.*apnum/(qr*qr*fbint(q1,foc,rleng,ikmax))))))
surfp  = SQRT(ABS(2./(q1*q1)*LOG(ABS((SIN(qr)/qr-COS(qr))  &
        *3.*apnum/(qr*qr*fbint(q1,fop,rleng,ikmax))))))
surfc0 = SQRT(ABS(2./(q1*q1)*LOG(ABS((SIN(qr)/qr-COS(qr))  &
        *3.*apnum/(qr*qr*fbint(q1,foc0,rleng,ikmax))))))

!WRITE(7,'(a,4(1pg13.5))') ' before formfactors: r_c,r_tot,r_prot,r_neut=',&
!    rmsc,rms,rmsp,rmsn
!dqt    = delk/2.0D0
!dqtd2  = 6.0D0/(dqt*dqt)
CALL rms_back(foc,rmsc,r4c)
IF(tpri_formf) THEN
  WRITE(176,'(a)') '# radial formfactor:',&
      '#    q          F_C           F_tot         F_p            F_n' 
  DO iq=1,ikmax
    q1     = iq*delk-delk
    WRITE(176,'(5(1pg13.5))')  q1,foc(iq),fo(iq),fop(iq),fon(iq)
  END DO
  WRITE(176,'(1x/1x/a)') '# radial charge distribution:',&
      '#  r    rho_charge' 
  DO iq=1,ikmax
    WRITE(176,'(5(1pg13.5))')  delrac*(iq-1),rhoback(iq)
  END DO
END IF
CALL rms_back(fo,rmsback,r4)
CALL rms_back(fop,rmspback,r4p)
CALL rms_back(fon,rmsnback,r4n)
CALL rms_back(foc0,rmsc0,r4c0)

WRITE(7,'(a,4(1pg13.5))') ' from formfactors(incl.c.m.):r_c,r_tot,r_prot,r_neutr=',&
    rmsc,rmsback,rmspback,rmsnback
WRITE(7,'(a,4(1pg13.5))') '                         r4_c,r4_tot,r4_prot,r4_neutr=',&
    r4c,r4,r4p,r4n
WRITE(7,'(a,4(1pg13.5))') ' from charge formfactor:           r,R,sigma,r_unsafe=',  &
    rmsc,rdmsc,surfc,r4c
WRITE(7,'(a,12f8.4)') 'FF:',rmsc,rdmsc,surfc,r4c,rmsc0,rdmsc0,surfc0,r4c0,rmspback,rdmsp,surfp,r4p


! the block to compute the nucleon corrections to rms radii directly
rlsmomp2=f1int(drhols(0:nz,0:nr,1)*rmsfield(0:nz,0:nr))       
rlsmomn2=f1int(drhols(0:nz,0:nr,2)*rmsfield(0:nz,0:nr))       
r2intrnact = (r2intrn*nneut)/nprot
darcor = 3D0*d4m
!       xmgne = -0.04197D0      ! magn.moment * 2 * darcor
!       xmgpr = 0.06131D0
r2lsp = -(xmgpr-isign_rso*d4m)*rlsmomp2/nprot 
r2lsn = -xmgne*rlsmomn2/nprot
r2cmcor = -6D0*cmwid
WRITE(6,'(a,9(1pg13.5))') 'rms_C^2:', &
   rmspback*rmspback+r2intrp+r2intrnact+darcor+r2lsp+r2lsn, &
   rmspback*rmspback-r2cmcor,r2intrp,r2intrnact,darcor,r2lsp,r2lsn,r2cmcor,rmsc**2
WRITE(6,'(a,7(1pg13.5))') '  rms_C:', &
   SQRT(rmspback*rmspback+r2intrp+r2intrnact+darcor+r2lsp+r2lsn),rmsc
!WRITE(7,'(a,9(1pg13.5))') 'rms_C^2:', &
WRITE(7,'(a,9f9.4)') 'rms_C^2:', &
   rmspback*rmspback+r2intrp+r2intrnact+darcor+r2lsp+r2lsn, &
   rmspback*rmspback-r2cmcor,r2intrp,r2intrnact,darcor,r2lsp,r2lsn,r2cmcor,rmsc**2
WRITE(7,'(a,7(1pg13.5))') '  rms_C:', &
   SQRT(rmspback*rmspback+r2intrp+r2intrnact+darcor+r2lsp+r2lsn),rmsc, &
   SQRT(rmspback*rmspback+r2intrp+r2intrnact+darcor)
!   SQRT(6D0*(foc(1)-foc(2))/(delk**2)/foc(1))


RETURN

CONTAINS
SUBROUTINE rms_back(foa,rmsact,r4act)
!     charge density by fourier back-transformation

REAL(8),INTENT(IN) :: foa(kfrm)   ! formfactor to be ttransformed
REAL(8),INTENT(OUT) :: rmsact,r4act  ! radial moments for given formfactors


LOGICAL,PARAMETER :: tprint=.FALSE.

REAL(8) :: dpi3,qa,acc,acc0,acc2,acc4,c,s,dc,ds,x,r,vol,r2
INTEGER :: ir

dpi3   = 0.5D0/(pi*pi)
!                                                         cut high q

DO iq=2,ikmax
  qa     = delk*iq-delk
  fob(iq)=foa(iq)/(1D0+EXP(((qa-cutk)/cutdel)))
END DO
!                                                         r=0.0 separately
IF(tprint) WRITE(78,'(a)') '# test density from back Fourier'
acc    = 0D0
DO iq=2,ikmax
  qa     = delk*iq-delk
  acc    = qa*qa*fob(iq)+acc
END DO
rhoback(1) = ABS(acc*delk*dpi3)
IF(tprint) WRITE(78,'(f6.2,1pg13.5)') 0D0,rhoback(1)
!                                                         now all other r
r      = 0D0
DO ir=2,ikmax
  r      = delrac+r
  
  dc     = COS(r*delk)
  ds     = SIN(r*delk)
  c      = 1D0
  s      = 0D0
  acc    = 0D0
  qa     = 0D0
  DO iq=2,ikmax
    qa     = delk+qa
    x      = c*dc-s*ds
    s      = s*dc+c*ds
    c      = x
    acc    = qa*s*fob(iq)+acc
  END DO
  rhoback(ir) = ABS(acc*delk*dpi3/r)
  IF(tprint) WRITE(78,'(f8.4,1pg13.5)') r,rhoback(ir)
END DO

! r.m.s. radius from backward trafo
acc0 = 0D0
acc2 = 0D0
acc4 = 0D0
r = 0D0
DO ir=2,ikmax
  r = delrac+r
  r2 = r*r
  vol = r2*rhoback(ir)
  acc0 = acc0 + vol
  acc2 = acc2 + vol*r2
  acc4 = acc4 + vol*r2*r2
END DO

rmsact = SQRT(acc2/acc0)
r4act = (acc4/acc0)**0.25D0

RETURN

END SUBROUTINE rms_back

END SUBROUTINE radius

!---------------------------------------------------------------------------  
! DESCRIPTION: fourf
!> @brief
!!Computes the spherical Fourier-Bessel transform of the proton and neutron 
!!densities. The formfactors thus obtained are augmented by the nucleon 
!!formfactor.

!> @details
!!The present version omits the contribution from the magnetic formfactor 
!!of the nucleon.
!>
!> @param[in] q1
!> REAL(8), takes the radial momentum at which the fourier-transformation
!! is to be evaluated.
!> @param[inout] cmwid
!> REAL(8), parameter for center of mass correction.
!> @param[out] foc0
!> REAL(8), returns the charge formfactor.
!> @param[out] fop
!> REAL(8), returns the proton formfactor.
!> @param[out] fon
!> REAL(8), returns the neutron formfactor.
!> @param[out] fo
!> REAL(8), returns the total formfactor.
!> @param[out] foc
!> REAL(8), returns the charge formfactor including spin-orbit.
!--------------------------------------------------------------------------- 
SUBROUTINE fourf(q1,cmwid,foc0,fop,fon,fo,foc)

USE fields
USE grid
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)                         :: q1
REAL(8), INTENT(IN OUT)                     :: cmwid
REAL(8), INTENT(OUT)                        :: foc0
REAL(8), INTENT(OUT)                        :: fop
REAL(8), INTENT(OUT)                        :: fon
REAL(8), INTENT(OUT)                        :: fo
REAL(8), INTENT(OUT)                        :: foc

REAL(8), PARAMETER ::  half=0.5D0

REAL(8) :: besfac,darfac,exphf,q2,qa,rad2,schsen,schsep,schses,sumc,&
           sumn,sump,schsev
INTEGER :: ir,iz

!      the parameters in data are for the sachs-formfactors:
!       sa..,sm.. are for the isoscalar electric formfactor
!       va..,vm.. are for the isovector electric formfactor
!       gp..,gm.. are for the magnetic formfactors of proton and neutron

!       isoscalar and isovector formfactors are recoupled to
!       the electric formfactors of proton and neutron.

!     the sachs-formfactors are then recoupled to the nucleon
!     vector formfactor (form1..) and the nucleon tensor
!     formfactor (form2..).

!     parameters for scalar and vector electric formfactor
!     (walther, priv.comm., sept. 86)

REAL(8),PARAMETER :: sa1=2.2907,sa2=-0.6777,sa3=-0.7923,sa4=0.1793,&
                     sm1=15.75, sm2=26.68,  sm3=41.01,  sm4=134.2
REAL(8),PARAMETER :: va1=0.3681,va2=1.2263,va3=-0.6316,va4=0.0372,&
                     vm1=5.00,  vm2=15.02,  vm3=44.08,  vm4=154.2

!     magnetic formfactor , simon et al n.p. a333 (1980) 381

REAL(8),PARAMETER :: gp1=0.694,gp2=0.719,gp3=-0.418,gp4=0.005,&
                     gm1=8.50, gm2=15.02,gm3=44.08, gm4=355.4

!     the anomanlous magnetic moments and the darwin factor 'd4m'
!REAL(8),PARAMETER :: xmgne=0.04197d0,xmgpr=-0.06131d0,d4m=0.010989d0  as in header, bit other sign

LOGICAL,PARAMETER :: ifrso=.TRUE.   ! include contribution from spin-orbit density

!REAL(8),PARAMETER :: apro0=0.74d0,apro11=0.0d0,apro12=0.0d0,bpro11=0.0d0
!REAL(8),PARAMETER :: bneu0=0.74d0,aneu11=0.0d0,aneu12=0.0d0,bneu11=0.0d0

REAL(8) :: sumsp,sumsn,schsm,schsmp,schsmn
REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
!     fourier-bessel transformation for proton- and neutron-density
wrk1=0D0;wrk2=0D0

q2 = q1*q1
DO ir=0,nr
  DO iz=0,nz
    rad2   = (z(iz)**2+r(ir)**2)*q2
    IF(rad2 < smal) THEN
      besfac = 1.0D0
    ELSE
      rad2 = SQRT(rad2)
      besfac = SIN(rad2)/rad2
    END IF
    wrk1(iz,ir) = rho(iz,ir,1)*besfac
    wrk2(iz,ir) = rho(iz,ir,2)*besfac
  END DO
END DO
sump  = f1int(wrk1)
sumn  = f1int(wrk2)

IF(ifrso) THEN
  DO ir=0,nr
    DO iz=0,nz
      rad2   = (z(iz)**2+r(ir)**2)*q2
      IF(rad2 < smal) THEN
        besfac = 1.0D0
      ELSE
        rad2 = SQRT(rad2)
        besfac = SIN(rad2)/rad2
      END IF
      wrk1(iz,ir) = drhols(iz,ir,1)*besfac
      wrk2(iz,ir) = drhols(iz,ir,2)*besfac
    END DO
  END DO
  sumsp  = f1int(wrk1)
  sumsn  = f1int(wrk2)
ELSE
  sumsp  = 0D0
  sumsn  = 0D0
END IF

!      isoscalar and isovector electrical formactors on 'schses'
!      and 'schsev', recoupled then to proton and neutron electric
!      formfactors  'schsep' and 'schsen'. proton and neutron
!      magnetic formfactor on 'schsmp' and 'schsmn'.

qa     = q1*q1
schses = sa1/(1D0+qa/sm1)+sa2/(1D0+qa/sm2)+sa3/(1D0+qa/sm3) +sa4/(1D0+qa/sm4)
schsev = va1/(1D0+qa/vm1)+va2/(1D0+qa/vm2)+va3/(1D0+qa/vm3) +va4/(1D0+qa/vm4)
schsep = half*(schses+schsev)
schsen = half*(schses-schsev)
sumc   = sump*schsep+sumn*schsen

schsm  = ( gp1/(one+qa/gm1)+gp2/(one+qa/gm2)+gp3/(one+qa/gm3)+gp4/(one+qa/gm4) )
schsmp = -schsm*xmgpr + isign_rso*schsep*d4m   ! +  for correct scheme
schsmn = -schsm*xmgne + isign_rso*schsen*d4m


!      factor for darwin correction 'darfac'
!      and the factor for the c.m. correction 'exphf'.

qa     = q1*q1
if(ifrso) then
  darfac = 1D0/sqrt(1D0+qa*d4m)
else
  darfac = 1D0
endif
exphf  = EXP(cmwid*qa)

!     finally assembling the formfactors

foc0   = exphf*sumc
fop    = exphf*sump
fon    = exphf*sumn
fo     = exphf*(sump+sumn)
foc    = exphf*(sumc*darfac + (schsmp*sumsp+schsmn*sumsn))

DEALLOCATE(wrk1,wrk2)
RETURN
END SUBROUTINE fourf

!---------------------------------------------------------------------------  
! DESCRIPTION: fbint
!> @brief
!!This routine interpolates the formfactor from the formfactor on a discrete 
!!grid at momentum \c qeff = \c q * \c rleng. 
!> @details
!!The fourier-bessel interpolation runs into problems with round-off
!!errors if \c qeff comes too close to a multiple of \f$ \pi \f$. We switch
!!to an explicit treatment of the pole-term in that case. \c precis defines 
!!the neighbourhood of a pole-term for which the routine switches to an explicit 
!!laurent expansion.
!>
!> @param[in] qin
!> REAL(8), takes the q value.
!> @param[in] formfa
!> REAL(8), array, takes the formfactor on the grid.
!> @param[in] delk
!> REAL(8), takes the grid spacing.
!> @param[in] rleng
!> REAL(8), scaling factor for effective q value.
!> @param[in] ikmax
!> REAL(8), grid size of \c formfa.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION fbint(qin,formfa,rleng,ikmax)

USE grid
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)         :: qin
REAL(8), INTENT(IN)         :: formfa(kfrm)
REAL(8), INTENT(IN)         :: rleng
INTEGER, INTENT(IN)         :: ikmax

REAL(8),PARAMETER :: precis=1D-9

REAL(8) :: qdif,qeff,signum,sinq,sum,sum2,xfx,xfx2
INTEGER :: i,ipole,ixfx

!-----------------------------------------------------------------------

IF(ikmax > kfrm) STOP ' too large grid in FBINT '

qeff   = qin*rleng
xfx    = (qeff/pi)
xfx2   = xfx*xfx
ixfx   = INT(xfx+precis)

IF(ABS(xfx-ixfx) < precis) THEN
!                                             pole term separately
  ipole  = ixfx+1
  signum = -one
  sum    =  zero
  DO i=2,ikmax
    signum = -signum
    IF(i /= ipole) sum  = formfa(i)*signum/(one-xfx2/((i-1)*(i-1))) + sum
  END DO
  signum = (MOD(ipole,2)*2-1)
  qdif   = (xfx-ixfx)
  sum2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
  qdif   = qdif*pi
  sinq   = (one-qdif*qdif*0.166666667D0)*signum
  fbint  = 2.0D0*(sum*qdif+sum2*pi)*sinq/qeff
ELSE
!                                             normal case
  signum = -one
  sum    = zero
  DO i=2,ikmax
    signum = -signum
    sum    = formfa(i)*signum/(one-xfx2/FLOAT((i-1)*(i-1))) + sum
  END DO
  fbint  = (sum+sum)*SIN(qeff)/qeff
END IF

RETURN
END FUNCTION fbint

!---------------------------------------------------------------------------  
! DESCRIPTION: qzero
!> @brief
!!Evaluates and returns the momentum at which the formfactor given on array 
!!\c foin has a zero.
!> @details
!!Scan for zero starts at the field position \c iqstrt. \c ikmax is the 
!!actual length of the array \c foin.\c delk is the stepsize in the momentum.
!>
!> @param[in] foin
!> REAL(8), array, takes the form factor.
!> @param[in] iqstart
!> REAL(8), takes the field position for start of scanning.
!> @param[in] delk
!> REAL(8), takes the grid spacing of \c foin.
!> @param[in] rleng
!> REAL(8), scaling factor for effective q value.
!> @param[in] ikmax
!> REAL(8), grid size of \c foin.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION qzero(foin,iqstrt,delk,rleng,ikmax)

USE grid
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)             :: foin(kfrm)
INTEGER, INTENT(IN)             :: iqstrt
REAL(8), INTENT(IN)             :: delk
REAL(8), INTENT(IN)             :: rleng
INTEGER, INTENT(IN)             :: ikmax

INTEGER,PARAMETER :: itq0mx=3
REAL(8),PARAMETER :: endq0=1D-4

REAL(8) :: c2,c3,ca,f1,f2,f3,qa1,qa2,q3old,q1,q2,q3,c11
INTEGER :: iq1,iq2,iscan,itq0,ikmaxm

!     scan changing sign in 'foin'.
ikmaxm = ikmax-1
DO  iscan=iqstrt,ikmaxm
  IF(foin(iscan)*foin(iscan+1) <= zero) EXIT
END DO

iq1   = iscan
iq2   = iq1+1
IF(foin(iq1) <= zero) iq2   = iq1-1
q1    = iq1*delk-delk
q2    = iq2*delk-delk
f1    = foin(iq1)
f2    = foin(iq2+1) !GH (+1)

!     come closer to 'foin'=0.0 by repeated linear, interpolation.
q3     = q1-f1*(q2-q1)/(f2-f1)
DO  itq0=1,itq0mx
  f3     = fbint(q3,foin,rleng,ikmax)
  IF(f1*f3 > zero) THEN
    q1     = q3
    f1     = f3
  ELSE
    q2     = q3
    f2     = f3
  END IF
  q3old  = q3
  q3     = q1-f1*(q2-q1)/(f2-f1)
  IF(ABS(q3-q3old) < endq0) EXIT
END DO

!     improve position of the zero by quadratic interpolation:

!                         the coefficients of the interpolating binomial
f1     = fbint(q1,foin,rleng,ikmax)
f2     = fbint(q2,foin,rleng,ikmax)
f3     = fbint(q3,foin,rleng,ikmax)
ca     = one/((q1-q2)*(q2-q3)*(q3-q1))
c11     = ca*(q3*(q3-q2)*q2*f1+q1*(q1-q3)*q3*f2+q2*(q2-q1)*q1*f3)
c2     = ca*((q2*q2-q3*q3)*f1+(q3*q3-q1*q1)*f2+(q1*q1-q2*q2)*f3)
c3     = ((q3-q2)*f1+(q1-q3)*f2+(q2-q1)*f3)*ca

!                         the two zeroes of the interpolating binomial
qa1    = ( SQRT(ABS(c2*c2-4.0D0*c11*c3))-c2)/(c3+c3)
qa2    = (-SQRT(ABS(c2*c2-4.0D0*c11*c3))-c2)/(c3+c3)

!                         select the zero which is closer to q3
q      = qa1
IF(ABS(qa1-q3) > ABS(qa2-q3)) q=qa2

qzero  = q

RETURN
END FUNCTION qzero

!---------------------------------------------------------------------------  
! DESCRIPTION: zpe_cm
!> @brief
!!Zero point energy for c.m. motion. as <P**2>/2mA
!--------------------------------------------------------------------------- 
SUBROUTINE zpe_cm

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid

IMPLICIT NONE

INTEGER,PARAMETER :: intest=0,itstpr=0

REAL(8),ALLOCATABLE :: psi1(:,:,:),psin(:,:,:),psip(:,:,:)

REAL(8) :: deltai,deltaj,edifi,edifj,eferm,equasi,equasj,pij,&
          sumdir,tstovl,ui,uj,vi,vj
INTEGER :: ir,is,ist,jst,jstmax,jstmin,mist,mm,nparr,npist

ALLOCATE(psi1(0:nz,0:nr,nsd),psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
psi1=0D0;psip=0D0;psin=0D0

!     1. double loop to accumulate x-y-width:
IF(intest == 1) WRITE(7,'(a)') ' Nabla_xy:'
sumdir = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
!     gradient operator Nabla+ on 'psip' and Nabla- on 'psin'
  CALL nablar_spin(psi(0:nz,0:nr,1:2,ist),psi1(0:nz,0:nr,1:2),mpsi(ist),1,one)
  
  DO is=1,nsd
    mm        = mist+is-1
    npist  = nparz(is,ist)
    nparr     = 1+MOD(mm,2)
    ir = 0
    IF(mm == 1) THEN
      psip(0:nz,ir,is) = 0D0
      psin(0:nz,ir,is) = 2D0*psi1(0:nz,ir,is)
    ELSE IF(mm == -1) then
      psin(0:nz,ir,is) = 0D0
      psip(0:nz,ir,is) = 2D0*psi1(0:nz,ir,is)
    ELSE
      psip(0:nz,ir,is) = 0D0 !psi1(0:nz,ir,is)
      psin(0:nz,ir,is) = 0D0 !psip(0:nz,ir,is)
    END IF
    DO ir=1,nr
      psip(0:nz,ir,is) = psi1(0:nz,ir,is)-mm/r(ir)*psi(0:nz,ir,is,ist)
      psin(0:nz,ir,is) = psi1(0:nz,ir,is)+mm/r(ir)*psi(0:nz,ir,is,ist)
    END DO
  END DO
!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  sumdir = 2.0*vi*vi*(tstovl+wfovlp(psin,psin))*wtpair_v2(ist) +sumdir
  IF(intest == 1) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <Nabla**2>=',ist,tstovl
  DO jst=jstmin,jstmax
    IF( ( npist == nparz(1,jst) .AND.  &
          (mist+1 == mpsi(jst) .OR. mist-1 == mpsi(jst)) )  &
        .OR.  &
        ( npist == nparz(2,jst) .AND.  &
          ( (mist == 0 .AND. mpsi(jst) == 0)  &
           .OR.(mist.eq.-1 .and. mpsi(jst).eq.-1)))      ) then
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      IF(isTsym.LE.0 .AND. mist == 0 .AND. mpsi(jst) == 0) THEN
        pij    = wfovoc(psi(0:nz,0:nr,1,jst),psin(0:nz,0:nr,2))  &
            + wfovoc(psi(0:nz,0:nr,2,jst),psin(0:nz,0:nr,1))
      ELSE IF(mist+1 == mpsi(jst)) THEN
        pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      ELSE IF(mist-1 == mpsi(jst)) THEN
        pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psin)
      ELSE
        pij = 0D0
      END IF
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
      IF(intest == 1) WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij, pij*(vj*ui-uj*vi)
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
!               pij    = pij*(vj*ui-uj*vi)
!               pij    = pij*pij
!               sumexy = pij/(equasi+equasj) + sumexy
!               sumwxy  = pij + sumwxy
    END IF
  END DO
END DO
IF(isTsym>0) sumdir=sumdir/2D0
widdxy = sumdir/2.0

!     2. loop to accumulate c.m. width along z
IF(intest == 1) WRITE(7,'(a)') ' Nabla_z:'
sumdir = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
!     gradient operator Nabla_z on 'psip'
  
  CALL nablaz_spin(psi(0:nz,0:nr,1:2,ist),psi1,nparz(1,ist),1,one)
  psip = psi1
!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  IF(intest /= 0) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <Delta_z**2>=',ist,tstovl
  sumdir = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir
  DO jst=jstmin,jstmax
    IF(npist == nparz(2,jst) .AND. mist == mpsi(jst)) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      IF(intest == 1)   WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij*pij,  &
          (pij*(vj*ui-uj*vi))**2
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
    END IF
  END DO
END DO
IF(isTsym>0) sumdir=sumdir/2D0
widdz  = sumdir


!anumb  = npart(1)+npart(2)
!IF(ifzpe > 0) THEN
!  ezp    = 17.3/anumb**0.2
!ELSE
!  ezp    = 0.0D0
!END IF

IF(ifzpe == 1) THEN
  ezp = (widdz+widdxy)*0.5D0*(h2m(1)+h2m(2))/(npart(1)+npart(2)) 
ELSE
  ezp = 0.0D0
END IF
IF(itest>=0) WRITE(7,'(a,3(1pg12.4))') ' total ZPE: ',ezp,widdz,widdxy
WRITE(7,'(a,3(1pg12.4))') ' total ZPE: ',ezp,widdz,widdxy

DEALLOCATE(psi1,psin,psip)

RETURN
END SUBROUTINE zpe_cm

!---------------------------------------------------------------------------  
! DESCRIPTION: localize
!> @brief
!!This routine computes the localization criterion. 
!> @details
!!First the routine writes the headers for density and localization files 
!!which can be viewed with plotmtv. Then the localization measure is 
!!constructed. Finally it is written to the plotmtv file.
!>
!> @param[in] outnam
!> CHARACTER, array, takes output name suffix.
!--------------------------------------------------------------------------- 
SUBROUTINE localize()

USE paramr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid

IMPLICIT NONE

LOGICAL, PARAMETER :: ttest=.false.

INTEGER, SAVE     :: numcall=0
CHARACTER (LEN=3) :: tnum
REAL(8)           :: chom,clocal,ctilde,curvt,fpair,rhoact,ttffac,&
                     xm2,xm2l,zpoint
INTEGER           :: is,iq,ir,ist,istlow,istup,iwunit,iz,nparr,nparrc,mm

REAL(8),ALLOCATABLE :: psin(:,:,:),psip(:,:,:)
REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:),wrk3(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr),wrk3(0:nz,0:nr))
ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
wrk1=0D0;wrk2=0D0;wrk3=0D0;psip=0D0;psin=0D0


numcall = 1 + numcall
OPEN(78,STATUS='scratch')
WRITE(78,'(i4)') 1000+numcall
REWIND(78)
READ(78,'(1x,a)') tnum
CLOSE(78)
iwunit = 78

!     outer loop over isospin

DO iq=1,2
  IF(iq == 1) THEN
    istlow = 1
    istup  = nst(1)
  ELSE
    istlow = nst(1)+1
    istup  = nst(1)+nst(2)
  END IF
  
!     output for 'plotmtv' of density
  
  IF(iq == 1) THEN
    OPEN(iwunit,FILE=TRIM(fpath('densip'))//'densip.'//tnum//'.'//outnam)
    WRITE(iwunit,'(a)') '# plot density for protons:'
  ELSE
    OPEN(iwunit,FILE=TRIM(fpath('densin'))//'densin.'//tnum//'.'//outnam)
    WRITE(iwunit,'(a)') '# plot density for neutrons:'
  END IF
  WRITE(iwunit,'(a,2i4,2f8.3)') '# nprot,nneut,beta2,beta3=',npart,beta2,beta3
  WRITE(iwunit,'(a/a,i3,a,f7.2/a,i3,a,i3/a,i3,a,f7.2/a,i3)')  &
      '$ data=contour', '% nx=',nz+1,' xmin= 0.0 xmax=',nz*dz,  &
      '% vxmin=',INT(nz*dz*(0.25-0.25*isym)),  &
      ' vxmax=',INT(nz*dz*(0.75-0.25*isym)),  &
      '% ny=',nr+1,' ymin= 0.0 ymax=',nr*dr, '% vymax=',INT(nr*dr/2.0)
  WRITE(iwunit,'(a)') '% contstyle=2',  &
      '% nsteps=50', '% interp=2',  &
      '% xlabel="z [fm]"', '% ylabel="r [fm]"'
  IF(iq == 1) THEN
    WRITE(iwunit,'(a,2i4,a/a,2f6.3,a)')  &
        '% toplabel ="proton density for Z,N=',npart,'" '
  ELSE
    WRITE(iwunit,'(a,2i4,a/a,2f6.3,a)')  &
        '% toplabel ="neutron density for Z,N=',npart,'" '
  END IF
  WRITE(iwunit,'(a,2f6.3,a)')  &
      '% subtitle ="deformation beta_2,beta_3=',beta2,beta3,'"'
  WRITE(iwunit,'(a)') '% grid=true',  &
      '% equalscale=true', '% fitpage=false',  &
      '% clip=true'
  WRITE(iwunit,'(a,f6.3)') '% xyratio=',1.0*nr/nz
  WRITE(iwunit,'(a)') '% cmin=0',  &
      '% cmax=0.1', '% cstep=0.01   '
  DO ir=0,nr
    WRITE(iwunit,'(10f8.5)') (rho(iz,ir,iq),iz=0,nz)
    WRITE(iwunit,'(a)') '     '
  END DO
  CLOSE(iwunit)
  
!     prepare header for 'plotmtv' file
  
  IF(iq == 1) THEN
    OPEN(iwunit,FILE=TRIM(fpath('localp'))//'localp.'//tnum//'.'//outnam)
    WRITE(iwunit,'(a)') '# plot localization for protons:'
  ELSE
    OPEN(iwunit,FILE=TRIM(fpath('localn'))//'localn.'//tnum//'.'//outnam)
    WRITE(iwunit,'(a)') '# plot localization for neutrons:'
  END IF
  WRITE(iwunit,'(a,2i4,2f8.3)') '# nprot,nneut,beta2,beta3=',npart,beta2,beta3
  WRITE(iwunit,'(a/a,i3,a,f7.2/a,i3,a,i3/a,i3,a,f7.2/a,i3)')  &
      '$ data=contour', '% nx=',nz+1,' xmin= 0.0 xmax=',nz*dz,  &
      '% vxmin=',INT(nz*dz*(0.25-0.25*isym)),  &
      ' vxmax=',INT(nz*dz*(0.75-0.25*isym)),  &
      '% ny=',nr+1,' ymin= 0.0 ymax=',nr*dr, '% vymax=',INT(nr*dr/2.0)
  WRITE(iwunit,'(a)') '% contstyle=2',  &
      '% nsteps=50', '% interp=2',  &
      '% xlabel="z [fm]"', '% ylabel="r [fm]"'
  IF(iq == 1) THEN
    WRITE(iwunit,'(a,2i4,a/a,2f6.3,a)')  &
        '% toplabel ="proton localization for Z,N=',npart,'" '
  ELSE
    WRITE(iwunit,'(a,2i4,a/a,2f6.3,a)')  &
        '% toplabel ="neutron localization for Z,N=',npart,'" '
  END IF
  WRITE(iwunit,'(a,2f6.3,a)')  &
      '% subtitle ="deformation beta_2,beta_3=',beta2,beta3,'"'
  WRITE(iwunit,'(a)') '% grid=true',  &
      '% equalscale=true', '% fitpage=false',  &
      '% clip=true'
  WRITE(iwunit,'(a,f6.3)') '% xyratio=',1.0*nr/nz
  WRITE(iwunit,'(a)') '% cmin=0',  &
      '% cmax=1', '% cstep=0.1   '
  
!     compute derivative of density consistently
  
  wrk1 = 0.0D0
  wrk2 = 0.0D0
  rho(0:nz,0:nr,iq) = 0.0D0
  tau(0:nz,0:nr,iq) = 0.0D0
  DO ist=istlow,istup
    fpair  = vsqf(ist)*wate(ist)
    mm     = mpsi(ist)
    xm2    = mm*mm
    xm2l   = (mm+1)*(mm+1)
    nparr  = 1+MOD(mm,2)
    nparrc = 3-nparr
    IF(ttest) WRITE(iwunit,'(a,3i5,1pg13.5)')  &
        ' localize: iq,ist,mm,fpair=',iq,ist,mm,fpair
    CALL nablaz_spin(psi(0:nz,0:nr,1:2,ist),psin(0:nz,0:nr,1:2),nparz(1,ist),1,one)
    CALL nablar_spin(psi(0:nz,0:nr,1:2,ist),psip(0:nz,0:nr,1:2),   mpsi(ist),1,one)
    CALL densit(psi(0:nz,0:nr,1:2,ist),vsqf(ist),ist,rho(0:nz,0:nr,iq),  &
        tau(0:nz,0:nr,iq),wrk3(0:nz,0:nr),wrk3(0:nz,0:nr),wrk3(0:nz,0:nr))
    DO is=1,2
      wrk1(0:nz,0:nr) = psip(0:nz,0:nr,is)*psi(0:nz,0:nr,is,ist)*fpair + wrk1(0:nz,0:nr)
      wrk2(0:nz,0:nr) = psin(0:nz,0:nr,is)*psi(0:nz,0:nr,is,ist)*fpair + wrk2(0:nz,0:nr)
    END DO
  END DO
  
!     interpolate to standard grid points and sum up
  
  IF(ttest) WRITE(6,'(a)')   'localization:',  &
      '  iz  ir    z     r      c      c_tilde   rho'
  WRITE(iwunit,'(a,i3/a)') '# for nucleon type iq=',iq,  &
      '#  iz  ir    z     r      c      c_tilde   rho'
  ttffac = 0.6*(3.0*pi**2)**0.6666666667
  DO ir=0,nr
    DO iz=0,nz
      rhoact = rho(iz,ir,iq)
      IF(rhoact > smal) THEN
        clocal = tau(iz,ir,iq) -(wrk1(iz,ir)**2+wrk2(iz,ir)**2)/rhoact
        chom   = ttffac*rhoact**1.666666667
        ctilde = 1.0D0/(1.0D0+(clocal/chom)**2)
      ELSE
        clocal = 1.0D-20
        chom   = 1.0D-10
        curvt  = 0.0D0
        ctilde = 0.0D0
      END IF
      zpoint  = ((iz)+sz-icm)*dz
      IF(ttest) WRITE(6,'(1x,2i5,2f8.3,3(1pg13.5))')  &
          iz,ir,zpoint,(ir-0.5D0)*dr,clocal,ctilde
      WRITE(iwunit,'(f8.5)') ctilde
    END DO
    WRITE(iwunit,'(5x,a)') ' '
  END DO
  CLOSE(78)
END DO
DEALLOCATE(psin,psip,wrk1,wrk2,wrk3)
RETURN
END SUBROUTINE localize

!---------------------------------------------------------------------------  
! DESCRIPTION: densit
!> @brief
!!Accumulates the densities \c rhod, \c taud, \c divlsd, calculated for the 
!!state in wave function \c psiw with the occupation probability \c weight.
!> @details
!!This state is supposed to be related to the global label \c ist.
!>
!> @param[in] psiw
!> REAL(8), array, takes the wave function of state \ist.
!> @param[in] weight
!> REAL(8), takes the occupation probability of the state \c ist (does not 
!!include the multiplicity \c wate.
!> @param[in] ist
!> INTEGER, takes the absolute adress of the state.
!> @param[out] rhod
!> REAL(8), array, takes and returns the density.
!> @param[out] taud
!> REAL(8), array, takes and returns the kinetic density.
!> @param[out] divlsd
!> REAL(8), array, takes and returns the divergence ls-density.
!> @param[out] rlsr
!> REAL(8), array, takes and returns the ls-density r-component.
!> @param[out] rlsz
!> REAL(8), array, takes and returns the ls-density z-component.
!--------------------------------------------------------------------------- 

SUBROUTINE densit(psiw,weight,ist,rhod,taud,divlsd,rlsr,rlsz)

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)                         :: psiw(0:nz,0:nr,nsd)
REAL(8), INTENT(IN)                         :: weight
INTEGER, INTENT(IN)                         :: ist
REAL(8), INTENT(INOUT)                      :: rhod(0:nz,0:nr)
REAL(8), INTENT(INOUT)                      :: taud(0:nz,0:nr)
REAL(8), INTENT(INOUT)                      :: divlsd(0:nz,0:nr)
REAL(8), INTENT(INOUT)                      :: rlsr(0:nz,0:nr)
REAL(8), INTENT(INOUT)                      :: rlsz(0:nz,0:nr)

REAL(8) :: fpair,xm2,xm2l,xmlow,xmup,zentfl,zentfu
INTEGER :: ir,mm,nparr,nparrc

REAL(8),ALLOCATABLE :: psip(:,:,:),psin(:,:,:)
REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
wrk1=0D0;wrk2=0D0;psip=0D0;psin=0D0
fpair  = weight*wate(ist)

!---- some combinations of <j_z>

mm     = mpsi(ist)
xm2    = mm*mm
xm2l   = (mm+1)*(mm+1)
nparr  = 1+MOD(mm,2)
nparrc = 3-nparr

!---- calculate derivatives of psiw

CALL nablaz_spin(psiw,psin,nparz(1,ist),1,one)
CALL nablar_spin(psiw,psip,   mpsi(ist),1,one)

!---- calculate rhod and taud --------------------------------------------

!---- value at 'r=0'

IF(mm == 0) THEN
  rhod(0:nz,0) = fpair * psiw(0:nz,0,1)*psiw(0:nz,0,1) + rhod(0:nz,0)
  wrk1(0:nz,0) = fpair * ( psin(0:nz,0,1)**2 + 2.0*psip(0:nz,0,2)**2 )    ! save that for spin-orbit
  taud(0:nz,0) = wrk1(0:nz,0)  + taud(0:nz,0)
ELSE IF(mm == 1) THEN
  wrk1(0:nz,0) = (fpair+fpair) * ( psip(0:nz,0,1)**2 )               ! save that for spin-orbit
  taud(0:nz,0) = wrk1(0:nz,0) + taud(0:nz,0)
!NEW ODD
ELSE IF(mm == (-1)) THEN
    rhod(0:nz,0) = fpair * psiw(0:nz,0,2)*psiw(0:nz,0,2) + rhod(0:nz,0)
    wrk1(0:nz,0) = fpair * ( psin(0:nz,0,2)*psin(0:nz,0,2)  &
        + 2.0*psip(0:nz,0,1)*psip(0:nz,0,1) )              ! save that for spin-orbit               
    taud(0:nz,0) = wrk1(0:nz,0)  + taud(0:nz,0)
ELSE IF(mm == (-2)) THEN
    wrk1(0:nz,0) = (fpair+fpair) * ( psip(0:nz,0,2)*psip(0:nz,0,2) )           ! save that for spin-orbit                
    taud(0:nz,0) = wrk1(0:nz,0) + taud(0:nz,0)
!END NEW ODD
END IF

!---- all other 'r'

DO ir=1,nr
  zentfu = xm2/(r(ir)**2)
  zentfl = xm2l/(r(ir)**2)
  rhod(0:nz,ir) = fpair*(psiw(0:nz,ir,1)**2+psiw(0:nz,ir,2)**2) + rhod(0:nz,ir)
  wrk1(0:nz,ir) = fpair*(psin(0:nz,ir,1)*psin(0:nz,ir,1)  &
      +psin(0:nz,ir,2)*psin(0:nz,ir,2) +psip(0:nz,ir,1)*psip(0:nz,ir,1)  &
      +psip(0:nz,ir,2)*psip(0:nz,ir,2)&
      +zentfu*psiw(0:nz,ir,1)**2+zentfl*psiw(0:nz,ir,2)**2) ! save that for spin-orbit
  taud(0:nz,ir) = wrk1(0:nz,ir) + taud(0:nz,ir)
END DO

!----- calculate spin-orbit density ----------------------------------

!     uses divlsd = sum_ist |sigma*nabla psiw|^2-|nabla psiw|^2

IF(t4 /= zero) THEN
  ir=0                        !  test???
  IF(mm == 0) THEN
    wrk2(0:nz,0) = -(fpair*(psin(0:nz,0,1)+ 2.0*psip(0:nz,0,2))**2-wrk1(0:nz,0))
    divlsd(0:nz,0) = wrk2(0:nz,0) + divlsd(0:nz,0)
    rlsr(0:nz,0) = fpair*psiw(0:nz,ir,1)*psin(0:nz,ir,2) +rlsr(0:nz,0)
    rlsz(0:nz,0) = -fpair*psiw(0:nz,ir,1)*2.0*psip(0:nz,ir,2) + rlsz(0:nz,0)
  ELSE IF(mm == 1) THEN
    wrk2(0:nz,0) = wrk1(0:nz,0)
    divlsd(0:nz,0) = wrk1(0:nz,0) + divlsd(0:nz,0)
  ELSE
    wrk2(0:nz,0) = 0.0D0
  END IF
  DO ir=1,nr
    xmup  = mm/(r(ir))
    xmlow = (mm+1)/(r(ir))
    wrk2(0:nz,ir)= -(fpair*( &
           (psin(0:nz,ir,1)+psip(0:nz,ir,2)+xmlow*psiw(0:nz,ir,2))**2 &
          +(-psin(0:nz,ir,2)+psip(0:nz,ir,1)-xmup*psiw(0:nz,ir,1))**2)&
          -wrk1(0:nz,ir))
    divlsd(0:nz,ir) = wrk2(0:nz,ir) + divlsd(0:nz,ir)
    rlsr(0:nz,ir) = fpair*(  &
           psiw(0:nz,ir,1)*(xmup*psiw(0:nz,ir,1)+psin(0:nz,ir,2))  &
          -psiw(0:nz,ir,2)*(psin(0:nz,ir,1)+xmlow*psiw(0:nz,ir,2))) +rlsr(0:nz,ir)
    rlsz(0:nz,ir) = fpair*( psiw(0:nz,ir,2)*  &
           (psip(0:nz,ir,1)-xmup*psiw(0:nz,ir,1) ) - psiw(0:nz,ir,1)*  &
           (psip(0:nz,ir,2)+xmlow*psiw(0:nz,ir,2) ) ) + rlsz(0:nz,ir)
  END DO
END IF
DEALLOCATE(psip,psin,wrk1,wrk2)

RETURN
END SUBROUTINE densit

!NEW ODD
SUBROUTINE densitodd(psiw,weight,ist,  &
    rrotj,zrotj,rotsig,rsig,zsig,slarsig,slazsig, curr,spdenkinr,spdenkinz)

! DESCRIPTION: densitodd
!> @brief
!!Accumulates the time-odd densities \c rotj, \c curr, \c sigma, calculated for the 
!!state in wave function \c psiw with the occupation probability \c weight.


!**********************************************************************
!                                                                     *
!        accumulates the time-odd densities
!                                                                     *
!**********************************************************************

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)     :: psiw(0:nz,0:nr,nsd)
REAL(8), INTENT(IN)     :: weight
INTEGER, INTENT(IN)     :: ist
REAL(8), INTENT(IN OUT)    :: rrotj(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: zrotj(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: rotsig(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: rsig(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: zsig(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: slarsig(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: slazsig(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: curr(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: spdenkinr(0:nz,0:nr)
REAL(8), INTENT(IN OUT)    :: spdenkinz(0:nz,0:nr)

REAL(8),ALLOCATABLE :: psip(:,:,:),psin(:,:,:)
REAL(8),ALLOCATABLE :: psipp(:,:,:),psinn(:,:,:)

REAL(8) :: fpair,fpair2,rinv,xmlow,fup,flow,xm2,xm2l,xmup
REAL(8) :: sigism,f,add
INTEGER :: iz,ir,mm,nparrc,nparr,ism

!-----------------------------------------------------------

ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
ALLOCATE(psipp(0:nz,0:nr,nsd),psinn(0:nz,0:nr,nsd))


!  in case of calculating one the half spectrum, set
!  time-odd densitys to zero

!WRITE(6,*) 'CUR in SUBR:',ist,SUM(curr**2)
IF(toddterms)THEN
  fpair  = weight*wate(ist)
ELSE
  fpair = 0D0
END IF


!---- some combinations of <j_z>

mm     = mpsi(ist)
xm2    = mm*mm
xm2l   = (mm+1)*(mm+1)
nparr  = 1+MOD(mm,2)
nparrc = 3-nparr

!---- calculate derivatives of psiw

CALL nablaz_spin(psiw,psin,nparz(1,ist),1,1D0)
CALL nablar_spin(psiw,psip,   mpsi(ist),1,1D0)
CALL nablaz_spin(psin,psinn,nparz(1,ist),1,1D0)
CALL nablar_spin(psip,psipp,   mpsi(ist),1,1D0)

!**********************************************************************
!       rrotj r-comp of rot j
!       zrotj z-comp of rot j
!**********************************************************************

!----- r = 0

!      fpair = (-1.0)*fpair
! only z-comp changes
IF(mm == 0) THEN
  DO iz=0,nz
    zrotj(iz,0)=zrotj(iz,0) + 2*fpair*psip(iz,0,2)*psip(iz,0,2)
  END DO
ELSE IF(mm == 1) THEN
  DO iz=0,nz
    zrotj(iz,0)=zrotj(iz,0) + 2*fpair*psip(iz,0,1)*psip(iz,0,1)
  END DO
ELSE IF(mm == (-1)) THEN
  DO iz=0,nz
    zrotj(iz,0)=zrotj(iz,0) - 2*fpair*psip(iz,0,1)*psip(iz,0,1)
  END DO
ELSE IF(mm == (-2)) THEN
  DO iz=0,nz
    zrotj(iz,0)=zrotj(iz,0) - 2*fpair*psip(iz,0,2)*psip(iz,0,2)
  END DO
END IF

!---- all other 'r'

DO ir=1,nr
  rinv  = 1D0/((sr+ir)*dr)
  xmup  = 2*mm*rinv
  xmlow = 2*( mm+1)*rinv
  DO iz=0,nz
    zrotj(iz,ir)=zrotj(iz,ir) + fpair*(xmup*psiw(iz,ir,1)*psip(iz,ir,1)+  &
        xmlow*psiw(iz,ir,2)*psip(iz,ir,2))
    rrotj(iz,ir)=rrotj(iz,ir) - fpair*(xmup*psiw(iz,ir,1)*psin(iz,ir,1)+  &
        xmlow*psiw(iz,ir,2)*psin(iz,ir,2))
  END DO
END DO

!      fpair = (-1.0)*fpair

!**********************************************************************
!                Current j
!  curr = phi-comp. Current
!**********************************************************************
!IF(ist==80 .OR. ist==23) THEN
!WRITE(6,*) SUM(curr**2)
!---- r=0
! nothing has to be done
DO iz=0,nz
  curr(iz,0)=zero
END DO

!---- all other r
DO ir=1,nr
  rinv  = 1D0/((sr+ir)*dr)
  fup  = mm*rinv
  flow = (mm+1D0)*rinv
  DO iz=0,nz
    add = fpair*(fup*psiw(iz,ir,1)*psiw(iz,ir,1) &
                 + flow*psiw(iz,ir,2)*psiw(iz,ir,2))
    curr(iz,ir)= curr(iz,ir)+add
  END DO
END DO
!WRITE(6,*) 'CURR:',ist,mm,itpartn(ist),mpsi(itpartn(ist)),mm+1,fpair,SUM(curr**2),SUM(psiw**2)
!END IF

!**********************************************************************

!              sigma
!             rsig = r-comp of sigma
!             zsig = z-comp of sigma

!**********************************************************************
fpair2= (1D0)*fpair/1D0
DO ir=0,nr
  DO iz=0,nz
    rsig(iz,ir)=rsig(iz,ir)+fpair2*2*psiw(iz,ir,1)*psiw(iz,ir,2)
    zsig(iz,ir)=zsig(iz,ir)+fpair2*(  &
        psiw(iz,ir,1)*psiw(iz,ir,1) - psiw(iz,ir,2)*psiw(iz,ir,2))
  END DO
END DO




!**********************************************************************

!             rot sigma
!             rotsig = phi-comp of rot sigma

!**********************************************************************

DO ir=0,nr
  DO iz=0,nz
    rotsig(iz,ir)= rotsig(iz,ir) + 2*fpair*(psiw(iz,ir,1)*  &
        psin(iz,ir,2) + psiw(iz,ir,2)*psin(iz,ir,1)  &
        -psiw(iz,ir,1)*psip(iz,ir,1) + psiw(iz,ir,2)*psip(iz,ir,2))
  END DO
END DO

!**********************************************************************

!             Laplace sigma
!             larsig = r-comp of Laplace sigma
!             lazsig = z-comp of Laplace sigma

!**********************************************************************

!----- r = 0


! only z-comp changes
ir=0
IF((mm == 0).OR.(mm == (-2))) THEN
  DO iz=0,nz
    slazsig(iz,ir)=slazsig(iz,ir) - 4*fpair*psip(iz,ir,2)*psip(iz,ir,2)
  END DO
END IF
IF((mm == 1).OR.(mm == (-1))) THEN
  DO iz=0,nz
    slazsig(iz,ir)=slazsig(iz,ir) + 4*fpair*psip(iz,ir,1)*psip(iz,ir,1)
  END DO
END IF

! now all other r
DO ir=1,nr
  rinv  = 1D0/((sr+ir)*dr)
  DO iz=0,nz
    slarsig(iz,ir)=slarsig(iz,ir)+fpair*(2*  &
        (psiw(iz,ir,1)*psipp(iz,ir,2) +psiw(iz,ir,2)*psipp(iz,ir,1)  &
        + 2*psip(iz,ir,1)*psip(iz,ir,2)) + rinv*2*  &
        (psip(iz,ir,1)*psiw(iz,ir,2) + psip(iz,ir,2)*psiw(iz,ir,1))  &
        - rinv*rinv*2*psiw(iz,ir,1)*psiw(iz,ir,2)  &
        + 2*(psiw(iz,ir,1)*psinn(iz,ir,2) + psiw(iz,ir,2)*psinn(iz,ir,1)  &
        + 2*psin(iz,ir,1)*psin(iz,ir,2)))
    
    slazsig(iz,ir)=slazsig(iz,ir) + 2*fpair*(  &
        psip(iz,ir,1)*psip(iz,ir,1) -psip(iz,ir,2)*psip(iz,ir,2)  &
        +psiw(iz,ir,1)*psipp(iz,ir,1)-psiw(iz,ir,2)*psipp(iz,ir,2)  &
        +rinv*(psiw(iz,ir,1)*psip(iz,ir,1)-psiw(iz,ir,2)*psip(iz,ir,2))  &
        +psinn(iz,ir,1)*psiw(iz,ir,1)-psinn(iz,ir,2)*psiw(iz,ir,2)  &
        +psin(iz,ir,1)*psin(iz,ir,1)-psin(iz,ir,2)*psin(iz,ir,2))
    
  END DO
END DO


!*********************************************************************
!                kinetic spin density
!********************************************************************

! r = 0
IF((mm == 0).OR.(mm == (-1)))  THEN
  ism = ABS(mm)+1
  sigism = 3 -2*ism
  DO iz=0,nz
! perhaps spdenkinr = 0
    spdenkinz(iz,0)=spdenkinz(iz,0)+sigism*fpair*(  &
        psin(iz,0,ism)*psin(iz,0,ism))
  END DO
END IF
IF((mm == 0).OR.(mm == (-2))) THEN
  DO iz=0,nz
    spdenkinz(iz,0)=spdenkinz(iz,0)-2*fpair* psin(iz,0,2)*psin(iz,0,2)
  END DO
ELSE IF((mm == 1).OR.(mm == (-1))) THEN
  DO iz=0,nz
    spdenkinz(iz,0)=spdenkinz(iz,0)+2*fpair* psin(iz,0,1)*psin(iz,0,1)
  END DO
END IF

! now all other r
DO ir=1,nr
  rinv  = 1D0/((sr+ir)*dr)
  f = mm*(mm+1)*(rinv*rinv)
  fup = mm*rinv
  flow = (mm+1)*rinv
  DO iz=0,nz
    spdenkinr(iz,ir)= spdenkinr(iz,ir)+2*fpair*(psip(iz,ir,1)*  &
        psip(iz,ir,2) + f*psiw(iz,ir,1)*psiw(iz,ir,2)+  &
        psin(iz,ir,2)*psin(iz,ir,1))
    spdenkinz(iz,ir)=spdenkinz(iz,ir)+fpair*(psip(iz,ir,1)*  &
        psip(iz,ir,1)+fup*fup*psiw(iz,ir,1)*psiw(iz,ir,1)  &
        + psin(iz,ir,1)*psin(iz,ir,1) - (  &
        psip(iz,ir,2)*psip(iz,ir,2)+flow*flow*psiw(iz,ir,2)*  &
        psiw(iz,ir,2) + psin(iz,ir,2)*psin(iz,ir,2) ) )
    
  END DO
END DO

DEALLOCATE(psip,psin,psipp,psinn)


RETURN
END SUBROUTINE densitodd
!END NEW ODD


!---------------------------------------------------------------------------  
! DESCRIPTION: sdenti
!> @brief
!!Accumulates the density \c rho.
!> @details
!!This state is supposed to be related to the global label \c ist. The fields 
!!have to be multiplied with (1.0- \c addrho) before the m-loop in the 
!!calling routine. This routine is used in connection with the iteration of 
!!quadrupole deformations. 
!>
!> @param[in] psiw
!> REAL(8), array, takes the wave function of state \ist.
!> @param[out] rhod
!> REAL(8), array, takes and returns the density.
!> @param[in] addwgt
!> REAL(8), takes the occupation probability of the state \c ist (does not 
!!include the multiplicity \c wate.
!> @param[in] ist
!> INTEGER, takes the absolute adress of the state.
!--------------------------------------------------------------------------- 

SUBROUTINE sdenti(psiw,rhoacc,addwgt,ist)

USE paramr
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)                         :: psiw(0:nz,0:nr,nsd)
REAL(8), INTENT(INOUT)                      :: rhoacc(0:nz,0:nr)
REAL(8), INTENT(IN)                         :: addwgt
INTEGER, INTENT(IN)                         :: ist

REAL(8),PARAMETER :: lowadd=1.0D-9

REAL(8) :: fpair

IF (addwgt < lowadd) RETURN
fpair  = wate(ist)*addwgt
rhoacc(0:nz,0:nr) = fpair * (psiw(:,:,1)**2+psiw(:,:,2)**2) + rhoacc(0:nz,0:nr)

RETURN
END SUBROUTINE sdenti


!NEW ODD
SUBROUTINE cpairpot(psiw1,psiw2,rhoacc,addwgt,ist)

!---------------------------------------------------------------------------  
! DESCRIPTION: cpairpot
!> @brief
!!Accumulates pairing density in case of broken Cramer' degeneracy

USE levels, ONLY: wate
USE grid, ONLY: nz,nr
USE wavefr, ONLY: wfovoc,wfovpc

IMPLICIT NONE

REAL(8), INTENT(IN)      :: psiw1(0:nz,0:nr,2)
REAL(8), INTENT(IN)      :: psiw2(0:nz,0:nr,2)
REAL(8), INTENT(OUT)     :: rhoacc(0:nz,0:nr)
REAL(8), INTENT(IN)      :: addwgt
INTEGER, INTENT(IN)      :: ist

REAL(8), PARAMETER :: lowadd=1D-9

IF(ABS(addwgt)<lowadd) RETURN

rhoacc(0:nz,0:nr) = wate(ist)*addwgt * &
      (psiw1(0:nz,0:nr,1)*psiw2(0:nz,0:nr,2)  &
     -psiw1(0:nz,0:nr,2)*psiw2(0:nz,0:nr,1)) + rhoacc(0:nz,0:nr)


RETURN
END SUBROUTINE cpairpot
!END NEW ODD


SUBROUTINE testconj_all(commentin)

!---------------------------------------------------------------------------  
! DESCRIPTION: cpairpot
!> @brief
!!Accumulates pairing density in case of broken Cramer' degeneracy

USE levels, ONLY: nst,mpsi,itpartn
USE grid, ONLY: nz,nr
USE wavefr, ONLY: psi


IMPLICIT NONE

CHARACTER(*),INTENT(IN) :: commentin

INTEGER :: ist
REAL(8) :: overconj,kval1,kval2

DO ist=1,nst(1)+nst(2)
  overconj = testconj(psi(:,:,:,ist),psi(:,:,:,itpartn(ist)))
  kval1 = mpsi(ist)+0.5D0
  kval2 = mpsi(itpartn(ist))+0.5D0
  IF(kval1*overconj<0D0) THEN
    WRITE(6,'(a,2i5,2f5.1,1pg13.5)') commentin//', sign CONJ repaired:',&
      ist,itpartn(ist),kval1,kval2,overconj
    psi(:,:,:,itpartn(ist)) = -psi(:,:,:,itpartn(ist))
  END IF
END DO

RETURN
END SUBROUTINE testconj_all

REAL(8) FUNCTION testconj(psiw1,psiw2)

!---------------------------------------------------------------------------  
! DESCRIPTION: testconj
!> @brief
!!Tests the overlap between conjugates states

USE grid, ONLY: nz,nr,f1int

IMPLICIT NONE

REAL(8), INTENT(IN) :: psiw1(0:nz,0:nr,2)
REAL(8), INTENT(IN) :: psiw2(0:nz,0:nr,2)

REAL(8), ALLOCATABLE :: rho1(:,:)

ALLOCATE(rho1(0:nz,0:nr))

  rho1(0:nz,0:nr) =  (psiw1(0:nz,0:nr,1)*psiw2(0:nz,0:nr,2)  &
                      -psiw1(0:nz,0:nr,2)*psiw2(0:nz,0:nr,1)) 
  testconj = f1int(rho1)

DEALLOCATE(rho1)

RETURN
END FUNCTION testconj

!---------------------------------------------------------------------------  
! DESCRIPTION: zcmobs
!> @brief
!!Computes the center-of-mass for the density distribution given in two 
!!components on \c rhoini.
!>
!> @param[in] psiin
!> REAL(8), array, takes the input density.
!--------------------------------------------------------------------------- 

FUNCTION zcmobs(rhoini)

USE paramr
USE dynpar
USE fields
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)    :: rhoini(0:nz,0:nr,kqm)

REAL(8) :: zcmobs
REAL(8) :: atnum,d
REAL(8) :: zfield(0:nz,0:nr)
INTEGER :: iz

!----------------------------------------------------------------------

!     compute actual center of mass

IF(.NOT.tsym) THEN
  atnum = f1int(rhoini(0:nz,0:nr,1)+rhoini(0:nz,0:nr,2))
  zfield=0.0d0
  DO iz=0,nz-1
    zfield(iz,:) = (rhoini(iz,0:nr,1)+rhoini(iz,0:nr,2))*z(iz)
  END DO
  d=f1int(zfield)
  zcmobs = d/atnum
ELSE
  zcmobs = zero
END IF

RETURN
END FUNCTION zcmobs

FUNCTION check_border_density()

USE paramr
USE grid
USE fields

IMPLICIT NONE
LOGICAL :: check_border_density
LOGICAL :: mask(0:nz,0:nr)
REAL(8) :: rhot(0:nz,0:nr)

mask=.FALSE.
mask(nz,0:nr)=.TRUE.
mask(0:nz,nr)=.TRUE.
rhot=rho(:,:,1)+rho(:,:,2)
check_border_density = (MAXVAL(rhot,MASK=mask)>1D-6)

END FUNCTION check_border_density

!---------------------------------------------------------------------------  
! DESCRIPTION: spradius
!> @brief
!!This routine computes the single-particle radius for the given
!!s.p. wavefunction. It serves for testing purposes. Printing is
!!done inside the routine.
!> @details
!!The properties are evaluated as integrals.
!>
!> @param[in] psiin
!> REAL(8), array, takes the input spinor wavefunction.
!!computed.
!--------------------------------------------------------------------------- 
!SUBROUTINE spradius(psiin)
REAL(8) FUNCTION spradius(psiin)

USE paramr
USE grid

IMPLICIT NONE

REAL(8), INTENT(IN)                         :: psiin(0:nz,0:nr,1:2)

REAL(8),ALLOCATABLE :: auxrho(:,:),auxz(:,:),auxr2(:,:)

INTEGER :: ir,iz,is
REAL(8) :: spnum,zcm,rms2

!Prepare helper fields if called for the first time
ALLOCATE(auxrho(0:nz,0:nr),auxz(0:nz,0:nr),auxr2(0:nz,0:nr))
 

!DO is=1,2
  auxrho(:,:) = psiin(:,:,1)**2+psiin(:,:,2)**2
!  auxrho(:,:) = psiin(:,:,is)
  DO ir=0,nr; DO iz=0,nz
    auxr2(iz,ir)=auxrho(iz,ir)*(r(ir)**2+z(iz)**2)
    auxz(iz,ir)=auxrho(iz,ir)*z(iz)
  END DO; END DO


!Execute integrals
  spnum= f1int(auxrho(0:nz,0:nr))
  IF(tsym) THEN
    zcm=0D0
  ELSE
    zcm  = f1int(auxz(0:nz,0:nr))/spnum
  END IF
  rms2 = f1int(auxr2(0:nz,0:nr))/spnum
  rms2 = rms2-zcm**2

!  WRITE(6,'(a,i3,4(1pg13.5))') 'SPRADIUS: is,norm,zcm,rms=',is,spnum,zcm,SQRT(rms2),rms2
!END DO
spradius = SQRT(rms2)

DEALLOCATE(auxrho,auxz,auxr2)
RETURN
END FUNCTION spradius
!END SUBROUTINE spradius

!------------calangularmomentum-------------

!  

!---------------------------------------------------------------------------  
! DESCRIPTION: angularmom
!> @brief
!!This routine computes the total angular momentum of the
!! s.p. wavefunction with label 'ist'.
!--------------------------------------------------------------------------- 
!SUBROUTINE calangularmomentum(ist,angularmom)
REAL(8) FUNCTION angularmom(ist)

USE paramr
USE grid, ONLY: nz,nr,sr,sz,icm,dr,dz,nablaz_spin,nablar_spin
USE wavefr
USE levels

IMPLICIT NONE

INTEGER, INTENT(IN) :: ist
!REAL(8), INTENT(OUT) :: angularmom
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'
!INCLUDE 'derivs.inc'

REAL(8),ALLOCATABLE :: psi1(:,:,:),psi2(:,:,:),psi3(:,:,:),psi4(:,:,:) !,psip(:,:,:)
REAL(8),ALLOCATABLE :: psin(:,:,:),psip(:,:,:)

INTEGER :: iz,ir
REAL(8) :: r,z,rinv,amom2


ALLOCATE(psi1(0:nz,0:nr,2),psi2(0:nz,0:nr,2),psi3(0:nz,0:nr,2),psi4(0:nz,0:nr,2))
ALLOCATE(psin(0:nz,0:nr,2),psip(0:nz,0:nr,2))

! calculate derivatives
CALL nablaz_spin(psi(:,:,:,ist),psin,nparz(1,ist),1,1D0)
CALL nablar_spin(psi(:,:,:,ist),psip,  mpsi(ist),1,1D0)

DO iz=0,nz
  DO ir=0,nr
    r = (ir+sr)*dr
    z = (iz+sz-icm)*dz
    IF(ir == 0)THEN
      rinv = zero
    ELSE
      rinv = 1.0/r
    END IF
! psi1 = L_r|psi>
    
    IF(ir > 0)THEN
      psi1(iz,ir,1)= rinv*z*mpsi(ist)*psi(iz,ir,1,ist)
      psi1(iz,ir,2)= rinv*z*(mpsi(ist)+1)*psi(iz,ir,2,ist)
    ELSE
      IF(mpsi(ist) == 1) THEN
        psi1(iz,ir,1)=z*psip(iz,ir,1)
        psi1(iz,ir,2)=0
      ELSE IF(mpsi(ist) == 0) THEN
        psi1(iz,ir,1)=0
        psi1(iz,ir,2)=z*psip(iz,ir,2)
      ELSE
        psi1(iz,ir,1)=0
        psi1(iz,ir,2)=0
      END IF
    END IF
! psi2 = L_phi|psi>
    psi2(iz,ir,1)=r*psin(iz,ir,1) - z*psip(iz,ir,1)
    psi2(iz,ir,2)=r*psin(iz,ir,2) - z*psip(iz,ir,2)
    
! psi3 = L_z|psi>
    psi3(iz,ir,1)=mpsi(ist)*psi(iz,ir,1,ist)
    psi3(iz,ir,2)=(mpsi(ist)+1)*psi(iz,ir,2,ist)
    
! psi4 = (L_r*s_r +L_phi*s_phi+L_z*s_z)|psi>
    psi4(iz,ir,1)= + r*psin(iz,ir,2) - z*psip(iz,ir,2)  &
        + mpsi(ist)*psi(iz,ir,1,ist)
    psi4(iz,ir,2)= - r*psin(iz,ir,1) + z*psip(iz,ir,1)  &
        - (mpsi(ist)+1)*psi(iz,ir,2,ist)
    IF(ir > 0) THEN
      psi4(iz,ir,1)=psi4(iz,ir,1) -z*rinv*psi(iz,ir,2,ist)*(mpsi(ist)+1)
      psi4(iz,ir,2)=psi4(iz,ir,2) -z*rinv*psi(iz,ir,1,ist)*(mpsi(ist))
    ELSE
      IF(mpsi(ist) == 0) THEN
        psi4(iz,ir,1)=psi4(iz,ir,1)-z*psip(iz,ir,2)
      ELSE IF(mpsi(ist) == 1)THEN
        psi4(iz,ir,2)=psi4(iz,ir,2)-z*psip(iz,ir,1)
      END IF
    END IF
  END DO
END DO
amom2 = wfovlp(psi1,psi1)+wfovlp(psi2,psi2) +wfovlp(psi3,psi3)  &
    +wfovlp(psi4,psi(:,:,:,ist))+0.75D0
DEALLOCATE(psi1,psi2,psi3,psi4,psin,psip)


angularmom = 0.5D0*(-1.0+ SQRT(1.0+4.0*(amom2)))


RETURN
END FUNCTION angularmom
!END SUBROUTINE calangularmomentum

!---------------------------------------------------------------------------  
! DESCRIPTION: Qvariance
!> @brief
!! Variance of quadrupole momentum, presently only for istsym>0
!--------------------------------------------------------------------------- 
SUBROUTINE Qvariance(DelQ2)

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
USE fields, ONLY: rho

IMPLICIT NONE

REAL(8), INTENT(OUT) :: DelQ2(-2:2)

INTEGER,PARAMETER :: itstpr=0
REAL(8),PARAMETER :: fac20=SQRT(5D0/(16D0*pi))
REAL(8),PARAMETER :: fac21=SQRT(15D0/(8D0*pi))
REAL(8),PARAMETER :: fac22=SQRT(15D0/(32D0*pi))
LOGICAL,PARAMETER :: tntest=.FALSE.

REAL(8),ALLOCATABLE :: psi1(:,:,:),psin(:,:,:),psip(:,:,:)

REAL(8) :: deltai,deltaj,edifi,edifj,eferm,equasi,equasj,pij,&
          sumdir,tstovl,ui,uj,vi,vj,r2,sumdir0
INTEGER :: iz,ir,is,ist,jst,jstmax,jstmin,mist,mm,nparr,npist


IF(istsym.LE.0) THEN
  DelQ2=0D0
END IF

ALLOCATE(psi1(0:nz,0:nr,nsd),psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
psi1=0D0;psip=0D0;psin=0D0

!     1. double loop to accumulate Q_20 width
IF(tntest) WRITE(7,'(a)') ' Q_20:'
sumdir = zero
sumdir0 = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
! apply Q_20
  DO ir=0,nr
    r2     = r(ir)**2 !(dr*(ir+sr))**2
    DO iz=0,nz
!      z     = ((iz)+sz-icm)*dz
      psip(iz,ir,:) = fac20*(2D0*z(iz)**2-r2)*psi(iz,ir,:,ist)
    END DO
  END DO
  
!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  sumdir = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir  
  sumdir0 =  2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir0  
  IF(itstpr /= 0) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <sp|Q**2|sp>=',ist,tstovl
  DO jst=jstmin,jstmax
    IF(  npist == nparz(1,jst) .AND. mist == mpsi(jst) ) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      IF(tntest .AND. vi*vj>0.5D0) &
        WRITE(7,'(a,3i4,3(1pg13.5))') 'pij=',mist,ist,jst,vi,vj,pij**2
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
      IF(itstpr /= 0) WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij, pij*(vj*ui-uj*vi)
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
!               pij    = pij*(vj*ui-uj*vi)
!               pij    = pij*pij
!               sumexy = pij/(equasi+equasj) + sumexy
!               sumwxy  = pij + sumwxy
    END IF
  END DO
END DO
!IF(isTsym>0) sumdir=sumdir/2D0
DelQ2(0) = sumdir/2.0
IF(tntest) WRITE(7,'(a,1pg13.5)') 'diag.sum 20:',sumdir0/2D0


!     2. loop to accumulate Q_21-variance
IF(tntest) WRITE(7,'(a)') ' Q_21:'
sumdir = zero
sumdir0 = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
!  action of Q_21 operator
  DO ir=0,nr
!    r     = (dr*(ir+sr))
    DO iz=0,nz
!      z     = ((iz)+sz-icm)*dz
      psip(iz,ir,:) = fac21*z(iz)*r(ir)*psi(iz,ir,:,ist)
    END DO
  END DO

!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  IF(itstpr /= 0) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <Delta_zr**2>=',ist,tstovl
  sumdir = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir
  sumdir0 = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir0
  DO jst=jstmin,jstmax
    IF(npist == nparz(2,jst) .AND. mist+1 == mpsi(jst)) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      IF(itstpr /= 0)   WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij*pij,  &
          (pij*(vj*ui-uj*vi))**2
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
    END IF
  END DO
END DO
!IF(isTsym>0) sumdir=sumdir/2D0
DelQ2(1) = sumdir/2.0
DelQ2(-1) = sumdir/2.0
IF(tntest) WRITE(7,'(a,1pg13.5)') 'diag.sum 21:',sumdir0/2D0



!     3. loop to accumulate Q_22-variance
IF(tntest) WRITE(7,'(a)') ' Q_22:'
sumdir = zero
sumdir0 = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
!  action of Q_21 operator
  DO ir=0,nr
    r2     = r(ir)**2 !(dr*(ir+sr))**2
    DO iz=0,nz
!      z     = ((iz)+sz-icm)*dz
      psip(iz,ir,:) = fac22*r2*psi(iz,ir,:,ist)
    END DO
  END DO

!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  IF(itstpr /= 0) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <Delta_zr**2>=',ist,tstovl
  sumdir = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir
  sumdir0 = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir0
  DO jst=jstmin,jstmax
    IF(npist == nparz(1,jst) .AND. mist+2 == mpsi(jst)) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      IF(itstpr /= 0)   WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij*pij,  &
          (pij*(vj*ui-uj*vi))**2
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
    END IF
  END DO
END DO
!IF(isTsym>0) sumdir=sumdir/2D0
DelQ2(2) = sumdir/2.0
DelQ2(-2) = sumdir/2.0
IF(tntest) WRITE(7,'(a,1pg13.5)') 'diag.sum 22:',sumdir0/2D0


! test versus r**4
IF(tntest) THEN
  DO ir=0,nr
    r2     = (dr*(ir+sr))**2
    DO iz=0,nz
      psip(iz,ir,1) = (z(iz)**2+r2)**2*(rho(iz,ir,1)+rho(iz,ir,2))
    END DO
  END DO
  WRITE(7,'(a,3(1pg13.5))') '<r^4>:',f1int(psip(:,:,1))/(4D0*pi) &
     ,f1int(rho(:,:,1)),f1int(rho(:,:,2))
END IF

DEALLOCATE(psi1,psin,psip)

RETURN
END SUBROUTINE Qvariance


END MODULE observ
