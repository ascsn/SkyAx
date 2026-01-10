!------------------------------------------------------------------------------
! MODULE: pairing
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routines to solve the pairing problem.
!------------------------------------------------------------------------------
MODULE pairin

USE dynpar, ONLY: pair_cutoff,tcutoff_WS,tcutoff_uvspace

IMPLICIT NONE


REAL(8), ALLOCATABLE :: wcoltt(:,:)   !<Workspace for Coulomb solver.

REAL(8), ALLOCATABLE :: rspln(:)      !<Workspace for Lipkin Nogami.
REAL(8), ALLOCATABLE :: suv2rs(:)     !<Workspace for Lipkin Nogami.

INTEGER :: ilnrsp=1       !<Activates mean-field contribution to Lipkin-Nogami

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: hfbpot
!> @brief
!! Computes the pair density, potential, energy and gaps.
!> @details
!!The pairing density is given by
!! \f[ \xi_q(\vec r) = \sum_{\alpha\in q} \sum_{s} u_{\alpha} v_{\alpha} 
!! \psi_{\overline{\alpha}}(\vec{r},s)\psi_{\alpha}(\vec{r},s) \f]
!! the pairing potential by
!! \f[ U_{\rm pair}(\vec r)=-V_\mathrm{pair,q}
!! \xi_q(\vec r) \left[1 -\frac{\rho}{\rho_{0,\mathrm{pair}}}\right] \f]
!! the pairing energy by
!! \f[ E_{\rm pair}=\frac{1}{4} \sum_{q\in\{p,n\}}V_\mathrm{pair,q}
!! \int \D^3r |\xi_q(\vec r)|^2 \left[1 -\frac{\rho(\vec r)}{\rho_{0,\mathrm{pair}}}\right]\f]
!! and the pairing gap by
!!\f[ \Delta_k=\int\D^3r V_P(\vec r)\left|\phi_k(\vec r)\right|^2\f]
!! \c ipair is used here to indicate the switch to surface pairing: 5,6: 
!! volume pairing; 7,8: surface pairing.
!--------------------------------------------------------------------------- 
SUBROUTINE hfbpot(iter)

USE paramr
USE dynpar
USE fields
USE observ, ONLY: sdenti,spradius,cpairpot
USE wavefr
USE levels
USE grid
USE fyfunc, ONLY: fy_hfbpot
!USE steputil, ONLY: cpairpot

IMPLICIT NONE

INTEGER, INTENT(IN) :: iter

REAL(8) :: smallp = 0.000001d0

REAL(8) ::  ecrit,npower
REAL(8) ::  epairtestp,epairtestn,epairtest
REAL(8) ::  buff(2),geff,rhofac

REAL(8), ALLOCATABLE :: work(:,:)
REAL(8), ALLOCATABLE :: paircrtest(:,:,:),pairpttest(:,:,:)

INTEGER :: i,iq,jst,iz,ir
REAL(8) :: uv,v0prot,v0neut,v2
REAL(8) :: densint0,densintfold,densintdiff,epairtestnfold
LOGICAL, PARAMETER :: testfolding=.FALSE.

ALLOCATE(work(0:nz,0:nr))
IF(testfolding) ALLOCATE(paircrtest(0:nz,0:nr,kqm),pairpttest(0:nz,0:nr,kqm))
work=0D0

!     zero field before accumulation
paircr(:,:,:) = 0D0
IF(testfolding) paircrtest(:,:,:) = 0D0

! reset rearrangement potential
IF(tregpair) regpairadd=1D0

!     accumulate new pair-density
!NEW ODD
IF(thfbfold) THEN
  DO jst=1,nst(1)+nst(2)
    iq    = icharg(jst)
    IF(jst==istblo(iq) .OR. jst==istcblo(iq)) CYCLE  
    v2    = vsqf(jst)
    uv    = wtpair_uv(jst)*wtpair_v2(jst)*SQRT(MAX(v2-v2*v2,smallp))*0.5
    CALL Gaussfolding_spin(psi(0:nz,0:nr,1:2,jst),psifold(0:nz,0:nr,1:2,jst),nparz(1,jst),mpsi(jst),1)
    IF(testfolding) CALL sdenti(psi(0:nz,0:nr,1:2,jst),paircrtest(0:nz,0:nr,iq),uv,jst)
    IF(toddterms) THEN
      IF(mpsi(jst).lt.0) uv = (-1.0)*uv
      CALL cpairpot(psifold(0:nz,0:nr,1:2,jst),psifold(0:nz,0:nr,1:2,itpartn(jst)),paircr(0,0,iq),uv,jst)
    ELSE
      CALL sdenti(psifold(0:nz,0:nr,1:2,jst),paircr(0:nz,0:nr,iq),uv,jst)
    END IF
  END DO
ELSE
  DO jst=1,nst(1)+nst(2)
    iq    = icharg(jst)
    IF(jst==istblo(iq) .OR. jst==istcblo(iq)) CYCLE  
    v2    = vsqf(jst)
    uv    = wtpair_uv(jst)*wtpair_v2(jst)*SQRT(MAX(v2-v2*v2,smallp))*0.5
    IF(toddterms) THEN
      IF(mpsi(jst).lt.0) uv = (-1.0)*uv
      CALL cpairpot(psi(0:nz,0:nr,1:2,jst),psi(0:nz,0:nr,1:2,itpartn(jst)),paircr(0,0,iq),uv,jst)
    ELSE
      CALL sdenti(psi(0:nz,0:nr,1:2,jst),paircr(0:nz,0:nr,iq),uv,jst)
    END IF
  END DO
END IF
!END NEW ODD


!     now multiply with strength to obtain local pair-potential
IF(tFayans) THEN
  CALL fy_hfbpot()
ELSE
  v0prot = -pairfc(1)
  v0neut = -pairfc(2)

  IF(ipair <= 6 .OR. ipair == 9) THEN
    IF(tregpair .AND. iter>50) THEN
      DO ir=0,nr; DO iz=0,nz
        geff = geff_regpair(fer(1),pair_cutoff,pot(iz,ir,1),bmass(iz,ir,1),v0prot) 
        pairpt(iz,ir,1)  = paircr(iz,ir,1) * geff
        regpairadd(iz,ir,1) = 1D0 !(geff/v0prot)**2   ???
        geff = geff_regpair(fer(2),pair_cutoff,pot(iz,ir,2),bmass(iz,ir,2),v0neut) 
        pairpt(iz,ir,2)  = paircr(iz,ir,2) * geff
        regpairadd(iz,ir,2) = 1D0 !(geff/v0neut)**2   ???
      END DO; END DO
    ELSE
      pairpt(0:nz,0:nr,1)  = v0prot * paircr(0:nz,0:nr,1)
      pairpt(0:nz,0:nr,2)  = v0neut * paircr(0:nz,0:nr,2)
      IF(testfolding) THEN
        pairpttest(:,:,1)  = v0prot * paircrtest(:,:,1)
        pairpttest(:,:,2)  = v0neut * paircrtest(:,:,2)
      END IF
    END IF
  ELSE
    IF(tregpair .AND. iter>50) THEN
      DO ir=0,nr; DO iz=0,nz
        rhofac = (1D0-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr)
        geff = geff_regpair(fer(1),pair_cutoff,pot(iz,ir,1),bmass(iz,ir,1),v0prot*rhofac) 
        pairpt(iz,ir,1)  = paircr(iz,ir,1) * geff
        regpairadd(iz,ir,1) = (geff/(v0prot*rhofac))**2
        geff = geff_regpair(fer(2),pair_cutoff,pot(iz,ir,2),bmass(iz,ir,2),v0neut) 
        pairpt(iz,ir,2)  = paircr(iz,ir,2) * geff
        regpairadd(iz,ir,2) = (geff/(v0neut*rhofac))**2
      END DO; END DO
    ELSE
      pairpt(:,:,1)  = v0prot * (1D0-(rho(:,:,1)+rho(:,:,2))/rho0pr) * paircr(:,:,1)
      pairpt(:,:,2)  = v0neut * (1D0-(rho(:,:,1)+rho(:,:,2))/rho0pr) * paircr(:,:,2)
      IF(testfolding) THEN
        pairpttest(:,:,1)  = v0prot * (1D0-(rho(:,:,1)+rho(:,:,2))/rho0pr) * paircrtest(:,:,1)
        pairpttest(:,:,2)  = v0neut * (1D0-(rho(:,:,1)+rho(:,:,2))/rho0pr) * paircrtest(:,:,2)
      END IF
    END IF
  END IF
END IF

!     compute pairing energy
work(:,:)=paircr(:,:,1)*pairpt(:,:,1)
epairtestp=f1int(work)  
  
work(:,:)=paircr(:,:,2)*pairpt(:,:,2)
epairtestn = f1int(work)  
epairtest= epairtestp+ epairtestn
!   communicate to calling routine
epair(1)=epairtestp
epair(2)=epairtestn

IF(testfolding) THEN
  densint0 = f1int(paircrtest(:,:,2))
  densintfold = f1int(paircr(:,:,2))
  work(:,:)=paircrtest(:,:,2)*pairpttest(:,:,2)
  epairtestnfold = f1int(work)  
  pairpttest = (pairpt-pairpttest)**2
  paircrtest = (paircr-paircrtest)**2
  densintdiff = SQRT(f1int(paircrtest(:,:,2)))
  WRITE(6,'(a,8(1pg15.7))') 'COMPARE densities:', &
     densint0,densintfold,densintdiff,densint0-densintfold, &
     epairtestn,epairtestnfold,epairtestn-epairtestnfold, &
     f1int(pairpttest(:,:,2))
  OPEN(317,file='pairpttest_obs')
    DO ir=0,nr
      DO iz=0,nz
        WRITE(317,'(2i3,2(1pg13.5))') iz,ir,pairpttest(iz,ir,2),paircrtest(iz,ir,2)
      END DO
      WRITE(317,'(1x)')
    END DO
  CLOSE(317)
  DEALLOCATE(paircrtest,pairpttest)
ELSE
  densint0 = f1int(paircr(:,:,2))
  work(:,:)=paircr(:,:,2)*pairpt(:,:,2)
  epairtestnfold = f1int(work)  
!  WRITE(6,'(a,8(1pg15.7))') 'TEST densities:', &
!     densint0,epairtestn,epairtestnfold,epairtestn-epairtestnfold
END IF



! finally compute the actual gaps as s.p. expectation values with
! the pair potential
IF(thfbfold) THEN
  DO jst=1,nst(1)+nst(2)
    iq     = icharg(jst)
    work(:,:)=(psifold(:,:,1,jst)**2+psifold(:,:,2,jst)**2)*pairpt(:,:,iq)
    deltaf(jst) = f1int(work)
  END DO
ELSE
  DO jst=1,nst(1)+nst(2)
    iq     = icharg(jst)
    work(:,:)=(psi(:,:,1,jst)**2+psi(:,:,2,jst)**2)*pairpt(:,:,iq)
    deltaf(jst) = f1int(work)
  END DO
END IF

IF(itest > 0) THEN
  WRITE(6,'(a)') ' dynamics gaps for protons, directly:'
  WRITE(6,'(5(1pg12.4))') (deltaf(i),i=1,nst(1))
  WRITE(6,'(a)') ' dynamics gaps for neutrons, directly:'
  WRITE(6,'(5(1pg12.4))') (deltaf(i),i=nst(1)+1,nst(1)+nst(2))
END IF

!     stabilized pairing gap
IF(pairstab.NE.0D0) THEN
  IF(tFayans) THEN
    buff(1)=ABS(pairstab/epairtestp)
    buff(2)=ABS(pairstab/epairtestn)
    DO jst=1,nst(1)+nst(2)
      iq     = icharg(jst)
      deltaf(jst)=deltaf(jst)*(1.0D0 + buff(iq)**2)
    END DO
  ELSE
    ecrit=1.0D0
    npower=2.0D0
    buff(1)=ABS(ecrit/epairtestp)
    buff(2)=ABS(ecrit/epairtestn)
    DO jst=1,nst(1)+nst(2)
      iq     = icharg(jst)
      deltaf(jst)=deltaf(jst)*(1.0D0 - (npower-1.0D0)*pairstab*buff(iq)**npower)
    END DO
  END IF
END IF

DEALLOCATE(work)

RETURN
END SUBROUTINE hfbpot


REAL(8) FUNCTION geff_regpair(efermi_a,ecut_a,pot_a,bmas_a,g_pair)
!
!  Computes the regularization kernel g_eff for
!  regularized pairing.
!  
IMPLICIT NONE
REAL(8),INTENT(IN) :: efermi_a     ! Fermi energy
REAL(8),INTENT(IN) :: ecut_a       ! Cutoff energy
REAL(8),INTENT(IN) :: pot_a        ! local potential U(r)
REAL(8),INTENT(IN) :: bmas_a       ! local inverse effective mass B*(r)=hbar**2/(2M*)
REAL(8),INTENT(IN) :: g_pair       ! local pairing strength from pairing model

  geff_regpair = g_pair/(1-g_pair*Keff_regpair())
!  geff_regpair = g_pair/(1-g_pair*Keff_regpair2())

  RETURN

CONTAINS

REAL(8) FUNCTION Keff_regpair()
!
!  Computes the regularization kernel K_eff for
!  regularized pairing.
!  

REAL(8),PARAMETER :: fourpi=12.56637061D0   ! 4*pi
REAL(8),PARAMETER :: fourpi2=39.4784176D0   ! 4*pi**2
REAL(8) :: kf2                     ! square of Fermi momentum
REAL(8) :: kf                      ! absolute value of Fermi momentum
REAL(8) :: kcut                    ! cutoff momentum, deduced from ecut_a
REAL(8) :: ecut_eff                ! cutoff energy relatice to cont. threshold

kf2= (efermi_a-pot_a)/bmas_a
kf = SQRT(ABS(kf2))
ecut_eff = ecut_a  !+efermi_a
kcut = (ecut_eff+efermi_a-pot_a)/bmas_a
IF(kcut<0D0) THEN
  WRITE(6,'(a,3(1pg13.5))') ' e_cut,e_fermi,pot=',ecut_a,efermi_a,pot_a
  STOP "'kcut' is imaginary"
END IF
kcut=SQRT(kcut)

IF(kf2.GE.0D0) THEN
  Keff_regpair = (kcut/(fourpi2*bmas_a)) * &
                 (1D0 - (kf/(kcut+kcut))*LOG((kcut+kf)/(kcut-kf))  )
ELSE
  Keff_regpair = (kcut/(fourpi2*bmas_a)) * &
                 (1D0 + (kf/(kcut))*ATAN((kf/(kcut)))  )
END IF

RETURN

END FUNCTION Keff_regpair


REAL(8) FUNCTION Keff_regpair2()
!
!  Computes the regularization kernel K_eff for
!  regularized pairing.
!  

REAL(8),PARAMETER :: fourpi=12.56637061D0   ! 4*pi
REAL(8),PARAMETER :: fourpi2=39.4784176D0   ! 4*pi**2
!REAL(8) :: kf2                     ! square of Fermi momentum
REAL(8) :: kf                      ! absolute value of Fermi momentum
REAL(8) :: kcut                    ! cutoff momentum, deduced from ecut_a
REAL(8) :: lcut                    ! lower cutoff momentum, deduced from ecut_a
REAL(8) :: ecut_eff                ! cutoff energy relatice to cont. threshold

REAL(8) :: a,b                     ! auxiliary

ecut_eff = ecut_a  !+efermi_a

IF ( pot_a + ecut_eff - efermi_a < 0D0 ) THEN
   !
   kf = SQRT( ( efermi_a - pot_a ) / bmas_a )
   kcut = SQRT( ( efermi_a + ecut_eff - pot_a ) / bmas_a )
   lcut = SQRT( ( efermi_a - ecut_eff - pot_a ) / bmas_a )
   a = kcut / fourpi2 / bmas_a &
      * ( 1 - kf / 2 / kcut * log( ( kcut + kf ) / ( kcut - kf ) ) )
   b = lcut / fourpi2 / bmas_a &
      * ( 1 - kf / 2 / lcut * log( ( lcut + kf ) / ( kf - lcut ) ) )
!         rega(p,it) = a
!         regb(p,it) = b
   !........................ geff(p,it) = 1 / ( 1 / g - a - b )
ELSE
   IF ( pot_a - efermi_a < 0D0 ) THEN
      !
      kf = SQRT( ( efermi_a - pot_a ) / bmas_a )
      kcut = SQRT( ( efermi_a + ecut_eff - pot_a ) / bmas_a )
      a = kcut / fourpi2 / bmas_a * &
         ( 1 - kf / 2 / kcut * log( ( kcut + kf ) / ( kcut - kf ) ) )
      b = 0D0
!            rega(p,it) = a
!            regb(p,it) = 0D0
      !..................... geff(p,it) = 1 / ( 1 / g - a )
   ELSE
      IF ( pot_a - ecut_eff - efermi_a < 0D0 ) THEN
         !
         kf = SQRT( - ( efermi_a - pot_a ) / bmas_a )
         kcut = SQRT( ( efermi_a + ecut_eff - pot_a ) / bmas_a )
         a = kcut / fourpi2 / bmas_a &
            * ( 1 + kf / kcut * atan( kf / kcut ) )
         b = 0D0
!               rega(p,it) = a
!               regb(p,it) = 0D0
         !.................. geff(p,it) = 1 / ( 1 / g - a )
      ELSE
         !
         a = 0D0
         b = 0D0
!               rega(p,it) = 0D0
!               regb(p,it) = 0D0
         !.................. geff(p,it) = g
      END IF
   END IF
END IF
Keff_regpair2 = a+b

RETURN

END FUNCTION Keff_regpair2

END FUNCTION geff_regpair


!---------------------------------------------------------------------------  
! DESCRIPTION: pairdn
!> @brief
!! Determines the fermi-energy \c eferm such that the "particle number 
!! condition" is fulfilled. This is achieved by iteration.
!>
!> @param[in] e
!> REAL(8), array, takes single particle energies.
!> @param[in] delta 
!> REAL(8), array, takes the pairing gap parameters.
!> @param[inout] gw
!> REAL(8), array, takes occupation probabilities, input is used as initial 
!! guess for the iteration. Returns the resulting occupation probabilities.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[inout] pairgw_uv
!> REAL(8), array, takes space for dynamical u*v cut-off factor.
!> @param[inout] pairgw_v2
!> REAL(8), array, takes space for dynamical v**2 cut-off factor.
!> @param[inout] partnm
!> REAL(8), takes the desired particle number and returns the particle number, 
!! from summing \c gw with degeneracy \c ph.
!> @param[out] eferm
!> REAL(8), returns Fermi energy, input is used as initial guess; if eferm=0.0 
!! is input then a suitable initial guess is computed from input \c gw.
!> @param[out] epair
!> REAL(8), returns the pairing energy.
!> @param[out] avdelt
!> REAL(8), returns the average pairing gap.
!> @param[out] efadd
!> REAL(8), pairing band above 'eferm', in case of dynamical pair space see 
!! routine  \c cutoff.
!> @param[in] nmax
!> INTEGER, takes maximum number of single particle states.
!> @param[in] ndim
!> INTEGER, takes the dimensioning of the arrays, particularly of \c vpairf,
!! in the calling routine.
!> @param[in] iter
!> INTEGER, takes number of pairing iterations.
!> @param[in] ipair
!> INTEGER, takes pairing type.
!> @param[in] iwfoff
!> INTEGER, takes offset for addressing wavefunctions.
!> @param[in] nblock
!> INTEGER, takes address of blocked state (none if 0).
!> @param[in] eps
!> INTEGER, takes end condition for pairing iteration; iteration terminates 
!! if error in particle number is smaller than \c eps and if \c eferm and/or 
!! \c delta are smaller than \c eps.
!> @param[in] iab
!> INTEGER, takes print level; significant values are -2, -1, 0 and 1;
!! larger value prints more.
!--------------------------------------------------------------------------- 
SUBROUTINE pairdn(e,delta,gw,ph,pairgw_uv,pairgw_v2,partnm,eferm,epair,avdelt,  &
    efadd,nmax,ndim,iter,ipair,iwfoff,nblock, eps,iab)

USE paramr
USE levels, ONLY: npsi
USE dynpar, ONLY: pair_cutoff,tcutoff_WS,tcutoff_uvspace
IMPLICIT NONE

REAL(8), INTENT(IN)       :: e(ndim)
REAL(8), INTENT(IN)       :: delta(ndim)
REAL(8), INTENT(INOUT)    :: gw(ndim)
REAL(8), INTENT(IN)       :: ph(ndim)
REAL(8), INTENT(INOUT)    :: pairgw_uv(ndim)
REAL(8), INTENT(INOUT)    :: pairgw_v2(ndim)
REAL(8), INTENT(IN)       :: partnm
REAL(8), INTENT(IN OUT)   :: eferm
REAL(8), INTENT(OUT)      :: epair
REAL(8), INTENT(OUT)      :: avdelt
REAL(8), INTENT(OUT)      :: efadd
INTEGER, INTENT(IN)       :: nmax
INTEGER, INTENT(IN)       :: ndim
INTEGER, INTENT(IN)       :: iter
INTEGER, INTENT(IN)       :: ipair
INTEGER, INTENT(IN)       :: iwfoff
INTEGER, INTENT(IN)       :: nblock
REAL(8), INTENT(IN)       :: eps
INTEGER, INTENT(IN)       :: iab

REAL(8),ALLOCATABLE :: vln(:)  ! workspace for LN

!     several numerical parameters, here for nuclei

REAL(8),PARAMETER :: delefe=20D0   ! difference of Fermi energy for derivative
INTEGER,PARAMETER :: ifcut=1       ! switch for flexible cutoff
!                                   ifcut = 0  --> no cutoff
!                                   ifcut = 1  --> cutoff 1.65 A**2/3
!                                   ifcut = 2  --> cutoff 5 MeV
REAL(8),SAVE      :: xlam2=0D0     ! to make LN parameter 'lambda_2' resident
REAL(8),PARAMETER :: xl2stp=1D0    ! trial step for LN parameter
REAL(8),PARAMETER :: dxllim=0.5D0  ! limiting 'eferm'-step in case of LN
REAL(8),PARAMETER :: dx2lim=0.1D0  ! limiting 'xlam2'-step in case of LN
REAL(8),PARAMETER :: detmin=1D-4   ! lower limit for determinant in LN step
REAL(8),PARAMETER :: dampin=0.5D0  ! damping of LN step
REAL(8),PARAMETER :: half=0.5D0, xsmallp=1D-20, fifty=50.0d0

INTEGER :: ifirst=0,i,iftest,il1,il2,inhanc,it,iterh,iterp,n,na,npart
REAL(8) :: actnum,arg,dchi1,dchi11,dchi12,dchi2,dchi22,delcut,deltaa,deltac,det,&
           dinv11,dinv12,dinv22,dl2f1,dl2f2,dlf1,dlf2,dxlam,dxlam2,edif,enhanc,&
           equasi,f1,f1xl,f1xl2,f2,f2xl,f2xl2,parnum,partsm,phact,quenf,quenf1,&
           quenf2,sumduv,sumepa,sumuv,tlam,tlam2,tstnum,v2,vol
REAL(8) :: avg
LOGICAL :: converged

ALLOCATE(vln(npsi(2)))
vln=0D0
IF(iab > -1) THEN
  WRITE(6,'(/a,i4,1(a,f9.4)/a,10(/1x,10f8.3))')  &
      ' *** pairing with dynamic gap entered with: iter=',iter,  &
      '  eferm=',eferm, '  single particle energies:',  &
      (e(i),i=1,nmax)
  WRITE(6,'(a,10(/1x,10f8.4))') '  delta:',(delta(i),i=1,nmax)
  WRITE(6,'(a,10(/1x,10f8.4))') '    ph:',(ph(i),i=1,nmax)
END IF

! preset pairing weights, use pairgw_uv as workspace
IF(eferm /= 0D0) THEN
  IF(tcutoff_WS) THEN
    efadd  = pair_cutoff
    deltac = efadd*0.1D0
    DO i=1,nmax
      arg    = (e(i)-(eferm+efadd))/deltac
      IF(arg > fifty) THEN
        pairgw_uv(i)   = 0D0
      ELSE IF(arg < -fifty) THEN
        pairgw_uv(i)   = 1D0
      ELSE
        pairgw_uv(i)  = 1D0/(1D0+EXP(arg))
      END IF
    END DO
  ELSE 
    actnum = partnm + pair_cutoff*partnm**0.6666667      ! active pairing space 
!              sharp cutoff with integer nr. of states for 3D writing
    delcut = 0.1
    efadd  = 20.0/partnm**0.33333333
    CALL cutoff(e,pairgw_uv,ph,actnum,nmax,eferm,efadd,delcut, tstnum,iter,eps,-2)
  END IF
ELSE
  pairgw_uv(1:nmax) = one
END IF
! place pairing weight into correct field
IF(tcutoff_uvspace) THEN
  pairgw_v2(1:nmax)=1D0
ELSE
  pairgw_v2(1:nmax)=pairgw_uv(1:nmax)
  pairgw_uv(1:nmax)=1D0
END IF

!      prepare fermi-energy 'eferm', if not yet initialised.
converged=.FALSE.
IF(eferm == 0D0) THEN
  partsm = 0D0
  DO it=1,nmax
    partsm = ph(it) + partsm
    IF(partsm-partnm >= 0D0) THEN
      converged=.TRUE.
      EXIT
    END IF
  END DO
  IF(.NOT.converged) STOP ' not enough states to reach particle number'
  eferm  = 0.5*(e(it)+e(it+1))
  IF(iab >= 1) WRITE(06,'(1x,a,2g11.4,a,g11.4)')  &
      'initialization: partsm,partnm=',partsm,partnm,' eferm=',eferm
  DO n=1,nmax
    IF(n <= it) THEN
      gw(n) = one
    ELSE
      gw(n) = 0D0
    END IF
  END DO
  xlam2 = 0.0
END IF

! the coupled gradient iteration for two cases, switched by 'ipair'
IF(ipair /= 9.AND.ipair /= 11) THEN
  
! standard case where p.n. is iterated, misusing 'gw' as workspace
  eferm = rbrent(e, delta, ph, gw, pairgw_uv,pairgw_v2, partnm, ndim, nmax, nblock, iab)
  parnum = parnmt(e, delta, ph, gw, pairgw_uv,pairgw_v2, eferm, ndim, nmax, nblock)
  it=0
  
!       the LN case where p.n. and p.n. fluctuations are iterated
ELSE
  iftest  = iab + 1
  iterh   = iter/2
  iterp   = iter/10
  converged=.FALSE.
  DO it=1,iter
    IF(iftest > 0) WRITE(6,'(a,i5,a)')  ' iteration nr.',it,' :'
    
!         initial call to pairing equations
    npart  = NINT(partnm)
    IF(npart == 12 .AND. ifirst == 1) THEN
      ifirst = 0
      WRITE(6,'(a)') '  lambda     lambda_2            f1           f2    '
      DO il2=0,20
        tlam2 = -0.5 + il2*0.05
        DO il1=0,20
          tlam  = -25.0 + il1*0.5
          CALL lm2def(e,delta,ph,pairgw_uv,pairgw_v2,1,nmax,iwfoff,npart,  &
              tlam,tlam2,gw,gw,iftest, f1,f2)
          WRITE(6,'(1x,2f8.2,2(1pg12.4))') tlam,tlam2,f1,f2,f1**2+5*f2**2
        END DO
      END DO
      STOP
    END IF
    CALL lm2def(e,delta,ph,pairgw_uv,pairgw_v2,1,nmax,iwfoff,npart,  &
        eferm,xlam2,gw,gw,iftest, f1,f2)
    IF(iftest > 0) WRITE(6,'(a,2g12.3,a,2g12.3)')  &
        ' xlam,xlam2=',eferm,xlam2,'   f1,f2=',f1,f2
    
! test termination
    IF(it > 1 .AND. ABS(f1) < eps .AND. ABS(f2) < eps) THEN
      converged=.TRUE.
      EXIT
    END IF
! derivatives by finite differences
    CALL lm2def(e,delta,ph,pairgw_uv,pairgw_v2,1,nmax,iwfoff,npart,  &
        eferm+delefe,xlam2,gw,vln,iftest, f1xl,f2xl)
    CALL lm2def(e,delta,ph,pairgw_uv,pairgw_v2,1,nmax,iwfoff,npart,  &
        eferm,xlam2+xl2stp,gw,vln,iftest, f1xl2,f2xl2)
    IF(iftest >= 0) WRITE(6,'(a,i4,6(1pg11.3))')  &
        ' LN: iter,xlam,xlam2,f1,f2=',it,eferm,xlam2,f1,f2
    dlf1   = (f1xl-f1)/delefe
    dlf2   = (f2xl-f2)/delefe
    dl2f1  = (f1xl2-f1)/xl2stp
    dl2f2  = (f2xl2-f2)/xl2stp
    IF(iftest >= 1) WRITE(6,'(a,4g12.4)')  &
        ' dlf1,dlf2,dl2f1,dl2f2=',dlf1,dlf2,dl2f1,dl2f2
    
! propagation to improved solution by inverse matrix
    IF(it < iterh) THEN
      det    = (dlf1*dl2f2-dlf2*dl2f1) !   - sign ?
      det    = SIGN(MAX(ABS(det),detmin),det) ! lower limit for 'det'
      dxlam  = (f2*dl2f1-f1*dl2f2)/det
      dxlam2 = (f1*dlf2-f2*dlf1)/det
      IF(iftest >= 0) WRITE(6,'(a,4g12.4)')  &
          ' dir: dxlam,dxlam2=',dxlam,dxlam2
      
! propagation by minimzation of least squares .
! this as such is equivalent to direct gradient,
! but here we can employ the "damping" by systematically
! enhancing the diagonal elements (see Bevington step).
    ELSE
      dchi1  = f1*dlf1 + 5.0*f2*dlf2
      dchi2  = f1*dl2f1 + 5.0*f2*dl2f2
      dchi11 = dlf1*dlf1 + 5.0*dlf2*dlf2
      dchi12 = dl2f1*dlf1 + 5.0*dl2f2*dlf2
      dchi22 = dl2f1*dl2f1 + 5.0*dl2f2*dl2f2
      inhanc = (it-iterh)/iterp+1
      enhanc = 2.0**inhanc
      dchi11 = dchi11*enhanc
      dchi22 = dchi22*enhanc
      det    = dchi11*dchi22-dchi12*dchi12
      det    = SIGN(MAX(ABS(det),detmin),det) ! lower limit for 'det'
      dinv11 = dchi22/det
      dinv22 = dchi11/det
      dinv12 = -dchi12/det
      dxlam  = -(dinv11*dchi1+dinv12*dchi2)
      dxlam2 = -(dinv12*dchi1+dinv22*dchi2)
      IF(iftest >= 0) WRITE(6,'(a,4g12.4)')  &
          ' chi: dxlam,dxlam2=',dxlam,dxlam2
    END IF
    
    IF(iftest >= 1) WRITE(6,'(a,4g12.4)')  &
        ' dxlam,dxlam2=',dxlam,dxlam2
    quenf1 = dxllim/MAX(ABS(dampin*dxlam ),dxllim)
    quenf2 = dx2lim/MAX(ABS(dampin*dxlam2),dx2lim)
    quenf  = MIN(quenf1,quenf2)
    dxlam  = quenf*dampin*dxlam
    dxlam2 = quenf*dampin*dxlam2
    IF(iftest >= 1) WRITE(6,'(a,4g12.4)')  &
        ' dxlam,dxlam2=',dxlam,dxlam2
    eferm  = eferm + dxlam
    xlam2  = xlam2 + dxlam2
    
  END DO
  
  IF(iftest >= 0.AND..NOT.converged) WRITE(6,'(a)')  &
      ' no convergence of pairing parameters in PAIRDN, LN-branch'

  parnum = f1+npart
  IF(iftest >= 0) WRITE(6,'(a,i4,6(1pg11.3))')  &
      ' PAIRDN-LN: iter,xlam,xlam2,f1,f2=',it,eferm,xlam2,f1,f2
END IF

! finally compute average gap and pairing energy
sumuv  = 0.0
sumduv = 0.0
sumepa = 0.0
DO na=1,nmax
  deltaa = delta(na)*pairgw_uv(na)
  edif   = e(na)-eferm
  equasi = SQRT(edif*edif+deltaa*deltaa)
  v2     = half-half*edif/equasi
  phact  = half*ph(na)   ! for axial
  vol    = phact*SQRT(MAX(v2-v2*v2,xsmallp))*pairgw_v2(na)
  sumuv  = vol        + sumuv
  sumduv = vol*deltaa + sumduv
  sumepa = phact*half*(deltaa*deltaa)/equasi + sumepa
END DO
sumuv  = MAX(sumuv,xsmallp)
avdelt = sumduv/sumuv
epair  = sumduv           ! or - ?
avg    = epair/sumuv**2
IF(iab >= -1) WRITE(6,'(a,i4,a,g13.5,a,2g13.5)')  &
    ' PAIRDN: iteration=',it,' part.nr.=',parnum, ' eferm,epair=',eferm,epair
IF(iab >= 0) WRITE(6,'(5x,a,3g12.4)')  &
    ' av.Delta,av.G-strength,pair.E=',avdelt,avg,epair
DEALLOCATE(vln)
RETURN
END SUBROUTINE pairdn

!---------------------------------------------------------------------------  
! DESCRIPTION: cutoff
!> @brief
!!This routine determines the effective Fermi energy and cut-off 
!!distribution for the pairing space.
!>
!> @param[in] e
!> REAL(8), array, takes single particle energies.
!> @param[inout] gw
!> REAL(8), array, takes occupation probabilities, input is used as initial 
!! guess for the iteration. Returns the resulting occupation probabilities.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[inout] anz
!> REAL(8), takes the desired effective particle number.
!> @param[in] nmax
!> INTEGER, takes maximum number of single particle states.
!> @param[inout] eferm
!> REAL(8), takes "true" Fermi energy as given by pairing routine.
!> @param[inout] efadd
!> REAL(8), pairing band above \c eferm, in case of dynamical pair space.
!> @param[inout] delta 
!> REAL(8), takes the relative width of the distribution.
!> @param[out] partnm 
!> REAL(8), returns particle number, from adding 'gw' with degeneracy 'ph'.
!> @param[in] iter
!> INTEGER, takes number of pairing iterations.
!> @param[in] eps
!> INTEGER, takes end condition for pairing iteration; iteration terminates 
!! if error in particle number is smaller than \c eps and if \c eferm and/or 
!! \c delta are smaller than \c eps.
!> @param[in] iab
!> INTEGER, takes print level; significant values are -2, -1, 0 and 1;
!! larger value prints more.
!--------------------------------------------------------------------------- 
SUBROUTINE cutoff(e,gw,ph,anz,nmax,eferm,efadd,delta, partnm,iter,eps,iab)

USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)         :: e(*)
REAL(8), INTENT(INOUT)      :: gw(*)
REAL(8), INTENT(IN)         :: ph(*)
REAL(8), INTENT(IN)         :: anz
INTEGER, INTENT(IN)         :: nmax
REAL(8), INTENT(INOUT)      :: eferm
REAL(8), INTENT(INOUT)      :: efadd
REAL(8), INTENT(INOUT)      :: delta
REAL(8), INTENT(OUT)        :: partnm
INTEGER, INTENT(IN)         :: iter
REAL(8), INTENT(IN)         :: eps
INTEGER, INTENT(IN)         :: iab

REAL(8) :: dparn(2)=(/0D0,0D0/) ! variations in the particle number as
!                                 consequence of 'deferm' and 'ddelta'.
REAL(8) :: gr(1,1)=0D0          ! double array of variation of particle number and of
!                                 gap equation, used for construction of the tangantial step

REAL(8),PARAMETER :: dampin=0.9D0 !damping factor for the tangential step

REAL(8),PARAMETER :: deferm=1D0   !step in the fermi-energy, to evaluate
!                                  numerically the derivatives.
! efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
REAL(8),PARAMETER :: efmax=3D0
REAL(8),PARAMETER :: detmin=1D-4
REAL(8),PARAMETER :: fifty=50D0
INTEGER,PARAMETER :: iftang=0

INTEGER :: k,ka,iwarng,it,itsum,isum,n,i
REAL(8) :: parnum,arg,deltac,det,dparn1,dparn2,dparn3,dparnm,eferm1,eferm2,eferm3,efstep,elam
LOGICAL :: converged
!-------------------------------------------------------------------------------
!      set switch variables
iwarng = 0

IF(iab >= 1) WRITE(6,'(/a,2i4,a,4f9.4)')  &
    ' *** cutoff entered with: iter,nmax=',iter,nmax,  &
    '  delta,efadd,eferm,partnm=',delta,efadd,eferm,anz
IF(iab >= 2) WRITE(6,'(a,10(/1x,10f8.3))')  &
    ' single particle energies:',(e(i),i=1,nmax)

!      prepare fermi-energy 'eferm', if not yet initialised.
converged=.FALSE.
IF(efadd == 0D0) THEN
  isum   = 0
  DO it=1,nmax
    itsum  = it
    isum   = INT(ph(it)+.000001) +isum
    IF(one*isum-anz >= 0) THEN
      converged=.TRUE.
      EXIT
    END IF
  END DO
  IF(.NOT.converged) STOP ' not enough states too reach particle number'
  efadd  = e(itsum)+deferm-eferm
  IF(iab >= 1) WRITE(06,'(5x,a,i5,a,g11.4)')  &
      'initialistion: isum=',isum,'  efadd=',efadd
END IF

!-------------------------------------------------------------------------------
!     iterate Fermi energy and the distribution
!-------------------------------------------------------------------------------

!     tangent step
IF(iftang == 1) THEN
  converged=.FALSE.
  DO it=1,iter
! particle number at 'eferm' (k=1) and 'eferm' varied (k=2).
    DO ka=1,2
! set actual fermi-energy 'elam'
      k      = 3-ka
      deltac  = efadd*delta
      IF(k == 2) THEN
        elam   = eferm+efadd+deferm*deltac
      ELSE
        elam   = eferm+efadd
      END IF
! compute actual particle number
      parnum = 0D0
      DO i=1,nmax
        arg    = (e(i)-elam)/deltac
        IF(arg > fifty) THEN
          gw(i)   = 0D0
        ELSE IF(arg < -fifty) THEN
          gw(i)   = one
        ELSE
          gw(i)  = 1.0/(1.0+EXP(arg))
        END IF
        parnum = gw(i)*ph(i)    + parnum
      END DO
! store particle number
      dparn(k) = parnum-anz                                 ! axial
    END DO
    
!      gradient step to next 'eferm'   (including step limiters)
    gr(1,1) = (dparn(2)-dparn(1))/(deferm*deltac)
    det    = SIGN(MAX(ABS(gr(1,1)),detmin),gr(1,1))
    efstep = dampin*dparn(1)/det
    efadd  = efadd-SIGN(MIN(ABS(efstep),efmax),efstep)
    IF(iab == 0) WRITE(6,'(a,i4,2(a,g11.3))')  &
        ' pairing-weight iteration it=',it,'  :  efadd=',efadd,  &
        '  dparn(1)=',dparn(1)
    IF(iab >= 1) WRITE(6,'(a,i4,2(a,g11.3)/a,2g11.3,a,g11.3/  &
        &a,10f7.4,10(/6x,10f7.4))') ' iteration at it=',it,' :   efadd=',efadd,  &
        '  det=',det, ' dparn=',dparn(1),dparn(2),' ;   gr=',gr(1,1),  &
        '   gw=',(gw(n),n=1,nmax)
    IF(ABS(dparn(1)) < eps) THEN
      converged=.TRUE.
      EXIT
    END IF
  END DO
  IF(.NOT.converged) iwarng = 1
ELSE
  
!     secant step
  eferm1 = 1.0D10
  eferm2 = -1.0D10
  DO i=1,nmax
    IF(eferm2 < e(i)) eferm2 = e(i)
    IF(eferm1 > e(i)) eferm1 = e(i)
  END DO

  DO it=1,iter
    IF(it == 1) THEN
      eferm3 = eferm1
    ELSE IF(it == 2) THEN
      eferm3 = eferm2
    ELSE
      eferm3 = eferm1 - dparn1*(eferm2-eferm1)/(dparn2-dparn1)
    END IF
    efadd   = eferm3-eferm
    deltac  = MAX(efadd*delta,delta)
    parnum  = 0D0
    DO i=1,nmax
      arg    = (e(i)-eferm3)/deltac
      IF(arg > fifty) THEN
        gw(i)   = 0D0
      ELSE IF(arg < -fifty) THEN
        gw(i)   = one
      ELSE
        gw(i)  = 1.0/(1.0+EXP(arg))
      END IF
      parnum = gw(i)*ph(i)    + parnum
    END DO
    dparnm = parnum-anz                                      ! axial
    IF(it == 1) THEN
      dparn1  = dparnm
      IF(dparn1 > 0D0) THEN
        WRITE(6,'(a/i4,5(1pg12.4))')  &
            ' it,deltac,eferm1,eferm2,dparn1,dparn2=',  &
            it,deltac,eferm1,eferm2,dparn1,dparn2
        WRITE(6,'(a,10(/1x,10f8.3))')  &
            ' single particle energies:',(e(i),i=1,nmax)
        STOP ' lower state too high in CUTOFF'
      END IF
    ELSE IF(it == 2) THEN
      dparn2  = dparnm
      IF(dparn2 < 0D0) THEN
        WRITE(6,'(a/i4,5(1pg12.4))')  &
            ' it,deltac,eferm1,eferm2,dparn1,dparn2=',  &
            it,deltac,eferm1,eferm2,dparn1,dparn2
        WRITE(6,'(a,10(/1x,10f8.3))')  &
            ' single particle energies:',(e(i),i=1,nmax)
        STOP ' upper state too low in CUTOFF'
      END IF
    ELSE
      dparn3  = dparnm
    END IF
    IF(iab >= 2) WRITE(6,'(a,i4,2(a,g11.3))')  &
        ' secant iteration it=',it,'  :  eferm=',eferm3, '  dparn3=',dparn3
    IF(it > 2) THEN
      IF(ABS(dparn3) < eps) THEN
        EXIT
      ELSE IF(dparn3 > 0D0) THEN
        eferm2 = eferm3
        dparn2 = dparn3
      ELSE
        eferm1 = eferm3
        dparn1 = dparn3
      END IF
    END IF
  END DO
END IF

!-------------------------------------------------------------------------------
!     finally the particle-number is checked and printed
!-------------------------------------------------------------------------------
partnm = 0D0
DO i=1,nmax
  partnm = gw(i)*ph(i) + partnm
END DO

IF(iab >= -1) WRITE(6,'(a,i4,2(1pg13.5))')  &
    ' CUTOFF finished with iter,part.nr.,efadd=',it,partnm,efadd
IF(iab >= 2) WRITE(6,'(a,10(/1x,10f8.3))')  &
    ' pairing-weights:',(gw(i),i=1,nmax)

!     in case of incomplete convergence print a warning

IF(iab >= -1 .AND. iwarng == 1) WRITE(6,'(2(a,g11.3),2(a,i3))')  &
    ' ---> no convergence in cutoff! err(part.nr)=',dparn(1)

RETURN
END SUBROUTINE cutoff

!---------------------------------------------------------------------------  
! DESCRIPTION: lm2def
!> @brief
!!Routine to compute the second order multiplies \c lambda_2 for the LN 
!!scheme in connection with the delta-force pairing.
!>
!> @param[in] spen
!> REAL(8), array, takes single particle energies.
!> @param[in] delta
!> REAL(8), array, takes pairing gaps.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[in] pairgw_uv
!> REAL(8), array, takes the u*v cutoff weights complementing gap \c delta.
!> @param[in] pairgw_uv
!> REAL(8), array, takes the v**2 cutoff weights complementing gap \c delta.
!> @param[in] nspd
!> INTEGER, takes lowest level in pairing.
!> @param[in] nspg
!> INTEGER, takes highest level in pairing.
!> @param[in] iwfoff
!> INTEGER, takes offset for addressing wave functions.
!> @param[in] npart
!> INTEGER, takes the particle number.
!> @param[in] xlam
!> REAL(8), takes the Fermi energy.
!> @param[in] xlam2
!> REAL(8), takes the initial value for fluctuation constraint.
!> @param[in] vv2in
!> REAL(8), array, takes the initial internal occupation probabilities.
!> @param[out] vv2out
!> REAL(8), array, returns the new occupation weights.
!> @param[in] iftest
!> INTEGER, switch for testprints.
!> @param[out] f1
!> REAL(8), returns the particle number condition: \f$ \sum(2v^2)-n \f$.
!> @param[out] f1
!> REAL(8), returns the Lipkin-Nogami condition: 
!> \f$ \texttt{xlam2}-\texttt{vpair}/4\cdot (\sum(u^3\cdot v)\sum(v^3\cdot u)-
!! \sum(u^4\cdot v^4))/((\sum(u^2\cdot v^2))^2-\sum(u^4\cdot v^4)) \f$.
!--------------------------------------------------------------------------- 
SUBROUTINE lm2def(spen,delta,ph,pairgw_uv,pairgw_v2,nspd,nspg,iwfoff,npart,  &
    xlam,xlam2,vv2in,vv2out,iftest,f1,f2)

USE paramr
USE dynpar
USE wavefr
USE grid
USE fields, ONLY: rho
USE observ, ONLY: sdenti

IMPLICIT NONE

REAL(8), INTENT(IN)                   :: spen(*)
REAL(8), INTENT(IN)                   :: delta(*)
REAL(8), INTENT(IN)                   :: ph(*)
REAL(8), INTENT(IN)                   :: pairgw_uv(*)
REAL(8), INTENT(IN)                   :: pairgw_v2(*)
INTEGER, INTENT(IN)                   :: nspd
INTEGER, INTENT(IN)                   :: nspg
INTEGER, INTENT(IN)                   :: iwfoff
INTEGER, INTENT(IN)                   :: npart
REAL(8), INTENT(IN)                   :: xlam
REAL(8), INTENT(IN)                   :: xlam2
REAL(8), INTENT(IN OUT)               :: vv2in(*)
REAL(8), INTENT(OUT)                  :: vv2out(*)
INTEGER, INTENT(IN)                   :: iftest
REAL(8), INTENT(OUT)                  :: f1
REAL(8), INTENT(OUT)                  :: f2

INTEGER, PARAMETER :: iterv2=50     !max.nr. iterations in 'vv2'
REAL(8), PARAMETER :: epsv2=1.0d-7  !precision on 'vv2'
REAL(8), PARAMETER :: epsdv2=1.0d-7 !limit on change of energy in determining 'vv2'
REAL(8), PARAMETER :: dampin=0.5    !general reduction of the step, top stabilize
REAL(8), PARAMETER :: half=0.5

INTEGER :: istate,it,k,iz,ir
REAL(8) :: alpha,alphd2,d2,denom,derf,dx,e,enume1,enume2,eq,erspln,f,fnn,su3v,&
           suv,suv2,suv2ra,suv3,suv4,sv2,u2,uv,v2,vpair,x,xmult
REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
wrk1=0D0;wrk2=0D0
! preparation
fnn    = npart
IF(iwfoff == 0) THEN
  erspln = rspln(1)
  suv2ra = suv2rs(1)
  vpair  = pairfc(1)
ELSE
  erspln = rspln(2)
  suv2ra = suv2rs(2)
  vpair  = pairfc(2)
END IF
IF(iftest >= 3) WRITE(6,'(a,i5,3g12.4)')  &
    ' iwfoff,vpair,erspln,suv2ra=',iwfoff,vpair,erspln,suv2ra

!     reset accumulators for modified pairing currents
wrk1 = 0.0d0
wrk2 = 0.0d0

!     reset accumulators for combinations of 'u' and 'v'.
sv2   = 0D0                                  ! sum v^2
suv   = 0D0                                  ! sum u*v
suv4  = 0D0                                  ! sum u^4 v^4
suv2  = 0D0                                  ! sum u^2 v^2
su3v  = 0D0                                  ! sum u*v u^2
suv3  = 0D0                                  ! sum u*v v^2

alpha  = 4.*xlam2

DO k=nspd,nspg
! iteration loop to compute occupation 'vv2'
  d2     = (delta(k)*pairgw_uv(k))**2
  alphd2 = alpha*half*d2
  e     = spen(k)-xlam
  v2    = vv2in(k)                                ! start from previous value
  x     = e+alpha*v2*pairgw_v2(k)
  IF(iftest >= 3) WRITE(7,'(a,i5,a,3(1pg12.4))')  &
      ' state nr.',k,' :  v2in,e,delta=',v2,e,delta(k)*pairgw_uv(k)
  
! iterate to s.c. occup.numbers and lambda_2
  DO it=1,iterv2
    EQ    = SQRT(d2+x*x)
    v2    = half-half*x/EQ
    f     = e+alpha*v2*pairgw_v2(k)-x
    derf  = -alphd2/(EQ*EQ*EQ)-one
    dx    = dampin*f/derf         ! optionally switch on if unstable
    IF(ABS(f) < epsv2 .AND. ABS(dx) < epsdv2) GO TO 7
    x     = x-dx
  END DO

  IF(iftest >= 3) WRITE(7,'(a)')  &
      ' maximal number of iterations exceeded in "LM2DEF"'
  
  7   CONTINUE
  v2   = MAX(SQRT(smal),v2)
  IF(iftest >= 1) WRITE(7,'(a,i5,a,i5,g12.3)')  &
      ' state nr.',k,' :  iterv,v2',it,v2
  vv2out(k)=v2
  
! accumulate sums with new 'v2'
  xmult = ph(k)*pairgw_v2(k)
  u2    = one-v2
  u2    = MAX(SQRT(smal),u2)
  uv    = SQRT(u2*v2)
  sv2   = sv2  + xmult*v2
  suv   = suv  + xmult*uv
  suv4  = suv4 + xmult*u2*u2*v2*v2
  suv2  = suv2 + xmult*u2*v2
  su3v  = su3v + xmult*uv*u2
  suv3  = suv3 + xmult*uv*v2
  
!       accumulate modified pair currents
  istate = k+iwfoff
  CALL sdenti(psi(0:nz,0:nr,1:2,istate),wrk1(0:nz,0:nr),uv*u2,istate)
  CALL sdenti(psi(0:nz,0:nr,1:2,istate),wrk2(0:nz,0:nr),uv*v2,istate)
  
END DO

!     first result: particle number
f1     = sv2-fnn

!     second result: the LN equation

!     but the numerator from Skyrme response is a bit special. It was 
!     computed with "v2" from previous iteration. The corresponding
!     "sum V2*U2" is stored in "suv2ra". In order to avoid
!     all to frequent recomputation of that expensive piece,
!     we rescale it here just by the new "suv2".
denom  = 32.0*(suv2*suv2-suv4)                   !  denominator
IF(ipair==11) THEN
  DO ir=0,nr; DO iz=0,nz
    wrk1(iz,ir)   = wrk1(iz,ir)*wrk2(iz,ir)*(1D0-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr) 
  END DO; END DO
ELSE
  wrk1 = wrk1*wrk2
END IF
enume1 = -vpair*2.0*f1int(wrk1)                  !  1. numerator
enume2 = suv2*suv2*erspln/(suv2ra*suv2ra)        !  2. numerator
f2     = xlam2+(enume1+enume2)/denom             ! p.n. fluctuation eq.

IF(iftest >= 3) WRITE(6,'(a,6g12.4)')  &
    ' suv2,suv4,suv2ra,erspln,vpair,echi=',  &
    suv2,suv4,suv2ra,erspln,vpair,f1int(wrk1(0:nz,0:nr))
DEALLOCATE(wrk1,wrk2)
RETURN
END SUBROUTINE lm2def


!---------------------------------------------------------------------------  
! DESCRIPTION: respax
!> @brief
!!Calculates the contribution from the Skyrme residual interaction to 
!!\c lambda_2 in Lipkin-Nogami Pairing.
!> @details
!!Computes the response-densities to the particle number operator using 
!!pairing weights from \c weight.
!!The response Hamiltonian is constructed from the Skyrme-Energy-functional, 
!!for the details see \cite Reinhard1992 \cite Reinhard1992a. 
!!The Lipkin-Nogami does only proton-proton pairing and neutron-neutron 
!!pairing. Thus only the Z^2 and N^2 response is required. 
!!The procedure is extensively commented for the first term \c xpotp and 
!!\c xpotn. The other terms run analogously.
!!Results are returned via
!!  - \c rspln = \f$ u^2v^2 \f$ response energy to be used in Lipk.Nogami.
!!  - \c suv2r = \f$ \sum u^2v^2 \f$, to be used in Lipkin-Nogami.
!>
!> @param[in] iter
!> INTEGER, takes the current iteration number. 
!--------------------------------------------------------------------------- 
SUBROUTINE respax(iter)

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
USE coulomb, ONLY: colrun,wcolt,irepc

IMPLICIT NONE

INTEGER, INTENT(IN) :: iter

REAL(8),PARAMETER :: uvlim=1D-8    ! lower limit for inclusion in sum

INTEGER :: iq,ir,irepcs,ist,iz,n
REAL(8) :: b3t1,b3t2,b3t3,b3t4,&
           rn,rp,rpow,rpowm1,rpowm2,rs,tmp1,twdthr,v2,xxx
REAL(8) :: xpotp,xkinp,xpd2p,xspop,xcol,xcex,xpotn,xkinn,xpd2n,xspon
           
REAL(8), ALLOCATABLE :: rrho(:,:,:)   !<Response density.
REAL(8), ALLOCATABLE :: rtau(:,:,:)   !<Response kinetic energy density.
REAL(8), ALLOCATABLE :: rrhols(:,:,:) !<Response spin-orbit current.
REAL(8), ALLOCATABLE :: rd2rho(:,:,:) !<Response laplacian on rho.
REAL(8), ALLOCATABLE :: rcoupt(:,:)   !<Response Coulomb potential.
REAL(8), ALLOCATABLE :: wgtres(:)     !<Workspace for response calculation.
REAL(8), ALLOCATABLE :: wrk1(:,:),wrk2(:,:)


ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
ALLOCATE(rrho(0:nz,0:nr,kqm), rtau(0:nz,0:nr,kqm), rrhols(0:nz,0:nr,kqm),  &
         rd2rho(0:nz,0:nr,kqm),rcoupt(0:nz,0:nr))
ALLOCATE(wgtres(npsi(2)))
wrk1=0D0;wrk2=0D0;rrho=0D0;rtau=0D0;rrhols=0D0;rd2rho=0D0;rcoupt=0D0

!---- compute the response weights
DO  n=1,nst(1)+nst(2)
  v2 = vsqf(n)
  wgtres(n) = MAX(2.0 * (v2-v2*v2),smal)
END DO

!---- calculate sum u^2 v^2 for use in bcslip
suv2rs(1)=sum(vsqf(1:nst(1))-vsqf(1:nst(1))**2)
suv2rs(2)=sum(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2)

suv2rs(1) = MAX(suv2rs(1),smal)
suv2rs(2) = MAX(suv2rs(2),smal)

!---- set response densities to 0D0 ------------------------------------
rrho  (0:nz,0:nr,1:2)=0D0
rtau  (0:nz,0:nr,1:2)=0D0
rrhols(0:nz,0:nr,1:2)=0D0
rd2rho(0:nz,0:nr,1:2)=0D0
rcoupt(0:nz,0:nr)    =0D0

!---- compute the response densities -------------------------------------
!     rrho   -->   < rho Delta N >
!     rtau   -->   < tau Delta N >
!     rrhols -->  < drhols Delta N >
DO ist = 1,nst(1)+nst(2),1
  iq = icharg(ist)
  IF(wgtres(ist) > uvlim) CALL densit(psi(0:nz,0:nr,1:2,ist),wgtres(ist),ist,  &
      rrho(0:nz,0:nr,iq),rtau(0:nz,0:nr,iq),rrhols(0:nz,0:nr,iq), &
      rholsr(0:nz,0:nr,iq),rholsz(0:nz,0:nr,iq) )
END DO

!---- rd2rho -->  Laplace -> < rho Delta N >
CALL laplace_scalar(rrho(0:nz,0:nr,1),rd2rho(0:nz,0:nr,1),+1,0,1,one)
CALL laplace_scalar(rrho(0:nz,0:nr,2),rd2rho(0:nz,0:nr,2),+1,0,1,one)

!---- terms depending on rrho * rrho
b3t1 = b3*(power+two)*(power+one)*edr
b3t2 = two*b3p*edr
b3t3 = two*power*b3t2
b3t4 = b3p*power*(power-one)*edr

DO ir=0,nr
  DO iz=0,nz
    rp = rho(iz,ir,1)         ! proton-density
    rn = rho(iz,ir,2)         ! neutron-density
    rs = rp+rn                ! total density
    rpowm2 = rs**(power-two)  !
    rpowm1 = rs**(power-one)
    rpow   = rs**(power)
    tmp1 = (b0-b0p) +rpow*(b3t1-b3t2)  &
        -rpowm2*b3t4*((rp**2)+(rn**2))
    wrk1(iz,ir) = rrho(iz,ir,1)*rrho(iz,ir,1) *(tmp1-rpowm1*b3t3*rp)
    wrk2(iz,ir) = rrho(iz,ir,2)*rrho(iz,ir,2) *(tmp1-rpowm1*b3t3*rn)
  END DO
END DO

IF (.NOT.tsym) THEN
  wrk1(0,0:nr) = 0D0
  wrk2(0,0:nr) = 0D0
END IF

xpotp = f1int(wrk1)
xpotn = f1int(wrk2)

!---- terms depending on rtau * rrho
wrk1(0:nz,0:nr) = rrho(0:nz,0:nr,1)*rtau(0:nz,0:nr,1)
wrk2(0:nz,0:nr) = rrho(0:nz,0:nr,2)*rtau(0:nz,0:nr,2)

xkinp = two*(b1-b1p)*f1int(wrk1(0:nz,0:nr))
xkinn = two*(b1-b1p)*f1int(wrk2(0:nz,0:nr))

!---- terms depending on rrho * (Laplace --> rrho)
wrk1(0:nz,0:nr) = rrho(0:nz,0:nr,1)*rd2rho(0:nz,0:nr,1)
wrk2(0:nz,0:nr) = rrho(0:nz,0:nr,2)*rd2rho(0:nz,0:nr,2)

xpd2p = (b2p-b2)*f1int(wrk1(0:nz,0:nr))
xpd2n = (b2p-b2)*f1int(wrk2(0:nz,0:nr))

!---- terms depending on rrho * rrhols
wrk1(0:nz,0:nr) = rrho(0:nz,0:nr,1)*rrhols(0:nz,0:nr,1)
wrk2(0:nz,0:nr) = rrho(0:nz,0:nr,2)*rrhols(0:nz,0:nr,2)

xspop = -two*(b4+b4p)*f1int(wrk1(0:nz,0:nr))
xspon = -two*(b4+b4p)*f1int(wrk2(0:nz,0:nr))

!---- Coulomb ---------------------------------------------------------------

irepcs = irepc ! has to be increased for response calculation.
irepc  = 200
!      icrflx = 1       ! flag for print of colrun properties

!---- save wcolt for use of colrun in gradstep
DO ir=0,nr
  DO iz=0,nz
    xxx = wcolt(iz,ir)
    wcolt(iz,ir) = wcoltt(iz,ir)
    wcoltt(iz,ir) = xxx
  END DO
END DO

!---- Coulomb-solver
CALL colrun(rrho(0:nz,0:nr,1),rcoupt(0:nz,0:nr),iter)

!---- save wcolt in wcoltt for use of colrun in respax
DO ir=0,nr
  DO iz=0,nz
    xxx = wcolt(iz,ir)
    wcolt(iz,ir) = wcoltt(iz,ir)
    wcoltt(iz,ir) = xxx
  END DO
END DO

wrk1(0:nz,0:nr) = rrho(0:nz,0:nr,1)*rcoupt(0:nz,0:nr)

xcol = f1int(wrk1(0:nz,0:nr))
irepc = irepcs ! reset irepc

!---- Coulomb-Exchange
twdthr = -(2.d0/3.d0)
wrk1(0:nz,0:nr) = rrho(0:nz,0:nr,1) *rrho(0:nz,0:nr,1) *(rho(0:nz,0:nr,1)**twdthr)

IF (.NOT.tsym) wrk1(0,0:nr) = 0.0D0

xcex = excfac*edr*f1int(wrk1(0:nz,0:nr))

!---- Summation of contributions -------------------------------------------
rspln(1) = xpotp+xkinp+xpd2p+xspop+xcol+xcex
rspln(2) = xpotn+xkinn+xpd2n+xspon

DEALLOCATE(wrk1,wrk2)
DEALLOCATE(rrho,rtau,rrhols,rd2rho,rcoupt)
DEALLOCATE(wgtres)
RETURN
END SUBROUTINE respax

!---------------------------------------------------------------------------  
! DESCRIPTION: rbrent
!> @brief
!!This subroutine is an adapted version of the <em> Van
!!Wijngaarden-Dekker-Brent</em> method for finding the root of a function
!!(see \cite Pre92aB ).
!>
!> @param[in] e
!> REAL(8), array, takes the single particle energies. 
!> @param[in] delta
!> REAL(8), array, takes the gap parameters. 
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[inout] gw
!> REAL(8), array, takes and returns occupation probabilities.
!> @param[in] pairgw_uv
!> REAL(8), array, takes the u*v cutoff weights complementing gap \c delta.
!> @param[in] pairgw_v2
!> REAL(8), array, takes the v**2 cutoff weights complementing gap \c delta.
!> @param[in] partnm
!> REAL(8), takes the desired particle number.
!> @param[in] ndim
!> INTEGER, takes actual number of s.p. states.
!> @param[in] nmax
!> INTEGER, takes maximum number of s.p. states.
!> @param[in] nblock
!> INTEGER, takes address of blocked state (none if 0).
!> @param[in] iab
!> INTEGER, takes print level; significant values are -2, -1, 0 and 1;
!! larger value prints more.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION rbrent(e,delta,ph,gw,pairgw_uv,pairgw_v2,partnm,ndim,nmax,nblock,iab)

USE paramr
IMPLICIT NONE

REAL(8), INTENT(IN)                  :: e(ndim)
REAL(8), INTENT(IN)                  :: delta(ndim)
REAL(8), INTENT(IN)                  :: ph(ndim)
REAL(8), INTENT(IN OUT)              :: gw(ndim)
REAL(8), INTENT(IN)                  :: pairgw_uv(ndim)
REAL(8), INTENT(IN)                  :: pairgw_v2(ndim)
REAL(8), INTENT(IN)                  :: partnm
INTEGER, INTENT(IN)                  :: ndim
INTEGER, INTENT(IN)                  :: nmax
INTEGER, INTENT(IN)                  :: nblock
INTEGER, INTENT(IN)                  :: iab

REAL(8),PARAMETER :: ra0=-1D2,rb0=1D2,reps=12D-14,rtolw=1D-16
INTEGER,PARAMETER :: iitmax=100

INTEGER :: ii
REAL(8) :: ra,rb,rc,rd,re,rfa,rfb,rfbo,rfc,rmin1,rmin2,rp,rq,rr,rs,rtol1,rxm

!     SETTING UP
ra  = ra0
rb  = rb0
rfa = partnm - parnmt(e, delta, ph, gw, pairgw_uv,pairgw_v2, ra, ndim, nmax,  &
    nblock)
rfb = partnm - parnmt(e, delta, ph, gw, pairgw_uv,pairgw_v2, rb, ndim, nmax,  &
    nblock)
rfc = rfb
rc  = rb
rd  = rb-ra
re  = rd

!     TESTING
IF ((rfa > 0.0).AND.(rfb > 0.0)) THEN
  IF (rfa > rfb) THEN
    rb = 1.0
  ELSE
    rb = 0.0
  END IF
ELSE IF ((rfa < 0.0).AND.(rfb < 0.0)) THEN
  IF (rfa > rfb) THEN
    rb = 0.0
  ELSE
    rb = 1.0
  END IF
ELSE
  
!     ITERATION
  iter: DO ii=1,iitmax
!     RENANME RA,RB,RC AND AJUST BOUNDING INTERVAL RD
    IF (((rfb > 0.0).AND.(rfc > 0.0)) .OR.((rfb < 0.0).AND.(rfc < 0.0))) THEN
      rc  = ra
      rfc = rfa
      rd  = rb-ra
      re  = rd
    END IF
    IF (ABS(rfc) < ABS(rfb)) THEN
      ra  = rb
      rb  = rc
      rc  = ra
      rfa = rfb
      rfb = rfc
      rfc = rfa
    END IF
!     CONVERGENCE CHECK
    rtol1 = 2.0*reps*ABS(rb) + 0.5*rtolw
    rxm   = 0.5*(rc-rb)
    IF ((ABS(rxm) <= rtol1).OR.((ABS(rfb) == 0.0))) EXIT iter
    IF ((ABS(re) >= rtol1).OR.(ABS(rfa) > ABS(rfb))) THEN
!     ATTEMPT INVERSE QUADRATIC INTERPOLATION
      rs=rfb / rfa
      IF (ra == rc) THEN
        rp = 2.0*rxm*rs
        rq = 1.0-rs
      ELSE
        rq = rfa / rfc
        rr = rfb / rfc
        rp = rs * (2.0*rxm*rq*(rq-rr)-(rb-ra)*(rr-1.0))
        rq = (rq-1.0)*(rr-1.0) *(rs-1.0)
      END IF
!     CHECK WHETHER IN BOUNDS
      IF (rp > 0.0) THEN
        rq = -rq
      END IF
      rp    = ABS(rp)
      rmin1 = 3.0*rxm*rq-ABS(rtol1*rq)
      rmin2 = ABS(rq*re)
      IF (2.0*rp < MIN(rmin1,rmin2)) THEN
!     ACCEPT INTERPOLATION
        re = rd
        rd = rp / rq
        
        IF(iab > 0) THEN
          WRITE(6,'(A)') ' QUADRATIC INTERPOLATION '
        END IF
        
      ELSE
!     INTERPOLATION FAILED -> BISECTION
        rd = rxm
        re = rd
        
        IF(iab > 0) THEN
          WRITE(6,'(A)') ' FAILED -> BISECTION '
        END IF
        
      END IF
!     BOUNDS DECREASING TOO SLOWLY -> BISECTION
    ELSE
      rd = rxm
      re = rd
      
      IF(iab > 0) THEN
        WRITE(6,'(A)') ' TOO SLOWLY -> BISECTION '
      END IF
      
    END IF
!     MOVE LAST BEST GUESS TO RA
    ra  = rb
    rfa = rfb
!     EVALUATE NEW TRIAL ROOT
    IF (ABS(rd) > rtol1) THEN
      rb = rb + rd
    ELSE
      rb = rb + SIGN(rtol1,rxm)
    END IF
    rfbo = rfb
    rfb  = partnm - parnmt( e, delta, ph, gw, pairgw_uv,pairgw_v2, rb,  &
        ndim,nmax,nblock)
    
    IF(iab > 0) THEN
      WRITE (6,'(I3,A,F28.18,A,F28.18,A,F28.18)')  &
          ii,' :  rb = ',rb,'  rfb = ',rfb,'  drfb = ', (rfbo-rfb)
    END IF
    IF (ii==iitmax) STOP 'no convergence in RBRENT'
  END DO iter
END IF
IF(iab > 0) THEN
  WRITE (6,'(A,F28.18)') 'EFERMI=NST= ', rb
END IF

rbrent = rb
RETURN
END FUNCTION rbrent

!---------------------------------------------------------------------------  
! DESCRIPTION: parnmt
!> @brief
!!Computes the particle number for given \c efermi 
!>
!> @param[in] e
!> REAL(8), array, takes the single particle energies. 
!> @param[in] delta
!> REAL(8), array, takes the gap parameters. 
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[inout] gw
!> REAL(8), array, takes and returns occupation probabilities.
!> @param[in] pairgw_uv
!> REAL(8), array, takes the u*v cutoff weights complementing gap \c delta.
!> @param[in] pairgw_v2
!> REAL(8), array, takes the v**2 cutoff weights complementing gap \c delta.
!> @param[in] efermi
!> REAL(8), takes the Fermi energy.
!> @param[in] ndim
!> INTEGER, takes actual number of s.p. states.
!> @param[in] nmax
!> INTEGER, takes maximum number of s.p. states.
!> @param[in] nblock
!> INTEGER, takes address of blocked state (none if 0).
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION parnmt(e, delta, ph, gw, pairgw_uv,pairgw_v2, efermi, ndim, nmax, nblock)

USE paramr
USE dynpar, ONLY: ipair

IMPLICIT NONE

REAL(8), INTENT(IN)            :: e(ndim)
REAL(8), INTENT(IN)            :: delta(ndim)
REAL(8), INTENT(IN)            :: ph(ndim)
REAL(8), INTENT(OUT)           :: gw(ndim)
REAL(8), INTENT(IN)            :: pairgw_uv(ndim)
REAL(8), INTENT(IN)            :: pairgw_v2(ndim)
REAL(8), INTENT(IN)            :: efermi
INTEGER, INTENT(IN)            :: ndim
INTEGER, INTENT(IN)            :: nmax
INTEGER, INTENT(IN)            :: nblock

REAL(8),PARAMETER :: mintemp=0.001D0
INTEGER :: na
REAL(8) :: delta2,edif,equasi

parnmt = 0.0

IF(ipair == 4) THEN            ! case of temperature
  DO  na=1,nmax
    equasi = (e(na)-efermi)/MAX(delta(na),mintemp)
    IF(equasi > 50.0d0) THEN
      gw(na) = 0.0
    ELSE IF(equasi < -50.0d0) THEN
      gw(na)  = 1.0
    ELSE
      gw(na)  = 1.0/(1.0+EXP(equasi))
    END IF
    parnmt    = gw(na)*ph(na)    + parnmt
  END DO
ELSE
  DO na = 1, nmax
    delta2 = delta(na)*pairgw_uv(na)
    delta2 = delta2*delta2
    edif   = e(na)-efermi
    equasi = SQRT(edif*edif+delta2)
    IF(na == nblock) THEN
      gw(na) = 0.5D0
    ELSE
      gw(na) = 0.5D0  * (1.0D0 - edif/equasi)
      gw(na) = MAX(gw(na),smal)
      gw(na) = MIN(gw(na),one-smal)
    END IF
    parnmt  = ph(na)*gw(na)*pairgw_v2(na) + parnmt
  END DO
END IF

RETURN
END FUNCTION parnmt


SUBROUTINE thermal_occ(t,e,gw,ph,pairgw_uv,pairgw_v2,efermi,parnum,ndim)
IMPLICIT NONE
REAL(8), INTENT(IN)                  :: t
REAL(8), INTENT(IN)                  :: e(ndim)
REAL(8), INTENT(IN OUT)              :: gw(ndim)
REAL(8), INTENT(IN)                  :: ph(ndim)
REAL(8), INTENT(IN)                  :: pairgw_uv(ndim)
REAL(8), INTENT(IN)                  :: pairgw_v2(ndim)
REAL(8), INTENT(INOUT)               :: efermi
REAL(8), INTENT(IN)                  :: parnum
INTEGER, INTENT(IN)                  :: ndim

REAL(8) :: tarr(ndim),pnum

  !eferm = rbrent(e, delta, ph, gw, pairgw_uv,pairgw_v2, partnm, ndim, nmax, nblock, iab)
  !parnum = parnmt(e, delta, ph, gw, pairgw_uv,pairgw_v2, eferm, ndim, nmax, nblock)
  tarr=t
  efermi = rbrent(e, tarr, ph, gw, pairgw_uv,pairgw_v2, parnum, ndim, ndim, 0, 0)
  pnum = parnmt(e, tarr, ph, gw, pairgw_uv,pairgw_v2, efermi, ndim, ndim, 0)

END SUBROUTINE thermal_occ

!NEW ODD
SUBROUTINE findconjpart(iq)

!> @brief
!!This subroutine determines the pairing conjugate partner states.

USE dynpar, ONLY: toddterms
USE wavefr, ONLY: psi
USE levels
USE grid, ONLY: nz,nr,f1int

IMPLICIT NONE

INTEGER, INTENT(IN)  :: iq


REAL(8),ALLOCATABLE :: wrk1(:,:)
INTEGER,ALLOCATABLE :: itpartnsav(:)
INTEGER :: istart,i1,i2
REAL(8) :: supnorm,a
LOGICAL :: tpartn

LOGICAL, PARAMETER :: toverlap=.true.  ! energy or overlap decision

IF(toddterms) THEN
  tpartn = SUM(itpartn)>0
  IF(tpartn) THEN
    ALLOCATE(itpartnsav(nst(1)+nst(2)))
    itpartnsav=itpartn
  END IF
  ALLOCATE(wrk1(0:nz,0:nr))
!     search the time-conjugate partners
  IF(iq == 1)THEN
    istart= 1
  ELSE
    istart=nst(1)+1
  END IF
  DO i1=istart,nst(iq)+istart-1
    IF(toverlap) THEN
      supnorm = 0D0
    ELSE
      supnorm = 1D30
    END IF
    DO i2=istart,nst(iq)+istart-1
      IF((mpsi(i1)+mpsi(i2)) == -1)THEN
        IF(toverlap) THEN
          wrk1(0:nz,0:nr)=psi(0:nz,0:nr,1,i1)*psi(0:nz,0:nr,2,i2)-  &
                          psi(0:nz,0:nr,2,i1)*psi(0:nz,0:nr,1,i2)
          a = ABS(f1int(wrk1))
!         IF(i1==istblo(iq)) WRITE(7,'(a,2i4,3(1pg13.5))')
!     &    'FCONJG:',i1,i2,spenrg(i1),spenrg(i2),a
          IF(supnorm < a) THEN
            supnorm = a
            itpartn(i1) = i2
            rtpartn(i1) = supnorm
          END IF
        ELSE
          a = ABS(spenrg(i1)-spenrg(i2))
          IF(supnorm > a) THEN
            supnorm = a
            itpartn(i1) = i2
            rtpartn(i1) = supnorm
          END IF
        END IF
      END IF
    END DO
  END DO
  DO i1=istart,nst(iq)+istart-1
    IF(i1 /= itpartn(itpartn(i1)).AND. (i1==istblo(1).OR.i1==istblo(2)))  &
        WRITE(7,'(a,3i5,2f5.1,2f9.3)')  &
        'FINDCONJ error: state=',i1,itpartn(i1),itpartn(itpartn(i1)),  &
        mpsi(i1)+0.5,mpsi(itpartn(i1))+0.5, spenrg(i1),spenrg(itpartn(i1))
    IF(tpartn .AND. itpartn(i1).NE.itpartnsav(i1)) WRITE(7,'(a,2i4)') &
        'FINDCONJ change: new,old=',itpartn(i1),itpartnsav(i1)
  END DO
  DEALLOCATE(wrk1)
  IF(tpartn) DEALLOCATE(itpartnsav)
ELSE
  itpartn = 0
END IF

RETURN
END SUBROUTINE findconjpart

SUBROUTINE combine_spenergies(spe_in,spe_out,nsttot)

!> @brief
!!This subroutine builds the average of s.p. energies between
!!conjugate states.

USE dynpar, ONLY: toddterms
USE levels, ONLY: itpartn

IMPLICIT NONE

INTEGER, INTENT(IN)   :: nsttot
REAL(8), INTENT(IN)   :: spe_in(nsttot)
REAL(8), INTENT(OUT)  :: spe_out(nsttot)

INTEGER :: ist

IF(.NOT.toddterms) THEN
  spe_out=spe_in
  RETURN
END IF

DO ist=1,nsttot
  spe_out(ist) = 0.5D0*(spe_in(ist)+spe_in(itpartn(ist)))
END DO

END SUBROUTINE combine_spenergies

!END NEW ODD


END MODULE


