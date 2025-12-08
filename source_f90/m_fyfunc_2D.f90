MODULE fyfunc

USE paramr, ONLY:pi,e2,xxxl,ixxxl
USE dynpar

IMPLICIT NONE

! Fayans model parameters are defined in module dynpar

! auxiliary variables deduced from given parameters
REAL(8) :: fxip
REAL(8) :: h1xip
REAL(8) :: hgradxip
REAL(8) :: fxin
REAL(8) :: h1xin
REAL(8) :: hgradxin
REAL(8) :: r0fy
REAL(8) :: r0fy2
REAL(8) :: rho0inv
REAL(8) :: efermfy
REAL(8) :: efacfy
REAL(8) :: ufacfy
REAL(8) :: pfacfy
REAL(8) :: clsp
REAL(8) :: clsm
REAL(8) :: ctensp
REAL(8) :: ctensm
REAL(8) :: ufacfypr
REAL(8) :: pgradfac
REAL(8) :: c0fac
REAL(8) :: xinfac
REAL(8) :: xipfac


REAL(8),ALLOCATABLE,PRIVATE :: psip(:,:,:),psin(:,:,:)
REAL(8),ALLOCATABLE,PRIVATE :: wrk1(:,:),wrk2(:,:),wrk3(:,:),wrk4(:,:)


CONTAINS

!-----ini_fyfunc------------ part of the hartree-fock package -----------------

SUBROUTINE ini_fyfunc()

!USE dynpar, ONLY:h2m,tpair

REAL(8) :: h2mav

!      IMPLICIT NONE
!     Initialize auxiliary parameters for Fayans functional.



!------------------------------------------------------------------------------

IF(.NOT.tFAYANS) STOP "INI_FYFUNC needs TFAYANS=.TRUE."

! resolve isovectormodel

IF(isovectormodel==0) THEN
  fxip=fxi
  h1xip=h1xi
  hgradxip=hgradxi
  fxin=fxi
  h1xin=h1xi
  hgradxin=hgradxi
ELSE IF(isovectormodel==2) THEN
  fxip=fxi+delfxi
  h1xip=h1xi
  hgradxip=hgradxi
  fxin=fxi-delfxi
  h1xin=h1xi
  hgradxin=hgradxi
ELSE IF(isovectormodel==3) THEN
  fxip=fxi+delfxi
  h1xip=h1xi+delh1xi
  hgradxip=hgradxi+delhgradxi
  fxin=fxi-delfxi
  h1xin=h1xi-delh1xi
  hgradxin=hgradxi-delhgradxi
ELSE IF(isovectormodel==1) THEN
  xipfac=(fxi+delfxi)/fxi
  fxip=fxi*xipfac
  h1xip=h1xi*xipfac
  hgradxip=hgradxi*xipfac
  xinfac=(fxi-delfxi)/fxi
  fxin=fxi*xinfac
  h1xin=h1xi*xinfac
  hgradxin=hgradxi*xinfac
END IF


!     auxiliary variables

!          general
h2mav=(h2m(1)+h2m(2))/2D0
r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
r0fy2 = r0fy*r0fy
efermfy = (1.125D0*pi)**(2D0/3D0)*h2mav/r0fy2
efacfy = efermfy*rho0fy*(2D0/3D0)
ufacfy = efermfy/3D0
pfacfy = 2D0*efermfy/(3D0*rho0fy)
rho0inv = 0.5D0/rho0fy
!         spin-orbit
c0fac = 2D0*efermfy/(3D0*rho0fy)*r0fy2
clsp = c0fac*c0nabj
clsm = c0fac*c1nabj
ctensp = c0fac*c0tens
ctensm = c0fac*c1tens
!         pairing
!      ufacfypr = pfacfy/(2D0*rho0xi)
ufacfypr = pfacfy/(rho0xi)     ! ???
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
pgradfac = -2D0*r0fy2*pfacfy/rho0xi
!      IF(ipair.NE.7 .AND. ipair.NE.8)
!     &   STOP 'Fy functional requires IPAIR=7 or 8'
!tpair = .true.

WRITE(6,'(a,8(1pg16.7))') 'rho0fy,efacfy,efermfy,pfacfy=',  &
    rho0fy,efacfy,efermfy,pfacfy,ctensp,ctensm,clsp,clsm
CALL flush(6)

RETURN
END SUBROUTINE ini_fyfunc

!-----fyfunc------------ part of the hartree-fock package ---------------------

SUBROUTINE fyfunc_local(rp,rn,rc,d1pr,d1pz,d1nr,d1nz,d2p,d2n,tp,tn,  &
                  sp,sn,ssp,ssn, potp,potn,bmasp,bmasn,  &
                  wlspr,wlspz,wlsnr,wlsnz,energd)

!     the basic Fayans energy density and derived quantities.

USE dynpar, ONLY:chgx,excfac

!      IMPLICIT NONE
!INCLUDE 'paramr.inc'

!INCLUDE 'dynpar.inc'
!INCLUDE 'wavefr.inc'


! local densities at a given point (input)
REAL(8), INTENT(IN)        :: rp  !< proton density
REAL(8), INTENT(IN)        :: rn  !< neutron density
REAL(8), INTENT(IN)        :: rc  !< charge density
REAL(8), INTENT(IN)        :: d1pr !< r-gradient of proton density
REAL(8), INTENT(IN)        :: d1pz !< z-gradient of proton density
REAL(8), INTENT(IN)        :: d1nr !< r-gradient of neutron density
REAL(8), INTENT(IN)        :: d1nz !< z-gradient of neutron density
REAL(8), INTENT(IN)        :: d2p !< Laplacian of proton density
REAL(8), INTENT(IN)        :: d2n !< Laplacian of neutron density
REAL(8), INTENT(IN)        :: tp  !< proton kinetic density
REAL(8), INTENT(IN)        :: tn  !< neutron kinetic density
REAL(8), INTENT(IN)        :: sp  !< proton spin-orbit density
REAL(8), INTENT(IN)        :: sn  !< neutron spin-orbit density
REAL(8), INTENT(IN)        :: ssp  !< DIV of proton spin-orbit density
REAL(8), INTENT(IN)        :: ssn  !< DIV of neutron spin-orbit density


! corresponding local potentials (output)
REAL(8), INTENT(OUT)      :: potp  !< proton potential
REAL(8), INTENT(OUT)      :: potn  !< neutron potential
REAL(8), INTENT(OUT)      :: bmasp  !< proton inverse mass
REAL(8), INTENT(OUT)      :: bmasn  !< neutron inverse mass
REAL(8), INTENT(OUT)      :: wlspr  !< proton l*s potential r-component
REAL(8), INTENT(OUT)      :: wlspz  !< proton l*s potential z-component
REAL(8), INTENT(OUT)      :: wlsnr  !< neutron l*s potential r-component
REAL(8), INTENT(OUT)      :: wlsnz  !< neutron l*s potential z-component
REAL(8), INTENT(OUT)      :: energd  !< energy density



REAL(8),PARAMETER :: dsmall=1D-30
LOGICAL, SAVE :: tfirst=.true.

REAL(8) :: rs,ts,ss,sss,xplus,xminus,xpsigma,xmsigma,d1sr2,d1mr2,add,add0
REAL(8) :: d1sr,d1sz,d2s,hxpsigmac

!------------------------------------------------------------------------------

IF(.NOT.tFAYANS) STOP "FYFUNC needs TFAYANS=.TRUE."


!     total densities and other auxiliary quantities

rs  = rp+rn
ts  = tp+tn
ss  = sp+sn
sss = ssp+ssn
d2s = d2p+d2n
d1sr = d1pr+d1nr
d1sz = d1pz+d1nz
xplus = rs*rho0inv
xminus = (rp-rn)*rho0inv
xpsigma = MAX(xplus,dsmall)**sigmavol
xmsigma = xminus  !**sigmavol
d1sr2 = (d1sr*d1sr+d1sz*d1sz)*r0fy2*rho0inv**2
d1mr2 = d1sr2 ! (d1p-d1n)**2*r0fy2*rho0inv**2
!      energ1 =  efacfy*(avp*(xplus*xplus)*
!     &   (1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma))

!     energy  (direct Coulomb separately, pairing part in HFBpot)

energd = efacfy*(avp*(xplus*xplus)* (1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    +avm*(xminus*xminus)* (1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    +asp*d1sr2/(1D0+hsp*xpsigma+hgradp*d1sr2)  &
    +asm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2) )  &
    +c0tau*rs*ts+c1tau*(rp-rn)*(tp-tn)

!     mean field potentials (surface term needs extra loop later)
!                           (direct Coulomb separately, pairing part in HFBpot)

potp = ufacfy*( 2D0*avp*xplus*(1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    -avp*sigmavol*xpsigma*xplus*(h1vp+h2vp)/ (1D0+h2vp*xpsigma)**2  &
    +2D0*avm*xminus*(1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    -avm*xminus*xminus*(h1vm+h2vm)/(1D0+h2vm*xplus)**2  &
    -asp*hsp*sigmavol*d1sr2*xpsigma/MAX(xplus,dsmall)/  &
    (1D0+hsp*xpsigma+hgradp*d1sr2)**2  &
    -asm*hsm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)**2 )  &
    +c0tau*ts+c1tau*(tp-tn)
!     &       +c0nabj*sss+c1nabj*(ssp-ssn)
potn = ufacfy*( 2D0*avp*xplus*(1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    -avp*sigmavol*xpsigma*xplus*(h1vp+h2vp)/ (1D0+h2vp*xpsigma)**2  &
    -2D0*avm*xminus*(1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    -avm*xminus*xminus*(h1vm+h2vm)/(1D0+h2vm*xplus)**2  &
    -asp*hsp*sigmavol*d1sr2*xpsigma/MAX(xplus,dsmall)/  &
    (1D0+hsp*xpsigma+hgradp*d1sr2)**2  &
    -asm*hsm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2)**2 )  &
    +c0tau*ts-c1tau*(tp-tn)
!     &       +c0nabj*sss-c1nabj*(ssp-ssn)
bmasp = c0tau*rs+c1tau*(rp-rn)
bmasn = c0tau*rs-c1tau*(rp-rn)

!     Coulomb exchange

IF(chgx /= 0D0) THEN
  xplus=MAX(xplus,1D-90)
  excfac = 1.4179979
  hxpsigmac = hcoul*xplus**sigmacoul
  add0 = -excfac*MAX(rc,1D-99)**0.333333333333D0
  add = add0*(1D0-hxpsigmac)
  energd = 0.75D0*add*rc + energd
  potp = add - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus + potp
  potn =  - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus + potn
END IF

!     spin-orbit

potp   = -clsp*sss - clsm*(ssp-ssn)   + potp
potn   = -clsp*sss + clsm*(ssp-ssn)   + potn
wlspr   = (clsp*d1sr+clsm*(d1pr-d1nr)) +ctensp*ss+ctensm*(sp-sn)
wlsnr   = (clsp*d1sr-clsm*(d1pr-d1nr)) +ctensp*ss-ctensm*(sp-sn)
wlspz   = (clsp*d1sz+clsm*(d1pz-d1nz)) +ctensp*ss+ctensm*(sp-sn)
wlsnz   = (clsp*d1sz-clsm*(d1pz-d1nz)) +ctensp*ss-ctensm*(sp-sn)
energd = -clsp*rs*sss - clsm*(rp-rn)*(ssp-ssn)  &
         +ctensp*ss*ss+ctensm*(sp-sn)**2 + energd


! a testprint, only once
IF(tfirst) THEN
  WRITE(6,'(a,6(1pg13.5))') 'r0fy,efermfy,C_0,ufacfy,rhoinv=',  &
      r0fy,efermfy,2D0*efermfy/(3D0*rho0fy),ufacfy,rho0inv,asp
  WRITE(6,'(a,4(1pg13.5))') 'potp etc=', rp,rn,xpsigma,potp,potn,energd
  WRITE(6,*) 'clsp,c0nabj,clsm,c1nabj=',clsp,c0nabj,clsm,c1nabj
  tfirst=.false.
END IF


RETURN
END SUBROUTINE fyfunc_local

!-----fyfunc_grad------------ part of the hartree-fock package -----------------

SUBROUTINE fyfunc_grad()

!IMPLICIT NONE

USE paramr, ONLY:kqm
USE grid, ONLY:nz,nr,nzd,nrd,nablaz_scalar,nablar_scalar
USE fields, ONLY:pot,rho,d2rho

!     The derivative contribution to the gradient term of the Fayans functional-
!     The result is accumulated onto the array 'pot'.
!     The routine must be called after 'fyfunc' having filled the array 'pot'.


REAL(8), ALLOCATABLE :: dpotaux(:,:,:),gpotadd(:,:,:,:)

REAL(8) :: rp,rn,rs,xminus,xmsigma,xplus,xpsigma,d1rp,d1zp,d1rn,d1zn,d1mr2,d1sr2
REAL(8) :: d1mr,d1sr,d1mz,d2mr,d1sz,d2sr,d2n,d2p
INTEGER :: iz,ir

!------------------------------------------------------------------------------

ALLOCATE(dpotaux(0:nz,0:nr,kqm),gpotadd(0:nz,0:nr,kqm,2))
ALLOCATE(psin(0:nz,0:nr,kqm),psip(0:nz,0:nr,kqm))

!     gradients of density, juse 'gpotadd' as workspace
CALL nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,1D0)  ! d1z proton
CALL nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,1D0)  ! d1z neutron
CALL nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,1D0)  ! d1r proton
CALL nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,1D0)  ! d1r neutron

DO iz=0,nz; DO ir=0,nr
  rp  = rho(iz,ir,1)
  rn  = rho(iz,ir,2)
  xminus = (rp-rn)*rho0inv
  xmsigma = xminus !**sigmavol
  rs  = rp+rn
  xplus = rs*rho0inv
  xpsigma = xplus**sigmavol
  d1rp = psip(iz,ir,1)
  d1zp = psin(iz,ir,1)
  d1rn = psip(iz,ir,2)
  d1zn = psin(iz,ir,2)
  d1mr2 = ((d1rp-d1rn)**2+(d1zp-d1zn)**2)*r0fy2*rho0inv**2
  d1sr2 = ((d1rp+d1rn)**2+(d1zp+d1zn)**2)*r0fy2*rho0inv**2
  dpotaux(iz,ir,1) = ufacfy*asp*r0fy2*(1D0+hsp*xpsigma)/  &
      (1D0+hsp*xpsigma+hgradp*d1sr2)**2
  dpotaux(iz,ir,2) = ufacfy*asm*r0fy2*(1D0+hsm*xplus)/  &
      (1D0+hsm*xplus+hgradm*d1sr2)**2
!        IF(iz+ir<5) WRITE(6,*) iz,ir,ufacfy,asp,r0fy2,dpotaux(iz,ir,1)
!
!          --> pairing contribution still missing
!
!        WRITE(6,*) iz,ir,dpotaux(iz,ir,1)
END DO; END DO

CALL nablaz_scalar(dpotaux(:,:,1),gpotadd(:,:,1,1), 1,1,1D0)
CALL nablaz_scalar(dpotaux(:,:,2),gpotadd(:,:,2,1), 1,1,1D0)
CALL nablar_scalar(dpotaux(:,:,1),gpotadd(:,:,1,2), 0,1,1D0)
CALL nablar_scalar(dpotaux(:,:,2),gpotadd(:,:,2,2), 0,1,1D0)

DO iz=0,nz; DO ir=0,nr
  
  d1rp = psip(iz,ir,1)
  d1zp = psin(iz,ir,1)
  d1rn = psip(iz,ir,2)
  d1zn = psin(iz,ir,2)
  d2p = d2rho(iz,ir,1)
  d2n = d2rho(iz,ir,2)
  
  d1mr = (d1rp-d1rn)*rho0inv
  d1sr = (d1rp+d1rn)*rho0inv
  d1mz = (d1zp-d1zn)*rho0inv
  d1sz = (d1zp+d1zn)*rho0inv
  d2mr = (d2p-d2n)*rho0inv
  d2sr = (d2p+d2n)*rho0inv
  
  pot(iz,ir,1) = pot(iz,ir,1)  &
      -2D0*(gpotadd(iz,ir,1,2)*d1sr+gpotadd(iz,ir,1,1)*d1sz  &
      +dpotaux(iz,ir,1)*d2sr)  &
      -2D0*(gpotadd(iz,ir,2,2)*d1mr+gpotadd(iz,ir,2,1)*d1mz  &
      +dpotaux(iz,ir,2)*d2mr)
!     &          -2D0*(gpotadd(iz,ir,1,2)*d1mr+gpotadd(iz,ir,1,1)*d1mz
!     &                +dpotaux(iz,ir,2)*d2mr)
  pot(iz,ir,2) = pot(iz,ir,2)  &
      -2D0*(gpotadd(iz,ir,1,2)*d1sr+gpotadd(iz,ir,1,1)*d1sz  &
      +dpotaux(iz,ir,1)*d2sr)  &
      +2D0*(gpotadd(iz,ir,2,2)*d1mr+gpotadd(iz,ir,2,1)*d1mz  &
      +dpotaux(iz,ir,2)*d2mr)
!     &          +2D0*(gpotadd(iz,ir,1,2)*d1mr+gpotadd(iz,ir,1,1)*d1mz
!     &                +dpotaux(iz,ir,2)*d2mr)
!        IF(iz+ir<5) WRITE(6,*) iz,ir,dpotaux(iz,ir,1),
!     &          -2D0*(gpotadd(iz,ir,1,2)*d1sr+gpotadd(iz,ir,1,1)*d1sz
!     &                +dpotaux(iz,ir,1)*d2sr)
END DO; END DO

DEALLOCATE(dpotaux,gpotadd,psin,psip)

RETURN
END SUBROUTINE fyfunc_grad

!-----fy_hfbpot------------ part of the hartree-fock package -----------------

SUBROUTINE fy_hfbpot()

USE paramr, ONLY:kqm
USE grid, ONLY:nr,nz,dr,dz,sz,sr,icm,nzd,nrd,f1int,nablaz_scalar,nablar_scalar,test_print
USE fields, ONLY:rho,paircr,pairpt
USE dynpar, ONLY:epair

!      IMPLICIT NONE
!INCLUDE 'paramr.inc'

!     The pairing potential of the Fayans functional.
!     Also used to compute pairing energy.

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

!INCLUDE 'dynpar.inc'
!      include 'numpar.inc'
!INCLUDE 'fields.inc'
!INCLUDE 'wavefr.inc'
!INCLUDE 'observ.inc'
!      include 'grid.inc'
!      include 'constn.inc'
REAL(8) :: d1sr,d1sz,d1sr2,xplus,xplusgamma,rhfacp,rhfacn,z,r
INTEGER :: ir,iz
LOGICAL,PARAMETER :: ttest=.FALSE.

!------------------------------------------------------------------------------

ALLOCATE(psin(0:nz,0:nr,kqm),psip(0:nz,0:nr,kqm))
ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))

IF(ttest) THEN
  DO ir=0,nr
    r     = dr*(ir+sr)
    DO iz=0,nz
      z     = ((iz)+sz-icm)*dz
      rho(iz,ir,1) = EXP(-(r**2+z**2)/9D0)
    END DO
  END DO
  rho(:,:,2) = rho(:,:,1)
  paircr(:,:,1)  = rho(:,:,1)
  paircr(:,:,2)  = rho(:,:,1)
  pairpt = 0D0
!  WRITE(6,*) 'FY_HFBPOT'
END IF
!        r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!        r0fy2 = r0fy*r0fy
!        efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      pfacfy = 2D0*efermfy/(3D0*rho0fy)
!      ufacfypr = pfacfy/(2D0*rho0xi)
!     gradients of density     <-- recycled from FY_PAIRPOT
CALL nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,1D0)  ! d1z proton
CALL nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,1D0)  ! d1z neutron
CALL nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,1D0)  ! d1r proton
CALL nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,1D0)  ! d1r neutron
!
IF(ttest) THEN
  ALLOCATE(wrk3(0:nz,0:nr),wrk4(0:nz,0:nr))
  wrk4=0D0
  wrk3=0D0
END IF
DO iz=0,nz; DO ir=0,nr
  d1sr = (psip(iz,ir,1)+psip(iz,ir,2))  !/rho0xi
  d1sz = (psin(iz,ir,1)+psin(iz,ir,2))  !/rho0xi
  d1sr2 = (d1sr*r0fy/rho0xi)**2+(d1sz*r0fy/rho0xi)**2
  xplus = (rho(iz,ir,1)+rho(iz,ir,2))/rho0xi
  xplusgamma = xplus**gammaxi
!        rhfac = pfacfy*(fxi+h1xi*xplusgamma+hgradxi*d1sr2)
  rhfacp = pfacfy*(fxip+h1xip*xplusgamma+hgradxip*d1sr2)
  rhfacn = pfacfy*(fxin+h1xin*xplusgamma+hgradxin*d1sr2)
  pairpt(iz,ir,1) = paircr(iz,ir,1) * rhfacp
  pairpt(iz,ir,2) = paircr(iz,ir,2) * rhfacn
!        paircrtot2 = paircr(i)**2+paircr(in)**2
!       epair = epair + 0.5D0*rhfac*paircrtot2*(ract*ract)
  wrk1(iz,ir) = (rhfacp*paircr(iz,ir,1)**2+rhfacn*paircr(iz,ir,2)**2)
  wrk2(iz,ir) =  rhfacn*paircr(iz,ir,2)**2
!  IF(ttest) THEN
!    wrk3(iz,ir) =  xplusgamma
!    wrk4(iz,ir) =  d1sr2
!  END IF
END DO; END DO
epair(2) = f1int(wrk2)
epair(1) = f1int(wrk1)-epair(2)
!      IF(nprot==20 .AND. nneut==20) THEN
!        WRITE(6,'(a,4(1pg13.5))') 'FYpair:',
!     &      pfacfy*fxip,pfacfy*h1xip,pfacfy*hgradxip,epair
!      END IF
IF(ttest) WRITE(6,'(a,5(1pg13.5))') 'FYpair:',epair, f1int(wrk3),f1int(wrk4),pfacfy
IF(ttest) THEN
  WRITE(7,'(a,4(1pg13.5))')'params:',pfacfy,fxin,h1xin,hgradxin
  CALL test_print(paircr(:,:,2),0,9,0,9,'paircr')
  CALL test_print(pairpt(:,:,2),0,9,0,9,'pairpt')
  CALL test_print(rho(:,:,1),0,9,0,9,'rho(:,:,1)')
  CALL test_print(wrk3(:,:),0,9,0,9,'xplusgamm')
  CALL test_print(wrk4,0,9,0,9,'d1sr2')
  CALL test_print(wrk2,0,9,0,9,'rhfacn')
  STOP
END IF
!        write(6,'(a)') ' Fy pair potential for neutrons:'
!        write(6,'(5(1pg12.4))') (pairpt(i),i=ngrid+1,2*ngrid)

DEALLOCATE(psin,psip)
DEALLOCATE(wrk1,wrk2)

RETURN
END SUBROUTINE fy_hfbpot
!-----fy_pairpot------------ part of the hartree-fock package -----------------

SUBROUTINE fy_pairpot()

!      IMPLICIT NONE
!INCLUDE 'paramr.inc'

USE paramr, ONLY:kqm
USE fields, ONLY:paircr,rho,pot,d2rho
USE grid, ONLY:nr,nz,dr,dz,sz,sr,icm,nzd,nrd,nablaz_scalar,nablar_scalar,test_print,laplace_scalar

!     Non-gradient part of the  pairing contribnution
!     to potential from the Fayans functional.

!     The result is accumulated onto the array 'pot'.

!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'
!INCLUDE 'wavefr.inc'
!INCLUDE 'derivs.inc'
!      include 'constn.inc'
REAL(8), ALLOCATABLE :: dpotauxp(:,:),gpotaddp(:,:,:)
REAL(8), ALLOCATABLE :: dpotauxn(:,:),gpotaddn(:,:,:)
REAL(8) :: d1sr,d1sz,d1sr2,d2sr,xplus,xplusgamma,add,r,z
INTEGER :: ir,iz
LOGICAL, PARAMETER :: ttest=.false.

!------------------------------------------------------------------------------

ALLOCATE(psin(0:nz,0:nr,kqm),psip(0:nz,0:nr,kqm))

IF(ttest) THEN
  DO ir=0,nr
    r     = dr*(ir+sr)
    DO iz=0,nz
      z     = ((iz)+sz-icm)*dz
      rho(iz,ir,1) = EXP(-(r**2+z**2)/9D0)
    END DO
  END DO
  rho(:,:,2) = rho(:,:,1)
  paircr(:,:,1)  = rho(:,:,1)
  paircr(:,:,2)  = rho(:,:,1)
  pot = 0D0
END IF
ALLOCATE(dpotauxp(0:nz,0:nr),gpotaddp(0:nz,0:nr,2))
ALLOCATE(dpotauxn(0:nz,0:nr),gpotaddn(0:nz,0:nr,2))
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
dpotauxp(0:nz,0:nr) = paircr(0:nz,0:nr,1)**2
dpotauxn(0:nz,0:nr) = paircr(0:nz,0:nr,2)**2
!      WRITE(6,*) 'paircr**2:',SUM(dpotaux(1:ngrid))
CALL nablaz_scalar(dpotauxp,gpotaddp(:,:,1), 1,1,1D0)
CALL nablaz_scalar(dpotauxn,gpotaddn(:,:,1), 1,1,1D0)
CALL nablar_scalar(dpotauxp,gpotaddp(:,:,2), 0,1,1D0)
CALL nablar_scalar(dpotauxn,gpotaddn(:,:,2), 0,1,1D0)
!     gradients of density
CALL nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,1D0)  ! d1z proton
CALL nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,1D0)  ! d1z neutron
CALL nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,1D0)  ! d1r proton
CALL nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,1D0)  ! d1z neutron
IF(ttest) THEN
  CALL laplace_scalar(rho(0,0,1),d2rho(0,0,1),+1,0,1,1D0)
  CALL laplace_scalar(rho(0,0,2),d2rho(0,0,2),+1,0,1,1D0)
END IF

!      WRITE(6,*) 'FY_PAIRPOT test pairing potential:',
!     &   ufacfypr*h1xi*gammaxi,pgradfac
DO iz=0,nz; DO ir=0,nr
  d1sr = (psip(iz,ir,1)+psip(iz,ir,2))/rho0xi
  d1sz = (psin(iz,ir,1)+psin(iz,ir,2))/rho0xi
  d1sr2 = (d1sr*r0fy)**2+(d1sz*r0fy)**2
  d2sr = (d2rho(iz,ir,1)+d2rho(iz,ir,2))/rho0xi
  xplus = (rho(iz,ir,1)+rho(iz,ir,2))/rho0xi
  IF(gammaxi==1D0) THEN
    xplusgamma = 1D0
  ELSE
!          WRITE(*,*) rho0xi,gammaxi,xplus
    xplusgamma = MAX(xplus,1D-90)**(gammaxi-1D0)
  END IF
  IF(ttest) THEN
    ALLOCATE(wrk2(0:nz,0:nr),wrk3(0:nz,0:nr),wrk4(0:nz,0:nr))
    wrk3(iz,ir) =  xplusgamma
    wrk4(iz,ir) =  d1sr
    wrk2(iz,ir) =  d2sr
  END IF
  add  = ufacfypr*gammaxi*xplusgamma*  &
      (dpotauxp(iz,ir)*h1xip+dpotauxn(iz,ir)*h1xin) +pgradfac*(  &
      hgradxip*(dpotauxp(iz,ir)*d2sr+gpotaddp(iz,ir,2)*d1sr  &
      +gpotaddp(iz,ir,1)*d1sz)  &
      +hgradxin*(dpotauxn(iz,ir)*d2sr+gpotaddn(iz,ir,2)*d1sr  &
      +gpotaddn(iz,ir,1)*d1sr) )
!         IF(5<iz+ir.AND.iz+ir<15) WRITE(*,*) 'add:',iz,ir,add
!        WRITE(6,*) i,add,pot(i),pot(j)
!        add=0D0
  pot(iz,ir,1) = pot(iz,ir,1) + add
  pot(iz,ir,2) = pot(iz,ir,2) + add
END DO; END DO
IF(ttest) THEN
  WRITE(7,'(a,4(1pg13.5))')'params:',ufacfypr,pgradfac,hgradxin
  CALL test_print(paircr(:,:,2),0,9,0,9,'paircr')
  CALL test_print(pot(:,:,2),0,9,0,9,'add pot')
  CALL test_print(wrk3,0,9,0,9,'xplusgamm')
  CALL test_print(wrk4,0,9,0,9,'d1sr')
  CALL test_print(wrk2,0,9,0,9,'d2sr')
  DEALLOCATE(wrk2,wrk3,wrk4)
!         STOP
!        WRITE(6,'(5(1pg13.5))') ppadd(1:ngrid)
END IF

deallocate(dpotauxp,gpotaddp,dpotauxn,gpotaddn)
DEALLOCATE(psin,psip)

RETURN
END SUBROUTINE fy_pairpot

END MODULE fyfunc
