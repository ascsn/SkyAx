MODULE fyfunc

IMPLICIT NONE

USE paramr, ONLY pi,e2

LOGICAL :: tpairFy
REAL(8) ::  avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,rho0fy, &
            asp,hsp,hgradp,asm,hsm,hgradm,hcoul,sigmacoul, &
            fxi,h1xi,hgradxi,gammaxi,rho0xi, &
            fxip,h1xip,hgradxip,fxin,h1xin,hgradxin, &
            c0tau,c1tau,c0nabj,c1nabj,c0tens,c1tens, &
            r0fy,r0fy2,rho0inv,efermfy,efacfy,ufacfy,pfacfy, &
            clsp,clsm,ctensp,ctensm,ufacfypr,pgradfac 
!INTEGER,PARAMETER :: isovectormodel=3
!INTEGER,PARAMETER :: isovectormodel=2
INTEGER,PARAMETER :: isovectormodel=1


CONTAINS

!-----ini_fyfunc------------ part of the hartree-fock package -----------------

SUBROUTINE ini_fyfunc()

!      IMPLICIT NONE
!INCLUDE 'paramr.inc'

!     Initialize auxiliary parameters for Fayans functional.

!     Parameters are input per common /dynpar/,

!INCLUDE 'dynpar.inc'
!      include 'constn.inc'


!------------------------------------------------------------------------------

!     switch to Fayans model

!      IF(.NOT.ioptsk.EQ.13) STOP "FYFUNC needs IOPTSK=13"

!     auxiliary variables

!          general
r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
r0fy2 = r0fy*r0fy
efermfy = (1.125D0*pi)**(2D0/3D0)*sum(h2m)*0.5D0/r0fy2
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
tpair = .true.

WRITE(7,'(a,8(1pg16.7))') 'rho0fy,efacfy,efermfy,pfacfy=',  &
    rho0fy,efacfy,efermfy,pfacfy,ctensp,ctensm,clsp,clsm
CALL flush(7)

RETURN
END SUBROUTINE ini_fyfunc

!-----fyfunc------------ part of the hartree-fock package ---------------------

SUBROUTINE fyfunc(rp,rn,rc,d1pr,d1pz,d1nr,d1nz,d2p,d2n,tp,tn,  &
    tlp,tln,sp,sn,ssp,ssn, potp,potn,bmasp,bmasn,bmaslp,bmasln,  &
    wlspr,wlspz,wlsnr,wlsnz,energd)

!      IMPLICIT NONE
INCLUDE 'paramr.inc'

!     the basic Fayans energy density and derived quantities.

!     Input per list:
!      rp,rn        = densities
!      rc           = charge density (used in Coulomb exchange)
!      d1p,d1n      = derivative of densities
!      d2p,d2n      = Laplacian of densities
!      tp,tn        = kinetic energy densities
!      tlp,tln      = L**2 densities
!      sp,sn        = spin orbit densities
!      ssp,ssn      = derivative of spin-orbit densities
!      frp,frn      = folded densities (for Yukawa term)

!     Forces are input per common /skypar/ und /tswitc/.
!     Ouput per list:
!      energd        = energy density
!      potp,potn     = local potentials
!      bmasp,bmasn   = inverse mass
!      bmaslp,bmasln = inverse L**2 mass
!      wlsp,wlsn     = spin-orbit potential

INCLUDE 'dynpar.inc'
!      include 'numpar.inc'
INCLUDE 'wavefr.inc'
!      include 'constn.inc'


REAL(:: 8),intent(in), INTENT(IN)        :: rp
REAL(:: 8),intent(in), INTENT(IN)        :: rn
REAL(:: 8),intent(in), INTENT(IN)        :: rc
REAL(:: 8),intent(in), INTENT(IN)        :: d1pr
REAL(:: 8),intent(in), INTENT(IN)        :: d1pz
REAL(:: 8),intent(in), INTENT(IN)        :: d1nr
REAL(:: 8),intent(in), INTENT(IN)        :: d1nz
REAL(:: 8),intent(in), INTENT(IN)        :: d2p
REAL(:: 8),intent(in), INTENT(IN)        :: d2n
REAL(:: 8),intent(in), INTENT(IN)        :: tp
REAL(:: 8),intent(in), INTENT(IN)        :: tn
REAL(:: 8),intent(in), INTENT(IN OUT)    :: tlp
REAL(:: 8),intent(in), INTENT(IN OUT)    :: tln
REAL(:: 8),intent(in), INTENT(IN)        :: sp
REAL(:: 8),intent(in), INTENT(IN)        :: sn
REAL(:: 8),intent(in), INTENT(IN)        :: ssp
REAL(:: 8),intent(in), INTENT(IN)        :: ssn
REAL(:: 8),intent(out), INTENT(OUT)      :: potp
REAL(:: 8),intent(out), INTENT(OUT)      :: potn
REAL(:: 8),intent(out), INTENT(OUT)      :: bmasp
REAL(:: 8),intent(out), INTENT(OUT)      :: bmasn
REAL(:: 8),intent(out), INTENT(IN OUT)   :: bmaslp
REAL(:: 8),intent(out), INTENT(IN OUT)   :: bmasln
REAL(:: 8),intent(out), INTENT(OUT)      :: wlspr
REAL(:: 8),intent(out), INTENT(OUT)      :: wlspz
REAL(:: 8),intent(out), INTENT(OUT)      :: wlsnr
REAL(:: 8),intent(out), INTENT(OUT)      :: wlsnz
REAL(:: 8),intent(out), INTENT(OUT)      :: energd



REAL(:: 8),PARAMETER :: dsmall=1D-30
LOGICAL,:: SAVE :: tfirst=.true.

REAL(:: 8) :: rs,ts,sss,xplus,xminus,xpsigma,xmsigma,d1sr2,d1mr2

!------------------------------------------------------------------------------

!     switch to Fayans model

IF(.NOT.ioptsk == 13) STOP "FYFUNC needs IOPTSK=13"

!     auxiliary variables

!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      efacfy = efermfy*rho0fy*(2D0/3D0)
!      ufacfy = efermfy/3D0
!      rho0inv = 0.5D0/rho0fy
!      clsp = 2D0*efermfy/(3D0*rho0fy)*r0fy2*c0nabj
!      clsm = 2D0*efermfy/(3D0*rho0fy)*r0fy2*c1nabj

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

!      IF(rc==zero) THEN
!        WRITE(911,'(a,8(1pg13.5))') 'FYFUNC:',
!     &    rs,1D0/rho0inv,rp-rn,avm,h1vm,h2vm,efacfy,
!     &    efacfy*avm*(xminus*xminus)*(1D0-h1vm*xplus)/
!     &    (1D0+h2vm*xplus)/rs
!      END IF
energd = efacfy*(avp*(xplus*xplus)* (1D0-h1vp*xpsigma)/(1D0+h2vp*xpsigma)  &
    +avm*(xminus*xminus)* (1D0-h1vm*xplus)/(1D0+h2vm*xplus)  &
    +asp*d1sr2/(1D0+hsp*xpsigma+hgradp*d1sr2)  &
    +asm*d1mr2/(1D0+hsm*xplus+hgradm*d1mr2) )  &
    +c0tau*rs*ts+c1tau*(rp-rn)*(tp-tn)
!      energ2=energd-energ1
!     &         +ctensp*ss*ss+ctensm*(sp-sn)**2
!     &         +clsp*rs*sss+clsm*(rp-rn)*(ssp-ssn)
!old     &                 (1D0-h1vm*xmsigma)/(1D0+h2vm*xmsigma)

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
  add0 = -excfac*rc**0.333333333333D0
  add = add0*(1D0-hxpsigmac)
  energd = 0.75D0*add*rc + energd
  potp = add - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus + potp
  potn =  - 0.75D0*sigmacoul*add0*hxpsigmac*rc*rho0inv/xplus + potn
END IF
!      energ3=energd-energ2

!     spin-orbit

!old      potp   = -clsp*sss - cslm*(ssp-ssn)   + potp
potp   = -clsp*sss - clsm*(ssp-ssn)   + potp
potn   = -clsp*sss + clsm*(ssp-ssn)   + potn
wlspr   = (clsp*d1sr+clsm*(d1pr-d1nr)) +ctensp*ss+ctensm*(sp-sn)
wlsnr   = (clsp*d1sr-clsm*(d1pr-d1nr)) +ctensp*ss-ctensm*(sp-sn)
wlspz   = (clsp*d1sz+clsm*(d1pz-d1nz)) +ctensp*ss+ctensm*(sp-sn)
wlsnz   = (clsp*d1sz-clsm*(d1pz-d1nz)) +ctensp*ss-ctensm*(sp-sn)
energd = -clsp*rs*sss - clsm*(rp-rn)*(ssp-ssn)  &
    +ctensp*ss*ss+ctensm*(sp-sn)**2 + energd

!      energ4=energd-energ3

!
IF(tfirst) THEN
  WRITE(6,'(a,6(1pg13.5))') 'r0fy,efermfy,C_0,ufacfy,rhoinv=',  &
      r0fy,efermfy,2D0*efermfy/(3D0*rho0fy),ufacfy,rho0inv,asp
  WRITE(6,'(a,4(1pg13.5))') 'potp etc=', rp,rn,xpsigma,potp,potn,energd
  WRITE(6,*) 'clsp,c0nabj,clsm,c1nabj=',clsp,c0nabj,clsm,c1nabj
  tfirst=.false.
END IF
!      WRITE(911,'(a,8(1pg16.7))')'rs,xplus,rho0inv,E=',rs,xplus,rho0inv,
!     &   energ1/rs,energ2/rs,energ3/rs,energ4/rs

RETURN
END SUBROUTINE fyfunc
!-----fyfunc_odd------------ part of the hartree-fock package ---------------------

SUBROUTINE fyfunc_odd(rotjrp,rotjrn,rotjzp,rotjzn,  &
    spdenrp,spdenrn,spdenzp,spdenzn,energd)

!      IMPLICIT NONE
INCLUDE 'paramr.inc'
INCLUDE 'dynpar.inc'
!     Energy contribution from minimal Galilean invariant  odd-odd terms.
!     ! presently only for l*s
!
!        rotjrp/n  =    r-component of rotation of current, protons and neutrons
!        rotjzp/n  =    z-component of rotation of current, protons and neutrons
!        spdenrp/n  =   r-component of spin density, protons and neutrons
!        spdenzp/n  =   r-component of spin density, protons and neutrons
!
rotjrs = rotjrp+rotjrn
rotjzs = rotjzp+rotjzn
rotjrd = rotjrp-rotjrn
rotjzd = rotjzp-rotjzn
spdenrs = spdenrp+spdenrn
spdenzs = spdenzp+spdenzn
spdenrd = spdenrp-spdenrn
spdenzd = spdenzp-spdenzn
energd = energd - clsp*(spdenrs*rotjrs+spdenzs*rotjzs)  &
    - clsm*(spdenrd*rotjrd+spdenzd*rotjzd)
RETURN
END SUBROUTINE fyfunc_odd

!-----fyfunc_grad------------ part of the hartree-fock package -----------------

SUBROUTINE fyfunc_grad()

!      IMPLICIT NONE
INCLUDE 'paramr.inc'

!     The derivative contribution to the gradient term of the Fayans functional-

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

INCLUDE 'dynpar.inc'
INCLUDE 'fields.inc'
INCLUDE 'wavefr.inc'
!      include 'grid.inc'
!      include 'constn.inc'

REAL(:: 8), allocatable :: dpotaux(:,:,:),gpotadd(:,:,:,:)

!------------------------------------------------------------------------------

allocate(dpotaux(0:nzd,0:nrd,kqm),gpotadd(0:nzd,0:nrd,kqm,2))
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      efacfy = efermfy*rho0fy*(2D0/3D0)
!      ufacfy = efermfy/3D0
!      rho0inv = 0.5D0/rho0fy
!      WRITE(6,*) 'FYFUNC_GRAD:',r0fy,r0fy2,efermfy,efacfy,
!     &         efacfy*asp*r0fy2

!     gradients of density
CALL nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,one)  ! d1z proton
CALL nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,one)  ! d1z neutron
CALL nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,one)  ! d1r proton
CALL nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,one)  ! d1r neutron

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

CALL nablaz_scalar(dpotaux(:,:,1),gpotadd(:,:,1,1), 1,1,one)
CALL nablaz_scalar(dpotaux(:,:,2),gpotadd(:,:,2,1), 1,1,one)
CALL nablar_scalar(dpotaux(:,:,1),gpotadd(:,:,1,2), 0,1,one)
CALL nablar_scalar(dpotaux(:,:,2),gpotadd(:,:,2,2), 0,1,one)

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

deallocate(dpotaux,gpotadd)
RETURN
END SUBROUTINE fyfunc_grad
!-----fy_hfbpot------------ part of the hartree-fock package -----------------

SUBROUTINE fy_hfbpot()

!      IMPLICIT NONE
INCLUDE 'paramr.inc'

!     The pairing potential of the Fayans functional.
!     Also used to compute pairing energy.

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

INCLUDE 'dynpar.inc'
!      include 'numpar.inc'
INCLUDE 'fields.inc'
INCLUDE 'wavefr.inc'
INCLUDE 'observ.inc'
!      include 'grid.inc'
!      include 'constn.inc'
REAL(:: 8) :: d1sr,d1sz,d1sr2,xplus,xplusgamma,rhfacp,rhfacn
LOGICAL,:: PARAMETER :: ttest=.false.

!------------------------------------------------------------------------------

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
  WRITE(6,*) 'FY_HFBPOT'
END IF
!        r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!        r0fy2 = r0fy*r0fy
!        efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      pfacfy = 2D0*efermfy/(3D0*rho0fy)
!      ufacfypr = pfacfy/(2D0*rho0xi)
!     gradients of density     <-- recycled from FY_PAIRPOT
CALL nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,one)  ! d1z proton
CALL nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,one)  ! d1z neutron
CALL nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,one)  ! d1r proton
CALL nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,one)  ! d1r neutron
!
wrk4=0D0
wrk3=0D0
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
  IF(ttest) THEN
    wrk3(iz,ir) =  xplusgamma
    wrk4(iz,ir) =  d1sr2
  END IF
END DO; END DO
epair(2) = f1int(wrk2)
epair(1) = f1int(wrk1)-epair(2)
!      IF(nprot==20 .AND. nneut==20) THEN
!        WRITE(6,'(a,4(1pg13.5))') 'FYpair:',
!     &      pfacfy*fxip,pfacfy*h1xip,pfacfy*hgradxip,epair
!      END IF
WRITE(6,'(a,5(1pg13.5))') 'FYpair:',epair, f1int(wrk3),f1int(wrk4),pfacfy
IF(ttest) THEN
  WRITE(7,'(a,4(1pg13.5))')'params:',pfacfy,fxin,h1xin,hgradxin
  CALL testpr(paircr(:,:,2),0,9,0,9,'paircr')
  CALL testpr(pairpt(:,:,2),0,9,0,9,'pairpt')
  CALL testpr(wrk3,0,9,0,9,'xplusgamm')
  CALL testpr(wrk4,0,9,0,9,'d1sr2')
  CALL testpr(wrk2,0,9,0,9,'rhfacn')
  STOP
END IF
!        write(6,'(a)') ' Fy pair potential for neutrons:'
!        write(6,'(5(1pg12.4))') (pairpt(i),i=ngrid+1,2*ngrid)

RETURN
END SUBROUTINE fy_hfbpot
!-----fy_pairpot------------ part of the hartree-fock package -----------------

SUBROUTINE fy_pairpot()

!      IMPLICIT NONE
INCLUDE 'paramr.inc'

!     Non-gradient part of the  pairing contribnution
!     to potential from the Fayans functional.

!     Parameters are input per common /skypar/,
!     potentials and densities per /wavfun/.

!     The result is accumulated onto the array 'pot'.

INCLUDE 'dynpar.inc'
INCLUDE 'fields.inc'
INCLUDE 'wavefr.inc'
INCLUDE 'derivs.inc'
!      include 'constn.inc'
REAL(:: 8), allocatable :: dpotauxp(:,:),gpotaddp(:,:,:)
REAL(:: 8), allocatable :: dpotauxn(:,:),gpotaddn(:,:,:)
REAL(:: 8) :: d1sr,d1sz,d1sr2,d2sr,xplus,xplusgamma,add
LOGICAL,:: PARAMETER :: ttest=.false.

!------------------------------------------------------------------------------

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
allocate(dpotauxp(0:nzd,0:nrd),gpotaddp(0:nzd,0:nrd,2))
allocate(dpotauxn(0:nzd,0:nrd),gpotaddn(0:nzd,0:nrd,2))
!      r0fy = (0.375D0/(pi*rho0fy))**(1D0/3D0)
!      r0fy2 = r0fy*r0fy
!      efermfy = (1.125D0*pi)**(2D0/3D0)*(amfc0(1)+amfc0(2))*0.5D0/r0fy2
!      ufacfypr = efermfy/(3D0*rho0fy*rho0xi)
dpotauxp(0:nz,0:nr) = paircr(0:nz,0:nr,1)**2
dpotauxn(0:nz,0:nr) = paircr(0:nz,0:nr,2)**2
!      WRITE(6,*) 'paircr**2:',SUM(dpotaux(1:ngrid))
CALL nablaz_scalar(dpotauxp,gpotaddp(:,:,1), 1,1,one)
CALL nablaz_scalar(dpotauxn,gpotaddn(:,:,1), 1,1,one)
CALL nablar_scalar(dpotauxp,gpotaddp(:,:,2), 0,1,one)
CALL nablar_scalar(dpotauxn,gpotaddn(:,:,2), 0,1,one)
!     gradients of density
CALL nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,one)  ! d1z proton
CALL nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,one)  ! d1z neutron
CALL nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,one)  ! d1r proton
CALL nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,one)  ! d1z neutron
IF(ttest) THEN
  CALL laplace_scalar(rho(0,0,1),d2rho(0,0,1),+1,0,1,one)
  CALL laplace_scalar(rho(0,0,2),d2rho(0,0,2),+1,0,1,one)
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
  CALL testpr(paircr(:,:,2),0,9,0,9,'paircr')
  CALL testpr(pot(:,:,2),0,9,0,9,'add pot')
  CALL testpr(wrk3,0,9,0,9,'xplusgamm')
  CALL testpr(wrk4,0,9,0,9,'d1sr')
  CALL testpr(wrk2,0,9,0,9,'d2sr')
!         STOP
!        WRITE(6,'(5(1pg13.5))') ppadd(1:ngrid)
END IF

deallocate(dpotauxp,gpotaddp,dpotauxn,gpotaddn)
RETURN
END SUBROUTINE fy_pairpot

END MODULE fyfunc
