c include block: dynpar
c
      character*80 title
      LOGICAL :: thfb
      LOGICAL, PARAMETER :: ttestconstr=.TRUE.
      REAL(8) :: v2lim,delmix         ! for damped step in 'lagpar'
      common  /iterpr/  cmono,cquad,chex,cdip,coct,
     &                  cmssqr,colerr,snl,snlhc,yukerr,dmsfac,
     &                  errenr,errvar,swdiag,    
     &                  pairenhance,v2limh,delmix,
     &                  iterhfb,iternat,modnat,iterbcs,
     &                  irepm,irepc,itpfix,
     &                  itert,      ! nr. of total iterations
     &                  thfb,title
*                                                                     *
***********************************************************************
*                                                                     *
*        number of iteration and acceleration terms for the solution  *
*        of meson equations with the yukawa solver and for the        *
*        coulomb fields.                                              *
*                                                                     *
***********************************************************************
c
      common  /gridpr/  dr,dz,sr,sz,drsqr,dzsqr,dvol,sym,
     &                  nr,nz,nrm1,nzm1,nro,nzo,isym,icm,
     &                  icutr,icutz,idrtyp
*                                                                     *
***********************************************************************
*                                                                     *
*        geometrical and mesh parameters:                             *
*          dr     = step size in r-direction                          *
*          dz     = step size in z-direction                          *
*          sr     = increment in r-direction  ( r=(ir+sr)*dr )        *
*          sz     = increment in z-direction  ( z=(iz+sz)*dz )        *
*          drsqr  = dr*dr                                             *
*          dzsqr  = dz*dz                                             *
*          dvol   = volume element for space integration              *
*          sym    = symmetry parameter,                               *
*                    1 -> z-reflection symmetry                       *
*                    0 -> no z-reflection symmetry                    *
*          nr     = grid size in r-direction                          *
*          nz     = grid size in z-direction                          *
*          nrm1   = nr-1                                              *
*          nzm1   = nz-1                                              *
*          isym   = ifix(sym)                                         *
*          icm    = point on z-grid for center of the nucleus         *
*                    icm=1      for sym=0.0                           *
*                    icm=nz/2+1 for sym=1.0                           *
*          icutr  = momentum comp. cutoff in r-direction              *
*          icutz  = momentum comp. cutoff in z-direction              *
*          idrtyp = switch for type of derivative                     *
*                                                                     *
***********************************************************************
c
      common /dampar/ x0inv,e0inv,addrho,addpot,addocc,x0fac,adyrho,
     &                adypot,mpmax,itsub
*
***********************************************************************
*                                                                     *
*        variables defining the damping operator for the gradient     *
*        iteration step:                                              *
*         x0inv  = overall factor in front of the damping operator    *
*         e0inv  = witdh of the damping operator                      *
*         addrho = fraction of new density to be added to the old,    *
*                  used to compute 'adyrho'                           *
*         addpot = fraction of new potential to be added to the old.  *
*         addocc = fraction of new occupation numbers 'VSQF'          *
*                  to be added to the old                             *
*         adyrho = dynamical  fraction of new density                 *
*         adypot = dynamical  fraction of new potential               *
*         mpmax  = max no. of z component of angular momentum         *
*         x0fac  = facor on 'x0inv' in sub-iterations                 *
*         itsub  = number of subiterations                            *
*        note: the damping operator is                                *
*                  x0inv/(1.+kinetic energy/e0inv) .                  *
*                                                                     *
***********************************************************************
c
      REAL(8) :: b_str_act(2)
      common  /forces/  t0,t1,t2,t3,t4,
     &                  x0,x1,x2,x3,b4p,power,h2m(2),h2meff(2),h2mav,
     &                  dglan,dglanp,b_strength,b_str_act,
     &                  chg,chgx,
     &                  iparam,ifzpe,ifls,ioptsk,icoulfold,
     &                  ifrso
      common /wspara/ ws0(2),wswid(2),wsb0(2),wsb2(2),wsb4(2),
     &                wskap(2),wsrad0(2)
*                                                                     *
***********************************************************************
*                                                                     *
*        skyrme -  parameters                                         *
*         t0-t4,x0-x3  --> usual convention                           *
*         power        --> power of density dependence in 't3'        *
*         h2m    = hbar**2/baryonmass                                 *
*         h2meff       --> effective mass (for 1-1/A correction)      *
*         dglan   = additional contribution to Landau parameter g_0   *
*         chg     = effective charge in units ?e2! , usually 1.0      *
*         chgx    = effective charge for exchange term                *
*         iparam  = number of preset force                            *
*         ifzpe   = 0 --> m-m/A recipe for c.m. correction            *
*                   1 --> no c.m. correction                          *
*                                                                     *
*         ifls    = switch for contribution from t1,t2 to             *
*                   spin-orbit force ('iftm' in spherical code)       *
*         icoulfold = activate folding of charge density              *
*         tpairspace = switches to v**2 pairing space
*                                                                     *
*         /wspara/ carries the Woods-Saxon parameters (for iparam=-1) *
*                                                                     *
***********************************************************************
c
      LOGICAL,PARAMETER :: tHFBvariat=.FALSE.   ! variational HFB Hamiltonian
      LOGICAL :: tpair
      LOGICAL :: tpairspace
      common  /pairpr/  ferp(2),ggg(2),gap(2),fer(2),fermx(2),ewin(2),
     x                  epair(2),fer2(2),pairfc(2),feradd(2),pairstab,
     x                  rho0pr,cutfac,addpar,epspar,wfac,errf,erra,
     &                  stwid,stslop,rhocut,ifstab,
     &                  npiter,ipair,ipairin,itpini,
     &                  tpairspace,tpair
*                                                                     *
***********************************************************************
*                                                                     *
*                      pairing parameters.                            *
*            and the parameters 'ifstab,stwid,stslop,rhocut' for the  *
*            stabilizing external parameters.                         *
*                                                                     *
***********************************************************************



c
c     Parameters for the Fayans functional
c
      LOGICAL :: tpairFy
      REAL(8) ::avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,rho0fy,
     &          asp,hsp,hgradp,asm,hsm,hgradm,hcoul,sigmacoul,
     &          fxi,h1xi,hgradxi,gammaxi,rho0xi,
     &          fxip,h1xip,hgradxip,fxin,h1xin,hgradxin,
     &          c0tau,c1tau,c0nabj,c1nabj,c0tens,c1tens,
     &          r0fy,r0fy2,rho0inv,efermfy,efacfy,ufacfy,pfacfy,
     &          clsp,clsm,ctensp,ctensm,ufacfypr,pgradfac
!      INTEGER,PARAMETER :: isovectormodel=3
      INTEGER,PARAMETER :: isovectormodel=2
!      INTEGER,PARAMETER :: isovectormodel=1
      common /faypar/avp,h1vp,h2vp,avm,h1vm,h2vm,sigmavol,rho0fy,
     &               asp,hsp,hgradp,asm,hsm,hgradm,hcoul,sigmacoul,
     &               fxi,h1xi,hgradxi,gammaxi,rho0xi,delfxi,
     &               fxip,h1xip,hgradxip,fxin,h1xin,hgradxin,
     &               c0tau,c1tau,c0nabj,c1nabj,c0tens,c1tens,
     &               r0fy,r0fy2,rho0inv,efermfy,efacfy,ufacfy,pfacfy,
     &               clsp,clsm,ctensp,ctensm,ufacfypr,pgradfac,
     &               tpairFy



c end include block
