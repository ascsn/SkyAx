! Wrapper routine to prepare input for the old F77 code from
! input read in the format of the F90 code.
PROGRAM wrapper_ZPE

USE paramr
USE constr
USE dynpar
!USE fields
USE observ
!USE wavefr
USE levels
USE grid
USE coulomb
USE staticstep
!USE init
USE fyfunc
!USE util, ONLY :plot_2d,write_skyax_header

IMPLICIT NONE


!REAL(8) :: devh2s,dnewoct,econs,oldoct,&
!           qact0,qact2,qact3,qact4,qconsav,&
!           devh2prev=100,econprev=100,x0invmin,pcutact
!REAL(8) :: pairfacact                        ! actual pairing enhancement factor
!REAL(8) :: pairfcsav(2),fxisav,h1xisav,hgradxisav
!INTEGER :: iplold,ipold,iter,itesto,&
!           itmaxact,iwu,nchg,nmass,number_threads=1
!INTEGER :: iqact0,iqact2,iqact3,iqact4,itemp,ntempmax
!INTEGER :: iqact0o,iqact2o,iqact3o,iqact4o    ! save previous values
!INTEGER :: ipairsav            ! saves initial setting of pairing model

INTEGER :: ispins=3    ! level of spin terms: 0 = degen, 1 = no spin terms, 2 = min Gal., 3 = full"
INTEGER :: ioptsk=10   ! way of Sykrma-force input
INTEGER :: iplot=99999 ! plot modulus -- obsolete
INTEGER :: ifold=0 
INTEGER :: ifrso=-1    ! new charge radius
!LOGICAL :: tvaryx0=.FALSE.     ! dynamical step -- not active in old code

! variables for constraint
!REAL(8) :: swicon0=0D0,betmin0,betmax0,betaw0,&
!           swicon2=0D0,betmin2,betmax2,betaw2,&
!           swicon3=0D0,betmin3,betmax3,betaw3,&
!           swicon4=0D0,betmin4,betmax4,betaw4,&
!           stslop=15D0,rhocut=0.03D0,stwid=10D0
!REAL(8) :: qerr=1E-4   
!REAL(8) :: c0=0.8      
!REAL(8) :: d0=0D0      
!REAL(8) :: qeps=-1D3   
!REAL(8) :: dalp=3.5D0  
!REAL(8) :: dgam=0.5D0   
!REAL(8) :: swifac=0.05D0

!INTEGER :: itcons0=99999,itcons2=99999,itcons3=99999,itcons4=99999,itypco
!LOGICAL :: tstab=.FALSE.,tlin=.FALSE.,tfresh=.FALSE.
INTEGER :: ifstab=0,iflin=0


NAMELIST /main/ iprint,itmax,itest,&
                iterbcs,iterhfb,pairenhance,modnat,iternat,tnobcs,&
                v2limh,delmix,&
                title,outnam,&
                ipri_fields,tpri_occup,tpri_conver,tpri_bubble,tpri_2d,&
                tpri_radii,tpri_formf
NAMELIST /grids/ tsym,nr,nz,dr,dz,idrtyp
NAMELIST /static/ nprot,nneut,betain,octini,deocc,mmax,ispins,&
                  dipin,e0inv,x0inv,addrho,addpot,addocc,&
                  errenr,errvar,erra,npiter,itpini,&
                  ndrun,nqrun,norun,nhrun,&
                  colerr,irepc,cmssqr,itpfix,inodiag,tvaryx0
NAMELIST /constraint/ swicon0,betmin0,betmax0,betaw0,&
                      swicon2,betmin2,betmax2,betaw2,&
                      swicon3,betmin3,betmax3,betaw3,&
                      swicon4,betmin4,betmax4,betaw4,&
                      itcons0,itcons2,itcons3,itcons4,&
                      itypco,qerr,c0,d0,qeps,tfresh,&
                      tstab,stwid,stslop,rhocut,&
                      dalp,dgam,tlin,swifac
IF(swicon0.NE.0D0 .OR. swicon3.NE.0D0 .OR. swicon4.NE.0D0) &
  STOP 'only quadrupole constraint implemented presently'

READ(5,main)

CALL read_force
IF(b4p == zero) b4p = t4*oh
IF(ipair<7) rho0pr=0.16D0
IF(tFayans) THEN
  ioptsk=13
  t0=0D0
  t1=0D0
  t2=0D0
  t3=0D0
  t4=0D0
  x0=0D0
  x1=0D0
  x2=0D0
  x3=0D0
  b4p=0D0
  pairfc=(0D0,0D0)
  rho0pr=0.16
ELSE
  isovectormodel = 2
END IF

READ(5,grids)
IF(tsym) THEN
  sym=one
  isym=1
ELSE
  sym=zero
  isym=0
END IF
IF(tcutoff_WS)  pair_cutoff=-pair_cutoff
IF(tcoulfold) ifold=1
IF(.NOT.tcutoff_uvspace) STOP 'SkyAx_odd not compiled for v2 pairing space'
IF(HFBwidth.NE.0D0) STOP 'wrapper not yte designed for finite pairing width'
IF(isovectormodel.NE.2) STOP 'SkyAx_odd compiled only for isovectormodel=2'

READ(5,static) 
IF(MOD(nprot,2)==0 .AND. MOD(nneut,2)==0) ispins=0

READ(5,constraint)

OPEN(31,file='skyinput_old.'//outnam)
WRITE(31,'(a)') " 100000          ITMXT   max. total iterations"
WRITE(31,'(a)') "'"//outnam//": odd code '"
WRITE(31,'(2i4,a)')  nprot,nneut,"  nprot,nneut"
WRITE(31,'(4(1pg20.12),a)')  t0,t1,t2,t3," t0,t1,t2,t3"
WRITE(31,'(4(1pg20.12),a)')  x0,x1,x2,x3,"  x0,x1,x2,x3"
WRITE(31,'(5(1pg20.12),a)')  t4/2,b4p,power,h2m(1),h2m(2),",,, b4,b4p,alpha,amfc,dglan"
WRITE(31,'(i3,f10.7,i3,2f5.1,3i3,a)') ifex,riftm,ifzpe,chg,0D0,ioptsk,ifold,ifrso," if:ex,ls,cm,CHG,optsk,fold"
!IF(ioptsk.eq.11) WRITE(31,'(3(1pg20.12),a/3(1pg20.12),a)') &
!           alpinv,betinv,riftm1,'    alpinv,betinv,riftm1',&
!           b7,b7p,pow7, '    b7,b7p,pow7'
IF(tFayans) THEN
  WRITE(31,'(3(1pg20.12),a)') avp,h1vp,h2vp," avp,h1vp,h2vp"
  WRITE(31,'(3(1pg20.12),a)') avm,h1vm,h2vm," avm,h1vm,h2vm"
  WRITE(31,'(2(1pg20.12),a)') sigmavol,rho0fy," sigmavol,rho0fy"
  WRITE(31,'(3(1pg20.12),a)') asp,hsp,hgradp ," asp,hsp,hgradp "
  WRITE(31,'(3(1pg20.12),a)') asm,hsm,hgradm ," asm,hsm,hgradm "
  WRITE(31,'(2(1pg20.12),a)') hcoul,sigmacoul," hcoul,sigmacoul"
  WRITE(31,'(5(1pg20.12),a)') fxi,h1xi,hgradxi,gammaxi,rho0xi," fxi,h1xi,hgradxi,gammaxi,rho0xi"
  WRITE(31,'(3(1pg20.12),a)') delfxi,-ABS(pairstab),pair_cutoff," delfxi,ecstab,cutfac"
  WRITE(31,'(2(1pg20.12),a)') c0tau,c1tau ," c0tau,c1tau "
  WRITE(31,'(4(1pg20.12),a)') c0nabj,c1nabj,c0tens,c1tens," c0nabj,c1nabj,c0tens,c1tens"
END IF
WRITE(31,'(i4,5(1pg17.9),a)')  ipair,pairfc,-ABS(pairstab),rho0pr,pair_cutoff,"  ipair,pairfp,pairfn"
WRITE(31,'(f4.1,2i4,2f6.2,a)') sym,nr,nz,dr,dz,"  0    0  0 0 0 0  "
WRITE(31,'(4f7.3,3i4,a)')  betain,octini,deocc,mmax,ispins,",,, init: BETAIN, OCTINI,  DEOCC(2),MMAX(2),ispins"
IF(swicon2>0D0) THEN
  WRITE(31,'(4f7.3,a)') swicon2,betmin2,betmax2,betaw2,&
  " 2  0 1,,, CHF: swicon,betmin,betmax,betaw,lqcon,ifresh,icrank"
ELSE
  WRITE(31,'(a)') "  0.0 0.0 0.0  0.025   2  0 0,,, CHF: swicon,betmin,betmax,betaw,lqcon,ifresh"
END IF
WRITE(31,'(a)')   "  0.0  0.0  1.0,,,,,   thermal: tempmn,tempmx,tempw"
WRITE(31,'(4i6,a)') iprint,iplot,itmax,itest,",,  print opt.: IPRINT,IPLOT,ITMAX,ITEST"
WRITE(31,'(5f8.3,a)') E0INV,X0INV,ADDrho,addpot,ADDOCC,"  step iter:E0INV,X0INV,ADDrho,adpot,ADDOCC "
WRITE(31,'(2(1pg13.5),a)') ERRENR,ERRVAR," limit iter:ERRENR,ERRVAR,"
WRITE(31,'(a)')   " 9999, 1,1.5,  extra step: IVS,ITSUB,X0FAC "
WRITE(31,'(a,2i5,a)') "1D-9,1D-9,3D3,3D3,",npiter,itpini,", pair iter:ERRF,ERRA,FERMIP,FERMIN,NPITER,ITPINI"
WRITE(31,'(a)') "0.5,1.0, -2.1,      pair iter:ADDPAR,EPSPAR,WFAC"
WRITE(31,'(a)') "1,1,1,1,   Coul.BC:  ndrun,NQRUN,norun,NHRUN"
WRITE(31,'(a)') "1D-9,100,6D-1,   Coul.Iter:COLERR,IREPC,CMSsQR"
WRITE(31,'(a)') "1D-7 ,-50.0,  20,  1.0, Yukawa Iter:YUKERR,SNL,IREPM,DMSFAC"
IF(swicon2>0D0) THEN
  WRITE(31,'(4(1pg13.5),i8,a)') qerr,c0,d0,qeps,itcons2,&
    "  1,  CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc"
  IF(tstab) ifstab=1
  IF(tlin) iflin=1
  WRITE(31,'(i3,3f7.3,a)') ifstab,stwid,stslop,rhocut," stab.CHF:ifstab,stwid,stslop,rhocut"
  WRITE(31,'(2f7.3,i3,f7.3,a)') dalp,dgam,iflin,swifac," 3.5,0.30,    0,   0.05,   damp.CHF:DALP,DGAM,IFLIN,SWIFAC,"
ELSE
  WRITE(31,'(a)') " 0.1E-3, 0.8, 0.0,  0.3,  999, 1,  CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc"
  WRITE(31,'(a)') " 0,  10.0,  15.0,  0.03,   stab.CHF:ifstab,stwid,stslop,rhocut"
  WRITE(31,'(a)') " 3.5,0.30,    0,   0.05,   damp.CHF:DALP,DGAM,IFLIN,SWIFAC,"
END IF
WRITE(31,'(a)') " 0,  IPARAM = switch type force  "
WRITE(31,'(a)') " 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ireg(2),kblock(2),mblock(2),isoblo,istblo(2)" 
WRITE(31,'(a)') " -1,,,,,,,,,,,,,,, occup.selection:  LOCC(15)  (inactive)  "
WRITE(31,'(a)') "-1, to terminate"
CLOSE(31)

IF(ipair==6 .OR. ipair==8) THEN
  OPEN(31,file='hfbparams')
  WRITE(31,'(3i8,f5.2,i8,2f5.2,l3)') &
     iterhfb,iternat,modnat,pairenhance,iterbcs,v2limh,delmix,tnobcs
  CLOSE(31)
END IF

! make sure that the acual SkyAx_old executable is linked to 'skyax.exe'
CALL system('./skyax.exe < skyinput_old.'//outnam//' >/dev/null ')
!CALL system('../../skyax_mkl.exe < skyinput_old.'//outnam//' >/dev/null ')
CALL system('mv for007.'//outnam//' out_details.'//outnam)
!CALL system('rm colp*')
CALL system('rm quad*')
CALL system('rm conv*')
CALL system('rm field*')

! postludium: combining entries from 'crankn.<outnam> to ZPE of topological GOA

OPEN(31,file='crankn.'//outnam,action='READ')
OPEN(32,file='zpetop.'//outnam,action='WRITE')

CALL zpe_topGOA()

CLOSE(31)
CLOSE(32)


STOP 

END PROGRAM wrapper_ZPE
!-----zpe_topGOA------------------------------------------------------------
 
SUBROUTINE zpe_topGOA()

IMPLICIT REAL(8) (a-h,o-z)

PARAMETER (maxdat=2000)

!     computes interpolated ZPE for beta_2 & rotation & triaxial vibr.
!     from file on unit 31 and writes results to file on uni 32

CHARACTER (LEN=80) :: a
REAL(8) :: beta(maxdat),xmsc(maxdat),energ(maxdat),  &
    xmcra(maxdat),amass(maxdat),widrot(maxdat),  &
    widbet(maxdat),cmcorr(maxdat),xmod(maxdat), zpetri(maxdat)


REAL(8) ::  swfac    !  switch factor of topological GOA
REAL(8) ::  zpesc    !  vibrational ZPE from self-consistent mass
REAL(8) ::  zpecra   !  vibrational ZPE from cranking mass
REAL(8) ::  zrot     !  rotational ZPE (axial rotor)
REAL(8) ::  zpe5     !  total kinetic ZPE triaxial 5D interpolation
REAL(8) ::  zpe3     !  total kinetic ZPE axial 3D grid
REAL(8) ::  zpepot   !  potential ZPE

LOGICAL :: tprintasgiven=.FALSE.  ! testprint of ZPE directly from 'crankn.*'

!-----------------------------------------------------------------------

READ(31,'(a)')  a
WRITE(32,'(a)') '#'//a(1:72)
READ(31,'(a)')  a
WRITE(32,'(a)') '#'//a(1:72)
READ(31,'(1x)')
!      read(31,'(1x)')

IF(tprintasgiven) THEN
  WRITE(32,'(a)') "# ZPE on beta-points from 'crankn.*' input"
  WRITE(32,'(a,a)') '#  beta_2   V_raw  ZPE_3    ZPE_5  switch fac. ',  &
                   '   zpe_vib   zpe_rot   zpe_tri'
END IF

DO i=1,maxdat-1
!old            read(31,*,err=98,end=99) b2,zb2s,zb2,xx,zrot,widrot,energ
  READ(31,*,ERR=98,END=98)  &
      beta(i),xmsc(i),xmcra(i),xmod(i),amass(i),widrot(i),  &
      widbet(i),energ(i),cmcorr(i),xd1,xd2,zpetri(i)
  maxact = i
  
  swfac  = swmon(widrot(i)*0.25,2)
  zpesc  = widbet(i)/(2.0*xmsc(i))
  zpecra = widbet(i)/(2.0*xmcra(i))
  zrot   = widrot(i)/(2.0*amass(i))  ! *2D0
  zpe5   = zpesc*5.0+(zrot+zpetri(i)-4.0*zpesc)*swfac
  zpe3   = zpesc*3.0+(zrot-2.0*zpesc)*swfac
  IF(tprintasgiven) WRITE(32,'(10(1pg13.5))') &
     beta(i),energ(i),zpe3,zpe5,swfac,zpesc,zrot,zpetri(i)
  
END DO
WRITE(6,'(a)') ' max.nr. of input lines exhausted, enhance MAXDAT!'

98   CONTINUE

IF(tprintasgiven) WRITE(32,'(//a)') '# ZPE on interpolated grid:'
WRITE(32,'(a)') '# V_raw   = raw potential energy surface from mean-field calculation', &
               '# ZPE_3   = total kinetic ZPE on axial 3D grid',&
               '# ZPE_5   = total kinetic ZPE on triaxial 5D grid',&
               '# ZPE_pot = potential contribution to ZPE',&
               '# switch fac. = switch factor for topological GOA',&
               '# zpe_vib,zpe_rot,zpi_tri = contributions to top. GOA'
WRITE(32,'(a,a)') '#  beta_2       V_raw        ZPE_3        ZPE_5        ZPE_pot ',  &
    '   switch fac.   zpe_vib       zpe_rot      zpe_tri'
!WRITE(32,'(//a/a,a)') ' corrected ZPE:',  &
!    '  beta_2    switch-fac.     zpe_SC   E   corr.E_SC    ',  &
!    ' zpe    corr.energy  corr.E_SC_pot    corr.E_pot  '

DO i=2,maxact
  betact = 0.5*(beta(i)+beta(i-1))
  CALL polint(beta(i-1),energ(i-1),3,betact,enact,error)
  CALL polint(beta(i-1),widrot(i-1),3,betact,widrac,error)
  CALL polint(beta(i-1),amass(i-1),3,betact,amasac,error)
  CALL polint(beta(i-1),zpetri(i-1),3,betact,ztriac,error)
  enrd2  = derv2(beta(i-1),energ(i-1))
  IF(i >= 3) THEN
    CALL polint(beta(i-2),energ(i-2),3,betact,enac2,error)
    CALL polint(beta(i-2),widrot(i-2),3,betact,widra2,error)
    CALL polint(beta(i-2),amass(i-2),3,betact,amasa2,error)
    CALL polint(beta(i-2),zpetri(i-2),3,betact,ztria2,error)
    IF(i < maxact) THEN
      enact  = 0.5*(enac2+enact)
      widrac = 0.5*(widrac+widra2)
      amasac = 0.5*(amasa2+amasac)
      ztriac = 0.5*(ztria2+ztriac)
      enrd2  = 0.5*(derv2(beta(i-2),energ(i-2))+enrd2)
    ELSE
      enact  = enac2
      widrac = widrac
      amasac = amasa2
      ztriac = ztria2
      enrd2  = derv2(beta(i-2),energ(i-2))
    END IF
  END IF
  widrac = MAX(0.0,widrac)
  cmcora = 0.5*(cmcorr(i)+cmcorr(i-1))
!vold        swfac  = swrot(widrac*0.25)
!old        swfac  = swgen(widrac*0.25,1)
!        ztriac=0D0
  swfac  = swmon(widrac*0.25D0,2)
  zpesc  = widbet(i)/(2D0*xmsc(i))
  zpecra = widbet(i)/(2D0*xmcra(i))
  zpeod  = widbet(i)/(2D0*xmod(i))
  zrot   = widrac/(2D0*amasac)  ! *2.0
  zpepot = enrd2/(8D0*widbet(i))
!
  zpe3   = zpesc*3.0+(zrot-2.0*zpesc)*swfac      ! 5D ZPE with SC mass
  zpe5   = zpesc*5.0+(zrot+ztriac-4.0*zpesc)*swfac      ! 5D ZPE with SC mass
  zpe    = zpecra*5.0+(zrot+ztriac-4.0*zpecra)*swfac    ! 5D ZPE with cranking mass
  encors = enact-zpe5!-cmcora
  encor  = enact-zpe!-cmcora
  WRITE(32,'(10(1pg13.5))') betact,enact,zpe3,zpe5,zpepot, &
                           swfac,zpesc,zrot,zpetri(i)
!  WRITE(32,'(1x,f7.3,f6.3,3(1pg12.4,1pg14.6,1pg12.4),1pg12.4, 4(1pg14.6))')  &
!      betact,swfac,zpe5,enact,encors,xmsc(i),  &
!      zpe,encor,xmcra(i),enrd2/xmcra(i),encors-zpepot,encor-zpepot ,zpesc,zrot
!     &       zpe,encor,xmcra(i),zpepot,encors-zpepot,encor-zpepot
END DO
99   STOP
END SUBROUTINE zpe_topGOA

!-----swrot ------------------------------------------------------------

FUNCTION swrot(alpha)

IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)             :: alpha

!     interpolation formula for rotational ZPE
!     see Z.Phys. A 285 (1978) 93 .

!     Formula is

!       ZPE      = 3*E_beta + (E_rot-E_beta) * swrot(<J_x**2>/2)

!                  int_0^1 a(x^2-1)exp[a(x^2-1)]
!       swrot(a) = -----------------------------
!                      int_0^1 exp[a(x^2-1)]

!     integration is done by brute force trapezoidal rule

!-----------------------------------------------------------------------

IF(alpha == 0.0D0) THEN
  swrot  = 0.0
  RETURN
ELSE IF(alpha < 16.0D0) THEN
  stepx  = 0.004/SQRT(alpha)
  ixmax  = 1.0/stepx
  stepx  = 1.0/ixmax
ELSE
  stepx  = 0.004/SQRT(alpha)
  ixmax  = 4.0/0.004
END IF

sumc   = 0.0
sumd   = 0.5
DO ix=1,ixmax-1
  x      = 1.0-ix*stepx
  arg    = alpha*(x*x-1.0)
  exa    = EXP(arg)
  sumc   = arg*exa    + sumc
  sumd   = exa        + sumd
END DO
arg    = -alpha
exa    = EXP(arg)
sumd    = (0.5*exa+sumd)*stepx
sumc    = (0.5*arg*exa+sumc)*stepx

swrot   = -sumc/sumd

RETURN
END FUNCTION swrot
!-----swgen ------------------------------------------------------------

FUNCTION swgen(alpha,norder)

IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)            :: alpha
INTEGER, INTENT(IN)            :: norder

!     interpolation formula (GCM switch factor) generalised to
!     various dimensions.
!     'norder=0' yields the 1D case as e.g. in particle-number projection
!     'norder=1' yields the 3D case as in the rotational projection
!                                      (axially symmetric)
!     'norder=2' yields the 5D case appropriate for interpolation in the
!                                   beta-gamma plane

!     The switch factor is defined as

!                   2*alpha*I_(norder+2)
!         swgen  =  --------------------
!                    (norder+1)*I_norder

!      I_norder  =  int_0^pi dphi exp(-a*sin(phi/2)**2) sin(phi/2)**norder

!     The limiting cases are

!         swgen --> 1                    for alpha --> infty

!         swgen --> alpha*2/(norder+2)   for alpha --> 0

!-----------------------------------------------------------------------

ixmax  = 200
stepx  = 3.1415926/ixmax
sum2   = 0.0
sum0   = 0.0
DO ix=0,ixmax
  x      = ix*stepx
  sin2   = SIN(x*0.5)**2
  fac    = EXP(-alpha*sin2)*SIN(x*0.5)**norder
  IF(ix == 0) fac = 0.5*fac
  sum0   = fac + sum0
  sum2   = fac*sin2 + sum2
!        write(6,'(1x,f8.2,4(1pg12.4))')
!     &   x,fac,sin2,sum0,sum2
END DO

swgen = (2*alpha/(norder+1))*sum2/sum0
RETURN
END FUNCTION swgen
!-----swmon ------------------------------------------------------------

FUNCTION swmon(alpha,norder)

IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)             :: alpha
INTEGER, INTENT(IN)             :: norder

!     Approximate GCM switch factor which avoids overshooting.
!     It interpolates the slope ate alpha --> 0 with the asymptotis
!     at alpha --> infty.

!-----------------------------------------------------------------------

slope  = 2.0D0/(norder+2.0D0)
swmon  = 1.0D0-EXP(-alpha*slope)

RETURN
END FUNCTION swmon

!-----polint------------------------------------------------------------

SUBROUTINE polint(xa,ya,n,x,y,dy)

IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)             :: xa(n)
REAL(8), INTENT(IN)             :: ya(n)
INTEGER, INTENT(IN)                      :: n
REAL(8), INTENT(IN)             :: x
REAL(8), INTENT(OUT)            :: y
REAL(8), INTENT(OUT)            :: dy

!     polint() polynomial interpolation thru the n points (xa,ya) -> pol(x),
!     returns y=pol(x) and est. error dy
!     (basically taken from "Numerical Recipes")
!     Input:
!     array xa = support points on x axis
!     array ya = function vaues at support points
!            n = order of the interpolating polynomial
!            x = x value at which the function y is wanted

!     Output:
!             y = resulting function value
!            dy = estimated error of interpolation

INTEGER, PARAMETER :: nmax=10
INTEGER :: i,m,ns
REAL(8) :: den,dif,dift,ho,hp,w,c(nmax),d(nmax)

ns=1
dif=ABS(x-xa(1))

DO i=1,n
  dift=ABS(x-xa(i))
  IF (dift < dif) THEN
    ns=i
    dif=dift
  END IF
  c(i)=ya(i)
  d(i)=ya(i)
END DO

y=ya(ns)
ns=ns-1

DO m=1,n-1
  DO i=1,n-m
    ho=xa(i)-x
    hp=xa(i+m)-x
    w=c(i+1)-d(i)
    den=ho-hp
    IF(den == 0.)STOP 'failure in polint'
    den=w/den
    d(i)=hp*den
    c(i)=ho*den
  END DO
  
  IF (2*ns < n-m)THEN
    dy=c(ns+1)
  ELSE
    dy=d(ns)
    ns=ns-1
  END IF
  y=y+dy
END DO

RETURN
END SUBROUTINE polint

!-----derv2 ------------------------------------------------------------

FUNCTION derv2(xin,yin)

IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)             :: xin(3)
REAL(8), INTENT(IN)             :: yin(3)

!     second derivate according to three-point finite difference
!      xin  = vector(1..3) of x-values
!      yin  = vector(1..3) of y-values
!     result is returned as 'derv2'



!-----------------------------------------------------------------------

derv2  = ( (yin(3)-yin(2))/(xin(3)-xin(2))  &
        - (yin(2)-yin(1))/(xin(2)-xin(1)) )/ ((xin(3)-xin(1))*0.5)
RETURN
END FUNCTION derv2

