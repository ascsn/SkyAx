! Wrapper routine to prepare input for the old F77 code from
! input read in the format of the F90 code.
PROGRAM wrapper

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
LOGICAL :: tvaryx0=.FALSE.     ! dynamical step -- not active in old code


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
  WRITE(31,'(3(1pg20.12),a)') delfxi,ABS(pairstab),pair_cutoff," delfxi,ecstab,cutfac"
  WRITE(31,'(2(1pg20.12),a)') c0tau,c1tau ," c0tau,c1tau "
  WRITE(31,'(4(1pg20.12),a)') c0nabj,c1nabj,c0tens,c1tens," c0nabj,c1nabj,c0tens,c1tens"
END IF
WRITE(31,'(i4,5(1pg17.9),a)')  ipair,pairfc,ABS(pairstab),rho0pr,pair_cutoff,"  ipair,pairfp,pairfn"
WRITE(31,'(f4.1,2i4,2f6.2,a)') sym,nr,nz,dr,dz,"  0    0  0 0 0 0  "
WRITE(31,'(4f7.3,3i4,a)')  betain,octini,deocc,mmax,ispins,",,, init: BETAIN, OCTINI,  DEOCC(2),MMAX(2),ispins"
WRITE(31,'(a)') "  0.0 0.0 0.0  0.025   2  0 0,,, CHF: swicon,betmin,betmax,betaw,lqcon,ifresh"
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
WRITE(31,'(a)') " 0.1E-3, 0.8, 0.0,  0.3,  999, 1,  CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc"
WRITE(31,'(a)') " 0,  10.0,  15.0,  0.03,   stab.CHF:ifstab,stwid,stslop,rhocut"
WRITE(31,'(a)') " 3.5,0.30,    0,   0.05,   damp.CHF:DALP,DGAM,IFLIN,SWIFAC,"
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


CALL system('../../skyax_mkl.exe < skyinput_old.'//outnam//' >/dev/null ')
CALL system('mv for007.'//outnam//' out_details.'//outnam)
CALL system('rm colp*')
CALL system('rm quad*')
CALL system('rm conv*')
CALL system('rm field*')

STOP 

END
