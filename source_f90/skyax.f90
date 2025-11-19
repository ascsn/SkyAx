!-----skyax -----------------------------------------------------------
PROGRAM skyax

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
USE coulomb
USE staticstep
USE init
USE fyfunc
USE util, ONLY :plot_2d,write_skyax_header

IMPLICIT NONE

!CHARACTER (LEN=20) :: set

LOGICAL,PARAMETER :: tafter=.true.           ! switch to "afterburn"

LOGICAL ::           tvaryx0=.FALSE.

REAL(8) :: devh2s,dnewoct,econs,oldoct,&
           qact0,qact2,qact3,qact4,qconsav,&
           devh2prev=100,econprev=100,x0invmin
INTEGER :: iplold,ipold,iter,itesto,&
           itmaxact,iwu,nchg,nmass,number_threads=1
INTEGER :: iqact0,iqact2,iqact3,iqact4,itemp,ntempmax
INTEGER :: iqact0o,iqact2o,iqact3o,iqact4o    ! save previous values

NAMELIST /main/ iprint,itmax,itest,title,outnam,&
                ipri_fields,tpri_occup,tpri_conver,tpri_bubble,tpri_2d,&
                tpri_radii,tpri_formf
NAMELIST /grids/ tsym,nr,nz,dr,dz,idrtyp
NAMELIST /static/ nprot,nneut,betain,octini,deocc,mmax,&
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
NAMELIST /thermal/ tempmn,tempmx,tempw

!**********************************************************************
!                                                                     *
!     read namelist main                                              *
!                                                                     *
!**********************************************************************
  READ(5,main)
  
  IF(tpri_bubble) THEN
    OPEN(173,FILE='bubbleness.'//outnam)
    WRITE(173,'(a)') '# protocol of central density',  &
        '# beta2   energy   rho_p(0,0)  rho_n(0,0)  bubbleness (total,prot,neuut)'
  END IF
  
  OPEN(7,STATUS='unknown',FORM='formatted',FILE='out_details.'//outnam)
  CALL write_skyax_header
  IF(ifoccs == 1) THEN
    WRITE(7,'(a)') ' WARNING: running IFOCCS=1 is a test version !'
    WRITE(6,'(a)') ' WARNING: running IFOCCS=1 is a test version !'
  END IF
!**********************************************************************
!                                                                     *
!     read namelist force                                             *
!                                                                     *
!**********************************************************************
  CALL read_force
  IF(tFayans) THEN
    CALL ini_fyfunc()
    tpair = .TRUE.
  END IF

!  ggp=pairfc(1); ggn=pairfc(2)              ! eliminate ggp,ggn later
  IF(riftm == 0D0) THEN                     ! riftm: Why not logical variable as input?
    ifls = 0
  ELSE
    ifls = 1
  END IF
    riftm1 = 0D0
  IF(tadsic) THEN
    chg   = chg*(1D0-1D0/nchg)
    chgx  = chgx*(1D0-(1D0/nchg)**0.3333333333D0)
  END IF
  IF(riftm1 == 0D0) riftm1 = riftm  !new
  IF(itest>0) WRITE(6,*) 'ipair,v0prot,v0neut,pairstab,rho0pr=',&
    ipair,pairfc,pairstab,rho0pr  
  
!**********************************************************************
!                                                                     *
!     read namelist grid                                              *
!                                                                     *
!**********************************************************************

  READ(5,grids)
  IF(tsym) THEN
    sym=one
    isym=1
  ELSE
    sym=zero
    isym=0
  END IF
  
  IF(itest>0) WRITE(6,*)  'Grids: ',tsym,sym,nr,nz,dr,dz,idrtyp

!**********************************************************************
!                                                                     *
!     read namelist static                                            *
!                                                                     *
!**********************************************************************
  READ(5,static) 
  nmass  = nprot+nneut
  nchg   = nprot
  npart(1) = nchg
  npart(2) = nmass-nchg
  anumb  = npart(1)+npart(2) 
  IF(tvaryx0) x0invmin=x0inv
  IF(tvaryx0) x0inv=x0inv*2.0
!**********************************************************************
!                                                                     *
!     read namelist constraint                                        *
!                                                                     *
!**********************************************************************

  READ(5,constraint)                         

!**********************************************************************
!                                                                     *
!     read namelist thermal                                           *
!                                                                     *
!**********************************************************************      

  READ(5,thermal)


!  kblock = 0                               

  IF(MAX(mmax(1),mmax(2)) > kmpsi+1)  &
      STOP ' MMAX too large. enhance KPSI in PARAMR.INC!'
  IF(itpfix==0) itpfix=itmax+1

                               ! force parameters
  istblo = 0
  IF(b4p == zero) b4p = t4*oh
  
  CALL prep_qsoll(nmass*1D0) !Also sets grid size
  CALL alloc_fields()
  
!     check consistency of input

  IF(tswicon3) THEN
    octini = 0D0
    WRITE(6,'(a)') ' constraint calc. --> initial octupole kick deactivated'
  END IF
  IF(dipin /= 0D0 .AND. tsym)  &
      STOP 'external dipole field requires TSYM false'
  IF(tswicon3 .AND. tsym) STOP ' octupole constraint requires TSYM false'

!**********************************************************************
!                                                                     *
!        write the parameters and options                             *
!                                                                     *
!**********************************************************************

  CALL write_headers()  

!**********************************************************************
!                                                                     *
!     initialisation of integration weights, pairing etc.             *
!                                                                     *
!**********************************************************************

  IF(ifzpe == 0) THEN
    h2meff(1)  = h2m(1)-h2m(1)/nmass                      ! optionally ZPE by m*
    h2meff(2)  = h2m(2)-h2m(2)/nmass                      ! optionally ZPE by m*
    WRITE(7,'(10X,a,2f10.3)') 'ZPE with m-m/A, h2m_eff now=',h2meff
  ELSE
    h2meff(1)  = h2m(1)
    h2meff(2)  = h2m(2)
    WRITE(7,'(10X,a,i3,a)') &
         ' ZPE with <P^2>/2M (evaluated every ',itzpe,' iterations)'
  END IF
  h2mav = 0.5D0*(h2meff(1)+h2meff(2))
  chg   = chg * e2        ! ?? changing input data is dangerous !
  eocc(1) = 1.5*npart(1)**edr + deocc(1)
  eocc(2) = 1.5*npart(2)**edr + deocc(2)
  IF(2*eocc(1) > mmax(1)+3.0)  &
       STOP ' MMAX too small for given DEOCC for protons'
  IF(2*eocc(2) > mmax(2)+3.0)  &
      STOP ' MMAX too small for given DEOCC for neutrons'
  CALL init_grid            ! grid and derivative matrices
  CALL idanuc               ! nucleon damping
  CALL idames               ! coulomb damping

  IF(.NOT.tfresh) THEN
    CALL ssetup                     ! the wavefunctions
    IF(ipair>0) CALL inipai(nmass)
  END IF

!**********************************************************************
!                                                                     *
!     loop over deformations                                          *
!                                                                     *
!**********************************************************************
  
  iqact0o=0;iqact2o=0;iqact3o=0;iqact4o=0;
  itert  = 0                         ! set counter for total iterations
  DO  iqact0=0,nqmax0; qact0=iqact0*qw0+qmin0
  DO  iqact2=0,nqmax2; qact2=iqact2*qw2+qmin2
  DO  iqact3=0,nqmax3; qact3=iqact3*qw3+qmin3
  DO  iqact4=0,nqmax4; qact4=iqact4*qw4+qmin4
    ! re-initialize state and constraints either always (tfresh==.TRUE.)
    ! or if constraint change too much (iqact...)
    IF(tfresh .OR. ABS(iqact0-iqact0o)>1 .OR. ABS(iqact2-iqact2o)>1 &
       .OR. ABS(iqact3-iqact3o)>1 .OR. ABS(iqact4-iqact4o)>1) THEN
      IF(tswicon3) octini = qact3/cfact3
      CALL ssetup                    
      CALL inipai(nmass)
      qlamb0=0D0;qlamb2=0D0;qlamb3=0D0;qlamb4=0D0
      qlambi0=0D0;qlambi2=0D0;qlambi3=0D0;qlambi4=0D0 
    END IF
    iqact0o=iqact0;iqact2o=iqact2;iqact3o=iqact3;iqact4o=iqact4;

    qsoll0=qact0; qsoll2=qact2; qsoll3=qact3; qsoll4=qact4
    IF(tfixcon0) qlamb0=qact0
    IF(tfixcon2) THEN
      qlamb2 = qact2
      qsoll2 = betain * cfact2  ! for 'inicms'
    END IF
    IF(tfixcon3) qlamb3=qact3
    IF(tfixcon4) qlamb4=qact4
    IF(itest >= 1) THEN
      WRITE(73,'(a,a)')  &
          '   qtilc      qcst     varqsc      qq       ratio    qdh/qq  ',  &
          '  qlcstp     ratio'
      WRITE(74,'(a,a)')  &
          '    qtil        q      varqsq      qq       ratio    qdh/qq  ',  &
          '  qlstep     ratio'
    END IF
    WRITE(*,*) 'qsoll:',iqact0,qact0,qsoll0,qlamb0
!     initialisation of damping mass for coulomb solver and
!     constraint potential
    
    CALL inicms
    CALL init_constraint
    DO iwu=6,7
      IF(tswicon0) THEN
        WRITE(iwu,'(a,i3,2g12.4,/)')  &
          ' now running constraint: L,qsoll,beta=',0,qsoll0,qsoll0/cfact0
      END IF
      IF(tswicon2) THEN
        WRITE(iwu,'(a,i3,2g12.4,/)')  &
          ' now running constraint: L,qsoll,beta=',2,qsoll2,qsoll2/cfact2
      END IF
      IF(tswicon3) THEN
        WRITE(iwu,'(a,i3,2g12.4,/)')  &
          ' now running constraint: L,qsoll,beta=',3,qsoll3,qsoll3/cfact3
      END IF
      IF(tswicon4) THEN
        WRITE(iwu,'(a,i3,2g12.4,/)')  &
          ' now running constraint: L,qsoll,beta=',4,qsoll4,qsoll4/cfact4
      END IF
    END DO
    
!     the loop over temperatures in case of IPAIR=4.
!     take the loop parameters from TEMPMX,TEMPMN and TEMPW
    IF(ipair /= 4) THEN
      tempmn = zero
      tempmx = zero
      tempw  = one
    END IF
    ntempmax=NINT((tempmx-tempmn)/tempw)
    DO  itemp=0,ntempmax
      tempac=itemp*tempw+tempmn
      IF(ipair == 4) THEN
!        ggg = tempac
!        ggg = tempac
        WRITE(6,'(/a,g12.4/)') ' now running T=',tempac
        WRITE(7,'(/a,g12.4/)') ' now running T=',tempac
      END IF
!**********************************************************************
!                                                                     *
!     open the necessary files                                        *
!                                                                     *
!**********************************************************************

    CALL ini_files()
      
!**********************************************************************
!                                                                     *
!     the static iteration                                            *
!                                                                     *
!**********************************************************************
      
      eold   = energy
      istblo = 0
      itmaxact=itmax
!      IF(iflmb0==1) itmaxact=itcons2     ! prepare for later direct iteration
      iterations: DO iter=0,itmaxact
        oldoct=beta3
        IF(MOD(iter,iprint)==1) THEN
          WRITE(6,*)' Iter, Energy,  beta2,   beta3,log(D^2E),log(dE/dIt),log(dQ/dIt), zcm x0inv'
          WRITE(7,*)' Iter, Energy,  beta2,   beta3,log(D^2E),log(dE/dIt),log(dQ/dIt), zcm  x0inv'
        END IF
        IF(iter /= 0) THEN
          WRITE(6,5050)  &
              iter,energy,beta2,beta3, LOG10(MAX(devh2,smal)),&
              LOG10(MAX(ABS(econ),smal)),LOG10(MAX(ABS(qcon),smal)), zcm,x0inv
          WRITE(7,5050)  &
              iter,energy,beta2,beta3, LOG10(MAX(devh2,smal)),&
              LOG10(MAX(ABS(econ),smal)),LOG10(MAX(ABS(qcon),smal)), zcm,x0inv
          IF(check_border_density()) WRITE(6,*)'***WARNING: Grid might be too small. Density at border >1E-6!'
          CALL flush(7)
          5050 FORMAT(i5,X,f9.2,2(X,f8.3),3(3X,f7.3),X,1pg12.4,0pf7.3)
        END IF
        
!     set addition-weights for mixing density and potential,
!     first iteration needs a refresh
        
        IF(iter <= 3) THEN
          adyrho = one
          adypot = one
        ELSE
          adyrho = addrho
          adypot = addpot

        END IF
        
        CALL sstep(iter)
        CALL sinfo(iter)
        
        dnewoct=beta3
        
        IF(econ < 10.0*errenr .AND. kblock > 0 .AND. istblo <= 0)  &
            CALL comp_block
        itert  = itert+1
        IF(iter > itpini+3 .AND. tswicon) THEN
          IF (econ < errenr .AND. devh2 < errvar .AND. &
              ABS(qcon) < qerr .AND. ABS(zcm) < 1E-3  .AND. &
              ABS(oldoct-dnewoct) < 1E-3) EXIT iterations
        ELSE IF (.NOT.tswicon) THEN
          IF(econ < errenr .AND. devh2 < errvar) EXIT iterations
        END IF
        
        IF(tvaryx0) THEN
          IF(energy<econprev .AND. devh2<(devh2prev*(1.0d0-1.0d-5))) THEN
             x0inv=x0inv*1.005
          ELSE
             x0inv=x0inv*0.9
          END IF
          IF(x0inv<x0invmin) x0inv=x0invmin
          econprev=energy
          devh2prev=devh2
       END IF
        
        IF(iter > itcons0) THEN
          qlamb0  = zero
          tswicon0 = .FALSE.
        END IF
        IF(iter > itcons2) THEN
          qlamb2  = zero
          tswicon2 = .FALSE.
        END IF
        IF(iter > itcons3) THEN
          qlamb3  = zero
          tswicon3 = .FALSE.
        END IF
        IF(iter > itcons4) THEN
          qlamb4  = zero
          tswicon4 = .FALSE.
        END IF
         
      END DO iterations
      
!**********************************************************************
!                                                                     *
!     end of static iteration                                         *
!                                                                     *
!**********************************************************************
      
!     one more iteration to put observables up to date
      
!      IF(tafter .AND. iflmb0 /= 1) THEN
      IF(tafter) THEN
        ipold  = iprint
        iprint = 1
        iplold = ipri_fields
        itesto = itest
        IF(ipri_fields>0) ipri_fields  = 1
        itest  = 0
        econs  = econ
        qconsav  = qcon
        devh2s = devh2
        adyrho = 1.0
        adypot = 1.0
        CALL sstep(iter)                       ! "afterburn"
        CALL sinfo(itmax)
        iprint = ipold
        ipri_fields  = iplold
        itest  = itesto
        econ   = econs
        qcon   = qconsav
        devh2  = devh2s
      END IF
      IF(tpri_2d)  CALL plot_2d()
      CALL occ_lz()
      
      CALL user_interface_1()        ! user defined corrections
      
!     optionally compute localization criterion
      
      IF(klocal == 1) CALL localize()
      
      CALL write_ensurf(iter)

    END DO       ! temperature loop
    
  END DO       ! L=4 loop
  END DO       ! L=3 loop
  END DO       ! L=2 loop
  END DO       ! L=0 loop
  
!**********************************************************************
!                                                                     *
!     end of loop over deformations                                   *
!                                                                     *
!**********************************************************************
  
  CALL user_interface_2()  ! user defined post processing

STOP 'normal termination '


CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: write_headers
!> @brief
!! Writes input parameters and conditions of actrual run on unit 7
!> @details
!!
!--------------------------------------------------------------------------- 
SUBROUTINE write_headers()

  IMPLICIT NONE
  INTEGER,  EXTERNAL   :: omp_get_num_threads
  
  !$OMP PARALLEL
  !$ number_threads=OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  
  WRITE(7,'(/,10X,a,i3,a)')'SkyAx runs with ',number_threads,' OMP threads.'

  WRITE(7,5001) title,outnam,nmass,nchg,nz,nzd-1,nr,nrd-1,dz,dr
  5001  FORMAT(//,10X,'Title: ',a,///, 10X,'qualifier of output files=',a,/,  &
      10X,'mass number=',i5,/, 10X,'total charge=',i5,//,  &
      10X,'mesh dimensions',/, 10X,'nz =',i5,'  maximum allowed =',i5,/,  &
      10X,'nr =',i5,'  maximum allowed =',i5,//, 10X,'mesh spacings  ',/,  &
      10X,'dz =',f5.3,'  fm',/, 10X,'dr =',f5.3,'  fm',/)
  IF(tsym) THEN
    WRITE(7,'(10x,a)') 'z-reflection symmetric calculation'
  ELSE
    WRITE(7,'(10x,a)') 'calculations without z-reflection symmetry'
  END IF
  IF(idrtyp == 0) THEN
    WRITE(7,'(10x,a,2i5)') 'd/dr & d/dz from fourier.'
  ELSE IF(idrtyp == 1) THEN
    WRITE(7,'(10x,a,2i5)') 'd/dr from finite difference & d/dz from fourier.'
  ELSE IF(idrtyp == 2) THEN
    WRITE(7,'(10x,a,2i5)') 'd/dr from fourier & d/dz from finite difference.'
  ELSE IF(idrtyp == 3) THEN
    WRITE(7,'(10x,a,2i5)') 'd/dr & d/dz from finite difference.'
  ELSE
    STOP ' invalid "idrtyp"'
  END IF
  IF(TRIM(name)=='') THEN 
    WRITE(7,'(//,10X,a)') 'Using new force defined explicitly'
  ELSE
    WRITE(7,'(//,10X,a,a)') 'Using predefined force: ',TRIM(name)
  END IF
  IF(tFayans) THEN
    WRITE(7,5003) avp,h1vp,h2vp,avm,h1vm,h2vm,asp,hsp,hgradp,asm,hsm,hgradm, &
         fxi,h1xi,hgradxi,delfxi,isovectormodel,c0nabj,c1nabj,c0tau,c1tau,hcoul, &
         e2,h2m,ifzpe,ifex
    5003  FORMAT(/, 10X,'avp    =',f10.3,10X,'h1vp   =',f10.3,10X,'h2vp   =',f10.3,/, &
      10X,'avm    =',f10.3,10X,'h1vm   =',f10.3,10X,'h2vm   =',f10.3,/, &
      10X,'asp    =',f10.6,10X,'hsp    =',f10.6,10X,'hgradp =',f10.6,/, &
      10X,'asm    =',f10.6,10X,'hsm    =',f10.3,10X,'hgradm =',f10.6,/, &
      10X,'fxi    =',f10.6,10X,'h1xi   =',f10.6,10X,'hgradxi=',f10.6,/, &
      10X,'delfxi =',f10.6,10X,'isovectormodel=',i12,/,  &
      10X,'c0nabj =',f10.6,10X,'c1nabj =',f10.6,/, &
      10X,'c0tau  =',f10.6,10X,'c1tau  =',f10.6,/, &
      10X,'hcoul  =',f10.6,10X,'e2     =',f10.3,/,  &
      10X,'h2m    =',2F12.6,/,10X,'ifzpe  =',i5,15X,'ifex   =',i5)
  ELSE
    WRITE(7,5002) t0,t1,t2,t3,t4,x0,x1,x2,x3,b4p,power,e2,h2m,ifzpe,riftm
    5002  FORMAT(/, 10X,'t0     =',f10.3,'  mev*fm^3',6X,/,  &
      10X,'t1     =',f10.3,'  mev*fm^5',6X,/,  &
      10X,'t2     =',f10.3,'  mev*fm^5',6X,/,  &
      10X,'t3     =',f10.3,'  mev*fm^3+3*power',6X,/,  &
      10X,'t4     =',f10.3,'  mev*fm^4',6X,/, 10X,'x0     =',f10.6,/,  &
      10X,'x1     =',f10.6,/, 10X,'x2     =',f10.6,/,  &
      10X,'x3     =',f10.6,/, 10X,'b4p    =',f10.6,/,  &
      10X,'power  =',f10.6,/, 10X,'e2     =',f10.3,/,  &
      10X,'h2m    =',2F10.3,/, 10X,'ifzpe  =',i5,/,  &
      10X,'riftm  =',f10.3)
  WRITE(7,5005) pairfc,rho0pr,pairstab
  5005 FORMAT(10X,'pairfc  =',f12.3,'  strength for p',' (mev)',/,  &
      10X,'pairfc  =',f12.3,'  strength for n',' (mev)',/,  &
      10X,'rho0pr  =',1pg12.4,'  density parameter in DDDI',/,  &
      10X,'pairstab=',f12.3,'  BCS stabilizing pairing energy'/)
  END IF
  
  WRITE(7,5017) colerr,irepc,ndrun,nqrun,norun,nhrun,cmssqr
  5017  FORMAT(//,10X,'laplace equation parameters',//,  &
      10X,'colerr = ',g10.2,'  error for col iter',/,  &
      10X,'irepc  = ',i10,'  no. of iterations  ',/,  &
      10X,'nqrun  = ',i10,'  dipole bound (y=1/n=0)  ',/,  &
      10X,'nqrun  = ',i10,'  quad bound (y=1/n=0)  ',/,  &
      10X,'norun  = ',i10,'  oct bound (y=1/n=0)  ',/,  &
      10X,'nhrun  = ',i10,'  hex bound (y=1/n=0)  ',/,  &
      10X,'cmssqr = ',g10.2,'  mass for col damping')
  IF(tsym) THEN
    ndrun = 0
    norun = 0
  ELSE
    ndrun = 1
    norun = 1
  END IF
  
  WRITE(7,5007) chg,chgx,e0inv,x0inv,addrho,addpot,addocc,tvaryX0
  5007  FORMAT(/,  &
      10X,'chg    =',f5.1,'  effective charge for direct field',/,  &
      10X,'chgx   =',f5.1,'  effective charge for exchange field',///,  &
      10X,'pade damping parameters',//,  &
      10X,'e0inv   =',f10.3,'  damping energy scale (mev)',/,  &
      10X,'x0inv   =',g10.3,'  expansion parameter',/,  &
      10X,'addrho  =',g10.3,'  weight for accumulating density',/,  &
      10X,'addpot  =',g10.3,'  weight for accumulating potential',/,  &
      10X,'addocc  =',g10.3,'  weight for accumulating occ.numbers',/, &
      10X,'tvaryX0 =',L6,4x,'  dynamical tuning of step size x0inv',/)
 
  WRITE(7,5004) ipair,npiter,itpini
  5004 FORMAT(//,10X,'bcs pairing model parameters',//,  &
      10X,'ipair   =',i12,'  ipair=-1 nopairing',/,  &
      31X,                '  ipair=4 temperature CHF',/,  &
      31X,                '  ipair=5 BCS delta-force',/,  &
      31X,                '  ipair=7 BCS dens.dep. delta-force',//,  &
      31X,                '  ipair=9 BCS-LN delta-force',/,  &
      31X,                '  ipair=11 BCS-LN dens.dep. delta-force',//,  &
      10X,'npiter  =',i12,'  no. of pairing iterations',//,  &
      10X,'itpini  =',i12,'  no. iterations before pairing')
  WRITE(7,5014)    addocc
  5014  FORMAT(  &
      10X,'addocc  =',f12.3,'  mixing of new and old occ. numbers')
  
  IF(nr > nrd-1) THEN
    WRITE(7,6001) nr,nrd
    6001     FORMAT(' insufficient r dimension - nr =',i5,' nrd = ',i5)
    STOP 'dimension in setup'
  END IF
  
  IF(nz > nzd-1) THEN
    WRITE(7,6002) nz,nzd
    6002     FORMAT(' insufficient z dimension - nz =',i5,' nzd = ',i5)
    STOP 'dimension in setup'
  END IF
  
  IF((ipair == 5.OR.ipair == 7)) THEN
    IF (pairstab > 0.0D0) THEN
      WRITE(7,6003) pairstab
      6003    FORMAT(' BCS stabilizing energy must not be positive - pairstab=',  &
          f10.3)
      STOP 'wrong setup for stabilized BCS'
    END IF
  ELSE IF (pairstab /= 0.0D0) THEN
    WRITE(7,6004)
    6004    FORMAT(' stabilized BCS only allowed for ipair=5,7')
    STOP 'wrong setup for stabilized BCS'
  END IF
  
! switch 'iflmb0' presently only for L=2
!  IF(itcons2 >= 0) THEN
!    iflmb0 = 0
!  ELSE
!    iflmb0 = 1
!  END IF
  itcons2 = ABS(itcons2)
  
  WRITE(7,5008) nmass,nchg,deocc,betain
  5008  FORMAT(///,10X,'static input parameters',//,  &
      10X,'nmass   =',i10,'  atomic mass number',/,  &
      10X,'nchg    =',i10,'  atomic charge number',/,  &
      10X,'deocc(1)=',f10.3,'  additional no. of quanta for proton ',  &
      '  shell ',/,  &
      10X,'deocc(2)=',f10.3,'  additional no. of quanta for neutron',  &
      '  shell ',/, &
      10X,'betain  =',f10.3,'  initialization deformation parameter ')
  WRITE(7,5009) ERRenr,errvar,erra
  5009  FORMAT(/10X,'errenr  =',1pg12.4,'  error in energy convergence',//,  &
      10X,'errvar  =',1pg12.4,'  error in s.p. variance',//,  &
      10X,'erra    =',1pg13.4,'  error in particle number (pairing)')
  WRITE(7,5010) itmax,iprint,ipri_fields,itest
  5010  FORMAT(/,  &
      10X,'itmax      =',i10,'  maximum no. of static iterations',/,  &
      10X,'iprint     =',i10,'  iteration modulus for printing',/,  &
      10X,'ipri_fields=',i10,'  iteration modulus for 1D density plot'/,  &
      10X,'itest      =',i10,'  level for testprints')
  IF(ABS(dipin)>0D0) WRITE(7,'(/a,1pg12.4/)') &
    ' external isovector dipole field dipin=',dipin
  IF (ABS(swicon0)+ABS(swicon2)+ABS(swicon3)+ABS(swicon4) > 0.5D0) THEN
    WRITE(7,'(/10x,a)') 'constraint calculation:'
    IF(itypco.GE.0) THEN
      IF(qeps<0D0) THEN         ! optionally set defaults
        IF(inodiag<itmax) THEN
          qeps = 0.05D0
        ELSE
          qeps = 0.3D0
        END IF
      END IF
      WRITE(7,5011)  qerr,c0,d0,qeps
      5011  FORMAT(/,  &
        10X,'qerr    =',g10.3,'  accuracy for quadrupol deformation',/,  &
        10X,'c0      =',f10.5,'  stepsize quadr.-correction',/,  &
        10X,'d0      =',f10.5,'  parameter d0 for constraint calc  '/,  &
        10X,'qeps    =',f10.5,'  stepsize qlambdacorrection ')
    ELSE
      WRITE(7,'(10x,a)') 'using augmented Lagrangian method' 
    END IF
    IF(ABS(swicon0)>0.5D0) &
    WRITE(7,5013)  'L=0 constraint:',betmin0,betmax0,betaw0,itcons0
    IF(swicon0<-0.5D0) WRITE(7,'(10x,a)') 'fixed constraining force'
    IF(ABS(swicon2)>0.5D0) &
    WRITE(7,5013)  'L=2 constraint:',betmin2,betmax2,betaw2,itcons2
    IF(swicon2<-0.5D0) WRITE(7,'(10x,a)') 'fixed constraining force'
    IF(ABS(swicon3)>0.5D0) &
    WRITE(7,5013)  'L=3 constraint:',betmin3,betmax3,betaw3,itcons3
    IF(swicon3<-0.5D0) WRITE(7,'(10x,a)') 'fixed constraining force'
    IF(ABS(swicon4)>0.5D0) &
    WRITE(7,5013)  'L=4 constraint:',betmin4,betmax4,betaw4,itcons4
    IF(swicon4<-0.5D0) WRITE(7,'(10x,a)') 'fixed constraining force'
    5013  FORMAT(//,10X,a,/  &
        10X,'betamin =',f10.4,'  minimal beta deformation value ',/,  &
        10X,'betamax =',f10.4,'  maximal beta deformation value ',/,  &
        10X,'betastep=',f10.4,'  beta step size '/,  &
        10X,'itcons  =',i10,   '  maximum number of iterations for constraint ')
    WRITE(7,5012) dalp,dgam,tfresh
    5012  FORMAT(/,  &
        10X,'dalp   =',f10.4,'  parameter det. the beginning of damping.',/,  &
        10X,'dgam   =',f10.4,'  parameter det. the width of damping  ',/,&
        10X,'tfresh =',L1,9x,'  initialize from scratch each deformation.',/)
  ELSE
    WRITE(7,'(/10x,a)') 'calculation without constraint'
  END IF
  WRITE(7,'(///)')
  RETURN

END SUBROUTINE write_headers


END PROGRAM skyax
