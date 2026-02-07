!------------------------------------------------------------------------------
! MODULE: staticstep
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the routines to perform one static iteration and several 
!!helper routines.
!------------------------------------------------------------------------------
MODULE staticstep

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
USE hdiag
USE pairin
USE steputil
!$ use omp_lib

IMPLICIT NONE
INTEGER, PUBLIC :: inodiag=500    !< Number of iterations without diagonalizing the Hamiltonian matrix
INTEGER, PUBLIC :: iterstat      !< to communicate actual iteration number
REAL(8), PUBLIC :: devhh2              !< counter for variance (h*h)
REAL(8), PUBLIC :: hdh                 !< counter for <hD h>
REAL(8), PUBLIC :: hddh                !< counter for <hD Dh>
REAL(8), PUBLIC :: hdhdh               !< counter for <hD h Dh>
REAL(8), PUBLIC :: qh0,qh1,qh2,qh3,qh4 !< counter for <Q h>
REAL(8), PUBLIC :: qhdig               !< counter for diagonal part of <Q h>, only L=2

LOGICAL, PUBLIC,PARAMETER :: tshiftcm=.TRUE.   !< shift c.m. to z=0


PRIVATE

!*************************************************************************
!  Parameter settings
!*************************************************************************
REAL(8) :: caugmlagrb=20D0  !50D0   !< Curvature parameter for augmented Lagrangian method, to be used in connection with 'itypco=-1'. Curvature here with respect to dimensionless deformation.

INTEGER,PARAMETER :: ifqqpr=1     !< switch for detailed protocol of <QQ> etc (experienced users only)
INTEGER,PARAMETER :: ioccmx=1     !< switch for propagation of occ.number
INTEGER,PARAMETER :: ifvarh=1     !< switch for full variance
INTEGER,PARAMETER :: ifopts=0     !< protocol of optimal stepsize
LOGICAL,PARAMETER :: tcstp=.FALSE. !< protocol of constraint step
REAL(8),PARAMETER :: v2lim=1D-7   !< lower limit for inclusion in 'rho' sum


LOGICAL,SAVE       :: tvarfile=.false.  !< switch to write out a file containing variances.

REAL(8),ALLOCATABLE :: hmatr(:,:)    !< Hamiltonian matrix



!*************************************************************************
!  Other variables
!*************************************************************************

LOGICAL :: tqq              !<switch for lambda update protocoll
LOGICAL :: tdiag            !< switch for diagonalization of Hamiltonian matrix
INTEGER :: jst              !< Variable for current wf

REAL(8) :: caugmlagr0=0D0   !< Parameter for augmented Lagrangian method  L=0
!REAL(8) :: caugmlagr1=0D0   !< Parameter for augmented Lagrangian method  L=1
REAL(8) :: caugmlagr2=0D0   !< Parameter for augmented Lagrangian method  L=2
REAL(8) :: caugmlagr3=0D0   !< Parameter for augmented Lagrangian method  L=3
REAL(8) :: caugmlagr4=0D0   !< Parameter for augmented Lagrangian method  L=4
REAL(8),SAVE :: qc          !< auxiliary storage for actual deformation

!REAL(8) :: ebin                !< offset in s.p. energy

! The following are already defined in 'm_observ.f90':
!    devh2               ! counter for variance
!    devoc2              ! counter for variance (occ.)
!    devmf2              ! counter for variance (mf)
!    devgap              ! counter for gap equation
!    qdh                 ! counter for <Q Dh>
!    qq                  ! counter for <Q Q>  


REAL(8),PRIVATE,ALLOCATABLE :: psin(:,:,:) !< temporal storage for wave functions
REAL(8),PRIVATE,ALLOCATABLE :: psip(:,:,:) !< temporal storage for wave functions

INTEGER, PRIVATE :: num_threads                !< number of OMP threads
INTEGER, PRIVATE :: mythread                   !< OMP thread number of current thread
INTEGER,ALLOCATABLE,PRIVATE :: load_m(:,:)     !< array for load balancing 
REAL(8),ALLOCATABLE,PRIVATE :: loadweight(:)   !< array to measure load of mblocks
PUBLIC sstep

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: sstep
!> @brief
!!Subroutine to propagate one gradient step 
!> @details
!! First opens protocoll files and sets several variables. The main step is 
!!performed for blocks of same quantum number m separately. The Hamiltonian 
!!is applied to the wave function and the Hamiltonian matrix is constructed. 
!!All variances are accumulated. After damping the iteration step is performed.
!! One extra subiteration can be performed. Then the wave functions are 
!!diagonalized and orthonormalized. The new densities are constructed (only 
!!\c rho for constraint iterations). In the next optional step, the 
!!constraints are applied and all densities are accumulated. 
!!Orthonormalization is performed again. Now the pairing weights and 
!!occupations are updated, optionally fields are printed and finally the new 
!!mean fields are evaluated. 
!>
!> @param[in] iter
!> INTEGER, takes the current iteration number.
!--------------------------------------------------------------------------- 
SUBROUTINE sstep(iter)

USE coulomb, ONLY: wcolt
USE util, ONLY: protfold
USE fields

IMPLICIT NONE

INTEGER, INTENT(IN)           :: iter

REAL(8) :: aopt,fac,qh,redfac,width,rmsb,v2act,symcond
INTEGER :: im,iq,is,kpsic,nadd
INTEGER :: nlo,nhi,i
!$ REAL(8) :: wtime1,wtime2
LOGICAL :: topened

REAL(8), ALLOCATABLE :: temp(:,:,:,:)
REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:),wrk4(:,:)
REAL(8),ALLOCATABLE :: rhot(:,:,:), taut(:,:,:), drholst(:,:,:),rholsrt(:,:,:),rholszt(:,:,:)


ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr),wrk4(0:nz,0:nr))
ALLOCATE(hmatr(1:npsi(2),1:npsi(2)))
ALLOCATE(rhot(0:nz,0:nr,nsd),taut(0:nz,0:nr,nsd),drholst(0:nz,0:nr,nsd),rholsrt(0:nz,0:nr,nsd),rholszt(0:nz,0:nr,nsd))
wrk1=0D0;wrk2=0D0;wrk4=0D0;hmatr=0D0;rhot=0D0;taut=0D0;drholst=0D0;rholsrt=0D0;rholszt=0D0
!psip=0D0;psin=0D0;
!IF(ifqq==1 .AND. inodiag<itmax) &
!  STOP "option IFQQ=1 incompatible with mean-field diagonalization"


IF(tcstp .AND. tswicon) THEN
  INQUIRE(unit=707,OPENED=topened)
  IF(.NOT.topened) THEN
    OPEN(707,file='protoc_constraint')
    WRITE(707,'(a)') ' iter q-deviat q-actual qcorr del_lambda  lambda'
  END IF
END IF


! scaling of parameter for augmented Lagrangian method (yet tentative)
IF(itypco==-1) THEN
  rmsb   = 0.929516003D0*(nprot+nneut)**0.333333333333D0
  IF(tswicon0)  STOP "augmented Lagrangian useless for L=0"
  IF(tswicon2)  caugmlagr2=caugmlagrb*pi**2/(5D0*rmsb*rmsb*(nprot+nneut))**2
  IF(tswicon3)  caugmlagr3=caugmlagrb*(1.946773765/((nprot+nneut)*rmsb**3))**2
  IF(tswicon4)  caugmlagr4=caugmlagrb*(1.507964474/((nprot+nneut)*rmsb**4))**2
!WRITE(6,*) 'caugmlagr:', caugmlagr0, caugmlagr2, caugmlagr3, caugmlagr4
END IF


IF(itypco==-1 .AND. caugmlagrb==0D0) &
  STOP 'CAUGMLAGRB>0 required for option ITYPCO=-1'

iterstat=iter
tdiag = (iter>inodiag)      ! pre-diagonalization in sp-space


tqq = (itypco == 1 .OR. ifqqpr == 1)    ! shorthand protocol lambda update
IF(tqq .AND. tcstp .AND. tswicon) THEN
  INQUIRE(unit=66,OPENED=topened)
  IF(.NOT.topened) THEN
    OPEN(66,file='protoc_more_constraint')
    WRITE(66,'(a)') ' iter  <Qh>  <Q_damp,h>  <QQ>  <Q, H_diag> '//&
        ' <Qh>/>QQ>  lambda'
  END IF
END IF

 ALLOCATE(temp(0:nz,0:nr,nsd,npsi(2)))
 temp=0D0

! open protocol files to track variances
IF(tvarfile) THEN
  OPEN(171,FILE='variances1')
  WRITE(171,'(a,600(i2,1pg12.4))')  &
      '#',(mpsi(jst),spenrg(jst),jst=1,nst(1)+nst(2))
  OPEN(172,FILE='variances2')
  WRITE(172,'(a,600(i2,1pg12.4))')  &
      '#',(mpsi(jst),spenrg(jst),jst=1,nst(1)+nst(2))
  tvarfile = .false.
END IF

IF(.NOT.tsym .AND. iter > 2) CALL zcm_adjust()    ! update c.m. (asymmetric case)

CALL user_interface_5()   ! user defined modifications on constraint etc

!     initialize the Coulomb solver for RESPAX
IF((ipair == 9 .OR. ipair == 11) .AND. ilnrsp == 1 .AND. iter == 0) &
      wcoltt(0:nz,0:nr) = wcolt(0:nz,0:nr)

!     initialize accumulators for density, currents, etc
redfac = one-adyrho
rho(0:nz,0:nr,1:2)     = redfac*rho(0:nz,0:nr,1:2)
tau(0:nz,0:nr,1:2)     = redfac*tau(0:nz,0:nr,1:2)
drhols(0:nz,0:nr,1:2)  = redfac*drhols(0:nz,0:nr,1:2)
rholsr(0:nz,0:nr,1:2)  = redfac*rholsr(0:nz,0:nr,1:2)
rholsz(0:nz,0:nr,1:2)  = redfac*rholsz(0:nz,0:nr,1:2)
!NEW ODD
IF(toddterms) THEN
  rotjr(0:nz,0:nr,1:2)  = redfac*rotjr(0:nz,0:nr,1:2)
  rotjz(0:nz,0:nr,1:2)  = redfac*rotjz(0:nz,0:nr,1:2)
  rotspden(0:nz,0:nr,1:2)  = redfac*rotspden(0:nz,0:nr,1:2)
  spdenr(0:nz,0:nr,1:2)  = redfac*spdenr(0:nz,0:nr,1:2)
  spdenz(0:nz,0:nr,1:2)  = redfac*spdenz(0:nz,0:nr,1:2)
  d2spdenr(0:nz,0:nr,1:2)  = redfac*d2spdenr(0:nz,0:nr,1:2)
  d2spdenz(0:nz,0:nr,1:2)  = redfac*d2spdenz(0:nz,0:nr,1:2)
  curphi(0:nz,0:nr,1:2)  = redfac*curphi(0:nz,0:nr,1:2)
  spdenkr(0:nz,0:nr,1:2)  = redfac*spdenkr(0:nz,0:nr,1:2)
  spdenkz(0:nz,0:nr,1:2)  = redfac*spdenkz(0:nz,0:nr,1:2)
END IF
!END NEW ODD
IF(tswicon) rhoti(0:nz,0:nr,1:2) = redfac*rhoti(0:nz,0:nr,1:2)  ! for CHF
IF(tstab) rhosur(0:nz,0:nr,1:2)=0D0   ! for stabilizing wall

!**********************************************************************
!                                                                     *
!        here starts the big loop over states, sorted in blocks       *
!        of azimuthal angular momenta 'm' and within over             *
!        the states in each block. the wavefunction indices are       *
!        computed and charge,m,z-parity, and occupation               *
!        numbers are indexed                                          *
!                                                                     *
!**********************************************************************

!  initialize accumulators for variances
devh2 = 0D0; devoc2= 0D0; devmf2= 0D0; devhh2= 0D0; devgap= 0D0; 
hdh   = 0D0; hddh  = 0D0; hdhdh = 0D0; 
qq0=0D0;qq1=0D0;qq2=0D0;qq3=0D0;qq4=0D0;
qh0=0D0;qh1=0D0;qh2=0D0;qh3=0D0;qh4=0D0;
qh = 0D0; qdh = 0D0; qhdig = 0D0
symcondiso = 0D0
IF(istsym>0) THEN
  CALL findconjpart(1)
  CALL findconjpart(2)
!  CALL comp_block()
END IF

IF(thfb .AND. 1==2) THEN
  nlo=nst(1)+1
  nhi=nst(1)+nst(2)
  WRITE(6,'(a,i6,1pg13.5)') 'SUMUV1:',iter,SUM(SQRT(vsqf(nlo:nhi)-vsqf(nlo:nhi)**2))
END IF

IF(1>2) WRITE(6,*) 'odd dens. before:' &
  ,SQRT(SUM(rotjr(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(rotjz(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(rotspden(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenr(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenz(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(curphi(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenkr(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenkz(:,:,1:2)**2)/(nr+1)/(nz+1)) 

IF(.NOT.ALLOCATED(load_m)) CALL init_load()
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(mythread,i,im,nhi,nlo,nadd,kpsic,iq,jst,psin,psip,rhot,taut,drholst,rholsrt,rholszt,wtime1,wtime2,v2act) &
!$OMP REDUCTION(+: devhh2,qq0,qq1,qq2,qq3,qq4,qh0,qh1,qh2,qh3,qh4,qhdig,devmf2,devgap,devh2,devoc2,qdh,hdh,hddh,hdhdh,symcond, &
!$OMP           rotjr,rotjz,rotspden,spdenr,spdenz,d2spdenr,d2spdenz,curphi,spdenkr,spdenkz)
!!!$ mythread=OMP_GET_THREAD_NUM()
ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
DO  i=1,nom         ! loop over all blocks
!$  mythread=OMP_GET_THREAD_NUM()
  im=load_m(mythread,i)
  IF(im==0) CYCLE
  nhi=mhi(im)        ! uppermost state for given 'im'
  nlo=mlo(im)        ! lowest state for given 'im''
  nadd=nhi-nlo+1     ! number of states in 'im' block
  kpsic=mim(im)      ! actual 'm' value for 'im,'
  iq=miq(im)         ! actual isospin 'iq' for given 'im'
!  IF(thfb) WRITE(6,'(6(a,i4)/)') 'mythread=',mythread,' mblock=',im,  &
!       ' nlo=',nlo,' nhi=',nhi, ' lz=',kpsic,' chg=',iq
  IF(itest >= 2)  WRITE(7,'(6(a,i4)/)') 'mblock=',im,  &
       '          nlo=',nlo,'   nhi=',nhi, '   lz=',kpsic,'   chg=',iq
  !$ wtime1 = OMP_GET_WTIME()
  IF(iter > 0) THEN
!    WRITE(7,'(a,i4,2(1pg13.5))') 'pairpt:',iter,pairpt(1,1,1:2)
    IF(thfb) THEN
      CALL hfbstep(temp,symcond,nlo,nhi,iter,.TRUE.)
!      WRITE(6,'(a,4i4,1pg13.5)') 'HFB:',mythread,i,nlo,nhi,symcond
      symcondiso(iq) = symcondiso(iq) + symcond
    ELSE
    IF(tdiag) temp(0:nz,0:nr,1:nsd,nlo:nhi) = psi(0:nz,0:nr,1:nsd,nlo:nhi)
    DO  jst=nlo,nhi
      
!**********************************************************************
!                                                                     *
!        here begins the central piece of the static step:            *
!                                                                     *
!**********************************************************************
      
      ! action of Hamiltonian on wf
      CALL hpsi(psi(0:nz,0:nr,1:2,jst),psin(0:nz,0:nr,1:2),jst)        !IN:psi,jst OUT: psin
      IF(tdiag) CALL collmatr(hmatr,temp(0:nz,0:nr,1:nsd,:),&
                              psin(0:nz,0:nr,1:2),nlo,nhi,jst,npsi(2))
!      IF(tdiag) CALL collmatr(hmatr(nlo:nhi,nlo:nhi),temp(0:nz,0:nr,1:nsd,nlo:nhi),&
!                              psin(0:nz,0:nr,1:2),nlo,jst,nadd)
      CALL accumulate_variances_1(psi(0:nz,0:nr,1:2,jst),psin,iq,jst,nlo,nhi,&
           spevar(jst),spenrg(jst),&     !s.p. properties
           qq0,qq1,qq2,qq3,qq4,qh0,qh1,qh2,qh3,&      !accumulated properties
                                  qh4,qhdig,devmf2,devgap,devh2,devoc2)
      ! apply damping
      CALL damnuc(psin(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2),jst)!INPUT:psin, OUTPUT:psip
        
      CALL accumulate_variances_2(jst,nlo,nhi,psip,qdh,hdh,hddh,hdhdh)!INPUT:jst,nlo,nhi OUTPUT:qdh,hdh,hddh,hdhdh

!      IF(tcstp) WRITE(*,*) 'QQ-test1:',qq2,qq3,qh2,qh3
      ! the key iteration step
      psi(0:nz,0:nr,1:nsd,jst) = psi(0:nz,0:nr,1:nsd,jst)-psip(0:nz,0:nr,1:nsd)

       
!        CALL hpsi(psi(0:nz,0:nr,1:2,jst),psin(0:nz,0:nr,1:2),jst)
!        IF(tdiag) CALL collmatr(hmatr(nlo:nhi,nlo:nhi),temp(0:nz,0:nr,1:nsd,nlo:nhi),&
!				psin(0:nz,0:nr,1:2),nlo,jst,nhi-nlo+1)
!          CALL damnuc(psin(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2),jst)
!          psi(0:nz,0:nr,1:nsd,jst) = psi(0:nz,0:nr,1:nsd,jst)  &
!                                    -1.5D0*psip(0:nz,0:nr,1:nsd)




!      end of inner loop over states within a block of same 'm'
    END DO
    

    IF(tdiag) CALL diag_wf(hmatr(nlo:nhi,nlo:nhi),nlo,nhi-nlo+1)

    END IF
    
!   Schmidt orthonormmalization of the states  ( needed after "diag_wf")
!    IF(thfb) WRITE(6,*) 'BEFORE SCHMIT:',mythread,SUM(psi(0:nz,0:nr,1:2,nlo:nhi)**2)
    CALL schmit(psi(0:nz,0:nr,1:2,nlo:nhi),nlo-1,nadd)
  END IF
  
!   begin loop to compute the new density and s.p.energies.
!   within a block of same 'm'
  rhot=0.0d0; taut=0.0d0; drholst=0.0d0; rholsrt=0.0d0; rholszt=0.0d0
  DO jst=nlo,nhi
    iq = icharg(jst)
    v2act=adyrho*vsqf(jst)*wtpair_v2(jst)
    IF(v2act > v2lim) THEN
      IF(tswicon .AND.iter > 0.AND.iter <= itcons) THEN
        CALL sdenti(psi(0:nz,0:nr,1:2,jst),rhot(0:nz,0:nr,iq),v2act,jst)
      ELSE
        CALL densit(psi(0:nz,0:nr,1:2,jst),v2act,jst,  &
                    rhot(0:nz,0:nr,iq),taut(0:nz,0:nr,iq),drholst(0:nz,0:nr,iq),  &
                    rholsrt(0:nz,0:nr,iq),rholszt(0:nz,0:nr,iq) )
!NEW ODD
        IF(toddterms) THEN
!          WRITE(6,*) 'CUR before call:',jst,SUM(curphi(:,:,iq)**2)
          CALL densitodd(psi(:,:,:,jst),v2act,jst,  &
              rotjr(:,:,iq), rotjz(:,:,iq), rotspden(:,:,iq),  &
              spdenr(:,:,iq),spdenz(:,:,iq),d2spdenr(:,:,iq),  &
              d2spdenz(:,:,iq),curphi(:,:,iq),spdenkr(:,:,iq), spdenkz(:,:,iq))
        END IF
!END NEW ODD
      END IF
    END IF
  END DO
  !$omp critical 
  IF(tswicon .AND.iter > 0.AND.iter <= itcons) THEN
    rhoti=rhoti+rhot
  ELSE
!    IF(thfb) WRITE(6,'(a,2i4,1pg13.5)') 'HFB 4:',mythread,jst,SUM(rhot)
    rho=rho+rhot
    tau=tau+taut
    drhols=drhols+drholst
    rholsr=rholsr+rholsrt
    rholsz=rholsz+rholszt
  END IF
  !$omp end critical 

!**********************************************************************
!                                                                     *
!      here ends the first loop over all blocks                       *
!                                                                     *
!**********************************************************************
  !$ wtime2 = OMP_GET_WTIME()
  !$ loadweight(im)=wtime2-wtime1
END DO
DEALLOCATE(psip,psin)
!$OMP END PARALLEL
IF(1>2) WRITE(6,*) 'odd dens. after:'  &
  ,SQRT(SUM(rotjr(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(rotjz(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(rotspden(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenr(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenz(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(curphi(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenkr(:,:,1:2)**2)/(nr+1)/(nz+1)) &
  ,SQRT(SUM(spdenkz(:,:,1:2)**2)/(nr+1)/(nz+1)) 

IF(loadweight(1)>0.0d0.AND.num_threads>1.AND.MOD(iter,10)==1) CALL loadbalancing()

! finalize variances
IF(ifvarh == 1) THEN
  devh2  = SQRT(ABS(devh2))
  devoc2 = SQRT(ABS(devoc2))
  devmf2 = SQRT(ABS(devmf2))
ELSE
  devmf2 = SQRT(ABS(devmf2))
  devh2  = devmf2
END IF
IF(thfb) devh2=SUM(symcondiso)
devhh2 = SQRT(devhh2)

! some protocol prints
IF(ifopts == 1) THEN
  aopt = (hdh-hdh*hddh)/MAX(hdhdh-hdhdh*hddh,smal)
  WRITE(7,'(a,2g12.4)') ' optimal stepsize=',aopt
  WRITE(6,'(a,2g12.4)') ' optimal stepsize=',aopt
END IF
!IF(thfb) WRITE(7,'(a,i5,2(1pg13.5))') '#HFB symmetry condition:',iter,symcondiso
! testprint, presently restricted to constraint L=2
IF(tswicon2 .AND. tqq .AND. tcstp) &
  WRITE(66,'(a,i4,6g12.4)') 'iter,qh,qdh,qq,qhdig,qh/qq,qlam=',  &
      iter,qh,qdh,qq2,qhdig,qh/MAX(qq2,1D-10),qlamb2

!**********************************************************************
!                                                                     *
!   (optional) part for iteration of constraint                       *
!                                                                     *
!**********************************************************************
IF (tswicon .AND. iter>0 .AND. iter<=itcons) THEN
! calculation of Q-corrective step and update of Lagrange parameter
  confield(0:nz,0:nr)=0D0
  wrk4(0:nz,0:nr)=1D0      ! storage for Q correction factor
  IF(tswicon0 .AND. .NOT.tfixcon0) &
       CALL update_lambda(qfac0(0:nz,0:nr),qcfac0(0:nz,0:nr),&
                          qlamb0,qlambi0,qsoll0,caugmlagr0,&
                          qq0,qh0,wrk4(0:nz,0:nr),iter,0)
  IF(tswicon0) confield(0:nz,0:nr) = qlamb0*qfac0(0:nz,0:nr) + confield(0:nz,0:nr)
  IF(tswicon2 .AND. .NOT.tfixcon2) &
       CALL update_lambda(qfac2(0:nz,0:nr),qcfac2(0:nz,0:nr),&
                          qlamb2,qlambi2,qsoll2,caugmlagr2,&
                          qq2,qh2,wrk4(0:nz,0:nr),iter,2)
  IF(tswicon2) confield(0:nz,0:nr) = qlamb2*qfac2(0:nz,0:nr) + confield(0:nz,0:nr)
!      IF(tcstp) WRITE(*,*) 'QQ-test2:',qq2,qq3,qh2,qh3
  IF(tswicon3 .AND. .NOT.tfixcon3) &
       CALL update_lambda(qfac3(0:nz,0:nr),qcfac3(0:nz,0:nr),&
                          qlamb3,qlambi3,qsoll3,caugmlagr3,&
                          qq3,qh3,wrk4(0:nz,0:nr),iter,3)
  IF(tswicon3) confield(0:nz,0:nr) = qlamb3*qfac3(0:nz,0:nr) + confield(0:nz,0:nr)
!      IF(tcstp) WRITE(*,*) 'QQ-test3:',qq2,qq3,qh2,qh3
  IF(tswicon4 .AND. .NOT.tfixcon4) & 
       CALL update_lambda(qfac4(0:nz,0:nr),qcfac4(0:nz,0:nr),&
                          qlamb4,qlambi4,qsoll4,caugmlagr4,&
                          qq4,qh4,wrk4(0:nz,0:nr),iter,4)
  IF(tswicon4) confield(0:nz,0:nr) = qlamb4*qfac4(0:nz,0:nr) + confield(0:nz,0:nr)
!IF(tcstp) WRITE(*,*) ' RHO: ',f1int(rho(0:nz,0:nr,1)),f1int(rho(0:nz,0:nr,2))

!**********************************************************************
!                                                                     *
!        second loop over the number of m blocks                      *
!                                                                     *
!**********************************************************************
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(i,mythread,im,nhi,nlo,nadd,iq,is,jst,psin,width,fac,rhot,taut,drholst,rholsrt,rholszt) &
!$OMP REDUCTION(+: rhosur,rotjr,rotjz,rotspden,spdenr,spdenz,d2spdenr,d2spdenz,curphi,spdenkr,spdenkz)  

!$ mythread=OMP_GET_THREAD_NUM()
  DO  i=1,nom         ! loop over all blocks
    im=load_m(mythread,i)
    IF(im==0) CYCLE
    nhi=mhi(im)
    nlo=mlo(im)
    nadd=nhi-nlo+1
    iq=miq(im)
    ! Multirupole correction and ortho-normalization
    DO jst=nlo,nhi
      DO is=1,2
        psi(0:nz,0:nr,is,jst) = wrk4(0:nz,0:nr)*psi(0:nz,0:nr,is,jst)
      END DO
    END DO
    CALL schmit(psi(0:nz,0:nr,1:2,nlo:nhi),nlo-1,nadd)
    ! new densities    
    rhot=0.0d0; taut=0.0d0; drholst=0.0d0; rholsrt=0.0d0; rholszt=0.0d0
    DO jst=nlo,nhi
      iq = icharg(jst)
      v2act=adyrho*vsqf(jst)*wtpair_v2(jst)
      IF(v2act > v2lim)  THEN
          CALL densit(psi(0:nz,0:nr,1:2,jst),v2act,jst,  &
          rhot(0:nz,0:nr,iq),taut(0:nz,0:nr,iq),drholst(0:nz,0:nr,iq), &
          rholsrt(0:nz,0:nr,iq),rholszt(0:nz,0:nr,iq) )
!NEW ODD
        IF(toddterms) THEN
!        IF(toddterms .AND. 1>2) THEN
          CALL densitodd(psi(:,:,1,jst),v2act,jst,  &
              rotjr(:,:,iq), rotjz(:,:,iq), rotspden(:,:,iq),  &
              spdenr(:,:,iq),spdenz(:,:,iq),d2spdenr(:,:,iq),  &
              d2spdenz(:,:,iq),curphi(:,:,iq),spdenkr(:,:,iq), spdenkz(:,:,iq))
        END IF
!END NEW ODD
      END IF
      IF(tstab) THEN
        width = stwid/SQRT(anumb)
        fac   = EXP(-((spenrg(jst)-fer(iq))/width)**2)
        IF(fac > v2lim) CALL sdenti(psi(0:nz,0:nr,1:2,jst),rhosur(0:nz,0:nr,iq),fac,jst)
      END IF
    END DO
    !$omp critical 
    rho=rho+rhot
    tau=tau+taut
    drhols=drhols+drholst
    rholsr=rholsr+rholsrt
    rholsz=rholsz+rholszt
    !$omp end critical 
  END DO
!$OMP END PARALLEL
!IF(tcstp) THEN
!  WRITE(*,*) ' RHOT:',f1int(rhot(0:nz,0:nr,1)),f1int(rhot(0:nz,0:nr,2))
!  WRITE(*,*) ' RHO: ',f1int(rho(0:nz,0:nr,1)),f1int(rho(0:nz,0:nr,2))
!END IF
!**********************************************************************
!                                                                     *
!        end of 2. loop over the number of m blocks                   *
!                                                                     *
!**********************************************************************

  IF(tswicon0) CALL protocol_constr(qfac0(0:nz,0:nr),qcfac0(0:nz,0:nr),qsoll0,qq0,1)
  IF(tswicon2) CALL protocol_constr(qfac2(0:nz,0:nr),qcfac2(0:nz,0:nr),qsoll2,qq2,2)
  IF(tswicon3) CALL protocol_constr(qfac3(0:nz,0:nr),qcfac3(0:nz,0:nr),qsoll3,qq3,3)
  IF(tswicon4) CALL protocol_constr(qfac4(0:nz,0:nr),qcfac4(0:nz,0:nr),qsoll4,qq4,4)
  
END IF
!**********************************************************************
!                                                                     *
!     here ends to block for constraint iteration                     *
!                                                                     *
!**********************************************************************


! WRITE(7,'(a,i4,2(1pg13.5))') 'pairpt before:',iter,pairpt(1,1,1:2)
CALL update_occupations(iter)
! WRITE(7,'(a,i4,2(1pg13.5))') 'pairpt after:',iter,pairpt(1,1,1:2)


!**********************************************************************
!                                                                     *
!     new mean-field hamiltonian                                      *
!                                                                     *
!**********************************************************************

!     save coulomb density
IF(tcoulfold) THEN
!  WRITE(6,*) 'before PROTFOLD:',f1int(rho)
  CALL protonfolding(rho(0:nz,0:nr,1),crho(0:nz,0:nr))
!  CALL protfold(rho(0:nz,0:nr,1),crho(0:nz,0:nr))
!  WRITE(6,*) 'after PROTFOLD:',f1int(rho),f1int(crho)
ELSE
  crho(0:nz,0:nr) = rho(0:nz,0:nr,1)
END IF

IF(itest >= 1 .AND. MOD(iter,iprint) == 0) THEN
  WRITE(7,'(a,i3,a)') ' densities in sstep (iter=',iter,'):'
  CALL test_print(rho(0:nz,0:nr,1),0,nz_testp,0,nr_testp,'  rhop')
  CALL test_print(rho(0:nz,0:nr,2),0,nz_testp,0,nr_testp,'  rhon')
  CALL test_print(tau(0:nz,0:nr,1),0,nz_testp,0,nr_testp,'  taup')
  CALL test_print(tau(0:nz,0:nr,2),0,nz_testp,0,nr_testp,'  taun')
  IF(t4 /= 0D0) THEN
    CALL test_print(drhols(0:nz,0:nr,1),0,nz_testp,0,nr_testp,' drholsp')
    CALL test_print(drhols(0:nz,0:nr,2),0,nz_testp,0,nr_testp,' drholsn')
  END IF
END IF

! optionally stabilizing potential
IF(tstab .AND. iter > 1) THEN
  CALL stapot(wrk1(0:nz,0:nr))                  
  IF(iter == itmax) THEN
    CALL test_print(wrk1(0:nz,0:nr),0,nz,0,nr,'from STAPOT')
  END IF
END IF

! the new mean field
CALL mfhmlt(iter)
!CALL user_interface_6()    ! to replace 'mfhmlt' for alternative mean-field

DEALLOCATE(hmatr)
IF(tdiag) DEALLOCATE(temp)
DEALLOCATE(wrk1,wrk2,wrk4)
DEALLOCATE(rhot,taut,drholst,rholsrt,rholszt)

RETURN
END SUBROUTINE sstep

include 'lagpar.f90'


!---------------------------------------------------------------------------  
! DESCRIPTION: zcm_adjust
!> @brief
!!Computes actual c.m. along z-direction and updates constraint field 
!!correspondingly. Optionally, it shifts the c.m. to z=0 before all that.
!> @details
!!
!>
!> @param[in] iter
!> INTEGER, takes the current iteration number.
!--------------------------------------------------------------------------- 

SUBROUTINE zcm_adjust()

IMPLICIT NONE
  IF (tshiftcm) THEN
    zcm = zcmobs(rho(0:nz,0:nr,:))
    DO jst=1,nst(1)+nst(2)
      CALL shiftcm(zcm,jst)
!      CALL shiftcm(zcm/2D0,jst)
    END DO
  END IF
!  zcm = zcmobs(rho(0:nz,0:nr,:))
!  CALL init_constraint

END SUBROUTINE zcm_adjust


!---------------------------------------------------------------------------  
! DESCRIPTION: accumulate_variances_1
!> @brief
!! Accumulation of various expectation values and variances part 1 : 
!! Accumulation before wavefunction damping
!> @details
!! 
!--------------------------------------------------------------------------- 

SUBROUTINE accumulate_variances_1(psiin,hpsiin,iq,jst,nlo,nhi,&
             spevar,spenrg,qq0,qq1,qq2,qq3,qq4,&
             qh0,qh1,qh2,qh3,qh4,qhdig,devmf2,devgap,devh2,devoc2)

IMPLICIT NONE

REAL(8),INTENT(IN)    :: psiin(0:nz,0:nr,1:2),hpsiin(0:nz,0:nr,1:2)
INTEGER,INTENT(IN)    :: iq,jst
REAL(8),INTENT(INOUT) :: spevar,spenrg,qq0,qq1,qq2,qq3,qq4,&
                         qh0,qh1,qh2,qh3,qh4,qhdig,devmf2,devgap,devh2,devoc2
INTEGER,INTENT(IN) :: nlo,nhi     

REAL(8) :: fpair    ! intermediate weight factor      
REAL(8) :: qab0,qab1,qab2,qab3,qab4
REAL(8) :: hab,offnor,offwgt,ua,ub,va,vb,vv2
INTEGER :: is,jstp

!REAL(8),ALLOCATABLE :: psi0(:,:,:),psi1(:,:,:),psi2(:,:,:),psi3(:,:,:),&
!                       psi4(:,:,:),psip(:,:,:)
REAL(8),ALLOCATABLE :: psi0(:,:,:),psi2(:,:,:),psi3(:,:,:),psi4(:,:,:),psip(:,:,:)

ALLOCATE(psip(0:nz,0:nr,nsd))

        
!  accumulation of <Qh> and <QQ>    
  IF(tqq)  THEN
    IF(tswicon0) ALLOCATE(psi0(0:nz,0:nr,2))
    IF(tswicon2) ALLOCATE(psi2(0:nz,0:nr,2))
    IF(tswicon3) ALLOCATE(psi3(0:nz,0:nr,2))
    IF(tswicon4) ALLOCATE(psi4(0:nz,0:nr,2))
    DO is=1,nsd
      IF(tswicon0) psi0(0:nz,0:nr,is) = qfac0(0:nz,0:nr)*psiin(0:nz,0:nr,is)
      IF(tswicon2) psi2(0:nz,0:nr,is) = qfac2(0:nz,0:nr)*psiin(0:nz,0:nr,is)
      IF(tswicon3) psi3(0:nz,0:nr,is) = qfac3(0:nz,0:nr)*psiin(0:nz,0:nr,is)
      IF(tswicon4) psi4(0:nz,0:nr,is) = qfac4(0:nz,0:nr)*psiin(0:nz,0:nr,is)
    END DO
    fpair = vsqf(jst)*wate(jst)
    IF(tswicon0) qq0 =  fpair*wfovlp(psi0(0:nz,0:nr,1:2),psi0(0:nz,0:nr,1:2)) + qq0
    IF(tswicon2) qq2 =  fpair*wfovlp(psi2(0:nz,0:nr,1:2),psi2(0:nz,0:nr,1:2)) + qq2
    IF(tswicon3) qq3 =  fpair*wfovlp(psi3(0:nz,0:nr,1:2),psi3(0:nz,0:nr,1:2)) + qq3
    IF(tswicon4) qq4 =  fpair*wfovlp(psi4(0:nz,0:nr,1:2),psi4(0:nz,0:nr,1:2)) + qq4
    IF(tswicon0) qh0 =  fpair*wfovlp(psi0(0:nz,0:nr,1:2),hpsiin(0:nz,0:nr,1:2)) + qh0
    IF(tswicon2) qh2 =  fpair*wfovlp(psi2(0:nz,0:nr,1:2),hpsiin(0:nz,0:nr,1:2)) + qh2
    IF(tswicon3) qh3 =  fpair*wfovlp(psi3(0:nz,0:nr,1:2),hpsiin(0:nz,0:nr,1:2)) + qh3
    IF(tswicon4) qh4 =  fpair*wfovlp(psi4(0:nz,0:nr,1:2),hpsiin(0:nz,0:nr,1:2)) + qh4
    va    = vsqf(jst)
    ua    = SQRT(one-va)
    va    = SQRT(va)
    qhdig = (ua*ua-va*va)*2.0*ua*va*deltaf(jst) + qhdig
    DO jstp=nlo,nhi
      IF(nparz(1,jst) == nparz(1,jstp)) THEN
        hab    = wfovlp(psi(0:nz,0:nr,1:2,jstp),hpsiin(0:nz,0:nr,1:2))
        IF(tswicon0) qab0 = wfovlp(psi(0:nz,0:nr,1:2,jstp),psi0(0:nz,0:nr,1:2))
        IF(tswicon2) qab2 = wfovlp(psi(0:nz,0:nr,1:2,jstp),psi2(0:nz,0:nr,1:2))
        IF(tswicon3) qab3 = wfovlp(psi(0:nz,0:nr,1:2,jstp),psi3(0:nz,0:nr,1:2))
        IF(tswicon4) qab4 = wfovlp(psi(0:nz,0:nr,1:2,jstp),psi4(0:nz,0:nr,1:2))
        vb     = vsqf(jstp)
        ub     = SQRT(one-vb)
        vb     = SQRT(vb)
        offwgt = va*vb*(ua*ub-va*vb)*wate(jst)
        IF(tswicon0) qh0 = offwgt*hab*qab0 + qh0
        IF(tswicon2) qh2 = offwgt*hab*qab2 + qh2
        IF(tswicon3) qh3 = offwgt*hab*qab3 + qh3
        IF(tswicon4) qh4 = offwgt*hab*qab4 + qh4
        IF(tswicon0) qq0 = offwgt*qab0*qab0 + qq0
        IF(tswicon2) qq2 = offwgt*qab2*qab2 + qq2
        IF(tswicon3) qq3 = offwgt*qab3*qab3 + qq3
        IF(tswicon4) qq4 = offwgt*qab4*qab4 + qq4
      END IF
    END DO
    IF(tswicon0) DEALLOCATE(psi0)
    IF(tswicon2) DEALLOCATE(psi2)
    IF(tswicon3) DEALLOCATE(psi3)
    IF(tswicon4) DEALLOCATE(psi4)
  END IF
        
!  accumulate variance of mean-field hamiltonian
  devmf2  = vsqf(jst)*wate(jst)*spevar**2 + devmf2
  vv2     = vsqf(jst)
  devgap  = wate(jst)* ((spenrg-fer(iq))*2.0*SQRT(vv2-vv2*vv2)  &
           -gap(iq)*(one-vv2-vv2) ) + devgap
  IF(ifvarh == 1) THEN
    psip(0:nz,0:nr,1:2) = hpsiin(0:nz,0:nr,1:2)
    CALL projec(psip(0:nz,0:nr,1:2),npar(jst),nlo-1,nhi-nlo+1,offnor)
    devh2  = vsqf(jst)*wate(jst)*wfovlp(psip(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2)) + devh2
    devoc2 = vsqf(jst)*wate(jst)*offnor + devoc2
  END IF
  DEALLOCATE(psip)

END SUBROUTINE accumulate_variances_1

!---------------------------------------------------------------------------  
! DESCRIPTION: accumulate_variances_2
!> @brief
!! Accumulation of various expectation values and variances part 2 : 
!! Accumulation after wavefunction damping.
!> @details
!! Note: presently only for L=2
!--------------------------------------------------------------------------- 

SUBROUTINE accumulate_variances_2(jst,nlo,nhi,psip,qdh,hdh,hddh,hdhdh)

IMPLICIT NONE

INTEGER,INTENT(IN)     :: jst
REAL(8),INTENT(IN OUT) :: qdh,hdh,hddh,hdhdh,psip(0:nz,0:nr,nsd)
INTEGER,INTENT(IN)     :: nlo,nhi     

INTEGER :: is,jstp
REAL(8) :: hab,offnor,offwgt,qab,ua,ub,va,vb
REAL(8) :: psit(0:nz,0:nr,nsd)

  IF(tqq .AND. tswicon2)  THEN
    DO is=1,nsd
      psit(0:nz,0:nr,is) = qfac2(0:nz,0:nr)*psi(0:nz,0:nr,is,jst)
    END DO
    qdh  = vsqf(jst)*wate(jst)*wfovlp(psip(0:nz,0:nr,1:2),psit(0:nz,0:nr,1:2)) + qdh
    va   = vsqf(jst)
    ua   = SQRT(1D0-va)
    va   = SQRT(va)
    DO jstp=nlo,nhi
      IF(nparz(1,jst) == nparz(1,jstp)) THEN
        vb     = vsqf(jstp)
        ub     = SQRT(one-vb)
        vb     = SQRT(vb)
        offwgt = va*vb*(ua*ub-va*vb)*wate(jst)
        hab    = wfovlp(psi(0:nz,0:nr,1:2,jstp),psip(0:nz,0:nr,1:2))
        qab    = wfovlp(psi(0:nz,0:nr,1:2,jstp),psit(0:nz,0:nr,1:2))
        qdh    = offwgt*hab*qab + qdh
      END IF
    END DO
  END IF
        
  IF(ifopts == 1) THEN                      ! test
    CALL projec(psip(0:nz,0:nr,1:2),npar(jst),nlo-1,nhi-nlo+1,offnor)
    hdh   = vsqf(jst)*wate(jst)*wfovlp(psit(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2)) + hdh
    hddh  = vsqf(jst)*wate(jst)*wfovlp(psip(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2)) + hddh
    CALL hpsi(psip(0:nz,0:nr,1:2),psit(0:nz,0:nr,1:2),jst)
    hdhdh = vsqf(jst)*wate(jst)*wfovlp(psit(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2)) + hdhdh
  END IF

END SUBROUTINE accumulate_variances_2

!---------------------------------------------------------------------------  
! DESCRIPTION: update_lambda
!> @brief
!! Calculation of Q-corrective step and update of Lagrange parameter \c qlambda.
!> @details
!! Various strategies to deduce 'qcorr' and 'qlambda' can be chosen: 
!!  - itypco = 0: difference of damped Q, local estimate of Delta Q (recommended)
!!  - itypco = 1: difference of damped Q, exact variances of damped Q
!!  - itypco = 2: from raw Q and local estimate of raw Delta Q  
!!  - itypco = -1: augmented Lagrangian, requires proper setting of 'caugmlagr', still in testing stage.
!! The routine updates the corrective array "wrk4".
!! Note: option 'itypco=1' requires smaller 'c0', typically around 0.2-0.3.
!>
!> @param[in] qfac
!> REAL(8), array, takes the actual damped constraint.
!> @param[in] qcfac
!> REAL(8), array, takes the actual raw constraint.
!> @param[in,out] qlamb
!> REAL(8), takes and returns the actual Lagrange parameter.
!> @param[in,out] qlambi
!> REAL(8), takes and returns the intrinsic Lagrange parameter (augmented Lagrangian).
!> @param[in] qsoll
!> REAL(8), takes the wanted value of deformation.
!> @param[in] qq
!> REAL(8), takes the variance of Q operator.
!> @param[in] qh
!> REAL(8), takes the variance <Qh>.
!> @param[in,out] corrfield
!> REAL(8), array, takes and returns the w.f. correction.
!> @param[in] iter
!> INTEGER, takes the current iteration number.
!--------------------------------------------------------------------------- 

SUBROUTINE update_lambda(qfac,qcfac,qlamb,qlambi,qsoll,caugmlagr,&
                         qq,qh,corrfield,iter,lcon)

IMPLICIT NONE

  REAL(8),INTENT(IN)    :: qfac(0:nz,0:nr)
  REAL(8),INTENT(IN)    :: qcfac(0:nz,0:nr)
  REAL(8),INTENT(INOUT) :: qlamb          
  REAL(8),INTENT(INOUT) :: qlambi
  REAL(8),INTENT(IN)    :: qsoll
  REAL(8),INTENT(IN)    :: caugmlagr
  REAL(8),INTENT(IN)    :: qq
  REAL(8),INTENT(IN)    :: qh
  REAL(8),INTENT(INOUT) :: corrfield(0:nz,0:nr)
  INTEGER,INTENT(IN)    :: iter
  INTEGER,INTENT(IN)    :: lcon

  REAL(8) :: qtil           !<
  REAL(8) :: qtilc          !<
  REAL(8) :: qcorr
  REAL(8),save :: qlambn,qlstep,qlcstp,qlambdif
  
  REAL(8) :: qtil2,qtilc2,x0eff

!WRITE(*,*) 'update_lambda:',lcon,qlamb,qlambi,qsoll,caugmlagr,qq,qh
!  IF(tcstp) WRITE(6,'(a,1pg13.5)')  'norm RHOTI:',f1int(rhoti(0:nz,0:nr))
  CALL obsti(rhoti(0:nz,0:nr,:),qfac(0:nz,0:nr),qcfac(0:nz,0:nr),qtil,qtilc,qtil2,qtilc2) ! new mulipole momenta
  IF(tcstp) WRITE(6,'(a,i2,6(1pg12.4))')  &
      ' OBSTI trial: lcon,q,qc,q2,qc2,qq,qh=', &
       lcon,qtil,qtilc,qtil2,qtilc2,qq,qh

! direct correction of expectation value of constraint
  x0eff  = x0inv/2.5D0       ! factor 1/2.5D0 empirical, may still be improved
  qlcstp = - (MAX(e0inv,one)/x0eff)  &
      *(qtilc-qcst)/(two*ABS(qtilc2-qtilc*qtilc/annum)+d0)
!     &          / adyrho
  qlstep = - (MAX(e0inv,one)/x0eff)  &
      *(qtil-q)/(two*ABS(qtil2-qtil*qtil/annum)+d0)

! update Lagrange parameter 'qlamb'
  IF(itypco == 1) THEN
    varqsq = qq
    qcorr  = c0*(qtilc-qsoll)/(two*varqsq+d0)
    qlambn = qlamb+qeps*qh/qq
    qlcstp = qh/qq
    qlambdif=qlambn-qlamb
    IF(tcstp) WRITE(707,'(i5,5(1pg13.5))') &
      iter,qtilc-qsoll,qtilc,qcorr,qlambdif,qlamb
  ELSE IF(itypco == 0) THEN
    varqsq = ABS(qtilc2-qtilc*qtilc/annum)
    qcorr  = c0*(qtilc-qsoll)/(two*varqsq+d0)
    qlambn = qlamb+qeps*qlcstp
    qlambdif=qlambn-qlamb
    IF(tcstp) WRITE(707,'(i5,5(1pg13.5))') &
      iter,qtilc-qsoll,qtilc,qcorr,qlambdif,qlamb
  ELSE IF(itypco == 2) THEN
    varqsq = ABS(qtil2-qtil*qtil/annum)
    qcorr  = c0*(qtil-qsoll)/(two*varqsq+d0)
    qlambn = qlamb+qeps*qlstep
  ELSE IF(itypco == -1) THEN
    varqsq = ABS(qtilc2-qtilc*qtilc/annum)
    qcorr  = c0*(qtilc-qsoll)/(two*varqsq+d0)
    qlambdif=caugmlagr*(qtilc-qsoll)
    qlamb  = qlambi-qlambdif
    qlambi = qlambi-c0*qlambdif   ! intrinsic 'lambda' damped  ?
    IF(tcstp) WRITE(707,'(2i5,5(1pg13.5))') &
      lcon,iter,qtilc-qsoll,qtilc,qlambdif,qlamb
  ELSE
    STOP ' invalid option ITYPCO'
  END IF
!                           limit step in 'qlamb'
  IF(itypco /= -1) THEN
    IF(itert /= iter .AND. iter < 10 .AND. qlamb /= 0D0) THEN
      qlstep = (qlambn-qlamb)/qlamb
      qlstep = SIGN(MIN(ABS(qlstep),qeps),qlstep)
      qlstep = qlamb*qlstep
      qlamb  = qlstep+qlamb
    ELSE
      qlamb  = qlambn
    END IF
  END IF
  IF(tcstp) WRITE(6,'(a,i2,4(1pg12.4))')  &
      ' steps: lcon,qcorr,qlstep,qlamb=',lcon,qcorr,qlstep,qlamb

!! accumulate total constraining potential
!  confield(0:nz,0:nr) = qlamb*qfac(0:nz,0:nr) + confield(0:nz,0:nr)
!     prepare stepping factor for quadrupole
  IF(itypco /= 3) corrfield(0:nz,0:nr) = &
            EXP(-qcorr*qfac(0:nz,0:nr))*corrfield(0:nz,0:nr)

  qcorrs(lcon)=qcorr
  IF(itypco == 2) THEN
    qlcstps(lcon)=qlstep
  ELSE
    qlcstps(lcon)=qlcstp
  END IF
  varqscs(lcon)=varqsq
  RETURN

END SUBROUTINE update_lambda

!---------------------------------------------------------------------------  
! DESCRIPTION: protocol_constr
!> @brief
!! Computes new multipole momenta and change thereof. Actual density comes 
!!through 'rhoti'
!> @details
!!
!>
!> @param[in] qfac
!> REAL(8), array, takes the actual damped constraint.
!> @param[in] qcfac
!> REAL(8), array, takes the actual raw constraint.
!> @param[in,out] qlamb
!> REAL(8), takes and returns the actual Lagrange parameter.
!> @param[in] qsoll
!> REAL(8), takes the wanted value of deformation.
!--------------------------------------------------------------------------- 

SUBROUTINE protocol_constr(qfac,qcfac,qsoll,qq,lcon)

IMPLICIT NONE

  REAL(8),INTENT(IN)    :: qfac(0:nz,0:nr)
  REAL(8),INTENT(IN)    :: qcfac(0:nz,0:nr)
  REAL(8),INTENT(IN)    :: qq
  REAL(8),INTENT(INOUT) :: qsoll
  INTEGER,INTENT(IN)    :: lcon

  REAL(8) :: qtil
  REAL(8) :: qtilc
  REAL(8) :: qtil2
  REAL(8) :: qtilc2 
  REAL(8) :: qsq
  REAL(8) :: qsqcst
  REAL(8) :: qlstep         
  
  REAL(8) :: qcorr,qlcstp,varqsc

  qcorr=qcorrs(lcon);qlcstp=qlcstps(lcon);varqsq=varqscs(lcon)
  
  CALL obsti(rhoti(0:nz,0:nr,:),qfac(0:nz,0:nr),qcfac(0:nz,0:nr),qtil,qtilc,qtil2,qtilc2) ! new multipole momenta
  IF(tcstp) WRITE(6,'(a,4(1pg12.4))')  &
      ' OBSTI trial: q,qc,q2,qc2=',qtil,qtilc,qtil2,qtilc2

  ! new multipole momenta
  CALL obsti(rho(0:nz,0:nr,:),qfac(0:nz,0:nr),qcfac(0:nz,0:nr),q,qcst,qsq,qsqcst)            
  IF(tcstp) WRITE(6,'(a,4(1pg12.4))')  &
      ' OBSTI  step: q,qc,q2,qc2=',q,qcst,qsq,qsqcst
  qold  = q
  qc  = qcst
  IF(itypco == 2) THEN
    qcon   = ABS((q-qsoll)/MAX(varqsq,smal))
  ELSE
    qcon   = ABS((qcst-qsoll)/MAX(varqsq,smal))
  END IF
  qcons(lcon)=qcon
  ! test prints
  IF(itest >= 1) THEN
    WRITE(7,'(a,5(1pg10.3))') ' qtil,q=',qtil,q
    WRITE(6,'(a,5(1pg10.3))') ' qtil,qold,q,varqsq,varex=',  &
        qtilc,qold,qc,qsqcst-qc*qc/annum  !,-0.5*(qc-qtilc)/qcorr
    IF(itypco == 1) THEN
      varqsc = ABS(qtilc2-qtilc*qtilc/annum)
      varqsq = ABS(qtil2-qtil*qtil/annum)
      WRITE(73,'(1x,8(1pg10.2))')  &
          qtilc,qcst,varqsc,qq,varqsc/qq,qdh/qq,qlcstp,qlcstp*qq/qdh
      WRITE(74,'(1x,8(1pg10.2))')  &
          qtil,q,varqsq,qq,varqsq/qq,qdh/qq,qlstep,qlstep*qq/qdh
      WRITE(6,'(a,5(1pg10.3))') ' qq,qdh,qdh/qq,qlstep,ratio=',  &
          qq,qdh,qdh/qq,qlstep,qlstep*qq/qdh
    ELSE
      WRITE(73,'(1x,8(1pg10.2))') qtilc,qcst,varqsc,0D0,0D0,0D0,qlcstp,0D0
      WRITE(74,'(1x,8(1pg10.2))') qtil,q,varqsq,0D0,0D0,0D0,qlstep,0D0
    END IF
  END IF

END SUBROUTINE protocol_constr

!---------------------------------------------------------------------------  
! DESCRIPTION: update_occupations
!> @brief
!! Calculate occupation numbers and ordering of the levels. \c iofk is the 
!! pointer to the 1st index of state q. 
!> @details
!!
!--------------------------------------------------------------------------- 

SUBROUTINE update_occupations(iter)

USE fyfunc, ONLY: fy_hfbpot

IMPLICIT NONE 

INTEGER, INTENT(IN) :: iter   ! actual iteration number

INTEGER :: nsttot,iofk,iq,ist2,na   
REAL(8) :: addact,conocc,deltaa,distmx,epairn,epairp,v2,vol

nsttot=nst(1)+nst(2)

! shift s.p. energies to initialize bubble shape
IF(einitbubble /= 0D0 .AND. iter <= itpini) THEN
  DO jst=1,nst(1)+nst(2)
    spenrg(jst) = spenrg(jst)-einitbubble*mpsi(jst)**2
  END DO
  WRITE(7,*) ' einitbubble=',einitbubble
END IF

! reorder and reoccupy states for output
!IF(thfb)  WRITE(6,'(a,i6,1pg13.5)') 'SUMUVb:',iter,SUM(SQRT(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2))
IF(ipair==-1 .OR. (ipair>=0.AND.itert<=itpini)) THEN
  DO iq=1,2
    iofk=iokk(iq)
    CALL vsqtr(vsqf(iofk),wate(iofk),nst(iq),  &
        spenrg(iofk),npart(iq),iorder(iofk))
  END DO
ELSE
  DO iq=1,2
    iofk=iokk(iq)
    CALL update_energy_order(nst(iq),spenrg(iofk),iorder(iofk))
  END DO
END IF

! save occupation numbers for all coming steps
IF(iter == itpfix) oc(1:nsttot) = vsqf(1:nsttot)  

IF(iter > itpfix) THEN
  vsqf(1:nsttot) = oc(1:nsttot)               ! restore frozen occupation
!IF(thfb)  WRITE(6,'(a,i6,1pg13.5)') 'SUMUVoc:',iter,SUM(SQRT(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2))
ELSE IF(ipair >= 0.AND.itert >= itpini) THEN
!IF(thfb)  WRITE(6,'(a,i6,1pg13.5)') 'SUMUVa:',iter,SUM(SQRT(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2))
  vsqfol(1:nsttot) = vsqf(1:nsttot)
  IF(ipair.GE.0 .AND. ipair.LE.3) THEN
    CALL user_interface_4()  ! switches to user defined pairing
  ELSE IF(ipair==4) THEN !Temperature
    CALL thermal_occ(tempac,spenrg,vsqf,wate,wtpair_uv,wtpair_v2,&
                     fer(1),one*npart(1),nst(1))
    ist2 = nst(1)+1
    CALL thermal_occ(tempac,spenrg(ist2),vsqf(ist2),wate(ist2),wtpair_uv(ist2),&
                     wtpair_v2(ist2),fer(2),one*npart(2),nst(2))
!  ELSE IF(ipair == 5 .OR. ipair == 7 .OR. ipair == 9  &
!        .OR. ipair == 11) THEN     ! HFB/BCS with zero-range pairing functional
  ELSE IF(ipair>4 .AND. ipair.NE.10) THEN  ! HFB/BCS with zero-range pairing functional
    CALL hfbpot(iter)
    epairp=epair(1)
    epairn=epair(2)
!    WRITE(7,'(a,2(1pg13.5))') 'after HFBPOT:',epair
    IF((ipair == 9 .OR. ipair == 11).AND. ilnrsp == 1 .AND. .NOT.tFayans) THEN
      CALL respax(iter)
      IF(itest > 0) WRITE(7,'(a,4g12.4)')  &
          ' RESPAX: rspln(1/2),suv2r2(1/2)=',rspln,suv2rs
    END IF
    itest = -1
!NEW ODD
    CALL combine_spenergies(spenrg,speold,nst(1)+nst(2))
    CALL pairdn(speold,deltaf,vsqf, wate,wtpair_uv,wtpair_v2,one*npart(1),  &
        fer(1),epair(1),gap(1),feradd(1),  &
        nst(1),nst(1),npiter,ipair,0,istblo(1),istcblo(1),erra,itest-1)
    ist2 = nst(1)+1
    CALL pairdn(speold(ist2),deltaf(ist2),vsqf(ist2),wate(ist2),  &
        wtpair_uv(ist2),wtpair_v2(ist2),one*npart(2), fer(2),epair(2),gap(2),feradd(2),  &
        nst(2),nst(2),npiter,ipair,nst(1),istblo(2)-nst(1),istcblo(2)-nst(1),erra,itest-1)
!END NEW ODD
!    WRITE(6,'(a,f9.5,200(/10f8.5))') 'neutron VSQF:',&
!         SUM(vsqf(ist2:nsttot)),(vsqf(jst),jst=ist2,nsttot)
    vsqobs(1:nsttot) = vsqf(1:nsttot)
!    IF(thfb .AND. iter>2*iterbcs .AND. MOD(iter,modnat)==0) THEN
    IF(thfb .AND. iter>iterbcs .AND. MOD(iter,modnat)==0) THEN
      IF(tFayans) THEN
        CALL fy_hfbpot()
        epairp=epair(1)
        epairn=epair(2)
      END IF
      CALL hfb_natural(1,fer(1),1D0*npart(1))
      CALL hfb_natural(2,fer(2),1D0*npart(2))
!      vsqobs(1:nsttot) = vsqf(1:nsttot)
    END IF
    IF(tFayans) THEN
      epair(1)=epairp
      epair(2)=epairn
    END IF
  ELSE
    STOP ' ipair outside -1..4 or 5, 7, 9 or 11 not active'
  END IF

  ! mixing old and new occupations
  IF(itert > itpini) THEN
  !                                      maximum distance in 'vsqf'
     IF(ioccmx == 1) THEN
      distmx = 0D0
      DO jst=1,nst(1)+nst(2)
        distmx = MAX(ABS(vsqf(jst)-vsqfol(jst)),distmx)
      END DO
      distmx = MAX(addocc,distmx)
      addact = addocc/distmx
    ELSE
      addact = addocc
    END IF
    conocc = 1D0-addact
!      WRITE(6,'(a,i6,1pg13.5)') 'SUMUVunmix:',iter,SUM(SQRT(vsqf(ist2:nst(1)+nst(2))-vsqf(ist2:nst(1)+nst(2))**2))
!      WRITE(6,'(a,i6,1pg13.5)') 'SUMUVold:',iter,SUM(SQRT(vsqfol(ist2:nst(1)+nst(2))-vsqfol(ist2:nst(1)+nst(2))**2))
    DO jst=1,nst(1)+nst(2)
      v2=vsqf(jst)
      vsqf(jst)   = addact*vsqf(jst)+conocc*vsqfol(jst)
!      IF(iter>800 .AND. jst>nst(1)) WRITE(6,'(2i4,4(1pg13.5))') jst,iorder(jst)+nst(1),spenrg(jst),v2,vsqfol(jst),vsqf(jst)
    END DO
!      WRITE(6,'(a,i6,1pg13.5)') 'SUMUVmix:',iter,SUM(SQRT(vsqf(ist2:nst(1)+nst(2))-vsqf(ist2:nst(1)+nst(2))**2))
    ! compute pairing energies from gaps and occupations
    ! (not best choice  --> replace by direct integration)
!    IF(ipair >= 5) THEN      ! disable: deltaf contains stabilization
    IF(ipair >= 12) THEN
      epairp = 0.0
      epairn = 0.0
      DO na=1,nst(1)+nst(2)
        v2     = vsqf(na)
        deltaa = deltaf(na)*wtpair_uv(na)
        vol    = 0.5*wate(na)*SQRT(MAX(v2-v2*v2,smal))*wtpair_v2(na)
        IF(na <= nst(1)) THEN
          epairp = vol*deltaa + epairp
        ELSE
          epairn = vol*deltaa + epairn
        END IF
      END DO
      epair(1) = epairp
      epair(2) = epairn
    END IF
  END IF
END IF

END SUBROUTINE update_occupations

!---------------------------------------------------------------------------  
! DESCRIPTION: init_load
!> @brief
!!Initialization of workload for OMP threads. Here every thread gets an equal 
!!number of loop iterations
!--------------------------------------------------------------------------- 
SUBROUTINE init_load()
IMPLICIT NONE
INTEGER :: i
mythread=0
num_threads=1
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mythread)
!$ num_threads=OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
ALLOCATE(load_m(0:num_threads-1,nom),loadweight(nom))
load_m=0
loadweight=-1.0d0
WRITE(*,*) 'INIT_LOAD:',num_threads,nom
DO i=1,nom
  CALL add_block(load_m,MOD(i-1,num_threads),i)
END DO
!WRITE(6,*) 'LOAD_M:'
!DO i=1,nom
!  WRITE(6,'(200i3)') load_m(:,i)
!END DO

END SUBROUTINE init_load

!---------------------------------------------------------------------------  
! DESCRIPTION: loadbalancing
!> @brief
!!Load balancing from actual loop iteration times. The algorithm starts with 
!!the most time consuming loop iteration and distributes it to the thread with 
!!the smallest load and repeats it until all tasks are distributed.
!--------------------------------------------------------------------------- 
SUBROUTINE loadbalancing()
IMPLICIT NONE
REAL(8) :: load_curr(num_threads)
INTEGER :: maxtime_ind,minload_ind

load_curr=0.0d0
load_m=0
IF(MINVAL(loadweight)<0.0d0) STOP 'ERROR IN LOADBALANCING'
DO WHILE(MAXVAL(loadweight)>0.0d0)
  maxtime_ind=MAXLOC(loadweight, DIM=1)
  minload_ind=MINLOC(load_curr, DIM=1)
  CALL add_block(load_m,minload_ind-1,maxtime_ind)
  load_curr(minload_ind)=load_curr(minload_ind)+loadweight(maxtime_ind)
  loadweight(maxtime_ind)=-1
END DO
END SUBROUTINE loadbalancing

!---------------------------------------------------------------------------  
! DESCRIPTION: add_block
!> @brief
!!Distributes a loop iteration to a certain thread.
!>
!> @param[in,out] load
!> INTEGER, array, takes variable which keeps track of loop iterations for each
!! thread.
!> @param[in] thread
!> INTEGER, takes thread number to which the loop iteration shall be distributed.
!> @param[in] mblock
!> INTEGER, takes the loop iteration which is to be distributed.
!--------------------------------------------------------------------------- 
SUBROUTINE add_block(load,thread,mblock)
INTEGER,INTENT(INOUT) :: load(0:num_threads-1,nom)
INTEGER,INTENT(IN) :: thread,mblock
INTEGER :: i
DO i=1,nom
  IF(load(thread,i)==0) THEN
    load(thread,i)=mblock
    RETURN
  END IF
END DO
STOP 'ERROR IN ADD_BLOCK'
END SUBROUTINE

END MODULE staticstep
