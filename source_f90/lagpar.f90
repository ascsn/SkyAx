SUBROUTINE hfbstep(hfbpsi,symcond,istmin,istmax,iter,tqq)

!
!     Subroutine to compute the action of the HFB Hamiltonian.
!     Computes the matrix of HFB Lagrange parameters and
!     with that the action of the HFB Hamiltonian on the
!     s.p. wavfunctions for both types of nucleons separately.
!
!      hfbpsi  = workspace for Hamiltonian on wavefunction
!      symcond = accumulator for symmetry condition
!      istmin  = minimal state nr.
!      istmax  = maximal state nr.
!      iter    = actual number of iteration
!      tqq     = switch to compute variances
!
USE paramr
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE observ, ONLY: qq0,qq1,qq2,qq3,qq4,qdh,devmf2,devgap,devh2,devoc2
!USE fyfunc, ONLY: fy_hfbpot
!USE init, ONLY: idanuc
!USE staticstep, ONLY: qh0,qh1,qh2,qh3,qh4,hdh,hddh,hdhdh,qhdig

!INCLUDE 'resden.inc'
!INCLUDE 'crankin.inc'


REAL(8), INTENT(OUT)    :: hfbpsi(0:nz,0:nr,nsd,npsi(2))
REAL(8), INTENT(OUT)    :: symcond
INTEGER, INTENT(IN)     :: istmin,istmax,iter
LOGICAL, INTENT(IN)     :: tqq


REAL(8), ALLOCATABLE :: hlagr(:,:)
REAL(8), ALLOCATABLE :: symatrix(:,:)

!REAL(8),PARAMETER :: xsmall=1D-99

LOGICAL, PARAMETER :: tbcs=.FALSE.    ! testswitch
LOGICAL, PARAMETER :: ttest=.FALSE.   ! testswitch
REAL(8) :: v2a,uva,fachfb
INTEGER :: mist,npzist,jst,ist,is,iq,mythread

REAL(8),ALLOCATABLE :: psip(:,:,:),psin(:,:,:)

IF(thfbfold) ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))

!  iq=OMP_GET_THREAD_NUM()
!  WRITE(6,*) 'HFBSTEP: mythread=',iq

! set isospin (protons or neutrons)
IF(istmin <= nst(1)) THEN
  iq = 1
ELSE
  iq = 2
END IF
!  CALL fy_hfbpot()
!DO ist=0,nr
!  WRITE(6,'(200(1pg13.5))') pairpt(0:nz,ist,iq)
!END DO
!STOP "testp"
!IF(iq==2) WRITE(6,'(a,1pg13.5)') 'PAIRPT:',pairpt(1,1,iq)
symcond = 0D0
!$  mythread=OMP_GET_THREAD_NUM()
!  WRITE(6,*) 'HFBSTEP 1: mythread=',mythread
!     compute and store action of HFB Hamiltonian
IF(iterbcs==0) fachfb = MIN(1D0*iter/iterhfb,1D0)
!WRITE(7,'(a,1pg13.5)') 'FACHFB=',fachfb
DO jst=istmin,istmax
!  WRITE(6,*) 'HFBSTEP 2: mythread,jst=',mythread,jst
  CALL hpsi(psi(0:nz,0:nr,1:2,jst),hfbpsi(0:nz,0:nr,1:2,jst),jst)
  v2a = vsqf(jst)
  uva = SQRT(MAX(v2a-v2a*v2a,xsmall))*wtpair_uv(jst)
! augment for BCS cooling in case of 'iterhfb=0'
  IF(iterbcs==0) uva = uva * fachfb
!  WRITE(6,*) 'HFBSTEP 3: mythread,jst=',mythread,jst,v2a,uva
  IF(tbcs) THEN
    v2a=1D0
    uva=0D0
    IF(ttest) WRITE(7,'(a,i4,1pg13.5)') 'ist,spenrg=',jst,spenrg(jst)
    hfbpsi(0:nz,0:nr,1:2,jst) = hfbpsi(0:nz,0:nr,1:2,jst)  &
                                -psi(0:nz,0:nr,1:2,jst)*spenrg(jst)
    IF(ttest) THEN
      WRITE(7,'(a,i4)') 'ortho-normality before: jst=',jst
      WRITE(7,'(100(1pg13.5))')  &
          (wfovlp(psi(0:nz,0:nr,1:2,ist),psi(0:nz,0:nr,1:2,jst)),ist=istmin,istmax)
    END IF
  ELSE
    IF(thfbfold) THEN
      psip(0:nz,0:nr,1)=pairpt(0:nz,0:nr,iq)*psifold(0:nz,0:nr,1,jst)
      psip(0:nz,0:nr,2)=pairpt(0:nz,0:nr,iq)*psifold(0:nz,0:nr,2,jst)
      CALL Gaussfolding_spin(psip,psin,nparz(1,jst),mpsi(jst),1)
      DO is=1,2
        hfbpsi(0:nz,0:nr,is,jst) = hfbpsi(0:nz,0:nr,is,jst)*v2a  &
            +psin(0:nz,0:nr,is)*uva
      END DO
    ELSE
      DO is=1,2
        hfbpsi(0:nz,0:nr,is,jst) = hfbpsi(0:nz,0:nr,is,jst)*v2a  &
            +pairpt(0:nz,0:nr,iq)*psi(0:nz,0:nr,is,jst)*uva
      END DO
    END IF
  END IF
!  WRITE(6,'(a,i4,4(1pg13.5))') 'HPSI:',jst,v2a,uva, &
!     wfovlp(psi(0:nz,0:nr,1:2,jst),hfbpsi(0:nz,0:nr,1:2,jst)),SUM(pairpt(0:nz,0:nr,iq))
  IF(tqq) CALL accumulate_variances_1(psi(0:nz,0:nr,1:2,jst),hfbpsi,iq,jst,istmin,istmax,&
                            spevar(jst),spenrg(jst),&     
                            qq0,qq1,qq2,qq3,qq4,qh0,qh1,qh2,qh3,&  
                            qh4,qhdig,devmf2,devgap,devh2,devoc2)
!  WRITE(6,*) 'HFBSTEP 4: mythread,jst=',mythread,jst,qq0,qq2
END DO
IF(thfbfold) DEALLOCATE(psip,psin)


!     compute matrix of Lagrange parameters
ALLOCATE(hlagr(istmin:istmax,istmin:istmax))
ALLOCATE(symatrix(istmin:istmax,istmin:istmax))
hlagr = 0D0
symatrix = 0D0
DO ist=istmin,istmax
  mist = mpsi(ist)
  npzist = nparz(1,ist)
!  WRITE(6,*) 'HFBSTEP 5: mythread,ist=',mythread,ist,mist,npzist
  DO jst=istmin,ist
    IF(mpsi(jst)==mist .AND. nparz(1,jst)==npzist) THEN
      hlagr(ist,jst) = wfovlp(psi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,jst))
      symatrix(ist,jst) = hlagr(ist,jst)
      IF(ist /= jst) THEN
        hlagr(jst,ist) = wfovlp(psi(0:nz,0:nr,1:2,jst),hfbpsi(0:nz,0:nr,1:2,ist))
        symatrix(ist,jst) = hlagr(ist,jst)-hlagr(jst,ist)
        symatrix(jst,ist) = symatrix(ist,jst)
        hlagr(ist,jst) = 0.5D0*(hlagr(ist,jst)+hlagr(jst,ist))
        hlagr(jst,ist) = hlagr(ist,jst)
      END IF
    ELSE
      hlagr(ist,jst) = 0D0
    END IF
!    WRITE(6,'(a,2i4,2(1pg13.5))') 'sym:',ist,jst,symatrix(ist,jst),hlagr(ist,jst)
  END DO
END DO
!STOP "test"

!     accumulate symmetry condition
DO ist=istmin,istmax
  DO jst=istmin,ist-1
    symcond = symcond + symatrix(ist,jst)**2
  END DO
END DO

!     subtract Lagrange matrix from 'hfbpsi'
IF(tbcs) hlagr=0D0
DO ist=istmin,istmax
  mist = mpsi(ist)
  npzist = nparz(1,ist)
  DO jst=istmin,istmax
    IF(mpsi(jst)==mist .AND. nparz(1,jst)==npzist)  &
        hfbpsi(0:nz,0:nr,1:2,ist) = hfbpsi(0:nz,0:nr,1:2,ist)  &
        -psi(0:nz,0:nr,1:2,jst)*hlagr(jst,ist)
  END DO
  IF(.NOT.tbcs) devmf2  = devmf2 + vsqf(ist)*wate(ist)*  &
                          wfovlp(hfbpsi(:,:,:,ist),hfbpsi(:,:,:,ist))
END DO
DEALLOCATE(hlagr,symatrix)
!  WRITE(6,*) 'HFBSTEP 6: mythread,devmf2=',mythread,devmf2


!     kinetic energy damping
DO ist=istmin,istmax
  IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
    WRITE(7,'(a,i4,3(1pg13.5))') 'HFB overlaps before:',  &
        ist,wfovlp(psi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(hfbpsi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(psi(0:nz,0:nr,1:2,ist),psi(0:nz,0:nr,1:2,ist))
  END IF
!  WRITE(6,*) 'HFBSTEP 7: mythread,norm=',mythread,SUM(hfbpsi(0:nz,0:nr,1:2,ist)**2)
  CALL damp_hfb(hfbpsi(:,:,:,ist),vsqf(ist),ist)
  IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
    WRITE(7,'(a,i4,3(1pg13.5))') 'HFB overlaps after:',  &
        ist,wfovlp(psi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(hfbpsi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(psi(0:nz,0:nr,1:2,ist),psi(0:nz,0:nr,1:2,ist))
  END IF
  IF(tqq) CALL accumulate_variances_2(ist,istmin,istmax,hfbpsi(:,:,:,ist),qdh,hdh,hddh,hdhdh)  
!  WRITE(6,*) 'HFBSTEP 8: mythread,v2=',mythread,vsqf(ist),&
!   SUM(psi(0:nz,0:nr,1:2,ist)**2),SUM(hfbpsi(0:nz,0:nr,1:2,ist)**2)
  psi(:,:,:,ist) = psi(:,:,:,ist)-hfbpsi(:,:,:,ist)
  IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
    WRITE(7,'(a,i4)') 'ortho-normality after: ist=',ist
    WRITE(7,'(100(1pg13.5))')  &
        (wfovlp(psi(0:nz,0:nr,1:2,jst),psi(0:nz,0:nr,1:2,ist)),jst=istmin,istmax)
  END IF
END DO



RETURN
END SUBROUTINE hfbstep

!-----damp_HFB
!

SUBROUTINE damp_hfb(psiact,v2act,istact)
!
!     Nucleon damping for HFB iteration.
!     The damping matrices are computed for each case anew.
!     Presently, the standard damping strategy is recycled.
!     This should be replaced later by direct inversion.
!
!     Input parameters are:
!      psiact    = s.p. wavefunction to be damped
!      v2act     = occupation of the actual state
!      istact    = number of actual state
!
USE paramr
USE dynpar
USE wavefr
USE levels
USE fields, ONLY: pairpt
USE init, ONLY: idanuc
!$ use omp_lib

REAL(8), INTENT(IN OUT) :: psiact(0:nz,0:nr,nsd)
REAL(8), INTENT(IN)     :: v2act
INTEGER, INTENT(IN)     :: istact

!REAL(8), PARAMETER :: xsmall=1D-99
LOGICAL, PARAMETER :: ttest=.FALSE.

REAL(8),ALLOCATABLE :: psip(:,:,:)
REAL(8) :: pptmax,e0sav,x0sav,uvact
INTEGER :: iq,is,mm

ALLOCATE(psip(0:nz,0:nr,nsd))
e0sav = e0inv
x0sav = x0inv
iq = icharg(istact)
uvact = SQRT(MAX(v2act-v2act*v2act,xsmall))

!     maximum gap on the grid
pptmax = maxval(ABS(pairpt(0:nz,0:nr,iq)))
IF(ttest) WRITE(7,'(a,i3,1pg13.5)') 'iq,pptmax=',iq,pptmax


!  step   = x0inv/e0inv
!IF(num_threads==1) THEN
IF(num_threads==-1) THEN
  e0inv = e0inv+pptmax*uvact/MAX(v2act,v2limh)*delmix
  x0inv = x0inv/MAX(v2act,v2limh)
  IF(ttest) WRITE(7,'(a,3(1pg13.5))')  &
    'IDANUC with v2act,e0inv,x0inv=',v2act,e0inv,x0inv
  CALL idanuc
  CALL damnuc(psiact,psip,istact)
ELSE
!!!!$  WRITE(6,*) 'OpenMP switch',OMP_GET_THREAD_NUM() !mythread
  CALL damnuc_HFB(psiact,psip,istact,e0inv+pptmax*uvact/MAX(v2act,v2limh)*delmix,x0inv/MAX(v2act,v2limh))
END IF
psiact = psip

!  set boundary zeroes if parity requires it.
DO is=1,2
  mm = mpsi(istact)+is-1
  IF (ABS(mm) > 0) psiact(0:nz,0,is)  =  0D0
  IF (nparz(is,istact) == 2) psiact(0,0:nr,is)  =  zero
END DO

e0inv = e0sav
x0inv = x0sav
IF(num_threads==-1) CALL idanuc

DEALLOCATE(psip)

RETURN
END SUBROUTINE damp_hfb


!---------------------------------------------------------------------------  
! DESCRIPTION: damnuc_HFB
!> @brief
!!Kinetic energy damping as in subroutinbe 'damnuc'. This version uses
!!dynamic damping arrays 'rdnuc_dyn' and 'zudnuc_dyn' to allow OpenMP
!in HFB calculations. 
!!Unlike 'damnuc', the routine embodies also the construction of the
!!damping arrays to accommodate changing damping operators. This part
!!is taken from the subroutine 'idanuc'.
!>
!> @param[in] psiin
!> REAL(8), array, takes the input wave function.
!> @param[out] psiout
!> REAL(8), array, returns the Hamiltonian applied to the input wave 
!!function.
!> @param[in] ist
!> INTEGER, absolute adress of the state.
!--------------------------------------------------------------------------- 

SUBROUTINE damnuc_HFB(psiin,psiout,ist,e0act,x0act)

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)  :: e0act,x0act
REAL(8), INTENT(IN)  :: psiin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT) :: psiout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)  :: ist

REAL(8),ALLOCATABLE :: zdnuc_dyn(:,:,:)        !<Damping operator for the nucleons in z-direction.
REAL(8),ALLOCATABLE :: rdnuc_dyn(:,:,:)        !<Damping operator for the nucleons in r-direction.


REAL(8) :: step,befac
INTEGER :: mm, is, i, impsi, ipi

ALLOCATE(zdnuc_dyn(0:nz,0:nz,4), rdnuc_dyn(0:nr,0:nr,0:kmpsi))

! first, the construction of the damping matrices

  befac  = h2mav/(two*e0act)
!     nucleon damping in z-direction
  zdnuc_dyn(0:nz,0:nz,:) = -befac * derz2f(0:nz,0:nz,:)
  DO i=0,nz
    zdnuc_dyn(i,i,:) = one + zdnuc_dyn(i,i,:)
  END DO
!     inversion 
  DO ipi=1,4    
    CALL matr_inv(zdnuc_dyn(0:nz,0:nz,ipi),nz+1,nz+1)
    IF (ipi == 2.OR.ipi == 4) THEN
      zdnuc_dyn(0,0:nz,ipi) = 0.0D0
      zdnuc_dyn(0:nz,0,ipi) = 0.0D0  ! addition PK
    END IF
  END DO

!     the nucleon damping in r-direction
  rdnuc_dyn(0:nr,0:nr,0:kmpsi) = -befac * derr2af(0:nr,0:nr,0:kmpsi,1)
  DO i=0,nr
    rdnuc_dyn(i,i,0:kmpsi) =  one + rdnuc_dyn(i,i,0:kmpsi)
  END DO

  mm = mpsi(ist)
!  IF(mm<0) STOP 'DAMNUC_HFB not yet prepared for negative MPSI'
  DO impsi=mm,mm+1
    CALL matr_inv(rdnuc_dyn(0:nr,0:nr,ABS(impsi)),nr+1,nr+1)
    IF (ABS(impsi) > 0) THEN
      rdnuc_dyn(0,0:nr,ABS(impsi)) = 0D0
      rdnuc_dyn(0:nr,0,ABS(impsi)) = 0D0 ! addition PK
    END IF
  END DO


!     the straightforward imaginary-time step in case of "e0act=0"

IF(e0act == 0D0) THEN
  psiout(0:nz,0:nr,1:2) = x0act*psiin(0:nz,0:nr,1:2)
ELSE
  
!     now the full kinetic energy damping
  CALL damp_nucleon(rdnuc_dyn(0:nr,0:nr,:),zdnuc_dyn(0:nz,0:nz,:),psiin(0:nz,0:nr,:),&
                    psiout(0:nz,0:nr,:),nparz(1,ist),mpsi(ist),x0act/e0act)

!       set boundary zeroes if parity requires it.
  DO is=1,2
    mm = mpsi(ist)+is-1
!    IF (mm > 0) THEN
    IF (ABS(mm) > 0) THEN
      psiout(0:nz,0,is)  =  0D0
    END IF
    IF (nparz(is,ist) == 2) THEN
      psiout(0,0:nr,is)  =  0D0
    END IF
  END DO
  
END IF

DEALLOCATE(zdnuc_dyn, rdnuc_dyn)

RETURN
END SUBROUTINE damnuc_HFB
