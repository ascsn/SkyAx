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

REAL(8) :: v2a,uva
INTEGER :: mist,npzist,jst,ist,is,iq


! set isospin (protons or neutrons)
IF(istmin <= nst(1)) THEN
  iq = 1
ELSE
  iq = 2
END IF
symcond = 0D0
!     compute and store action of HFB Hamiltonian
DO jst=istmin,istmax
  CALL hpsi(psi(0:nz,0:nr,1:2,jst),hfbpsi(0:nz,0:nr,1:2,jst),jst)
  v2a = vsqf(jst)
  uva = SQRT(MAX(v2a-v2a*v2a,xsmall))
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
    DO is=1,2
      hfbpsi(0:nz,0:nr,is,jst) = hfbpsi(0:nz,0:nr,is,jst)*v2a  &
          +pairpt(0:nz,0:nr,iq)*psi(0:nz,0:nr,is,jst)*uva
    END DO
  END IF
  IF(tqq) CALL accumulate_variances_1(psi(0:nz,0:nr,1:2,jst),psin,iq,jst,istmin,istmax,&
                            spevar(jst),spenrg(jst),&     
                            qq0,qq1,qq2,qq3,qq4,qh0,qh1,qh2,qh3,&  
                            qh4,qhdig,devmf2,devgap,devh2,devoc2)
END DO


!     compute matrix of Lagrange parameters
ALLOCATE(hlagr(istmin:istmax,istmin:istmax))
ALLOCATE(symatrix(istmin:istmax,istmin:istmax))
hlagr = 0D0
symatrix = 0D0
DO ist=istmin,istmax
  mist = mpsi(ist)
  npzist = nparz(1,ist)
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
  END DO
END DO

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


!     kinetic energy damping
DO ist=istmin,istmax
  IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
    WRITE(7,'(a,i4,3(1pg13.5))') 'HFB overlaps before:',  &
        ist,wfovlp(psi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(hfbpsi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(psi(0:nz,0:nr,1:2,ist),psi(0:nz,0:nr,1:2,ist))
  END IF
  CALL damp_hfb(hfbpsi(:,:,:,ist),vsqf(ist),ist)
  IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
    WRITE(7,'(a,i4,3(1pg13.5))') 'HFB overlaps after:',  &
        ist,wfovlp(psi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(hfbpsi(0:nz,0:nr,1:2,ist),hfbpsi(0:nz,0:nr,1:2,ist)),  &
        wfovlp(psi(0:nz,0:nr,1:2,ist),psi(0:nz,0:nr,1:2,ist))
  END IF
  IF(tqq) CALL accumulate_variances_2(ist,istmin,istmax,hfbpsi(:,:,:,ist),qdh,hdh,hddh,hdhdh)  
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

REAL(8), INTENT(IN OUT) :: psiact(0:nz,0:nr,nsd)
REAL(8), INTENT(IN)     :: v2act
INTEGER, INTENT(IN)     :: istact

REAL(8), PARAMETER :: v2limh=1D-6,delmix=0.8D0  ! flexible damping parameters
!REAL(8), PARAMETER :: xsmall=1D-99
LOGICAL, PARAMETER :: ttest=.FALSE.

REAL(8) :: pptmax,e0sav,x0sav,uvact
INTEGER :: iq,is,mm

e0sav = e0inv
x0sav = x0inv
iq = icharg(istact)
uvact = SQRT(MAX(v2act-v2act*v2act,xsmall))

!     maximum gap on the grid
pptmax = maxval(ABS(pairpt(0:nz,0:nr,iq)))
IF(ttest) WRITE(7,'(a,i3,1pg13.5)') 'iq,pptmax=',iq,pptmax

e0inv = e0inv+pptmax*uvact/MAX(v2act,v2limh)*delmix
x0inv = x0inv/MAX(v2act,v2limh)
IF(ttest) WRITE(7,'(a,3(1pg13.5))')  &
    'IDANUC with v2act,e0inv,x0inv=',v2act,e0inv,x0inv
CALL idanuc

!  step   = x0inv/e0inv
CALL damnuc(psiact,psip,istact)
psiact = psip

!  set boundary zeroes if parity requires it.
DO is=1,2
  mm = mpsi(istact)+is-1
  IF (ABS(mm) > 0) psiact(0:nz,0,is)  =  0D0
  IF (nparz(is,istact) == 2) psiact(0,0:nr,is)  =  zero
END DO

e0inv = e0sav
x0inv = x0sav
CALL idanuc


RETURN
END SUBROUTINE damp_hfb
