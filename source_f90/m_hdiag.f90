!------------------------------------------------------------------------------
! MODULE: hdiag
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routines to compute Hamiltonian matrices and to 
!!diagonalize them. The matrices are always assumed to be in blocks for groups 
!!with equal m quantum number.
!------------------------------------------------------------------------------
MODULE hdiag
IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: collmatr
!> @brief
!!Computes one line of a Hamiltonian matrix submatrix with equal quantum number m.
!> @details
!!The subroutine takes the Hamiltonian applied on a wave function with index \c jst and 
!!calculates the overlap with wave functions with the same quantum number m. If the 
!!parity in z-direction is different, the overlap is set to 0.
!>
!> @param[out] hmatr
!> REAL(8), array, returns the Hamiltonian matrix.
!> @param[in] psi
!> REAL(8), array, takes wave functions in current m-block.
!> @param[in] hpsi
!> REAL(8), array, takes Hamiltonian applied to the wave function with index jst.
!> @param[in] nlo
!> INTEGER, array, takes index of first wave function in the m-block.
!> @param[in] jst
!> INTEGER, takes index of wave function for which overlaps are calculated.
!> @param[in] ndimact
!> INTEGER, takes size of current m-block.
!--------------------------------------------------------------------------- 
SUBROUTINE collmatr(hmatr,psi,hpsi,nlo,nhi,jst,ndimact)

USE wavefr  , ONLY : wfovlp
USE levels  , ONLY : nparz,spenrg
USE paramr  , ONLY : nsd
USE grid    , ONLY : nr,nz

IMPLICIT NONE

INTEGER, INTENT(IN)      :: ndimact
REAL(8), INTENT(OUT)     :: hmatr(ndimact,ndimact)
REAL(8), INTENT(IN)      :: hpsi(0:nz,0:nr,nsd)
REAL(8), INTENT(IN)      :: psi(0:nz,0:nr,nsd,ndimact)
INTEGER, INTENT(IN)      :: nlo,nhi
INTEGER, INTENT(IN)      :: jst

INTEGER                  :: nact,kst,nhst
!----------------------------------------------------------------------

DO kst=nlo,nhi
  nact = jst-nlo+1

  IF (nparz(1,jst) == nparz(1,kst)) THEN
    hmatr(kst,jst) = wfovlp(psi(0:nz,0:nr,1:2,kst),hpsi(0:nz,0:nr,1:2))
  ELSE
    hmatr(kst,jst) = 0.0D0
  END IF
  IF(kst == jst) hmatr(kst,jst)= hmatr(kst,jst)+spenrg(jst)

END DO

RETURN
END SUBROUTINE collmatr

!---------------------------------------------------------------------------  
! DESCRIPTION: diag_wf
!> @brief
!!Diagonalizes 'hmatr' and transforms s.p. wavefunctions.
!> @details
!!The subroutine takes the Hamiltonian applied on a wave function with index jst and 
!!calculates the overlap with wave functions with the same quantum number m. If the 
!!parity in z-direction is different, the overlap is set to 0.
!>
!> @param[inout] hmatr
!> REAL(8), array, takes the Hamiltonian matrix.
!> @param[in] nlo
!> INTEGER, array, takes index of first wave function in the m-block.
!> @param[in] ndimact
!> INTEGER, takes size of current m-block.
!--------------------------------------------------------------------------- 

SUBROUTINE diag_wf(hmatr,nlo,ndimact)

USE paramr  , ONLY : one,zero,nsd
USE wavefr  , ONLY : psi
USE grid    , ONLY : nr,nz

IMPLICIT NONE

INTEGER, INTENT(IN)               :: nlo,ndimact
REAL(8), INTENT(IN OUT)           :: hmatr(ndimact,ndimact)

LOGICAL,PARAMETER                 :: tsimpdiag=.true.

REAL(8),ALLOCATABLE               :: eigvec(:,:),eigval(:),work(:),temp(:,:),psitemp(:,:,:,:)
INTEGER,ALLOCATABLE               :: iwork(:),j(:)
INTEGER                           :: griddim,lwork,liwork,info,i,k

!----------------------------------------------------------------------

IF(ndimact == 1) RETURN                   ! nothing to do

lwork=1+6*ndimact+2*ndimact**2
liwork=3+5*ndimact
ALLOCATE(eigvec(ndimact,ndimact),eigval(ndimact),work(lwork),iwork(liwork),&
         temp(ndimact,ndimact),j(1))
eigvec=0D0;eigval=0D0;work=0D0;iwork=0;
temp=0D0;j=0
IF(tsimpdiag) THEN
  CALL diagsimp(hmatr,eigvec,ndimact)
ELSE
  temp=hmatr
  eigvec=temp
  CALL DSYEVD('V','U',ndimact,eigvec,ndimact,eigval,work,lwork,iwork,liwork,info)
  DO k=1,ndimact
  DO i=k,ndimact
    j=MAXLOC(abs(eigvec(:,i)))
    temp(:,1)=eigvec(:,i)
    eigvec(:,i)=eigvec(:,j(1))
    eigvec(:,j(1))=temp(:,1)
  END DO
  END DO
END IF
!     transform wavefunctions
griddim=(nz+1)*(nr+1)*nsd

ALLOCATE(psitemp(0:nz,0:nr,1:nsd,nlo:nlo+ndimact))
psitemp=0D0

CALL DGEMM('N','N',griddim,ndimact,ndimact,one,psi(0:nz,0:nr,1:nsd,nlo:nlo+ndimact),&
           griddim,eigvec,ndimact,zero,psitemp(0:nz,0:nr,1:nsd,nlo:nlo+ndimact),griddim)

psi(0:nz,0:nr,1:nsd,nlo:nlo+ndimact-1) = psitemp(0:nz,0:nr,1:nsd,nlo:nlo+ndimact-1)

DEALLOCATE(eigvec,eigval,work,iwork,temp,j,psitemp)
!count=count+1

RETURN
END SUBROUTINE diag_wf

!---------------------------------------------------------------------------  
! DESCRIPTION: diagsimp
!> @brief
!!Employs a simple gradient step method for approximate diagonalization of 
!!the Hamiltonian matrix.
!>
!> @param[inout] hmatr
!> REAL(8), array, takes the Hamiltonian matrix. It is destroyed and used as workspace.
!> @param[out] eigvec
!> REAL(8), array, returns the eigenvectors of the Hamiltonian Matrix.
!> @param[in] ndimact
!> INTEGER, takes size of current m-block.
!--------------------------------------------------------------------------- 

SUBROUTINE diagsimp(hmatr,eigvec,ndimact)

USE paramr, ONLY :smal
IMPLICIT NONE

INTEGER, INTENT(IN)      :: ndimact
REAL(8), INTENT(IN OUT)  :: hmatr(ndimact,ndimact)
REAL(8), INTENT(OUT)     :: eigvec(ndimact,ndimact)

INTEGER, PARAMETER       :: itdiagmax=200
REAL(8), PARAMETER       :: diagstep=0.01D0
REAL(8), PARAMETER       :: diagprecis=1D-6
LOGICAL, PARAMETER       :: tprint=.false.

REAL(8),ALLOCATABLE      :: workv1(:,:),workv2(:,:)
REAL(8)                  :: variance,sumeig,hmin,sumeig2
INTEGER                  :: na,nb,itdiag

!----------------------------------------------------------------------

IF(ndimact == 1) THEN
  eigvec(1,1) = 1D0
  RETURN
END IF

ALLOCATE(workv1(ndimact,ndimact),workv2(ndimact,ndimact))
workv1=0D0
workv2=0D0

IF(tprint) THEN
  WRITE(6,'(a)') 'H matrix:'
  DO na=1,ndimact
    WRITE(6,'(10(1pg13.5))') (hmatr(na,nb),nb=1,ndimact)
  END DO
END IF
!***********************************************************
!     prepare stepping matrix and eigenvectors
!***********************************************************
!search for minimum diagonal entry
hmin = 1D99
DO na=1,ndimact
  hmin = MIN(hmatr(na,na),hmin)
END DO

hmatr = -diagstep*hmatr
eigvec=0.0d0
DO na=1,ndimact
  hmatr(na,na)=1.0D0+hmatr(na,na)+diagstep*hmin 
  eigvec(na,na)=1.0D0
END DO
!     iterate
iterations: DO itdiag=1,itdiagmax  
!       analysis of variance
  workv1 = MATMUL(hmatr(:,:),eigvec(:,:))      
  workv2 = MATMUL(hmatr(:,:),workv1(:,:))
  DO na=1,ndimact
    sumeig = sum(eigvec(:,na)*workv1(:,na))
    sumeig2 = sum(eigvec(:,na)*workv2(:,na))
    variance = sumeig2-sumeig**2
  END DO
  variance = SQRT(MAX(ABS(variance),smal))
  IF(tprint) WRITE(6,'(a,2i5,1pg12.4)')  &
      'diag: block,iter,variance:',ndimact,itdiag, variance
  IF(variance < diagprecis) EXIT iterations
!       stepping
  eigvec = workv1
  CALL orthmatr(eigvec,ndimact)
END DO iterations

DEALLOCATE(workv1,workv2)
RETURN
END SUBROUTINE diagsimp

!---------------------------------------------------------------------------  
! DESCRIPTION: orthmatr
!> @brief
!!Employs the Gram-Schmidt algorithm to orthonormalize eigenvectors.
!>
!> @param[inout] eigvec
!> REAL(8), array, takes the non-orthonormal eigenvectors and returns the 
!! orthonormal ones.
!> @param[in] ndimact
!> INTEGER, takes size of current m-block.
!--------------------------------------------------------------------------- 

SUBROUTINE orthmatr(eigvec,ndimact)

IMPLICIT NONE

INTEGER,INTENT(IN)             :: ndimact
REAL(8),INTENT(IN OUT)         :: eigvec(ndimact,ndimact)

LOGICAL, PARAMETER             :: ttest=.false.

INTEGER                        :: na,nb
REAL(8)                        :: ov
!----------------------------------------------------------------------

DO na=1,ndimact
  DO nb=1,na-1
    ov = sum(eigvec(:,nb)*eigvec(:,na))
    eigvec(:,na) = eigvec(:,na)-ov*eigvec(:,nb)
  END DO
  ov = 1.0D0/SQRT( sum(eigvec(:,na)**2) )
  eigvec(:,na) = ov*eigvec(:,na)
END DO

!     test
IF(.NOT.ttest) RETURN

WRITE(6,'(a)') ' test orthonormalization:'
DO na=1,ndimact
  DO nb=1,ndimact
    ov = sum(eigvec(:,nb)*eigvec(:,na))
    WRITE(6,'(2i5,1pg13.5)') na,nb,ov
  END DO
END DO

RETURN
END SUBROUTINE orthmatr

END MODULE hdiag
