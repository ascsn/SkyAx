!---------------------------------------------------------------------------  
! DESCRIPTION: init_Gaussfolding
!> @brief
!! Initialization of matrices 'gaufolz' and 'gaufolr' for Gaussian folding.
!> @details
!! Input are the matrices for 2. derivatives for which the eigenstates
!! on the grid are computed by diagonalization. The Gaussian folding
!! matrices are composed from eigenstates and eigenvalues.
!! Input and output are communicated through the header of the module.
!--------------------------------------------------------------------------- 
SUBROUTINE init_Gaussfolding()

USE grid

IMPLICIT NONE

REAL(8),ALLOCATABLE:: eigvec(:,:),eigval(:),work(:)
INTEGER,ALLOCATABLE:: iwork(:)

!REAL(8)            :: a(0:2*imaxdrvrng,0:2*imaxdrvrng)
!REAL(8)            :: b(0:2*imaxdrvrng)

INTEGER            :: iparmax, ibuf,in, ndimact
INTEGER            :: inx,k,nx,i,im,ix,j,jx
INTEGER            :: lwork,liwork,info
REAL(8)            :: buf,buf1,buf2


ALLOCATE(gaufolz(0:nz,0:nz,4),gaufolr(0:nr,0:nr,0:kmpsi,2))

!     ===============================================================
!     initialize folding matrix in z direction
!     ===============================================================
ndimact=nz
lwork=1+6*ndimact+2*ndimact**2
liwork=3+5*ndimact
ALLOCATE(eigvec(ndimact,ndimact),eigval(ndimact),work(lwork),iwork(liwork))

eigvec=0D0;eigval=0D0;work=0D0;iwork=0;j=0
DO in=1,4
  eigvec=derz2f(:,:,in)
  CALL DSYEVD('V','U',ndimact,eigvec,ndimact,eigval,work,lwork,iwork,liwork,info)
  eigval = EXP(eigval*(0.5D0*hfbwidth)**2)
  DO i=1,ndimact; DO j=1,ndimact
!   sum over "Bessel" eigenstates                                                                                                             
    gaufolz(i,j,in)=0D0
    DO k = 1, ndimact
      gaufolz(i,j,in) = gaufolz(i,j,in) + eigvec(i,k)*eigval(k)*eigvec (j,k)
    END DO
  END DO; END DO
END DO

DEALLOCATE(eigvec,eigval,work,iwork)
WRITE(6,*) ' Setup of Gaussian folding in z direction completed'

!     ===============================================================
!     initialize folding matrix in r direction
!     ===============================================================
ndimact=nr
lwork=1+6*ndimact+2*ndimact**2
liwork=3+5*ndimact
ALLOCATE(eigvec(ndimact,ndimact),eigval(ndimact),work(lwork),iwork(liwork))

DO in=1,2; DO im=0,kmpsi
  
  eigvec=0D0;eigval=0D0;work=0D0;iwork=0;j=0
  eigvec=derr2af(:,:,im,in)
  CALL DSYEVD('V','U',ndimact,eigvec,ndimact,eigval,work,lwork,iwork,liwork,info)
  eigval = EXP(eigval*(0.5D0*hfbwidth)**2)
  DO i=1,ndimact; DO j=1,ndimact
!   sum over "Bessel" eigenstates                                                                                                              
    gaufolr(i,j,im,in)=0D0
    DO k = 1, ndimact
      gaufolr(i,j,im,in) = gaufolr(i,j,im,in) + eigvec(i,k)*eigval(k)*eigvec (j,k)
    END DO
  END DO; END DO

END DO; END DO
DEALLOCATE(eigvec,eigval,work,iwork)

WRITE(6,*) ' Setup of Gaussian folding in r direction completed'

RETURN
END SUBROUTINE init_Gaussfolding


!---------------------------------------------------------------------------  
! DESCRIPTION: Gaussfolding_spin
!> @brief
!! Calculate the action of the Gaussian folding on a spinor.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] izpfld
!> INTEGER, takes the z-parity of upper component (1 for positive, 2 for negative,
!!1 and isym=0 for no good z-parity)
!> @param[in] imfld
!> INTEGER, takes the z-component of angular momentum m of the upper component.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE Gaussfolding_spin(fldin,fldout,izpfld,imfld,ibpfld,factor)

USE grid

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)        :: izpfld
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

REAL(8), ALLOCATABLE :: work(:,:)

INTEGER :: imact,izpact,ibpidx,is,izpidx

ALLOCATE(work(0:nz,0:nr))
DO is=1,nsd
  imact  = ABS(imfld+is-1)
  izpact = 3*(izpfld-1) - (2*izpfld-3)*is
  izpidx = izpact*isym  + (1-isym)
  ibpidx = izpidx + 2*ibpfld - 2
  
  CALL mulz(gaufolz(:,:,ibpidx),fldin(:,:,is),work,factor,0)
  CALL mulr(gaufolr(:,:,imact,ibpfld),work,fldout(:,:,is),factor,0)

END DO
DEALLOCATE(work)

RETURN
END SUBROUTINE Gaussfolding_spin
