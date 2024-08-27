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


ALLOCATE(gaufolz(0:nz,0:nz),gaufolr(0:nr,0:nr,0:kmpsi))

!     ===============================================================
!     initialize folding matrix in z direction
!     ===============================================================
ndimact=nz
lwork=1+6*ndimact+2*ndimact**2
liwork=3+5*ndimact
ALLOCATE(eigvec(ndimact,ndimact),eigval(ndimact),work(lwork),iwork(liwork))

eigvec=0D0;eigval=0D0;work=0D0;iwork=0;j=0
eigvec=derz2f(:,:,1)
CALL DSYEVD('V','U',ndimact,eigvec,ndimact,eigval,work,lwork,iwork,liwork,info)
eigval = EXP(-eigval*(0.5D0*foldwid)**2)
DO i=1,ndimact; DO j=1,ndimact
!   sum over "Bessel" eigenstates                                                                                                             
  gaufolz(i,j)=0D0
  DO k = 2, ndimact
    gaufolz(i,j) = gaufolz(i,j) + eigvec(i,k)*eigval(k)*eigvec (j,k)
  END DO
END DO; END DO

DEALLOCATE(eigvec,eigval,work,iwork)
WRITE(6,*) ' Setup of Gaussian folding in z direction completed'

!     ===============================================================
!     initialize folding matrix in r direction
!     ===============================================================
ndimact=nr
lwork=1+6*ndimact+2*ndimact**2
liwork=3+5*ndimact
ALLOCATE(eigvec(ndimact,ndimact),eigval(ndimact),work(lwork),iwork(liwork))

DO im=0,kmpsi
  
  eigvec=0D0;eigval=0D0;work=0D0;iwork=0;j=0
  eigvec=derr2af(:,:,im,1)
  CALL DSYEVD('V','U',ndimact,eigvec,ndimact,eigval,work,lwork,iwork,liwork,info)
  eigval = EXP(-eigval*(0.5D0*foldwid)**2)
  DO i=1,ndimact; DO j=1,ndimact
!   sum over "Bessel" eigenstates                                                                                                              
    gaufolr(i,j,im)=0D0
    DO k = 2, ndimact
      gaufolr(i,j,im) = gaufolr(i,j,im) + eigvec(i,k)*eigval(k)*eigvec (j,k)
    END DO
  END DO; END DO

END DO
DEALLOCATE(eigvec,eigval,work,iwork)

WRITE(6,*) ' Setup of Gaussian folding in r direction completed'

RETURN
END SUBROUTINE init_Gaussfolding
