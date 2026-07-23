!---------------------------------------------------------------------------  
! DESCRIPTION: mulz
!> @brief
!!Matrix-matrix multiplication for operations in z-direction. 
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix \c a.
!> @param[inout] x
!> REAL(8), array, takes the matrix \c x.
!> @param[inout] y
!> REAL(8), array, takes the matrix \c y.
!> @param[in] factor
!> REAL(8), takes the value \c factor.
!> @param[in] iaccu
!> INTEGER, 0 for no accumulation, 1 for accumulation of \c y.

!--------------------------------------------------------------------------- 

SUBROUTINE mulz(a,x,y,factor,iaccu)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nz,0:nz)
REAL(8), INTENT(IN)             :: x(0:nz,0:nr)
REAL(8), INTENT(INOUT)          :: y(0:nz,0:nr)
REAL(8), INTENT(IN)             :: factor
INTEGER, INTENT(IN)             :: iaccu

REAL(8) :: beta

!WRITE(6,*) 'new MULZ'
beta=one
IF (iaccu == 0) beta = zero

IF (ABS(factor) > smal) CALL DGEMM('N','N',nz+1,nr+1,nz+1,factor,a,nz+1,x,nz+1,beta,y,nz+1) 

RETURN
END SUBROUTINE mulz

!---------------------------------------------------------------------------  
! DESCRIPTION: mulr
!> @brief
!!Matrix-matrix multiplication for operations in r-direction. 
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix \c a.
!> @param[inout] x
!> REAL(8), array, takes the matrix \c x.
!> @param[inout] y
!> REAL(8), array, takes the matrix \c y.
!> @param[in] factor
!> REAL(8), takes the value \c factor.
!> @param[in] iaccu
!> INTEGER, 0 for no accumulation, 1 for accumulation of \c y.
!--------------------------------------------------------------------------- 

SUBROUTINE mulr(a,x,y,factor,iaccu)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nr,0:nr)
REAL(8), INTENT(IN)             :: x(0:nz,0:nr)
REAL(8), INTENT(INOUT)          :: y(0:nz,0:nr)
REAL(8), INTENT(IN)             :: factor
INTEGER, INTENT(IN)             :: iaccu

REAL(8) :: beta

beta=one
IF (iaccu == 0) beta = zero

IF (ABS(factor) > smal) CALL DGEMM('N','T',nz+1,nr+1,nr+1,factor,x,nz+1,a,nr+1,beta,y,nz+1) 

RETURN
END SUBROUTINE mulr

!---------------------------------------------------------------------------  
! DESCRIPTION: mulr0
!> @brief
!!Matrix-matrix multiplication for operations in r-direction (only r=0). 
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix \c a.
!> @param[inout] x
!> REAL(8), array, takes the matrix \c x.
!> @param[inout] y
!> REAL(8), array, takes the matrix \c y.
!> @param[in] factor
!> REAL(8), takes the value \c factor.
!> @param[in] iaccu
!> INTEGER, 0 for no accumulation, 1 for accumulation of \c y.
!--------------------------------------------------------------------------- 

SUBROUTINE mulr0(a,x,y,factor,iaccu,iband)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nr,0:nr)
REAL(8), INTENT(IN)             :: x(0:nz,0:nr)
REAL(8), INTENT(INOUT)          :: y(0:nz,0:nr)
REAL(8), INTENT(IN)             :: factor
INTEGER, INTENT(IN)             :: iaccu
INTEGER, INTENT(IN)             :: iband      ! bandwidth, presently only -1 accepted

REAL(8) :: beta

IF(iband.NE.-1) STOP 'MULR0: finite bandwidth not implemented'
beta=one
IF (iaccu == 0) beta = zero

IF (ABS(factor) > smal) CALL DGEMV('N',nz+1,nr+1,factor,x,nz+1,a,nz+1,beta,y,1)

RETURN
END SUBROUTINE mulr0
