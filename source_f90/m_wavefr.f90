!------------------------------------------------------------------------------
! MODULE: wavefr
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the main variable for wave functions and other variables 
!!for auxiliary storaged of s.p. wave functions. Furthermore it contains routines 
!!for overlaps and projections of s.p.states and time-reversed s.p. states.
!------------------------------------------------------------------------------

MODULE wavefr

USE paramr
USE grid, ONLY: nz,nr,int_weight,init_int_weight,f2int

IMPLICIT NONE
 
REAL(8),POINTER :: psi(:,:,:,:)                  !< The main storage for all s.p. states.
REAL(8),POINTER :: psifold(:,:,:,:)              !< Storage for all folded s.p. states.

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_wf
!> @brief
!!Allocates the storage for wave functions.
!>
!> @param[in] npsi
!> INTEGER, takes the total number of wave function.
!--------------------------------------------------------------------------- 
SUBROUTINE alloc_wf(npsi)
  USE dynpar, ONLY: thfbfold
  IMPLICIT NONE
  INTEGER,INTENT(IN) ::npsi
  ALLOCATE(psi(0:nz,0:nr,nsd,npsi))
  psi=0.0d0
  IF(thfbfold) THEN
    ALLOCATE(psifold(0:nz,0:nr,nsd,npsi))
    psifold=0D0
  END IF
END SUBROUTINE alloc_wf

!---------------------------------------------------------------------------  
! DESCRIPTION: dealloc_wf
!> @brief
!!Deallocates the storage for wave functions.
!--------------------------------------------------------------------------- 
SUBROUTINE dealloc_wf()
  IMPLICIT NONE
  DEALLOCATE(psi)
END SUBROUTINE dealloc_wf

!---------------------------------------------------------------------------  
! DESCRIPTION: wfovlp
!> @brief
!!Computes the overlap between two wave functions.
!>
!> @param[in] psia
!> REAL(8), array, takes the first wave function.
!> @param[in] psib
!> REAL(8), array, takes the second wave function.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION wfovlp(psia,psib)

  IMPLICIT NONE
  
  REAL(8), INTENT(IN)         :: psia(0:nz,0:nr,nsd)
  REAL(8), INTENT(IN)         :: psib(0:nz,0:nr,nsd)

  wfovlp=f2int(psia(0:nz,0:nr,1),psib(0:nz,0:nr,1))+f2int(psia(0:nz,0:nr,2),psib(0:nz,0:nr,2))

  RETURN
  
END FUNCTION wfovlp

!---------------------------------------------------------------------------  
! DESCRIPTION: conjg_wfovlp
!> @brief
!!Computes the overlap between a wave function and a the time-conjugated
!!other wavefauntion.
!>
!> @param[in] psia
!> REAL(8), array, takes the first wave function.
!> @param[in] psib
!> REAL(8), array, takes the second, time-conjugated  wave function.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION conjg_wfovlp(psia,psib)

  IMPLICIT NONE
  
  REAL(8), INTENT(IN)         :: psia(0:nz,0:nr,nsd)
  REAL(8), INTENT(IN)         :: psib(0:nz,0:nr,nsd)

  conjg_wfovlp=-f2int(psia(0:nz,0:nr,1),psib(0:nz,0:nr,2))+f2int(psia(0:nz,0:nr,2),psib(0:nz,0:nr,1))

  RETURN
  
END FUNCTION conjg_wfovlp

!---------------------------------------------------------------------------  
! DESCRIPTION: time_reversal
!> @brief
!!Returns the time-reversed state of \c psia
!>
!> @param[in] psia
!> REAL(8), array, takes the first wave function.
!--------------------------------------------------------------------------- 
PURE FUNCTION time_reversal(psia)

  IMPLICIT NONE

  REAL(8), INTENT(IN)         :: psia(0:nz,0:nr,nsd)
  REAL(8)                     :: time_reversal(0:nz,0:nr,nsd)

  time_reversal(:,:,1)=-psia(:,:,2)
  time_reversal(:,:,2)=psia(:,:,1)
  RETURN
  
END FUNCTION time_reversal


!---------------------------------------------------------------------------  
! DESCRIPTION: wfproj
!> @brief
!!projects components of spinor \c psia out of \c psib:
!!\f[  |b\rangle = |b\rangle - fac \times |a\rangle  \f]
!>
!> @param[in] psia
!> REAL(8), array, takes the first wave function.
!> @param[inout] psib
!> REAL(8), array, takes the second wave function.
!> @param[in] fac
!> REAL(8), takes the factor of wave function \c psia.
!--------------------------------------------------------------------------- 

SUBROUTINE wfproj(psia,psib,fac)

IMPLICIT NONE

REAL(8), INTENT(IN)       :: psia(0:nz,0:nr,nsd)
REAL(8), INTENT(INOUT)    :: psib(0:nz,0:nr,nsd)
REAL(8), INTENT(IN)       :: fac

  psib = psib - fac * psia

RETURN
END SUBROUTINE wfproj


!---------------------------------------------------------------------------  
! DESCRIPTION: wfovpc
!> @brief
!!Computes the norm of a wave functions, without spin.
!>
!> @param[in] psia
!> REAL(8), array, takes the wave function.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION wfovpc(psia)

IMPLICIT NONE

REAL(8), INTENT(IN)            :: psia(0:nz,0:nr)

wfovpc = wfovoc(psia,psia)

RETURN
END FUNCTION wfovpc

!---------------------------------------------------------------------------  
! DESCRIPTION: wfovoc
!> @brief
!!Computes the overlap between two wave functions, without spin.
!>
!> @param[in] psia
!> REAL(8), array, takes the first wave function.
!> @param[in] psib
!> REAL(8), array, takes the second wave function.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION wfovoc(psia,psib)

  IMPLICIT NONE

  REAL(8), INTENT(IN)         :: psia(0:nz,0:nr)
  REAL(8), INTENT(IN)         :: psib(0:nz,0:nr) 

  wfovoc=f2int(psia(0:nz,0:nr),psib(0:nz,0:nr))
RETURN
END FUNCTION wfovoc

END MODULE wavefr
