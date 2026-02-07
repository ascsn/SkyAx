!------------------------------------------------------------------------------
! MODULE: comp_block
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the blocking variables and routines
!------------------------------------------------------------------------------
MODULE blocking

USE paramr
USE levels


INTEGER, PARAMETER :: kblist=20   !<Maximum possible number of blocked states in each species.
INTEGER :: akblock(kblist,2)      !<Storage of 2*K values of blocked s.p. state.
INTEGER :: aktcblock(kblist,2)    !<Storage of 2*K values of conjugate to blocked s.p. state. 
INTEGER :: amblock(kblist,2)      !Storage of M values of blocked s.p. state. 
INTEGER :: aistblo(kblist,2)      !<Storage of number of blocked s.p. state in level sequence.
INTEGER :: aistcblo(kblist,2)     !<Storag of number of conjugate to blocked s.p. state in level sequence.
INTEGER :: aisoblo(kblist)        !<Storage of isospin of blocked s.p. state. To be set at compile time.


CONTAINS

!NEW ODD
!---------------------------------------------------------------------------  
! DESCRIPTION: comp_block
!> @brief
!! Computes the blocked state after some pre-iterations. Blocked is a state 
!! with \c kblock, \c mblock and \c isoblock as given in m_levels.f90. The 
!! routine finds out what number this state has in the normal ordering 
!! sequence of the code.
!> @details
!! Present version only for given state-number of blocked states. To be extended.
!
!--------------------------------------------------------------------------- 
SUBROUTINE comp_block()
USE dynpar
USE wavefr
USE grid
IMPLICIT NONE

!REAL(8) :: signup
INTEGER :: iq,iunit !,mbact,nfin,ninit

!IF(2>1) RETURN

DO iq=1,2
  
  IF(MOD(npart(iq),2).EQ.0 .AND. istblo(iq).NE.0) &
      STOP 'ISTBLO>0 not possible for even particle number'
  IF(istblo(iq).NE.0) THEN
    istcblo(iq) = itpartn(istblo(iq))
    mblock(iq) = mpsi(istblo(iq))
    kblock(iq) = 2*mblock(iq)+1
    vsqf(istblo(iq))=1D0
    vsqf(istcblo(iq))=0D0
  END IF

END DO

DO iunit=6,7
  WRITE(iunit,'(a,20i5)') ' blocking set: kblock,mblock,istblo,istcblo=',  &
      kblock,mblock,istblo,istcblo
END DO

RETURN
END SUBROUTINE comp_block
!END NEW ODD


END MODULE blocking
