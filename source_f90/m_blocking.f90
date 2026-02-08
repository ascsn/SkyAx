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
USE dynpar
USE grid
USE constr, ONLY: tfresh

INTEGER, PARAMETER :: kblist=20   !<Maximum possible number of blocked states in each species.
INTEGER :: akblock(kblist,2)      !<Storage of 2*K values of blocked s.p. state.
INTEGER :: aktcblock(kblist,2)    !<Storage of 2*K values of conjugate to blocked s.p. state. 
INTEGER :: amblock(kblist,2)      !Storage of M values of blocked s.p. state. 
INTEGER :: aistblo(kblist,2)      !<Storage of number of blocked s.p. state in level sequence.
INTEGER :: aistcblo(kblist,2)     !<Storag of number of conjugate to blocked s.p. state in level sequence.
INTEGER :: aisoblo(kblist,2)      !<Storage of isospin of blocked s.p. state. To be set at compile time.
REAL(8) :: edifblo(kblist,2)      !<Storage for quasi-particle energies of blocked states
REAL(8) :: updo(kblist,2)         !<Storage for difference of upper to lower component
! The basic blocking fields are defined and allocated in module 'levels'

INTEGER :: maxblist(2)
LOGICAL :: blnprot,blnneut        !<Switches to indicate odd nuclein number
REAL(8) :: absdiff=2D0            !<Maximum s.p. excitation energy for blocking configurations

LOGICAL, PARAMETER :: ttrack=.TRUE.  !<Switch for test output


CONTAINS


!---------------------------------------------------------------------------  
! DESCRIPTION: scan_block
!> @brief
!! Determines the s.p. states closest to the Fermi energy and finds out what 
!! number this state has in the normal ordering sequence of the code. Then
!! it prepares the appropriate blocking parameters \c kblock, \c mblock and \c isoblock.
!
!--------------------------------------------------------------------------- 
SUBROUTINE scan_block()
!USE wavefr
IMPLICIT NONE

IF(blnprot) CALL search_block(1)  
IF(blnneut) CALL search_block(2) 





END SUBROUTINE scan_block


!---------------------------------------------------------------------------  
! DESCRIPTION: block1
!> @brief
!! Computes the SHF solution for one blocked nucleon configuration.
SUBROUTINE block1(iq)

IMPLICIT NONE

INTEGER, INTENT(IN) :: iq

INTEGER :: blist

DO  blist = 1,maxblist(iq)
  
  mblock(iq) = amblock(blist,iq)
  kblock(iq) = akblock(blist,iq)
  isoblo = aisoblo(blist,iq)
!  ediffblock(iq) = edifblo(blist,iq)
   
 
  
  istblo(iq) = aistblo(blist,iq)
  istcblo(iq) = itpartn(istblo(iq)) ! time conj. partner
  
  
!  iqtest = MOD(iq,2)+1              ! blocking only neutrons or protons
!! necessary for odd-odd nuclei
!  kblock(iqtest) = 0
  
  
  WRITE(6,'(a,4(i4),1pg13.5)')'Search the minimal energy:  &
      mblock kblock isoblock istblo',  &
      mblock(iq),kblock(iq),isoblo,istblo(iq),edifblo(blist,iq)
  WRITE(7,'(a,5i4,1pg13.5)')'Search the minimal energy:  &
      mblock kblock isoblock istblo',  &
      mblock(iq),kblock(iq),isoblo,istblo(iq),istcblo(iq),edifblo(blist,iq)
  CALL flush(7)

  IF(tfresh) THEN
    CALL ssetup                    
    CALL inipai(npart(1)+npart(2))
  END IF

  CALL iteration_loop()
  
  
!  ang = 0
!  CALL calangularmomentum(istblo(iq),ang)
!  blenergy(blist,iq) = energy       ! binding energy for the particular blocking state
!  iparity = npar(istblo(iq))*nrpar(istblo(iq))
!  WRITE(17,'(a)') 'SEARCH_ENERGY: blist,m,k,iso,maxb,energ,ist,par,ang'
!  WRITE(17,'(5i4,1pg13.5,2i4,2(1pg13.5))')  &
!      blist,mblock(iq),kblock(iq),isoblo,maxblist(iq)  &
!      ,blenergy(blist,iq),istblo, iparity,ang
  
  
END DO

END SUBROUTINE block1

!---------------------------------------------------------------------------  
! DESCRIPTION: block2
!> @brief
!! Computes the SHF solution for one proton and one neutron blocked.
SUBROUTINE block2()

IMPLICIT NONE

INTEGER :: iq,maxblist_total,ip,ip1,ip2,blist1,blist2,ipsav(2)
REAL(8) :: energmin,energmax
REAL(8), ALLOCATABLE ::  total_energy(:),ipoint(:,:)


maxblist_total=maxblist(1)*maxblist(2)
allocate(total_energy(maxblist_total),ipoint(1:2,maxblist_total))

ipoint = 0
total_energy =0D0
ip=0
DO blist1 = 1,maxblist(1)
  DO blist2 = 1,maxblist(2)
    ip=ip+1
    total_energy(ip)=edifblo(blist1,1)+edifblo(blist2,2)
    ipoint(1,ip)=blist1
    ipoint(2,ip)=blist2
  END DO
END DO
IF(ttrack) THEN
  WRITE(7,'(a,4i4)') 'nr. of states:',maxblist(:),maxblist_total,ip
  WRITE(7,'(a)') 'energy before sorting of 2qp states:'
  WRITE(7,'(400(2i5,3(1pg13.5)/))') (ipoint(:,ip1),total_energy(ip1),  &
    edifblo(ipoint(1,ip1),1),edifblo(ipoint(2,ip1),2), ip1=1,maxblist_total)
END IF

!     sort on total energy
DO ip1=1,maxblist_total
  energmin = total_energy(ip1)
  DO ip2=ip1+1,maxblist_total
    IF(total_energy(ip2)<energmin) THEN
      total_energy(ip1)=total_energy(ip2)
      total_energy(ip2)=energmin
      energmin = total_energy(ip1)
      ipsav = ipoint(:,ip1)
      ipoint(:,ip1) = ipoint(:,ip2)
      ipoint(:,ip2) = ipsav
    END IF
  END DO
END DO
IF(ttrack) THEN
  WRITE(7,'(a)') 'final energy sorting of 2qp states:'
  WRITE(7,'(2i5,1pg13.5)')  &
    (ipoint(:,ip1),total_energy(ip1),ip1=1,maxblist_total)
END IF
energmin = minval(total_energy(:))
energmax = maxval(total_energy(:))
energmax = energmin+(energmax-energmin)/2D0

!     now scan the states
DO ip1=1,maxblist_total
  blist1 = ipoint(1,ip1)
  blist2 = ipoint(2,ip1)
  IF(total_energy(ip1)>energmax) EXIT
  
  istblo = 0
  mblock(1) = amblock(blist1,1)
  kblock(1) = akblock(blist1,1)
  mblock(2) = amblock(blist2,2)
  kblock(2) = akblock(blist2,2)
  istblo(1) = aistblo(blist1,1)
  istcblo(1) = itpartn(istblo(1)) ! time conj. partner
  istblo(2) = aistblo(blist2,2)
  istcblo(2) = itpartn(istblo(2)) ! time conj. partner
!  ediffblock(1) = edifblo(blist1,1)
!  ediffblock(2) = edifblo(blist2,2)
  
  
!        iqtest = mod(iq,2)+1              ! blocking only neutrons or protons
! necessary for odd-odd nuclei
!        kblock(iqtest) = 0
  
  IF(ttrack) THEN
    WRITE(6,'(a,8(i4))')'Search the minimal energy:  &
        mblock kblock istblo',mblock,kblock,istblo
    WRITE(7,'(a,8(i4))')'Search the minimal energy:  &
        mblock kblock istblo',mblock,kblock,istblo
    WRITE(*,*) 'before iteration loop in BLOCK2'
    WRITE(7,*) 'before iteration loop in BLOCK2'
  END IF
  
  CALL flush(7)

  IF(tfresh) THEN
    CALL ssetup                    
    CALL inipai(npart(1)+npart(2))
  END IF

  CALL iteration_loop()
  
  
!  ang = 0
!  CALL calangularmomentum(istblo(1),angularmom(1))
!  CALL calangularmomentum(istblo(2),angularmom(2))
!  blenergy(blist1,1) = energy       ! binding energy for the particular blocking state
!  iparity = npar(istblo(iq))*nrpar(istblo(iq))
!  WRITE(17,'(a)') 'SEARCH_ENERGY2: blist,m,k,iso,maxb,energ,ist,par,ang,rmsc'
!  WRITE(17,'(6i4,1pg13.5,2i4,3(1pg13.5))')  &
!      blist1,blist2,mblock,kblock,energy,istblo,angularmom,rmsc
!  CALL flush(17)
  
END DO




!mblock(iq) = amblock(blist,iq)
!kblock(iq) = akblock(blist,iq)
!isoblo = aisoblo(blist,iq)
!
!istblo(iq) = aistblo(blist,iq) !test
!istcblo (iq) = itpartn(istblo(iq))
!
!IF(updo(blist,iq) > 0)THEN
!  spinblo(iq) = 1
!ELSE
!  spinblo(iq) = 0
!END IF
!
!
!testblock(iq)=kblock(iq) ! for odd-odd nuclei
!testblock2(iq)=mblock(iq)!

DEALLOCATE(total_energy,ipoint)

RETURN

END SUBROUTINE block2


!---------------------------------------------------------------------------  
! DESCRIPTION: search_block
!> @brief
!! Searches for the s.p. configurations with the smallest distance to the
!! Fermi energy, stores them in book-keeping arrays, and sort on excitaton
!! energy.
SUBROUTINE search_block(iq)

IMPLICIT NONE

INTEGER, INTENT(IN)  :: iq   ! isospin for which levels are scanned

INTEGER :: blist,ine,ind,isav,ib2,ist,nstm
REAL(8) :: equasi              ! quasi-particle energy
REAL(8) :: emin,sav





IF(ttrack) WRITE(17,'(a,i2)') 'SEARCH_BLOCK: iq=',iq

blist   = 0
nstm=nst(iq)
DO ist=1,nstm
  ine    = iokk(iq)+ist-1
  ind    = iorder(ine)+iokk(iq)-1

  equasi=SQRT((spenrg(ind)-fer(iq))**2+deltaf(ind)**2)  
  IF(equasi < absdiff) THEN  

    ! exclude a couple of unwanted cases
    IF(( ((blnprot.AND.blnneut).AND.iq==1)  &
        .OR.(.NOT.(blnprot.AND.blnneut))   ) .AND. mpsi(ind)<0 ) CYCLE
!    IF(ABS(mpsi(ind)+0.5D0-angular(ind))>0.3D0) CYCLE    ! select maximum K orientation --> postpone

    blist = blist +1                              ! counter for the blocking states
    amblock(blist,iq) = mpsi(ind)    
    akblock(blist,iq) = 2*mpsi(ind)+1    !   2*K to store it as integer
    aisoblo(blist,iq) = icharg(ind)
    edifblo(blist,iq) = equasi
    updo(blist,iq)    = (cmpn(1,ind)-cmpn(2,ind)) ! to characterize the state in subroutine comp_block
    aistblo(blist,iq) = ind
    
    IF(ttrack) THEN
      WRITE(7,'(a,i4)') 'success: blist=',blist
      WRITE(7,'(a)') 'ind,charge,m,3*parity,v2,comps,spenerg'
      WRITE(7,'(i5,i4,f6.1,3i5,f7.4,2f5.2,f9.3)')  &
          ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind), nparz(1,ind),  &
          vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind)
      WRITE(7,'(a)') 'amblock,akblock,aisoblock,diff,updo:'
      WRITE(7,'(5(1pg13.5))') amblock(blist,iq),akblock(blist,iq),  &
        aisoblo(blist,iq),equasi,updo(blist,iq)
    END IF
    
    
  END IF
  
END DO

maxblist(iq) = blist

WRITE(7,'(a,i1,a,i5)') 'nr. of blocked states of type',iq,'=', maxblist(iq)


! sort on energy
DO blist=1,maxblist(iq)
  emin = edifblo(blist,iq)
  DO ib2=blist+1,maxblist(iq)
    IF(edifblo(ib2,iq)<emin) THEN
      emin = edifblo(ib2,iq)
      isav = amblock(blist,iq)
      amblock(blist,iq) = amblock(ib2,iq)
      amblock(ib2,iq) = isav
      isav = akblock(blist,iq)
      akblock(blist,iq) = akblock(ib2,iq)
      akblock(ib2,iq) = isav
      isav = aisoblo(blist,iq)
      aisoblo(blist,iq) = aisoblo(ib2,iq)
      aisoblo(ib2,iq) = isav
      isav = aistblo(blist,iq)
      aistblo(blist,iq) = aistblo(ib2,iq)
      aistblo(ib2,iq) = isav
      sav = edifblo(blist,iq)
      edifblo(blist,iq) = edifblo(ib2,iq)
      edifblo(ib2,iq) = sav
      sav = updo(blist,iq)
      updo(blist,iq) = updo(ib2,iq)
      updo(ib2,iq) = sav
    END IF
  END DO
END DO


RETURN
END SUBROUTINE search_block





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
!USE wavefr
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


END MODULE blocking
