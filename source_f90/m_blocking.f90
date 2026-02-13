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
USE init, ONLY: ssetup,inipai
USE staticstep, ONLY: sstep

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
REAL(8) :: maxexenerg=2D0         !<Maximum s.p. excitation energy for blocking configurations
REAL(8) :: maxamomdiff=1D3        !<Maximum difference between atzimuthal and total angular momentum

LOGICAL, PARAMETER :: ttrack=.TRUE.  !<Switch for test output


CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: iteration_loop
!> @brief
!!Subroutine for one full static solution 
SUBROUTINE iteration_loop()

USE observ !, ONLY: itmax
USE constr
USE pairin, ONLY: findconjpart
USE util, ONLY: plot_2d

IMPLICIT NONE

INTEGER,PARAMETER :: ntermcount=3 ! nr. of termination counts required   
LOGICAL,PARAMETER :: tafter=.true.           ! switch to "afterburn"

REAL(8) :: oldoct,dnewoct,econprev,devh2prev,x0invmin,econs,qconsav,devh2s
REAL(8) :: pairfcsav(2),fxisav,h1xisav,hgradxisav,pairfacact
INTEGER :: iplold,ipold,iter,itesto,itmaxact,itermcount !,iwu,nchg,nmass,number_threads=1
INTEGER :: ipairsav,iteract
LOGICAL :: tfirst=.TRUE.,thfbsav
      
  thfbsav = ipair==6 .OR. ipair==8
  ipairsav = ipair
  IF(thfbsav) THEN
    fxisav = fxi
    h1xisav = h1xi
    hgradxisav = hgradxi
    pairfcsav = pairfc
  END IF

      eold   = energy
      itmaxact=itmax
!      tfirst = iqact0==0.AND.iqact2==0.AND.iqact3==0.AND.iqact4==0
!      IF(iflmb0==1) itmaxact=itcons2     ! prepare for later direct iteration
      itermcount=0                        ! set termination counter
      iterations: DO iter=0,itmaxact
        oldoct=beta3

        ! set actual HFB/BCS option
        IF(iter>iterbcs .OR. (tnobcs .AND. .NOT.tfirst) .OR. .NOT.thfbsav) THEN
          thfb = thfbsav
          ipair = ipairsav
        ELSE
          thfb = .FALSE.
          ipair = ipairsav-1
        END IF

        ! initialize gradual pairing enhancment 
        IF(thfbsav .AND. (tfresh .OR. tfirst)) THEN
          IF(iter<iterhfb) THEN
            pairfacact = 1D0+pairenhance*MAX(0,iterhfb-iter)/iterhfb
            fxi = fxisav*pairfacact
            h1xi = h1xisav*pairfacact
            hgradxi = hgradxisav*pairfacact
            IF(.NOT.tFayans) pairfc = pairfcsav*pairfacact
          ELSE
            fxi = fxisav
            h1xi = h1xisav
            hgradxi = hgradxisav
            IF(.NOT.tFayans) pairfc = pairfcsav
          END IF
        ELSE
        END IF

!NEW ODD
        ! update book-keeping of conjugated states
!        IF(toddterms .AND. MOD(iter,mod_upd_conjstates)==0 .AND. iter>iterhfb) THEN
        IF((MOD(iter,mod_upd_conjstates)==1 .AND. iter>iterbcs) .OR. &
           iter==iterbcs .OR. iter==iterhfb) THEN
          CALL findconjpart(1)
          CALL findconjpart(2)
          CALL comp_block()            ! ???
        END IF
!END NEW ODD


        IF(MOD(iter,iprint)==1) THEN
          WRITE(6,*)' Iter, Energy,  beta2,   beta3,log(D^2E),log(dE/dIt),log(dQ/dIt), zcm x0inv'
          WRITE(7,*)' Iter, Energy,  beta2,   beta3,log(D^2E),log(dE/dIt),log(dQ/dIt), zcm  x0inv'
        END IF
!        IF(iter /= 0) THEN
        IF(iter >2) THEN
          WRITE(6,5050)  &
              iter,energy,beta2,beta3, LOG10(MAX(devh2,smal)),&
              LOG10(MAX(ABS(econ),smal)),LOG10(MAX(ABS(qcon),smal)), zcm,x0inv
          WRITE(7,5050)  &
              iter,energy,beta2,beta3, LOG10(MAX(devh2,smal)),&
              LOG10(MAX(ABS(econ),smal)),LOG10(MAX(ABS(qcon),smal)), zcm,x0inv
          IF(check_border_density()) WRITE(6,*)'***WARNING: Grid might be too small. Density at border >1E-6!'
          CALL flush(7)
          5050 FORMAT(i5,X,f9.2,2(X,f8.3),3(3X,f7.3),X,1pg12.4,0pf7.3)
        END IF
        
!     set addition-weights for mixing density and potential,
!     first iteration needs a refresh
        
        IF(iter <= 3) THEN
          adyrho = 1D0
          adypot = 1D0
        ELSE
          adyrho = addrho
          adypot = addpot

        END IF
        
!IF(thfb)  WRITE(6,'(a,i6,1pg13.5)') 'SUMUV10:',iter,SUM(SQRT(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2))
        CALL sstep(iter)
!IF(thfb)  WRITE(6,'(a,i6,1pg13.5)') 'SUMUV11:',iter,SUM(SQRT(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2))
        CALL sinfo(iter)
!IF(thfb)  WRITE(6,'(a,i6,1pg13.5)') 'SUMUV12:',iter,SUM(SQRT(vsqf(nst(1)+1:nst(1)+nst(2))-vsqf(nst(1)+1:nst(1)+nst(2))**2))
        
        dnewoct=beta3
        iteract = iter
        
!        IF(econ < 10.0*errenr .AND. kblock > 0 .AND. SUM(istblo) <= 0)  &
!            CALL comp_block                           !NEW ODD: activate with KBLOCK-entry
        itert  = itert+1
        IF(iter>iterhfb) THEN
          IF(iter > itpini+3 .AND. tswicon) THEN
            IF (econ < errenr .AND. devh2 < errvar .AND. &
                ABS(qcon) < qerr .AND. ABS(zcm) < 1E-3  .AND. &
                ABS(oldoct-dnewoct) < 1E-3) itermcount=1+itermcount
          ELSE IF (.NOT.tswicon) THEN
            IF(econ < errenr .AND. devh2 < errvar) THEN
              itermcount=1+itermcount
              WRITE(71,'(a,i4)') 'itermcount enhanced:',itermcount
            ELSE
              IF(itermcount>0) WRITE(71,'(a,i4)') 'itermcount reset'
              itermcount = 0
            END IF
          END IF
          IF(itermcount>0) THEN
            WRITE(71,'(a,i3)') 'itermcount=',itermcount
            WRITE(7,'(a,i3)') 'itermcount=',itermcount
            WRITE(6,'(a,i3)') 'itermcount=',itermcount
          END IF
          IF( itermcount==ntermcount) EXIT iterations
        END IF
        
        IF(tvaryx0) THEN
          IF(energy<econprev .AND. devh2<(devh2prev*(1.0d0-1.0d-5))) THEN
             x0inv=x0inv*1.005
          ELSE
             x0inv=x0inv*0.9
          END IF
          IF(x0inv<x0invmin) x0inv=x0invmin
          econprev=energy
          devh2prev=devh2
       END IF
        
        IF(iter > itcons0) THEN
          qlamb0  = zero
          tswicon0 = .FALSE.
        END IF
        IF(iter > itcons2) THEN
          qlamb2  = zero
          tswicon2 = .FALSE.
        END IF
        IF(iter > itcons3) THEN
          qlamb3  = zero
          tswicon3 = .FALSE.
        END IF
        IF(iter > itcons4) THEN
          qlamb4  = zero
          tswicon4 = .FALSE.
        END IF
         
      END DO iterations
      
!**********************************************************************
!                                                                     *
!     end of static iteration                                         *
!                                                                     *
!**********************************************************************
      
!     one more iteration to put observables up to date
      
!      IF(tafter .AND. iflmb0 /= 1) THEN
      IF(tafter) THEN
        iter = iteract
        ipold  = iprint
        iprint = 1
        iplold = ipri_fields
        itesto = itest
        IF(ipri_fields>0) ipri_fields  = 1
        itest  = 0
        econs  = econ
        qconsav  = qcon
        devh2s = devh2
        adyrho = 1.0
        adypot = 1.0
        CALL sstep(iteract)                       ! "afterburn"
        CALL sinfo(itmax)
        iprint = ipold
        ipri_fields  = iplold
        itest  = itesto
        econ   = econs
        qcon   = qconsav
        devh2  = devh2s
      END IF
      IF(tpri_2d)  CALL plot_2d()
      CALL occ_lz()
      
      CALL user_interface_1()        ! user defined corrections
      
!     optionally compute localization criterion
      
      IF(klocal == 1) CALL localize()
      


      angblock=0D0
      IF(blnprot .AND. .NOT.blnneut .AND. SUM(istblo)>0) angblock(1)=angularmom(istblo(1))  
      IF(.NOT.blnprot .AND. blnneut .AND. SUM(istblo)>0) angblock(2)=angularmom(istblo(2))  
      IF(blnprot .AND. blnneut .AND. SUM(istblo)>0) THEN
        angblock(1)=angularmom(istblo(1))  
        angblock(2)=angularmom(istblo(2))  
      END IF
      CALL write_ensurf(iter)

      tfirst=.FALSE.

      IF(tpri_conver) WRITE(71,'(1x/1x)')

      RETURN

END SUBROUTINE iteration_loop



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


  
IF(blnprot .AND. .NOT.blnneut) CALL block1(1)
IF(.NOT.blnprot .AND. blnneut) CALL block1(2)
IF(blnprot .AND. blnneut)  CALL block2()


END SUBROUTINE scan_block


!---------------------------------------------------------------------------  
! DESCRIPTION: block1
!> @brief
!! Computes the SHF solution for one blocked nucleon configuration.
SUBROUTINE block1(iq)

USE constr, ONLY: tfresh

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
  
  
  WRITE(6,'(a,4(i4),1pg13.5)')'Search the ground state for:  &
      mblock kblock isoblock istblo',  &
      mblock(iq),kblock(iq),isoblo,istblo(iq),edifblo(blist,iq)
  WRITE(7,'(a,5i4,1pg13.5)')'Search the ground state for:  &
      mblock kblock isoblock istblo',  &
      mblock(iq),kblock(iq),isoblo,istblo(iq),istcblo(iq),edifblo(blist,iq)
  CALL flush(7)
  OPEN(90,STATUS='unknown',FORM='formatted',FILE='ensurf.'//outnam,  &
          POSITION='APPEND')
  WRITE(90,'(a,4i5,f8.3)') 'Blocking  isoblo,mblock,kblock,istblo,equasi=',&
    isoblo,mblock(iq),kblock(iq),istblo(iq),edifblo(blist,iq)
  CLOSE(90)

  IF(tfresh) THEN
    CALL ssetup                    
    CALL inipai(npart(1)+npart(2))
  END IF

  CALL iteration_loop()
  
  
  
END DO

END SUBROUTINE block1

!---------------------------------------------------------------------------  
! DESCRIPTION: block2
!> @brief
!! Computes the SHF solution for one proton and one neutron blocked.
SUBROUTINE block2()

USE constr, ONLY: tfresh

IMPLICIT NONE

INTEGER :: iq,maxblist_total,ip,ip1,ip2,blist1,blist2,ipsav(2)
REAL(8) :: energmin,energmax
REAL(8), ALLOCATABLE :: total_energy(:)
INTEGER, ALLOCATABLE :: ipoint(:,:)

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
  
  WRITE(7,'(a,8(i4))') 'Ground state for:  &
        mblock kblock istblo',mblock,kblock,istblo
  IF(ttrack) THEN
    WRITE(6,'(a,8(i4))')'Search the minimal energy:  &
        mblock kblock istblo',mblock,kblock,istblo
    WRITE(*,*) 'before iteration loop in BLOCK2'
    WRITE(7,*) 'before iteration loop in BLOCK2'
  END IF
  OPEN(90,STATUS='unknown',FORM='formatted',FILE='ensurf.'//outnam,  &
          POSITION='APPEND')
  WRITE(90,'(a,6i5,2f8.3)') 'Double blocking mblock,kblock,istblo,equasi=',&
    mblock,kblock,istblo,edifblo(blist1,1),edifblo(blist2,2)
  CLOSE(90)
  
  CALL flush(7)

  IF(tfresh) THEN
    CALL ssetup                    
    CALL inipai(npart(1)+npart(2))
  END IF

  CALL iteration_loop()
  
  
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

USE observ, ONLY: angularmom

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
  IF(equasi < maxexenerg) THEN  

    ! exclude a couple of unwanted cases
    IF(( ((blnprot.AND.blnneut).AND.iq==1)  &
        .OR.(.NOT.(blnprot.AND.blnneut))   ) .AND. mpsi(ind)<0 ) CYCLE
    IF(ABS(mpsi(ind)+0.5D0-angularmom(ind))>maxamomdiff) CYCLE   

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
