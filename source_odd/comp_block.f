c----block-------------------------------------------------------
c
      subroutine block(blnprot,blnneut,outnam,x0sav)
*
***********************************************************************
*    
*     subroutine to calculate the optimal blocked state for an odd 
*     nucleus by minimizing the binding energy
*    
*                                                                     *
***********************************************************************
      
      include 'paramr.inc'

      include 'constr.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'inipar.inc'
      include 'observ.inc'
      include 'regime.inc'

      include 'blocking.inc'



      character*3  outnam

      integer blist   ! blocking

      logical blnprot       ! for blocking 
      logical blnneut

c--------------------------------------------------------------------------



      IF(ttrack) THEN
        WRITE(17,'(a,2L2,2i2)') 'BLOCK entered:',
     &    blnprot,blnneut,kblock
        WRITE(*,'(a,2L2,2i2)') 'BLOCK entered:',
     &    blnprot,blnneut,kblock
      END IF
      if((blnneut.or.blnprot)
     & .and.(kblock(1).eq.0).and.(kblock(2).eq.0))then   ! search-branch

c     make sure that conjugate states are known

c
c     If the nucleus is even-odd or odd-even 
c     and kblock(iq)=0
c     --> Blocking

c      
c     static iteration to get the single particle energies
c
      call inicms
      if (swicon.gt.oh) call inicon
      call findconjpart(1)
      call findconjpart(2)
      IF(ttrack) WRITE(17,'(a,2L2,2i2)') 'BLOCK: before pre-iteration'
      WRITE(7,'(a,200i4)') 
     &    'conj. states:',itpartn(nst(1)+1:nst(1)+nst(2))
      CALL flush(17)
      IF(ttrack) THEN
         WRITE(*,*) '1. iteration loop in COMP_BLOCK'
         WRITE(7,*) '1. iteration loop in COMP_BLOCK'
         WRITE(17,*) '1. iteration loop in COMP_BLOCK'
      END IF
      call iteration_loop(outnam,x0sav,qact)
      IF(ttrack) WRITE(17,*) 'BLOCK: preparatory static iteration '
      CALL flush(17)

c
c     searching the levels which were next to the fermi energy 
c   

      amblock = 0          
      akblock = 0
      aisoblo = 0

      if(blnprot.and.blnneut)then
        call search_block(blnprot,blnneut,1) 
        call search_block(blnprot,blnneut,2) 
      elseif(blnprot)then   
        call search_block(blnprot,blnneut,1)  !only for protons
      elseif(blnneut)then
        call search_block(blnprot,blnneut,2)  !only for neutrons
      endif



c
c     write the common blocks for restart                            
c                                                                     

      call write_com(outnam)

      IF(ttrack) WRITE(17,*) 'BLOCK: after search, testblock=',
     &    testblock,testblock2

c
c     Here the different states are blocked to search the one 
c     with the minimal energy
c

      if(blnprot.and.blnneut) then
!      do 9898  iq=1,2
        call search_energy2(iter,outnam,x0sav,qact)
! 9898 continue 
        IF(ttrack) WRITE(17,*) 'odd-odd branch: k,m=',
     &     testblock,testblock2
        do iq=1,2
          kblock(iq) = testblock(iq)      ! necessary because in the previous loop we have overwritten kblock 
          mblock(iq) = testblock2(iq)
        enddo   
      elseif(blnprot)then
        IF(ttrack) WRITE(17,'(a)') 'odd-even branch'
        call search_energy(1,iter,outnam,x0sav,qact)
      elseif(blnneut)then
        IF(ttrack) WRITE(17,'(a)') 'even-odd branch'
        call search_energy(2,iter,outnam,x0sav,qact)
      endif
      WRITE(7,'(a,2i5)') 'after SEARCH_ENERGY: istblo,itscblo=',
     &    istblo(2),istcblock(2)
      IF(ttrack) THEN
        WRITE(17,'(a,4i4)') 'after SEARCH_ENERGY: kblock,istblo=',
     &    kblock,istblo          
        CALL flush(17)
      END IF

      else                 ! alternative branch

      istcblock(1)=itpartn(istblo(1))
      istcblock(2)=itpartn(istblo(2))
      IF(ttrack) THEN
         WRITE(17,*) 'BLOCK: directly to iteration, params=',
     &       kblock,mblock,istblo,istcblock,isoblo
         WRITE(*,*) 'BLOCK: directly to iteration, params=',
     &       kblock,mblock,istblo,istcblock,isoblo
      END IF
      IF(ttrack) THEN
         WRITE(*,*) '2. iteration loop in COMP_BLOCK'
         WRITE(7,*) '2. iteration loop in COMP_BLOCK'
         WRITE(17,*) '2. iteration loop in COMP_BLOCK'
      END IF
      call iteration_loop(outnam,x0sav,qact)
!        call search_energy(isoblo,iter,outnam,x0sav,qact)

      endif    

      WRITE(*,*) 'BLOCK terminated. return.'

      return 
      
      end

c
c----------------------------------------------------------------------
c

c-----comp_block-------------------------------------------------------
c
      subroutine comp_block(iq)
*
***********************************************************************
*                                                                     *
*     Computes the blocked state after some pre-iterations.           *
*     Blocked is a state with 'kblock', 'mblock' and 'isoblock'       *
*     as given in 'levels.inv'.                                       *
*     The routine finds out what number this state has in the         *
*     normal ordering sequence of the code.                           *
*                                                                     *
***********************************************************************
      include 'paramr.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      

      double precision bldist(20) 
      integer          bljst
      integer          maxbljst
      integer          jstblo(20)
      double precision distmin

      integer ud

c
c----------------------------------------------------------------------
c


      bldist = 1000

      jstblo = 0

      isoblo = iq 


 
      if(isoblo.eq.1) then
        ninit = 1
        nfin  = nst(1)
      elseif(isoblo.eq.2) then
        ninit = nst(1)+1
        nfin  = nst(1)+nst(2)
      else
        stop ' invalid isospin in blocking routine COMP_BLOCK'
      endif


 
      write(17,'(a,5i3)') 'COMP_BLOCK: m,k,iso=',
     &    mblock(1),mblock(2),kblock(1),kblock(2),isoblo
      write(*,'(a,5i3)') 'COMP_BLOCK: m,k,iso=',
     &    mblock(1),mblock(2),kblock(1),kblock(2),isoblo
!      IF(ttrack) write(17,'(a,5i3)') 'COMP_BLOCK: m,k,iso=',
!     &    mblock(1),mblock(2),kblock(1),kblock(2),isoblo

      bljst = 1
      




      do ist=ninit,nfin


         if((cmpn(1,ist)-cmpn(2,ist)).gt.0)then
            ud = 1
         else   
            ud = 0
         endif   


        if(mpsi(ist).eq.mblock(iq)) then
           if(ud.eq.spinblo(iq))then
             if(abs(spenrg(ist)-fer(isoblo)).lt.8.0D0) then



            bldist(bljst) = abs(spenrg(ist)-fer(isoblo)) 

            
            

            jstblo(bljst) = ist

            bljst = bljst + 1 

            endif
         endif
        endif
      enddo
      
      maxbljst = bljst - 1



      do 2145 bljst = 1,maxbljst                               ! searching the state with the demanded K
                                                               ! next to the fermi level
      distmin = minval(bldist)  
         



      if(bldist(bljst).eq.distmin) goto 2045


 2145 continue

 2045 continue 




      istblo(iq) = jstblo(bljst)

      istcblock(iq)=itpartn(istblo(iq)) !search the time conjugate state of the blocked state




      do iunit=6,7
       write(iunit,'(a,6i5)')
     &  ' blocking found: iq,kblock(iq),mblock(iq),isoblo,istblo(iq)=',
     &    iq,kblock(iq),mblock(iq),istblo(iq)
      enddo
      IF(ttrack) write(17,'(a,6i5)')
     &  ' blocking found: iq,kblock(iq),mblock(iq),isoblo,istblo(iq)=',
     &    iq,kblock(iq),mblock(iq),istblo(iq)

c
      return
      end






c------ search_block
c
c     subroutine to find the energy levels which were next to the
c     fermi energy
c
cccccccccccccccccccccccccccccccccc

      subroutine search_block(blnprot,blnneut,iq)
      

      include 'paramr.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      


      logical blnprot
      logical blnneut

      dimension absdiff(1:2)
      integer blist
      integer iq

      

c     
c     searching the energy levels next to the fermi energy 
c
c--------------------------------------------------------------------
c
      
!     set the energy interval for the search

!      absdiff(1) = 2.0
!      absdiff(2) = 2.2
!      absdiff(1) = 6.0
!      absdiff(2) = 6.2

!      absdiff(1) = 3.2
!      absdiff(2) = 3.2
      absdiff(1) = 5.2
      absdiff(2) = 5.2
!      absdiff(1) = 0.5
!      absdiff(2) = 0.5


      blist   = 1


      write(17,*)'SEARCH BLOCK: iq=',iq
      IF(ttrack) WRITE(17,'(a,i2)') 'SEARCH_BLOCK: iq=',iq      

      nstm=nst(iq)

c      print levels

      do ist=1,nstm
        ine    = iokk(iq)+ist-1
        ind    = iorder(ine)+iokk(iq)-1      


!old      diff(ind)  =  spenrg(ind)-fer(iq)           ! distance to the fermi energy
        diff(ind)=SQRT((spenrg(ind)-fer(iq))**2+deltaf(ind)**2)  ! QP energy


      
      

!      if((abs(diff(ind)).lt.absdiff(iq)).and.mpsi(ind).ge.0)then  ! only one of the time conjugate partners is blocked
      if((abs(diff(ind)).lt.absdiff(iq)))then  ! both partners can be blocked
!      IF( mpsi(ind).GE.0 ) WRITE(7,'(a,2i4,f7.2)') 
!     &   'SEARCH: ind,mpsi,j=',ind,mpsi(IND),angular(ind)
      IF(( ((blnprot.AND.blnneut).AND.iq==1)
     &    .OR.(.NOT.(blnprot.AND.blnneut))   )
     &   .AND. mpsi(ind)<0 ) CYCLE
      IF(ABS(mpsi(IND)+0.5D0-angular(ind))>0.6D0) CYCLE
!      IF(ABS(mpsi(IND)+0.5D0-angular(ind))>0.3D0) CYCLE
!      IF(mpsi(ind).NE.0 ) CYCLE
!      IF(mpsi(ind).NE.1 ) CYCLE
!      IF(mpsi(ind).NE.2 ) CYCLE
!      IF(mpsi(ind).NE.5 ) CYCLE
!      IF(mpsi(ind).NE.4 ) CYCLE
!     &   .AND. mpsi(ind)<0 ) CYCLE
!      IF(( ((blnprot.AND.blnneut).AND.iq==1)
!     &    .OR.(.NOT.(blnprot.AND.blnneut))   )
!     &   .AND. (mpsi(ind)<0. 
!     &          OR. mpsi(ind).NE.4 )) CYCLE
!      IF(( ((blnprot.AND.blnneut).AND.iq==1)
!     &    .OR.(.NOT.(blnprot.AND.blnneut))   )
!     &   .AND. (mpsi(ind)<0.OR.mpsi(ind).NE.4)) CYCLE
!      IF(mpsi(ind)<0 .OR. mpsi(ind).NE.5) CYCLE
!      IF(mpsi(ind)<0 .OR. mpsi(ind).NE.4) CYCLE
!      IF(mpsi(ind)<0 .OR. mpsi(ind).NE.1) CYCLE
!      IF(mpsi(ind)<0 .OR. (mpsi(ind).NE.3 .AND. iq==1)) CYCLE
      amblock(blist,iq) = mpsi(ind)               ! states which are in the vicinity of the fermi energy
      akblock(blist,iq) = 2*(mpsi(ind)+0.5)
      aisoblo(blist,iq) = icharg(ind)
      edifblo(blist,iq) = ABS(diff(ind))



      updo(blist,iq)    = (cmpn(1,ind)-cmpn(2,ind)) ! to characterize the state in subroutine comp_block


      aistblo(blist,iq) = ind

      WRITE(17,'(a,i4)') 'success: blist=',blist
      WRITE(17,'(a)') 'ind,charge,m,3*parity,v2,comps,spenerg'
      write(17,'(i5,i4,f6.1,3i5,f7.4,2f5.2,f9.3)')
     &       ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),
     &       nparz(1,ind),
     &       vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind)

      WRITE(17,'(a)') 'amblock,akblock,aisoblock,diff,updo:'
      write(17,'(5(1pg13.5))') amblock(blist,iq),akblock(blist,iq),
     &   aisoblo(blist,iq),diff(ind),updo(blist,iq)

      
      blist = blist +1                              ! counter for the blocking states
      
      else

!      write(17,'(a,i5)')  
!     &        "no success: nr. state=",ind
      endif
      
      
      enddo

      
      maxblist(iq) = blist - 1
      
      WRITE(17,'(a,i5)') 'maxblist=', maxblist(iq)


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


      return 
      end


     
c
c------ search_energy
c
c     subroutine to find the minimal energy for given isospin 'iq'
c    
c
cccccccccccccccccccccccccccccccccc

      subroutine search_energy(iq,iter,outnam,x0sav,qact)
      

      include 'paramr.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'constr.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'inipar.inc'
      include 'regime.inc'


      CHARACTER(3) :: outnam
      logical blnprot
      logical blnneut

      double precision minenergy


      integer blist

      dimension blenergy(1:nblint,1:2)
      

      integer iqtest

      character         outnump*1
      character         outnumn*1


c
c     searching the next energy level above the fermi energy 
c
c--------------------------------------------------------------------
c

!      WRITE(7,'(a,200i4)') 
!     &    'conj. before loop:',itpartn(nst(1)+1:nst(1)+nst(2))
      blenergy = 0

      do 1345 blist = 1,maxblist(iq)

!      IF(blist==1) WRITE(7,'(a,200i4)') 
!     &    'conj. before read:',itpartn(nst(1)+1:nst(1)+nst(2))

      blenergy = 0
      istblo = 0   

c
c      read the common blocks                                          
c     

      x0sav  = x0inv 
      call read_com(iter,outnam)
      eold   = energy
     
!      WRITE(7,'(a,i4/200i4)') 
!     &    'conj. in loop:',blist,itpartn(nst(1)+1:nst(1)+nst(2))
      
      mblock(iq) = amblock(blist,iq)
      kblock(iq) = akblock(blist,iq) 
      isoblo = aisoblo(blist,iq) 
      ediffblock(iq) = edifblo(blist,iq)

      

      istblo(iq) = aistblo(blist,iq)
      istcblock(iq) = itpartn(istblo(iq)) ! time conj. partner
      

      iqtest = mod(iq,2)+1              ! blocking only neutrons or protons
                                        ! necessary for odd-odd nuclei 
      kblock(iqtest) = 0


      write(6,'(a,4(i4),1pg13.5)')'Search the minimal energy: 
     &       mblock kblock isoblock istblo',
     &     mblock(iq),kblock(iq),isoblo,istblo(iq),ediffblock(iq)
      
      

      write(7,'(a,5i4,1pg13.5)')'Search the minimal energy: 
     & mblock kblock isoblock istblo',
     &     mblock(iq),kblock(iq),isoblo,istblo(iq),istcblock(iq),
     &     ediffblock(iq)
      IF(ttrack) write(17,'(a,4(i4))')'Search the minimal energy: 
     & mblock kblock isoblock istblo',
     &     mblock(iq),kblock(iq),isoblo,istblo(iq)
      CALL flush(17)

      IF(ttrack) THEN
         WRITE(*,*) '3. iteration loop in COMP_BLOCK'
         WRITE(7,*) '3. iteration loop in COMP_BLOCK'
         WRITE(17,*) '3. iteration loop in COMP_BLOCK'
      END IF
      call iteration_loop(outnam,x0sav,qact)
      

       ang = 0
       call calangularmomentum(istblo(iq),ang)
       blenergy(blist,iq) = energy       ! binding energy for the particular blocking state
       iparity = npar(istblo(iq))*nrpar(istblo(iq))
       WRITE(17,'(a)') 
     &    'SEARCH_ENERGY: blist,m,k,iso,maxb,energ,ist,par,ang'
       write(17,'(5i4,1pg13.5,2i4,2(1pg13.5))') 
     &    blist,mblock(iq),kblock(iq),isoblo,maxblist(iq)
     &    ,blenergy(blist,iq),istblo, iparity,ang


 1345 continue



       minenergy = minval(blenergy)



      do 1145 blist = 1,maxblist(iq)                                    ! searching the miminal binding energy

      if(blenergy(blist,iq).eq.minenergy) goto 1045

 1145 continue

 1045 continue 

      
      mblock(iq) = amblock(blist,iq)
      kblock(iq) = akblock(blist,iq) 
      isoblo = aisoblo(blist,iq)

      istblo(iq) = aistblo(blist,iq) !test
      istcblock(iq) = itpartn(istblo(iq))

       if(updo(blist,iq).gt.0)then
       spinblo(iq) = 1
       else
       spinblo(iq) = 0 
       endif   


      testblock(iq)=kblock(iq) ! for odd-odd nuclei
      testblock2(iq)=mblock(iq)!


      return 

      end


c
c------ search_energy
c
c     subroutine to find the minimal energy for doubly blocked odd-odd nuclei
c    
c
cccccccccccccccccccccccccccccccccc

      subroutine search_energy2(iter,outnam,x0sav,qact)
      

      include 'paramr.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'constr.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'inipar.inc'
      include 'regime.inc'


      CHARACTER(3) :: outnam
      logical blnprot
      logical blnneut

      double precision minenergy


      integer blist1,blist2

      dimension blenergy(1:nblint,1:2)
      

      integer iqtest

      character         outnump*1
      character         outnumn*1

      REAL(8) :: angularmom(2)
      
      REAL(8), ALLOCATABLE :: total_energy(:)
      INTEGER, ALLOCATABLE :: ipoint(:,:)
      INTEGER :: ipsav(1:2)
c
c     searching the next energy level above the fermi energy 
c
c--------------------------------------------------------------------
c

! new structure:
!     combine into one linear field
      maxblist_total=maxblist(1)*maxblist(2)
      ALLOCATE(total_energy(maxblist_total),ipoint(1:2,maxblist_total))
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
      WRITE(17,'(a,4i4)') 'nr. of states:',maxblist(:),maxblist_total,ip
      WRITE(17,'(a)') 'energy before sorting of 2qp states:'
      WRITE(17,'(400(2i5,3(1pg13.5)/))') 
     &     (ipoint(:,ip1),total_energy(ip1),
     &     edifblo(ipoint(1,ip1),1),edifblo(ipoint(2,ip1),2),
     &     ip1=1,maxblist_total)
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
      WRITE(17,'(a)') 'final energy sorting of 2qp states:'
      WRITE(17,'(2i5,1pg13.5)') 
     &  (ipoint(:,ip1),total_energy(ip1),ip1=1,maxblist_total)
      energmin = MINVAL(total_energy(:))
      energmax = MAXVAL(total_energy(:))
      energmax = energmin+(energmax-energmin)/2D0
      
!     now scan the states
      DO ip1=1,maxblist_total
        blist1 = ipoint(1,ip1)
        blist2 = ipoint(2,ip1)
        IF(total_energy(ip1)>energmax) EXIT

!  old loop
!      DO isum=0,maxblist(2)-1
!      DO idif=-isum,isum,2
!      blist1=(isum+idif)/2  +1
!      IF(blist1>maxblist(1)) CYCLE
!      blist2=(isum-idif)/2  +1
!      IF(blist2>maxblist(2)) CYCLE



        istblo = 0   

c
c      read the common blocks                                          
c     

        x0sav  = x0inv 
        call read_com(iter,outnam)
        eold   = energy
     
      
        mblock(1) = amblock(blist1,1)
        kblock(1) = akblock(blist1,1) 
        mblock(2) = amblock(blist2,2)
        kblock(2) = akblock(blist2,2) 
!        isoblo = aisoblo(blist,iq) 
     
        istblo(1) = aistblo(blist1,1)
        istcblock(1) = itpartn(istblo(1)) ! time conj. partner
        istblo(2) = aistblo(blist2,2)
        istcblock(2) = itpartn(istblo(2)) ! time conj. partner
        ediffblock(1) = edifblo(blist1,1)      
        ediffblock(2) = edifblo(blist2,2)      


!        iqtest = mod(iq,2)+1              ! blocking only neutrons or protons
                                        ! necessary for odd-odd nuclei 
!        kblock(iqtest) = 0


        write(6,'(a,8(i4))')'Search the minimal energy: 
     &     mblock kblock istblo',mblock,kblock,istblo
        write(7,'(a,8(i4))')'Search the minimal energy: 
     &     mblock kblock istblo',mblock,kblock,istblo
        IF(ttrack) WRITE(17,'(a,8(i4))')
     &     'SEARCH_ENERGY2: minimal energy:'//
     &     ' mblock kblock istblo',mblock,kblock,istblo
      
     

        IF(ttrack) THEN
           WRITE(*,*) '4. iteration loop in COMP_BLOCK'
           WRITE(7,*) '4. iteration loop in COMP_BLOCK'
           WRITE(17,*) '4. iteration loop in COMP_BLOCK'
        END IF
        call iteration_loop(outnam,x0sav,qact)
      

         ang = 0
         call calangularmomentum(istblo(1),angularmom(1))
         call calangularmomentum(istblo(2),angularmom(2))
         blenergy(blist1,1) = energy       ! binding energy for the particular blocking state
         iparity = npar(istblo(iq))*nrpar(istblo(iq))
         WRITE(17,'(a)') 
     &    'SEARCH_ENERGY2: blist,m,k,iso,maxb,energ,ist,par,ang,rmsc'
         write(17,'(6i4,1pg13.5,2i4,3(1pg13.5))') 
     &    blist1,blist2,mblock,kblock,energy,istblo,angularmom,rmsc
         CALL flush(17)

       END DO
!       END DO     ! only for old loop



       minenergy = minval(blenergy)



      do 1145 blist = 1,maxblist(iq)                                    ! searching the miminal binding energy

      if(blenergy(blist,iq).eq.minenergy) goto 1045

 1145 continue

 1045 continue 

      
      mblock(iq) = amblock(blist,iq)
      kblock(iq) = akblock(blist,iq) 
      isoblo = aisoblo(blist,iq)

      istblo(iq) = aistblo(blist,iq) !test
      istcblock(iq) = itpartn(istblo(iq))

       if(updo(blist,iq).gt.0)then
       spinblo(iq) = 1
       else
       spinblo(iq) = 0 
       endif   


      testblock(iq)=kblock(iq) ! for odd-odd nuclei
      testblock2(iq)=mblock(iq)!

      DEALLOCATE(total_energy,ipoint)

      return 

      end
