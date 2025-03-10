c-----block-------------------------------------------------------
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

      include 'beta.inc'


      character*3  outnam

      integer blist   ! blocking

      logical blnprot       ! for blocking 
      logical blnneut

c--------------------------------------------------------------------------




       write(99,*)'if1',blnneut,blnprot,kblock(1),kblock(2)

      if((blnneut.or.blnprot)
     & .and.(kblock(1).eq.0).and.(kblock(2).eq.0))then
c
c     If the nucleus is even-odd or odd-even 
c     and kblock(iq)=0
c     --> Blocking
          write(99,*)'if2'
c      
c     static iteration to get the single particle energies
c

      call iteration_loop(outnam,x0sav,qact)

c
c     searching the levels which were next to the fermi energy 
c   
    


      amblock = 0         
      akblock = 0
      aisoblo = 0


      if(blnprot.and.blnneut)then
      do 7878  iq=1,2
      call search_block(blnprot,blnneut,iq) !for neutrons and protons
 7878 continue 
      elseif(blnprot)then   
      call search_block(blnprot,blnneut,1)  !only for protons
      elseif(blnneut)then
      call search_block(blnprot,blnneut,2)  !only for neutrons
      endif



c
c     write the common blocks for restart                            
c                                                                     

      call write_com(outnam)


c
c     Here the different states are blocked to search the one 
c     with the minimal energy
c

      if(blnprot.and.blnneut)then
      do 9898  iq=1,2
      call search_energy(iq,iter,outnam,x0sav,qact)
 9898 continue 
      do iq=1,2
      kblock(iq) = testblock(iq)          ! necessary because in the previous loop we have overwritten kblock 
      mblock(iq) = testblock2(iq)
      istblo(iq)= testblock3(iq) !test
      enddo   
      elseif(blnprot)then
      call search_energy(1,iter,outnam,x0sav,qact)
      elseif(blnneut)then
      call search_energy(2,iter,outnam,x0sav,qact)
      endif

      endif    


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

c
c----------------------------------------------------------------------
c


      bldist = 1000

      jstblo = 0

      isoblo = iq 


 






      if(kblock(iq).gt.2*mblock(iq)) then
        signup = 1.0D0
        mbact  = mblock(iq)
      else
        signup = -1.0D0
        mbact  = mblock(iq)-1
      endif
      if(isoblo.eq.1) then
        ninit = 1
        nfin  = nst(1)
      elseif(isoblo.eq.2) then
        ninit = nst(1)+1
        nfin  = nst(1)+nst(2)
      else
        stop ' invalid isospin in blocking routine COMP_BLOCK'
      endif


 
      
      write(15,*)'iq',iq

      bljst = 1
      

      do ist=ninit,nfin
        if(mpsi(ist).eq.mbact) then
          if(signup*(cmpn(1,ist)-cmpn(2,ist)).ge.0.0D0) then
             if(abs(spenrg(ist)-fer(isoblo)).lt.8.0D0) then



            bldist(bljst) = abs(spenrg(ist)-fer(isoblo)) 

            
            

            jstblo(bljst) = ist

            write(12,*) iq,fer(isoblo),bldist(bljst),jstblo(bljst)
            

            bljst = bljst + 1 

            endif
         endif
        endif
      enddo
      
      maxbljst = bljst - 1



      do 2145 bljst = 1,maxbljst                               ! searching the state with the demanded K
                                                               ! next to the fermi level
      distmin = min(bldist(1),bldist(2),bldist(3),
     & bldist(4),bldist(5),bldist(6),bldist(7),
     & bldist(8),bldist(9),bldist(10),bldist(11),
     & bldist(12),bldist(13),bldist(14),bldist(15))  
         



      if(bldist(bljst).eq.distmin) goto 2045


 2145 continue

 2045 continue 

      istblo(iq) = jstblo(bljst)

c      write(12,*) istblo(iq),jstblo(bljst),distmin



      do iunit=6,7
       write(iunit,'(a,6i5)')
     &  ' blocking found: iq,kblock(iq),mblock(iq),isoblo,istblo(iq)=',
     &    iq,kblock(iq),mblock(iq),istblo(iq)
      enddo

c
      return
      end





c
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
      
      absdiff(1) = 8.2
      absdiff(2) = 8.2




      blist   = 1


      write(17,*)'iq',iq
      

      nstm=nst(iq)

c      print levels

      do ist=1,nstm
          ine    = iokk(iq)+ist-1
          ind    = iorder(ine)+iokk(iq)-1      


      diff(ind)  =  spenrg(ind)-fer(iq)           ! distance to the fermi energy

      
      

      if(abs(diff(ind)).lt.absdiff(iq))then


      amblock(blist,iq) = mpsi(ind)               ! states which are in the vicinity of the fermi energy
      akblock(blist,iq) = 2*(mpsi(ind)+0.5)
      aisoblo(blist,iq) = icharg(ind)
      
      aistblo(blist,iq) = ind

      write(17,'(i5,i4,f6.1,3i5,f7.4,2f5.2,f9.3)')
     &       ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),
     &       nparz(1,ind),
     &       vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind)

      write(17,*) amblock(blist,iq),akblock(blist,iq),aisoblo(blist,iq),
     &  diff(ind)

      
      blist = blist +1                              ! counter for the blocking states
      
      else

      write(17,'(i5)')  
     &        ind
      endif
      
      
      enddo

      
      maxblist(iq) = blist - 1
      
      write(14,*) absdiff(1),maxblist(1),absdiff(2),maxblist(2)


      return 

      end


     
c
c------ search_energy
c
c     subroutine to find the minimal energy 
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
      include 'beta.inc'


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

      blenergy = 0



      do 1345 blist = 1,maxblist(iq),1

      istblo = 0   

c
c      read the common blocks                                          
c     
      x0sav  = x0inv 
      call read_com(iter,outnam)
      eold   = energy
     
      
      mblock(iq) = amblock(blist,iq)
      kblock(iq) = akblock(blist,iq) 
      isoblo = aisoblo(blist,iq) 

      istblo(iq) = aistblo(blist,iq)

     
      iqtest = mod(iq,2)+1              ! blocking only neutrons or protons
                                        ! necessary for odd-odd nuclei 
      kblock(iqtest) = 0





      write(19,*)'iq,iqtest,k',iq,iqtest,kblock(1),kblock(2)

   

 
      
      write(6,*)'Search the minimal energy' ,
     &       'mblock kblock isoblock istblo',
     &     mblock(iq),kblock(iq),isoblo,istblo(iq)
      
      

      write(7,*)'Search the minimal energy', 
     & 'mblock kblock isoblock istblo',
     &     mblock(iq),kblock(iq),isoblo,istblo(iq)
      

      call iteration_loop(outnam,x0sav,qact)
      


       blenergy(blist,iq) = energy       ! binding energy for the particular blocking state

       write(17,*) blist,mblock(iq),kblock(iq),isoblo,maxblist(iq)
     & ,blenergy(blist,iq)


 1345 continue



      do 1145 blist = 1,maxblist(iq),1                                    ! searching the miminal binding energy

      minenergy = min(blenergy(1,iq),blenergy(2,iq),blenergy(3,iq),
     & blenergy(4,iq),blenergy(5,iq),blenergy(6,iq),blenergy(7,iq),
     & blenergy(8,iq),blenergy(9,iq),blenergy(10,iq),blenergy(11,iq),
     & blenergy(12,iq),blenergy(13,iq),blenergy(14,iq),blenergy(15,iq))  
         
     

      if(blenergy(blist,iq).eq.minenergy) goto 1045





 1145 continue

 1045 continue 

      



      write(11,*) blist,iq,blenergy(blist,iq),beta2,amblock(blist,iq),
     &            akblock(blist,iq)


      mblock(iq) = amblock(blist,iq)
      kblock(iq) = akblock(blist,iq) 
      isoblo = aisoblo(blist,iq)

      istblo(iq) = aistblo(blist,iq) !test
      testblock3(iq)=istblo(iq) !test

      testblock(iq)=kblock(iq)
      testblock2(iq)=mblock(iq)
     
      write(11,*)'iq,klock',iq,kblock(iq),mblock(iq)


      return 

      end
