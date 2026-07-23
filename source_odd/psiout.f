c---------------------
      subroutine psiout(outnam)  
*
***********************************************************************
*        V.O. Nesterenko, 24.07.2002, San Paulo                       *                                                         *
*        Subroutine psiout:                                           *
*        -- assigns Nillson quantum numbers to the s-p levels,        *
*        -- computes the complete space parity ntpar,                 *                 
*        -- writes s-p w.f., energies, quantum numbers and occupation *
*           numbers to the output file SPOUT.                         *
*        -- writes the local s-p potential to the output file         *
*           'sppoten'                                                 *
*        ksev - max. number of s-p states of EVEN-parity and of fixed *
*               isospin and  magnetic number Lambda                   * 
*        ksod - max. number of s-p states of ODD-parity and of fixed  *
*               isospin and  magnetic number Lambda                   *  
*        klev - max. number of Lambda's for EVEN-parity states        *
*               of a fixed isospin                                    *
*        klod - max. number of Lambda's for ODD-parity states         *
*               of a fixed isospin                                    *        
*        nqn(ist) - Nillson quantum numbers of the state 'ist'        *                                      
*  iorder (i) - numbers of the ordered levels (protons from 1 to      *
*               nst(1) and neutrons from 1 to to nst(2))              *
*               For example, if there are 20 proton and 30 neutron    *
*               levels, then iorder changes in the intervals 1-20 for *
*               protons and 1-30 for neutrons. At the same time 'ist' *
*               changes from 1 to 50.                                 *
***********************************************************************
      include 'paramr.inc'
      include 'levels.inc'
      include 'observ.inc'
      include 'wavefr.inc'
      include 'regime.inc' 
      include 'dynpar.inc'
      include 'fields.inc' 
c      common/ord/iord2(npsid),spen2(npsid)  
      character*3  outnam
c                
      parameter (nsemax=25,nlemax=8,nsomax=30,nlomax=10)
c      dimension ntpar(npsid),
      dimension ak(npsid) 
      dimension niprev(nsemax,nlemax),niprod(nsomax,nlomax)
      dimension niobev(nsemax,nlemax),niobod(nsomax,nlomax)
      dimension icount(2,10),nqn(npsid)
c
c     nsemax,nsomax - numbers of [N n_z \Lambda] for a given K^{\pi}
c     nlemax,nlomax - numbers of K^{\pi} for a given parity. 
c
c     Nillson quantim numbers:
c     prolate, even parity:  
      data niprev/0,220,211,200,440,431,420,411,400, ! K=1/2+
     & 660,651,640,631,620,611,600,
     & 880,871,860,851,840,831,820,811,800,  
     & 211,202,431,422,411,402,                      ! K=3/2+
     & 651,642,631,622,611,602,
     & 871,862,851,842,831,822,811,802,
     & 0,0,0,0,0,           
     & 202,422,413,402,642,633,622,613,602,          ! K=5/2+
     & 862,853,842,833,822,813,802,
     & 0,0,0,0,0,0,0,0,0,
     & 413,404,633,624,613,604,                      ! K=7/2+
     & 853,844,833,824,813,804,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,
     & 404,624,615,604,844,835,824,815,804,          ! K=9/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,    
     & 615,606,835,826,815,806,                      ! K=11/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,
     & 606,826,817,806,                              ! K=13/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,
     & 817,808,                                      ! K=15/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0/
c
c    oblate, even parity:   
       data niobev/0,200,211,220,400,411,420,431,440, ! K=1/2+
     & 600,611,620,631,640,651,660,
     & 800,811,820,831,840,851,860,871,880,  
     & 202,211,402,411,422,431,                      ! K=3/2+
     & 602,611,622,631,642,651,
     & 802,811,822,831,842,851,862,871,
     & 0,0,0,0,0,           
     & 202,402,413,422,602,613,622,633,642,          ! K=5/2+
     & 802,813,822,833,842,853,862,
     & 0,0,0,0,0,0,0,0,0,
     & 404,413,604,613,624,633,                      ! K=7/2+
     & 804,813,824,833,844,853,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,
     & 404,604,615,624,804,815,824,835,844,          ! K=9/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,    
     & 606,615,806,815,826,835,                      ! K=11/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,
     & 606,806,817,826,                              ! K=13/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,
     & 808,817,                                      ! K=15/2+
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0/      
c
c     prolate, odd parity: 
      data niprod/110,101,330,321,310,301,          ! K=1/2-
     & 550,541,530,521,510,501,
     & 770,761,750,741,730,721,710,701,
     & 990,981,970,961,950,941,930,921,910,901,  
     & 101,321,312,301,541,532,521,512,501,         ! K=3/2-
     & 761,752,741,732,721,712,701,
     & 981,972,961,952,941,932,921,912,901,
     & 0,0,0,0,0,          
     & 312,303,532,523,512,503,                     ! K=5/2-
     & 752,743,732,723,712,703,
     & 972,963,952,943,932,923,912,903,
     & 0,0,0,0,0,0,0,0,0,0,
     & 303,523,514,503,743,734,723,714,703,         ! K=7/2-
     & 963,954,943,934,923,914,903,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,               
     & 514,505,734,725,714,705,                     ! K=9/2-
     & 954,945,934,925,914,905,          
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,    
     & 505,725,716,705,945,936,925,916,905,         ! K=11/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,
     & 716,707,936,927,916,907,                     ! K=13/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,
     & 707,927,918,907,                             ! K=15/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,
     & 918,909,                                     ! K=17/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,
     & 909,                                         ! K=19/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0/     
c      
c     oblate, odd parity: 
      data niobod/101,110,301,310,321,330,          ! K=1/2-
     & 501,510,521,530,541,550,
     & 701,710,721,730,741,750,761,770,
     & 901,910,921,930,941,950,961,970,981,990,  
     & 101,301,312,321,501,512,521,532,541,         ! K=3/2-
     & 701,712,721,732,741,752,761,
     & 901,912,921,932,941,952,961,972,981,
     & 0,0,0,0,0,           
     & 303,312,503,512,523,532,                     ! K=5/2-
     & 703,712,723,732,743,752,
     & 903,912,923,932,943,952,963,972,
     & 0,0,0,0,0,0,0,0,0,0,
     & 303,503,514,523,703,714,723,734,743,         ! K=7/2-
     & 903,914,923,934,943,954,963,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,               
     & 505,514,705,714,725,734,                     ! K=9/2-
     & 905,914,925,934,945,954,          
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,    
     & 505,705,716,725,905,916,925,936,945,         ! K=11/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,
     & 707,716,907,916,927,936,                     ! K=13/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,
     & 707,907,918,927,                             ! K=15/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,
     & 909,918,                                     ! K=17/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,
     & 909,                                         ! K=19/2-
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0,0,
     & 0,0,0,0,0,0,0,0,0/   
c------write single-particle potential:     
       open(81,status='unknown',form='formatted',
     &     file='sppoten.'//outnam)
       do iz=0,nzd
         zdist=dz*iz
         do ir=0,nrd
           rdist=dr*ir
       write(81,'(3e13.5)') zdist,rdist,pot(iz,ir,1) 
         enddo
       enddo  
c------end write single-particle potential
       open(80,status='unknown',form='unformatted',
     &     file='nillson.'//outnam) 
       if (ireg(1).GT.0) then
       write(80) atnum,apnum,annum,beta2,beta2p,beta2n,
     &  beta4,beta4p,beta4n,nst(1),nst(2),nzd,nrd,dz,dr
c     &  gap(1),gap(2),fer(1),fer(2) 
       endif 
       do 1 np=1,2      
c------ protons (np=1) and neutrons (np=2)       
        if(np.eq.1) then 
           write(7,'(10x,a)') 'PROTONS:'
	else
	   write(7,'(10x,a)') 'NEUTRONS:' 
	endif
	 write(7,'(10x,a)') 'Reinhards quasiordering'           
       init=1+nst(1)*(np-1)
       ifin=nst(1)+nst(2)*(np-1)
       do i=init,ifin
          iord2(i)=iorder(i)
          ist=iord2(i)+nst(1)*(np-1)
          spen2(i)=spenrg(ist)
           write(7,'(2x,2i3,2x,e11.5)') i,iord2(i),spen2(i)
        enddo    
c----------------------------------------
c        Ordering of Reinhard  is approximate. 
c        So, before the asignment, we repeat the ordering.
c----------------------------------------    
c********************        
        call order(init,ifin) 
c******************** 
c       write(7,'(2x,a)')  'after ordering 1:'
c       do i=init,ifin
c           write(7,'(2x,i3,2x,e11.5)') iord2(i),spen2(i)
c        enddo          
       write(7,'(2x,a)') '  i  ist [Nn_zl]  k   occ    energy' 	  
c--------assignment of asypt. Nillson quantum numbers:  
       do i1=1,2
c         i1 means parity but not isospin
          do i2=1,10
c            i2 means number of K for a given parity.           
            icount(i1,i2)=0              
          enddo 
       enddo        
       do 2 i=init,ifin      
         ist=iord2(i)+nst(1)*(np-1)
	 lam=mpsi(ist) 
	 ipar=npar(ist)*nrpar(ist)
	 if(ipar.eq.1) then
	   i1=1
	   i2=lam+1  
c   i2=mpsi+1=(mpsi+1/2)+1/2=K+1/2
c   K   i2
c   1/2  1
c   3/2  2
c   ...  ...	   
	   icount(i1,i2)=icount(i1,i2)+1
	   i3=icount(i1,i2) 
c   i3 number of [Nn_z|Lambda] for a given K^{\pi}  	   
	   if(beta2.ge.0.0) then 
	     nqn(ist)=niprev(i3,i2)    	     
	   else
	     nqn(ist)=niobev(i3,i2)
	   endif  
	 else
	   i1=2
	   i2=lam+1  
	   icount(i1,i2)=icount(i1,i2)+1
	   i3=icount(i1,i2) 
	   if(beta2.ge.0.0) then 
	     nqn(ist)=niprod(i3,i2)    
	   else
	     nqn(ist)=niobod(i3,i2)
	   endif  	   
	 endif  
	 ak(ist)=mpsi(ist)+0.5
	 write(7,'(2x,3(i3,2x),f5.1,f8.5,2x,e11.5)') i,
     &	      ist,nqn(ist),ak(ist),vsqf(ist),spenrg(ist)        
         if(.not.(np.eq.1.and.vsqf(ist).gt.0d0.and.ireg(1).eq.1)) 
     &                                       goto 2
c----write w.f.             
         write(80) (((psi(iz,ir,is,ist),
     &               iz=0,nzd),ir=0,nrd),is=1,nsd)
        write(80) spenrg(ist),mpsi(ist),npar(ist),nrpar(ist),
     &	          ipar,icharg(ist),vsqf(ist),iorder(ist),
     &            nqn(ist),spn(ist)
c----end write w.f.
  2      continue  
  1      continue 
      close(80) 
      close(81)      
      return 
      end 
c***************************************
c****************************************
      subroutine order(imin,imax)
      include 'paramr.inc'
      include 'regime.inc'
c      common/ord/iord2(npsid),spen2(npsid)   
      if(imin.ge.imax) goto 8
      max1=imax-1
      min1=imin
c       write(7,'(2x,3i5)'), imin,imax,npsid
    1 ip=1
      j1=ip
      j=j1
      do 7 i=min1,max1
        if(spen2(i).le.spen2(i+1))
     *                                  GOTO    7
        x1=spen2(i)
        spen2(i)=spen2(i+1)
        spen2(i+1)=x1
        k=iord2(i)
        iord2(i)=iord2(i+1)
        iord2(i+1)=k
        if(.not.(ip.eq.1))
     *                         GOTO    6
        j=i-1
        ip=0
    6   j1=i
    7 CONTINUE
      max1=j1
      min1=max(j,imin)
      if(max1.gt.min1)
     *                         GOTO    1
    8 continue
      return
      end
