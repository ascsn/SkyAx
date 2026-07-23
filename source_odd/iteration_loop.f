c---------- iteration loop
c
      subroutine iteration_loop(outnam,x0sav,qact)



      include 'paramr.inc'
*
!      character*80 title
      character*20 set
      character*3  outnam
c
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

!     HFB driving parameters
      REAL(8) :: pairenhance=1.0
      INTEGER :: iterhfb=600,iternat=0,modnat=100000
      LOGICAL :: thfbparams=.FALSE.
      REAL(8) :: pairfcsav(2),ang(2)
c
cMPI      include 'mpir.inc'
c
      data iallop/0/                    ! protocol of opening files
   
      logical tafter                    ! switch to "afterburn"
      data tafter/.FALSE./
*
      LOGICAL :: forceblo,thfbsav
      LOGICAL :: tfirst=.TRUE.               ! print headers in first call
      LOGICAL :: tnobcs=.FALSE.     ! suppress BCS and pairenhance after first beta
      LOGICAL,PARAMETER :: tprecis=.TRUE.
      INTEGER :: iterbcs=300 


      INQUIRE(FILE='hfbparams',EXIST=thfbparams)
      v2limh=1D-6; delmix=0.8D0            ! default for flexible damping in 'lagpar'
                                          ! switch off with 'delmix=0'&'v2limh=1'
      IF(thfbparams) THEN
        OPEN(111,file='hfbparams')
        READ(111,*) iterhfb,iternat,modnat,pairenhance,iterbcs,
     &              v2limh,delmix,tnobcs
        WRITE(96,'(1x/1x)')
        WRITE(96,'(a,4i5,f7.2,2(1pg13.5),L1)')
     &   '# iterhfb,iternat,modnat,iterbcs,pairenhance,v2limh,delmix=',
     &   iterhfb,iternat,modnat,iterbcs,pairenhance,v2limh,delmix,tnobcs
        WRITE(6,'(a,4i5,f7.2,2(1pg13.5),L1)')
     &   '# iterhfb,iternat,modnat,iterbcs,pairenhance,v2limh,delmix=',
     &   iterhfb,iternat,modnat,iterbcs,pairenhance,v2limh,delmix,tnobcs
        CLOSE(111)
      END IF
      WRITE(7,'(a,4i5,f7.2,2(1pg13.5),L1)')
     &   '# iterhfb,iternat,modnat,iterbcs,pairenhance,v2limh,delmix=',
     &   iterhfb,iternat,modnat,iterbcs,pairenhance,v2limh,delmix,tnobcs
!      thfb = (iterhfb > 0)

!      IF(.NOT.tfirst .AND. SUM(kblock)==0) iterhfb=0
!      IF(.NOT.tfirst) iterhfb=0
      IF(thfb) THEN
        fxisav = fxi
        h1xisav = h1xi
        hgradxisav = hgradxi
        fxipsav = fxip
        h1xipsav = h1xip
        hgradxipsav = hgradxip
        fxinsav = fxin
        h1xinsav = h1xin
        hgradxinsav = hgradxin
        pairfcsav = pairfc
      END IF

      thfbsav = thfb
      IF(tnobcs.AND..NOT.tfirst) pairenhance=0D0
      suv_lambda=0D0
      do 1000 iter=0,itmax


        IF(iter>iterbcs .OR. (tnobcs.AND..NOT.tfirst)) THEN
          thfb = thfbsav
          ipair=ipairin
        ELSE
          thfb = .FALSE.
        END IF
        IF(thfbsav) THEN
          IF(iter<iterhfb) THEN
            pairfacact = 1D0+pairenhance*MAX(0,iterhfb-iter)/iterhfb
            fxi = fxisav*pairfacact
            h1xi = h1xisav*pairfacact
            hgradxi = hgradxisav*pairfacact
            fxip = fxipsav*pairfacact
            h1xip = h1xipsav*pairfacact
            hgradxip = hgradxipsav*pairfacact
            fxin = fxinsav*pairfacact
            h1xin = h1xinsav*pairfacact
            hgradxin = hgradxinsav*pairfacact
            pairfc = pairfcsav*pairfacact
          ELSE
            fxi = fxisav
            h1xi = h1xisav
            hgradxi = hgradxisav
            fxip = fxipsav
            h1xip = h1xipsav
            hgradxip = hgradxipsav
            fxin = fxinsav
            h1xin = h1xinsav
            hgradxin = hgradxinsav
            pairfc = pairfcsav
          END IF
        END IF
        IF(ioptsk==13) WRITE(7,'(a,2i5,3(1pg13.5),L2)') 
     &       'iter,ipair,fxi,fxip,fxin,thfb:',
     &        iter,ipair,fxi,fxip,fxin,thfb

!          write(13,*) betmin,betmax,betaw,qact,qmin,qmax,qw   


c
         oldoct=beta3
c
         if(iter.ne.0) then
         IF(thfb) THEN
           variance=SUM(symmiso) !devmf2 !SUM(symmiso)
         ELSE
           variance=devh2
         END IF
      write(6,'(a,i5,f9.2,2f8.3,3f5.1,1pg12.4)')
     &    ' It,E,b2,b3,D^2E,dE/dIt,dQ/dIt,zcm=',
     &      iter,energy,
     &      beta2,beta3,
     &      log10(max(variance,smal)),log10(max(econ,smal)),
     &      log10(max(qcon,smal)),
     &      zcm
c     write(6,'(a,1pg12.4)') ' zcm=',zcm
        IF(isym.eq.0) THEN
          write(7,'(a,i5,f9.2,2f8.3,3f6.1,1pg12.4)')
     &           ' It,E,b2,b3,D^2E,dE/dIt,dQ/dIt=',
     &         iter,energy,
     &         beta2,beta3,
     &         variance,log10(max(econ,smal)),
     &         log10(max(qcon,smal)),
     &         zcm
        ELSE
          write(7,'(a,i5,f9.2,2f8.3,3f6.1,2f9.3)')
     &           ' It,E,b2,b3,D^2E,dE/dIt,dQ/dIt=',
     &         iter,energy,
     &         beta2,beta3,
     &         variance,log10(max(econ,smal)),
     &         log10(max(qcon,smal)),
     &         ggg !epair
        END IF
            CALL flush(7)
         endif
c      write(7,'(a,1pg12.4)') ' zcm=',zcm
c
c     set addition-weights for mixing density and potential,
c     first iteration needs a refresh
c
      if(iter.le.3) then
        adyrho = one
        adypot = one
      else
        adyrho = addrho
        adypot = addpot
c        adyrho = max(addrho,one-exp(-iter*one/30.0))
      endif
ctest      if(iter.gt.100) addocc=0.01
ctest      if(iter.gt.100) qeps = 0.94 * qeps
ctest        if(iter.eq.101) then
ctest          e0inv = zero
ctest          x0inv = 0.0005
ctest        endif
        if(iter.le.15) then
          x0inv = x0sav  !  0.5*x0sav
        else
          x0inv = x0sav
        endif
c
        IF(iter>100 .AND. MOD(iter,1)==0) THEN
          swdiag = 1D0
        ELSE
          swdiag = 0D0
        END IF
        call sstep(iter)
        call sinfo(iter,outnam)
        FLUSH(6)
c
        dnewoct=beta3
c
c here we have to distinguish between odd-odd, even-odd and odd-even


        write(678,*)'compblock:',mblock(1),kblock(1),mblock(2),kblock(2)


        if(iter.ge.(0.8*itmax))then
        forceblo = 1
        else
        forceblo = 0   
        endif   


        if(iter.gt.itpini+10)then                             ! blocking can't start before pairing is switched on 
        if((forceblo.or.econ.lt.10.0*errenr )
     &   .and. kblock(1).ne.0    ! blocking of neutrons and protons 
     &   .and. kblock(2).ne.0
     &   .and. istblo(1).le.0 
     &   .and. istblo(2).le.0)then
              call comp_block(1)
              call comp_block(2)
        elseif((forceblo.or.econ.lt.10.0*errenr) 
     &      .and. kblock(1).ne.0 ! blocking of protons
     &   .and. kblock(2).eq.0
     &   .and. istblo(1).le.0)then 
c     &   .and. istblo(2).le.0)then
           call comp_block(1)
        elseif((forceblo.or.econ.lt.10.0*errenr) 
     &      .and. kblock(1).eq.0 ! blocking of neutrons 
     &   .and. kblock(2).ne.0
c     &   .and. istblo(1).le.0 
     &   .and. istblo(2).le.0)then
           call comp_block(2)
        endif    
        endif 




        itert  = itert+1

        

        do iq =1,2
         if(istblo(iq).ne.0) then   
!            wate(istblo(iq)) = one  ! for the blocked state the degeneracy wate is one
         endif                      ! -> contribution to epair is 0
        enddo                       ! -> degeneracy for calculation of the energy is 0 
                                    !    for this state
       

        
!      WRITE(7,'(a,3i5,2L3)') 'iter,iterhfb,iterbcs,thfb,thfbsav=',
!     &     iter,iterhfb,iterbcs,thfb,thfbsav
c        if(iter.gt.2 .and. swicon.gt.oh) then
c         if(iter.gt.itpini+3 .and. swicon.gt.oh) then 
        if(iter.gt.itpini+13 .and. swicon.gt.oh) then             
!          if (((econ.le.errenr .and. variance.lt.errvar).or.
!     &         (econ.le.1E-10 .and. variance.lt.10.0D0*errvar)).and.
          if ((econ.le.errenr .and. variance.lt.errvar).and.
     &          qcon.le.qerr.and.
c     &          iter.gt.50.and.
     &          abs(zcm).lt.1E-3.and. 
     &          (iter>iterhfb+1000 .OR. .NOT.thfbsav) .AND.
     &          abs(oldoct-dnewoct).lt.1E-2) goto 900
        else if ((swicon.lt.oh).and.(iter.gt.itpini+13)) then  ! blocking
          if(econ.lt.errenr .and.
     &          (iter>iterhfb+1000 .OR. .NOT.thfbsav) .AND.
     &       (variance.lt.errvar .OR. thfb)) goto 900
        endif
c

        if(iter.gt.itcons .and. iflmb0.eq.1) then
          qlamb  = zero
          swicon = zero
        endif
c
 1000 continue

      do iq = 1,2
         IF(istblo(iq)>0) write(11,*)'wate(1/2)',iq,wate(istblo(iq))
      enddo

*
***********************************************************************
*                                                                     *
*     end of static iteration                                         *
*                                                                     *
***********************************************************************
*
  900 continue
c
c     one more iteration to put observables up to date
c
      if(tafter) then
        ipold  = iprint
        iprint = 1
        iplold = iplot
        itesto = itest
        iplot  = 1
        itest  = 0
        econs  = econ
        qcons  = qcon
        devh2s = devh2
        devmf2s = devmf2
        adyrho = 1.0
        adypot = 1.0
        call sstep(iter)
        call zpe_cm()
        call sinfo(itmax+1,outnam)                       ! "afterburn"
        iprint = ipold
        iplot  = iplold
        itest  = itesto
        econ   = econs
        qcon   = qcons
        devh2  = devh2s
        devmf2 = devmf2s
      else
         call sinfo(iprint,outnam) ! "afterburn"
      endif
c
c      iter   = iter +1
c      tmass = swicon.gt.oh .and. econ.le.20.0*errenr 
      IF(thfb) THEN
        variance=SUM(symmiso) !devmf2 !SUM(symmiso)
      ELSE
        variance=devh2
      END IF
      tmass = econ.le.20.0*errenr                                    !!!!
     &        .and. (variance.lt.20.0*errvar .OR. thfb)
     &        .and. ifcran.eq.1 .and. qcon.le.10.0*qerr
!      tmass=.FALSE.
c
      if(ifcran.eq.1) then 
ctest        write(25,'(a,f8.3)') 'LINRESP called',beta2 !!!!
        call linresp(outnam)
        call crankin()
        call crnkam()
        call zpe_cm()
        call trizpe
        if(swicon.lt.oh) call linresp(outnam)
      elseif(ifcran.eq.2) then
        call seprpa
      else
        call zpe_cm()
        call crinfo(outnam)
      endif
c
c      call wrt_spprop(outnam,0)

     



c
c     optionally print potentials on separate file
c
c      call pripot(title,outnam)
c      if(klocal.eq.1) call localize(outnam)
c
c     printout on file 'ensurf.???' is done here.
c     check first call this routine and compose a header.
c
      open(90,status='unknown',form='formatted',file='ensurf.'//outnam,
     &        position='append')
!      do ispool=1,9999
!       read(90,'(1x)',end=9999)
!      enddo
! 9999 continue
!      WRITE(90,'(a,2(1pg13.5),2i4)') 'test90:',qact,qmin,istblo
!      if (qact.eq.qmin .and. tempac.eq.tempmn) then
       if(ipair.ne.4) then          
          IF(istblo(1)+istblo(2)==0) THEN
          if(tfirst) then
          tfirst = .FALSE.
            write(90,'(a)')
     &       '#'//title(4:75),
     &       '#X:&P25b','Y:&P16PES',
     &      '#S: DIST(4) SCAL(0.8)'
            write(90,1098) title,npart(1)+npart(2),npart(1),iparam,
     &               dz,dr,nz,nr,
     &               addrho,addpot,addocc,c0,qeps,
     &               x0inv,e0inv,itsub,x0fac
 1098 format('# Plot potential energy surface for:',/,1x,a,/,
     &       '# Nmass=',i5,', Ncharge=',i5,',  force from IPARAM=',i5,/,
     &       '#  dz,dr: ',2f6.3,' nz,nr:',2i4,/,
     &       '# addrho,addpot,addocc: ',3f6.3,'  c0,qeps:',2f6.3,/,
     x       '# x0inv(fm),e0inv(mev): ',2f7.2,' itsub,x0fac:',i3,f6.3)
            WRITE(90,'(a,L3,f8.2)') 
     &           '# pairing: v**2 space, cutoff=',tpairspace,cutfac
            WRITE(90,'(a,5f7.2)')
     &           '# sum_uv constraint:',tempmn,tempmx,tempw,suv_mue
          endif  
          END IF
        DO iq = 1,2 
          IF(istblo(iq).gt.0) THEN
            WRITE(90,'(a,4i5,1pg13.5)') 
     &        '# blocking: kblock,mblock,isoblo,istblo,ediff=',
     &        kblock(iq),mblock(iq),iq,istblo(iq),ediffblock(iq)
            IF(tmass) THEN
              OPEN(76,status='unknown',form='formatted',
     &        file='zpecor.'//outnam,position='APPEND')
              WRITE(76,'(a,4i5)') 
     &          '# blocking: kblock,mblock,isoblo,istblo=',
     &          kblock(iq),mblock(iq),iq,istblo(iq),ediffblock(iq)
              CLOSE(76)
              OPEN(75,status='unknown',form='formatted',
     &             file='crankn.'//outnam,position='APPEND')
              WRITE(75,'(a,4i5)') 
     &          '# blocking: kblock,mblock,isoblo,istblo=',
     &          kblock(iq),mblock(iq),iq,istblo(iq),ediffblock(iq)
              CLOSE(75)
            END IF
          END IF
        END DO
          
          IF(istblo(1)+istblo(2)==0) THEN
            if(isym.eq.0) then
              write(90,'(3a)')
     &       '#  bet1 bet1p bet1n  bet2 bet2p bet2n  bet3 bet3p bet3n',
     &       '  bet4 bet4p bet4n   energy     epairp  epairn    iter',
     &       '   rms    rmsp   rmsn'
!              write(90,'(a)') 'H: N N N X  N N N N N N N N Y1,L ' 
            else
              write(90,'(3a)')
     &       '#  bet2  bet2p  bet2n',
     &       '  bet4   bet4p  bet4n   energy       epairp        ',
     &       'epairn    iter  rms    rmsp   rmsn  jz_p jz_n   j'
!              write(90,'(a)') 'H: X  N N N N N Y1,L ' 
            endif
          END IF
	else
          write(90,1092) title,npart(1)+npart(2),npart(1),iparam,
     &               tempmn,tempmx,tempw,
     &               dz,dr,nz,nr,
     &               addrho,addpot,addocc,c0,qeps,
     &               x0inv,e0inv,itsub,x0fac
 1092 format('# Plot potential energy surface for:',/,1x,a,/,
     &       '# Nmass=',i5,', Ncharge=',i5,',  force from IPARAM=',i5,/,
     &       '# tempmn,tempmx,tempw=',3f6.3,/,
     &       '#  dz,dr: ',2f6.3,' nz,nr:',2i4,/,
     &       '# addrho,addpot,addocc: ',3f6.3,'  c0,qeps:',2f6.3,/,
     x       '# x0inv(fm),e0inv(mev): ',2f7.2,' itsub,x0fac:',i3,f6.3,/,
     &       '#   T    beta2   beta2p    beta4   beta4p   ',
     &       '#   energy     efree         T*S   level-d.')
	endif
!      endif
      and=0D0
      if(ipair.ne.4) then
        if(isoblo.ne.0) then
      	  call calangularmomentum(istblo(isoblo),ang(isoblo))
          WRITE(*,*) 'BLOCK: isoblo,ist,ang=',isoblo,istblo(isoblo),
     &                ang(isoblo)
        ELSE IF(kblock(1).NE.0 .AND. kblock(2).NE.0) THEN
      	  call calangularmomentum(istblo(1),ang(1))
      	  call calangularmomentum(istblo(2),ang(2))
        endif
        if(isym.eq.0) then
          write(90,'(3f6.2,f7.3,2f6.2,f7.3,5f6.2,f10.3,4f7.1,i5,4f8.4,
     &               4f5.2,3f10.4)')
     &      beta1,beta1p,beta1n,beta2,beta2p,beta2n,
     &      beta3,beta3p,beta3n,beta4,beta4p,beta4n,
     &      energy,epair,fer,iter,rms,rmsp,rmsn,rmsc,
     &      0.5*kblock(1),ang(1),0.5*kblock(2),ang(2),r4mom,qp,
     &      sum_uv(2)
        else
         IF(tprecis) THEN
           write(90,'(3f16.12,3(1pg20.12),0pf15.8,2f8.2,2(1pg20.12),
     &               i5,3(0pf15.11),4f7.2,3(1pg20.12),2(0pf10.4))')
     &      beta2,beta2p,beta2n,beta4,beta4p,beta4n,
     &      energy,epair,gap,iter,rms,rmsp,rmsn,!fer,iter,rms,rmsp,rmsn,
     &      0.5*kblock(1),ang(1),0.5*kblock(2),ang(2),rmsc,r4mom,qp,
     &      sum_uv(2),suv_wanted(2)
         ELSE
           write(90,'(f7.3,5f7.3,f10.3,4f7.1,i5,3f8.4,4f7.2,
     &              2f7.4,3f8.2)')
     &      beta2,beta2p,beta2n,beta4,beta4p,beta4n,
     &      energy,epair,gap,iter,rms,rmsp,rmsn,!fer,iter,rms,rmsp,rmsn,
     &      0.5*kblock(1),ang(1),0.5*kblock(2),ang(2),rmsc,r4mom,qp,
     &      sum_uv(2),suv_wanted(2)
         ENDIF
        endif
        if(abs(apnum-npart(1)).gt.1.0E-5 
     &     .or. abs(annum-npart(2)).gt.1.0E-5 )
     &    write(90,'(35x,a,2(1pg14.6))')   
     &       'part.nr.=',apnum,annum
      else
	efree  = energy  
	energy = efree-etpair  
        write(90,'(f6.2,4f8.3,3f12.3,f7.2)')
     &     tempac,beta2,beta2p,beta4,beta4p,energy,efree,-etpair,
     &     -3.0*etpair/(tempac*pi)**2
      endif
      close(90)

      CALL pri_2D(outnam)
c      if(Min(NX3,NY3,NZ3).gt.1) stop ' only one run with RESTART'


      return

      end
