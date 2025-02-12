c-----sinfo -----------------------------------------------------------
c
      subroutine sinfo(iter,outnam)
*
***********************************************************************
*                                                                     *
*     print of densities if iter matches print modulus.               *
*     evaluate energies, mass number , deformation, radii, etc.       *
*     prints those quantities if iter matches print modulus.          *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
*
      character*3  outnam
c
      character*4 itype(2,2)
c

      dimension  work(0:nzd,0:nrd)
!      dimension  angular(1:nst(1)+nst(2))
      dimension  parity(1:nst(1)+nst(2))
      REAL(8) :: totjz(2),totlz(2),totsz(2)
      REAL(8),SAVE :: econold,econnew=1D3
!      REAL(8) :: sumuv(kqm),sumuv2(kqm),sumuve(kqm)
c
c      data rho0pr/0.16/               ! preliminary standard for surface pairing
      data half,xsmall/0.5,1.0E-20/



c     some auxiliary fields for printing single particle properties
c     ?? note: one may change the horizontal printout to a vertical
c              printout which may look nicer and requires less
c              bookkeeping
c
c
      data itype /'prot','on  ','neut','ron '/
      data iq2old/0/               ! switch r_rms <--> R_dms in beta_i
      LOGICAL,PARAMETER :: testobsti=.FALSE.
!
!      WRITE(6,'(a)')'rho(0:4,0:4):'
!      DO ir=0,4
!        WRITE(6,'(i2,5(1pg13.5))') ir,rho(0:4,ir,1)
!      END DO
!      CALL test_intpol2D(rho(:,:,1))
!      STOP
!
!     test Q
!
!      nsttot=nst(1)+nst(2)
!      WRITE(6,'(a,5(1pg13.5))') 'enter info:   PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wate(1:nsttot)*wtpair(1:nsttot))
!      DO i=1,nsttot
!        WRITE(6,'(i4,3(1pg13.5))') i,wate(i),wtpair(i),vsqf(i)
!      END DO
      IF(testobsti) THEN
        CALL fdenti(psi,psip)
        CALL obsti(psip,qtil,qtilc,qtil2,qtilc2)
!        wrk1=psip(:,:,1)+psip(:,:,2)
!        CALL obs(wrk1,atnum,rms,d,q,o,h,qsq,zcm)
        WRITE(65,'(a,i6,4(1pg12.4))')
     &   ' INFO-0: iter,q,qc=',iter,qtil,qtilc
!     &   ' INFO-0: iter,A,rms,q,qc=',iter,atnum,rms,qtil,qtilc
      END IF


*
***********************************************************************
*                                                                     *
*         block to compute the shell energy, the field energy         *
*         and the rearrangement energy                                *
*                                                                     *
***********************************************************************
*
c      if(mod(iter,50).eq.0) call testls
c                                                constraint energy
      do ir=0,nr
      do iz=0,nz
        wrk1(iz,ir) = (rho(iz,ir,1)+rho(iz,ir,2))*qfac(iz,ir)
      enddo
      enddo
      equad  = qlamb*f1int(wrk1)
c
c     the zero point energy in a rough harmonic oscillator approximation
c
cold  ezp    = 30.75*(npart(1)+npart(2))**(-edr)
c
c     the zero point energy as estimated from an average over spherical results
c     ?? !! note: this recipe is different from the relativistic axial code.
c
      anumb  = npart(1)+npart(2)
      if(ifzpe>0) then
!        IF(iter.LE.itmax .OR. istsym>0) THEN
        IF(iter.LE.itmax+10 .OR. istsym>0) THEN  
          ezp    = 17.3/anumb**0.2
        ELSE
!         replace old estimate by correct Z.P.E. ! uncertain for odd nuclei
          CALL zpe_cm()
          ezp = (widdz+widdxy)*0.5D0*(h2m(1)+h2m(2))/(npart(1)+npart(2))
        END IF
      else
       ezp    = 0.0D0
      endif
*
***********************************************************************
*                                                                     *
*     pairing energy:                                                 *
*        ipair=0 for constant gap model                               *
*        ipair=1 for constant force model                             *
*        ipair=2 for constant force model  (BCS)                      *
*        ipair=3 for Lipkin-Nogami                                    *
*        ipair=4 free energy for thermal CHF                          *
*                                                                     *
***********************************************************************
*
      IF(ioptsk.NE.13) THEN
      if(ipair.eq.0 .or. ipair.eq.1) then
         epair(1)=-gap(1)*gap(1)/max(ggg(1),smal)
         epair(2)=-gap(2)*gap(2)/max(ggg(2),smal)
      elseif(ipair.eq.4) then                         ! compute entropy
	 epair(1) = zero
	 epair(2) = zero
	 do ist=1,lst
	   iq   = icharg(ist)
	   occ  = max(vsqf(ist),smal)
	   occd = max(one-occ,smal)
           epair(iq) = wate(ist)*(occ*log(occ)+occd*log(occd))
     &                +epair(iq)
         enddo
         epair(1) = gap(1)*epair(1)                 ! -T*S
         epair(2) = gap(2)*epair(2)
      end if
      END IF
c      etpair=epair(1)+epair(2)

cJE-20-7-07
ccccccccccccccccccccccnew calculation of the pairing energycccccccccccccccccccccccc
c
c
c

c
c--------------------------------------------------------------------
c
c


!      WRITE(6,'(a,5(1pg13.5))') 'before paircr:PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wtpair(1:nsttot))
      IF(ioptsk.NE.13) THEN
c
c     zero field before accumulation
c
      do iq=1,2
      do ir=0,nr
      do iz=0,nz
       paircr(iz,ir,iq)  = zero
      enddo
      enddo
      enddo

c     accumulate new pair-density
cc
cc      smallp = 0.000001
      smallp = 1D-40
      do jst=1,nst(1)+nst(2)
        iq    = icharg(jst)
        IF(jst==istblo(iq) .OR. jst==istcblock(iq)) CYCLE
!        IF(jst==istblo(iq)) CYCLE
        v2    = vsqf(jst)
        uv    = wtpair(jst)*sqrt(max(v2-v2*v2,smallp))*0.5

        if(isTSym.eq.0) then
         call sdenti(psi(0,0,1,jst),paircr(0,0,iq),uv,jst)
        else
         if(mpsi(jst).lt.0)then
          uv = (-1.0)*uv
         endif
          call cpairpot(psi(0,0,1,jst),psi(0,0,1,itpartn(jst)),
     &     paircr(0,0,iq),uv,jst)
        endif

      enddo

cc
cc     now multiply with strength to obtain local pair-potential
cc
      v0prot = -pairfc(1)
      v0neut = -pairfc(2)
      if(ipair.le.6 .or. ipair.eq.9) then
        do ir=0,nr
        do iz=0,nz
          pairpt(iz,ir,1)  = v0prot * paircr(iz,ir,1)
          pairpt(iz,ir,2)  = v0neut * paircr(iz,ir,2)
        enddo
        enddo
      else
        do ir=0,nr
        do iz=0,nz
          rhofac = one-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr
          pairpt(iz,ir,1)  = v0prot * rhofac * paircr(iz,ir,1)
         pairpt(iz,ir,2)  = v0neut * rhofac * paircr(iz,ir,2)
        enddo
        enddo
      endif
      END IF

cc
cc     compute pairing energy 
cc
c      do ir=0,nr
c      do iz=0,nz     
c      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)
c      enddo
c      enddo
c      
c      epairtestp=f1int(work)  !siehe oben
c      
c      do ir=0,nr
c      do iz=0,nz     
c      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)
c      enddo
c      enddo
c      
c      epairtestn=f1int(work)  !siehe oben
c
ccccccccccccccccccccccc Nachbesserung DDDI 18.02.08 ccccccccccccc
c
c     compute pairing energy 
c

      IF(ioptsk.NE.13) THEN
      if(ipair.le.6 .or. ipair.eq.9) then
      do ir=0,nr
      do iz=0,nz     
      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)
      enddo
      enddo
      
      epairtestp=f1int(work)  !0.25 oder 0.5,was ist v0 ?
      
      do ir=0,nr
      do iz=0,nz     
      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)
      enddo
      enddo
      
      epairtestn = f1int(work)  !0.25 oder 0.5,was ist v0 ?

      

      else
      
      

      do ir=0,nr
      do iz=0,nz     
      rhofac = one-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr
      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)* rhofac
      enddo
      enddo
      
      epairtestp=f1int(work)  !0.25 oder 0.5,was ist v0 ?
      
      do ir=0,nr
      do iz=0,nz     
      rhofac = one-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr
      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)* rhofac
      enddo
      enddo
      
      epairtestn = f1int(work)  !0.25 oder 0.5,was ist v0 ?
      
      endif
      


      etpair= epairtestp+ epairtestn

      epair(1) = epairtestp
      epair(2) = epairtestn
      ELSE
        epairtestp = 0D0
        epairtestn = 0D0
        do na=1,nst(1)+nst(2)
          iq    = icharg(na)
          IF(tpairspace) THEN
             deltaa = deltaf(na)
           ELSE
             deltaa = deltaf(na)*wtpair(na)
           END IF
	   v2  = vsqf(na)
           phact  = 0.5D0*wate(na)   ! for axial
           IF(tpairspace) THEN
             vol    = phact*sqrt(max(v2-v2*v2,xsmall))*wtpair(na)
           ELSE
             vol    = phact*sqrt(max(v2-v2*v2,xsmall))
           END IF
           IF(na==istblo(iq) .OR. na==istcblock(iq)) vol=0D0
           sumuv  = vol        + sumuv
           IF(iq==2) THEN           
             epairtestn = epairtestn + vol*deltaa 
           ELSE
             epairtestp = epairtestp + vol*deltaa 
           END IF
        enddo
!        WRITE(7,'(a,4(1pg13.5))') 'EPAIR:',epairtestp,epairtestn,epair
        etpair= epairtestp+ epairtestn
      END IF

c      
c      etpair=epair(1)+epair(2)
c
*
***********************************************************************
*                                                                     *
*     finally the total energy                                        *
*                                                                     *
***********************************************************************
*
!      WRITE(6,'(a,5(1pg13.5))') 'before E_tot: PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wtpair(1:nsttot))
      if((itest.ge.1 .or. iparam.eq.3) .AND. ioptsk.NE.13) then
       call energ2()         ! energy using rearrangement energy
       energr = energr - ezp + etpair
      endif

c
      IF(ioptsk==13) THEN
        do ir=0,nr
        do iz=0,nz
           tp      = tau(iz,ir,1)
           tn      = tau(iz,ir,2)
           wrk2(iz,ir) = h2meff(1)*tp+h2meff(2)*tn
        enddo
        enddo
        ekin   = f1int(wrk2)
!        CALL  fy_hfbpot()
        etot1  = ekin+energyFy+ecoul/2D0
        energy = ekin+energyFy+ecoul/2D0-ezp+SUM(epair)
        IF(itest>0) WRITE(7,*) '2. energyFy=',etot1,epair,ezp
      ELSE
        call energ()           ! energy directly from density
        etot1  = energy
        energy = energy-ezp+etpair
      END IF
c
      if(itest.ge.1) then
       write(6,'(a,i5,a,2f12.4)')
     &  ' iter=',iter,
     &  ', energies: from rearr. vs. from densities=',energr,energy
       write(7,'(a,2f12.4/a,4(1pg12.5))')
     &  '  energies: from rearr. vs. from densities=',energr,energy,
     &  '  ekin,epot,etot=',ekin,epot,etot1,energy
      endif
c
c
      varmax=-1.0d5
      varmin= 1.0d5
      do ist=1,lst
        varmax = max(varmax,spevar(ist))
        varmin = min(varmin,spevar(ist))
      enddo
      do ist=1,lst
        speold(ist) = spenrg(ist)
      enddo
c
c
c     the spatial observables derived from the densities.
c
c       a(t/n/p)num ---> mass number for the nucleus(t),
c                        for protons(p) and neutrons(n)
c
c       rms( /n/p)  ---> rms radius for the nucleus( ),
c                        for protons(p) and neutrons(n)
c
c         q( /n/p)  ---> quadrupole moment for the nucleus( ),
c                        for protons(p) and neutrons(n)
c
c         h( /n/p)  ---> hexadecupole moment for the nucleus( ),
c                        for protons(p) and neutrons(n)
c
c     are computed in the subroutine obs and retruned via common.
c
c
c      zcm = zcmobs(rho)
      IF(ttestconstr) THEN
        call obsti(rho,qtil,qtilc,qtil2,qtilc2) ! new quadrupole momenta
        write(65,'(a,i6,4(1pg12.4))')
     &   ' INFO: iter,q,qc,q2,qc2=',iter,qtil,qtilc,qtil2,qtilc2
      END IF
      call obs(rho(0,0,1),apnum,rmsp,dp,qp,op,hp,qsqp,zcm)
      call obs(rho(0,0,2),annum,rmsn,dn,qn,on,hn,qsqn,zcm)
      atnum = apnum+annum
      rms   = sqrt((rmsp*rmsp*apnum+rmsn*rmsn*annum)/atnum)
      if(itest.gt.0) write(6,'(a,7(1pg12.4))')
     &   'pr: Z,rms,d,q,o,h=',apnum,rmsp,dp,qp,op,hp,
     &     SUM(vsqf(1:nst(1))*wate(1:nst(1)))
      if(itest.gt.0) write(6,'(a,7(1pg12.4))')
     &   'ne: N,rms,d,q,o,h=',annum,rmsn,dn,qn,on,hn,
     &    SUM(vsqf(nst(1)+1:nst(1)+nst(2))*wate(nst(1)+1:nst(1)+nst(2)))
      d     = dp+dn
      q     = qp+qn
      o     = op+on
      h     = hp+hn
      qsq   = qsqp+qsqn
ctest      if (swicon.gt.oh.and.iter.le.itcons) then
ctest        call quadc
ctest      endif
c
c     recouple to the dimensionless deformations
c     !! note: this uses only the r.m.s. radius, not the <r^l> moments !!
c
CR11OCT98      fac    = 0.2*sqrt(5.0*pi)
      fac     = sqrt(0.2*pi)
c      fac1   = 1.25/pi
      if(iq2old.eq.1) then
        rmsb   = rms
        rmsbp  = rmsp
        rmsbn  = rmsn
      else
        rmsb   = 0.93*atnum**0.3333333                ! 0.93 = sqrt(0.6)*1.2
        rmsbp  = rmsb
        rmsbn  = rmsb
      endif
      beta2  = fac*q /max(rmsb*rmsb*atnum,smal)
      beta2p = fac*qp/max(rmsbp*rmsbp*apnum,smal)
      beta2n = fac*qn/max(rmsbn*rmsbn*annum,smal)
      beta4  = h * 1.507964474/(max(atnum*rmsb**4,smal))    ! should be 2.5132=4*pi/5
      beta4p = hp* 1.507964474/(max(apnum*rmsbp**4,smal))   ! is 3/5 smaller
      beta4n = hn* 1.507964474/(max(annum*rmsbn**4,smal))
      if(isym.eq.0) then
        beta1  = d * 3.244622941/max(atnum*rmsb,smal)
        beta1p = dp* 3.244622941/max(apnum*rmsbp,smal)
        beta1n = dn* 3.244622941/max(annum*rmsbn,smal)
        beta3  = o * 1.946773765/max(atnum*rmsb**3,smal)
        beta3p = op* 1.946773765/max(apnum*rmsbp**3,smal)
        beta3n = on* 1.946773765/max(annum*rmsbn**3,smal)
      endif
c      beta0  = beta0 *(1.0 + fac1*beta0*beta0)
c      betap  = betap *(1.0 + fac1*betap*betap)
c      betan  = betan *(1.0 + fac1*betan*betan)
c
c      momentum of inertia for rotation about the x-axis
c
      do ir=0,nr
      r     = dr*(ir+sr)
      r2h   = r*r*0.5
        do iz=0,nz
        z     = ((iz)+sz-icm)*dz
        r2t   = z*z+r2h
          wrk1(iz,ir) = (rho(iz,ir,1)+rho(iz,ir,2))*r2t
        enddo
      enddo
cold      xinert  = f1int(wrk1)/h2meff
      xinert  = f1int(wrk1)/h2mav
c
c     collective quadrupole mass in scaling approach
c
CR11OCT98      bmassc = 2.66666667*h2meff*(2.0*rms*rms*atnum+5.0*q)   ! cartesian
c      bmassc = 8.0*h2meff*(2.0*rms*rms*atnum+q)              ! cartesian
      bmassc = 8.0*h2mav*(2.0*rms*rms*atnum+q)              ! cartesian
      bmassr = bmassc*(fac/(rms*rms*atnum))**2               ! r-quadr.
*
***********************************************************************
*                                                                     *
*        begin block to write out the observables                     *
*        econ is the energy convergence parameter                     *
*        eold is the energy form the last iteration                   *
*                                                                     *
***********************************************************************
*

!      WRITE(6,'(a,5(1pg13.5))') 'before print: PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wtpair(1:nsttot))
!       call testpr(paircr(:,:,2),0,9,0,9,'paircr')
!       call testpr(pairpt(:,:,2),0,9,0,9,'pairpt')
      IF(iter==itmax) THEN
        CALL testplot(pairpt(:,:,1),0,nz,0,nr,'pairpt_p')
        CALL testplot(pairpt(:,:,2),0,nz,0,nr,'pairpt_n')
      END IF
      if(iter.eq.0) then
       do iwrin=71,72
        write(iwrin,1010) npart(1)+npart(2),npart(1),iparam,
     &               dz,dr,nz,nr,
     &               addrho,addpot,addocc,c0,qeps,
     &               x0inv,e0inv,itsub,x0fac
 1010   format(
     &       '     '/'     '/,
     &       ' Nmass=',i5,', Ncharge=',i5,',  force from IPARAM=',i5,/,
     &       '  dz,dr: ',2f6.3,'   nz,nr=',2i5,/,
     &       ' addrho,addpot,addocc: ',3f6.3,'   c0,qeps=',2f6.3,/,
     &       ' x0inv(fm),e0inv(mev): ',2f7.2,' itsub,x0fac:',i3,f6.3)
       enddo
       if(swicon.lt.oh) then
         IF(thfb) THEN
           write(71,'(a,7x,a)')
     &       ' iter   energy    variance  difference    qtotal',
     &       '   devoc2    symmetry conditions (p/n) '
         ELSE
           write(71,'(a,7x,a)')
     &       ' iter   energy    variance  difference    qtotal',
     &       '   devoc2    devmf2     devgap '
         END IF
       else
         write(71,'(a,7x,a)')
     &       ' iter   energy    variance  difference    qtotal',
     &       '   qcorr     qlambd    qcon'
       endif
       write(72,'(a)')
     &       ' iter   q_total     q_prot.     q_damp.    hexadec.'
      else
       if(swicon.lt.oh) then
        IF(thfb) THEN
          write(71,'(i5,f11.4,6(1pg13.5))')
     &        iter,energy,devh2,econ,q,devmf2,symmiso
        ELSE
          write(71,'(i5,f11.4,6(1pg13.5))')
     &        iter,energy,devh2,econ,q,devmf2,devoc2,devgap
        END IF
       else
        IF(thfb) THEN
          write(71,'(i5,f11.4,2(1pg13.5),5x,6(1pg13.5))')
     &        iter,energy,devh2,econ,q,qcorr,qlamb,qcon,symmiso
        ELSE
          write(71,'(i5,f11.4,2(1pg13.5),5x,4(1pg13.5))')
     &        iter,energy,devh2,econ,q,qcorr,qlamb,qcon
        END IF
       endif
       write(72,'(i5,4(1pg12.3))') iter,q,qp,qcst,h
      end if
      IF(energy.ne.zero) THEN
        econold= econnew
        econnew   = abs((energy-eold)/energy)
        econ = 0.5D0*(econold+econnew)    ! average to avoid cusps
      END IF
c
      eold   = energy
c      qold   = q
      if (iter.eq.0) return
      if(iprint.eq.0) return
      if(mod(iter,iprint).ne.0) return
      ebin   = energy
      ebinpn = ebin/max(atnum,smal)
      write(7,25010) iter,ebin,ebinpn,eshell,ecoul/2D0,ecolex,erearr,
     x               ezp,epot,ekin,etpair,equad,econ
25010 format('1'/,' iteration no. = ',i10,
     x           /,' binding energy (mev) = ',
     x           f10.3,' binding energy per nucleon (mev) = ',f10.3,/,
     x             '   eshell = ',f10.3,
     x             '   ecoul  = ',f10.3,
     x             '   ecolex = ',f10.3,/,
     x             '   erearr = ',f10.3,
     x             '   ezp    = ',f10.3,
     x             '   epot   = ',f10.3,/,
     x             '   ekin   = ',f10.3,
     x             '   epair  = ',f10.3,
     x             '   equad  = ',f10.3,/,
     x             '   de/e   = ',1pg10.4)
      WRITE(7,'(a,3(1pg15.7))') ' ekin_p/n=',
     &  h2meff(1)*f1int(tau(:,:,1)),h2meff(2)*f1int(tau(:,:,2))
     &  ,h2meff(2)*f1int(tau(:,:,2))+etpair
      if(ipair.ge.0)
     &  write(7,'(3(a,2(1pg10.3)))')
     &             '   gap(1/2)=',gap,'   ggg(1/2)=',ggg,
     &             '   fer(1/2)=',fer
      if (swicon.gt.oh) then
         write(7,25700) qlamb
25700     format(/,' parameter determing the deformation:',
     x       /,' qlamb=',f12.9,/)
      endif
      WRITE(7,'(a,4(1pg13.5))') 'deviations: mf,full,gap=',
     &   devmf2,devh2,devgap
c
c        cartesian quadrupole moments are transformed
c        in spherical ones (qs,qsp,qsn).
c

         csf = sqrt(5.0/pi)/4.0
         qs  = q*csf
         qsp = qp*csf
         qsn = qn*csf
         write(7,'(/a,a)')
     &    '         part.no.  radius   dipol       quad        octu  ',
     &    '     hexadeca     quad_cart'
         write(7,'(a,2f8.3,5(1pg12.4))')
     &    '  proton:',apnum,rmsp,dp,qsp,op,hp,qp
         write(7,'(a,2f8.3,5(1pg12.4))')
     &    ' neutron:',annum,rmsn,dn,qsn,on,hn,qn
         write(7,'(a,2f8.3,5(1pg12.4))')
     &    '   total:',atnum,rms,d,qs,o,h,q
         write(7,'(/a,1pg12.4)') 
     &    '            center of mass=',zcm
         write(7,'(/a)')
     &    '            beta_1      beta_2      beta_3      beta_4  '
         write(7,'(a,4(1pg12.4))')
     &    '  proton:',beta1p,beta2p,beta3p,beta4p
         write(7,'(a,4(1pg12.4))')
     &    ' neutron:',beta1n,beta2n,beta3n,beta4n
         write(7,'(a,4(1pg12.4))')
     &    '   total:',beta1,beta2,beta3,beta4
      if(ipair.ne.-1.and.iter.gt.itpini) then
!         CALL pairobs(sumuv,sumuv2,sumuve)
         write(7,'(/11x,2a)')
     &     'e-pair      e-fermi     v2-gap  ',
     &     '   uv-gap      sumuv   '
         write(7,'(a,7(1pg12.4))')
     &     '  proton:',epair(1),fer(1),ggg(1),gap(1),sum_uv(1)!,ewin(1)
         write(7,'(a,7(1pg12.4))')
     &     ' neutron:',epair(2),fer(2),ggg(2),gap(2),sum_uv(2)!,ewin(2)
         write(7,'(a,5(1pg12.4))')
     &     '   total:',etpair
      endif
      call printenerg()
!     total angular momenta in case of B field
      IF(b_strength.NE.0D0) THEN
        totlz=0D0
        totsz=0D0
        DO ist=1,nst(1)+nst(2)
          iq = icharg(ist)
          DO is=1,2
            isi       = 3-is
            mm        = mpsi(ist)+is-1
            sigis     = mod(1D0*(is),2D0)*2D0 - 1D0
            WRITE(*,*) 'compute a.m.:=',ist,iq,is,mm,sigis,vsqf(ist)
            totlz(iq) = totlz(iq)+mm*vsqf(ist)*cmpn(is,ist)
            totsz(iq) = totsz(iq)+0.5D0*sigis*vsqf(ist)*cmpn(is,ist)
          END DO
        END DO
        totjz=totlz+totsz
        energmag = -b_str_act(1)*(totlz(1)+g_magnetic(1)*totsz(1))
     &             -b_str_act(2)*g_magnetic(2)*totsz(2)

        WRITE(7,'(a,f10.3)') 'B field: total magnetic energy=',energmag
        WRITE(7,'(a,3f8.3)') 'a.m. proton:',totlz(1),totsz(1),totjz(1)
        WRITE(7,'(a,3f8.3)') 'a.m.neutron:',totlz(2),totsz(2),totjz(2)
        WRITE(7,'(a,3f8.3)') 'a.m.  total:',SUM(totlz(:)),
     &        SUM(totsz(:)),SUM(totjz(:))
        WRITE(37,'(a,f10.3)') 'B field: total magnetic energy=',energmag
        WRITE(37,'(a,3f8.3)') 'a.m. proton:',totlz(1),totsz(1),totjz(1)
        WRITE(37,'(a,3f8.3)') 'a.m.neutron:',totlz(2),totsz(2),totjz(2)
        WRITE(37,'(a,3f8.3)') 'a.m.  total:',SUM(totlz(:)),
     &        SUM(totsz(:)),SUM(totjz(:))
      END IF
*
***********************************************************************
*                                                                     *
*        begin block to write out the single particle observables     *
*        note that the ordering in energy is not the natural ordering *
*        of the levels                                                *
*                                                                     *
***********************************************************************
*


!      WRITE(6,'(a,5(1pg13.5))') 'before j**2:  PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wtpair(1:nsttot))
      j2matrix=0D0   ! reset J**2 matrix for new sampling
      do 50000 iq=1,2
      nstm=nst(iq)
      write(7,25030) itype(1,iq),itype(2,iq),nst(iq)
25030 format(//2a4,' levels:',5x,i4,'  maximum',//)

codd     sort single particle states by energy
      call sortenerg(nstm,spenrg(iokk(iq)),iorder(iokk(iq)))
      call findconjpart(iq)

c calculate effective parity
      CALL calc_parity(parity)

c Calculate angular momentum
      do iaa=1,nstm
          ine    = iokk(iq)+iaa-1
          ia    = iorder(ine)+iokk(iq)-1
!          CALL calc_angularmomentum(ia,amom,avar)
          call calangularmomentum(ia, angular(ia))
      enddo


cende

codd
c
c     print levels
c
      nsttot=nst(1)+nst(2)
!      WRITE(6,'(a,5(1pg13.5))') 'before levels:PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wtpair(1:nsttot))
!      WRITE(7,'(a,5(1pg13.5))') 'before levels:PARTNM=',
!     &   SUM(vsqf(1:nsttot)*wtpair(1:nsttot))
      if(isTSym.ge.1)then
        write(7,'(/18x,a,9x,a/a,a)')
     &    'parities','components',
     &    '   Nr Iso   j_z total   r    z     v^2  1.   2.',
     &    '   energy partner |<i|T|j>|'
      else
        write(7,'(/18x,a,9x,a/a,a)')
     &    'parities','components',
     &    '   Nr Iso   j_z total   r    z   v^2  uv  delta',
     &    '   energy'
      endif

      if(ipair.lt.6) then

!        WRITE(6,'(a,200i4)') 'ITPARTN-inf:',itpartn(1:nst(1)+nst(2))
        do ist=1,nstm
          ine    = iokk(iq)+ist-1
          ind    = iorder(ine)+iokk(iq)-1
          equasi =SQRT((spenrg(ind)-fer(icharg(ind)))**2+deltaf(ind)**2)
          if(isTSym.ge.1) then
          write(7,'(i5,i4,f6.1,3i5,2f7.4,2f5.2,f9.4,i8.1, 4f6.3)')
     &       ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),
     &       nparz(1,ind),wtpair(ind),
     &       vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind),itpartn(ind),
     &       rtpartn(ind), angular(ind), equasi,parity(ind)
          else
          write(7,'(i5,i4,f6.1,3i5,2f7.4,2f5.2,f9.5, 4f6.3)')
     &       ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),
     &       nparz(1,ind),wtpair(ind),
     &       vsqf(ind),cmpn(1,ind),cmpn(2,ind),spenrg(ind),angular(ind),
     &       equasi,parity(ind)
          endif
        enddo
        write(7,'(a)') ' '
      else

!        accpartn=0D0
        do ist=1,nstm
          ine    = iokk(iq)+ist-1
          ind    = iorder(ine)+iokk(iq)-1
          equasi =SQRT((spenrg(ind)-fer(icharg(ind)))**2+deltaf(ind)**2)
!          accpartn=vsqf(ind)*wtpair(ind)+accpartn
          CALL calc_angularmomentum(ind,amom,avar)
        if(isTSym.ge.1)then
          write(7,'(i5,i4,f6.1,3i5,3f7.4,f7.3,f15.10,i4,2f6.3,
     &              3(1pg12.4),0pf6.3,f7.3)')
     &       ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),
     &       nparz(1,ind),wtpair(ind),
     &       vsqf(ind),sqrt(max(abs(vsqf(ind)-vsqf(ind)**2),smal)),
     &       deltaf(ind)*wtpair(ind),spenrg(ind),itpartn(ind),
     &       rtpartn(ind),angular(ind),amom,avar,equasi,parity(ind),
     &       calcspquad(ind)
        else
          write(7,'(i5,i4,f6.1,3i5,3f7.4,f7.3,f15.10,f6.3,3(1pg12.4),
     &              0pf6.3,f7.3)')
     &       ind,icharg(ind),mpsi(ind)+0.5,npar(ind),nrpar(ind),
     &       nparz(1,ind),wtpair(ind),
     &       vsqf(ind),sqrt(max(abs(vsqf(ind)-vsqf(ind)**2),smal)),
     &       deltaf(ind)*wtpair(ind),spenrg(ind),angular(ind),amom,avar,
     &       equasi,parity(ind),calcspquad(ind)
          endif
          angular(ind)=amom
        enddo
!        WRITE(7,'(a,1pg13.5)') '# test particle number=',accpartn
        write(7,'(a)') ' '
      endif
c
50000 continue

!     temporarily print deviation from gap equation
      WRITE(7,'(a)') 'Deviation from gap equation:',
     &  ' state  isospin  s.p. energy  gap   dev-gap'
      devgap=0D0
      DO jst=1,nst(1)+nst(2)
        iq = icharg(jst)
        vv2     = vsqf(jst)
        IF(tpairspace) THEN
          deltact=deltaf(jst)
        ELSE
          deltact=deltaf(jst)*wtpair(jst)
        END IF
        equasi = sqrt((spenrg(jst)-fer(iq))**2+deltact**2)
        vv2new=0.5D0-0.5D0*(spenrg(jst)-fer(iq))/equasi
        devgap_state=((spenrg(jst)-fer(iq))*2.0*sqrt(vv2-vv2*vv2)
     &               +deltact*(one-vv2-vv2) )
        WRITE(7,'(i4,i2,8(1pg13.5))') jst,iq,spenrg(jst)-fer(iq),
     &            deltact,devgap_state,
     &            vv2,sqrt(vv2-vv2*vv2),one-vv2-vv2,vv2new,
     &            vv2new-vv2
        devgap=devgap+wate(jst)*devgap_state
      END DO
      WRITE(7,'(a,1pg13.5)') 'total deviation:',devgap
!

c
c     compute variance of total angular momentum
c
!      averj2 = 0D0
!      varyj2 = 0D0
!      DO ist=1,nst(1)+nst(2)
!        averj2 = averj2 + j2squar(ist)*vsqf(ist)
!        averj2 = averj2 + angular(ist)*vsqf(ist)
!        DO jst=1,nst(1)+nst(2)
!        END DO
!      END DO
c
c     optionally save occupation numbers
c
      if(ifoccs.eq.1 .and. ifocsv.eq.0) then
        do iq=1,2
          nstm=nst(iq)
          do ist=1,nstm
            ine     = iokk(iq)+ist-1
            oc(ine) = vsqf(ine)
          enddo
        enddo
        ifocsv = 1
        write(6,'(a)') ' set flag IFOCSV=1'
      endif
c
      return
      end
c
c-----energ -----------------------------------------------------------
c
      subroutine energ()
*
***********************************************************************
*                                                                     *
*        compute the the total energy directly from densities         *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'levels.inc'
c
c
      if(iparam.eq.3) return                   ! exclude case of external field
c
c   prepare derivatives, 'd2rho' is overriden because it should have
c   been done already in 'mfhmlt' and stored.
c
c     prepare 'd2rho' and first derivatives on 'psip' and 'psin'
c
c      call mulr(derr2a(0,0,0),rho(0,0,1),d2rho(0,0,1),one,0)
c      call mulr(derr2a(0,0,0),rho(0,0,2),d2rho(0,0,2),one,0)
c      call mulz(derz2(0,0,1),rho(0,0,1),d2rho(0,0,1),one,1)
c      call mulz(derz2(0,0,1),rho(0,0,2),d2rho(0,0,2),one,1)
c
c     compose combination of skyrme parameters
c
      b0     = t0*(one+oh*x0)
      b0p    = t0*(oh+x0)
      b3     = t3*(one+oh*x3)/4.0
      b3p    = t3*(oh+x3)/4.0
      b1     = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      b1p    = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      b2     = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
      b2p    = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
      b4     = oh*t4


codd
      b6     = oh*oh*t0*x0
      b6p    = oh*oh*t0
      b7     = oh*oh*oh*t3*(x3/3)
      b7p    = oh*oh*oh*(t3/3)
      b8     = oh*oh*oh*oh*(3*t1*x1-t2*x2)
      b8p    = oh*oh*oh*oh*(3*t1+t2)

      if(isTSym.le.2) then
        b6  = zero
        b6p = zero
        b7  = zero
        b7p = zero
        b8  = zero
        b8p = zero
      elseif(isTSym.eq.3) then
        b8  = zero
        b8p = zero
      endif

      excfac = -1.418033*chgx               ! coulomb exchange in Slater appr.
      vdr    = 4.0*edr
c                                             for t_1,t_2 contr. to spin-orbit
      if(ifls.eq.1) then
        t1m20  = (t1-t2)/8.0
        t1p2x  = (x1*t1+x2*t2)/8.0
      else
        t1m20  = 0.0
        t1p2x  = 0.0
      endif
c      write(6,'(a,6(1pg12.4))')
c     &  ' test ls: t1xx,int(rhols)=',
c     &  t1m20,t1p2x,f1int(rholsr),f1int(rholsz),
c     &  f1int(rholsr(0,0,2)),f1int(rholsz(0,0,2))
c
c     compose integrand for energy on 'wrk1'
c
coddtest
*        call nablar_scalar(curphi(0,0,1),rotjz(0,0,1),0,1,one)
*        call nablar_scalar(curphi(0,0,2),rotjz(0,0,2),0,1,one)
*        call nablaz_scalar(curphi(0,0,1),rotjr(0,0,1),1,1,one)
*        call nablaz_scalar(curphi(0,0,2),rotjr(0,0,2),1,1,one)

coddtest
      do ir=0,nr
      do iz=0,nz
        rp      = rho(iz,ir,1)
        rn      = rho(iz,ir,2)
        rs      = max(rn+rp,smal)
        rspow   = rs**power
        tp      = tau(iz,ir,1)
        tn      = tau(iz,ir,2)
        ts      = tp+tn
        d2p     = d2rho(iz,ir,1)
        d2n     = d2rho(iz,ir,2)
        d2s     = d2p+d2n
        ssp     = drhols(iz,ir,1)
        ssn     = drhols(iz,ir,2)
        sss     = ssp+ssn
        rlsrp   = rholsr(iz,ir,1)
        rlszp   = rholsz(iz,ir,1)
        rlsrn   = rholsr(iz,ir,2)
        rlszn   = rholsz(iz,ir,2)
codd

        curp     = curphi(iz,ir,1)
        curn     = curphi(iz,ir,2)
        tc       = curp + curn
        rsigp    = spdenr(iz,ir,1)
        rsign    = spdenr(iz,ir,2)
        rsigs    = spdenr(iz,ir,1) + spdenr(iz,ir,2)
        zsigp    = spdenz(iz,ir,1)
        zsign    = spdenz(iz,ir,2) 
        zsigs    = spdenz(iz,ir,1) + spdenz(iz,ir,2) 
        d2rsigp    = d2spdenr(iz,ir,1)
        d2rsign    = d2spdenr(iz,ir,2)
        d2rsigs    = d2rsigp + larsign
        d2zsigp    = d2spdenz(iz,ir,1)
        d2zsign    = d2spdenz(iz,ir,2) 
        d2zsigs    = d2zsigp + d2zsign
        rotjrp   = rotjr(iz,ir,1)
        rotjrn   = rotjr(iz,ir,2)
        rotjrs   = rotjrp + rotjrn
        rotjzp   = rotjz(iz,ir,1)
        rotjzn   = rotjz(iz,ir,2)
        rotjzs   = rotjzp + rotjzn

        sdrp     = spdenkr(iz,ir,1)
        sdrn     = spdenkr(iz,ir,2)
        sdrs     = sdrp + sdrn
        sdzp     = spdenkz(iz,ir,1)
        sdzn     = spdenkz(iz,ir,2)
        sdzs     = sdzn+sdzp

codd
c
cold        wrk2(iz,ir) = h2meff*ts
        wrk2(iz,ir) = h2meff(1)*tp+h2meff(2)*tn
        wrk1(iz,ir) =
     &                + 0.5*(b0*rs*rs-b0p*(rp*rp+rn*rn))
     &                + rspow/3.0*(b3*rs*rs-b3p*(rp*rp+rn*rn))
     &                + b1*(rs*ts-tc*tc)-b1p*(rp*tp+rn*tn
     &                 - curp*curp - curn*curn)
     &                - 0.5*(b2*rs*d2s-b2p*(rp*d2p+rn*d2n))
     &                - b4*(rs*sss+rsigs*rotjrs + zsigs*rotjzs)
     &              -b4p*(rp*ssp+rn*ssn+ rsigp*rotjrp + rsign*rotjrn
     &                 + zsigp*rotjzp + zsign*rotjzn)
     &                - 0.5*t1p2x*((rlsrp+rlsrn)*(rlsrp+rlsrn)
     &                            +(rlszp+rlszn)*(rlszp+rlszn)
     &                  -rsigs*sdrs-zsigs*sdzs)
     &                + 0.5*t1m20*(rlsrp*rlsrp+rlsrn*rlsrn
     &                            +rlszp*rlszp+rlszn*rlszn
     &            - rsigp*sdrp-zsigp*sdzp - rsign*sdrn-zsign*sdzn)
     &       +(b6+b7*rspow+dglan-dglanp)*(rsigs*rsigs+zsigs*zsigs)
     &       - (b6p+b7p*rspow-2D0*dglanp)*(rsigp*rsigp+zsigp*zsigp+
     &       rsign*rsign+zsign*zsign)
     &      - oh*b8*(rsigs*d2rsigs+zsigs*d2zsigs)
     &      +oh*b8p*(rsigp*d2rsigp+zsigp*d2zsigp+rsign*d2rsign
     &       +zsign*d2zsign)
        if (chg.ne.0.0)
     &     wrk1(iz,ir) =  0.5*rp*coulpt(iz,ir) + wrk1(iz,ir)
        if(chgx.ne.zero)
     &     wrk1(iz,ir) =  0.75*excfac*max(rp,smal)**vdr + wrk1(iz,ir)
           wrk3(iz,ir) =  0.75*excfac*max(rp,smal)**vdr
!        wrk3(iz,ir) = - b4*(rsigs*rotjrs + zsigs*rotjzs)
!     &              -b4p*(rsigp*rotjrp + rsign*rotjrn
!     &                 + zsigp*rotjzp + zsign*rotjzn)
c
      enddo
      enddo
      ekin   = f1int(wrk2)
      epot   = f1int(wrk1)
!      WRITE(7,'(a,1pg15.7)') ' Coulomb exchange=',f1int(wrk3)
      ecolex=f1int(wrk3)
c  details
!      do ir=0,nr
!      do iz=0,nz
!        rp      = rho(iz,ir,1)
!        rn      = rho(iz,ir,2)
!        rs      = max(rn+rp,smal)
!        rspow   = rs**power
!        tp      = tau(iz,ir,1)
!        tn      = tau(iz,ir,2)
!        ts      = tp+tn
!        d2p     = d2rho(iz,ir,1)
!        d2n     = d2rho(iz,ir,2)
!        d2s     = d2p+d2n
!        ssp     = drhols(iz,ir,1)
!        ssn     = drhols(iz,ir,2)
!        sss     = ssp+ssn
!        rlsrp   = rholsr(iz,ir,1)
!        rlszp   = rholsz(iz,ir,1)
!        rlsrn   = rholsr(iz,ir,2)
!        rlszn   = rholsz(iz,ir,2)
!c
!        wrk1(iz,ir) =  0.5*(b0*rs*rs-b0p*(rp*rp+rn*rn))
!     &                + rspow/3.0*(b3*rs*rs-b3p*(rp*rp+rn*rn))
!        wrk2(iz,ir) =  b1*(rs*ts)-b1p*(rp*tp+rn*tn)
!        wrk3(iz,ir) = - 0.5*(b2*rs*d2s-b2p*(rp*d2p+rn*d2n))
!        wrk4(iz,ir) = - b4*rs*sss - b4p*rp*ssp+rn*ssn
!      enddo
!      enddo
!      WRITE(7,'(a,4(1pg15.7))') 'energs:',
!     &    f1int(wrk1),f1int(wrk2),f1int(wrk3),f1int(wrk4)
!      espin  = f1int(wrk3)
!      WRITE(7,'(a,1pg13.5)') 'espin=',espin
c
c     2. radial derivative for compensation of the error in
c     "partial integration on the grid".
c
      call laplacer_scalar(rho(0,0,1),wrk1, 0, 1, one)
      call laplacer_scalar(rho(0,0,2),wrk2, 0, 1, one)
      ekincr = -0.5*(h2meff(1)*f1int(wrk1)+h2meff(2)*f1int(wrk2))
cPK      call mulr(derr2a(0,0,0),rho(0,0,1),wrk3(0,0),-0.5*h2meff,0)
cPK      call mulr(derr2a(0,0,0),rho(0,0,2),wrk3(0,0),-0.5*h2meff,1)
cPK      ekincr = f1int(wrk3)
      do ir=0,nr
      do iz=0,nz
        rp = rho(iz,ir,1)
        rn = rho(iz,ir,2)
        wrk1(iz,ir) = (rp+rn)*(rp+rn)
        wrk2(iz,ir) = rp*rp+rn*rn
      enddo
      enddo
      call laplacer_scalar(wrk1,wrk3, 0, 1, one)
      call laplacer_scalar(wrk2,wrk4, 0, 1, one)
cPK      call mulr(derr2a(0,0,0),wrk1(0,0),wrk3(0,0),one,0)
cPK      call mulr(derr2a(0,0,0),wrk2(0,0),wrk4(0,0),one,0)
      epotcr = -0.125*(b1*f1int(wrk3)-b1p*f1int(wrk4))       ! sign ???

c
c     final composition
c
      epot    = epot+epotcr
      ekin    = ekin+ekincr
      energy = epot+ekin

c
      if(itest.ge.1) then
       write(7,'(a,2(1pg12.5))')
     &  '  ekincr,epotcr=',ekincr,epotcr
      endif
c
      return
      end
c
c-----energ2-----------------------------------------------------------
c
      subroutine energ2()
*
***********************************************************************
*                                                                     *
*        compute the the total energy from sumj of s.p. energies      *
*        and rearrangement energy.                                    *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'levels.inc'
c
c                                                sum of s.p. energies
c
      sumspe = 0.0
      do ist=1,nst(1)+nst(2)
        if (vsqf(ist).gt.1.0d-9) 
     &    sumspe = spenrg(ist)*wate(ist)*vsqf(ist) + sumspe
      enddo
c
c
c                                                two-body potential energy
      do ir=0,nr
      do iz=0,nz
        wrk1(iz,ir) = rho(iz,ir,1)*pot(iz,ir,1)
     &               +rho(iz,ir,2)*pot(iz,ir,2)
     &               +tau(iz,ir,1)*(bmass(iz,ir,1)-h2meff(1))
     &               +tau(iz,ir,2)*(bmass(iz,ir,2)-h2meff(2))
c        wrk2(iz,ir) = (tau(iz,ir,1)+tau(iz,ir,2))
      enddo
      enddo
      epotr  = oh*f1int(wrk1)         ! ?? note: misses half spin-orbit part
c      ekin   = h2meff*f1int(wrk2)
      ekin   = h2meff(1)*f1int(tau(0,0,1))+h2meff(2)*f1int(tau(0,0,2))
c
c     2. radial derivative for compensation of the error in
c     "partial integration on the grid".
c
      call laplacer_scalar(rho(0,0,1),wrk1, 0, 1, one)
      call laplacer_scalar(rho(0,0,2),wrk2, 0, 1, one)
      ekincr = -0.5*(h2meff(1)*f1int(wrk1)+h2meff(2)*f1int(wrk2))
cPK      call mulr(derr2a(0,0,0),rho(0,0,1),wrk3(0,0),one,0)
cPK      call mulr(derr2a(0,0,0),rho(0,0,2),wrk3(0,0),one,1)
cPK      ekincr = -0.5*h2meff*f1int(wrk3)
      ekin   = ekincr+ekin
      if(itest.ge.1) write(7,'(a,2(1pg12.5))')
     &   ' ekin,ekincor=',ekin,ekincr
c                                                compose
      eshell = oh*(sumspe+ekin+equad)
c
c     the rearrangement energy
c
      if(iparam.eq.3) then
       erearr = zero
      else
c                                                   prepare some parameters
       t3a    = t3*(one+oh*x3)/6.0
       t3b    = t3*(oh+x3)/6.0
       t3as   = oh*power*t3a
       t3bs   = oh*power*t3b
       excfac = -1.418033*chgx
       exrfac = -excfac*oh
       vdr    = 4.0*edr

c                                                   prepare integrand
       do ir=0,nr
       do iz=0,nz
         rp      = rho(iz,ir,1)
         rm      = rho(iz,ir,2)
         rs      = rm+rp
         rs      = max(rm+rp,smal)
         rspow   = rs**power
         wrk1(iz,ir) = rspow*(t3as*rs*rs-t3bs*(rp*rp+rm*rm))
         if(chgx.ne.zero)
     &     wrk1(iz,ir) = exrfac*max(rp,smal)**vdr + wrk1(iz,ir)
       enddo
       enddo
c                                                  finally integrate
       erearr = -f1int(wrk1)*oh
       epotr  = epotr+erearr
      endif
c
c     final composition
c
      energr = (eshell+erearr)
c
      return
      end

c
c-----obs--------------------------------------------------------------
c
      subroutine obs(rhoini,atnum,rms,d,q,o,h,qsq,zcm)
c
c     computes the observables
c       atnum  = mass number,
c       rms    = r.m.s. radius,
cbm     d      = spherical dipole moment
c       q      = cartesian quadrupole moment
cbm     o      = spherical octupole moment
c       h      = spherical hexadecupole moment and
c       qsq    = squared quadrupole moment, as estimate of variance
c       zcm    = input: center-of-mass relative to which the above
c                       moments are computed
c     for the density given on 'rhoini'
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
!      common /r4m/r4mom    ! quick and dirty transfer
      dimension rhoini(0:nzd,0:nrd)
c
      data three,seven,eight,twfour /3.0d0,7.0,8.0,24.0/
c
c----------------------------------------------------------------------
c
c     summing up the wavefunction
c
      drezw  = one/12.0D0
      rhot   = 0.0D0
c
cbm   at1    = oh*rhoini(0,0)
      at1    = rhoini(0,0)
      if (isym.eq.1) at1 = at1 * oh
!      r4mom1 = zero
      rms1   = zero
      d1     = zero
      q1     = zero
      qsq1   = zero
      o1     = zero
      h1     = zero
c
      do iz=1,nz
       z     = ((iz)+sz-icm)*dz-zcm
       z2    = z*z
       r2t   = z2
       qt    = z2 + z2
       q2t   = qt*qt
       ot    = two * z * z2
       ht    = eight*z2*z2
       rhot  = rhoini(iz,0)
       at1   = rhot     + at1
       rms1  = rhot*r2t + rms1
!       r4mom1  = rhot*r2t*r2t + r4mom1
       d1    = rhot*z   + d1
       q1    = rhot*qt  + q1
       qsq1  = rhot*q2t  + qsq1
       o1    = rhot*ot  + o1
       h1    = rhot*ht  + h1
      enddo
      d1     = d1 - rhot*z
      o1    =  o1 - rhot*ot
      at1    = drezw * at1
      rms1   = drezw * rms1
      d1     = drezw * d1
      q1     = drezw * q1
      qsq1   = drezw * qsq1
      o1     = drezw * o1
      h1     = drezw * h1
c
      at2    = zero
!      r4mom2  = zero
      rms2   = zero
      d2     = zero
      q2     = zero
      qsq2   = zero
      o2     = zero
      h2     = zero
      z      = (sz) * dz - zcm
      do ir=1,nr
       dir   = (ir)+sr
       z2    = z*z
       r     = dir*dr
       r2    = r*r
       r2t   = r2 + z2
       qt    = z2 + z2 - r2
       q2t   = qt*qt
       ot    = two*z*z2 - three*z*r2
       ht    = eight*z2*z2-twfour*z2*r2+three*r2*r2
       rhot  = rhoini(0,ir)*dir
       at2   = rhot     + at2
       rms2  = rhot*r2t + rms2
!       r4mom2 = rhot*r2t*r2t + r4mom2
       d2    = rhot*z   + d2
       q2    = rhot*qt  + q2
       qsq2  = rhot*q2t  + qsq2
       o2    = rhot*ot  + o2
       h2    = rhot*ht  + h2
      enddo
      if (isym.eq.1) then
        at2    = oh * at2
        rms2   = oh * rms2
        d2     = oh * d2
        q2     = oh * q2
        qsq2   = oh * qsq2
        o2     = oh * o2
        h2     = oh * h2
      endif
c
      at3    = zero
!      r4mom3 = zero
      rms3   = zero
      d3     = zero
      q3     = zero
      qsq3   = zero
      o3     = zero
      h3     = zero
      do ir=1,nr
       dir   = (ir)+sr
       r     = dir*dr
       r2    = r*r
       do iz=1,nz
         z     = ((iz)+sz-icm)*dz-zcm
         z2    = z*z
         r2t   = z2+r2
         qt    = z2+z2-r2
         q2t   = qt*qt
         ot    = two*z*z2 - three*z*r2
         ht    = eight*z2*z2-twfour*z2*r2+three*r2*r2
         rhot  = rhoini(iz,ir)*dir
         at3   = rhot     + at3
         rms3  = rhot*r2t + rms3
!         r4mom3 = rhot*r2t*r2t + r4mom3
         d3    = rhot*z   + d3
         q3    = rhot*qt  + q3
         qsq3  = rhot*q2t  + qsq3
         o3    = rhot*ot  + o3
         h3    = rhot*ht  + h3
       enddo
       d3     = d3 - rhot*z
       o3    =  o3 - rhot*ot
      enddo
c
c     the different parts are summed up:
c
      atnum   =  at1+at2+at3
      rms     =  rms1 +rms2 +rms3
!      r4mom   =  r4mom1 + r4mom2 + r4mom3
      d       =  d1 +d2 +d3
      q       =  q1 +q2 +q3
      qsq     =  qsq1 +qsq2 +qsq3
      o       =  o1 +o2 +o3
      h       =  h1 +h2 +h3
c
c     the final result
c
      rms   = sqrt(abs(rms/max(atnum,smal)))
!      r4mom = abs(r4mom/max(atnum,smal))**0.25D0
c
      vol   = (sym+one)*dvol
c
      atnum = atnum*vol
c
      q     = q  * vol     !!!!!!!!!!!!!
      qsq   = qsq* vol
      h     = h  * vol * three/(two*eight*sqrt(pi))
c
      if (isym.eq.1) then
        d = zero
        o = zero
      else
        d = d * vol * oh * sqrt(three/pi)
        o = o * vol * oh * oh * sqrt(seven/pi)
      endif
c
      return
      end
c-----obsti------------------------------------------------------------
c
      subroutine obsti(rhoin,qa,qc,qa2,qc2)
c
c     computes
c       qa   = constraining momentum from 'qcfac'
c       qc   = momentum of 'qfac'
c       qa2  = squared momentum from 'qcfac^2'
c       qc2  = momentum of 'qfac^2'
c     for density 'rhoin'.
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
c
      dimension rhoin(0:nzd,0:nrd,kqm)
c
c
c     summing up the wavefunction
c
      drezw  = one/12.0
c
      sum1   = zero
      sumd1  = zero
      sums1  = zero
      sumsd1 = zero
      do iz=1,nz
c       z     = ((iz)+sz-icm)*dz
c       z2    = z*z
c       qt    = two*z2
       qt    = qcfac(iz,0)
       qtd   = qfac(iz,0)
       rhotot= rhoin(iz,0,1)+rhoin(iz,0,2)
       sum1  = qt *rhotot+sum1
       sumd1 = qtd*rhotot+sumd1
       sums1 = qt*qt  *rhotot+sums1
       sumsd1= qtd*qtd*rhotot+sumsd1
      enddo
      sum1   =  drezw*sum1
      sumd1  =  drezw*sumd1
      sums1  =  drezw*sums1
      sumsd1 =  drezw*sumsd1
c
      sum2   = zero
      sumd2  = zero
      sums2  = zero
      sumsd2 = zero
      do ir=1,nr
       dir   = (ir)+sr
c       r     = dir*dr
c       r2    = r*r
c       qt    =-r2
       qt    = qcfac(0,ir)
       qtd   = qfac(0,ir)
       rhotot= (rhoin(0,ir,1)+rhoin(0,ir,2))*dir
       sum2  = qt *rhotot+sum2
       sumd2 = qtd*rhotot+sumd2
       sums2 = qt*qt  *rhotot+sums2
       sumsd2= qtd*qtd*rhotot+sumsd2
      enddo
      sum2   = oh*sum2
      sumd2  = oh*sumd2
      sums2  = oh*sums2
      sumsd2 = oh*sumsd2
c
      sum3   = zero
      sumd3  = zero
      sums3   = zero
      sumsd3  = zero
      do ir=1,nr
      dir   = (ir)+sr
c      r     = dir*dr
c      r2    = r*r
      do iz=1,nz
c       z     = ((iz)+sz-icm)*dz
c       z2    = z*z
c       qt    = two*z2-r2
       qt    = qcfac(iz,ir)
       qtd   = qfac(iz,ir)
       rhotot= (rhoin(iz,ir,1)+rhoin(iz,ir,2))*dir
       sum3  = qt *rhotot + sum3
       sumd3 = qtd*rhotot + sumd3
       sums3  = qt*qt  *rhotot + sums3
       sumsd3 = qtd*qtd*rhotot + sumsd3
      enddo
      enddo
c
c
c     the final result
c
      qa  = (sum1+sum2+sum3)   *(sym+one)*dvol
      qc  = (sumd1+sumd2+sumd3)*(sym+one)*dvol
      qa2 = (sums1+sums2+sums3)   *(sym+one)*dvol
      qc2 = (sumsd1+sumsd2+sumsd3)*(sym+one)*dvol
c
      return
      end
c



c-----crinfo -----------------------------------------------------------
c
      subroutine crinfo(outnam)
*
***********************************************************************
*                                                                     *
*     print of quantities from cranking part                          *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
!      common /r4m/r4mom    ! quick and dirty transfer
*
      character*3  outnam
c
c     formfactor and related quantities
c
      rchrg2 = rmsp*rmsp + 0.64  - 2.25/(widdz+widdxy)
      call radius(.TRUE.)
!      IF(ttrack) write(17,'(a,4(1pg13.5))')
!     &  ' testprint r_est,r,R,sigma=',rchrg2,rmsc,rdmsc,surfc
      write(7,'(a,4(1pg13.5))')
     &  ' testprint r_est,r,R,sigma=',rchrg2,rmsc,rdmsc,surfc
!      open(90,status='unknown',form='formatted',file='ensurf.'//outnam,
!     &        position='append')
!      write(90,'(a,4(1pg13.5))')
!     &  ' testprint r,R,sigma=',rmsc,rdmsc,surfc
!      close(90)
c
      anumb  = npart(1)+npart(2)
      if(ifzpe.ne.0) then
       ezp    = 17.3/anumb**0.2
      else
       ezp    = 0.0D0
      endif
c
c     print cranking results
c
*
c      if(swicon.gt.oh) then
      if(swicon.gt.-oh) then
        if(lqcon.eq.3) then
          betawr = beta3
          betawp = beta3p
        elseif(lqcon.eq.2) then
          betawr = beta2
          betawp = beta2p
        else
          betawr = log(rms)
          betawp = log(rmsp)
c          write(6,'(a,2g12.4)') 
c     &      ' INFO: betawr,log(rms)=',betawr,sqrt(log(rms))
        endif     
       write(7,'(a,2(1pg13.5))') 
     &   ' rigid body momentum of inertia and inverse   :',
     &   0.5*xinert,one/xinert
       write(7,'(a,2(1pg13.5))') 
     &   ' inv. scaling mass for quadrupole (cart., b_2):',
     &   bmassc,bmassr
       if(crmasp.ne.zero) write(7,'(a,2(1pg13.5))') 
     &   ' Inglis mass for quadrupole (direct,inverse)  :',
     &   crmasp,one/crmasp
       if(crmaso.ne.zero) write(7,'(a,2(1pg13.5))') 
     &   ' off.d. mass for quadrupole (direct,inverse)  :',
     &   crmaso,one/crmaso
       if(crmasq.ne.zero) write(7,'(a,2(1pg13.5))') 
     &   ' 3.Inglis mass for quadrupole (direct,inverse):',
     &   crmasq,one/crmasq
       if(pwid.ne.zero) write(7,'(a,2(1pg13.5))') 
     &   ' GCM width,  ZPE correction                   :',
     &   pwid,pwid/(2.0*crmasp)
       if(pwido.ne.zero) write(7,'(a,2(1pg13.5))') 
     &   ' off diag. GCM width,  ZPE correction         :',
     &   pwido,pwido/(2.0*crmaso)
       if(ammass.ne.zero) write(7,'(a,3(1pg13.5))') 
     &   ' a.m. mass   a.m. width  a.m. ZPE             :',
     &   ammass,amwid,amwid/(2.0*ammass)
       if(widcmz.ne.zero) write(7,'(a,4(1pg13.5))') 
     &   ' widcm_z&xy, zpe_cm_z&xy                      :',
     &   widcmz,widcmx,zpecmz,zpecmx
!       WRITE(7,*) 'before CRANKN:',widdz,tmass
       if(widdz.ne.zero .and. tmass) then
         open(75,status='unknown',form='formatted',
     &        file='crankn.'//outnam,position='APPEND')
         if(ifzpe.eq.1) then
           zpemod = 
     &        (widdz+widdxy)*0.5D0*(h2m(1)+h2m(2))/(npart(1)+npart(2))
     &        -ezp 
         else
           zpemod = 0.0D0
         endif
!         WRITE(7,*) 'write CRANKN:'
         write(75,'(1x,f7.4,6(1pg12.4),0pf13.5,0pf7.3,4(1pg12.4))')
     &    betawr,scmass,crmasp,crmaso,ammass,amwid,pwid,
     &    energy,zpemod,pnmass,pnmas0,zpetri,ampnms
         close(75)
       endif
c       if(crmasp.ne.zero .and. tmass) then
       if(tmass) then
         open(76,status='unknown',form='formatted',
     &        file='zpecor.'//outnam,position='APPEND')
         rchrg2 = rmsp*rmsp + 0.64  - 2.25/(widdz+widdxy)
         write(7,'(a,3f10.3)') 
     &       ' rmsp2,rchrg2,c.m. r**2=',
     &         rmsp*rmsp,rchrg2,2.25/(widdz+widdxy)
         write(76,'(1x,f6.3,5(1pg12.4),1pg14.6,2(1pg12.4),0pf10.4,
     &              4f9.5)')
     &    betawr,pwid/(2.0*scmass),
     &    pwid/(2.0*crmasp),pwido/(2.0*crmaso),amwid/(2.0*ammass),
     &    amwid,energy-zpemod,zpetri,omgtri,rchrg2,betawp,
     &    rdmsc,surfc,zpemod
         close(76)
         if(kintpl.eq.1) then
           write(77) 
     &       betawr,scmass,crmasp,crmaso,ammass,amwid,pwido,
     &       energy,zpemod,fer,
     &    (spenrg(ist)-fer(1),vsqf(ist),deltaf(ist)*wtpair(ist),
     &     ist=1,nst(1)),
     &    (spenrg(ist)-fer(2),vsqf(ist),deltaf(ist)*wtpair(ist),
     &     ist=nst(1)+1,nst(1)+nst(2))
         endif
       endif
      endif
*
      return
      end
c
c-----printenerg -----------------------------------------------------------
c
      subroutine printenerg()
*
***********************************************************************
*                                                                     *
*        compute the the total energy directly from densities         *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'levels.inc'
c
c
      dimension pot1(0:nzd,0:nrd),pot2(0:nzd,0:nrd),pot3(0:nzd,0:nrd)
      dimension pot4(0:nzd,0:nrd),pot5(0:nzd,0:nrd),pot6(0:nzd,0:nrd)
      dimension pot7(0:nzd,0:nrd),pot8(0:nzd,0:nrd),pot9(0:nzd,0:nrd)
      dimension pot10(0:nzd,0:nrd),pot11(0:nzd,0:nrd),pot12(0:nzd,0:nrd)
      dimension pot13(0:nzd,0:nrd),pot14(0:nzd,0:nrd),pot15(0:nzd,0:nrd)
      dimension pot16(0:nzd,0:nrd),pot17(0:nzd,0:nrd),pot18(0:nzd,0:nrd)
      dimension pot19(0:nzd,0:nrd),pot20(0:nzd,0:nrd),pot21(0:nzd,0:nrd)
      dimension pot22(0:nzd,0:nrd),pot23(0:nzd,0:nrd),pot24(0:nzd,0:nrd)
      dimension pot25(0:nzd,0:nrd)

      if(iparam.eq.3) return                   ! exclude case of external field
c
c   prepare derivatives, 'd2rho' is overriden because it should have
c   been done already in 'mfhmlt' and stored.
c
c     prepare 'd2rho' and first derivatives on 'psip' and 'psin'
c
c      call mulr(derr2a(0,0,0),rho(0,0,1),d2rho(0,0,1),one,0)
c      call mulr(derr2a(0,0,0),rho(0,0,2),d2rho(0,0,2),one,0)
c      call mulz(derz2(0,0,1),rho(0,0,1),d2rho(0,0,1),one,1)
c      call mulz(derz2(0,0,1),rho(0,0,2),d2rho(0,0,2),one,1)
c
c     compose combination of skyrme parameters
c
    
      b0     = t0*(one+oh*x0)
      b0p    = t0*(oh+x0)
      b3     = t3*(one+oh*x3)/4.0
      b3p    = t3*(oh+x3)/4.0
      b1     = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      b1p    = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      b2     = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
      b2p    = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
      b4     = oh*t4



codd      

       b6     = oh*oh*t0*x0
       b6p    = oh*oh*t0
       b7     = oh*oh*oh*t3*x3/3
       b7p    = oh*oh*oh*t3/3
       b8     = oh*oh*oh*oh*oh*(3*t1*x1-t2*x2)
       b8p    = oh*oh*oh*oh*oh*(3*t1+t2)
      if(isTSym.le.2) then
        b6  = zero
        b6p = zero
        b7  = zero
        b7p = zero
        b8  = zero
        b8p = zero
      elseif(isTSym.eq.3) then
        b8  = zero
        b8p = zero
      endif
codd

      excfac = -1.418033*chgx               ! coulomb exchange in Slater appr.
      vdr    = 4.0*edr
c                                             for t_1,t_2 contr. to spin-orbit
      if(ifls.eq.1) then
        t1m20  = (t1-t2)/8.0
        t1p2x  = (x1*t1+x2*t2)/8.0
      else
        t1m20  = 0.0
        t1p2x  = 0.0
      endif
      do iz=0,nz
       do ir=0,nr
        rp      = rho(iz,ir,1)
        rn      = rho(iz,ir,2)
        rs      = max(rn+rp,smal)
        rspow   = rs**power
        tp      = tau(iz,ir,1)
        tn      = tau(iz,ir,2)
        ts      = tp+tn
        d2p     = d2rho(iz,ir,1)
        d2n     = d2rho(iz,ir,2)
        d2s     = d2p+d2n
        ssp     = drhols(iz,ir,1)
        ssn     = drhols(iz,ir,2)
        sss     = ssp+ssn
        rlsrp   = rholsr(iz,ir,1)
        rlszp   = rholsz(iz,ir,1)
        rlsrn   = rholsr(iz,ir,2)
        rlszn   = rholsz(iz,ir,2)
codd

        curp     = curphi(iz,ir,1)
        curn     = curphi(iz,ir,2)
        curs       = curp + curn
        rsigp    = spdenr(iz,ir,1)
        rsign    = spdenr(iz,ir,2)
        rsigs    = spdenr(iz,ir,1) + spdenr(iz,ir,2)
        zsigp    = spdenz(iz,ir,1)
        zsign    = spdenz(iz,ir,2) 
        zsigs    = spdenz(iz,ir,1) + spdenz(iz,ir,2) 
        d2rsigp    = d2spdenr(iz,ir,1)
        d2rsign    = d2spdenr(iz,ir,2)
        d2rsigs    = d2rsigp + larsign
        d2zsigp    = d2spdenz(iz,ir,1)
        d2zsign    = d2spdenz(iz,ir,2) 
        d2zsigs    = d2zsigp + d2zsign
        rotjrp   = rotjr(iz,ir,1)
        rotjrn   = rotjr(iz,ir,2)
        rotjrs   = rotjrp + rotjrn
        rotjzp   = rotjz(iz,ir,1)
        rotjzn   = rotjz(iz,ir,2)
        rotjzs   = rotjzp + rotjzn


        sdrp     = spdenkr(iz,ir,1)
        sdrn     = spdenkr(iz,ir,2)
        sdrs     = sdrp + sdrn
        sdzp     = spdenkz(iz,ir,1)
        sdzn     = spdenkz(iz,ir,2)
        sdzs     = sdzn+sdzp

*    calculate all the densities
        pot1(iz,ir) = 0.5*(b0*rs*rs)
        pot2(iz,ir) = -0.5*b0p*(rp*rp+rn*rn)
        pot3(iz,ir) = b1*rs*ts
        pot4(iz,ir) = -b1*curs*curs
        pot5(iz,ir) = -b1p*(rp*tp+rn*tn)
        pot6(iz,ir) = b1p*(curn*curn+curp*curp)
        pot7(iz,ir) = -0.5*b2*rs*d2s
        pot8(iz,ir) = 0.5*b2p*(rp*d2p +rn*d2n)
        pot9(iz,ir) = (b3/3)*rspow*rs*rs
        pot10(iz,ir) = -(b3p/3)*rspow*(rp*rp+rn*rn)
        pot11(iz,ir) = - b4*rs*sss
        pot12(iz,ir) = - b4*(rsigs*rotjrs + zsigs*rotjzs)
        pot13(iz,ir) = -b4p*(rp*ssp+rn*ssn)
        pot14(iz,ir) = -b4p*(rsigp*rotjrp + rsign*rotjrn
     &                 + zsigp*rotjzp + zsign*rotjzn)

        pot15(iz,ir) = - 0.5*t1p2x*((rlsrp+rlsrn)*(rlsrp+rlsrn)
     &                            +(rlszp+rlszn)*(rlszp+rlszn))
        pot16(iz,ir) = 0.5*t1p2x*(rsigs*sdrs+zsigs*sdzs)
        pot17(iz,ir) =  + 0.5*t1m20*(rlsrp*rlsrp+rlsrn*rlsrn
     &                            +rlszp*rlszp+rlszn*rlszn)
        pot18(iz,ir) = -0.5*t1m20*( rsigp*sdrp+zsigp*sdzp 
     &            +rsign*sdrn+zsign*sdzn)
        pot19(iz,ir) = b6*(rsigs*rsigs + zsigs*zsigs)
        pot20(iz,ir) = b7*rspow*(rsigs*rsigs+zsigs*zsigs)
        pot21(iz,ir) = - b6p*(rsigp*rsigp+zsigp*zsigp+
     &       rsign*rsign+zsign*zsign)
        pot22(iz,ir) = - b7p*rspow*(rsigp*rsigp+zsigp*zsigp+
     &       rsign*rsign+zsign*zsign)
        pot23(iz,ir) = - b8*(rsigs*d2rsigs+zsigs*d2zsigs)
        pot24(iz,ir) =  +b8p*(rsigp*d2rsigp+zsigp*d2zsigp+rsign*d2rsign
     &       +zsign*d2zsign)
         wrk1(iz,ir)=zero
        if (chg.ne.0.0)
     &     wrk1(iz,ir) =  0.5*rp*coulpt(iz,ir) 
        if(chgx.ne.zero)
     &     wrk1(iz,ir) =  0.75*excfac*max(rp,smal)**vdr +wrk1(iz,ir)
c
       enddo
      enddo
* calculate integral
      epot1 = f1int(pot1)
      epot2 = f1int(pot2)
      epot3 = f1int(pot3)
      epot4 = f1int(pot4)
      epot5 = f1int(pot5)
      epot6 = f1int(pot6)
      epot7 = f1int(pot7)
      epot8 = f1int(pot8)
      epot9 = f1int(pot9)
      epot10 = f1int(pot10)
      epot11 = f1int(pot11)
      epot12 = f1int(pot12)
      epot13 = f1int(pot13)
      epot14 = f1int(pot14)
      epot15 = f1int(pot15)
      epot16 = f1int(pot16)
      epot17 = f1int(pot17)
      epot18 = f1int(pot18)
      epot19 = f1int(pot19)
      epot20 = f1int(pot20)
      epot21 = f1int(pot21)
      epot22 = f1int(pot22)
      epot23 = f1int(pot23)
      epot24 = f1int(pot24)
      epotcr1 = f1int(wrk1)
      do ir=0,nr
      do iz=0,nz
        rp = rho(iz,ir,1)
        rn = rho(iz,ir,2)
        wrk1(iz,ir) = (rp+rn)*(rp+rn)
        wrk2(iz,ir) = rp*rp+rn*rn
      enddo
      enddo
      call laplacer_scalar(wrk1,wrk3, 0, 1, one)
      call laplacer_scalar(wrk2,wrk4, 0, 1, one)
      epotcr2 = -0.125*(b1*f1int(wrk3)-b1p*f1int(wrk4))       ! sign ???

      epotsum = epot1+epot2+epot3+epot5+
     & epot4+epot6+epot7+epot8+epot9+epot10+
     & epot11+epot13+epot12+epot14+epot15+epot17+
     & epot16+epot18+epot19+epot21+epot20+epot22+
     & epot23+epot24
      epotcorr = epotcr1 + epotcr2
      anumb  = npart(1)+npart(2)
      if(ifzpe.ne.0) then
       ezp    = 17.3/anumb**0.2
      else
       ezp    = 0.0D0
      endif
      esum = epotsum + epotcorr+ekin+etpair-ezp


*  write the energies in the file


      write(7,25011) epot1,epot2,epot3,epot5,
     & epot4,epot6,epot7,epot8,epot9,epot10,
     & epot11,epot13,epot12,epot14,epot15,epot17,
     & epot16,epot18,epot19,epot21,epot20,epot22,
     & epot23,epot24,epotsum,epotcorr,ekin,esum

25011 format(/' energies (mev)'/,'   rho^2           = ',f11.4,
     x           ' rho_q^2           = ',
     x           f11.4,/,'   rho*tau         = ',
     x  f11.4, ' rho_q*tau_q       = ', f11.4,/,
     x '   j^2             = ', f11.4,' j_q^2             = ',f11.4,/,
     x'   rho*D(rho)      = ',f11.4,' rho_q*D(rho_q)    = ',f11.4
     x ,/,'   rho^(alpha+2)   = ',f11.4,
     x ' rho^alpha*rho_q   = ',f11.4,/,
     x '   rho(div.J)      = ', f11.4, ' rho_q(div.J_q)    = ', f11.4,/,
     x '   sig.(rot j)     = ', f11.4, ' sig_q.(rot j_q)   = ',f11.4,/,
     x '   J^2             = ', f11.4, ' J_q^2             = ', f11.4,/,
     x '   sig.tau         = ',f11.4,  ' sig_q.tau_q       = ', f11.4,/,
     x '   sig^2           = ', f11.4, ' sig_q^2           = ', f11.4,/,
     x '   rho^alpha*sig^2 = ',f11.4,  ' rho^alpha*sig_q^2 = ',f11.4
     x ,/,'   sig*D(sig)      = ',f11.4, ' sig_q*D(sig_q)    = ',f11.4,/,
     x '   E_pot           = ',f11.4,  ' E_pot corr        = ',f11.4, /,
     x '   E_kin           = ', f11.4, ' E_binding         = ', f11.4,/
     x  )
      return
      end
c
c-----printenerg -----------------------------------------------------------
c
      subroutine findconjpart(iq)
*
***********************************************************************
*                                                                     *
*        compute the the total energy directly from densities         *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'levels.inc'
      LOGICAL,PARAMETER :: toverlap=.TRUE.  ! energy or overlap decision


      if(isTSym.ge.1) then
c     search the time-conjugate partners
      if(iq.eq.1)then
       istart= 1
      else
       istart=nst(1)+1
      endif
      do i1=istart,nst(iq)+istart-1
       IF(toverlap) THEN 
         supnorm = zero
       ELSE
         supnorm = 1D30
       END IF
       do i2=istart,nst(iq)+istart-1
       if((mpsi(i1)+mpsi(i2)).eq.-1)then
        IF(toverlap) THEN
         do iz=0,nz
          do ir=0,nr
            wrk1(iz,ir)=psi(iz,ir,1,i1)*psi(iz,ir,2,i2)-
     &      psi(iz,ir,2,i1)*psi(iz,ir,1,i2)
          enddo
         enddo
         a = abs(f1int(wrk1))
!         IF(i1==istblo(iq)) WRITE(7,'(a,2i4,3(1pg13.5))')
!     &    'FCONJG:',i1,i2,spenrg(i1),spenrg(i2),a
         if(supnorm.lt.a) then
           supnorm = a
           itpartn(i1) = i2
           rtpartn(i1) = supnorm
         endif
        ELSE
         a = abs(spenrg(i1)-spenrg(i2))
         if(supnorm.gt.a) then
           supnorm = a
           itpartn(i1) = i2
           rtpartn(i1) = supnorm
         endif
        END IF
       endif
       enddo 
      enddo
      do i1=istart,nst(iq)+istart-1
        IF(i1.NE.itpartn(itpartn(i1)).AND.
     &      (i1==istblo(1).OR.i1==istblo(2))) 
     &    WRITE(7,'(a,3i5,2f5.1,2f9.3)')
     &   'FINDCONJ error: state=',i1,itpartn(i1),itpartn(itpartn(i1)),
     &    mpsi(i1)+0.5,mpsi(itpartn(i1))+0.5,
     &    spenrg(i1),spenrg(itpartn(i1))
      enddo
      ELSE
        itpartn = 0
      endif
      return
      end
c-----calc_parity -----------------------------------------------------------
c
      subroutine calc_parity(parity)
*
***********************************************************************
*                                                                     *
*        compute the the total energy directly from densities         *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'levels.inc'

      real(8) :: parity(1:nst(1)+nst(2))

      IF(isym==1) THEN
        parity(1:nst(1)+nst(2)) = 0D0
        RETURN
      END IF
      do i1=1,nst(1)+nst(2)
        psip = 0D0
        do is=1,2
          isig = 3-2*is
          do iz=1,nz-1
            izc = 2*icm-iz
!            WRITE(*,*) 'iz,izc=',iz,izc
            do ir=0,nr
              psip(iz,ir,is)=isig*psi(izc,ir,is,i1)
            enddo
          enddo
        enddo
        sig = 1D0-2D0*MOD(ABS(mpsi(i1)),2)
        parity(i1) = sig*wfovlp(psip,psi(:,:,:,i1))
      enddo
      return
      end

c------------calcspquad-------------
c
c  caclulates dimensionless quadrupole deformation of an s.p. state
c
      REAL(8) FUNCTION calcspquad(ist)

      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      REAL(8) :: totalbeta = 0D0, totalpartnum=0D0

      wrk1=0D0
      call sdenti(psi(0,0,1,ist),wrk1,1D0,ist)
      call obs(wrk1,spnum,srms,sd,sq,so,sh,sqsq,szcm)
      fac    = SQRT(0.2D0*pi)
      rmsb   = 0.93*atnum**0.3333333                ! 0.93 = sqrt(0.6)*1.2
      spbeta  = fac*sq /max(rmsb*rmsb*spnum,smal)
!      totalbeta = totalbeta + spbeta*vsqf(ist)*wate(ist)
!      totalpartnum = totalpartnum+vsqf(ist)*wate(ist)
!      WRITE(*,*) ' testprint ATNUM,totalbeta=',atnum,
!     &    totalbeta/totalpartnum,totalpartnum
      calcspquad  = spbeta

      RETURN 
      END

c
c------------calangularmomentum-------------
c
c  caclulates total angular momentum of wavefunction
c
      subroutine calangularmomentum(ist,angularmom)

      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'

c calculate derivatives
      call nablaz_spin(psi(0,0,1,ist),psin,nparz(1,ist),1,one)
      call nablar_spin(psi(0,0,1,ist),psip,  mpsi(ist),1,one)

      do iz=0,nz
       do ir=0,nr
         r = (ir+sr)*dr
	 z = (iz+sz-icm)*dz
         if(ir.eq.0)then
          rinv = zero
         else
          rinv = 1.0/r
         endif
c psi1 = L_r|psi>

         if(ir.gt.0)then
	  psi1(iz,ir,1)= rinv*z*mpsi(ist)*psi(iz,ir,1,ist)
	  psi1(iz,ir,2)= rinv*z*(mpsi(ist)+1)*psi(iz,ir,2,ist)
         else
           if(mpsi(ist).eq.1) then
             psi1(iz,ir,1)=z*psip(iz,ir,1) 
             psi1(iz,ir,2)=0
           elseif(mpsi(ist).eq.0) then
             psi1(iz,ir,1)=0
             psi1(iz,ir,2)=z*psip(iz,ir,2) 
           else
             psi1(iz,ir,1)=0
             psi1(iz,ir,2)=0
           endif
         endif
c psi2 = L_phi|psi>
         psi2(iz,ir,1)=r*psin(iz,ir,1) - z*psip(iz,ir,1)
         psi2(iz,ir,2)=r*psin(iz,ir,2) - z*psip(iz,ir,2)

c psi3 = L_z|psi>
         psi3(iz,ir,1)=mpsi(ist)*psi(iz,ir,1,ist)
         psi3(iz,ir,2)=(mpsi(ist)+1)*psi(iz,ir,2,ist)

c psi4 = (L_r*s_r +L_phi*s_phi+L_z*s_z)|psi>
          psi4(iz,ir,1)= 
     & + r*psin(iz,ir,2) - z*psip(iz,ir,2)
     & + mpsi(ist)*psi(iz,ir,1,ist)
         psi4(iz,ir,2)=
     & - r*psin(iz,ir,1) + z*psip(iz,ir,1)
     & - (mpsi(ist)+1)*psi(iz,ir,2,ist)
         if(ir.gt.0) then
         psi4(iz,ir,1)=psi4(iz,ir,1)
     &            -z*rinv*psi(iz,ir,2,ist)*(mpsi(ist)+1) 
         psi4(iz,ir,2)=psi4(iz,ir,2)
     &            -z*rinv*psi(iz,ir,1,ist)*(mpsi(ist)) 
         else
          if(mpsi(ist).eq.0) then
		psi4(iz,ir,1)=psi4(iz,ir,1)-z*psip(iz,ir,2)
          elseif(mpsi(ist).eq.1)then
		psi4(iz,ir,2)=psi4(iz,ir,2)-z*psip(iz,ir,1)
          endif
         endif
       enddo
      enddo
c calculate wavefuntion overlap
!        WRITE(7,*) ' old: ist,ovl=',ist,wfovlp(psi2,psi2)
!        averL2 = wfovlp(psi1,psi1)+wfovlp(psi2,psi2)
!     &         +wfovlp(psi3,psi3)
!        aver4 = wfovlp(psi4,psi(0,0,1,ist))
        eigvalJ = wfovlp(psi1,psi1)+wfovlp(psi2,psi2)
     &         +wfovlp(psi3,psi3)
     &         +wfovlp(psi4,psi(0,0,1,ist))+0.75
        angularmom = oh*(-1.0+ sqrt(1.0+4.0*(eigvalJ)))
!      WRITE(7,'(a,2i4,7(1pg12.4))') ' indirect J: ist,',ist,mpsi(ist),
!     &  averL2,aver4,oh*(-1D0+ sqrt(1D0+4D0*(averL2))),angularmom,
!     &  averJ2b
      
      return
      end


c------------calangularmomentum-------------
c
c  caclulates total angular momentum of wavefunction
c  2. version with L**2 directly
c
      subroutine calc_angularmomentum(ist,angularmom,angularvar)

      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'

c 1. derivatives - 1. round
      call nablaz_spin(psi(0,0,1,ist),psin,nparz(1,ist),1,one)
      call nablar_spin(psi(0,0,1,ist),psip,  mpsi(ist),1,one)

      do iz=0,nz
       do ir=0,nr
         r = (ir+sr)*dr
	 z = (iz+sz)*dz
         psi1(iz,ir,1)= r*psin(iz,ir,1)-z*psip(iz,ir,1)
         psi1(iz,ir,2)= r*psin(iz,ir,2)-z*psip(iz,ir,2)
       enddo
      enddo
      test1=wfovlp(psi1,psi1)
!      call nablar_spin(psip,psin, mpsi(ist),1,one)
!      do iz=0,nz
!       do ir=0,nr
!         r = (ir+sr)*dr
!	 z = (iz+sz)*dz
!         if(ir.eq.0)then
!          rinv = zero
!         else
!          rinv = 1.0/r
!         endif
!         psi2(iz,ir,1)= psin(iz,ir,1)+rinv*psip(iz,ir,1)
!         psi2(iz,ir,2)= psin(iz,ir,2)+rinv*psip(iz,ir,2)
!       enddo
!      enddo
!      WRITE(7,'(a,i4,3(1pg12.4))') ' compare: ist,ovls=',
!     &   ist,wfovlp(psip,psip),-wfovlp(psi(0,0,1,ist),psin),
!     &   -wfovlp(psi(0,0,1,ist),psi2)

c contribution from spin-orbit:
c  psi4 = (L_r*s_r +L_phi*s_phi+L_z*s_z)|psi>
      do iz=0,nz
       do ir=0,nr
         r = (ir+sr)*dr
	 z = (iz+sz)*dz
         if(ir.eq.0)then
          rinv = zero
         else
          rinv = 1.0/r
         endif
         psi4(iz,ir,1)= 
     & + r*psin(iz,ir,2) - z*psip(iz,ir,2)
     & + mpsi(ist)*psi(iz,ir,1,ist)
         psi4(iz,ir,2)=
     & - r*psin(iz,ir,1) + z*psip(iz,ir,1)
     & - (mpsi(ist)+1)*psi(iz,ir,2,ist)
         if(ir.gt.0) then
          psi4(iz,ir,1)=psi4(iz,ir,1)
     &            -z*rinv*psi(iz,ir,2,ist)*(mpsi(ist)+1) 
          psi4(iz,ir,2)=psi4(iz,ir,2)
     &            -z*rinv*psi(iz,ir,1,ist)*(mpsi(ist)) 
         else
          if(mpsi(ist).eq.0) then
		psi4(iz,ir,1)=psi4(iz,ir,1)-z*psip(iz,ir,2)
          elseif(mpsi(ist).eq.1)then
		psi4(iz,ir,2)=psi4(iz,ir,2)-z*psip(iz,ir,1)
          endif
         endif
       enddo
      enddo

c 1. derivatives - 2. round
      nparzco = 3-nparz(1,ist)
      call nablaz_spin(psi1(0,0,1),psin,nparzco,1,one)
      call nablar_spin(psi1(0,0,1),psip,  mpsi(ist),1,one)

      do iz=0,nz
       do ir=0,nr
         r = (ir+sr)*dr
	 z = (iz+sz)*dz
         if(ir.eq.0)then
          rinv = zero
         else
          rinv = 1.0/r
         endif
         psi1(iz,ir,1)= r*psin(iz,ir,1)-z*psip(iz,ir,1)
     &                  -z*rinv*psi1(iz,ir,1)
         psi1(iz,ir,2)= r*psin(iz,ir,2)-z*psip(iz,ir,2)
     &                  -z*rinv*psi1(iz,ir,2)
       enddo
      enddo
!      WRITE(7,*) ' new: ist,ovl=',ist,-wfovlp(psi(0,0,1,ist),psi1),
!     &    test1

c contribution l_0
      mm=mpsi(ist)
      mml=mpsi(ist)+1
      do iz=0,nz
       do ir=0,nr
         r = (ir+sr)*dr
	 z = (iz+sz)*dz
         if(ir.eq.0)then
          rinv = zero
         else
          rinv = 1.0/r
         endif
         psi1(iz,ir,1)=-psi1(iz,ir,1)+mm**2*psi(iz,ir,1,ist)
         psi1(iz,ir,2)=-psi1(iz,ir,2)+mml**2*psi(iz,ir,2,ist)
!         IF(mm>0) THEN
           rad2=z*z*rinv**2
           psi1(iz,ir,1)=psi1(iz,ir,1)+mm**2*rad2*psi(iz,ir,1,ist)
           psi1(iz,ir,2)=psi1(iz,ir,2)+mml**2*rad2*psi(iz,ir,2,ist)
!         END IF
       enddo
      enddo



c compose to full action J**2
!      averL2 = wfovlp(psi(0,0,1,ist),psi1)
      psi1=psi1+psi4+0.75D0*psi(:,:,:,ist)

c calculate wavefuntion overlap
        averJ2 = wfovlp(psi(0,0,1,ist),psi1)
!        averJ2 = averL2 + wfovlp(psi(0,0,1,ist),psi4) +0.75D0
!        aver4 = wfovlp(psi4,psi(0,0,1,ist))
        varyJ2 = wfovlp(psi1,psi1)-averJ2**2
        varyJ2 = MAX(varyJ2,1D-100)
        angularvar = varyJ2**0.25D0
        angularmom = oh*(-1.0+ sqrt(1.0+4.0*(averJ2)))

!      WRITE(7,'(a,2i4,7(1pg12.4))') '   direct J: ist,',ist,mpsi(ist),
!     &  averL2,aver4,oh*(-1D0+ sqrt(1D0+4D0*(averL2))),angularmom

!     scan J**2 matrix
!      j2squar(ist) = varyJ2
!      iq1=icharg(ist)
!      mm=mpsi(ist)
!      nzp=nparz(1,ist)
!      DO jst=1,nst(1)+nst(2)
!        IF(icharg(jst)==iq1 .AND. mpsi(jst)==mm 
!     &     .AND. nzp==nparz(1,jst)) THEN
!           j2matrix(jst,ist) = wfovlp(psi(0,0,1,jst),psi1)
!        END IF
!      END DO

      return
      end
