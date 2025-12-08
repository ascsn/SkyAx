c-----hcommQ-----------------------------------------------------------
c
      subroutine hcommQ(qopin,popin,betain,popout,icase)
c
c     Commutator [H,Q] of Hamiltonian with Q-response matrix,
c     matrix elements of Q are entered on 'qopin',
c     result is returned on 'popout'.
c     The collective mass <[Q,[H,Q]]> is also accumulated and
c     returned on 'scmass'.
c     The particle-number mass is accumulated on
c     'pnmass' = with linear response.
c     The routine requires that the derivatives of the wave functions
c     are already computed and stored in 'dpsir' and 'dpsiz'.
c
c     icase: for test-operatores, icase=0 normal use, for details on
c            icase see source of Oabtst() and hcOT()
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
      dimension qopin(ktriac),popin(ktriac)
      dimension popout(ktriac),betain(ktriac)
c
      data ifhcts/0/                 ! switch for testing
      data ifb4ad/1/                 ! switch for b4&b4p terms
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
      if(KTESTM) stop ' old HCOMMQ needs KTESTM=.FALSE.' 
      if(ifls.ne.0) write(6,*) 
     &   ' Warning: IFLS=0 needed for linear response'
c      if(ifls.ne.0) stop ' IFLS=0 needed for linear response'
c
c     prepare derivatives for later use
c
      do ist=1,nst(1)+nst(2)
        call nablaz_spin(psi(0,0,1,ist),dpsiz(0,0,1,ist),
     &                   nparz(1,ist),1,one)
        call nablar_spin(psi(0,0,1,ist),dpsir(0,0,1,ist),
     &                      mpsi(ist),1,one)
cPK        do is=1,2
cPK          isi    = 3-is
cPK          mm     = mpsi(ist)+is-1
cPK          nparr  = 1+mod(mm,2)
cPK          call mulr(derr1(0,0,nparr),psi(0,0,is,ist),
cPK     &              dpsir(0,0,is,ist),one,0)
cPK          call mulz(derz1(0,0,nparz(is,ist)),psi(0,0,is,ist),
cPK     &              dpsiz(0,0,is,ist),one,0)
cPK        enddo
      enddo
c      write(7,'(a)') ' HCOMMQ: derivatives are prepared'
c
c     1. step: response currents are collected
c
      do ir=0,nr
         do iz=0,nz
            curqr(iz,ir,1)  = 0.0
            curqr(iz,ir,2)  = 0.0
            curqz(iz,ir,1)  = 0.0
            curqz(iz,ir,2)  = 0.0
            spinph(iz,ir,1) = 0.0
            spinph(iz,ir,2) = 0.0
         enddo
      enddo
c
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         vi     = sqrtb(vsqf(ist))
         ui     = sqrtb(one-vi*vi)
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            eferm  = fer(1)
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
            eferm  = fer(2)
            jstmax = nst(1)+nst(2)
         endif
         do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               vj     = sqrtb(vsqf(jst))
               uj     = sqrtb(one-vj*vj)
CR030399 check
cold            factor = (ui*ui*vj*vj+vi*vi*uj*uj)*qopin(ipostr)
               if(ist.ne.jst) then ! exlude diagonal element ? yes, should be zero.
                  factor = (ui*ui*vj*vj-vi*vi*uj*uj)*qopin(ipostr) 
                  do ir=0,nr
                     do iz=0,nz
C
C factor 1/2 ?
C
CR030399 check factor ?
                        curqr(iz,ir,iq)  = -2.0        !0.5 ! j_Q,q;r
     &                       *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                       + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                       - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)
     &                       - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist))
     &                       *factor  +  curqr(iz,ir,iq)
                     enddo
                     do iz=0,nz
CR030399 check
                        curqz(iz,ir,iq)  = -2.0        !0.5 ! j_Q,q;z
     &                       *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)
     &                       + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)
     &                       - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)
     &                       - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist))
     &                       *factor  +  curqz(iz,ir,iq)
c                        spinph(iz,ir,iq) = 0.5 ! sigma_Q,q;z
c     &                       *(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
c     &                       - psi(iz,ir,1,ist)*psi(iz,ir,2,jst)
c     &                       - psi(iz,ir,2,jst)*psi(iz,ir,1,ist)
c     &                       + psi(iz,ir,1,jst)*psi(iz,ir,2,ist))
c     &                       *factor  +  spinph(iz,ir,iq)
                     enddo
                     do iz=0,nz
CR030399 check
                        spinph(iz,ir,iq) = 2.0 ! sigma_Q,q;z 
     &                       *(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
     &                       - psi(iz,ir,1,ist)*psi(iz,ir,2,jst))
     &                       *factor  +  spinph(iz,ir,iq)
                     enddo
                  enddo
               endif
            endif
         enddo
      enddo
      iposmx = ipostr
c     
crold      iact = min(kcranc,1)
c############################################################################
CR250599
c test output of currents
c
      if(icase.ne.0) then
         call ptstcur(icase)
      endif
CR250599 end of currents test output
c############################################################################
C
c
c     2. step: response potentials are composed
c
      b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      b1c = b1-b1p
      if(ifb4ad.ne.1) then
        b4  = 0.0   ! t4*0.5
        b4c = 0.0   ! b4+b4p
      else
        b4  = t4*0.5
        b4c = b4+b4p
      endif
c     
      if(isym.eq.1)then
         nparzs = 2
      else
         nparzs = 1
      endif
C     ! rot(sigma_Q,q)
cnonsense      call mulr(derr1(0,0,2),spinph(iz,ir,1),wrk3,one,0)
cnonsense      call mulr(derr1(0,0,2),spinph(iz,ir,2),wrk4,one,0)
cnonsense      call mulz(derz1(0,0,nparzs),spinph(iz,ir,1),wrk1,one,0)
cnonsense      call mulz(derz1(0,0,nparzs),spinph(iz,ir,2),wrk2,one,0)
      call nablaz_scalar(spinph(0,0,1),wrk1,nparzs,1,one)
      call nablaz_scalar(spinph(0,0,2),wrk2,nparzs,1,one)
      call nablar_scalar(spinph(0,0,1),wrk3,     1,1,one)
      call nablar_scalar(spinph(0,0,2),wrk4,     1,1,one)
cPK      call mulr(derr1(0,0,2),spinph(0,0,1),wrk3,one,0)
cPK      call mulr(derr1(0,0,2),spinph(0,0,2),wrk4,one,0)
cPK      call mulz(derz1(0,0,nparzs),spinph(0,0,1),wrk1,one,0)
cPK      call mulz(derz1(0,0,nparzs),spinph(0,0,2),wrk2,one,0)
c     ! A_Q,q
      do ir=0,nr
         do iz=0,nz
            dynpor(iz,ir,1)  = 
     &           -2.0*(b1c*curqr(iz,ir,1) + b1*curqr(iz,ir,2))
     &           +     b4c*wrk1(iz,ir)    + b4*wrk2(iz,ir)
            dynpor(iz,ir,2)  = 
     &           -2.0*(b1c*curqr(iz,ir,2) + b1*curqr(iz,ir,1))
     &           +     b4c*wrk2(iz,ir)    + b4*wrk1(iz,ir)
            dynpoz(iz,ir,1)  = 
     &           -2.0*(b1c*curqz(iz,ir,1) + b1*curqz(iz,ir,2))
     &           -     b4c*wrk3(iz,ir)    - b4*wrk4(iz,ir)
            dynpoz(iz,ir,2)  = 
     &           -2.0*(b1c*curqz(iz,ir,2) + b1*curqz(iz,ir,1))
     &           -     b4c*wrk4(iz,ir)    - b4*wrk3(iz,ir)
         enddo
      enddo
c     ! rot(j_Q,q)
cnonsense      call mulr(derr1(0,0,1),curqz(iz,ir,1),wrk3,one,0)
cnonsense      call mulr(derr1(0,0,1),curqz(iz,ir,2),wrk4,one,0)
cnonsense      call mulz(derz1(0,0,1),curqr(iz,ir,1),wrk1,one,0)
cnonsense      call mulz(derz1(0,0,1),curqr(iz,ir,2),wrk2,one,0)
      call nablaz_scalar(curqr(0,0,1),wrk1, 1,1,one)
      call nablaz_scalar(curqr(0,0,2),wrk2, 1,1,one)
      call nablar_scalar(curqz(0,0,1),wrk3, 0,1,one)
      call nablar_scalar(curqz(0,0,2),wrk4, 0,1,one)
cPK      call mulr(derr1(0,0,1),curqz(0,0,1),wrk3,one,0)
cPK      call mulr(derr1(0,0,1),curqz(0,0,2),wrk4,one,0)
cPK      call mulz(derz1(0,0,1),curqr(0,0,1),wrk1,one,0)
cPK      call mulz(derz1(0,0,1),curqr(0,0,2),wrk2,one,0)
c     ! S_Q,q
      do ir=0,nr
         do iz=0,nz
            spipop(iz,ir,1)  = -b4c*(wrk1(iz,ir)-wrk3(iz,ir))
     &                         - b4*(wrk2(iz,ir)-wrk4(iz,ir))
            spipop(iz,ir,2)  = -b4c*(wrk2(iz,ir)-wrk4(iz,ir))
     &                         - b4*(wrk1(iz,ir)-wrk3(iz,ir))
         enddo
      enddo
c
c     3. step: accumulate matrix elements by spatial integration with
c              transition currents
c
      pnmass = 0.0
      pqnorm = 0.0
      scmass = 0.0
      tstmas = 0.0              ! Inglis mass as test
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            eferm  = fer(1)
            jstmax = nst(1)      
         else                   ! neutrons
            iq     = 2
            eferm  = fer(2)
            jstmax = nst(1)+nst(2)
         endif
         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-eferm
         equasi = sqrt(edifi*edifi+deltai*deltai)
         vi     = sqrtb(vsqf(ist))
         ui     = sqrtb(one-vi*vi)
c
         do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               vj     = sqrtb(vsqf(jst))
               uj     = sqrtb(one-vj*vj)
               deltaj = deltaf(jst)*wtpair(jst)
               edifj  = spenrg(jst)-eferm
               equasj = sqrt(edifj*edifj+deltaj*deltaj)
c     ! h_Q,alpha,beta
               do ir=0,nr
                  do iz=0,nz
                     wrk1(iz,ir) = 
     &                    -0.5                                     ! sign?
     &                    *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                    + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                    - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)
     &                    - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist))
     &                    * dynpor(iz,ir,iq)
c                         
     &                    -0.5                                    ! sign ????
     &                    *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)
     &                    + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)
     &                    - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)
     &                    - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist))
     &                    * dynpoz(iz,ir,iq)
c     
c     &                    +0.5
c     &                    *(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
c     &                    - psi(iz,ir,1,ist)*psi(iz,ir,2,jst)
c     &                    - psi(iz,ir,2,jst)*psi(iz,ir,1,ist)
c     &                    + psi(iz,ir,1,jst)*psi(iz,ir,2,ist))
c     &                    * spipop(iz,ir,iq)
     &                    +(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
     &                    - psi(iz,ir,1,ist)*psi(iz,ir,2,jst))
     &                    * spipop(iz,ir,iq)
                  enddo
               enddo
               popout(ipostr) =  qopin(ipostr)*(equasi+equasj)
     &              + (ui*vj-vi*uj)*f1int(wrk1)
               if(ist.eq.jst) then
                  scmass = 2.0 * popout(ipostr)*qopin(ipostr)   + scmass
                  pqnorm = 2.0 * popin(ipostr)*betain(ipostr)   + pqnorm
ctest                  tstmas = popin(ipostr)**2/(equasi+equasj)     + tstmas
                  if(ifhcts.eq.1)
     &            tstmas = 2.0 * qopin(ipostr)**2*(equasi+equasj)+tstmas
                  pnmass = 2.0 * qopin(ipostr)**2
     &                         *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
               else
                  scmass = 4.0*popout(ipostr)*qopin(ipostr)     + scmass
                  pqnorm = 4.0*popin(ipostr)*betain(ipostr)     + pqnorm
ctest                  tstmas = 2.0*popin(ipostr)**2/(equasi+equasj) + tstmas
                  if(ifhcts.eq.1)
     &            tstmas = 4.0*qopin(ipostr)**2*(equasi+equasj) + tstmas
                  pnmass = 4.0 * qopin(ipostr)**2
     &                         *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
               endif
            endif
         enddo
      enddo
cPK      if(ifhcts.eq.1) then
      write(7,'(a,4(1pg12.4))')
     &     ' hcommQ: scmass,tstmas,crmasp,pqnorm=',
     &               scmass,tstmas,crmasp,pqnorm
cPK   endif
c      write(6,'(a,i6)') ' HCOMMQ: ipostr=',ipostr
      if(ifb4ad.ne.1) then
        b4  = t4*0.5
        b4c = b4+b4p
      endif
c
c     all done
c
      return
      end 
