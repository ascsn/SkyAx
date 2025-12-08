c-----hcommQ-----------------------------------------------------------
c
      subroutine hcommQ(qopin,popin,betaop,popout,icase)
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
      dimension qopin(ktriac),popin(ktriac),popout(ktriac),
     &          betaop(ktriac)
c
      data ifhcts/0/                 ! switch for testing
      data ifb4ad/1/                 ! switch for b4&b4p terms
c
c     buffered SQRT as internal statement function
c
c      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
      if(ifls.ne.0) write(6,*) 
     &   ' Warning: IFLS=0 needed for linear response'
c
c     present settings
c
      mcomp  = 0
      npcomp = 1
c
      call Qcurre(qopin,mcomp,npcomp)
      call Qpoten(mcomp,npcomp)
      call Qhact(qopin,popout,betaop,mcomp,npcomp)
c
      return
      end
c-----Qcurre-----------------------------------------------------------
c
      subroutine Qcurre(qopin,mcomp,npcomp)
c
c     Computes the response currents and spin densities for the Q-like
c     operator given on 'qopin'.
c     The m-component of the operator is given on 'mcomp'.
c     It should satisfy 'mcomp>=0'. The case 'mcomp<0' is trivially
c     obtained by time reversal.
c     Note that the case 'mcomp=0' is much simpler and uses a special 
c     switch within that subroutine.
c     The z-parity of the operator is 'npcomp'.
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
      logical KTEST
      parameter (KTEST=.FALSE.)
c
      dimension qopin(ktriac)
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
c      if(ifls.ne.0) stop ' IFLS=0 needed for linear response'
c
c     prepare derivatives for later use
c
      do ist=1,nst(1)+nst(2)
        do is=1,2
          isi    = 3-is
          mm     = mpsi(ist)+is-1
          nparr  = 1+mod(mm,2)
          call mulr(derr1(0,0,nparr),psi(0,0,is,ist),
     &              dpsir(0,0,is,ist),one,0)
          call mulz(derz1(0,0,nparz(is,ist)),psi(0,0,is,ist),
     &              dpsiz(0,0,is,ist),one,0)
        enddo
      enddo
      write(7,'(a)') ' QCURRE: derivatives are prepared'
c
c     response currents are collected
c
      do iq=1,2
      do ir=0,nr
      do iz=0,nz
            curqr(iz,ir,iq)  = 0.0
            curqz(iz,ir,iq)  = 0.0
            curqph(iz,ir,iq)  = 0.0
            spinr(iz,ir,iq) = 0.0
            spinz(iz,ir,iq) = 0.0
            spinph(iz,ir,iq) = 0.0
      enddo
      enddo
      enddo
c
c     the case 'mcomp=0'
c
      if(mcomp.eq.0 .and. .NOT.KTESTM) then
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
CR0303  99 check
cold              factor = (ui*ui*vj*vj+vi*vi*uj*uj)*qopin(ipostr)
               if(ist.ne.jst) then ! exlude diagonal element ? yes, should be zero.
                 factor = (ui*ui*vj*vj-vi*vi*uj*uj)*qopin(ipostr) 
                 do ir=0,nr
                 do iz=0,nz
C  
                   curqr(iz,ir,iq)  = -2.0             ! j_Q,q;r
     &                       *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                       + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                       - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)
     &                       - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist))
     &                       *factor  +  curqr(iz,ir,iq)
                   curqz(iz,ir,iq)  = -2.0             ! j_Q,q;z
     &                       *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)
     &                       + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)
     &                       - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)
     &                       - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist))
     &                       *factor  +  curqz(iz,ir,iq)
                   spinph(iz,ir,iq) = 2.0              ! sigma_Q,q;z 
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
c     the case 'mcomp>0'
c
      elseif(mcomp.gt.0 .or. KTESTM) then
        rmcomp = mcomp*1.0D0
        ipostr = 0
        do ksign=-1,1,2                ! to span all 'k' values 
        iup   = (3-ksign)/2            ! effective upper component
        idw   = (3+ksign)/2            ! effective lower component
        do ist=1,nst(1)+nst(2)
           vi     = sqrtb(vsqf(ist))
           ui     = sqrtb(one-vi*vi)
           npist  = nparz(1,ist)*ksign
           kist2  = (2*mpsi(ist)+1)*ksign     ! 2*'k'
           kgoal  = kist2+2*mcomp
           if(kgoal.gt.0) then
             mgoal  = (kgoal-1)/2
             ksig2  = 1
             npgoal = npist*npcomp
             jup    = 1
             jdw    = 2
           else
             mgoal  = (-kgoal-1)/2
             ksig2  = -1
             npgoal = -npist*npcomp
             jup    = 2
             jdw    = 1
           endif
           if(ist.le.nst(1)) then ! protons
              iq     = 1
              eferm  = fer(1)
              jstmin = 1
              jstmax = nst(1)
           else                   ! neutrons
              iq     = 2
              eferm  = fer(2)
              jstmin = nst(1)+1
              jstmax = nst(1)+nst(2)
           endif
           do jst=jstmin,jstmax
           if(mpsi(jst).eq.mgoal .and. nparz(1,jst).eq.npgoal) then
              ipostr = 1 + ipostr
              vj     = sqrtb(vsqf(jst))
              uj     = sqrtb(one-vj*vj)
              factor = (ui*ui*vj*vj-vi*vi*uj*uj)*qopin(ipostr) 
c              write(6,'(a,3i5,2(1pg12.4))')
c     &        'ist,jst,ipostr,factor,qopin(ipostr)=',
c     &         ist,jst,ipostr,factor,qopin(ipostr)

              if(KTEST)  write(6,'(a,6i5,2f5.1)')
     &        'QCURRE: ist,jst,mpsi(ist),mgoal,npist,npgoal,ki,kj=',
     &                ist,jst,mpsi(ist),mgoal,npist,npgoal,
     &                kist2/2.0,kgoal/2.0
              do ir=0,nr
              do iz=0,nz
                curqr(iz,ir,iq)  = -1.0                    ! j_Q,q;r
     &               *( psi(iz,ir,iup,ist)*dpsir(iz,ir,jup,jst)
     &                 +psi(iz,ir,idw,ist)*dpsir(iz,ir,jdw,jst)
     &                 -psi(iz,ir,jup,jst)*dpsir(iz,ir,iup,ist)
     &                 -psi(iz,ir,jdw,jst)*dpsir(iz,ir,idw,ist))
     &                        *factor  +  curqr(iz,ir,iq)
                curqz(iz,ir,iq)  = -1.0                    ! j_Q,q;z
     &               *( psi(iz,ir,iup,ist)*dpsiz(iz,ir,jup,jst)
     &                 +psi(iz,ir,idw,ist)*dpsiz(iz,ir,jdw,jst)
     &                 -psi(iz,ir,jup,jst)*dpsiz(iz,ir,iup,ist)
     &                 -psi(iz,ir,jdw,jst)*dpsiz(iz,ir,idw,ist))
     &                       *factor  +  curqz(iz,ir,iq)
                curqph(iz,ir,iq) = rcomp                   ! j_Q,q;phi
     &               *( psi(iz,ir,idw,ist)*psi(iz,ir,jup,jst)
     &                 +psi(iz,ir,iup,ist)*psi(iz,ir,jdw,jst))
     &                       *factor  +  curqph(iz,ir,iq)
                spinr(iz,ir,iq)  = 1.0                     ! sigma_Q,q;r
     &               *( psi(iz,ir,idw,ist)*psi(iz,ir,jup,jst)
     &                 +psi(iz,ir,iup,ist)*psi(iz,ir,jdw,jst))
     &                       *factor  +  spinr(iz,ir,iq)
                spinz(iz,ir,iq)  = 1.0                     ! sigma_Q,q;z 
     &               *( psi(iz,ir,iup,ist)*psi(iz,ir,jup,jst)
     &                 -psi(iz,ir,idw,ist)*psi(iz,ir,jdw,jst))
     &                       *factor  +  spinz(iz,ir,iq)
                spinph(iz,ir,iq) = 1.0                     ! sigma_Q,q;phi
     &               *( psi(iz,ir,idw,ist)*psi(iz,ir,jup,jst)
     &                 -psi(iz,ir,iup,ist)*psi(iz,ir,jdw,jst))
     &                       *factor  +  spinph(iz,ir,iq)
              enddo
              enddo
           endif
           enddo
         enddo
         enddo
c
c     the unwanted 'mcomp<0'
c
      else
        stop "QCURRE: mcomp < 0 not yet implemented "
      endif
c
      write(6,'(a,6(1pg12.4))')
     &  ' norms:',f1int(curqr),f1int(curqz),f1int(curqph),
     &           f1int(spinr),f1int(spinz),f1int(spinph)
c     
      return
      end

c-----Qpoten-----------------------------------------------------------
c
      subroutine Qpoten(mcomp,npcomp)
c
c     Computes the time-odd response potentials 'dynpoz'...
c     for given currents and spin densities 'jcurz'...
c     The m-component of the operator is given on 'mcomp'.
c     It should satisfy 'mcomp>=0'. The case 'mcomp<0' is trivially
c     obtained by time reversal.
c     Note that the case 'mcomp=0' is simpler and uses a special 
c     switch within that subroutine.
c     The z-parity of the operator is 'npcomp'.
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
      data ifb4ad/1/                 ! switch for b4 & b4p terms
c
c-----------------------------------------------------------------------
c
      if(ifls.ne.0) write(6,*) 
     &   ' Warning: IFLS=0 needed for linear response'
      if(ifls.ne.0) write(7,*) 
     &   ' Warning: IFLS=0 needed for linear response'
c
c     prepare parameters
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
      mparop = 1-2*mod(mcomp,2)
      npc    = (3-npcomp)/2
      npco   = (3+npcomp)/2
c
c     case 'mcomp=0'
c
      if(mcomp.eq.0) then     
c
C       rot(sigma_Q,q)
c
        call mulr(derr1(0,0,2),spinph(0,0,1),wrk3,one,0)
        call mulr(derr1(0,0,2),spinph(0,0,2),wrk4,one,0)
        call mulz(derz1(0,0,nparzs),spinph(0,0,1),wrk1,one,0)
        call mulz(derz1(0,0,nparzs),spinph(0,0,2),wrk2,one,0)
c     ! A_Q,q
        do ir=0,nr
           do iz=0,nz
              dynpor(iz,ir,1)  = 
     &             -2.0*(b1c*curqr(iz,ir,1) + b1*curqr(iz,ir,2))
     &             +     b4c*wrk1(iz,ir)    + b4*wrk2(iz,ir)
              dynpor(iz,ir,2)  = 
     &             -2.0*(b1c*curqr(iz,ir,2) + b1*curqr(iz,ir,1))
     &             +     b4c*wrk2(iz,ir)    + b4*wrk1(iz,ir)
              dynpoz(iz,ir,1)  = 
     &             -2.0*(b1c*curqz(iz,ir,1) + b1*curqz(iz,ir,2))
     &             -     b4c*wrk3(iz,ir)    - b4*wrk4(iz,ir)
              dynpoz(iz,ir,2)  = 
     &             -2.0*(b1c*curqz(iz,ir,2) + b1*curqz(iz,ir,1))
     &             -     b4c*wrk4(iz,ir)    - b4*wrk3(iz,ir)
           enddo
        enddo
c
c       rot(j_Q,q)
c
        call mulr(derr1(0,0,1),curqz(0,0,1),wrk3,one,0)
        call mulr(derr1(0,0,1),curqz(0,0,2),wrk4,one,0)
        call mulz(derz1(0,0,1),curqr(0,0,1),wrk1,one,0)
        call mulz(derz1(0,0,1),curqr(0,0,2),wrk2,one,0)
c       ! S_Q,q
        do ir=0,nr
        do iz=0,nz
          spipop(iz,ir,1)  = -b4c*(wrk1(iz,ir)-wrk3(iz,ir))
     &                       - b4*(wrk2(iz,ir)-wrk4(iz,ir))
          spipop(iz,ir,2)  = -b4c*(wrk2(iz,ir)-wrk4(iz,ir))
     &                       - b4*(wrk1(iz,ir)-wrk3(iz,ir))
        enddo
        enddo
c
c     case 'mcomp>0'
c
      elseif(mcomp.gt.0) then     
c
c       (rot j)_phi  and  S_phi
c
        call mulr(derr1(0,0,mparop),curqz(0,0,1),wrk1,one,0)
        call mulz(derz1(0,0,npc),curqr(0,0,1),wrk1,-one,1)
        call mulr(derr1(0,0,mparop),curqz(0,0,2),wrk2,one,0)
        call mulz(derz1(0,0,npc),curqr(0,0,2),wrk2,-one,1)
        do ir=0,nr
        do iz=0,nz
          spipop(iz,ir,1)  = -b4c*wrk1(iz,ir) - b4*wrk2(iz,ir)
          spipop(iz,ir,2)  = -b4c*wrk2(iz,ir) - b4*wrk1(iz,ir)
        enddo
        enddo
c
c       (rot j)_z  and  S_z
c
        do ir=0,nr
        do iz=0,nz
          wrk1(iz,ir)  = mcomp*curqr(iz,ir,1)       ! imaginary ?
          wrk2(iz,ir)  = mcomp*curqr(iz,ir,2)
        enddo
        enddo
        call mulr(derr1(0,0,mparop),curqph(0,0,1),wrk1,-one,1)
        call mulr(derr1(0,0,mparop),curqph(0,0,2),wrk2,-one,1)
        do ir=0,nr
        do iz=0,nz
          spipoz(iz,ir,1)  = -b4c*wrk1(iz,ir) - b4*wrk2(iz,ir)
          spipoz(iz,ir,2)  = -b4c*wrk2(iz,ir) - b4*wrk1(iz,ir)
        enddo
        enddo
c
c       (rot j)_r  and  S_r
c
        do ir=0,nr
        do iz=0,nz
          wrk1(iz,ir)  = -mcomp*curqz(iz,ir,1)       ! imaginary ?
          wrk2(iz,ir)  = -mcomp*curqz(iz,ir,2)
        enddo
        enddo
        call mulz(derz1(0,0,npc),curqph(0,0,1),wrk1,one,1)
        call mulz(derz1(0,0,npc),curqph(0,0,2),wrk2,one,1)
        do ir=0,nr
        do iz=0,nz
          spipor(iz,ir,1)  = -b4c*wrk1(iz,ir) - b4*wrk2(iz,ir)
          spipor(iz,ir,2)  = -b4c*wrk2(iz,ir) - b4*wrk1(iz,ir)
        enddo
        enddo
c
c       (rot sigma)_phi  and  A_phi
c
        call mulr(derr1(0,0,mparop),spinz(0,0,1),wrk1,one,0)
        call mulz(derz1(0,0,npc),spinr(0,0,1),wrk1,-one,1)
        call mulr(derr1(0,0,mparop),spinz(0,0,2),wrk2,one,0)
        call mulz(derz1(0,0,npc),spinr(0,0,2),wrk2,-one,1)
        do ir=0,nr
        do iz=0,nz
          dynpop(iz,ir,1)  = -b4c*wrk1(iz,ir) - b4*wrk2(iz,ir)
     &             -2.0*(b1c*curqph(iz,ir,1) + b1*curqph(iz,ir,2))
          dynpop(iz,ir,2)  = -b4c*wrk2(iz,ir) - b4*wrk1(iz,ir)
     &             -2.0*(b1c*curqph(iz,ir,2) + b1*curqph(iz,ir,1))
        enddo
        enddo
c
c       (rot sigma)_z  and  A_z
c
        do ir=0,nr
        do iz=0,nz
          wrk1(iz,ir)  = -mcomp*spinz(iz,ir,1)       ! imaginary ?
          wrk2(iz,ir)  = -mcomp*spinz(iz,ir,2)
        enddo
        enddo
        call mulr(derr1(0,0,mparop),spinph(0,0,1),wrk1,-one,1)
        call mulr(derr1(0,0,mparop),spinph(0,0,2),wrk2,-one,1)
        do ir=0,nr
        do iz=0,nz
          dynpoz(iz,ir,1)  = -b4c*wrk1(iz,ir) - b4*wrk2(iz,ir)
     &             -2.0*(b1c*curqz(iz,ir,1) + b1*curqz(iz,ir,2))
          dynpoz(iz,ir,2)  = -b4c*wrk2(iz,ir) - b4*wrk1(iz,ir)
     &             -2.0*(b1c*curqz(iz,ir,2) + b1*curqz(iz,ir,1))
        enddo
        enddo
c
c       (rot sigma)_r  and  A_r
c
        do ir=0,nr
        do iz=0,nz
          wrk1(iz,ir)  = mcomp*spinr(iz,ir,1)       ! imaginary ?
          wrk2(iz,ir)  = mcomp*spinr(iz,ir,2)
        enddo
        enddo
        call mulz(derz1(0,0,npco),spinph(0,0,1),wrk1,one,1)
        call mulz(derz1(0,0,npco),spinph(0,0,2),wrk2,one,1)
        do ir=0,nr
        do iz=0,nz
          dynpor(iz,ir,1)  = -b4c*wrk1(iz,ir) - b4*wrk2(iz,ir)
     &             -2.0*(b1c*curqz(iz,ir,1) + b1*curqz(iz,ir,2))
          dynpor(iz,ir,2)  = -b4c*wrk2(iz,ir) - b4*wrk1(iz,ir)
     &             -2.0*(b1c*curqr(iz,ir,2) + b1*curqr(iz,ir,1))
        enddo
        enddo
c
c     unwanted case 'mcomp<0'
c
      else
        stop "QPOTEN: mcomp < 0"
      endif
c
      return
      end


c-----Qhact-----------------------------------------------------------
c
      subroutine Qhact(qopin,popout,betaop,mcomp,npcomp)
c
c     Action of time-odd response Hamiltonian.
c     'qopin' is the matrix of the Q-like operator on which H acts and
c     #popout' is the matrix of resulting P-like operator.
c     The m-component of the operator is given on 'mcomp'.
c     It should satisfy 'mcomp>=0'. The case 'mcomp<0' is trivially
c     obtained by time reversal.
c     Note that the case 'mcomp=0' is simpler and uses a special 
c     switch within that subroutine.
c     The z-parity of the operator is 'npcomp'.
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
      dimension qopin(ktriac),popout(ktriac),betaop(ktriac)
c
      data ifb4ad/1/                 ! switch for b4 & b4p terms
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
c
c     prepare derivatives for later use
c
      do ist=1,nst(1)+nst(2)
        do is=1,2
          isi    = 3-is
          mm     = mpsi(ist)+is-1
          nparr  = 1+mod(mm,2)
          call mulr(derr1(0,0,nparr),psi(0,0,is,ist),
     &              dpsir(0,0,is,ist),one,0)
          call mulz(derz1(0,0,nparz(is,ist)),psi(0,0,is,ist),
     &              dpsiz(0,0,is,ist),one,0)
        enddo
      enddo
      write(7,'(a)') ' HCOMMQ: derivatives are prepared'
c
c     accumulate matrix elements by spatial integration with
c     transition currents
c
      if(mcomp.eq.0 .and. .NOT.KTESTM) then
        pnmass = 0.0
        scmass = 0.0
        tstmas = 0.0              ! Inglis mass as test
        pqnorm = 0.0
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
     &                 -0.5                                     ! sign?
     &                 *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                 + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                 - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)
     &                 - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist))
     &                 * dynpor(iz,ir,iq)
c                         
     &                 -0.5                                    ! sign ????
     &                 *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)
     &                 + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)
     &                 - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)
     &                 - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist))
     &                 * dynpoz(iz,ir,iq)
c     
     &                 +(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
     &                 - psi(iz,ir,1,ist)*psi(iz,ir,2,jst))
     &                 * spipop(iz,ir,iq)
               enddo
             enddo
             popout(ipostr) =  qopin(ipostr)*(equasi+equasj)
     &              + (ui*vj-vi*uj)*f1int(wrk1)
             if(ist.eq.jst) then
               scmass = 2.0 * popout(ipostr)*qopin(ipostr)   + scmass
               pqnorm = 2.0 * popout(ipostr)*betaop(ipostr)  + pqnorm
               if(ifhcts.eq.1)
     &           tstmas = 2.0 * qopin(ipostr)**2*(equasi+equasj)+tstmas
                 pnmass = 2.0 * qopin(ipostr)**2
     &                        *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
               else
                 scmass = 4.0*popout(ipostr)*qopin(ipostr)     + scmass
                 pqnorm = 4.0*popout(ipostr)*betaop(ipostr)    + pqnorm
ctest                  tstmas = 2.0*popin(ipostr)**2/(equasi+equasj) + tstmas
                 if(ifhcts.eq.1)
     &             tstmas = 4.0*qopin(ipostr)**2*(equasi+equasj) + tstmas
                 pnmass = 4.0 * qopin(ipostr)**2
     &                        *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
               endif
             endif
           enddo
         enddo
c         if(ifhcts.eq.1) then
         write(7,'(a,4(1pg12.4))')
     &     ' hcommQ: scmass,tstmas,crmasp,norm=',
     &       scmass,tstmas,crmasp,pqnorm
c         endif
         write(6,'(a,i6)') ' HCOMMQ: ipostr=',ipostr
c
c     mcomp>0
c
      elseif(mcomp.gt.0 .or. KTESTM) then
        pnmass = 0.0
        scmass = 0.0
        tstmas = 0.0              ! Inglis mass as test
        tstres = 0.0
        ipostr = 0
        rmcomp = mcomp*1.0D0
        ipostr = 0
        do ksign=-1,1,2                ! to span all 'k' values 
        iup   = (3-ksign)/2            ! effective upper component
        idw   = (3+ksign)/2            ! effective lower component
        do ist=1,nst(1)+nst(2)
          vi     = sqrtb(vsqf(ist))
          ui     = sqrtb(one-vi*vi)
          npist  = nparz(1,ist)*ksign
          kist2  = (2*mpsi(ist)+1)*ksign     ! 2*'k'
          kgoal  = kist2+2*mcomp
          if(kgoal.gt.0) then
            mgoal  = (kgoal-1)/2
            ksig2  = 1
            npgoal = npist*npcomp
            jup    = 1
            jdw    = 2
          else
            mgoal  = (-kgoal-1)/2
            ksig2  = -1
            npgoal = -npist*npcomp
            jup    = 2
            jdw    = 1
          endif
          if(ist.le.nst(1)) then ! protons
             iq     = 1
             eferm  = fer(1)
             jstmin = 1
             jstmax = nst(1)
          else                   ! neutrons
             iq     = 2
             eferm  = fer(2)
             jstmin = nst(1)+1
             jstmax = nst(1)+nst(2)
          endif
          deltai = deltaf(ist)*wtpair(ist)
          edifi  = spenrg(ist)-eferm
          equasi = sqrt(edifi*edifi+deltai*deltai)
          do jst=jstmin,jstmax
          if(mpsi(jst).eq.mgoal .and. nparz(1,jst).eq.npgoal) then
            ipostr = 1 + ipostr
            vj     = sqrtb(vsqf(jst))
            uj     = sqrtb(one-vj*vj)
            deltaj = deltaf(jst)*wtpair(jst)
            edifj  = spenrg(jst)-eferm
            equasj = sqrt(edifj*edifj+deltaj*deltaj)
            xmup   = ((kist2+kgoal-2)/2)*0.5D0
            xmdw   = ((kist2+kgoal+2)/2)*0.5D0
c     ! h_Q,alpha,beta
            do ir=0,nr
              do iz=0,nz
                wrk1(iz,ir) = 
     &                -0.5                                     ! sign?
     &                *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)
     &                - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist))
     &                * dynpor(iz,ir,iq)
c                        
     &                -0.5                                    ! sign ????
     &                *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)
     &                + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)
     &                - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)
     &                - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist))
     &                * dynpoz(iz,ir,iq)
c
     &                +0.5                                     ! sign?
     &                *(xmup*psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                + xmdw*psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                - xmup*psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)
     &                - xmdw*psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist))
     &                * dynpop(iz,ir,iq)
c     
     &                +(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
     &                + psi(iz,ir,1,ist)*psi(iz,ir,2,jst))
     &                * spipor(iz,ir,iq)
c     
     &                +(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                - psi(iz,ir,2,ist)*psi(iz,ir,2,jst))
     &                * spipoz(iz,ir,iq)
c     
     &                +(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
     &                - psi(iz,ir,1,ist)*psi(iz,ir,2,jst))
     &                * spipop(iz,ir,iq)
              enddo
            enddo
            popout(ipostr) =  qopin(ipostr)*(equasi+equasj)
     &             + (ui*vj-vi*uj)*f1int(wrk1)
            tstres = tstres  + (ui*vj-vi*uj)*f1int(wrk1)
c            write(6,'(a,2i5,1pg12.4)')
c     &        ' ist,jst,resint=',ist,jst,f1int(wrk1)
            if(ist.eq.jst) then
              scmass = 1.0 * popout(ipostr)*qopin(ipostr)   + scmass
              if(ifhcts.eq.1)
     &          tstmas = 1.0 * qopin(ipostr)**2*(equasi+equasj)+tstmas
              pnmass = 1.0 * qopin(ipostr)**2
     &                       *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
            else
              scmass = 1.0*popout(ipostr)*qopin(ipostr)     + scmass
ctest                  tstmas = 2.0*popin(ipostr)**2/(equasi+equasj) + tstmas
              if(ifhcts.eq.1)
     &          tstmas = 1.0*qopin(ipostr)**2*(equasi+equasj) + tstmas
              pnmass = 1.0 * qopin(ipostr)**2
     &                       *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
            endif
            endif
          enddo
        enddo
        enddo
        if(ifhcts.eq.1) write(7,'(a,3(1pg12.4))')
     &    ' hcommQ: scmass,tstmas,crmasp=',scmass,tstmas,crmasp
        write(6,'(a,i6,1pg12.4)') 
     &    ' HCOMMQ: ipostr,tstres=',ipostr,tstres
c
c     the unwanted mcomp<0
c
      else
        stop "QHACT: mcomp < 0"
      endif
c
      if(ifb4ad.ne.1) then
        b4  = t4*0.5
        b4c = b4+b4p
      endif
c
      return
      end
