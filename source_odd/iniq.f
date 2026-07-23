c-----inie2q-----------------------------------------------------------
c
      subroutine inie2q(mcomp,npcomp)
c
c     Initializes the two-quasi-particle energies on 'e2quas' in
c     'common /crankin/'.
c     The storage sequence depends on angular momentum 'mcomp'
c     and z-parity of the mode 'npcomp'.
c     This routine is to be called first for each mode seperately.
c     It also defines the storage space 'iposmx'.
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'crankin.inc'
c
      logical KTEST
      parameter (KTEST=.FALSE.)
c
c
c-----------------------------------------------------------------------
c
c
      if(mcomp.eq.0 .and. .NOT.KTESTM) then
        ipostr = 0             ! set position for lower triangular storage
        do ist=1,nst(1)+nst(2)
          mist   = mpsi(ist)
          npist  = nparz(1,ist)
          if(ist.le.nst(1)) then
            jstmax = nst(1)
            eferm  = fer(1)
          else
            jstmax = nst(1)+nst(2)
            eferm  = fer(2)
          endif
          deltai = deltaf(ist)*wtpair(ist)
          edifi  = spenrg(ist)-eferm
          equasi = sqrt(edifi*edifi+deltai*deltai)
c
          do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
              deltaj = deltaf(jst)*wtpair(jst)
              edifj  = spenrg(jst)-eferm
              equasj = sqrt(edifj*edifj+deltaj*deltaj)
              ipostr = 1 + ipostr
              if(ipostr.gt.ktriac) 
     &           stop "INIE2Q: tridiagonal storage out of range"
c
              e2quas(ipostr)  = equasi+equasj
c
              if(KTEST)  write(6,'(a,6i5,2f5.1,g12.4)')
     &           'INIE2Q: ist,jst,mpsi(ist),mgoal,npist,npgoal,ki,kj=',
     &               ist,jst,mpsi(ist),mpsi(ist),npist,npist,
     &               (2*mpsi(ist)+1)/2.0,(2*mpsi(ist)+1)/2.0,
     &               e2quas(ipostr)
c
            endif
          enddo
        enddo
        write(6,'(a,i6)')  ' INIE2Q: initialized, IPOSTR=',ipostr
        iposmx = ipostr
c
      elseif(mcomp.gt.0 .or. KTESTM) then
        rmcomp = mcomp*1.0D0
        ipostr = 0
        do ksign=-1,1,2                ! to span all 'k' values 
        iup   = (3-ksign)/2            ! effective upper component
        idw   = (3+ksign)/2            ! effective lower component
        do ist=1,nst(1)+nst(2)
          npist  = nparz(1,ist)*ksign
          kist2  = (2*mpsi(ist)+1)*ksign     ! 2*'k'
          kgoal  = kist2+2*mcomp             ! double 'k' for final state
          if(kgoal.gt.0) then
            mgoal  = (kgoal-1)/2             ! spatial a.m. for final state
            ksig2  = 1                       ! sign of final 'k'
            npgoal = npist*npcomp            ! parity of final state
            jup    = 1                       ! no time reversal
            jdw    = 2
          else
            mgoal  = (-kgoal-1)/2
            ksig2  = -1
            npgoal = -npist*npcomp
            jup    = 2                       ! time reversal
            jdw    = 1
          endif
          if(ist.le.nst(1)) then ! protons
            iq     = 1
            eferm  = fer(1)
            jstmin = 1                       ! lower bound of loop for final
            jstmax = nst(1)                  ! upper bound of loop for final
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
c
            if(ipostr.gt.ktriac) 
     &           stop "INE2Q: tridiagonal storage out of range"
c
            deltaj = deltaf(jst)*wtpair(jst)
            edifj  = spenrg(jst)-eferm
            equasj = sqrt(edifj*edifj+deltaj*deltaj)
            e2quas(ipostr)  = equasi+equasj
c
            if(KTEST)  write(6,'(a,6i5,2f5.1,g12.4)')
     &      'INIE2Q: ist,jst,mpsi(ist),mgoal,npist,npgoal,ki,kj=',
     &               ist,jst,mpsi(ist),mgoal,npist,npgoal,
     &               kist2/2.0,kgoal/2.0,e2quas(ipostr)
c
          endif
          enddo
        enddo
        enddo
        write(6,'(a,i6)')  ' INIE2Q: initialized, IPOSTR=',ipostr
        iposmx = ipostr
c
      else
      endif
c
      return
      end

c-----quadms-----------------------------------------------------------
c
      subroutine quadms(qopout,mcomp,npcomp)
c
c     To compute the simple mass using qadrupole momenta and equasi**-3
c     results is 'crmasq' stored on common .
c     Q operator is returned on 'qopout'.
c     The z-component of angular momentum is 'mcomp' whereby only
c     'mcomp=0' or 'mcomp=2' make sense. 
c     The z-parity of the mode 'npcomp'; it must be +1.
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
      dimension qopout(ktriac)
c
c      data itstpr/0/
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
      if(npcomp.ne.1) stop 'QUADMS: parity NPCOMP must not be negative' 
c
c PK  qufac  = sqrt(0.2*pi)/(rms*rms*atnum)
      qufac  = sqrt(0.2D0*pi)/(0.6D0*1.2D0**2*atnum**(5.0D0/3.0D0))
      write(6,'(a,2f8.3)') ' QUADMS: rms,A=',rms,atnum
c
c     the component 'mcomp=0'
c
      if(mcomp.eq.0 .and. .NOT.KTESTM) then
        sumcrs = zero
        sumnor = zero
c  --- test PK ---
        sumbeta= zero
c  --- test PK ---
        ipostr = 0             ! set position for lower triangular storage
        do ist=1,nst(1)+nst(2)
          mist   = mpsi(ist)
          npist  = nparz(1,ist)
          if(ist.le.nst(1)) then
            jstmax = nst(1)
          else
            jstmax = nst(1)+nst(2)
          endif
          vi    = sqrtb(vsqf(ist))
          ui    = sqrtb(one-vi*vi)
          do is=1,nsd
          do ir=0,nr
            r     = dr*(ir+sr)
            r2    = r*r
            do iz=0,nz
              z     = ((iz)+sz-icm)*dz
              z2    = z*z
              psip(iz,ir,is) = (z2+z2-r2)*psi(iz,ir,is,ist)
            enddo
          enddo
          enddo
C
          do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
              vj    = sqrtb(vsqf(jst))
              uj    = sqrtb(one-vj*vj)
              qij    = wfovlp(psi(0,0,1,jst),psip)*qufac
              ipostr = 1 + ipostr
              if(ipostr.gt.ktriac) 
     &           stop "QUADMS: tridiagonal storage out of range"
c --- test PK ---
              if (ist.eq.jst) then
                 sumbeta = 2.0D0*vsqf(ist)*qij + sumbeta
              endif
c --- test PK ---
                
              qij    = qij*(vj*ui+uj*vi)
              qopout(ipostr)  = qij
              if(ist.eq.jst) then
                 qij    = qij*qij
              else
                 qij    = qij*qij*2.0
              endif
              sumcrs  = qij/e2quas(ipostr)**3 + sumcrs
              sumnor  = qij/e2quas(ipostr) + sumnor
c                sumcrs  = qij/(equasi+equasj)**3 + sumcrs
c                sumnor  = qij/(equasi+equasj) + sumnor
            endif
          enddo
        enddo
c
c     the component 'mcomp=2'
c
      elseif(mcomp.eq.2 .or. KTESTM) then
        rmcomp = mcomp*1.0D0
        ipostr = 0
        do ksign=-1,1,2                ! to span all 'k' values 
        iup   = (3-ksign)/2            ! effective upper component
        idw   = (3+ksign)/2            ! effective lower component
        do ist=1,nst(1)+nst(2)
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
          vi    = sqrtb(vsqf(ist))
          ui    = sqrtb(one-vi*vi)
          deltai = deltaf(ist)*wtpair(ist)
          edifi  = spenrg(ist)-eferm
          equasi = sqrt(edifi*edifi+deltai*deltai)
          do is=1,nsd
          do ir=0,nr
            r     = dr*(ir+sr)
            r2    = r*r
            do iz=0,nz
              z     = ((iz)+sz-icm)*dz
              z2    = z*z
              if(KTESTM .and. mcomp.eq.0) then
                psip(iz,ir,is) = (z2+z2-r2)*psi(iz,ir,is,ist)
              else
                psip(iz,ir,is) = (r2)*psi(iz,ir,is,ist)
              endif
            enddo
          enddo
          enddo
          do jst=jstmin,jstmax
          if(mpsi(jst).eq.mgoal .and. nparz(1,jst).eq.npgoal) then
            ipostr = 1 + ipostr
            vj    = sqrtb(vsqf(jst))
            uj    = sqrtb(one-vj*vj)
            qij    = wfovlp(psi(0,0,1,jst),psip)*qufac
            qij    = qij*(vj*ui+uj*vi)
            qopout(ipostr)  = qij
            qij    = qij*qij
            sumcrs  = qij/e2quas(ipostr)**3 + sumcrs
            sumnor  = qij/e2quas(ipostr) + sumnor
          endif
          enddo
        enddo
        enddo
c
c     everything else
c
      else
        stop 'QUADMS: m-component other than 0 or 2 not possible'
      endif
c
      crmasq = sumcrs/(sumnor*sumnor)
      if(ipostr.ne.iposmx) 
     &   stop 'in QUADMS: inconsistent nr. of matrix elements'
      write(7,'(a,4(1pg12.4),i6)')
     &  'QUADMS: sumcrs,sumnor,crmasq,beta,ipostr=',
     &   sumcrs,sumnor,crmasq,sumbeta,ipostr
c
      return
      end
C-----cringl-----------------------------------------------------------
c
      subroutine cringl(popout)
c
c     Prepares collective momenum operator on 'popout' and
c     computes Inglis cranking mass for collectiev (quadrupole) motion
c     by building derivative from previously saved wavefunction.
c     Actual wavefunctions are saved on scratch for subsequent call.
c
c     This routine contains some test switches e.g. for testing
c     'nabla' operator (use only without z-symmetry).
c     The switch 'KTESTM' invokes rectangular storage instead
c     of triangular one; this is provided for testing other
c     routines as e.g. QCURRE etc.
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
      dimension popout(ktriac)
c
      data ifirst/0/            
      data itstcm/0/            ! switch to test case c.m. motion
      data itstpr/0/            ! switch to testprints of Inglis Cranking
ctest      data isavsp/1/            ! save s.p. properties on file for interpol.
c
cold      data e2quas/ktriac*1.0/
      data delb/0.0/
c      data dummy/0.0/
      save beto2,beto2p,beto2n,beto3,ferolp,feroln
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
c     double loop over states with twofold intention:
c      1. to compute and store the matrix of collective derivative 'popout'
c      2. to accumulate Inglis cranking mass from derivatives
c
      if(itstcm.eq.1) then
         ifirst = 1
         betao2=0.0
         betao3=0.0
         if(sym.gt.zero) stop ' test c.m. requires SYM=0'
      endif
c
      if(ifirst.ne.0) then
         if(lqcon.eq.2) then
           delb   = beta2-beto2
         elseif(lqcon.eq.3) then
           delb   = beta3-beto3
         endif
         if(itstpr.ne.0) 
     &        write(6,'(a,3(1pg12.4))')
     &        ' BETAs=',delb,beta2,beto2
         sumcrp = zero
         sumdig = zero
         sumwdp = zero
         sumwdd = zero
         ipostr = 0             ! set position for lower triangular storage
         isigmin = 1
         if(KTESTM) isigmin=-1
         do isigm=isigmin,1,2
         do ist=1,nst(1)+nst(2)
            mist   = mpsi(ist)
            npist  = nparz(1,ist)
            if(ist.le.nst(1)) then
               if(KTESTM) jstmin = 1 
               jstmax = nst(1)
c               eferm  = fer(1)
CR               efermo = ferolp
            else
               if(KTESTM) jstmin = nst(1)+1 
               jstmax = nst(1)+nst(2)
c               eferm  = fer(2)
CR               efermo = feroln
            endif
            if(.NOT.KTESTM) jstmin = ist
            if(itstcm.ne.1) then
               read(98)
     &              (((psip(iz,ir,is),iz=0,nz),ir=0,nr),is=1,nsd),
     &              speol,vsqol,delol,wtpol
               if(itstpr.gt.0) write(7,'(a,i5,4(1pg12.4))')
     &              ' ist,speol,vsqol,delol,wtpol=',ist,speol,vsqol,
     &              delol,wtpol
            endif
c                                         test overlaps
            if(itstpr.ne.0) then
               do jst=1,jstmax          
                  if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) 
     &                 write(6,'(a,2i4,2(1pg12.4))')
     &                 ' test overlaps: ist,jst,ov=',
     &                 ist,jst,wfovlp(psip,psi(0,0,1,jst)),
     &                 wfovlp(psi(0,0,1,ist),psi(0,0,1,jst))
               enddo
            endif
c                                         derivative
            if(itstcm.ne.1) then
               do is=1,nsd
                  do ir=0,nr
                     do iz=0,nz
                        psip(iz,ir,is) 
     &                       = (psi(iz,ir,is,ist)-psip(iz,ir,is))
     &                       /delb
                     enddo
                  enddo
               enddo
c
c     gradient operator Nabla_z on 'psip' 
c
            else        
               do is=1,nsd
                  mm        = mist+is-1
                  npis     = nparz(is,ist)
                  call mulz(derz1(0,0,npis),psi(0,0,is,ist),wrk4,one,0)
                  do ir=0,nr
                     do iz=0,nz
                        psip(iz,ir,is) = wrk4(iz,ir)
                     enddo
                  enddo
               enddo
            endif
c                                         auxiliary variables
            vi    = sqrtb(vsqf(ist))
            ui    = sqrtb(one-vi*vi)
c            deltai = deltaf(ist)*wtpair(ist)
c            edifi  = spenrg(ist)-eferm
c            equasi = sqrt(edifi*edifi+deltai*deltai)
c                                         diagonal part of cranking mass
            vold   = sqrtb(vsqol)
            uold   = sqrtb(one-vsqol)
            vnew   = vi
            unew   = ui
            uaver  = 0.5*(unew+uold)
            vaver  = 0.5*(vnew+vold)
            if(vaver.gt.0.5*one) then
               pii    = -(unew-uold)/(vaver*delb)
            else
               pii    = (vnew-vold)/(uaver*delb)
            endif
cold            ipostr = 1 + ipostr
            if(ipostr.gt.ktriac) 
     &         stop "CRINGL: tridiagonal storage out of range"
cold            if(kcrank.gt.0) then
cold              popout(ipostr)  = pii     
ccold              e2quas(ipostr) = equasi+equasi
cold            endif
c            sumdig = pii*pii/(equasi+equasi) + sumdig
            sumwdd = pii*pii + sumwdd
ctest            write(7,'(a,5(1pg12.4))')
ctest     &        '     uold,vold,unew,vnew,delb=',uold,vold,unew,vnew,delb
c                                         off-diagonal elements
c
c
            do jst=jstmin,jstmax
               if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
                  ipostr = 1 + ipostr
                  vj    = sqrtb(vsqf(jst))
                  uj    = sqrtb(one-vj*vj)
c                  deltaj = deltaf(jst)*wtpair(jst)
c                  edifj  = spenrg(jst)-eferm
c                  equasj = sqrt(edifj*edifj+deltaj*deltaj)
                  pij    = wfovlp(psi(0,0,1,jst),psip)
                  if(itstpr.ne.0) 
     &               write(6,'(a,2i4,2f8.2,4(1pg12.4))')
     &                 ' i,j,equasi/j,phase,pij=',
     &                 ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
                  pij    = pij*(vj*ui-uj*vi)
                  if(ist.eq.jst) then
                     if(kcrank.gt.0) popout(ipostr)  = pii  
                     sumdig = pii*pii/e2quas(ipostr) + sumdig
                     pij     = pij*pij
                     if(itstpr.eq.1) write(7,'(a,2i5,2(1pg12.4))')
     &               ' ist,ipostr,pii,equasi=',ist,ipostr,pii,equasi
                  else
                     if(kcrank.gt.0) popout(ipostr)  = pij
                     if(KTESTM) then
                       pij     = pij*pij
                     else
                       pij     = 2.0*pij*pij
                     endif
                  endif
                  sumcrp = pij/e2quas(ipostr) + sumcrp
                  sumwdp = pij + sumwdp
               endif
            enddo
         enddo
         rewind 98
         enddo
         if(KTESTM) then
           crmaso = sumcrp
           crmasp = (sumcrp+sumdig)
           pwido  = sumwdp*0.5D0
           pwid   = (sumwdp+sumwdd)*0.5D0
         else
           crmaso = sumcrp*2.0
           crmasp = (sumcrp+sumdig)*2.0
           pwido  = sumwdp
           pwid   = sumwdp+sumwdd
         endif
c
         write(6,'(a,4g12.4,i5)')
     &   ' CRINGL: crmaso,crmasp,pwido,pwid,ipostr=',
     &             crmaso,crmasp,pwido,pwid,ipostr
         if(ipostr.ne.iposmx) 
     &      stop 'CRINGL: inconsistent nr. of matrix elements '
c
      endif
      if(KTESTM)  write(7,'(a,i10)') ' CRINGL: IPOSTR=',ipostr
c
c     save actual wavefunction and associated s.p. properties
c
      if(ifirst.eq.0) then
         open(98,form='unformatted',status='scratch')
         ifirst = 1
      else
        rewind 98
      endif
      beto2  = beta2
      beto3  = beta3
      beto2p = beta2p
      beto2n = beta2n
      ferolp = fer(1)
      feroln = fer(2)
      do ist=1,nst(1)+nst(2)
         write(98)
     &        (((psi(iz,ir,is,ist),iz=0,nz),ir=0,nr),is=1,nsd),
     &        spenrg(ist),vsqf(ist),deltaf(ist),wtpair(ist)
         write(6,'(a,i4)') ' scratch file written. record nr.',ist
      enddo
      write(98)   ifirst        ! dummy record
      rewind 98
c
c     optionally save s.p. properties for later interpolation
c
ctest      if(isavsp.eq.1) then
ctest        if(ifirst.eq.0) then
ctest           open(97,form='unformatted',file='testint')
ctest           write(97) nst
ctest        endif
ctest        write(97) 
ctest     &    beta2,crmaso,pwido,fer,
ctest     &    (spenrg(ist)-fer(1),vsqf(ist),deltaf(ist)*wtpair(ist),
ctest     &     ist=1,nst(1)),
ctest     &    (spenrg(ist)-fer(2),vsqf(ist),deltaf(ist)*wtpair(ist),
ctest     &     ist=nst(1)+1,nst(1)+nst(2))
ctest      endif
c
c
      return
      end
