c-----hcommP-----------------------------------------------------------
c
      subroutine hcommP(popin,qopout)
c
c     Commutator [H,P] of Hamiltonian with P-response matrix,
c     matrix elements of P are entered on 'popin',
c     result is returned on 'qopout'.
c     These matrices carry the matrix elements of the operator as
c     two-quasiparticle excitations. They are stored in upper
c     triangular form.
c
c     The routine uses the array 'dpsir' in 'common /crankin/' as
c     workspace. 
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
      dimension popin(ktriac),qopout(ktriac)
c
      data stepq0/0.000001/          ! determines step size --> 'stepr'
c
c     internal statement function for extrapolation to r=0
      data trhalf,trfift,oneten/1.5,0.6,0.1/
      extrap(f1,f2,f3) = trhalf*f1-trfift*f2+oneten*f3

c
c-----------------------------------------------------------------------
c
      if(KTESTM) stop ' old HCOMMQ needs KTESTM=.FALSE.' 
      if(ifls.ne.0) stop ' IFLS=0 needed for linear response'
      if(ipair.ne.5 .and. ipair.ne.9) 
     &    stop ' IPAIR=5 needed for linear P response'
c
c     1. step: accumulate response currents
c
      do iq=1,2
      do ir=0,nr
      do iz=0,nz
        rhop(iz,ir,iq) = 0.0
        taup(iz,ir,iq) = 0.0
        drhlsp(iz,ir,iq) = 0.0
        paircp(iz,ir,iq) = 0.0
      enddo
      enddo
      enddo
c
      stepr = 0.01/sqrt(opovlp(popin,popin))
c
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         vi     = sqrt(vsqf(ist))
         ui     = sqrt(one-vi*vi)
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
         addfac = 1.0
         do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               xm2    = mist*mist
               xm2l   = (mist+1)*(mist+1)
               ipostr = 1 + ipostr
               vj     = sqrt(vsqf(jst))
               uj     = sqrt(one-vj*vj)
               factor = addfac*(ui*vj+vi*uj)*popin(ipostr)
               factpa = addfac*(ui*uj-vi*vj)*popin(ipostr)
               do ir=1,nr
                 r      = (sr+ir)*dr
                 zentfu = xm2/(r*r)
                 zentfl = xm2l/(r*r)
                 do iz=0,nz
                   rhoacu =  psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
                   rhoacl =  psi(iz,ir,2,ist)*psi(iz,ir,2,jst)
                   rhop(iz,ir,iq) = factor*(rhoacu+rhoacl) 
     &                             + rhop(iz,ir,iq)
                   paircp(iz,ir,iq) = factpa*(rhoacu+rhoacl) 
     &                               + paircp(iz,ir,iq)
                   addtau = factor*
     &                       (dpsir(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                        +dpsir(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                        +dpsiz(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)
     &                        +dpsiz(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)
     &                        +zentfu*rhoacu+zentfl*rhoacl         )
                   wrk1(iz,ir) = addtau          ! save that for spin-orbit
                   taup(iz,ir,iq) = addtau + taup(iz,ir,iq)
                 enddo
               enddo
               if(t4.ne.zero .OR. ioptsk==13) then
                 do ir=1,nr
                   rinv  = one/((sr+ir)*dr)
                   xmup  = mist*rinv
                   xmlow = (mist+1)*rinv
                   do iz=0,nz
                     spsiui = dpsiz(iz,ir,1,ist)+dpsir(iz,ir,2,ist)
     &                     +xmlow*psi(iz,ir,2,ist)
                     spsili = -dpsiz(iz,ir,2,ist)+dpsir(iz,ir,1,ist)
     &                     -xmup*psi(iz,ir,1,ist)
                     spsiuj = dpsiz(iz,ir,1,jst)+dpsir(iz,ir,2,jst)
     &                     +xmlow*psi(iz,ir,2,jst)
                     spsilj = -dpsiz(iz,ir,2,jst)+dpsir(iz,ir,1,jst)
     &                     -xmup*psi(iz,ir,1,jst)
                     addls  = factor*(spsiui*spsiuj+spsili*spsilj)
     &                     -wrk1(iz,ir)
                     drhlsp(iz,ir,iq) = -addls + drhlsp(iz,ir,iq)
                   enddo
                 enddo
               endif      
            endif
            addfac = 2.0
         enddo
      enddo
c                                                 extrapolate to ir=0
      do iz=0,nz
        rhop(iz,0,1) = extrap(rhop(iz,1,1),rhop(iz,2,1),rhop(iz,3,1))
        rhop(iz,0,2) = extrap(rhop(iz,1,2),rhop(iz,2,2),rhop(iz,3,2))
        taup(iz,0,1) = extrap(taup(iz,1,1),taup(iz,2,1),taup(iz,3,1))
        taup(iz,0,2) = extrap(taup(iz,1,2),taup(iz,2,2),taup(iz,3,2))
        drhlsp(iz,0,1) = extrap(drhlsp(iz,1,1),drhlsp(iz,2,1),
     &                          drhlsp(iz,3,1))
        drhlsp(iz,0,2) = extrap(drhlsp(iz,1,2),drhlsp(iz,2,2),
     &                          drhlsp(iz,3,2))
      enddo
c
c     2. step: shift the main density by the response densities,
c              save the main densities on the response fields as
c              workspace
c
      rhopno = wfovlp(rhop,rhop)
      rhopov = wfovlp(rhop,rho)
      rhono  = wfovlp(rho,rho)
      varrho = rhopno/rhono
c      varrho = wfovlp(rhop,rhop)/wfovlp(rho,rho)
      taupno = wfovlp(taup,taup)
      taupov = wfovlp(taup,tau)
      tauno  = wfovlp(tau,tau)
      vartau = taupno/tauno
c      vartau = wfovlp(taup,taup)/wfovlp(tau,tau)
      stepr  = sqrt(min(stepq0/varrho,stepq0/vartau))
      write(7,'(a,3(1pg12.4))')
     &   ' HCOMMP: stepr,varrho,vartau=',stepr,varrho,vartau
      write(7,'(a/5(1pg12.4))')
     &  ' N(rho)  N(rhop)  Overlap   varrho  difference ',
     &   rhono,rhopno,rhopov,varrho,(rhono+rhopno-2.0*rhopov)/rhono
      write(7,'(a/5(1pg12.4))')
     &  ' N(tau)  N(taup)  Overlap   vartau  difference ',
     &   tauno,taupno,taupov,vartau,(tauno+taupno-2.0*taupov)/tauno
      steph = 0.5*stepr
      do iq=1,2
      do ir=0,nr
      do iz=0,nz
        xsav   = rho(iz,ir,iq)
        x      = sqrt(xsav)
c        rho(iz,ir,iq) = (x+steph*rhop(iz,ir,iq)/x)**2
        rho(iz,ir,iq) = xsav+stepr*rhop(iz,ir,iq)
        rhop(iz,ir,iq) = xsav
        xsav   = tau(iz,ir,iq)
        tau(iz,ir,iq) = xsav+stepr*taup(iz,ir,iq)
        taup(iz,ir,iq) = xsav
        xsav   = drhols(iz,ir,iq)
        drhols(iz,ir,iq) = xsav+stepr*drhlsp(iz,ir,iq)
        drhlsp(iz,ir,iq) = xsav
      enddo
      enddo
      enddo
      do ir=0,nr
      do iz=0,nz
        crho(iz,ir) = rho(iz,ir,1)
      enddo
      enddo
c
c     test new density
c
      rhono  = wfovlp(rhop,rhop)
      rhopov = wfovlp(rhop,rho)
      rhopno = wfovlp(rho,rho)
      tauno  = wfovlp(taup,taup)
      taupov = wfovlp(taup,tau)
      taupno = wfovlp(tau,tau)
c      stepr  = sqrt(min(1.0E-4/varrho,1.0E-4/vartau))
      write(7,'(a/a/5(1pg12.4))')
     &  ' HCOMMP:',
     &  ' N(rho)  N(rhop)  Overlap   angle ',
     &   rhono,rhopno,rhopov,one-rhopov*rhopov/(rhono*rhopno)
      write(7,'(a/5(1pg12.4))')
     &  ' N(tau)  N(taup)  Overlap   angle ',
     &   tauno,taupno,taupov,one-taupov*taupov/(tauno*taupno)
c
c     3. step: compute Hamiltonian for shifted densities
c              and deduce matrix elements of response from that,
c              result is stored on 'qopout'
c
ctest      call testpr(rho,0,10,0,4,' rho')
ctest      call testpr(tau,0,10,0,4,' tau')
ctest      call testpr(pot,0,10,0,4,' pot')
ctest      call testpr(bmass,0,10,0,4,' bmass')
      call mfhmlt(0)
ctest      call testpr(pot,0,10,0,4,' pot')
ctest      call testpr(bmass,0,10,0,4,' bmass')
c      write(7,'(a)') ' after 1. MFHMLT in HCOMMP: ist, spenrg'
      do ist=1,nst(1)+nst(2)
        sphc(ist) = spenrg(ist)                        ! save s.p. energy
c        write(7,'(1x,i5,1pg12.4)') ist,spenrg(ist)
      enddo
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         vi     = sqrt(vsqf(ist))
         ui     = sqrt(one-vi*vi)
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-eferm
         equasi = sqrt(edifi*edifi+deltai*deltai)
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            eferm  = fer(1)
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
            eferm  = fer(2)
            jstmax = nst(1)+nst(2)
         endif
         call hpsi(psi(0,0,1,ist),psin,ist)
         pairfa = -pairfc(iq)
         do is=1,2
         do ir=0,nr
         do iz=0,nz
          psip(iz,ir,is) = pairfa*paircp(iz,ir,iq)*psi(iz,ir,is,ist)
         enddo
         enddo
         enddo
         do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               vj     = sqrt(vsqf(jst))
               uj     = sqrt(one-vj*vj)
               deltaj = deltaf(jst)*wtpair(jst)
               edifj  = spenrg(jst)-eferm
               equasj = sqrt(edifj*edifj+deltaj*deltaj)
c               qopout(ipostr) = popin(ipostr)*(equasi+equasj) 
               qopout(ipostr) = popin(ipostr)*e2quas(ipostr)
     &              +(ui*vj+vi*uj)*wfovlp(psi(0,0,1,jst),psin)/stepr
c     &              +(ui*uj-vi*vj)*wfovlp(psi(0,0,1,jst),psip)/stepr
c              if(jst.eq.ist) then
c                qopout(ipostr) = (spenrg(ist)-sphc(ist))/stepr
c                qopout(ipostr) = (           -sphc(ist))/stepr
c    &                           + qopout(ipostr)
c              endif
            endif
         enddo
      enddo
c
c     4. step: restore densities and mean field Hamiltonian
c
c      write(7,'(a)') ' after 3.step in HCOMMP: ist, spenrg'
      do ist=1,nst(1)+nst(2)
        spenrg(ist) = sphc(ist)                        ! restore s.p. energy
c        write(7,'(1x,i5,1pg12.4)') ist,spenrg(ist)
      enddo
      do iq=1,2
      do ir=0,nr
      do iz=0,nz
        rho(iz,ir,iq) = rhop(iz,ir,iq)
        tau(iz,ir,iq) = taup(iz,ir,iq)
        drhols(iz,ir,iq) = drhlsp(iz,ir,iq)
      enddo
      enddo
      enddo
      do ir=0,nr
      do iz=0,nz
        crho(iz,ir) = rho(iz,ir,1)
      enddo
      enddo
      call mfhmlt(0)
c
      return
      end


