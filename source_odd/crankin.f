C-----crankin-----------------------------------------------------------
c
      subroutine crankin
c
c     saves wavefunctions
c     and computes Inglis cranking mass for quadrupole motion
c     by building derivative from previously save wavefunction
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
      logical tderv              ! checks if derivative is available
c
c      real*4 e2quas(ktriac)     ! storage for 2-quasiparticle energies
      double precision betain2(ktriac)
c
      data iqtest/0/            ! switch to test given Q-operators
      data iptest/0/            ! switch to test 'hcomP'
c      data itstcm/0/            ! switch to test case c.m. motion
c      data itstpr/1/            ! switch to testprints of Inglis Cranking
c      data isavsp/1/            ! save s.p. properties on file for interpol.
c
c      data delb/0.0/
c      save beto2,beto2p,beto2n,beto3,ferolp,feroln
c
c     buffered SQRT as internal statement function
c
c      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
      if(kcrank.le.0) return
c
      call inie2q(0,1)                         ! 2QP energies on 'e2quas'
c
      call quadms(qresp,0,1)                   ! simple quadrupole mass
      call quadms(betain2,0,1)
c
      inquire(UNIT=98,OPENED=tderv)         ! checks status
      call cringl(presp)                    ! prepare 'presp'
      if(.NOT.tderv) return
c
c     first a few tests
c
      if(iqtest.eq.1) then
         do icase=1,6
c
c     with cases 
c     icase = 1...6
c     1: Q = 2z^2
c     2: Q = -r^2
c     3: Q = 2z^2-r^2
c     4: Q = z
c     5: Q = sigma_r

c     6: Q = p_r sigma_r
c           
            b4=0.5*t4
            call hcOT(icase)
            t1sav = t1
            t2sav = t2
            x1sav = x1
            x2sav = x2
            b4sav = b4
            b4psav = b4p
            call Oabtst(icase)
            t1 = 0.0
            t2 = 0.0
            x1 = 0.0
            x2 = 0.0
            b4 = 0.0
            b4p = 0.0
            call Oabtst(icase)
            t1 = t1sav
            t2 = t2sav
            x1 = x1sav
            x2 = x2sav
            b4 = b4sav
            b4p = b4psav
         enddo
         stop 'cranking test-stop'
      endif
      write(7,'(a,2(1pg12.4))') 
     &   ' CRANKN: P-response initialized  <P**2>=',opovlp(presp,presp)
c
c     optionally test 'hcommP'
c
      if(iptest.eq.1 .and. tderv) then
         call hcPT
c         stop ' after HCPT'
      endif
c
c     initialize Q-response with PRESP
c      
c      write(7,'(a)') ' ipostr   ist   jst    equasij  deltai  deltaj',
c     &   '  eferm   spei  spej'    
ctest      ipostr = 0
ctest      do ist=1,nst(1)+nst(2)
ctest        mist   = mpsi(ist)
ctest        npist  = nparz(1,ist)
ctest        if(ist.le.nst(1)) then
ctest          jstmax = nst(1)
ctest        else
ctest          jstmax = nst(1)+nst(2)
ctest        endif
ctest        do jst=ist,jstmax
ctest          if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
ctest            ipostr = 1 + ipostr
      do ipostr=1,iposmx
        qresp(ipostr) = presp(ipostr)/e2quas(ipostr)
      enddo
ctest          endif
ctest        enddo
ctest      enddo     
      write(7,'(a,2(1pg12.4))') 
     &   ' CRANKN: Q-response initialized  <Q**2>=',opovlp(qresp,qresp)
c      iposmx = ipostr
c      write(7,'(a,i10)') ' second IPOSTR=',ipostr
c
c     iterate    Q_new = Q - stepQ * ( [H,Q] - P ),
c     the commutator [H,Q] is computed in 'hcommQ' and returned on 'ptest'.
c     The collective mass <[Q,[H,Q]]> is also accumulated in 'hcommQ' and
c     returned on 'scmass'. It is printed as protocol of convergence.
c
      stepQ  =   0.05                !  this needs to become an input variable
      itmaxQ =   150                !  this needs to become an input variable
      epsQ   = 0.0001
      scmaso = 0.0
      do itq=1,itmaxQ
        call hcommQ(qresp,presp,betain2,ptest,0)
        if(scmass.gt.zero) pnmass = pnmass/scmass**2
        if(itq.eq.1) pnmas0 = pnmass                      ! save Inglis level
        do ipostr=1,iposmx
          qresp(ipostr) = qresp(ipostr) 
     &                   - stepQ*(ptest(ipostr)-presp(ipostr))
     &                     /e2quas(ipostr)                  ! to damp step
        enddo
        write(7,'(a,i4,a,4(1pg12.4))')
     &    ' LRiter=',itq,
     &    '  SCmass,diff.,PNmasses=',scmass,scmass-scmaso,pnmass,pnmas0
        if(itq.gt.10 .and. scmass.ne.zero) then
         if(abs(scmaso/scmass-one).lt.epsQ .or. scmass.lt.zero
     &       .or. scmass.gt.1.0E20) goto 999
c          continue
        elseif(scmass.eq.zero) then
          goto 999
        endif
        scmaso = scmass
      enddo
 999  continue
      write(7,'(a,i5)')
     &  ' mass iteration ended at itq=',itq
c
      return
      end




c-----crnkam-----------------------------------------------------------
c
      subroutine crnkam
c
c     and computes Inglis cranking mass for angular momentum
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
c
      data itstpr/0/            ! switch to testprints of Inglis Cranking
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
c
c     double loop to accumulate cranking mass for angular momentum
c
      sumj2  = zero
      sumj2e = zero
      sumtst = 0.0
      ampnms = 0.0
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         if(ist.le.nst(1)) then
            jstmax = nst(1)
            jstmin = 1
            eferm  = fer(1)
CR            efermo = ferolp
         else
            jstmax = nst(1)+nst(2)
            jstmin = nst(1)+1
            eferm  = fer(2)
CR            efermo = feroln
         endif
c
c     ang.mom. operator J+ on 'psip' and J- on 'psin'
c
c                                     first, the L+ and L-
         call nablaz_spin(psi(0,0,1,ist),psi1,nparz(1,ist),1,one)
         call nablar_spin(psi(0,0,1,ist),psi2,   mpsi(ist),1,one)

         do is=1,nsd
            mm        = mist+is-1
            npist  = nparz(is,ist)
            nparr     = 1+mod(mm,2)
cPK         call mulr(derr1(0,0,nparr),psi(0,0,is,ist),psi2(0,0,is),one,0)
cPK         call mulz(derz1(0,0,npist),psi(0,0,is,ist),psi1(0,0,is),one,0)
            ir = 0
            if(mm.eq.1) then
               do iz=0,nz
                  z     = (iz+sz-icm)*dz
                  psip(iz,ir,is) = 0.0
c              psin(iz,ir,is) = (z+z)*psi2(iz,ir,is)
                  psin(iz,ir,is) = -(z+z)*psi2(iz,ir,is)
               enddo
            else
               do iz=0,nz
                  z     = (iz+sz-icm)*dz
c              psip(iz,ir,is) = zero
                  psip(iz,ir,is) = z*psi2(iz,ir,is)
c              psin(iz,ir,is) = psip(iz,ir,is)
                  psin(iz,ir,is) = -psip(iz,ir,is)
               enddo
            endif
            do ir=1,nr
               r    = (ir+sr)*dr
               rinv = mm/r
               do iz=0,nz
                  z     = (iz+sz-icm)*dz
                  psip(iz,ir,is) = z*(psi2(iz,ir,is)-rinv
     &                 *psi(iz,ir,is,ist))
     &                 -r*psi1(iz,ir,is)
c              psin(iz,ir,is) = z*(wrk2(iz,ir)+rinv*psi(iz,ir,is,ist))
c     &                        -r*wrk4(iz,ir)
                  psin(iz,ir,is) = -z*(psi2(iz,ir,is)+rinv
     &                 *psi(iz,ir,is,ist))
     &                 +r*psi1(iz,ir,is)
               enddo
            enddo
         enddo
c                                          the S+ and S-
         do ir=0,nr
            do iz=0,nz
               psip(iz,ir,1) = psip(iz,ir,1)+psi(iz,ir,2,ist)
               psin(iz,ir,2) = psin(iz,ir,2)+psi(iz,ir,1,ist)
            enddo
         enddo
c                                   test       the S+ and S-
ct        do ir=0,nr
ct        do iz=0,nz
ct          psip(iz,ir,1) = -psi(iz,ir,2,ist)
ct          psip(iz,ir,2) = 0.0
ct          psin(iz,ir,1) = 0.0
ct          psin(iz,ir,2) = psi(iz,ir,1,ist)
ct        enddo
ct        enddo
c                                       auxiliary variables
         npist  = nparz(1,ist)
         vi    = sqrtb(vsqf(ist))
         ui    = sqrtb(one-vi*vi)
         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-eferm
         equasi = sqrt(edifi*edifi+deltai*deltai)
c                                       the overlaps
         tstovl = wfovlp(psip,psip)
ctest         write(7,'(a,i3,1pg12.4)') ' test: ist, <J**2>=',ist,tstovl
         sumtst = tstovl*vi*vi + sumtst
         do jst=jstmin,jstmax
            if((npist.eq.nparz(2,jst) .and.
     &           (mist+1.eq.mpsi(jst) .or. mist-1.eq.mpsi(jst)) )
     &           .or.
     &           (npist.eq.nparz(1,jst) .and.
     &           (mist.eq.0 .and. mpsi(jst).eq.0))            ) then
               vj    = sqrtb(vsqf(jst))
               uj    = sqrtb(one-vj*vj)
               deltaj = deltaf(jst)*wtpair(jst)
               edifj  = spenrg(jst)-eferm
               equasj = sqrt(edifj*edifj+deltaj*deltaj)
               if(mist.eq.0 .and. mpsi(jst).eq.0) then
ctest            if((ist.eq.2 .and. jst.eq.5) .or.
ctest     &         (ist.eq.5 .and. jst.eq.2)) 
ctest     &         write(7,'(a,2(2x,2(1pg12.4)))')
ctest     &        ' 2&5: upper,lower=',
ctest     &        wfovoc(psi(0,0,1,jst),psin(0,0,2)),
ctest     &        wfovoc(psi(0,0,2,jst),psin(0,0,1)),
ctest     &        wfovoc(psi(0,0,1,jst),psi(0,0,1,ist)),
ctest     &        wfovoc(psi(0,0,2,jst),psi(0,0,2,ist))
c            do ir=0,nr
c            do iz=0,nz
c              wrk2(iz,ir)   = psin(iz,ir,1)
c              psin(iz,ir,1) = psin(iz,ir,2)
c              psin(iz,ir,2) = wrk2(iz,ir)
c              psin(iz,ir,2) = -psav
c            enddo
c            enddo
                  pij    = +wfovoc(psi(0,0,2,jst),psin(0,0,1))
     &                     -wfovoc(psi(0,0,1,jst),psin(0,0,2))
                  pij    = pij*(vj*ui-uj*vi)
c            write(7,'(a,3(2i2,1x),2(1pg12.4))')
c     &     ' 0:ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
c     &       ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
c     &       pij*(vj*ui-uj*vi)
               elseif(mist+1.eq.mpsi(jst)) then
                  pij    = wfovlp(psi(0,0,1,jst),psip)
                  pij    = pij*(vj*ui-uj*vi)
c            write(7,'(a,3(2i2,1x),2(1pg12.4))')
c     &     ' +:ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
c     &       ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
c     &       pij*(vj*ui-uj*vi)
               elseif(mist-1.eq.mpsi(jst)) then
                  pij    = wfovlp(psi(0,0,1,jst),psin)
                  pij    = pij*(vj*ui-uj*vi)
c            write(7,'(a,3(2i2,1x),2(1pg12.4))')
c     &     ' -:ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
c     &       ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
c     &       pij*(vj*ui-uj*vi)
               endif
               if(itstpr.ne.0) 
     &              write(6,'(a,2i4,2f8.2,4(1pg12.4))')
     &              ' i,j,equasi/j,phase,pij=',
     &              ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
cPK               pij    = pij*(vj*ui-uj*vi)
               pij    = pij*pij
               sumj2e = pij/(equasi+equasj) + sumj2e
               sumj2  = pij + sumj2
               ampnms = pij/(equasi+equasj)**2
     &                  *(ui*ui-vi*vi+uj*uj-vj*vj)
     &                  + ampnms
            endif
         enddo
      enddo
      ammass = sumj2e/2.0
      amwid  = sumj2/2.0
      ampnms = ampnms/2.0/ammass**2
      write(7,'(a,1pg12.4)') ' test <J**2>=',sumtst
c
      return
      end








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
c     This version computes the variation of the response Hamiltonian
c     directly from the potential by finite differences.
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
      dimension potp(0:nzdc,0:nrdc,kqm),
     x          coulptp(0:nzdc,0:nrdc),
     x          bmassp(0:nzdc,0:nrdc,kqm),
     x          potlsp(0:nzdc,0:nrdc,kqm,2)
      equivalence (potp,dpsir(0,0,1,1))
      equivalence (coulptp,dpsir(0,0,1,3))
      equivalence (bmassp,dpsir(0,0,1,5))
      equivalence (potlsp,dpsiz)
c
      data stepq0/0.000001/          ! determines step size --> 'stepr'
c
c     internal statement function for extrapolation to r=0
      data trhalf,trfift,oneten/1.5,0.6,0.1/
      extrap(f1,f2,f3) = trhalf*f1-trfift*f2+oneten*f3
c
c-----------------------------------------------------------------------
c
      if(ifls.ne.0) write(*,*) 
     &   ' IFLS=0 needed for linear response'
c      if(ifls.ne.0) stop ' IFLS=0 needed for linear response'
      if(ipair.ne.5 .and. ipair.ne.9 .and. ipair.ne.11) 
     &    write(*,*) ' IPAIR=5 or 9 needed for linear P response'
c     &    stop ' IPAIR=5 or 9 needed for linear P response'
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
         addfac = 2.0d0   
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
            addfac = 4.0 ! 2.0
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
c        x      = sqrt(xsav)
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
c              and deduce response Hamiltonian from that.
c              Response fields reside on potential fields.
c
      do ist=1,nst(1)+nst(2)
        sphc(ist) = spenrg(ist)                        ! save s.p. energy
        spenrg(ist) = 0.0d0
      enddo
      do ir=0,nr
      do iz=0,nz
        coulptp(iz,ir) = coulpt(iz,ir)
        potp(iz,ir,1) = pot(iz,ir,1)
        potp(iz,ir,2) = pot(iz,ir,2)
        bmassp(iz,ir,1) = bmass(iz,ir,1)
        bmassp(iz,ir,2) = bmass(iz,ir,2)
        potlsp(iz,ir,1,1) = potls(iz,ir,1,1)
        potlsp(iz,ir,2,1) = potls(iz,ir,2,1)
        potlsp(iz,ir,1,2) = potls(iz,ir,1,2)
        potlsp(iz,ir,2,2) = potls(iz,ir,2,2)
      enddo
      enddo
      call mfhmlt(0)
      sd = -1.0d0/stepr
      do ir=0,nr
      do iz=0,nz
        coulpt(iz,ir) = (coulptp(iz,ir)-coulpt(iz,ir))*sd
        pot(iz,ir,1) = (potp(iz,ir,1)-pot(iz,ir,1))*sd
        pot(iz,ir,2) = (potp(iz,ir,2)-pot(iz,ir,2))*sd
        bmass(iz,ir,1) = (bmassp(iz,ir,1)-bmass(iz,ir,1))*sd
        bmass(iz,ir,2) = (bmassp(iz,ir,2)-bmass(iz,ir,2))*sd
        potls(iz,ir,1,1) = (potlsp(iz,ir,1,1)-potls(iz,ir,1,1))*sd
        potls(iz,ir,2,1) = (potlsp(iz,ir,2,1)-potls(iz,ir,2,1))*sd
        potls(iz,ir,1,2) = (potlsp(iz,ir,1,2)-potls(iz,ir,1,2))*sd
        potls(iz,ir,2,2) = (potlsp(iz,ir,2,2)-potls(iz,ir,2,2))*sd
      enddo
      enddo

      call nablaz_scalar(bmass(0,0,1),dzbmas(0,0,1), 1,1,one)
      call nablaz_scalar(bmass(0,0,2),dzbmas(0,0,2), 1,1,one)
      call nablar_scalar(bmass(0,0,1),drbmas(0,0,1), 0,1,one)
      call nablar_scalar(bmass(0,0,2),drbmas(0,0,2), 0,1,one)

      call testpr(pot(0,0,1),0,9,0,9,'resp pot')
      call testpr(pot(0,0,2),0,9,0,9,'resp pot')
      call testpr(bmass(0,0,1),0,9,0,9,'resp pot')
      call testpr(bmass(0,0,2),0,9,0,9,'resp pot')
c
c     4. step: compute action of response Hamiltonian and
c              compose new Q operator
c
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         vi     = sqrt(vsqf(ist))
         ui     = sqrt(one-vi*vi)
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
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
               qopout(ipostr) = popin(ipostr)*e2quas(ipostr)
     &              +(ui*vj+vi*uj)*wfovlp(psi(0,0,1,jst),psin)
c     &              +(ui*uj-vi*vj)*wfovlp(psi(0,0,1,jst),psip)/stepr
            endif
         enddo
      enddo
c      write(6,'(a,i6)') ' HCOMMP: ipostr=',ipostr
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
        coulpt(iz,ir) = coulptp(iz,ir)
        pot(iz,ir,1) =  potp(iz,ir,1)
        pot(iz,ir,2) =  potp(iz,ir,2)
        bmass(iz,ir,1) =  bmassp(iz,ir,1)
        bmass(iz,ir,2) =  bmassp(iz,ir,2)
        potls(iz,ir,1,1) =  potlsp(iz,ir,1,1)
        potls(iz,ir,2,1) =  potlsp(iz,ir,2,1)
        potls(iz,ir,1,2) =  potlsp(iz,ir,1,2)
        potls(iz,ir,2,2) =  potlsp(iz,ir,2,2)
      enddo
      enddo
c      call mfhmlt(0)
      do ist=1,nst(1)+nst(2)
        spenrg(ist) = sphc(ist)
      enddo
c
      return
      end
c-----opovlp--------------------------------------------------------------
c
      function opovlp(op1in,op2in)
c
c     Overlap of 1ph-operators 'op1in' and 'op2in'.
c
      include 'paramr.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'crankin.inc'
c
      dimension op1in(ktriac),op2in(ktriac)
c
c-----------------------------------------------------------------------
c
c     accumulate matrix elements
c
      sum    = 0.0        
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         if(ist.le.nst(1)) then ! protons
            jstmax = nst(1)      
         else                   ! neutrons
            jstmax = nst(1)+nst(2)
         endif
         do jst=ist,jstmax
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               if(ist.eq.jst) then
                  sum = 2.0 * op1in(ipostr)*op2in(ipostr)   + sum
               else
                  sum = 4.0 * op1in(ipostr)*op2in(ipostr)    + sum
               endif
            endif
         enddo
      enddo
      opovlp = sum
c      
      return
      end 

