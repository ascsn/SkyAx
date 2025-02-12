c-----pairol------------ part of the hartree-fock package ---------------------
c                        changed to fit into the axial H.F.code
c
      subroutine pairol(e,gw,ph,nz,nmax,gp,eferm,delta,partnm,
     &                iter,ipair,eps,iab)
c
c     the include file 'params.inc' communicates the dimensioning of
c     the fields 'e' and 'gw' with the parameter 'ksl2'. it allows
c     furthermore to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c     Furthermore it carries the common block /CONSTN/ from which
c     the pseudo-constants ONE=1.0 and ZERO=0.0 are taken here.
c
      include 'paramr.inc'
c
c     * * * * * * * * * *     pairing iteration      * * * * * * * * * * * *
c
c     determines the fermi-energy 'eferm' and the gap 'delta'
c     such that the "gap-equation" and the "particle number condition"
c     are fulfilled simultanously. this is achieved by iteration
c     with a gradient scheme (newton's tangential formula)
c     where the gradient is evaluated by finite differences.
c
c     there are two schemes possible:
c     1. the "constant force" approach, where the pairing force is given
c                                       in 'gp'.
c     2. the "constant gap" approach,   where the gap is kept fixed and
c                                       is given in 'delta'.
c     there is only 'eferm' to be determined in the "constant gap"
c     approach. the gap-equation is used at the end to determine an
c     according force 'gp'.
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initila guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c      nz      = desired particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'npsid'
c                 and smaller than 46, for dimensioning reasons).
c      gp      = size of pairing-force, is input in case of
c                "constant force" approach.
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      delta   = pairing gap, is input in case of "constant gap" approach.
c      iter    = number of pairing iterations,
c                abs(iter) is the number of iterations
c      ipair   = switch for pairing type
c                 ipair=1 -->  "constant force" approach, using
c                              pairing force 'gp'
c                 ipair=0 -->  "constant gap" approach, using
c                              given gap 'delta'
c                 ipair=4 -->  thermal occupation with temperature
c                              given by 'delta'
c      eps     = end condition for pairing iteration;
c                iteration terminates if error in particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      gp      = size of pairing-force, is output in case of
c                "constant gap" approach.
c      delta   = pairing gap, is output in case of
c                "constant force" approach.
c      eferm   = fermi-energy.
c      sum     = particle number, from adding 'gw' with degeneracy 'ph'.
c
c     some auxiliary fields are:
c      ph      = array(1..46) of the multiplicities of the single
c                particle levels. array is ordered according to
c                the harmonic oscillator scheme, with rising shell
c                number and lowering angular momentum within each shell.
c      dparn   = array(1..3) of variation in the particle number as
c                consequence of 'deferm' and 'ddelta'.
c      dgap    = array(1..3) of variation in the gap equation as
c                consequence of 'deferm' and 'ddelta'.
c      gr      = double array of variation of particle number and of
c                gap equation, used for construction of the tangantial
c                step
c
      dimension e(npsid),gw(npsid),ph(npsid)
      dimension dparn(3),dgap(3),gr(2,2)
c
      data gr,dparn,dgap/10*0./
c
c     further numerical variables are:
c      deferm  = step in the fermi-energy, to evaluate
c                numerically the derivatives.
c      ddelta  = step in the gap, to evaluate numerically the
c                derivatives; used only in the "constant force" approach.
c      dampin  = damping factor for the tangential step
c      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
c
c     note: these parameters are adjusted for an application with
c           the nuclear shell model. it is assumed that the single
c           particle energies are given in units of [mev] and display
c           a reasonable nuclear spectrum.
c
      data dampin,deferm,ddelta/0.9,1.0,0.1/
      data efmax,dlmax,detmin,defmin,ddlmin/3.,.6,.1e-3,.1e-3,.1e-3/
cclust      data dampin,deferm,ddelta/0.9,0.01,0.1/
cclust      data efmax,dlmax,detmin,defmin,ddlmin/0.03,.6,.1e-3,.1e-3,.1e-3/
c
      data xmaxlw,x10lw/1.0e-30,1.0e-10/
c      data one,zero/1.0,0.0/
c
c-------------------------------------------------------------------------------
c
c      set switch variables
c
      iwarng = 0
c
      if(iab.ge.1)
     & write(6,'(/a,i4,3(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** pair entered with: iter=',iter,'  delta=',delta,
     &   '  eferm=',eferm,'  gp=',gp,
     &   '                        single particle energies:',
     &   (e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       isum   = 0
       do 120 it=1,nmax
       itsum  = it
       isum   = int(ph(it)+.000001) +isum
       if(isum+isum-nz.ge.0) goto 110
  120  continue
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = e(itsum)+deferm
       if(iab.ge.1) write(06,'(1x,a,i5,(a,g11.4),a,i5)')
     &   'initialistion: isum=',isum,'  eferm=',eferm,'  nz=',nz
c
      endif
c
      if(ipair.eq.1) then
c
c-------------------------------------------------------------------------------
c
c     case of constant pairing-force, given by 'gp'.
c      gap 'delta', fermi-energy 'eferm', and occupations 'gw'
c      are computed.
c
c-------------------------------------------------------------------------------
c
       delti  = delta
       if(abs(delti).lt.x10lw) delti = one
       do 200 it=1,abs(iter)
c
c       gap-equation & particle-number at 'eferm' & 'delti' (for k=1)
c       and at 'eferm' varied (k=2) and at 'delta' varied (k=3).
c
        do 210 ka=1,3
c                                          actual fermi-energy and gap
        k      = 4-ka
        if(k.eq.1) then
         elam   = eferm
         delact = delti
        elseif(k.eq.2) then
         elam   = eferm+deferm
         delact = delti
        else
         elam   = eferm
         delact = delti+ddelta
        endif
c                                         compute actual quasiparticle energies
c                                         and occupation-weights,
c                                         accumulate gap-equation and
c                                         particle number
        parnum = zero
        gapeq  = zero
         do 220 i=1,nmax
         equasi = sqrt((e(i)-elam)*(e(i)-elam)+delact*delact)
         gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
         gapeq  = ph(i)/equasi   + gapeq
         parnum = gw(i)*ph(i)    + parnum
  220    continue
c                                         store actual part.number and
c                                         gap-equation
        dgap(k)  = 0.25*gp*gapeq-one          !?? extra 1/2 because of WATE=2
        dparn(k) = parnum-nz
  210   continue
c
c       construction of the gradient-matrix and its inverse
c        ( determinant 'det' is limited in value to avoid singularities )
c
        gr(1,1) = ( dgap(2)- dgap(1))/deferm
        gr(1,2) = (dparn(2)-dparn(1))/deferm
        gr(2,1) = ( dgap(3)- dgap(1))/ddelta
        gr(2,2) = (dparn(3)-dparn(1))/ddelta
        det     = gr(1,1)*gr(2,2)-gr(1,2)*gr(2,1)
        det     = sign(max(abs(det),detmin),det)
        if(abs(det).lt.xmaxlw)  det    = sign(xmaxlw,det)
c
        if(iab.ge.1) then
         write(6,'(a,i4,3(a,g11.3)/a,3g11.3,a,3g11.3/a,4g11.3)')
     &    ' iteration at it=',it,' :  delti=',delti,'  eferm=',eferm,
     &    '  det=',det,
     &    ' dgap=',dgap,' ,  dparn=',dparn,
     &    '   gr=',gr
         write(6,'(a,10f7.4,10(/6x,10f7.4))')
     &    '   gw=',gw
        endif
c
        grsav   =  gr(1,1)
        gr(1,1) =  gr(2,2)/det
        gr(1,2) = -gr(1,2)/det
        gr(2,1) = -gr(2,1)/det
        gr(2,2) = grsav/det
c
c       tangential step  (including some step limiters)
c
        efold  = eferm
        delold = delti
        efstep = dampin*(gr(1,1)*dgap(1)+gr(2,1)*dparn(1))
        dlstep = dampin*(gr(1,2)*dgap(1)+gr(2,2)*dparn(1))
        eferm  = eferm-sign(min(abs(efstep),efmax),efstep)
        delti  = abs(delti-sign(min(abs(dlstep),dlmax),dlstep) )
        deferm = efold-eferm
        ddelta = delold-delti
        if(abs(deferm).lt.defmin) deferm = defmin
        if(abs(ddelta).lt.ddlmin) ddelta = ddlmin
c
        if(iab.eq.0)
     &   write(6,'(a,i4,2(a,g11.3))')
     &    ' pairing iteration it=',it,'  :  delti=',delti,
     &    '  eferm=',eferm
c
        if(abs(dgap(1)).lt.eps .and. abs(dparn(1)).lt.eps)   goto 299
  200   continue
        iwarng = 1
c
  299   continue
        delta  = delti
c
      else
c
c-------------------------------------------------------------------------------
c
c      case of constant paring-gap, given by 'delta'.
c      using now the safer secant step.
c
c-------------------------------------------------------------------------------
c
c      determine upper and lower bound for secant step
c
       itmaxp = abs(iter)
       eferm1 = eferm
       eferm2 = eferm
       parn1  = parnm(e,gw,ph,nmax,delta,eferm,ipair)
       parn2  = parn1
       If(parn1-nz.gt.zero) then
        emin = +1.0E30
        do 311 i=1,nmax
        emin = min(e(i),emin)
  311   continue
        deferm = (eferm-emin)/itmaxp
        if(iab.ge.1)
     &   write(6,'(a,3(1pg12.4))')
     &   ' search brackets: parn2,eferm,emin=',parn2,eferm,emin
        do 313 it=1,itmaxp
        eferm1 = eferm1-deferm
        parn1  = parnm(e,gw,ph,nmax,delta,eferm1,ipair)
        if(iab.ge.1)
     &   write(6,'(a,i5,3(1pg12.4))')
     &   ' search brackets: it,parn1,eferm1=',it,parn1,eferm1
        if(parn1-nz.lt.zero) goto 319
        eferm2 = eferm1
  313   continue
        stop ' particle number not embraced in PAIR'
       elseif(parn1-nz.lt.zero) then
        emax = -1.0E30
        do 312 i=1,nmax
        emax = max(e(i),emax)
  312   continue
        deferm = (emax-eferm)/itmaxp
        if(iab.ge.1)
     &   write(6,'(a,3(1pg12.4))')
     &   ' search brackets: parn1,eferm,emax=',parn1,eferm,emax
        do 314 it=1,itmaxp
        eferm2 = eferm2+deferm
        parn2  = parnm(e,gw,ph,nmax,delta,eferm2,ipair)
        if(iab.ge.1)
     &   write(6,'(a,i5,3(1pg12.4))')
     &   ' search brackets: it,parn2,eferm2=',it,parn2,eferm2
        if(parn2-nz.gt.zero) goto 319
        eferm1 = eferm2
  314   continue
        stop ' particle number not embraced in PAIR'
       else
       if(iab.ge.1) write(6,'(a)') ' solution immediately found '
        goto 399            ! solution immediately found. jump to end
       endif
  319  continue             !  particle number 'NZ' between 'eferm1 <--> 2'
       if(iab.ge.1) write(6,'(a)') ' brackets found '
c
c      now perform the secant search
c
       do it=1,itmaxp
         eferm3 = eferm1 + (nz-parn1)*(eferm2-eferm1)/(parn2-parn1)
         parn3  = parnm(e,gw,ph,nmax,delta,eferm3,ipair)
        if(iab.ge.2)
     &    write(6,'(a,i4,2(a,g11.3))')
     &     ' secant iteration it=',it,'  :  eferm=',eferm3,
     &     '  dparn(1)=',parn3-nz
         if(abs(parn3-nz).lt.eps) then
           goto 399
         elseif(parn3-nz.gt.zero) then
           eferm2 = eferm3
           parn2  = parn3
         else
           eferm1 = eferm3
           parn1  = parn3
         endif
       enddo
c
c      the secant step did not succeed. use bisection
c
       do it=1,itmaxp
         eferm3 = 0.5*(eferm2+eferm1)
         parn3  = parnm(e,gw,ph,nmax,delta,eferm3,ipair)
         if(iab.ge.0)
     &     write(6,'(a,i4,2(a,g11.3))')
     &      ' bisection it=',it,'  :  eferm=',eferm3,
     &      '  dparn(1)=',parn3-nz
         if(abs(parn3-nz).lt.eps) then
           goto 399
         elseif(parn3-nz.lt.zero) then
           eferm1 = eferm3
           parn1  = parn3
         elseif(parn3-nz.gt.zero) then
           eferm2 = eferm3
           parn2  = parn3
         endif
       enddo
       iwarng = 1             ! still not converged. set a flag
c
  399  continue
       eferm = eferm3
       dparn(1) = parn3 - nz
       dgap(1)  = zero
       partnm   = parn3
c
c      determine force 'gp' according to given 'delta' and the
c      equilibrium reached.
c
       if(ipair.eq.0) then
         gapeq  = zero
         do i=1,nmax
           gapeq = ph(i)/sqrt((e(i)-eferm)*(e(i)-eferm)+delta*delta)
     *            + gapeq
         enddo 
         gp     = 4.0/gapeq                       !?? account for WATE=2.0
       endif
c
c     end of big switch between cases
c
      endif
c
c-------------------------------------------------------------------------------
c
c     finally the gap-equation and particle-number are checked and printed
c
c-------------------------------------------------------------------------------
c
c
c     in case of incomplete convergence print a warning
c
      if(iab.ge.-1 .and. iwarng.eq.1)
     & write(6,'(2(a,g11.3),2(a,i3))')
     &  ' ---> no convergence in pair! err(gap-eq)=',dgap(1),
     &  '  err(part.nr)=',dparn(1),'  nz=',nz,'  nmzx=',nmax
c
c     check gap equation
c
      if(iab.ge.0) then
        if(ipair.ne.4) then
          gphalf = gp*0.5
          del2   = max(xmaxlw,delta*delta)
          gapeq  = zero
          partnm = zero
          do i=1,nmax
            partnm = gw(i)*ph(i) + partnm
            gapeq  = gphalf*ph(i)/sqrt(del2+(e(i)-eferm)**2) + gapeq
          enddo
        else
          gapeq = zero
        endif
c
        write(6,'(a,i4,2(a,g13.5))')
     &  ' pair finished with ',it,' iterations: gap-eq.=',gapeq,
     &  '  part.nr=',partnm
      endif
c
      return
      end
c-----parnm------------ part of the hartree-fock package ----------------------c                        changed to fit into the axial H.F.code
c
      function parnm(e,gw,ph,nmax,delta,elam,ipair)
c
c     the include file 'params.inc' communicates the precision
c     transfering a line 'implicit double precision (a-h,o-z)' or not.
c     Furthermore it carries the common block /CONSTN/ from which
c     the pseudo-constant ONE=1.0 is taken here.
c
      include 'paramr.inc'
c
c     the function computes the particle number in a pairing or
c     temperature distribution. The parameters are:
c       e      = array of s.p. energies
c       gw     = array of occupation weights (is workspace and output!)
c       ph     = array of degeneracy factors
c       nmax   = number of states in the above three arrays
c       delta  = pairing gap or temperature
c       elam   = Fermi energy
c       ipair  = switch to pairing (0) or temperature (4).
c
      dimension e(npsid),gw(npsid),ph(npsid)
c
c----------------------------------------------------------------------
c

        sum    = zero
        if(ipair.eq.0) then
         do 330 i=1,nmax
         equasi = sqrt( (e(i)-elam)*(e(i)-elam) + delta*delta )
         gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
         sum    = gw(i)*ph(i)    + sum
  330    continue
        elseif(ipair.eq.4) then            ! case of temperature
         do 331 i=1,nmax
         equasi = (e(i)-elam)/delta
         if(equasi.gt.one*50.0) then
          gw(i) = 0.0
         elseif(equasi.lt.-one*50.0) then
          gw(i)  = 1.0
         else
          gw(i)  = 1.0/(1.0+exp(equasi))
         endif
         sum    = gw(i)*ph(i)    + sum
  331    continue
        endif
c
        parnm = sum
c
        return
        end

c-----hfbpot------------------------------------------------------
c
      subroutine hfbpot(work)
c
      include 'paramr.inc'
c
c     computes the pair density
c
c       sum_n  u_n*v_n*|phi_n(r)|^2  ,
c
c     subsequently the pair potential on 'pairpt' from /wavfun/,
c     and finally the dynamic gaps on 'delatf' from /shells/
c     as expectation value of 'pairpt'.
c
c     'ipair' is used here to indicate the switch to surface pairing
c             5,6 --> volume pairing
c             7,8 --> surface pairing
c
c
      include 'wavefr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'levels.inc'
      include 'observ.inc'
      include 'constr.inc'
c
      
      REAL(8) ::  Ecrit,cfac,npower
      REAL(8) ::  epairtestp,epairtestn,epairtest
      REAL(8) ::  buff(2)
      REAL(8) ::  add_delta(2)


      dimension  work(0:nzd,0:nrd)
      dimension  workn(0:nzd,0:nrd)
      dimension  workp(0:nzd,0:nrd)


c
c      data rho0pr/0.16/               ! preliminary standard for surface pairing
c
c--------------------------------------------------------------------
c
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
c
c     accumulate new pair-density
c
c      smallp = 0.000001
!      WRITE(6,'(a,200i4)') 'ITPARTN-HFB:',itpartn(1:nst(1)+nst(2))
      smallp = 1D-40
      do jst=1,nst(1)+nst(2)
        iq    = icharg(jst)
        IF(jst==istblo(iq) .OR. jst==istcblock(iq)) CYCLE
c        if(jst.eq.istblo)then
c        v2  = 2*vsqf(jst)
c        write(55,'(1(1pg12.4))')
c     &  v2  
c        else
        v2  = vsqf(jst)
c        endif
        uv  = wtpair(jst)*sqrt(max(v2-v2*v2,smallp))*0.5

coddAnfang  JE
        if(isTSym.eq.0) then
         call sdenti(psi(0,0,1,jst),paircr(0,0,iq),uv,jst)
        else
         if(mpsi(jst).lt.0)then
          uv = (-1.0)*uv
         endif
          call cpairpot(psi(0,0,1,jst),psi(0,0,1,itpartn(jst)),
     &     paircr(0,0,iq),uv,jst)
        endif
coddEnde
      enddo


c
c     now multiply with strength to obtain local pair-potential
c
      IF(ioptsk==13) THEN
        CALL fy_hfbpot()
      ELSE
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



c
c     compute pairing energy 
c
c      do ir=0,nr
c      do iz=0,nz     
c      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)
c      enddo
c      enddo
c      
c      epairtestp=f1int(work)  !0.25 oder 0.5,was ist v0 ?
c      
c      do ir=0,nr
c      do iz=0,nz     
c      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)
c      enddo
c      enddo
c      
c      epairtestn = f1int(work)  !0.25 oder 0.5,was ist v0 ?c
c
c      
c      epairtest= epairtestp+ epairtestn
cc
c
c
ccccccccccccccccccccccccccNachbesserung fuer DDDI 18.02.08ccccccccccccccccccc


c
c     compute pairing energy 
c

      if(ipair.le.6 .or. ipair.eq.9) then
      do ir=0,nr
      do iz=0,nz     
      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)
      enddo
      enddo
      
      epairtestp=f1int(work)  
      
      do ir=0,nr
      do iz=0,nz     
      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)
      enddo
      enddo
      
      epairtestn = f1int(work)  

      

      else
      
      

      do ir=0,nr
      do iz=0,nz     
      rhofac = one-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr
      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)* rhofac
      enddo
      enddo
      
      epairtestp=f1int(work)  
      
      do ir=0,nr
      do iz=0,nz     
      rhofac = one-(rho(iz,ir,1)+rho(iz,ir,2))/rho0pr
      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)* rhofac
      enddo
      enddo
      
      epairtestn = f1int(work)  
      
      endif

      epairtest= epairtestp+ epairtestn

      END IF      ! end 'ioptsk' switch
c

c     finally compute the actual gaps as s.p. expectation values with
c     the pair potential
c
      IF(suv_mue(1)==0D0) THEN
        add_delta(1)=0D0
      ELSE
        add_delta(1) = -suv_lambda(1)
     &                 +2D0*suv_mue(1)*(sum_uv(1)-suv_wanted(1))
        suv_lambda(1) = suv_lambda(1)
     &                 -2D0*c0_suv*suv_mue(1)*(sum_uv(1)-suv_wanted(1))
      END IF
      IF(suv_mue(2)==0D0) THEN
        add_delta(2)=0D0
      ELSE
        add_delta(2) = -suv_lambda(2)
     &                 +2D0*suv_mue(2)*(sum_uv(2)-suv_wanted(2))
        suv_lambda(2) = suv_lambda(2)
     &                 -2D0*c0_suv*suv_mue(2)*(sum_uv(2)-suv_wanted(2))
      END IF
      do jst=1,nst(1)+nst(2)
        iq     = icharg(jst)
        do ir=0,nr
        do iz=0,nz
          if(isTSym.eq.0)then
          wfa    = psi(iz,ir,1,jst)
          wfb    = psi(iz,ir,2,jst)
          work(iz,ir)  = (wfa*wfa+wfb*wfb)*pairpt(iz,ir,iq)
         else
codd
          wfa1    = psi(iz,ir,1,jst)
          wfb1    = psi(iz,ir,2,jst)
          wfa2    = psi(iz,ir,1,itpartn(jst))
          wfb2    = psi(iz,ir,2,itpartn(jst))
          work(iz,ir) = (wfa1*wfb2-wfb1*wfa2)*pairpt(iz,ir,iq)
          endif
codd
        enddo
        enddo
        IF(jst==istblo(iq) .OR. jst==istcblock(iq)) THEN
          deltaf(jst) = 0D0
        ELSE
          deltaf(jst) = f1int(work)
          if(mpsi(jst).lt.0) deltaf(jst) = -deltaf(jst)
          deltaf(jst) = deltaf(jst)+add_delta(iq)
        END IF

      enddo
c
      if(itest.gt.0) then
        write(6,'(a)') ' dynamics gaps for protons, directly:'
        write(6,'(5(1pg12.4))') (deltaf(i),i=1,nst(1))
        write(6,'(a)') ' dynamics gaps for neutrons, directly:'
        write(6,'(5(1pg12.4))') (deltaf(i),i=nst(1)+1,nst(1)+nst(2))
      endif
c


c
c     stabilized pairing gap
c
c      cfac=-0.5D0
      IF(ioptsk==13 .AND. ABS(pairstab)>0D0) THEN
        buff(1)=abs(pairstab/epair(1))
        buff(2)=abs(pairstab/epair(2))
!        WRITE(*,*) 'pairstab,epair:',pairstab,epair
        do jst=1,nst(1)+nst(2)
           iq     = icharg(jst)
           deltaf(jst)=
     &     SIGN(MIN(ABS(deltaf(jst)*(1.0D0 + buff(iq)**2)),5D0),
     &          deltaf(jst))
        enddo
!        write(6,'(a)') ' dynamics gaps for protons, after:'
!        write(6,'(5(1pg12.4))') (deltaf(i),i=1,nst(1))
      ELSE If(ABS(pairstab)>0D0) THEN

      Ecrit=1.0D0
      npower=2.0D0
      
      buff(1)=abs(Ecrit/epairtestp)
      buff(2)=abs(Ecrit/epairtestn)

     
      do jst=1,nst(1)+nst(2)
         iq     = icharg(jst)
         deltaf(jst)=deltaf(jst)*(1.0D0 - 
     &    (npower-1.0D0)*pairstab*buff(iq)**npower)

      enddo      

      END IF

      IF(ABS(pairstab)==0D0) deltaf(jst) = MIN(deltaf(jst),-1D-6)

      return
      end
c-----pairdn------------ variant for the axial SHF package ---------------------
c
      subroutine pairdn(e,delta,gw,ph,pairgw,partnm,eferm,epair,avdelt,
     &                  avdeltv2,efadd,cutfac,nmax,ndim,iter,ipair,
     &                  iwfoff,nblock,ntcblock,eps,iab,iq,tpairspace)
!     &                  avg,efadd,cutfac,nmax,ndim,iter,ipair,iwfoff,
!     &                  nblock,ntcblock,eps,iab,iq)
c
C
CR 31-10-1998 Brent algo. for determination of eferm
C
      include 'paramr.inc'
c
c     the include file 'params.inc' allows
c     to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c     Furthermore it carries the common block /CONSTN/ from which
c     the pseudo-constants ONE=1.0 and ZERO=0.0 are taken here.
c
c
c     * * * * * * * * * *   pairing with full interaction   * * * * * * * * * *
c
c     determines the fermi-energy 'eferm'
c     such that the "particle number condition"
c     is fulfilled. this is achieved by iteration.
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      delta   = array(1..nmax) of gap parameters.
c      gw      = array(1..nmax) of occupation probabilities,
c                input is used as initial guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c                This multiplicities INCLUDES here the spin degeneracy,
c                e.g. a 1s1/2 state has ph = 2 .
c      pairgw  = space for dynamical cut-off factor
c      partnm  = desired particle number
c      eferm   = fermi-energy, input is used as initial guess;
c                if eferm=0.0 is input then a suitable initial guess
c                is computed from input 'gw'.
c      epair   = pairing energy
c      avdelt  = average uv gap
c      avdeltv2  = average v**2 gap
c      avg     = average pairing force
c      efadd   = pairing band above 'eferm', in case of dynamical pair space
c                see routine  CUTOFF
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iter    = number of pairing iterations,
c                abs(iter) is the number of iterations
c      iwfoff  = offset for addressing wavefunctions
c      nblock  = address of blocked state (none if 0)
c      ntcblock= address of conjugated blocked state (none if 0)
c      ipair   = pairing case, switches to LN if 'ipair.eq.9'
c      eps     = end condition for pairing iteration;
c                iteration terminates if error in particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      eferm   = fermi-energy.
c
      LOGICAL,INTENT(IN) :: tpairspace
      dimension e(ndim),gw(ndim),ph(ndim),delta(ndim),pairgw(ndim)
      dimension vln(npsid)                                   ! workspace for LN
c
c     several numerical parameters, here for nuclei
c
C      data stepl /0.2/          ! stepsize for iteration of Fermi energy
      data delefe/20.0/          ! difference of Fermi energy for derivative
      data ifcut/1/             ! switch for flexibl cutoff
C      data slpmin/2.0/          ! minimum slope in dN/dE
      data xlam2/0.0/           ! to make LN parameter 'lambda_2' resident
      data xl2stp/1.0/         ! trial step for LN parameter
      data dxllim/0.5/          ! limiting 'eferm'-step in case of LN
      data dx2lim/0.1/          ! limiting 'xlam2'-step in case of LN
      data detmin/1.0E-4/       ! lower limit for determinant in LN step
      data dampin/0.5/          ! damping of LN step
c
      data half,xsmall/0.5,1.0E-20/
c
      data ifirst/0/
C
c-----------------------------------------------------------------------
c
      RTMP=0.0 ! tbcs is not yet used
C
      if(iab.gt.-1) then
        write(6,'(/a,i4,1(a,f9.4)/a,10(/1x,10f8.3))')
     &   ' *** pairing with dynamic gap entered with: iter=',iter,
     &   '  eferm=',eferm,
     &   '  single particle energies:',
     &   (e(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') '  delta:',(delta(i),i=1,nmax)
        write(6,'(a,10(/1x,10f8.4))') '    ph:',(ph(i),i=1,nmax)
      endif

c
c                               preset pairing weights
c
!      WRITE(6,*) 'cutfac,eferm=',cutfac,eferm
      if(eferm.ne.zero .and. ifcut.eq.1 .AND. cutfac>0D0) then
       actnum = partnm + cutfac*partnm**0.6666667      ! active pairing space  1.65
       delcut = 0.1
       efadd  = 20.0/partnm**0.33333333
       call cutoff(e,pairgw,ph,actnum,nmax,eferm,efadd,delcut,
     &             tstnum,iter,eps,-2,iq)
c     &             tstnum,iter,eps,iab+1)
!         WRITE(6,'(a,5(1pg13.5))') 'CUTOFF: actnum,efadd,tstnum=',
!     &        actnum,efadd,tstnum,cutfac,partnm
      else if(eferm.ne.0D0 .and. ifcut.eq.1 .AND. cutfac<0D0) then
        efadd=-cutfac
        efcut=eferm+efadd
        acc = 0D0
        do n=1,nmax
          pairgw(n) = 1D0/(1D0+EXP((e(n)-efcut)/(0.100D0*efadd)))
          acc = pairgw(n)*ph(n) + acc
        enddo
!        WRITE(6,'(a,2i4,1pg13.5)') 'pairing space=',iter,nmax,acc
!        WRITE(6,*) 's.p.energies:'
!        WRITE(6,'(5(1pg13.5))') e(1:nmax)
!        WRITE(6,*) 'pairing space weights:'
!       WRITE(6,'(5(1pg13.5))') pairgw(1:nmax)
      else
       do n=1,nmax
        pairgw(n) = 1.0
       enddo
      endif
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(eferm.eq.zero) then
       partsm = zero
       do it=1,nmax
        itsum  = it
cspher        partsm = 2.0*ph(it) + partsm
        partsm = ph(it) + partsm
        if(partsm-partnm.ge.zero) goto 110
       enddo
       stop ' not enough states too reach particle number'
  110  continue
       eferm  = 0.5*(e(itsum)+e(itsum+1))
       if(iab.ge.1) write(06,'(1x,a,2g11.4,a,g11.4)')
     &   'initialization: partsm,partnm=',partsm,partnm,' eferm=',eferm
       do n=1,nmax
         if(n.le.itsum) then
           gw(n) = one
         else
           gw(n) = zero
         endif
       enddo
       xlam2 = 0.0
c
      endif
c
c     the coupled gradient iteration for two cases, switched by 'ipair'
c
      if(ipair.ne.9.and.ipair.ne.11) then
c
c       standard case where p.n. is iterated, misusing 'gw' as workspace
c
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CR 03-04-1998
C     CR 31-10-1998
C
C

         EFERM = RBRENT(E, DELTA, PH, GW, PAIRGW, RTMP, PARTNM, 
     &        NDIM, NMAX, nblock,ntcblock, IAB,iq,tpairspace)
         PARNUM = PARNMT(E, DELTA, PH, GW, PAIRGW, RTMP, 
     &        EFERM, NDIM, NMAX, nblock, ntcblock,iq,tpairspace)    

         IT=0.0

C     
c
c       the LN case where p.n. and p.n. fluctuations are iterated
c     
      else
         iftest  = iab     !+ 1
         iterh   = iter/2
         iterp   = iter/10
         do it=1,iter
            if(iftest.gt.0)
     &           write(6,'(a,i5,a)')  ' iteration nr.',it,' :'
c     
c         initial call to pairing equations
c
            npart  = partnm+0.000001
            if(npart.eq.12 .and. ifirst.eq.1) then
              ifirst = 0
              write(6,'(a)')  
     &          '  lambda     lambda_2            f1           f2    '
              do il2=0,20
                tlam2 = -0.5 + il2*0.05
                do il1=0,20
                  tlam  = -25.0 + il1*0.5
                    call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &                          tlam,tlam2,gw,gw,iftest, f1,f2)
                  write(6,'(1x,2f8.2,2(1pg12.4))')
     &              tlam,tlam2,f1,f2,f1**2+5*f2**2
                enddo
              enddo
              stop
            endif
            call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &           eferm,xlam2,gw,gw,iftest, f1,f2)
            if(iftest.gt.0)
     &           write(6,'(a,2g12.3,a,2g12.3)')
     &           ' xlam,xlam2=',eferm,xlam2,'   f1,f2=',f1,f2
c
c         test termination
c
            if(it.gt.1 .and. abs(f1).lt.eps 
     &           .and. abs(f2).lt.eps) goto 1099
c     
c     derivatives by finite differences
c
            call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &           eferm+delefe,xlam2,gw,vln,iftest, f1xl,f2xl)
            call lm2def(e,delta,ph,pairgw,1,nmax,iwfoff,npart,
     &           eferm,xlam2+xl2stp,gw,vln,iftest, f1xl2,f2xl2)
         if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     &        ' LN: iter,xlam,xlam2,f1,f2=',it,eferm,xlam2,f1,f2
            dlf1   = (f1xl-f1)/delefe
            dlf2   = (f2xl-f2)/delefe
            dl2f1  = (f1xl2-f1)/xl2stp
            dl2f2  = (f2xl2-f2)/xl2stp
            if(iftest.ge.1) 
     &           write(6,'(a,4g12.4)')
     &           ' dlf1,dlf2,dl2f1,dl2f2=',dlf1,dlf2,dl2f1,dl2f2
c
c         propagation to improved solution by inverse matrix
c
            if(it.lt.iterh) then
              det    = (dlf1*dl2f2-dlf2*dl2f1) !   - sign ?
              det    = sign(max(abs(det),detmin),det) ! lower limit for 'det'
              dxlam  = (f2*dl2f1-f1*dl2f2)/det
              dxlam2 = (f1*dlf2-f2*dlf1)/det
              if(iftest.ge.0) 
     &           write(6,'(a,4g12.4)')
     &           ' dir: dxlam,dxlam2=',dxlam,dxlam2
c
c       propagation by minimzation of least squares .
c       this as such is equivalent to direct gradient,
c       but here we can employ the "damping" by systematically
c       enhancing the diagonal elements (see Bevington step).
c
            else
              dchi1  = f1*dlf1 + 5.0*f2*dlf2
              dchi2  = f1*dl2f1 + 5.0*f2*dl2f2
              dchi11 = dlf1*dlf1 + 5.0*dlf2*dlf2
              dchi12 = dl2f1*dlf1 + 5.0*dl2f2*dlf2
              dchi22 = dl2f1*dl2f1 + 5.0*dl2f2*dl2f2
              inhanc = (it-iterh)/iterp+1
              enhanc = 2.0**inhanc
              dchi11 = dchi11*enhanc
              dchi22 = dchi22*enhanc
              det    = dchi11*dchi22-dchi12*dchi12
              det    = sign(max(abs(det),detmin),det) ! lower limit for 'det'
              dinv11 = dchi22/det
              dinv22 = dchi11/det
              dinv12 = -dchi12/det
              dxlam  = -(dinv11*dchi1+dinv12*dchi2)
              dxlam2 = -(dinv12*dchi1+dinv22*dchi2)
              if(iftest.ge.0) 
     &           write(6,'(a,4g12.4)')
     &           ' chi: dxlam,dxlam2=',dxlam,dxlam2
            endif
c
            if(iftest.ge.1) 
     &           write(6,'(a,4g12.4)')
     &           ' dxlam,dxlam2=',dxlam,dxlam2
            quenf1 = dxllim/max(abs(dampin*dxlam ),dxllim)
            quenf2 = dx2lim/max(abs(dampin*dxlam2),dx2lim)
            quenf  = min(quenf1,quenf2)
            dxlam  = quenf*dampin*dxlam
            dxlam2 = quenf*dampin*dxlam2
cold            dxlam  = sign(min(abs(dampin*dxlam ),dxllim),dxlam )
cold            dxlam2 = sign(min(abs(dampin*dxlam2),dx2lim),dxlam2)
            if(iftest.ge.1) 
     &           write(6,'(a,4g12.4)')
     &           ' dxlam,dxlam2=',dxlam,dxlam2
            eferm  = eferm + dxlam
            xlam2  = xlam2 + dxlam2
c     
         enddo
c     
         if(iftest.ge.0) write(6,'(a)')
     &      ' no convergence of pairing parameters in PAIRDN, LN-branch'
c
 1099    continue
         parnum = f1+npart
         if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     &        ' PAIRDN-LN: iter,xlam,xlam2,f1,f2=',it,eferm,xlam2,f1,f2
      endif

!     print result
!      WRITE(6,*) 'pairing weights:'
!      WRITE(6,'(5(1pg13.5))') gw(1:nmax)

c
c     finally compute average gap and pairing energy
c
      sumuv  = 0.0
      sumduv = 0.0
      sumdv2 = 0.0
      sumv2  = 0.0
      sumepa = 0.0
      do na=1,nmax
         IF(tpairspace) THEN
           deltaa = delta(na)
         ELSE
           deltaa = delta(na)*pairgw(na)
         END IF
         edif   = e(na)-eferm
         equasi = sqrt(edif*edif+deltaa*deltaa)
         v2     = MAX(half-half*edif/equasi,v2lowlim)
         phact  = half*ph(na)   ! for axial
         IF(tpairspace) THEN
           vol    = phact*sqrt(max(v2-v2*v2,xsmall))*pairgw(na)
         ELSE
           vol    = phact*sqrt(max(v2-v2*v2,xsmall))
         END IF
         IF(na==nblock .OR. na==ntcblock) vol=0D0
         sumuv  = vol        + sumuv
         sumduv = vol*deltaa + sumduv
         sumv2  = phact*v2 + sumv2
         sumdv2 = phact*v2*deltaa + sumdv2
         sumepa = phact*half*(deltaa*deltaa)/equasi + sumepa
      enddo
      sumuv  = max(sumuv,xsmall)
      avdelt = sumduv/sumuv
      avdeltv2 = sumdv2/sumv2
      epair  = sumduv           ! or - ?
      avg    = epair/sumuv**2
!      if(iab.ge.-1)
      if(iab.ge.-5)
     &     write(6,'(a,i4,a,1(1pg13.5),a,2(1pg13.5))')
     &     ' PAIRDN: iteration=',it,' part.nr.=',parnum-partnm,
     &     ' eferm,epair,avdelta=',eferm,epair,avdelt,avdeltv2
      if(iab.ge.0) write(6,'(5x,a,3g12.4)')
     &     ' av.Delta,av.G,pair.E=',avdelt,avg,epair
c     
      return
      end
c-----cutoff------------ variant for the axial SHF package ---------------------
c
      subroutine cutoff(e,gw,ph,anz,nmax,eferm,efadd,delta,
     &                  partnm,iter,eps,iab,iq)
c
c     this routine determines the effective "Fermi-energy" and
c     cut-off distribution for the pairing space.
c
c     the include file 'params.inc' communicates the dimensioning of
c     the fields 'e' and 'gw' with the parameter 'ksl2'. it allows
c     furthermore to make the whole subroutine double precision by
c     transfering a line 'implicit double precision (a-h,o-z)'.
c
      include 'paramr.inc'
      include 'levels.inc'
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      gw      = array(1..nmax) of cut-off weights,
c                input is used as initial guess for the iteration.
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels, to be provided by the calling routine.
c      anz     = desired effective particle number
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ksl2'
c                 and smaller than 46, for dimensioning reasons).
c      eferm   = true fermi-energy as given by the pairing routine.
c      efadd   = additional "fermi-energy" for the pairing-space weights
c      delta   = relative width of the distribution
c      iter    = number of pairing iterations
c      eps     = end condition for iteration;
c                iteration terminates if error in effective particle number
c                is smaller than 'eps' and if 'eferm' and/or 'delta'
c                are smaller than 'eps'.
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c
c     output via list:
c      gw      = array(1..nmax) of the occupation probabilities
c                finally achieved.
c      efadd   = additional "fermi-energy" for the pairing-space weights
c      partnm  = particle number, from adding 'gw' with degeneracy 'ph'.
c
c     some auxiliary fields are:
c      dparn   = array(1..3) of variation in the particle number as
c                consequence of 'deferm' and 'ddelta'.
c      dgap    = array(1..3) of variation in the gap equation as
c                consequence of 'deferm' and 'ddelta'.
c      gr      = double array of variation of particle number and of
c                gap equation, used for construction of the tangantial
c                step
c
      dimension e(npsid),gw(npsid),ph(npsid)
      dimension dparn(2),gr(1,1)
c
cold      data ph/1.,2.,1.,3.,2.,1.,4.,3.,2.,1.,5.,4.,3.,2.,1.,6.,5.,
cold     &        4.,3.,2.,1.,7.,6.,5.,4.,3.,2.,1.,8.,7.,6.,5.,4.,3.,2.,1.,
cold     &        9.,8.,7.,6.,5.,4.,3.,2.,1.,10./
      data gr,dparn/3*0./
c
c     further numerical variables are:
c      deferm  = step in the fermi-energy, to evaluate
c                numerically the derivatives.
c      ddelta  = step in the gap, to evaluate numerically the
c                derivatives; used only in the "constant force" approach.
c      dampin  = damping factor for the tangential step
c      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
c
c     note: these parameters are adjusted for an application with
c           the nuclear shell model. it is assumed that the single
c           particle energies are given in units of [mev] and display
c           a reasonable nuclear spectrum.
c
      data dampin,deferm/0.9,1.0/
      data efmax,detmin/3.,.1e-3/
c      data defmin/.1e-3/
c
c      data xmaxlw,x10lw/1.0e-30,1.0e-10/
c      data one,zero/1.0,0.0/
      data fifty/50.0/
      data iftang/0/
c
c-------------------------------------------------------------------------------
c
c      set switch variables
c
      iwarng = 0
c
      if(iab.ge.1)
     & write(6,'(/a,2i4,a,4f9.4)')
     &   ' *** cutoff entered with: iter,nmax=',iter,nmax,
     &   '  delta,efadd,eferm,partnm=',delta,efadd,eferm,anz
      if(iab.ge.2)
     & write(6,'(a,10(/1x,10f8.3))')
     &   ' single particle energies:',(e(i),i=1,nmax)
c
c      prepare fermi-energy 'eferm', if not yet initialised.
c
      if(efadd.eq.zero) then
       isum   = 0
       do it=1,nmax
        itsum  = it
        isum   = int(ph(it)+.000001) +isum
        if(one*isum-anz.ge.0) goto 110                        ! axial
       enddo
       stop ' not enough states too reach particle number'
  110  continue
       efadd  = e(itsum)+deferm-eferm
       if(iab.ge.1) write(06,'(5x,a,i5,a,g11.4)')
     &   'initialistion: isum=',isum,'  efadd=',efadd
c
      endif
c
c
c-------------------------------------------------------------------------------
c
c     iterate Fermi energy and the distribution
c
c-------------------------------------------------------------------------------
c
c     tangent step
c
      if(iftang.eq.1) then
        do it=1,iter
c
c       particle number at 'eferm' (k=1) and 'eferm' varied (k=2).
c
         do ka=1,2
c                                        set actual fermi-energy 'elam'
          k      = 3-ka
          deltac  = efadd*delta
          if(k.eq.2) then
           elam   = eferm+efadd+deferm*deltac
          else
           elam   = eferm+efadd
          endif
c                                         compute actual particle number
          parnum = zero
          do i=1,nmax
codd          if(isTSym.eq.1) then
codd             ip = i
codd          else
codd          ip=itpartn(i+(iq-1)*nst(1))-(iq-1)*nst(1)
codd          endif
           arg    = (e(i)-elam)/deltac
codd          arg = (0.5D0*(e(i)+e(ip))-elam)/deltac

           if(arg.gt.fifty) then
             gw(i)   = zero
           elseif(arg.lt.-fifty) then
             gw(i)   = one
           else
             gw(i)  = 1.0/(1.0+exp(arg))
           endif
           parnum = gw(i)*ph(i)    + parnum
          enddo
c                                         store particle number
          dparn(k) = parnum-anz                                 ! axial
         enddo
c
c      gradient step to next 'eferm'   (including step limiters)
c
         gr(1,1) = (dparn(2)-dparn(1))/(deferm*deltac)
         det    = sign(max(abs(gr(1,1)),detmin),gr(1,1))
c         efold  = efadd
         efstep = dampin*dparn(1)/det
         efadd  = efadd-sign(min(abs(efstep),efmax),efstep)
c         deferm = efold-efadd
c         deferm = sign(max(abs(deferm),defmin),deferm)
         if(iab.eq.0)
     &    write(6,'(a,i4,2(a,g11.3))')
     &     ' pairing-weight iteration it=',it,'  :  efadd=',efadd,
     &     '  dparn(1)=',dparn(1)
         if(iab.ge.1)
     &    write(6,'(a,i4,2(a,g11.3)/a,2g11.3,a,g11.3/
     &              a,10f7.4,10(/6x,10f7.4))')
     &     ' iteration at it=',it,' :   efadd=',efadd,
     &     '  det=',det,
     &     ' dparn=',dparn(1),dparn(2),' ;   gr=',gr(1,1),
     &     '   gw=',(gw(n),n=1,nmax)
         if(abs(dparn(1)).lt.eps)   goto 350
        enddo
        iwarng = 1
  350   continue
      else
c
c     secant step
c
       eferm1 = 1.0D10
       eferm2 = -1.0D10
       do i=1,nmax
         if(eferm2.lt.e(i)) eferm2 = e(i)
         if(eferm1.gt.e(i)) eferm1 = e(i)
       enddo
c       eferm2 = eferm2+2.0*max((eferm2-eferm)*delta,delta)
       do it=1,iter
         if(it.eq.1) then
           eferm3 = eferm1
         elseif(it.eq.2) then
           eferm3 = eferm2
         else
           eferm3 = eferm1 - dparn1*(eferm2-eferm1)/(dparn2-dparn1)
         endif
         efadd   = eferm3-eferm
         deltac  = max(efadd*delta,delta)
         parnum  = 0.0
         do i=1,nmax
           if(isTSym.eq.0)then
            ip = i
           else
            ip = itpartn(i+(iq-1)*nst(1))-(iq-1)*nst(1)
           endif
codd          arg    = (e(i)-eferm3)/deltac
          arg    = (0.5D0*(e(i)+e(ip))-eferm3)/deltac
          if(arg.gt.fifty) then
            gw(i)   = zero
          elseif(arg.lt.-fifty) then
            gw(i)   = one
          else
            gw(i)  = 1.0/(1.0+exp(arg))
          endif
          parnum = gw(i)*ph(i)    + parnum
         enddo
         dparnm = parnum-anz                                      ! axial
         if(it.eq.1) then
           dparn1  = dparnm
           if(dparn1.gt.zero) then
             write(6,'(a/i4,5(1pg12.4))')
     &        ' it,deltac,eferm1,eferm2,dparn1,dparn2=',
     &          it,deltac,eferm1,eferm2,dparn1,dparn2
             write(6,'(a,10(/1x,10f8.3))')
     &        ' single particle energies:',(e(i),i=1,nmax)
             stop ' lower state too high in CUTOFF'
           endif
         elseif(it.eq.2) then
           dparn2  = dparnm
           if(dparn2.lt.zero) then
             write(6,'(a/i4,5(1pg12.4))')
     &        ' it,deltac,eferm1,eferm2,dparn1,dparn2=',
     &          it,deltac,eferm1,eferm2,dparn1,dparn2
             write(6,'(a,10(/1x,10f8.3))')
     &        ' single particle energies:',(e(i),i=1,nmax)
             stop ' upper state too low in CUTOFF'
           endif
         else
           dparn3  = dparnm
         endif
         if(iab.ge.2)
     &    write(6,'(a,i4,2(a,g11.3))')
     &    ' secant iteration it=',it,'  :  eferm=',eferm3,
     &    '  dparn3=',dparn3
         if(it.gt.2) then
           if(abs(dparn3).lt.eps) then
             goto 399
           elseif(dparn3.gt.zero) then
             eferm2 = eferm3
             dparn2 = dparn3
           else
             eferm1 = eferm3
             dparn1 = dparn3
           endif
         endif
       enddo
  399  continue
      endif
c
c-------------------------------------------------------------------------------
c
c     finally the particle-number is checked and printed
c
c-------------------------------------------------------------------------------
c
      partnm = zero
      do i=1,nmax
       partnm = gw(i)*ph(i) + partnm
      enddo
c
      if(iab.ge.-1)
     &  write(6,'(a,i4,2(1pg13.5))')
     &   ' CUTOFF finished with iter,part.nr.,efadd=',it,partnm,efadd
      if(iab.ge.2)
     & write(6,'(a,10(/1x,10f8.3))')
     &   ' pairing-weights:',(gw(i),i=1,nmax)
c
c     in case of incomplete convergence print a warning
c
      if(iab.ge.-1 .and. iwarng.eq.1)
     & write(6,'(2(a,g11.3),2(a,i3))')
     &  ' ---> no convergence in cutoff! err(part.nr)=',dparn(1)
c
      return
      end
c-----lm2def-------------------------------------------------------------------
c
      subroutine lm2def(spen,delta,ph,pairgw,nspd,nspg,iwfoff,npart,
     &                  xlam,xlam2,vv2in,vv2out,iftest,  f1,f2)
c
c---------------------------------------------------------------------
c
      include 'paramr.inc'
      parameter (kmaxln=npsid)
c
      include 'resden.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
c
c---------------------------------------------------------------------
c
c     routine to compute the second order multiplies 'lambda_2`
c     for the LN scheme in connection with the delta-force pairing.
c
c     !!!! attention: version for v**2 pairing needs yet checking
c     !!!! attention: not yet adapted to blocking
c
c     input:
c       spen   = array of s.p. energies
c       delta  = array of pairing gaps
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       pairgw = cutoff weight complementing gap 'delta'
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       iwfoff = offset for addressing wavefunctions,
c                note that wf's require appropriate
c                absolute addresses because information
c                on wavefunctions is entered via common /wavefr/ .
c       npart  = particle number
c       vpair  = strength of pairing force
c       xlam   = Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2in  = array of initial internal occupation probabilities
c       iftest = switch for testprints
c       ilipng = switch between BCS and Lipkin-Nogami,
c                 1 --> f3=0
c     the following two parameters are entered via common /resden/ or /dynpar/
c     whereby proton <--> neutron is switched by 'iwfoff'
c       erspln = Skyrme response energy for 'lambda_2'
c       suv2ra =   u^2*v^2 from response routine, for normalization
c       vpair  = pairing strength
c     output:
c       f1     = sum(2v**2)-n         =  particle number condition
c       f2     = xlam2-vpair/4*(sum(u**3*v)sum(v**3*u)-
c                               sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                               sum(u**4*v**4))
c              = LN condition
c       vv2out = new occupation weights
c
c---------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2in(kmaxln),vv2out(kmaxln),ph(kmaxln),
     &          delta(kmaxln),pairgw(kmaxln)      
c
c       iterv2 = max.nr. iterations in 'vv2'
c       epsv2  = precision on 'vv2'
c       epsdv2 = limit on change of energy in determining 'vv2'
c       dampin = general reduction of the step, top stabilize
      data iterv2,epsv2,epsdv2,dampin
     &    /    50, 2*1.0D-7   ,  0.5 /
c
      data half/0.5/
c
c---------------------------------------------------------------------
c
c     preparation
c
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      fnn    = npart
      if(iwfoff.eq.0) then
        erspln = rspln(1)
        suv2ra = suv2rs(1)
        vpair  = pairfc(1)
      else
        erspln = rspln(2)
        suv2ra = suv2rs(2)
        vpair  = pairfc(2)
      endif
      if(iftest.ge.3) 
     &  write(6,'(a,i5,3g12.4)')
     &   ' iwfoff,vpair,erspln,suv2ra=',iwfoff,vpair,erspln,suv2ra
c
c     reset accumulators for modified pairing currents
c
      do ir=0,nr
      do iz=0,nz
        wrk1(iz,ir)  = 0.0
        wrk2(iz,ir)  = 0.0
      enddo
      enddo      
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = zero                                  ! sum v^2
      suv   = zero                                  ! sum u*v
      suv4  = zero                                  ! sum u^4 v^4
      suv2  = zero                                  ! sum u^2 v^2
      su3v  = zero                                  ! sum u*v u^2
      suv3  = zero                                  ! sum u*v v^2
c
      alpha  = 4.*xlam2
c
      do k=nspd,nspg       
c
c     iteration loop to compute occupation 'vv2'
c
        IF(tpairspace) THEN
          d2     = delta(k)**2
        ELSE
          d2     = (delta(k)*pairgw(k))**2
        END IF
        alphd2 = alpha*half*d2
        e     = spen(k)-xlam
        v2    = vv2in(k)                                ! start from previous value
        IF(tpairspace) THEN
          x     = e+alpha*v2*pairgw(k)
        ELSE
          x     = e+alpha*v2
        END IF
        if(iftest.ge.3)
     &   write(7,'(a,i5,a,3(1pg12.4))')
     &    ' state nr.',k,' :  v2in,e,delta=',v2,e,delta(k)*pairgw(k)
c
c                            iterate to s.c. occup.numbers and lambda_2
c
        do it=1,iterv2
          eq    = sqrt(d2+x*x)
          v2    = half-half*x/eq
          IF(tpairspace) THEN
            f     = e+alpha*v2*pairgw(k)-x
          ELSE
            f     = e+alpha*v2-x
          END IF
          derf  = -alphd2/(eq*eq*eq)-one
          dx    = dampin*f/derf         ! optionally switch on if unstable
CTEST          dx    = f/derf
          if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
          x     = x-dx
        enddo
c
ctst        write(6,'(a)')
ctst     &    ' maximal number of iterations exceeded in "LM2DEF"'
        if(iftest.ge.3) write(7,'(a)')
     &    ' maximal number of iterations exceeded in "LM2DEF"'
c
    7   continue
        v2   = max(sqrt(smal),v2)
        if(iftest.ge.1)
     &   write(7,'(a,i5,a,i5,g12.3)')
     &    ' state nr.',k,' :  iterv,v2',it,v2
        vv2out(k)=v2
c
c       accumulate sums with new 'v2'
c
        IF(tpairspace) THEN
          xmult = ph(k)*pairgw(k)
        ELSE
          xmult = ph(k)
        END IF
        u2    = one-v2
        u2    = max(sqrt(smal),u2)
        uv    = sqrt(u2*v2)
        sv2   = sv2  + xmult*v2
        suv   = suv  + xmult*uv
        suv4  = suv4 + xmult*u2*u2*v2*v2
        suv2  = suv2 + xmult*u2*v2
        su3v  = su3v + xmult*uv*u2
        suv3  = suv3 + xmult*uv*v2
c
c       accumulate modified pair currents
c
        istate = k+iwfoff
        call sdenti(psi(0,0,1,istate),wrk1,uv*u2,istate)
        call sdenti(psi(0,0,1,istate),wrk2,uv*v2,istate)
c
      enddo
c
c     first result: particle number 
c
      f1     = sv2-fnn          
c
c     second result: the LN equation
c     
c     but the numerator from Skyrme response is a bit special.
c     it was computed with an older "v2". The corresponding
c     "sum V2*U2" is stored in "suv2ra". In order to avoid
c     all to frequent recomputation of that expensive piece,
c     we rescale it here just by the new "suv2".
c
      denom  = 32.0*(suv2*suv2-suv4)                   !  denominator
      do ir=0,nr
      do iz=0,nz
        wrk1(iz,ir)  = wrk1(iz,ir)*wrk2(iz,ir)
      enddo
      enddo      
      enume1 = -vpair*2.0*f1int(wrk1)                  !  1. numerator
      enume2 = suv2*suv2*erspln/(suv2ra*suv2ra)        !  2. numerator
      f2     = xlam2+(enume1+enume2)/denom             ! p.n. fluctuation eq.
c
      if(iftest.ge.3) 
     &  write(6,'(a,6g12.4)')
     &   ' suv2,suv4,suv2ra,erspln,vpair,echi=',
     &   suv2,suv4,suv2ra,erspln,vpair,f1int(wrk1)
c
      return
      end

c========================================================================
c========================================================================
c    
c     Here follows the Pairing and Lipkin-Nogami package of
c      -->   G.A.Leander and W.Nazarewicz, march 1985   <--
c     revised by P.-G. Reinhard, june 1994
c
c     BCS is standard, but includes correction 'gpair*v^2' to
c     the s.p. energies.
c
c     The method of Lipkin and Nogami and coworkers is explained,
c     e.g. in Nucl. Phys. A201 (1973) 357 .
c
c     This version is adapted for use in the axial Skyrme code "SKYAX".
c
c========================================================================
c-----bcslip-------------------------------------------------------------------
c
      subroutine bcslip(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                  vv2,vln,epair,ecorr,entot,erspln,suv2rs,
     &                  ilipng,iftest)
c
ctest      implicit double precision (a-h,o-z)
c
ctest      parameter (kmaxln=46)
      include 'paramr.inc'                       ! activate that for
      parameter (kmaxln=npsid)                    ! SKYAX package
c
c------------------------------------------------------------------------------
c
c     Program to solve the Lipkin-Nogami, or BCS, equations
c       gap equation:      0 = f1 = delta-gpair*sum(u*v)
c       particle number:   0 = f2 = sum(2v**2)-n
c       p.n. fluctuation:  0 = f3 = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                                   sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                                   sum(u**4*v**4))
c     the f3 is omitted for BCS.
c
c     The 3 equations are solved by Newton's method.
c     The matrix of derivatives is computed by finite differences.
c     The internal numerical parameters are explained below.
c
c------------------------------------------------------------------------------
c
c     The input parameters are:
c
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       npart  = particle number
c       gpair  = strength of pairing force
c       delta  = initial value for pairing gap
c       xlam   = initial value for Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2    = array of internal occupation probabilities (initial)
c       erspln = Skyrme response energy for 'lambda_2', needed for
c                Skyrme H.F. codes, to be provided in a separte routine.
c                In all other cases set to 0.0 .
c       suv2rs = u^2*v^2 from Skyrme response routine, for normalization.
c                In all other cases set to 1.0 .
c       ilipng = switch between BCS and Lipkin-Nogami
c                 1 --> Lipkin-Nogami
c                 0 --> BCS
c       iftest = switch for testprints
c
c     The output parameters are:
c       delta  = result for pairing gap
c       xlam   = result for Fermi energy
c       xlam2  = result for fluctuation constraint
c       vv2    = array of internal occupation probabilities
c       vln    = array of occupation probabilities for external observables
c       epair  = naive pairing energy , -delta^2/gpair
c       ecorr  = consistent pairing energy as correction for the hf energy:
c                =-delta**2/gpair
c                 -gpair*sum{v(i)**4-v(i)**2}
c                 -4*xlam2*sum_(u(i)*v(i))**2
c                = -(4*xlam2-gpair)*sum_(u(i)*v(i))**2
c       entot  = total s.q.p. energy (= s.p. sum - g * npart/2 for normal state)
c
c------------------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2(kmaxln),vln(kmaxln),ph(kmaxln)
c
ctest      dimension derv(3,3)
c
c------------------------------------------------------------------------------
c
c     some useful constants, allow easy switch single <--> double precision
c
cold      data zero,one/0.0,1.0/
c
c------------------------------------------------------------------------------
c
c     iteration parameters:
c      itbcsl  = max.nr. iterations in 3 pairing parameters
c      iterv2  = max.nr. iterations in 'vv2'
c      epsf1   = precision on pairing equation 'f1'
c      epsf2   = precision on pairing equation 'f2'
c      epsf3   = precision on pairing equation 'f3'
c      epsv2   = precision on 'vv2'
c      epsdv2  = limit on change of energy in determining 'vv2'
c      detmin  = limiting minimum value for the determinant in linear eq.
c      ddllim  = limit on stepsize in 'delta'
c      dxllim  = limit on stepsize in 'xlam'
c      dx2lim  = limit on stepsize in 'xlam2'
c      dampin  = general reduction of the step, top stabilize
c      delstp  = step in 'delta' for first derivative in finite difference
c      xlstp   = step in 'xlam' for first derivative in finite difference
c      xl2stp  = step in 'xlam2' for first derivative in finite difference
c
      data itbcsl,iterv2,epsf1,epsf2,epsf3,epsv2,epsdv2
     &    /    50,    50,     5*1.0D-8                 /
      data detmin,ddllim,dxllim,dx2lim,dampin
     &    /1.0D-1,  0.2,    3.0,   1.0, 0.8/
      data delstp,xlstp,xl2stp
     &    /  0.01, 0.01,  0.01/
c
c---------------------------------------------------------------------
c
c     preparation and tests
c
      if(nspg.gt.kmaxln) stop ' last level > "KMAXLN" in BCSLIP'
      if(ilipng.ne.1) then
       xlam2 = zero
       dxlam2 = zero
      endif
c
c    iteration of the three pairing parameters 'delta', 'xlam', and 'xlam2' .
c
      do 1000 itbcs=1,itbcsl
      if(iftest.gt.0)
     &   write(6,'(a,i5,a)')  ' iteration nr.',itbcs,' :'
c
c     initial call to pairing equations
c
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &           vv2,vv2,erspln,suv2rs,iftest,ilipng,
     &           iterv2, epsv2,epsdv2,dampin,
     &           f1,f2,f3)
      if(iftest.gt.0) then
        write(6,'(a,3g12.3/a,3g12.3)')
     &   ' delta,xlam,xlam2=',delta,xlam,xlam2,
     &   '         f1,f2,f3=',f1,f2,f3
      endif
c
c     test termination
c
      if(itbcs.gt.1 .and. abs(f1).lt.epsf1 .and. abs(f2).lt.epsf2
     &   .and. abs(f3).lt.epsf3) goto 1099
c
c     derivatives by finite differences
c
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta*(one+delstp),
     &            xlam,xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,  ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1del,f2del,f3del)
      call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam+xlstp,
     &            xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,       ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1xl,f2xl,f3xl)
      if(ilipng.eq.1)
     & call eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2+xl2stp,
     &            vv2,vln,erspln,suv2rs,iftest,ilipng,             ! vln = aux
     &            iterv2, epsv2,epsdv2,dampin,
     &            f1xl2,f2xl2,f3xl2)
      ddelac = delta*delstp
      ddf1   = (f1del-f1)/ddelac
      ddf2   = (f2del-f2)/ddelac
      dlf1   = (f1xl-f1)/xlstp
      dlf2   = (f2xl-f2)/xlstp
      if(ilipng.eq.1) then
        dl2f1  = (f1xl2-f1)/xl2stp
        dl2f2  = (f2xl2-f2)/xl2stp
        ddf3   = (f3del-f3)/ddelac
        dlf3   = (f3xl-f3)/xlstp
        dl2f3  = (f3xl2-f3)/xl2stp
c
c     solve the linear equations for the case of 3 steps
c
c   define matrix in orders
c      derv(1,1) = ddf1
c      derv(2,1) = ddf2
c      derv(3,1) = ddf3
c      derv(1,2) = dlf1
c      derv(2,2) = dlf2
c      derv(3,2) = dlf3
c      derv(1,3) = dl2f1
c      derv(2,3) = dl2f2
c      derv(3,3) = dl2f3
c
c    matrix of inverse is computed explicitely and stored on bij
c
        det    =(ddf1*dlf2*dl2f3-ddf1*dl2f2*dlf3-ddf2*dlf1*dl2f3
     #          +ddf2*dl2f1*dlf3+ddf3*dlf1*dl2f2-ddf3*dl2f1*dlf2)
        det    = sign(max(abs(det),detmin),det)        ! lower limit for 'det'
        detinv = one/det
        b11 = (dlf2*dl2f3-dl2f2*dlf3)*detinv
        b12 = -(dlf1*dl2f3-dl2f1*dlf3)*detinv
        b13 = (dlf1*dl2f2-dl2f1*dlf2)*detinv
        b21 = -(ddf2*dl2f3-dl2f2*ddf3)*detinv
        b22 = (ddf1*dl2f3-dl2f1*ddf3)*detinv
        b23 = -(ddf1*dl2f2-dl2f1*ddf2)*detinv
        b31 = (ddf2*dlf3-dlf2*ddf3)*detinv
        b32 = -(ddf1*dlf3-dlf1*ddf3)*detinv
        b33 = (ddf1*dlf2-dlf1*ddf2)*detinv
ctest      if(iftest.ge.1) write(6,'(a,3g12.4,2(/19x,3g12.4))')
ctest     &   ' derivative matrix:',derv
ctest      call  MATIN2 (derv,3,3,NERROR,DET)
ctest      if(iftest.ge.1)
ctest     &  write(6,'(a,9g12.4)')
ctest     &   'd:',derv,
ctest     &   'b:',b11,b21,b31,b12,b22,b32,b13,b23,b33
ctest      ddelta = -derv(1,1)*f1-derv(1,2)*f2-derv(1,3)*f3
ctest      dxlam  = -derv(2,1)*f1-derv(2,2)*f2-derv(2,3)*f3
ctest      dxlam2 = -derv(3,1)*f1-derv(3,2)*f2-derv(3,3)*f3
        ddelta = -b11*f1-b12*f2-b13*f3
        dxlam  = -b21*f1-b22*f2-b23*f3
        dxlam2 = -b31*f1-b32*f2-b33*f3
        if(iftest.gt.0)
     &   write(6,'(a,4g12.3)')
     &    ' det,ddelta,dxlam,dxlam2=',det,ddelta,dxlam,dxlam2
      else
c
c     case of BCS, step in two quantities
c
        det    = -(dlf1*ddf2-dlf2*ddf1)
        det    = sign(max(abs(det),detmin),det)      ! lower limit for 'det'
        ddelta = (f2*dlf1-f1*dlf2)/det
        dxlam  = (f1*ddf2-f2*ddf1)/det
      endif
c
c     limit stepsize if necessary
c
      ddelta = dampin*sign(min(abs(ddelta/delta),ddllim),ddelta)
      dxlam  = dampin*sign(min(abs(dxlam ),dxllim),dxlam )
      if(ilipng.eq.1)
     &   dxlam2 = dampin*sign(min(abs(dxlam2),dx2lim),dxlam2)
c
c     step forward
c
      delta  = delta*exp(ddelta)
      xlam   = xlam  + dxlam
      if(ilipng.eq.1)
     &    xlam2  = xlam2 + dxlam2
c
 1000 continue
c
      if(iftest.ge.0) write(6,'(a)')
     & ' no convergence of pairing parameters in BCSLIP'
c
 1099 continue
      if(iftest.ge.0) write(6,'(a,i4,6(1pg11.3))')
     & ' BCSLIP: iter,delta,xlam,xlam2,f1,f2,f3=',
     &   itbcs,delta,xlam,xlam2,f1,f2,f3
c
c     compute some observables after convergence
c
      suv    = zero
      sv2    = zero
      sv4    = zero
      suv2   = zero
      sspe   = zero
      do 200 k=nspd,nspg
      v2     = vv2(k)
      uv     = sqrt(v2-v2*v2)
      vln(k) = uv                                     ! auxiliary storage
      xmult  = ph(k)
      suv    = suv    + xmult*uv
      sv2    = sv2    + xmult*v2
      sv4    = sv4    + xmult*v2*v2
      suv2   = suv2   + xmult*(v2 - v2*v2)
      sspe   = sspe   + xmult*spen(k)*v2
  200 continue
cold      ecorr  = -delta*delta/gpair-gpair*sv4-4.0*xlam2*suv2+gpair*sv2
      ecorr  = -delta*delta/gpair-(4.0*xlam2-gpair)*suv2
      entot  =  sspe - gpair * npart/2
      epair1 = -gpair*sv2
      efluc  = -4.*xlam2*suv2
      epair  = -gpair*suv
cold      etot   = 2.0*sspe+epair1+efluc+epair
cold      epairr = epair
      epair1 = epair1+efluc+epair
c
c           corrected coefficients, vln(i)=<a^+_i a_i>_ln
c
      if(ilipng.eq.1) then
        sum2=0.
        sum4=0.
        sumc=0.
        sumc2=0.
c
        do 120 k=nspd,nspg
        v2    = vv2(k)
        xmult = ph(k)
cold        uv    = sqrt(v2-v2*v2)
        uv    = vln(k)
        s     = uv+uv
        c     = v2+v2-one
        sum2  = sum2  + xmult*s**2
        sum4  = sum4  + xmult*s**4
        sumc  = sumc  + xmult*s*s*c
        sumc2 = sumc2 + xmult*(s*c)**2
  120   continue
c
        od=sum2*0.25/(sumc2*sum2-sumc*sumc+0.25*sum2*(sum2*sum2-sum4))
c
        do 130 k=nspd,nspg
        v2    = vv2(k)
        uv    = vln(k)
cold        uv    = sqrt(v2-v2*v2)
        s     = uv+uv
        c     = v2+v2-one
        vln(k)= v2+od*s*s*(c*sum2-sumc)
  130   continue
      else
        do 135 k=nspd,nspg
        vln(k)= vv2(k)
  135   continue
      endif
c
      return
      end
c-----eqdef-------------------------------------------------------------------
c
      subroutine eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,
     &                 vv2in,vv2out,erspln,suv2rs,iftest,ilipng,
     &                 iterv2, epsv2,epsdv2,dampin,
     &                 f1,f2,f3)
c
c---------------------------------------------------------------------
c
ctest      implicit double precision (a-h,o-z)
c
ctest      parameter (kmaxln=46)
      include 'paramr.inc'
      parameter (kmaxln=npsid)
c
c---------------------------------------------------------------------
c
c     routine to compute the three pairing equations, as explained
c     in the header of 'BCSLIP'.
c
c     input:
c       spen   = array of s.p. energies
c       ph     = array of multiplicities, e.g. a j=7/2 state has ph = 4 .
c       nspd   = lowest level in pairing
c       nspg   = highest level in pairing
c       npart  = particle number
c       gpair  = strength of pairing force
c       delta  = pairing gap
c       xlam   = initial value for Fermi energy
c       xlam2  = initial value for fluctuation constraint
c       vv2in  = array of initial internal occupation probabilities
c       erspln - Skyrme response energy for 'lambda_2'
c       suv2rs -   u^2*v^2 from response routine, for normalization
c       iftest = switch for testprints
c       ilipng = switch between BCS and Lipkin-Nogami,
c                 1 --> f3=0
c       iterv2 = max.nr. iterations in 'vv2'
c       epsv2  = precision on 'vv2'
c       epsdv2 = limit on change of energy in determining 'vv2'
c       dampin = general reduction of the step, top stabilize
c     output:
c       f1     = delta-gpair*sum(u*v) =  gap equation
c       f2     = sum(2v**2)-n         =  particle number condition
c       f3     = xlam2-gpair/4*(sum(u**3*v)sum(v**3*u)-
c                               sum(u**4*v**4))/((sum(u**2*v**2))**2 -
c                               sum(u**4*v**4))
c              = p.n. fluctuation
c       vv2out = new occupation weights
c
c---------------------------------------------------------------------
c
c
      dimension spen(kmaxln),vv2in(kmaxln),vv2out(kmaxln),ph(kmaxln)
c
cold      data zero,half,one/0.0,0.5,1.0/
      data half/0.5/
c
c---------------------------------------------------------------------
c
c     preparation
c
      fnn    = npart
      if(nspg.gt.kmaxln) stop ' last level > "kmaxln"'
      d2     = delta*delta
      alpha  = 4.*xlam2-gpair
      alphd2 = alpha*half*d2
c
c     reset accumulators for combinations of 'u' and 'v'.
c
      sv2   = zero                                  ! sum v^2
      suv   = zero                                  ! sum u*v
c      if(ilipng.eq.1) then
       suv4  = zero                                  ! sum u^4 v^4
       suv2  = zero                                  ! sum u^2 v^2
       su3v  = zero                                  ! sum u*v u^2
       suv3  = zero                                  ! sum u*v v^2
c      endif
c
      do 500 k=nspd,nspg                            ! each state separately
c
c     iteration loop to compute occupation 'vv2'
c
      e     = spen(k)-xlam
      v2    = vv2in(k)                                ! start from previous value
      x     = e+alpha*v2
c
      do 6 it=1,iterv2
      eq    = sqrt(d2+x*x)
      v2    = half-half*x/eq
      f     = e+alpha*v2-x
      derf  = -alphd2/(eq*eq*eq)-one
ctest      dx    = dampin*f/derf         ! optionally switch on if unstable
      dx    = f/derf
      if(abs(f).lt.epsv2 .and. abs(dx).lt.epsdv2) go to 7
      x     = x-dx
    6 continue
c
      write(6,'(a)')
     &  ' maximal number of iterations exceeded in "EQDEF"'
c
    7 continue
      if(iftest.gt.1)
     & write(6,'(a,i5,a,i5,g12.3)')
     &  ' state nr.',k,' :  iterv,v2',it,v2
      vv2out(k)=v2
c
c     accumulate sums with new 'v2'
c
      xmult = ph(k)
      u2    = one-v2
      uv    = sqrt(u2*v2)
      sv2   = sv2  + xmult*v2
      suv   = suv  + xmult*uv
      if(ilipng.eq.1) then
       suv4  = suv4 + xmult*u2*u2*v2*v2
       suv2  = suv2 + xmult*u2*v2
       su3v  = su3v + xmult*uv*u2
       suv3  = suv3 + xmult*uv*v2
      endif
c
  500 continue
c
c     final assembling of the three pairing equations
c
      f1     = delta-gpair*suv                    ! gap equation
      f2     = sv2+sv2-fnn                        ! particle number
c     
c     the numerator from Skyrme response is a bit special.
c     it was computed with an older "v2". The corresponding
c     "sum V2*U2" is stored in "suv2rs". In order to avoid
c     all to frequent recomputation of that expensive piece,
c     we rescale it here just by the new "suv2".
c
      if(ilipng.eq.1) then
       denom  = (suv2*suv2-suv4)                   !  denominator/32
       enume1 = -gpair/4.0*(su3v*suv3-suv4)        !  1. numerator/32
       enume2 = suv2*suv2*erspln/(suv2rs*suv2rs*32.0)   !  2. numerator/32
       f3     = xlam2+(enume1+enume2)/denom        ! p.n. fluctuations
cold   f3=xlam2-gpair/4.*(su3v*suv3-suv4)/(suv2*suv2-suv4)
       endif
c
      return
      end
c================================================================
      function gavrg(ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0)
c
c         gavrg calculates the strength of the monopole pairing
c         interaction using the average gap method. the average
c         gap model is based on moller and nix, nucl. phys.a536
c         (1992) 20. the average pair density is computed using
c         the semi-classical approximation  for   the  harmonic
c         oscillator model modified  by  the  isospin-dependent
c         factor derived based on the wigner-kirkwood expansion
c         for the woods-saxon potential (+ s.o. + coulomb).
c
c         input:    iz:    proton  number
c                   in:    neutron number
c                   iprot: 0 (neutrons) or 1 (protons)
c                   n1:    number of first h.f. level used
c                          in the bcs (ln) calculations.
c                          (preferably n1=1)
c                   n2:    number of last h.f. level used
c                          in the bcs (ln) calculations.
c                          (preferably n2=npart, i.e., iz or in)
c                   xmeff: effective mass (=1 in the woods
c                          saxon model)
c                   ibcs:  1 (bsc) or .ne.1 (ln)
c                   a0:    surface parameter in a Saxon Woods shape,
c                          is 0.58 * surface width from formfactor
c
c        output:    gavrg: pairing strength
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                   c
c      written by w. nazarewicz, 1994                               c
c                                                                   c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c      implicit real*8 (a-h,o-z)
      include 'paramr.inc'
c
      dimension alphap(2),bethap(2)

      data alphap /5.13,-2.35/,
     :     bethap /-1.28e-2,-2.38e-4/
c      data zero/0.0/
c
c----------------------------------------------------------------------
c
      a=iz+in
c
      ntot=iz
      if (iprot.eq.0) ntot=in
c
      if(rhosc.eq.zero) then
c
c          pair density at the fermi level
c
ctest      write(6,'(a/4i4,1f8.3,2i4,2f8.3)')
ctest     &   ' ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0',
ctest     &     ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0
       alpha=alphap(iprot+1)
       betha=bethap(iprot+1)
       fac=(2.*ntot/a)**(1.d0/3.)
       part3=(3.*ntot)**(2.0/3.)
       hom=41/a**(1.0/3)*fac
       tz=(in-iz)/a
       tz2=tz*tz
cold      rhosc=part3/hom/(1.-0.25/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
c4      rhosc=part3/hom/(1.-4.0/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
       rhosc=part3/hom/(1.-1.0/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
c
c                here: correction for the effective mass
c
       r=1.14*a**(1.d0/3)
       eta=a0/r
c
       xm=log(xmeff)
       pieta=(3.1415926*eta)**2
       xmm=xmeff*(1.-eta*xm)**3*(1.+2*eta*xm*pieta/(1.+pieta))
       rhosc=rhosc*xmm**1.5
c
      endif
      robar=rhosc
c

      if (ibcs.eq.1) then
c
c              bcs variant
c
         if (iprot.eq.0) then
            delavr=13.3/sqrt(a)
         elseif (iprot.eq.1) then
            delavr=13.9/sqrt(a)
         else
            delavr = 0.0D0
         endif

      else
c
c              ln variant
c
         if (iprot.eq.0) then
            delavr=8.6/sqrt(a)
         elseif (iprot.eq.1) then
            delavr=9.9/sqrt(a)
         else
            delavr = 0.0D0
         endif
c
      endif
c
      np=ntot/2
c
      y1=(-np+n1-1)/robar
      y2=(-np+n2  )/robar
c
      gavrg=2./robar/(log(sqrt(y2**2+delavr**2)+y2)-
     *                log(sqrt(y1**2+delavr**2)+y1))
c
      return
      end

c-----semrho------------ part of the hartree-fock package ----------------------
c
      function semrho(eref,potmf,bmasmf)
c
      include 'paramr.inc'
c
c     semiclassical level density
c     for the potential 'potmf'
c     with effective mass 'bmasmf'
c     taken at the reference energy 'eref'.
c     grid parameters enter via 'grid.inc'
c
c
      include 'dynpar.inc'
      include 'wavefr.inc'
c
      dimension potmf(0:nzd,0:nrd),bmasmf(0:nzd,0:nrd)
c
c
c-------------------------------------------------------------------------------
c
      do ir=0,nr
      do iz=0,nz
        bb     = one/bmasmf(iz,ir)
        eover  = eref-potmf(iz,ir)
        if(eover.gt.zero) then
          wrk4(iz,ir) = bb*sqrt(bb*eover)
        else
          wrk4(iz,ir) = zero
        endif
      enddo
      enddo
      semrho = f1int(wrk4)/(4.0*pi*pi)
c
      return
      end

c-----nucmat------- part of the Hartree-Fock package ------------------
c
      subroutine nucmat(rhonm,eovera,compr,effmas,asymme)
c
c     evaluates the stationary state for nuclear matter from
c     minimising the energy as function of the density.
c     
c     the force parameters are communicated via common/skypar/.
c
c     then a few nuclear matter parameters at the ninimum are
c     computed:
c       rhonm   = density
c       eovera  = energy per particle
c       compr   = compressibility
c       effmas  = effective mass
c
c     original version from joerg friedrich, rewritten 20.6.87
c     adapted for the SKYAX package 16.9.94
c
      include 'paramr.inc'
c
c      include 'skypar.inc'
      include 'dynpar.inc'
c      include 'grid.inc'
c      include 'wavfun.inc'
c     include 'constn.inc'
c
c      data pi/3.141592654/
cold      data rhoini,drhoit,epsil,itrmax
cold     *    / 0.16 , 0.01 ,1.e-5,    10/
      data rhoini,drhoit,itrmax
     *    / 0.16 , 0.01 ,    20/
      data dabe/1.0E-6/
c
c-----------------------------------------------------------------------
c
c     preparation of recoupled parameters
c     t12p = cm(krewald), t45p = cm1(krewald), diss. p.64, gl.(22,23)
c
      powp1  = power+1.0
      betab  = 3./5.*(3.*pi*pi/2.)**(2./3.)
      hqqd2m = 0.5D0*(h2m(1)+h2m(2))
c     hqqd2m = 197.3**2/2./936.
      t12p   = 3.*t1+5.*t2*(1.+4./5.*x2)
c      t45p   = (3.*t5+5.*t6*(1.+4./5.*x6))/6.0
      t45p   = zero
c
c     iteration to find the minimum energy
c
      eold   = 1.e10
      rhonm0 = rhoini
      do 20 iter=1,itrmax
      rhonmp = rhonm0+drhoit
      eplus  = hqqd2m*betab*rhonmp**(2./3.)+3./8.*t0*rhonmp
     &        +betab/16.*(t12p+t45p*rhonmp)*rhonmp**(5./3.)
     &        +t3/16.*rhonmp**(powp1)
      energy = hqqd2m*betab*rhonm0**(2./3.)+3./8.*t0*rhonm0
     &        +betab/16.*(t12p+t45p*rhonm0)*rhonm0**(5./3.)
     &        +t3/16.*rhonm0**(powp1)
      rhonmm = rhonm0-drhoit
      eminus = hqqd2m*betab*rhonmm**(2./3.)+3./8.*t0*rhonmm
     &        +betab/16.*(t12p+t45p*rhonmm)*rhonmm**(5./3.)
     &        +t3/16.*rhonmm**(powp1)
      rhonm0 = -0.5*drhoit*(eplus-eminus)/(eplus-energy-energy+eminus)
     *         +rhonm0
c      if(ifprnt.ge.0)
c     & write(02,'(a,i3,a,2g11.4,a,3g11.4)') ' it=',iter,' rho=',
c     &      rhonmm+drhoit,rhonm0,' e+0-=',eplus,energy,eminus
      if(abs(eold-energy).lt.dabe) goto 29
      eold   = energy
   20 continue
      write(06,'(a,i5,a)') ' no convergency after ',itrmax,
     *            ' steps in subroutine nucmat !'
c
   29 continue
c
c     now evaluate the nuclear matter parameters at the minimum
c     (e/a nach ergaenzt gemaess krewald, diss., p. 65, gl (25) )
c
      rhonm  = rhonm0
      fermp2 = (3.*pi*pi/2.*rhonm0)**(0.66666667)
      eovera = hqqd2m*betab*rhonm0**(2./3.)+3./8.*t0*rhonm0
     &        +betab/16.*(t12p+t45p*rhonm0)*rhonm0**(5./3.)
     &        +t3/16.*rhonm0**(powp1)
      compr  = 2.0*hqqd2m*betab*rhonm0**(2./3.)+9./4.*t0*rhonm0
     &        +5./4.*betab*(t12p+14./5.*t45p*rhonm0)*rhonm0**(5./3.)
     &        +9./16.*(powp1)*(power+2./3.)*t3*rhonm0**(powp1)
      effmas = 1./(1.+rhonm0/16./hqqd2m*(t12p+t45p*rhonm0))
      asymme = hqqd2m*fermp2/3.0 -t0/4.0*rhonm0*(0.5+x0)
     &         -t3/24.0*rhonm0**powp1*(0.5+x3)
     &         -t1*x1/8.0*rhonm0*fermp2
     &         +t2/6.0*rhonm0*fermp2*(1.0+1.25*x2)
c     &         -t5*x5/8.0*rhonm0*rhonm0*fermp2
c     &         +t6/6.0*rhonm0*rhonm0*fermp2*(0.5+x6)
c
      return
      end


*
* $RCSfile: liraxtest.f,v $
*
* $Revision: 1.10 $
* $Date: 1995/06/23 15:13:58 $
* $Author: bender $
* $Locker: bender $
* $State: Exp $
*
**********************************************************
*
* $Log: liraxtest.f,v $
* Revision 1.10  1995/06/23  15:13:58  bender
* Aufraeumarbeiten
*
* Revision 1.9  1995/05/11  18:31:11  bender
* Calculate Response-Densities with
* subroutine density
*
* Revision 1.8  1995/05/08  12:25:29  bender
* Neue Routine zur Berechnung der Dichten
*
* Revision 1.6  1995/04/20  18:20:48  bender
* Aufraeumarbeiten
*
* Revision 1.5  1995/04/12  13:18:05  bender
* Ausgabe ver"andert
*
* Revision 1.4  1995/04/07  18:27:21  bender
* Aufraeumarbeiten
*
* Revision 1.3  1995/04/03  11:54:22  bender
* Aufraeumarbeiten
*
* Revision 1.2  1995/03/31  16:42:25  bender
* write to coulomb.plt
*
* Revision 1.1  1995/03/31  13:03:41  bender
* Initial revision
*
*
**********************************************************
*
c      include 'density.f'
*
*-----respax----------------------------------------------------------*
*                                                                     *
      subroutine respax(iter)
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*     the contribution from the Skyrme residual interaction to        *
*     the "lambda_2" in Lipkin-Nogami Pairing.                        *
*                                                                     *
*                     -- axial symmetry --                            *
*                                                                     *
*     computes the response-densities to the particle number operator *
*     using pairing weights from 'weight'.                            *
*                                                                     *
*     The response Hamiltonian is constructed from the Skyrme-Energy- *
*     functional, for the details see                                 *
*     'Annalen der Physik, 1 (1992) 598 & 632'.                       *
*                                                                     *
*     the Lipkin-Nogami does only proton-proton pairing and           *
*     neutron-neutron pairing. thus only the Z^2 and N^2 response     *
*     is required.                                                    *
*     The procedure is extensively commented for the first term       *
*     'xpotp' and 'xpotn'. The other terms run analogously.           *
*                                                                     *
*---------------------------------------------------------------------*
*                                                                     *
*      Last modification: 03.04.1995                                  *
*                                                                     *
*----------------------------------------------------------------------
*                                                                     *
*     information is handled via common:                              *
*                                                                     *
*     results are returned via common                                 *
*                                                                     *
*      rspln    = U^2*V^2 response energy to be used in Lipk.Nogami   *
*      suv2r    = sum U^2*V^2, to be used in Lipkin-Nogami            *
*                                                                     *
***********************************************************************
*                                                                     *
*     constants and variables defined in SKYAX                        *
*     ----------------------------------------                        *
*                                                                     *
*-----mathem. and phys. constants (common /constv/ in paramr.inc)     *
*                                                                     *
*     oh = 0.5d0, one = 1.0d0, two = 2.0d0, edr = 1/3 d0              *
*                                                                     *
*-----parameter statments for setting dimensions:                     *
*     ( in common /????  in ???.inc )
*                                                                     *
*                  nzd    -   no. of mesh points in z direction.      *
*                  nrd    -   no. of mesh points in r direction.      *
*                  nsd    -   no. of spinor components.               *
*                  kqm    -   no. of isospins  (2)                    *
*                  npsid  -   maximum no. of wavefunctions            *
*                                                                     *
*-----geometrical and mesh parameters: ( in common /????/ )           *
*                                                                     *
*     dr     = step size in r-direction                               *
*     dz     = step size in z-direction                               *
*     sr     = increment in r-direction  ( r=(ir+sr)*dr )             *
*     sz     = increment in z-direction  ( z=(iz+sz)*dz )             *
*     drsqr  = dr*dr                                                  *
*     dzsqr  = dz*dz                                                  *
*     dvol   = 2*pi*drsqr*dz volume element for space integration     *
*     sym    = symmetry parameter,                                    *
*                    1 -> z-reflection symmetry                       *
*                    0 -> no z-reflection symmetry                    *
*     nz     - no. of z direction mesh points used in the
*              calculation. note nz < nzd.
*     nr     - no. of r direction mesh points used in the
*              calculation. note nr < nrd.
*     isym   = ifix(sym)                                              *
*     icm    = point on z-grid for center of the nucleus              *
*                    icm=1      for sym=0.0                           *
*                    icm=nz/2+1 for sym=1.0                           *
*                                                                     *
*        arrays defining the nucleus (running over iso=1,2):          *
*          nst    = total no. of states                               *
*          npart  = particle number                                   *
*          nlevel = number of occupied levels (=phase space for       *
*                   pairing                                           *
*          iokk   = pointer to the 1 st index of state with iso       *
*          ival   = number of valence particles                       *
*          lst    = total no. of states                               *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'inipar.inc'
      include 'observ.inc'
      include 'resden.inc'
*
      data uvlim/1.0D-8/              ! lower limit for inclusion in sum
*
*--------------------------------------------------------------------
*
*     compose combinations of skyrme parameters t_i, x_i
*
      b0     = t0*(one+oh*x0)
      b0p    = t0*(oh+x0)
      b3     = 0.25*t3*(one+oh*x3)
      b3p    = 0.25*t3*(oh+x3)
      b1     = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      b1p    = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      b2     = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
      b2p    = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
      b4     = oh*t4
      excfac = -1.418033*chgx    ! coulomb exchange in Slater appr.
*
*--------------------------------------------------------------------
*
*---- compute the response weights
*
      do 100 n=1,nst(1)+nst(2)
         v2 = vsqf(n)
         wgtres(n) = max(2.0 * (v2-v2*v2),smal)
  100 continue
*
*---- calculate sum u^2 v^2 for use in bcslip
*
      suv2rs(1) = zero
      suv2rs(2) = zero
      do n = 1, nst(1), 1
         v2 = vsqf(n)
         suv2rs(1) = suv2rs(1) + (v2-v2*v2)
      enddo
      do n = nst(1)+1, nst(1)+nst(2), 1
         v2 = vsqf(n)
         suv2rs(2) = suv2rs(2) + (v2-v2*v2)
      enddo
      suv2rs(1) = max(suv2rs(1),smal)
      suv2rs(2) = max(suv2rs(2),smal)
*
*     if(iprint.ge.1) then
*       write(8,113) suv2rs(1),suv2rs(2)
*113    format(1x,'suv2rs(p,n) :',2F12.4)
*     endif
*
*---- set response densities to zero ------------------------------------
*
      do 11 iz = 0, nz
         do 12 ir = 0, nr
            do 13 iq = 1, 2
               rrho(iz,ir,iq)   = zero
               rtau(iz,ir,iq)   = zero
               rrhols(iz,ir,iq) = zero
               rd2rho(iz,ir,iq) = zero
   13       continue
            rcoupt(iz,ir) = zero
   12    continue
   11 continue
*
*---- compute the response densities -------------------------------------
*
*     rrho   -->   < rho Delta N >
*     rtau   -->   < tau Delta N >
*     rrhols -->  < drhols Delta N >
*
*
      do ist = 1,nst(1)+nst(2),1
         iq = icharg(ist)
         if(wgtres(ist).gt.uvlim)
     &    call densit(psi(0,0,1,ist),wgtres(ist),ist,
     &                rrho(0,0,iq),rtau(0,0,iq),rrhols(0,0,iq),
     &                rholsr(0,0,iq),rholsz(0,0,iq) )
      enddo
*
*---- rd2rho -->  Laplace -> < rho Delta N >
*
      call laplace_scalar(rrho(0,0,1),rd2rho(0,0,1),+1,0,1,one)
      call laplace_scalar(rrho(0,0,2),rd2rho(0,0,2),+1,0,1,one)
cPK      call mulr(derr2a(0,0,0),rrho(0,0,1),rd2rho(0,0,1),one,0)
cPK      call mulr(derr2a(0,0,0),rrho(0,0,2),rd2rho(0,0,2),one,0)
cPK      call mulz(derz2(0,0,1),rrho(0,0,1),rd2rho(0,0,1),one,1)
cPK      call mulz(derz2(0,0,1),rrho(0,0,2),rd2rho(0,0,2),one,1)
*
*-------------------------------------------------------------------------
*
*---- terms depending on rrho * rrho
*
      b3t1 = b3*(power+two)*(power+one)*edr
      b3t2 = two*b3p*edr
      b3t3 = two*power*b3t2
      b3t4 = b3p*power*(power-one)*edr
*
      do ir=0,nr
         do iz=0,nz
            rp = rho(iz,ir,1)         ! proton-density
            rn = rho(iz,ir,2)         ! neutron-density
            rs = rp+rn                ! total density
            rpowm2 = rs**(power-two)  !
            rpowm1 = rs**(power-one)
            rpow   = rs**(power)
            tmp1 = (b0-b0p)
     &             +rpow*(b3t1-b3t2)
     &             -rpowm2*b3t4*((rp**2)+(rn**2))
*
c           wrk1(iz,ir) = rrho(iz,ir,1)*rrho(iz,ir,1)
c    &                  *( (b0-b0p)
c    &                   +rpowm2*( (b3t1-b3t2)*rs*rs
c    &                            -b3t3*rs*rp
c    &                            -b3t4*((rp**2)+(rn**2))))
*           wrk2(iz,ir) = rrho(iz,ir,2)*rrho(iz,ir,2)
*    &                  *( (b0-b0p)
*    &                   +rpowm2*( (b3t1-b3t2)*rs*rs
*    &                            -b3t3*rs*rn
*    &                            -b3t4*((rp**2)+(rn**2))))
*
            wrk1(iz,ir) = rrho(iz,ir,1)*rrho(iz,ir,1)
     &                  *(tmp1-rpowm1*b3t3*rp)
            wrk2(iz,ir) = rrho(iz,ir,2)*rrho(iz,ir,2)
     &                  *(tmp1-rpowm1*b3t3*rn)
*
         enddo
      enddo
c
      if (isym.eq.0) then
         do ir=0,nr
            wrk1(0,ir) = 0.0D0
            wrk2(0,ir) = 0.0D0
         enddo
      endif
c
      xpotp = f1int(wrk1)
      xpotn = f1int(wrk2)
*
*---- terms depending on rtau * rrho
*
      do ir=0,nr
         do iz=0,nz
            wrk1(iz,ir) = rrho(iz,ir,1)*rtau(iz,ir,1)
            wrk2(iz,ir) = rrho(iz,ir,2)*rtau(iz,ir,2)
         enddo
      enddo
      xkinp = two*(b1-b1p)*f1int(wrk1)
      xkinn = two*(b1-b1p)*f1int(wrk2)
*
*---- terms depending on rrho * (Laplace --> rrho)
*
      do ir=0,nr
         do iz=0,nz
            wrk1(iz,ir) = rrho(iz,ir,1)*rd2rho(iz,ir,1)
            wrk2(iz,ir) = rrho(iz,ir,2)*rd2rho(iz,ir,2)
         enddo
      enddo
      xpd2p = (b2p-b2)*f1int(wrk1)
      xpd2n = (b2p-b2)*f1int(wrk2)
*
*---- terms depending on rrho * rrhols
*
      do ir=0,nr
         do iz=0,nz
            wrk1(iz,ir) = rrho(iz,ir,1)*rrhols(iz,ir,1)
            wrk2(iz,ir) = rrho(iz,ir,2)*rrhols(iz,ir,2)
         enddo
      enddo
      xspop = -two*(b4+b4p)*f1int(wrk1)
      xspon = -two*(b4+b4p)*f1int(wrk2)
*
*---- Coulomb ---------------------------------------------------------------
*
      irepcs = irepc ! one should not manipulate input data
*                      like this, but this number has to be increased
*                      for this calculation
      irepc  = 200
c      icrflx = 1       ! flag for print of colrun properties
*
*---- save wcolt for use of colrun in gradstep
*
      do ir=0,nr
         do iz=0,nz
            xxx = wcolt(iz,ir)
            wcolt(iz,ir) = wcoltt(iz,ir)
            wcoltt(iz,ir) = xxx
         enddo
      enddo
*
*---- Coulomb-solver
*
      mangmo = 0
      call colrun(rrho(0,0,1),rcoupt(0,0),iter,mangmo)
*
*---- save wcolt in wcoltt for use of colrun in respax
*
      do ir=0,nr
         do iz=0,nz
            xxx = wcolt(iz,ir)
            wcolt(iz,ir) = wcoltt(iz,ir)
            wcoltt(iz,ir) = xxx
         enddo
      enddo
*
      do ir=0,nr
         do iz=0,nz
            wrk1(iz,ir) = rrho(iz,ir,1)*rcoupt(iz,ir)
         enddo
      enddo
      xcol = f1int(wrk1)
      irepc = irepcs ! reset irepc
c      icrflx = 0     ! reset icrflx
*
*---- Coulomb-Exchange
*
      twdthr = -(2.d0/3.d0)
      do ir=0,nr
         do iz=0,nz
            wrk1(iz,ir) = rrho(iz,ir,1) *rrho(iz,ir,1)
     &                    *(rho(iz,ir,1)**twdthr)
         enddo
      enddo
c
      if (isym.eq.0) then
         do ir=0,nr
            wrk1(0,ir) = 0.0D0
         enddo
      endif
      xcex = excfac*edr*f1int(wrk1)
*
*---- Summation der Beitraege -------------------------------------------
*
      rspln(1) = xpotp+xkinp+xpd2p+xspop+xcol+xcex
      rspln(2) = xpotn+xkinn+xpd2n+xspon
*
*---- write to response.plt
*
c      write(20,998) itert,xpotp,xkinp,xpd2p,xspop,
c     &              xcol,xcex,rspln(1),
c     &              xpotn,xkinn,xpd2n,xspon,rspln(2)
c  998 format(1x,1I4,12F9.3)
*
*
      return
      end
C#########################
C     CR 03-04-1998
C     CR 31-10-1998 now for skyax but first for rtmp=0
C
C     rbrent calculates efermi in pairdn for finite temperature using
C     zerofkt
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C Van Wijngaarden-Dekker-Brent Method for Root Finding
C
C W.H. Press et. al., Numerical Recipes in C, Sec. Ed, Cambridge 1995
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      FUNCTION RBRENT(E, DELTA, PH, GW, PAIRGW, RTMP, PARTNM, 
     &     NDIM, NMAX, nblock,ntcblock, IAB,IQ,tpairspace)
C
      include 'paramr.inc'
C
C----------------------------------------------------------------------------
C 
C      rtmp     = real temperature in MeV
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C
C      iab      = print level; significant values are -2, -1, 0 and 1;
C                 larger value prints more.
C
C      rxa,rxb,rxc,rlna,rlnb,rda,rdb are constants for finding roots
C                  of the zerofkt 
C
C
C      DATA FOR Wijngaarden-Dekker-Brent Method for Root
C
      LOGICAL,INTENT(IN) :: tpairspace
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
      DATA     RA0,   RB0,    REPS, IITMAX,    RTOL
     &     /-100.0, 100.0, 1.0e-14,    100, 1.0e-16/
C
C-----------------------------------------------------------------------------
C
C     SETTING UP
C


      RA  = RA0
      RB  = RB0
      RFA = partnm 
     &     - parnmt(E, DELTA, PH, GW, PAIRGW, RTMP, RA, NDIM, NMAX,
     &              nblock,ntcblock,IQ,tpairspace)
      RFB = partnm 
     &     - parnmt(E, DELTA, PH, GW, PAIRGW, RTMP, RB, NDIM, NMAX,
     &              nblock,ntcblock,IQ,tpairspace)
      RFC = RFB
      RC  = RB
      RD  = RB-RA
      RE  = RD
C     
C     TESTING
C
      IF ((RFA .GT. 0.0).AND.(RFB .GT. 0.0)) THEN
         IF (RFA.GT.RFB) THEN
            RB = 1.0
         ELSE
            RB = 0.0
         ENDIF
      ELSE IF ((RFA .LT. 0.0).AND.(RFB .LT. 0.0)) THEN
         IF (RFA.GT.RFB) THEN
            RB = 0.0
         ELSE
            RB = 1.0
         ENDIF
      ELSE
C
C     ITERATION
C
         DO II=1,IITMAX
C     RENANME RA,RB,RC AND AJUST BOUNDING INTERVAL RD
            IF (((RFB .GT. 0.0).AND.(RFC .GT. 0.0))
     &           .OR.((RFB .LT. 0.0).AND.(RFC .LT. 0.0))) THEN
               RC  = RA
               RFC = RFA
               RD  = RB-RA
               RE  = RD
            ENDIF
            IF (ABS(RFC).LT.ABS(RFB)) THEN
               RA  = RB
               RB  = RC
               RC  = RA
               RFA = RFB
               RFB = RFC
               RFC = RFA
            ENDIF
C     CONVERGENCE CHECK
            RTOL1 = 2.0*REPS*ABS(RB) + 0.5*RTOL
            RXM   = 0.5*(RC-RB)
            IF ((ABS(RXM).LE.RTOL1).OR.((ABS(RFB) .EQ. 0.0))) GOTO 2999
            IF ((ABS(RE).GE.RTOL1).OR.(ABS(RFA).GT.ABS(RFB))) THEN
C     ATTEMPT INVERSE QUADRATIC INTERPOLATION
               RS=RFB / RFA
               IF (RA .EQ. RC) THEN
                  RP = 2.0*RXM*RS
                  RQ = 1.0-RS
               ELSE
                  RQ = RFA / RFC
                  RR = RFB / RFC
                  RP = RS * (2.0*RXM*RQ*(RQ-RR)-(RB-RA)*(RR-1.0))
                  RQ = (RQ-1.0)*(RR-1.0) *(RS-1.0)
               ENDIF
C     CHECK WHETHER IN BOUNDS
               IF (RP .GT. 0.0) THEN
                  RQ = -RQ
               ENDIF
               RP    = ABS(RP)
               RMIN1 = 3.0*RXM*RQ-ABS(RTOL1*RQ)
               RMIN2 = ABS(RQ*RE)
               IF (2.0*RP .LT. MIN(RMIN1,RMIN2)) THEN
C     ACCEPT INTERPOLATION
                  RE = RD
                  RD = RP / RQ
C
                  IF(IAB .GT. 0) THEN
                     WRITE(6,'(A)') ' QUADRATIC INTERPOLATION '
                  ENDIF
C
               ELSE
C     INTERPOLATION FAILED -> BISECTION
                  RD = RXM
                  RE = RD
C
                  IF(IAB .GT. 0) THEN
                     WRITE(6,'(A)') ' FAILED -> BISECTION '
                  ENDIF
C
               ENDIF     
C     BOUNDS DECREASING TOO SLOWLY -> BISECTION
            ELSE
               RD = RXM
               RE = RD
C     
               IF(IAB .GT. 0) THEN
                  WRITE(6,'(A)') ' TOO SLOWLY -> BISECTION '
               ENDIF
C
            ENDIF
C     MOVE LAST BEST GUESS TO RA
            RA  = RB
            RFA = RFB
C     EVALUATE NEW TRIAL ROOT
            IF (ABS(RD) .GT. RTOL1) THEN
               RB = RB + RD
            ELSE
               RB = RB + SIGN(RTOL1,RXM)
            ENDIF
            RFBO = RFB
            RFB  = partnm 
     &            - parnmt( E, DELTA, PH, GW, PAIRGW, RTMP, RB, 
     &                      NDIM,NMAX,nblock,ntcblock,IQ,tpairspace)
C
            IF(IAB .GT. 0) THEN
               WRITE (6,'(I3,A,F28.18,A,F28.18,A,F28.18)')
     &              II,' :  rb = ',RB,'  rfb = ',RFB,'  drfb = ',
     &              (RFBO-RFB) 
            ENDIF
C
         ENDDO
         STOP ':( LOST'
 2999 ENDIF
C
      IF(IAB .GT. 0) THEN
         WRITE (6,'(A,F28.18)') 'EFERMI=NST= ', RB
      ENDIF
C
      RBRENT = RB
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     CR 03-04-1998
C     CR 31-10-1998
C
C     'parnmt' computes the particle number for given 'efermi'
C              in case of thermal BCS at temperature 'RTMP'
C
c
      FUNCTION parnmt(E, DELTA, PH, GW, PAIRGW, RTMP, EFERMI, 
     &                NDIM, NMAX, nblock, ntcblock,iq,tpairspace)
C
      include 'paramr.inc'
      include 'levels.inc'
C
C----------------------------------------------------------------------------
C 
C      rtmp     = real temperature in MeV
C      e        = array(1..nmax) of single particle energies.
C      delta    = array(1..nmax) of gap parameters.
C      ph       = array(1..nmax) of the multiplicities of the single
C                 particle levels.
C                 This multiplicities excludes (!) the spin degeneracy,
C                 e.g. a 1s1/2 state has ph = 1 .
C      pairgw   = space for dynamical cut-off factor
C      partnm   = desired particle number
C      efermi   = fermi-energy
c      nblock   = address of blocked state
C
C
      LOGICAL,INTENT(IN) :: tpairspace
      DIMENSION E(NDIM),PH(NDIM),GW(NDIM),DELTA(NDIM),PAIRGW(NDIM)
C
C-----------------------------------------------------------------------------

codd
!      call findconjpart(1)
!      call findconjpart(2)

codd
C
      PARTN = 0.0
C
      DO NA = 1, NMAX
C
         IF(tpairspace) THEN
           DELTA2 = DELTA(NA)
         ELSE
           DELTA2 = DELTA(NA)*PAIRGW(NA)
         END IF
         DELTA2 = DELTA2*DELTA2
coddAnfang
codd         if(isTSym.eq.1)then
          EDIF   = E(NA)-EFERMI 
codd         else
codd          NAP = itpartn(NA+(iq-1)*nst(1))-(iq-1)*nst(1)
codd          EDIF   = 0.5D0*(E(NA)+E(NAP)) - EFERMI
codd         endif
coddEnde
         EQUASI = SQRT(EDIF*EDIF+DELTA2)

C
c         IF (RTMP .EQ. 0.0) THEN
            if(NA.eq.nblock .AND. istsym==0) then
              GW(NA) = 0.5D0
            elseif(NA.eq.nblock) then
              GW(NA) = 1.0D0
            elseif(NA.eq.ntcblock) then
              GW(NA) = 1.0D-30
            else
              GW(NA) = MAX(0.5D0  * (1.0D0 - EDIF/EQUASI),v2lowlim)
              gw(na) = max(gw(na),smal)
              gw(na) = min(gw(na),one-smal)
            endif
c         ELSE
c            TEXP   = 1.0  / (EXP(EQUASI/RTMP) + 1.0)
c            GW(NA) = TEXP + (0.5 - TEXP) * (1.0 - EDIF/EQUASI)
c         ENDIF
C
codd        write(87,'(i3,a,i3,a,f7.3,a,f7.3,a,f8.3,a,
c     &   f8.3,a,f8.3,a,f8.3,a,f8.3)') 
c     & NA, ' ', NAP, ' ' , E(NA),' ', E(NAP), ' '
c     &,EDIF,' ',DELTA2,' ', GW(NA),' ', pairgw(NA), ' ', pairgw(NAP)



         IF(tpairspace) THEN
           PARTN  = 0.5*PH(NA)*GW(NA)*pairgw(na) + PARTN 
         ELSE
           PARTN  = 0.5*PH(NA)*GW(NA) + PARTN 
         END IF
      ENDDO
C
      parnmt  = PARTN+PARTN
!      IF(iq==2) WRITE(7,'(a,4i5,2f9.4)') 'CONJG:',
!     &   istblo(2),nblock+nst(1),itpartn(istblo(2)),ntcblock+nst(1),
!     &   gw(nblock),gw(ntcblock)
c      GW(nblock) = 1.0D0
      
C
      RETURN
      END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c-----pairinfo------------------------------------------------------
c
      subroutine pairinfo
c
      include 'paramr.inc'
c
c     output of beta2, gap,eff_gap,E_pair,E_eff,E_ges,E_eff_ges
c

c
c     computes the pair density
c
c       sum_n  u_n*v_n*|phi_n(r)|^2  ,
c
c     subsequently the pair potential on 'pairpt' from /wavfun/,
c     and finally the dynamic gaps on 'delatf' from /shells/
c     as expectation value of 'pairpt'.
c
c     'ipair' is used here to indicate the switch to surface pairing
c             5,6 --> volume pairing
c             7,8 --> surface pairing
c
c
      include 'wavefr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'levels.inc'
      include 'observ.inc'
c
      dimension  work(0:nzd,0:nrd)
c
      data rho0pr/0.16/               ! preliminary standard for surface pairing
      data half,xsmall/0.5,1.0E-20/

      

c
c--------------------------------------------------------------------
c
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
c
c     accumulate new pair-density
c
c      smallp = 0.000001
      smallp = 1D-40
      do jst=1,nst(1)+nst(2)
        iq    = icharg(jst)
        IF(jst==istblo(iq) .OR. jst==istcblock(iq)) CYCLE
        v2    = vsqf(jst)
        uv    = wtpair(jst)*sqrt(max(v2-v2*v2,smallp))*0.5
        call sdenti(psi(0,0,1,jst),paircr(0,0,iq),uv,jst)
      enddo
c
c     now multiply with strength to obtain local pair-potential
c
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
c
c     compute pairing energy 
c
      do ir=0,nr
      do iz=0,nz     
      work(iz,ir)=v0prot* paircr(iz,ir,1)* paircr(iz,ir,1)
      enddo
      enddo
      
      epairtestp=f1int(work)  !siehe oben
      
      do ir=0,nr
      do iz=0,nz     
      work(iz,ir)=v0neut* paircr(iz,ir,2)* paircr(iz,ir,2)
      enddo
      enddo
      
      epairtestn=f1int(work)  !siehe oben

      
      epairtest= epairtestp+ epairtestn



c
c     finally compute the actual gaps as s.p. expectation values with
c     the pair potential
c
      do jst=1,nst(1)+nst(2)
        iq     = icharg(jst)
        IF(jst==istblo(iq) .OR. jst==istcblock(iq)) THEN
          deltaf(jst) = 0D0
        ELSE
          do ir=0,nr
          do iz=0,nz
            wfa    = psi(iz,ir,1,jst)
            wfb    = psi(iz,ir,2,jst)
            work(iz,ir)  = (wfa*wfa+wfb*wfb)*pairpt(iz,ir,iq)
          enddo
          enddo
          deltaf(jst) = f1int(work)
          if(mpsi(jst).lt.0) deltaf(jst) = -deltaf(jst)
          deltaf(jst) = MIN(deltaf(jst),-1D-6)
        END IF
      enddo

c
c     compute average gap for protons
c
      sumuv  = 0.0
      sumduv = 0.0
      do na=1,nst(1)
        IF(jst==istblo(1) .OR. jst==istcblock(1)) CYCLE
        v2    = vsqf(na)
        uv    =sqrt(max(v2-v2*v2,smallp))*0.5
        sumuv  = uv*wate(na)                 + sumuv
        sumduv =wtpair(na)*uv*deltaf(na)*wate(na)      + sumduv
      enddo
      
      sumuv  = max(sumuv,xsmall)
      avdeltp =sumduv/sumuv
      epair2p =sumduv


c
c     compute average gap for neutrons 
c      
      sumuv  = 0.0
      sumduv = 0.0
      do na=nst(1)+1,nst(2)+nst(1)
        IF(jst==istblo(2) .OR. jst==istcblock(2)) CYCLE
        v2    = vsqf(na)
        uv    = sqrt(max(v2-v2*v2,smallp))*0.5
        sumuv  =uv*wate(na)                 + sumuv
        sumduv = wtpair(na)*uv*deltaf(na)*wate(na)      + sumduv
      enddo
      
      sumuv  = max(sumuv,xsmall)
      avdeltn =sumduv/sumuv
      
      
      write(987,'(7(1pg12.4))') 
     &         beta2,avdeltp,avdeltn,epairtestp,epairtestn,epair2p,
     &         sumduv
      
            


      etpair=epairtest

c      call energ()           ! energy directly from density
c      etot1  = energy
c      energy = energy-ezp+etpair

c      call energ2()         ! energy using rearrangement energy
c       energr = energr - ezp + etpair


      anumb  = npart(1)+npart(2)
      if(ifzpe.ne.0) then
       ezp    = 17.3/anumb**0.2
      else
       ezp    = 0.0
      endif

c      call energ2()             ! energy using rearrangement energy
c       energr = energr - ezp + etpair
      

      call energ()              ! energy directly from density
      etot1  = energy
      energy = energy+etpair-ezp



      write(88,'(1pg12.5,0pf13.5,0pf13.5)')
     &     beta2,energy






      return
      end
