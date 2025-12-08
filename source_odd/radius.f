c-----radius -----------------------------------------------------------
c
      subroutine radius()
c
c     Compute formfactor and deduced radius-observables               
c
      include 'paramr.inc'
c      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
c      include 'wavefr.inc'
c      include 'derivs.inc'
c
c
      parameter (kfrm=999)        ! maximum dimension of formfactor fields
c

      LOGICAL,PARAMETER :: trad=.TRUE.  ! switch to compute radii  and print

      dimension foc(kfrm),fop(kfrm),fon(kfrm),fo(kfrm),foc0(kfrm)
      dimension rhochr(kfrm),rhop(kfrm),rhon(kfrm)

c
cPK   parameter (delrad=0.15D0)       ! assumed spherical grid spacing
      parameter (delrad=0.3D0)        ! assumed spherical grid spacing
      REAL(8),PARAMETER :: cutk=4D0,cutdel=0.4D0    ! soft cutoff in k space
c
c----------------------------------------------------------------------
c
c     determine grid in spherical Fourier-Bessel space
c
      zmax   = ((nz)+sz-icm)*dz-zcm
      rmax   = (nr+sr)*dr
      rleng  = min(zmax,rmax)
      ikmax  = rleng/delrad+1
      if(ikmax.gt.kfrm) stop ' in RADIUS: KFRM too small'
      delrac = rleng/ikmax
      delk   =  PI/rleng
c      write(7,'(a,i5,3(1pg12.4))')
c     &  ' ikmax,rleng,delrac,delk=',ikmax,rleng,delrac,delk
c
c     width parameter for c.m. correction
c
      if(ifzpe.eq.0) then
       cmwid  = 1.58/(6.0*atnum**0.66666667)
      elseif(ifzpe<0) then
       cmwid=0D0
      else
        ezp = 17.3/anumb**0.2
        IF(trad) THEN
         zpeact=(widdz+widdxy)*0.5D0*(h2m(1)+h2m(2))/(npart(1)+npart(2))
         cmwid  = 0.375*h2mav/(atnum*zpeact)
        ELSE
         cmwid  = 0.375*h2mav/(atnum*ezp)
        END IF
      endif
      IF(trad) write(7,'(a,3(1pg12.4))') 
     &    ' zpeact(2),cmwid=',zpeact,ezp,cmwid
c
c     compute all formfactors
c
      write(7,'(a)') ' iq   q    formfactors '
      do iq=1,ikmax
        q1     = iq*delk-delk
        call fourf(q1,cmwid,foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq))
c        write(7,'(i5,f8.3,5(1pg12.4))')
c     &    iq,q1,foc0(iq),fop(iq),fon(iq),fo(iq),foc(iq)
      enddo
c      write(7,'(a)') ' test FBINT:'
c      do iq=1,4*ikmax
c        q1 = iq*delk*0.25
c        write(7,'(f8.4,1pg12.4)')
c     &    q1,fbint(q1,foc,delk,rleng,ikmax)
c      enddo
c
      IF(trad) THEN
c
c     finding first diffraction radii 'rdms..' at first zero
c     of the formfactor,
c
      rdmsc  = 4.493/qzero(foc,1,delk,rleng,ikmax)
c
c      evaluating surface width from height of first maximum.
c      the 'q1'="momentum_of_first_maximum" is estimated from
c      'rdmsc' as in a square-well density-distribution.
c
      q1     = 5.600/rdmsc
      qr     = q1*rdmsc
      surfc  = sqrt(abs(2./(q1*q1)*log(abs((sin(qr)/qr-cos(qr))
     &         *3.*apnum/(qr*qr*fbint(q1,foc,delk,rleng,ikmax))))))
c
c     rms-radii from curvature of formfactor at q=0
c     ---> extremely dangerous procedure due to insufficient interpolations
c     ----> values will be overwritten later by direct integration
c
c      dqt    = delk/10.0     !  delk/4.0
      dqt    = delk/2.0D0
      dqtd2  = 6.0D0/(dqt*dqt)
      rmsc   = sqrt(abs(dqtd2*(1.25D0*foc(1)
     &         -1.33333333D0*fbint(dqt,foc,delk,rleng,ikmax)
     &         +fbint((dqt+dqt),foc,delk,rleng,ikmax)/12.0D0)
     &         /foc(1)))
      END IF
c
c     charge density by fourier back-transformation
c
      dpi3   = 0.5D0/(pi*pi)
c                                                         cut high q
      WRITE(78,'(a)') '# formfactor'
      write(78,'(f8.3,5(1pg12.4))')  0D0,fop(1),fon(1),foc(1)
      do iq=2,ikmax
        qa     = delk*iq-delk
        IF(trad) write(78,'(f8.3,5(1pg12.4))')qa,fop(iq),fon(iq),foc(iq)
        cutbas = (qa-cutk)/cutdel
        cutexp = EXP(SIGN(cutbas**2,cutbas))
        foc(iq)=foc(iq)/(1D0+cutexp)
        fop(iq)=fop(iq)/(1D0+cutexp)
        fon(iq)=fon(iq)/(1D0+cutexp)
      enddo
      IF(trad) THEN
        WRITE(78,'(1x)')
        do iq=2,ikmax
          qa     = delk*iq-delk
          IF(trad) write(78,'(f8.3,5(1pg12.4))')
     &      qa,REAL(fop(iq),4),REAL(fon(iq),4),REAL(foc(iq),4)
        END DO
      END IF
c
c     first charge radius without magnetic contributions from 'foc0'
c     need only r>0
c     
      acc0=0D0
      accr=0D0
      r      = 0D0
      do ir=2,ikmax
        r      = delrac+r
c
        dc     = cos(r*delk)
        ds     = sin(r*delk)
        c      = 1D0
        s      = 0D0
        acc    = 0D0
        qa     = 0D0
        do iq=2,ikmax
          qa     = delk+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          acc    = qa*s*foc0(iq)+acc
        enddo
        rhochr(ir) = acc*delk*dpi3/r
        vol=r*r
        acc0=acc0+vol*rhochr(ir)
        accr=accr+vol*vol*rhochr(ir) 
      enddo
      rmsc0=SQRT(accr/acc0)
c
c     now full charge density and radius
c
      IF(trad) WRITE(78,'(1x/1x)') 
c                                                         r=0.0 separately
      IF(trad) WRITE(78,'(a)') 
     &    '#L=0: r  rho_proton   rho_neutron   rho_charge'
      acc    = 0D0
      accp   = 0D0
      accn   = 0D0
      do iq=2,ikmax
        qa     = delk*iq-delk
        acc    = qa*qa*foc(iq)+acc
        accp   = qa*qa*fop(iq)+accp
        accn   = qa*qa*fon(iq)+accn
      enddo
      rhochr(1) = (acc*delk*dpi3)
      rhop(1) = (accp*delk*dpi3)
      rhon(1) = (accn*delk*dpi3)
      IF(trad) WRITE(78,'(f6.2,3(1pg13.5))') 
     &     0D0,rhop(1),rhon(1),rhochr(1)
!      WRITE(7,'(l)') TDENSITY
!      IF(TDENSITY) WRITE(7,'(a)') '# r rhochr'
!      IF(TDENSITY) WRITE(7,'(f6.2,1pg13.5)') 0D0,rhochr(1)
c                                                         now all other r
      acc0=0D0
      accr=0D0
      accr4=0D0
      r      = 0D0
      do ir=2,ikmax
        r      = delrac+r
c
        dc     = cos(r*delk)
        ds     = sin(r*delk)
        c      = 1D0
        s      = 0D0
        acc    = 0D0
        accp   = 0D0
        accn   = 0D0
        qa     = 0D0
        do iq=2,ikmax
          qa     = delk+qa
          x      = c*dc-s*ds
          s      = s*dc+c*ds
          c      = x
          acc    = qa*s*foc(iq)+acc
          accp   = qa*s*fop(iq)+accp
          accn   = qa*s*fon(iq)+accn
        enddo
        rhochr(ir) = acc*delk*dpi3/r
        rhop(ir) = accp*delk*dpi3/r
        rhon(ir) = accn*delk*dpi3/r
        vol=r*r
        acc0=acc0+vol*rhochr(ir)
        accr=accr+vol*vol*rhochr(ir) 
        accr4=accr4+vol*vol*vol*rhochr(ir) 
        IF(trad) WRITE(78,'(f6.2,3(1pg13.5))') 
     &             r,rhop(ir),rhon(ir),rhochr(ir)
      enddo
      IF(trad) THEN
        CALL flush(78)
      ELSE
        RETURN
      END IF
      IF(trad) THEN
        WRITE(78,'(1x/1x/a)') '#L=2: formfactor'
        write(78,'(f8.3,5(1pg12.4))')  0D0,0D0,0D0
        fop(1)=0D0
        fon(1)=0D0
        do iq=2,ikmax
          qa     = delk*iq-delk
          CALL fourfquad(qa,cmwid,fop(iq),fon(iq))
          cutbas = (qa-cutk)/cutdel
          cutexp = EXP(SIGN(cutbas**2,cutbas))
          fop(iq)=fop(iq)/(1D0+cutexp)
          fon(iq)=fon(iq)/(1D0+cutexp)
!          fop(iq)=fop(iq)/(1D0+EXP(((qa-cutk)/cutdel)))
!          fon(iq)=fon(iq)/(1D0+EXP(((qa-cutk)/cutdel)))
          write(78,'(f8.3,5(1pg12.4))')  qa,fop(iq),fon(iq)
        enddo
!        CALL fourbquad(foc,rhochr,delk,delrac,ikmax)
        CALL fourbquad(fop,rhop,delk,delrac,ikmax)
        CALL fourbquad(fon,rhon,delk,delrac,ikmax)
        WRITE(78,'(1x/1x/a)') 
     &    '#L=2: r  rho_charge rho_proton   rho_neutron'
        WRITE(78,'(f6.2,3(1pg13.5))') 0D0,0D0,0D0
        r=0D0
        do ir=2,ikmax
          r      = delrac+r
          WRITE(78,'(f6.2,3(1pg13.5))') 
     &             r,rhop(ir),rhon(ir)
        enddo
        CALL flush(78)
      END IF
c
      rmscunsafe=rmsc
      rmsc=SQRT(accr/acc0)
      r4mom=(accr4/acc0)**0.25D0
      write(7,'(a,7(1pg13.5))')
     &  ' from charge formfactor: R,sigma,r,r0,r_unsafe,r4=',
     &   rdmsc,surfc,rmsc,rmsc0,rmscunsafe,r4mom,zpeact
      CALL flush(7)
c
      OPEN(71,file='densrad.res')
      WRITE(71,'(a/a)') '# radial proton densities',
     &   '#  r    L=0     L=2    L=4'
!      CALL test_intpol2D(rho(:,:,1))
      CALL radial_dens(rho(:,:,1))
      WRITE(71,'(/1x/a)') 'radial neutron densities'
      CALL radial_dens(rho(:,:,2))
      CLOSE(71)
c
      return
      end



c-----fourf ------------ part of the skyax package ----------------------
c
      subroutine fourf(q1,cmwid,foc0,fop,fon,fo,foc)
c
      include 'paramr.inc'
c
c     Computes the spherical Fourier-Bessel transform of the
c     proton and neutron densities. The formfactors thus obtained
c     are augmented by the nucleon formfactor.
c     The present version omits the contribution from the 
c     magnetic formfactor of the nucleon.
c     Input variables are:
c      q1    : radial momentum at which the fourier-transformation
c              is to be evaluated
c      cmwid : parameter for center of mass correction
c     The densities and the grid parameter are entered via common.
c     output paramaters:
c      foc0 :   charge formfactor
c      fop  :   proton "
c      fon  :   neutron    "
c      fo   :   total      "
c      foc  :   charge     "   including spin-orbit current
c
c     Further input and output is communicated through common blocks 
c     included via the '*.inc' files. Most important is the switch
c     for the l*s contribution to the formfactor (set in 'skyax.f'):
c      ifrso   =  1 --> charge density with contribution from ls current - old
c                -1 --> charge density with contribution from ls current - new
c                 0 --> charge density without contribution from ls current 
c
c     
c
      include 'fields.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  half/0.5D0/
c
c      the parameters in data are for the sachs-formfactors:
c       sa..,sm.. are for the isoscalar electric formfactor
c       va..,vm.. are for the isovector electric formfactor
c       gp..,gm.. are for the magnetic formfactors of proton and neutron
c
c       isoscalar and isovector formfactors are recoupled to
c       the electric formfactors of proton and neutron.
c
c     the sachs-formfactors are then recoupled to the nucleon
c     vector formfactor (form1..) and the nucleon tensor
c     formfactor (form2..).
c
c
c     parameters for scalar and vector electric formfactor
c     (walther, priv.comm., sept. 86)
c
      data    sa1,    sa2,    sa3,   sa4,   sm1,  sm2,  sm3,  sm4
     *    /2.2907,-0.6777,-0.7923,0.1793, 15.75,26.68,41.01,134.2/
      data    va1,    va2,    va3,   va4,   vm1,  vm2,  vm3,  vm4
     *    /0.3681, 1.2263,-0.6316,0.0372,  5.00,15.02,44.08,154.2/
c
c     magnetic formfactor , simon et al n.p. a333 (1980) 381
c
      data    gp1,    gp2,    gp3,   gp4,   gm1,  gm2,  gm3,  gm4
     *    / 0.694,  0.719,-0.418 ,0.005 ,  8.50,15.02,44.08,355.4/
c
c     old version of the
c     parameters for the magnetic formfactor (arenhoevel 1980)
c
c      data gp1,gp2,gp3,gp4,gm1,gm2,gm3,gm4/.794,.594,-.393,.005,
c     *                      8.73,15.02,44.08,355.40/
c
c     the anomanlous magnetic moments and the darwin factor 'd4m'
c
      data apro0,apro11,apro12,bpro11
     *    / 0.74, 0.0  ,  0.0 ,  0.0  /
      data aneu0,aneu11,aneu12,bneu11
     *    / 0.74, 0.0  ,  0.0 ,  0.0  /
      data xmgne,xmgpr,d4m/.04197,-.06131,.010989/
c
      LOGICAL :: trso    ! activte spin-orbit contribution
c
c----------------------------------------------------------------
c
      IF(abs(ifrso)>1) STOP 'IFRSO too large'  ! must be -1,0, or 1.
      trso = ABS(ifrso)==1
c
c     fourier-bessel transformation for proton- and neutron-density
c
      q2 = q1*q1
      do ir=0,nr
        r     = dr*(ir+sr)
        r2    = r*r
        do iz=0,nz
          z      = ((iz)+sz-icm)*dz
          rad2   = (z*z+r2)*q2
          if(rad2.LT.smal) then
            besfac = 1.0D0
          else
            rad2 = sqrt(rad2)
            besfac = sin(rad2)/rad2
          endif
          wrk1(iz,ir) = rho(iz,ir,1)*besfac
          wrk2(iz,ir) = rho(iz,ir,2)*besfac
c          write(7,'(2i5,6(1pg12.4))')
c     &      iz,ir,z,r,q1,rad2,besfac,wrk1(iz,ir)
        enddo
      enddo
      sump  = f1int(wrk1)
      sumn  = f1int(wrk2)
c
c     fourier-bessel transformation for proton- and neutron spin-orbit density
c
      IF(trso) THEN
        q2 = q1*q1
        do ir=0,nr
          r     = dr*(ir+sr)
          r2    = r*r
          do iz=0,nz
            z      = ((iz)+sz-icm)*dz
            rad2   = (z*z+r2)*q2
            if(rad2.LT.smal) then
              besfac = 1.0D0
            else
              rad2 = sqrt(rad2)
              besfac = sin(rad2)/rad2
            endif
            wrk1(iz,ir) = drhols(iz,ir,1)*besfac
            wrk2(iz,ir) = drhols(iz,ir,2)*besfac
c          write(7,'(2i5,6(1pg12.4))')
c     &      iz,ir,z,r,q1,rad2,besfac,wrk1(iz,ir)
          enddo
        enddo
        sumsp  = f1int(wrk1)
        sumsn  = f1int(wrk2)
      ELSE
        sumsp  = 0D0
        sumsn  = 0D0
      END IF
c
c      isoscalar and isovector electrical formactors on 'schses'
c      and 'schsev', recoupled then to proton and neutron electric
c      formfactors  'schsep' and 'schsen'. proton and neutron
c      magnetic formfactor on 'schsmp' and 'schsmn'.
c
       qa     = q1*q1
       schses = sa1/(1.+qa/sm1)+sa2/(1.+qa/sm2)+sa3/(1.+qa/sm3)
     &         +sa4/(1.+qa/sm4)
       schsev = va1/(1.+qa/vm1)+va2/(1.+qa/vm2)+va3/(1.+qa/vm3)
     &         +va4/(1.+qa/vm4)
       schsep = half*(schses+schsev)
       schsen = half*(schses-schsev)
       schsm  = ( gp1/(one+qa/gm1)+gp2/(one+qa/gm2)+gp3/(one+qa/gm3)
     &           +gp4/(one+qa/gm4) )
       IF(ifrso==1) THEN
         schsmp = schsm*xmgpr - schsep*d4m     ! old: yields mu+1/2
         schsmn = schsm*xmgne - schsen*d4m
       ELSE 
         schsmp = schsm*xmgpr + schsep*d4m     ! new: yields mu-1/2
         schsmn = schsm*xmgne + schsen*d4m
       END IF
c
       sumc   = sump*schsep+sumn*schsen
c
c      factor for darwin correction 'darfac'
c      and the factor for the c.m. correction 'exphf'.
c
       qa     = q1*q1
       if(trso) then
        darfac = one/sqrt(one+qa*d4m)
       else
        darfac = one
       endif
       exphf  = exp(cmwid*qa)
c
c     finally assembling the formfactors
c
      foc0   = exphf*sumc
      fop    = exphf*sump
      fon    = exphf*sumn
      fo     = exphf*(sump+sumn)
      foc    = exphf*(sumc*darfac + (schsmp*sumsp+schsmn*sumsn))
c
      return
      end

c-----fbint ------------ part of the hartree-fock package ----------------------
c
      double precision function fbint(qin,formfa,delk,rleng,ikmax)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns formfactor at the momentum 'qeff=q*rleng' by
c     fourier-bessel interpolation
c     from formfactors at fourier-grid points given on array 'formfa'
c     in the span '1..ikmax'.
c
      parameter (kfrm=999)        ! maximum dimension of formfactor fields
      dimension formfa(kfrm)
c
c      common /radpar/delk,rleng,ikmax           ! to communicate to subr.
c
c     'precis' defines the neighbourhood of a pole-term for which the
c     routine switches to an explicit laurent expansion.
c
      data precis/1.d-9/          ! /1.e-2/
      data  zero,one,pi/0.0d0,1.0d0,3.14159265358979d0/
c
c-----------------------------------------------------------------------
c
      if(ikmax.gt.kfrm) stop ' too large grid in FBINT '
c
      qeff   = qin*rleng
      xfx    = (qeff/pi)
      xfx2   = xfx*xfx
      ixfx   = (xfx+precis)
c
c     the fourier-bessel interpolation runs into problems with round-off
c     errors if 'qeff' comes too close to a multiple of 'pi'. we switch
c     to an explicit treatment of the pole-term in that case (i.e. first
c     case in the if-closure).
c
      if(abs(xfx-ixfx) .lt. precis) then
c                                             pole term separately
       ipole  = ixfx+1
       signum = -one
       sum    =  zero
       do i=2,ikmax
         signum = -signum
         if(i.ne.ipole)
     *   sum  = formfa(i)*signum/(one-xfx2/((i-1)*(i-1))) + sum
       enddo
       signum = (mod(ipole,2)*2-1)
       qdif   = (xfx-ixfx)
       sum2   = signum*formfa(ipole)*(ixfx*ixfx)/((ixfx+ixfx)+qdif)
       qdif   = qdif*pi
       sinq   = (one-qdif*qdif*0.166666667D0)*signum
       fbint  = 2.0D0*(sum*qdif+sum2*pi)*sinq/qeff
      else
c                                             normal case
       signum = -one
       sum    = zero
       do i=2,ikmax
         signum = -signum
         sum    = formfa(i)*signum/(one-xfx2/float((i-1)*(i-1))) + sum
       enddo
       fbint  = (sum+sum)*sin(qeff)/qeff
      endif
c
      return
      end
c-----qzero ------------ part of the hartree-fock package ----------------------
c
      double precision function qzero(foin,iqstrt,delk,rleng,ikmax)
c
      implicit double precision (a-h,o-z)
c
c     evaluates and returns the momentum at which the formfactor
c     given on array 'foin' has a zero.
c     scan for zero starts at the field position 'iqstrt'.
c     'ikmax' is the actual length of the array 'foin'.
c     'delk' is the stepsize in the momentum.
c
      parameter (kfrm=999)        ! maximum dimension of formfactor fields
      dimension foin(kfrm)
c
c      common /radpar/delk,rleng,ikmax           ! to communicate to subr.
c
c     some numerical parameters for termination of the search
c
      data itq0mx,endq0/3,0.0001D0/
      data  zero,one/0.0,1.0/
c
c---------------------------------------------------------------------
c
c     scan changing sign in 'foin'.
c      'iq1' is position before the zero and 'iq1' after the zero.
c
      ikmaxm = ikmax-1
      do 10 iscan=iqstrt,ikmaxm
      if(foin(iscan)*foin(iscan+1).le.zero) goto 19
   10 continue
   19 continue
      iq1   = iscan
      iq2   = iq1+1
      if(foin(iq1).le.zero) iq2   = iq1-1
      q1    = iq1*delk-delk
      q2    = iq2*delk-delk
      f1    = foin(iq1)
      f2    = foin(iq2+1) !GH (+1)
c      write(7,'(a,2i5,4(1pg12.4))')
c     &   ' iq1,iq2,q1,q2,f1,f2=',iq1,iq2,q1,q2,f1,f2
c
c     come closer to 'foin'=0.0 by repeated linear, interpolation.
c
      q3     = q1-f1*(q2-q1)/(f2-f1)
      do 20 itq0=1,itq0mx
      f3     = fbint(q3,foin,delk,rleng,ikmax)
c      write(7,'(a,2(1pg12.4))') 
c     &  ' q3,f3=',q3,f3
      if(f1*f3.gt.zero) then
       q1     = q3
       f1     = f3
      else
       q2     = q3
       f2     = f3
      endif
      q3old  = q3
      q3     = q1-f1*(q2-q1)/(f2-f1)
      if(abs(q3-q3old).lt.endq0) goto 29
   20 continue
   29 continue
c
c     improve position of the zero by quadratic interpolation:
c
c                         the coefficients of the interpolating binomial
c
      f1     = fbint(q1,foin,delk,rleng,ikmax)
      f2     = fbint(q2,foin,delk,rleng,ikmax)
      f3     = fbint(q3,foin,delk,rleng,ikmax)
      ca     = one/((q1-q2)*(q2-q3)*(q3-q1))
      c11     = ca*(q3*(q3-q2)*q2*f1+q1*(q1-q3)*q3*f2+q2*(q2-q1)*q1*f3)
      c2     = ca*((q2*q2-q3*q3)*f1+(q3*q3-q1*q1)*f2+(q1*q1-q2*q2)*f3)
      c3     = ((q3-q2)*f1+(q1-q3)*f2+(q2-q1)*f3)*ca
c      write(7,'(a,7(1pg12.4))')
c     &  ' f1,f2,f3,ca,c11,c2,c3=',f1,f2,f3,ca,c11,c2,c3
c
c                         the two zeroes of the interpolating binomial
c
      qa1    = ( sqrt(abs(c2*c2-4.0D0*c11*c3))-c2)/(c3+c3)
      qa2    = (-sqrt(abs(c2*c2-4.0D0*c11*c3))-c2)/(c3+c3)
c
c                         select the zero which is closer to q3
c
      q      = qa1
      if(abs(qa1-q3).gt.abs(qa2-q3)) q=qa2
c      write(7,'(a,2(1pg12.4))') ' qa1,qa2,q=',qa1,qa2,q
c
c     return zero thus found.
c
      qzero  = q
c
      return
      end


c-----fourquad ------------ part of the skyax package ----------------------
c
      subroutine fourfquad(q1,cmwid,fop,fon)
c
      include 'paramr.inc'
c
c     Computes the spherical Fourier-Bessel transform of the
c     proton and neutron densities for L=2 (quadrupole formfactor).
c     Input variables are:
c      q1    : radial momentum at which the fourier-transformation
c              is to be evaluated
c      cmwid : parameter for center of mass correction
c     The densities and the grid parameter are entered via common.
c     output paramaters:
c      fop  :   proton "
c      fon  :   neutron    "
c
      include 'fields.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
c
c      data  zero,one,pi/0.0,1.0,3.1415926/
      data  half/0.5D0/
      REAL(8),PARAMETER :: Y20fac=SQRT(5D0/(16D0*3.14159036D0))
c
c
c----------------------------------------------------------------
c
c     fourier-bessel transformation for proton- and neutron-density
c
      q2 = q1*q1
      do ir=0,nr
        r     = dr*(ir+sr)
        r2    = r*r
        do iz=0,nz
          z      = ((iz)+sz-icm)*dz
          rad2   = (z*z+r2)*q2
          if(rad2.LT.smal) then
            besfac = 0D0
          else
            besfac = bessel(rad2,2)*(2D0*z*z-r2)*Y20fac
          endif
          wrk1(iz,ir) = rho(iz,ir,1)*besfac
          wrk2(iz,ir) = rho(iz,ir,2)*besfac
c          write(7,'(2i5,6(1pg12.4))')
c     &      iz,ir,z,r,q1,rad2,besfac,wrk1(iz,ir)
        enddo
      enddo
      sump  = f1int(wrk1)
      sumn  = f1int(wrk2)
c
c     c.m. correction
c
      qa     = q1*q1
      exphf  = 1D0 !exp(cmwid*qa)
c
c     finally assembling the formfactors
c
      fop    = exphf*sump
      fon    = exphf*sumn
c
      return
      end

c-----bessel------------------------------------------------------
c
      function bessel(x,l)
c
      include 'paramr.inc'
c
c     spherical bessel function of angular momentum 'l' at
c     argument 'x'.
c
c     uses series expansion for x<l=1 and asymptotic expression
c     for x>l-1.
c
c     precision is tested only up to 5 digits.
c
!      include 'constn.inc'
c
      data xargmx/1.0D10/
c
c------------------------------------------------------------------------------
c
cc       if(l.gt.7) stop ' l>7 in "bessel"'
c
c      case of zero argument
c
       if(x.eq.0D0) then
        if(l.eq.0) then
         b      = 1D0
        else
         b      = 0D0
        endif
c
c      very large argument  --> set result to 0D0
c
       elseif(abs(x).gt.xargmx) then
         b      = 0D0
c
c      small arguments  --> use series expansion
c
       elseif(abs(x).le.l-1D0) then
        add   = 1D0
        do i=1,l
          add   = x/(2*i+1)*add
        enddo
        arg2  = -0.25D0*x*x
        b     = add
        do i=1,l+1
          add   = arg2*add/(i*(l+i+0.5D0))
          b     = add + b
        enddo
c
c      large arguments  --> use expansion in sin and cos
c
       else
        s      = sin(x)/x
        c      = cos(x)/x
        if(l.eq.0) then
         b      = s
        elseif(l.eq.1) then
         b      = s/x-c
        elseif(l.ge.2) then
         xinv   = 1D0/x
         bjm    = s
         b      = s/x-c
         do lrec=2,l
           bjp    = ((lrec+lrec)*xinv-xinv)*b-bjm
           bjm    = b
           b      = bjp
         enddo
         b      = b
        endif
       endif
c
       bessel = b
       return
       end



c-----fourbquad ------------ part of the skyax package ----------------------
c
      subroutine fourbquad(fourin,rhoout,delk,delrac,ikmax)
c
      include 'paramr.inc'
c
c     Fourier-Bessel backward for L=2
c
      parameter (kfrm=999)        ! maximum dimension of formfactor fields (align!)

      REAL(8),INTENT(IN) :: fourin(kfrm),delk,delrac
      INTEGER,INTENT(IN) :: ikmax
      REAL(8),INTENT(OUT) :: rhoout(kfrm)
c
      REAL(8),PARAMETER :: cutk=4D0,cutdel=0.4D0    ! soft cutoff in k space
      LOGICAL,PARAMETER :: trad=.FALSE.
c
      dpi3   = 0.5D0/(pi*pi)
c                                                         cut high q
      IF(trad) WRITE(78,'(a)') '# formfactor L=2'
      IF(trad) write(78,'(f8.3,5(1pg12.4))')  0D0,fourin(1)
!
!  damp high-q components  <--> to be done before
!      do iq=2,ikmax
!        qa     = delk*iq-delk
!        fourin(iq)=fourin(iq)/(1D0+EXP(((qa-cutk)/cutdel)))
!        IF(trad) write(78,'(f8.3,5(1pg12.4))')  qa,fourin(iq)
!      enddo
!
      IF(trad) WRITE(78,'(1x/1x)') 
c                                                         r=0.0 separately
      IF(trad) WRITE(78,'(a)') '# density L=2'
      rhoout(1) = 0D0
      IF(trad) WRITE(78,'(f6.2,1pg13.5)') 0D0,rhoout(1)
!      WRITE(7,'(l)') TDENSITY
!      IF(TDENSITY) WRITE(7,'(a)') '# r rhoout'
!      IF(TDENSITY) WRITE(7,'(f6.2,1pg13.5)') 0D0,rhoout(1)
c                                                         now all other r
      r      = 0D0
      do ir=2,ikmax
        r      = delrac+r
        acc    = 0D0
        qa     = 0D0
        do iq=2,ikmax
          qa     = delk+qa
          acc    = qa**2*fourin(iq)*bessel(r*qa,2)+acc
        enddo
        rhoout(ir) = acc*delk*dpi3/r
        IF(trad) WRITE(78,'(f6.2,1pg13.5)') r,rhoout(ir)
      enddo
c
      RETURN
      END

      include 'intpol2D.f'
