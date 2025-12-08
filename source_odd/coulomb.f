c=======================================================================
c
c
c     block #1: the package for the coulomb solver
c
c
c=======================================================================
c-----coulmb----------------------------------------------------------
c
      subroutine coulmb(iter,mangmo)
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
      include 'wavefr.inc'
c
      if(mangmo.GT.KANGMO) STOP ' COULMB: too large MANGMO'
      if (iter.eq.0) then
c       call colin1(crho,coulpt)
         if(mangmo.EQ.0) then
           call colin2(crho,coulpt)
         else
           do minit=0,KANGMO
             call colin3(crho,coulpt,mangmo)
           enddo
         endif
c         call testpr(coulpt,0,15,0,9,'coul.ini')
      endif
      call colrun(crho,coulpt,iter,mangmo)
      if(itest.gt.1 .and. mod(iter,iprint).eq.0)
     &     call testpr(coulpt,0,15,0,9,' coulpt')
c
      return
      end
c-----inicm2-----------------------------------------------------------
c
      subroutine inicms
c
c     initialises damping mass for coulomb potential
c     this variant is oriented towards the boundaries of the
c     box and carries the damping 2fm inside
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
*
*
      integer ir, iz
c      double precision drmin, drmin2, drr, dzadd, dradd,
c     &                 dtermz, dtermr, dz, dr, dqr, dz2, dr2,
c     &                 dzmax, drmax, ddmass, ddw, df
*
      data ddw/0.4/          ! width of damping
      data ddmass/1.0/       ! factor on damping mass
      data d2/2.0/           ! thickness of surface zone
      data dqr/0.0/          ! to make this variable static
*
      zmax    = ((nz-icm)+sz)*dz
      rmax    = ((nr)+sr)*dr
      drmin = min(zmax,rmax) - d2
      drmin2 = drmin * drmin
      dzadd = zmax - drmin - d2
      dradd = rmax - drmin - d2
      dtermr = drmin2 - dradd*dradd
      dtermz = drmin2 - dzadd*dzadd
*
      do iz = 0, nz
        z  = (iz+sz-icm)*dz
        z2 = z * z
        if( dzadd.ge.dradd ) then
          if (z.ge.dzadd) then
            dqr = sqrt(dtermz + d2*dzadd*z)
          else
            dqr = sqrt(z2 + drmin2)
          endif
        endif
        do ir = 0, nr
          r  = (ir+sr)*dr
          r2 = r * r
          drr = dsqrt(z2 + r2)
          if( dradd.gt.dzadd ) then
            if (r.ge.dradd) then
              dqr = sqrt(dtermr + d2*dradd*r)
            else
              dqr = sqrt(r2 + drmin2)
            endif
          endif
          sqrm(iz,ir) = ddmass / (one + exp((dqr-drr)/ddw))
        enddo
      enddo
      if(itest.gt.0) call testpr(sqrm,0,nz,0,nr,'sqrm')
*
      return
      end
c-----colrun------------------------------------------------------------
c
      subroutine colrun(rhoin,wcol,iter,mangmo)
c
c     this subroutine solves the poisson equation of
c     charged particles.
c
c
      include 'paramr.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'observ.inc'
c
      dimension wcol(0:nzd,0:nrd),rhoin(0:nzd,0:nrd),rhoa(0:nzd,0:nrd)
c
      data ephimx/0.0/     ! to make it static
c
      if (iter.le.10) then
        irepmx = 10*irepc
      else
        irepmx = irepc
      endif
      alpha   = 4.0*pi*chg
      dfac    = one/cmssqr
ctest      write(6,'(a,2g11.3)') ' colrun: cmssqr,dfac=',cmssqr,dfac
      call obs(rhoin,atnump,rmsp,dp,qp,op,hp,qsqp,zero)
      if(itest.ge.1) write(6,'(a,/4g16.6)')
     &  ' colrun with dp,qp,op,hp=',dp,qp,op,hp
        if(nhrun.eq.1) fach   = hp/chex
        if(nqrun.eq.1) facq   = qp/cquad
cbm
      if(isym.eq.0.and.ndrun.eq.1) then
        facd   = dp/cdip
      else
        facd   = zero
      endif
      if(isym.eq.0.and.nhrun.eq.1) then
        faco   = op/coct
      else
        faco   = zero
      endif
      facm   = atnump/cmono
c
      if (iter.eq.0) then
       do ir=0,nr
       do iz=0,nz
          rhoa(iz,ir)  = alpha*(rhoin(iz,ir)
     &                   -facm*crhom(iz,ir,mangmo))
       enddo
       enddo
      else if (isym.eq.1) then
        if (nqrun.eq.0 .and. nhrun.eq.0) then
         do ir=0,nr
         do iz=0,nz
            rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                   -facm*crhom(iz,ir,mangmo))
            wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else if (nqrun.eq.1 .and. nhrun.eq.0) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                       -facm*crhom(iz,ir,mangmo)
     &                       -facq*crhoq(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else if (nqrun.eq.1 .and. nhrun.eq.1) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                     -facm*crhom(iz,ir,mangmo)
     &                     -facq*crhoq(iz,ir,mangmo)
     &                     -fach*crhoh(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else
         stop ' no such case for "nqrun" or "nhrun"'
        endif
      else if (isym.eq.0) then
        if (ndrun.eq.0) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                      -facm*crhom(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else if (ndrun.eq.1 .and. nqrun.eq.0) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                      -facm*crhom(iz,ir,mangmo)
     &                      -facd*crhod(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else if (nqrun.eq.1 .and. norun.eq.0) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                      -facm*crhom(iz,ir,mangmo)
     &                      -facd*crhod(iz,ir,mangmo)
     &                      -facq*crhoq(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else if (norun.eq.1 .and. nhrun.eq.0) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                   -facm*crhom(iz,ir,mangmo)
     &                   - facd*crhod(iz,ir,mangmo)
     &                   - facq*crhoq(iz,ir,mangmo)
     &                   - faco*crhoo(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        else if (nhrun.eq.1) then
         do ir=0,nr
         do iz=0,nz
           rhoa(iz,ir) = alpha*(rhoin(iz,ir)
     &                   - facm*crhom(iz,ir,mangmo)
     &                   - facd*crhod(iz,ir,mangmo)
     &                   - facq*crhoq(iz,ir,mangmo)
     &                   - faco*crhoo(iz,ir,mangmo)
     &                   - fach*crhoh(iz,ir,mangmo))
           wcol(iz,ir) = wcolt(iz,ir)
         enddo
         enddo
        endif
      endif
c
      do 40 irep=1,irepmx
c     if (mod(irep,20).eq.0) then
c     write(6,'(a,i3)') ' iteration nr=',irep-1
c     call testpr(wcol,0,9,0,9,' wcol  ')
c     call testpr(wcol,nz-9,nz,nr-9,nr,' wcol  ')
c     endif
c
cPK      do ir=0,nr
cPK      do iz=0,nz
cPK        wcolt(iz,ir) = sqrm(iz,ir)*wcol(iz,ir)
cPK      enddo
cPK      enddo
c

      call laplace_scalar(wcol,wcolt,+1,mangmo,1,-one)
!      WRITE(*,*) 'LAPLACE:',f1int(wcol),f1int(wcolt)
cPK      call mulz(derz2(0,0,1),wcol,wcolt,-one,1,iband_derz)
cPK      call mulr(derr2a(0,0,mangmo),wcol,wcolt,-one,1,iband_derr)

c
      do ir=0,nr
      do iz=0,nz
         wcolt(iz,ir) = wcolt(iz,ir) +
     &                  sqrm(iz,ir)*wcol(iz,ir) - rhoa(iz,ir)
      enddo
      enddo
c
      call damp_meson(wcolt,wcolt,one)
!      WRITE(*,*) 'DAMP:',f1int(wcolt),f1int(zdmes(0,0,1)),
!     &    f1int(rdmes(0,0,1))
cPK   call mulz(zdmes(0,0,1),wcolt,wrk1,one,0,-1)
cPK   call mulr(rdmes(0,0,1),wrk1,wcolt,one,0,-1)
c
      ephimx = -one
      do ir=0,nr
      do iz=0,nz
       ephi    = wcolt(iz,ir)
       if (abs(ephi).gt.ephimx) then
         ephimx = abs(ephi)
       endif
       wcol(iz,ir) = wcol(iz,ir) - dfac*ephi
      enddo
      enddo

c     if (mod(irep,1).eq.0) then
c     write(6,'(/a,g15.5,3(a,i3)/)') ' max error=',ephimx,' rep=',irep,
c    &                               ' ir=',irmx,' iz=',izmx
c     call testpr(wcol,0,9,0,9,' wcol  ')
c     endif
c
c     if (abs(ephimx).le.colerr.or.irep.eq.irepc) then
      if (abs(ephimx).le.colerr) then
c     write(6,'(/a,g15.5,3(a,i3)/)') ' max error=',ephimx,' rep=',irep,
c    &                               ' ir=',irmx,' iz=',izmx
c     call testpr(wcol, 0, 9, 0, 9,' wcol  ')
c     call testpr(wcol, 0, 9,10,19,' wcol  ')
c     call testpr(wcol,10,19, 0, 9,' wcol  ')
c     call testpr(wcol,10,19,10,19,' wcol  ')
      goto 70
      endif
 40   continue
c
 70   continue
ccc
c      write(6,*) irep
cbm
      if (isym.eq.1) then
        if (nqrun.eq.0.or.iter.eq.0) then
          do ir=0,nr
          do iz=0,nz
           wcolt(iz,ir)= wcol(iz,ir)
           wcol (iz,ir)= wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
          enddo
          enddo
        else if (nqrun.eq.1.and.nhrun.eq.0) then
          do ir=0,nr
          do iz=0,nz
           wcolt(iz,ir)= wcol(iz,ir)
           wcol(iz,ir) = wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
     &              + facq*coulq(iz,ir,mangmo)
          enddo
          enddo
        else if (nqrun.eq.1.and.nhrun.eq.1) then
          do ir=0,nr
          do iz=0,nz
            wcolt(iz,ir)= wcol(iz,ir)
            wcol(iz,ir) = wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
     &                   + facq*coulq(iz,ir,mangmo) 
     &                   + fach*coulh(iz,ir,mangmo)
          enddo
          enddo
        endif
      else if (isym.eq.0) then
        if (ndrun.eq.0.or.iter.eq.0) then
          do ir=0,nr
          do iz=0,nz
            wcolt(iz,ir)= wcol(iz,ir)
            wcol (iz,ir)= wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
          enddo
          enddo
        else if (ndrun.eq.1.and.nqrun.eq.0) then
          do ir=0,nr
          do iz=0,nz
            wcolt(iz,ir)= wcol(iz,ir)
            wcol (iz,ir)= wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
     &                   + facd*could(iz,ir,mangmo)
          enddo
          enddo
        else if (nqrun.eq.1.and.norun.eq.0) then
          do ir=0,nr
          do iz=0,nz
            wcolt(iz,ir)= wcol(iz,ir)
            wcol (iz,ir)= wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
     &                   + facd*could(iz,ir,mangmo)
     &                   + facq*coulq(iz,ir,mangmo)
          enddo
          enddo
        else if (norun.eq.1.and.nhrun.eq.0) then
          do ir=0,nr
          do iz=0,nz
            wcolt(iz,ir)= wcol(iz,ir)
            wcol (iz,ir)= wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
     &                   + facd*could(iz,ir,mangmo)
     &                   + facq*coulq(iz,ir,mangmo)
     &                   + faco*coulo(iz,ir,mangmo)
          enddo
          enddo
        else if (nhrun.eq.1) then
          do ir=0,nr
          do iz=0,nz
            wcolt(iz,ir)= wcol(iz,ir)
            wcol (iz,ir)= wcol(iz,ir) + facm*coulm(iz,ir,mangmo)
     &                   + facd*could(iz,ir,mangmo)
     &                   + facq*coulq(iz,ir,mangmo)
     &                   + faco*coulo(iz,ir,mangmo)
     &                   + fach*coulh(iz,ir,mangmo)
          enddo
          enddo
        endif
        if (itest.ge.1) then
           write(6,'(a,5g15.5)') 'fm,fd,fq,fo,fh =',facm,facd,facq,faco,
     &          fach
        endif
      endif
c
c      if (abs(ephimx).ge.colerr.and.iter.gt.30) then
      if (abs(ephimx).ge.colerr) then
c       write(6,*)colerr,abs(ephimx),irepmx
       write(6,'(/2(a,g15.5),i3,a/)') ' colrun: max error=',
     &                   colerr,' actual error=',abs(ephimx),
     &                   irepmx,' iterations'
c     call testpr(wcol,0,9,0,9,' wcol  ')
      endif
      if(itest.ge.1)
     &  call coldif(wcol,rhoin,wcolt,' full solution')
c
      return
      end
c-------coldif---------------------------------------------------------
c
      subroutine coldif(colin,rhoin,work,text)
c
c     tests error in poisson equation
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
c
      character*(*) text
c
      dimension  colin(0:nzd,0:nrd),rhoin(0:nzd,0:nrd)
      dimension  work(0:nzd,0:nrd)
c
      write(6,'(2a)')
     &   ' test poisson equation for field:',text
c
      call laplace_scalar(colin,work,+1,0,1,-one)
cPK   call mulz(derz2,colin,work,-one,0)
cPK   call mulr(derr2a(0,0,0),colin,work,-one,1)
      alpha   = 4.0*pi*chg
      errav   = zero
      errwgt  = zero
      errmax  = zero
      irmax   = 0
      izmax   = 0
      do ir=0,nr
      do iz=0,nz
        rhowgt = rhoin(iz,ir)
        differ  = abs(work(iz,ir)-alpha*rhowgt)
        work(iz,ir) = differ
        if(differ.gt.errmax) then
         irmax  = ir
         izmax  = iz
         errmax = differ
        endif
        errwgt = rhowgt + errwgt
        errav  = differ*differ*rhowgt+errav
      enddo
      enddo
      errav  = errav/errwgt
c
      write(6,'(a,g12.4,a,g12.4,a,2i5)')
     &   ' average error=',sqrt(abs(errav)),
     &   ' max. error=',errmax,' at ir,iz=',irmax,izmax
      if(itest.ge.2) then
       write(6,'(a,20(/1x,10(1pg11.3)))')
     &  ' error along r-axis:',(work(0,ir),ir=0,nr)
       write(6,'(a,20(/1x,10(1pg11.3)))')
     &  ' error along z-axis:',(work(iz,0),iz=0,nz)
      endif
c
      return
      end
c-------coldf2--------------------------------------------------------
c
      subroutine coldf2(colin,rhoin,work,text)
c
c     tests error in poisson equation
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'levels.inc'
c
      character*(*) text
c
      dimension  colin(0:nzd,0:nrd),rhoin(0:nzd,0:nrd)
      dimension  work(0:nzd,0:nrd)
c
      open(16,file='coldif')
      write(16,'(2a)')
     &   ' test poisson equation for field:',text
      write(16,*) ' z      r       Lapl.Phi        4PiRho       Differ'
c
      call laplace_scalar(colin,work,+1,0,1,-one)
cPK      call mulz(derz2,colin,work,-one,0)
cPK      call mulr(derr2a(0,0,0),colin,work,-one,1)
      alpha   = 4.0*pi*chg
      do ir=0,0 !nr
         r     = (ir+sr)*dr
         do iz=0,nz
            z     = (iz+sz-icm)*dz
            differ  = abs(work(iz,ir)-alpha*rhoin(iz,ir))
            write(16,'(2f6.2,3g15.6)') z,r,
     &           work(iz,ir),alpha*rhoin(iz,ir),differ
         enddo
      enddo
      close(16)
c
      return
      end
c-----colin2----------------------------------------------------------
c
      subroutine colin2(rhoini,wcol)
c
c     initialises the monpole and quadrupole background fields for
c     the subroutine colrun, here using solution by radial integration
c     of the green's function
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
      common /colinp/ am2,rad02,aq2,rad04,ah2,radpow,          ! for sym=1
     &                ad2,rad0,ao2,rad03                       ! for sym=0
c
      dimension  wcol(0:nzd,0:nrd),rhoini(0:nzd,0:nrd)
      external fm1,fq1,fh1,fd1,fo1
c
c     initialisation of wcol
c
      do ir=0,nr
      do iz=0,nz
        wcol(iz,ir) = rhoini(iz,ir)
        wcol(iz,ir) = 0.0
      enddo
      enddo
c
c     general parameters for background fields
c
      zmax    = ((nz-icm)+sz)*dz
      rmax    = ((nr)+sr)*dr
c      radmax  = sqrt(zmax*zmax+rmax*rmax)
      zrmin   = min(zmax,rmax)
      zrmax   = max(zmax,rmax)
      binup   = zrmax+5.0d0
      fzrm    = zrmin/5.5d0
      rad0    = zrmax
      rad02   = rad0 *rad0
c
c     monopole field and density
c
cold      pro     = (npart-njell)
      pro     = one
      am      = 6.0*sqrt(-log(oh))/zrmin
      am2     = am*am
      am3     = am2*am
      bm      = two*am3*pro/pi
      pifac   = two*sqrt(pi)
      dm0     = bm/pifac
      fm0     = bm*pifac*chg
      do ir=0,nr
      r = (ir+sr)*dr
      do iz=0,nz
         z     = (iz+sz-icm)*dz
         rad2  = z*z+r*r
         rad   = sqrt(rad2)
         radpow = zero                              ! setting for 'fm1'
         crhom(iz,ir,0) = dm0*fm1(rad)
c
c        the  subroutine dglnq(f,x1,x2,o) computes the integral
c           int_{x1}^{x2} f(t) dt
c        with gauss legendre of order o.
c        the function f has to be provided as subroutine.
c
         bin1 = am * rad
         if (bin1.lt.1.d-2) then
          qt1 = zero
         else
          radpow = two                               ! setting for 'fm1'
          call qgaus (fm1, zero, rad, qt1)
          qt1 = qt1/(rad)
         endif
         radpow = one
         call qgaus(fm1, rad , binup, qt2)
         coulm(iz,ir,0) = fm0*(qt1+qt2)
      enddo
      enddo
      call obs(crhom,cmono,rmsm,cmdip,cmquad,cmoct,cmhex,qsqcm,zero)
           write(6,'(a/3x,7(1pg11.3))')
     &  ' from "crhom": cmono,rmsm,cmdip,cmquad,cmoct,cmhex,qsqcm=',
     &  cmono,rmsm,cmdip,cmquad,cmoct,cmhex,qsqcm
      if(itest.ge.1) then
        call coldif(coulm,crhom,wcolt,' monopole background')
        call coldf2(coulm,crhom,wcolt,' monopole background')
      endif
c     do 10 ir=0,nr
c     do 10 iz=0,nz
c       wcol(iz,ir) = coulm(iz,ir)
c10   continue
c
c     call testpr(coulm,0,9,0,9,' coulm')
c     call testpr(crhom,0,9,0,9,' crhom')
c
c     optionally initialize the quadrupole fields
c
      if(nqrun.eq.1) then
        qp      = one
        aq      = sqrt( -log( oh ) )/fzrm
        aq2     = aq*aq
c        aq3     = aq2*aq
c        aq7     = aq3*aq3*aq
        radpow = 6.0
        call qgaus(fq1, zero, binup, qinteg)
        pifac   = 0.25 * sqrt(5.0/pi)
        bq      = qp / qinteg * pifac
        dq0     = bq
        fq0     = bq*0.8d0*pi*chg
        pfy     = pifac
        do ir=0,nr
        r = (ir+sr)*dr
        do iz=0,nz
          z     = (iz+sz-icm)*dz
          rad2  = z*z+r*r
          rad   = sqrt(rad2)
          rad5  = rad2*rad2*rad
          y20p  = pfy*(two*z*z-r*r)
          radpow = zero
          crhoq(iz,ir,0) = dq0*fq1(rad)*y20p
c
c        the  subroutine dglnq(f,x1,x2,o) computes the integral
c           int_{x1}^{x2} f(t) dt
c        with gauss legendre of order o.
c        the function f has to be provided as subroutine.
c
          if (rad.lt.1.d-2) then
           qt1 = zero
          else
           radpow = 6.0
           call qgaus(fq1, zero, rad, qt1)
           qt1 = qt1/rad5
          endif
          radpow = one
          call qgaus(fq1, rad , binup, qqt2)
          coulq(iz,ir,0) = fq0*(qt1+qqt2)*y20p
        enddo
        enddo
        call obs(crhoq,cqmono,rmsq,cqdip,cquad,cqoct,cqhex,qsqcq,
     &           zero)
        write(6,'(a/3x,7(1pg11.3))')
     &   ' from "crhoq": cqmono,rmsq,cqdip,cquad,cqoct,cqhex,qsqcq=',
     &    cqmono,rmsq,cqdip,cquad,cqoct,cqhex,qsqcq
        if(itest.ge.1)
     &   call coldif(coulq,crhoq,wcolt,' quadrupole background')
      else
        do ir=0,nr
        do iz=0,nz
          coulq(iz,ir,0) = 0.0
          crhoq(iz,ir,0) = 0.0
        enddo
        enddo
       endif
c
c     optionally initialize the hexadecupole fields
c
      if(nhrun.eq.1) then
        hp      = one
        rad04   = rad02*rad02
        ah      = sqrt( -log( oh ) )/fzrm
        ah2     = ah*ah
c        ah3     = ah2*ah
c        ah7     = ah3*ah3*ah
        radpow = 10.0
        call qgaus(fh1, zero, binup, hint)
        bh      = hp / hint
        dh0     = bh
        fh0     = bh*4.0d0/9.0d0*pi*chg
        pfy     = 3.0 / 16.0 / sqrt(pi)
        do ir=0,nr
        r = (ir+sr)*dr
        do iz=0,nz
         z     = (iz+sz-icm)*dz
         rad2  = z*z+r*r
         rad4  = rad2*rad2
         rad   = sqrt(rad2)
         rad5  = rad4*rad
         rad9  = rad4*rad5
         y40p  = pfy*(8.0d0*z*z*z*z + 3.0d0 *r*r*r*r - 24.0d0 *z*z*r*r)
         radpow = zero
         crhoh(iz,ir,0) = dh0*fh1(rad)*y40p
c
c        the  subroutine dglnq(f,x1,x2,o) computes the integral
c           int_{x1}^{x2} f(t) dt
c        with gauss legendre of order o.
c        the function f has to be provided as subroutine.
c
         if (rad.lt.1.d-2) then
          qt1 = zero
         else
          radpow = 10.0
          call qgaus(fh1, zero, rad, qt1)
          qt1 = qt1/rad9
         endif
         radpow = one
         call qgaus(fh1, rad , binup, qqt2)
         coulh(iz,ir,0) = fh0*(qt1+qqt2)*y40p
        enddo
        enddo
        call obs(crhoh,chmono,rmsh,chdip,chquad,choct,chex,qsqch,
     &           zero)
        write(6,'(a/3x,7(1pg11.3))')
     &  ' from "crhoh": chmono,rmsh,chdip,chquad,choct,chex,qsqch=',
     &    chmono,rmsh,chdip,chquad,choct,chex,qsqch
        if(itest.ge.1)
     &    call coldif(coulh,crhoh,wcolt,' hexadecupole background')
      else
        do ir=0,nr
        do iz=0,nz
          coulh(iz,ir,0) = 0.0
          crhoh(iz,ir,0) = 0.0
        enddo
        enddo
      endif
cbm
cbm  initialize dipole and octupole fields in case of isym = 0
cbm
      if (isym.eq.0) then
cbm
cbm   optionally initialize the dipole fields
cbm
        if (ndrun.eq.1) then
          dp      = one
          ad      = 1.35*sqrt( -log( oh ) )/fzrm
          ad2     = ad * ad
          radpow  = 4.0
          call qgaus(fd1, zero, binup, dinteg)
          dh      = dp / dinteg
          fd0     = dh*4.0d0/3.0d0*pi*chg
          pfy     = oh * sqrt(3./pi)
          do ir=0,nr
          r = (ir+sr)*dr
          do iz=0,nz
             z     = (iz+sz-icm)*dz
             rad2  = z*z+r*r
             rad   = sqrt(rad2)
             rad3  = rad2 * rad
             y10p  = pfy * z
             radpow = zero
             crhod(iz,ir,0) = dh*fd1(rad)*y10p
c   mind periodic boundaries:
             if (iz.eq.nz) crhod(iz,ir,0) = zero
             if (iz.eq.0.and.isym.eq.0) crhod(iz,ir,0) = zero
c
c            the  subroutine dglnq(f,x1,x2,o) computes the integral
c               int_{x1}^{x2} f(t) dt
c            with gauss legendre of order o.
c            the function f has to be provided as subroutine.
c
             if (rad.lt.1.d-2) then
              qt1 = zero
             else
              radpow = 4.0
              call qgaus(fd1, zero, rad, qt1)
              qt1 = qt1/rad3
             endif
             radpow = one
             call qgaus(fd1, rad , binup, qt2)
             could(iz,ir,0) = fd0*(qt1+qt2)*y10p
c   mind periodic boundaries:
             if (iz.eq.nz) could(iz,ir,0) = zero
             if (iz.eq.0.and.isym.eq.0) could(iz,ir,0) = zero
          enddo
          enddo
          call obs(crhod,cdmono,rmsd,cdip,cdquad,cdoct,cdhex,qsqcd,
     &             zero)
          write(6,'(a/3x,7(1pg11.3))')
     &    ' from "crhod": cdmono,rmsd,cdip,cdquad,cdoct,cdhex,qsqcd=',
     &    cdmono,rmsd,cdip,cdquad,cdoct,cdhex,qsqcd
        if(itest.ge.1) then
           call coldif(could,crhod,wcolt,' dipole background')
c          call coldf2(could,crhod,wcolt,' dipole background')
         endif
        else
          do ir=0,nr
          do iz=0,nz
            could(iz,ir,0) = 0.0
            crhod(iz,ir,0) = 0.0
          enddo
          enddo
        endif
cbm
cbm   optionally initialize the octupole fields
cbm
        if (norun.eq.1) then
          rad03   = rad0 * rad02
          op      = one
          ao      = sqrt( -log( oh ) )/fzrm
          ao2     = ao * ao
          radpow  = 8.0
          call qgaus(fo1, zero, binup, oint)
          do      = op / oint
          fo0     = do * 4.0d0/7.0d0 *pi*chg
          pfy     = oh * oh * sqrt(7./pi)
          do ir=0,nr
          r = (ir+sr)*dr
          do iz=0,nz
             z     = (iz+sz-icm)*dz
             rad2  = z*z+r*r
             rad   = sqrt(rad2)
             rad3  = rad2 * rad
             rad7  = rad3 * rad3 * rad
             y30p  = pfy * (two*z*z*z - 3.*z*r*r)
             radpow = zero
             crhoo(iz,ir,0) = do*fo1(rad)*y30p
c   mind periodic boundaries:
             if (iz.eq.nz) crhoo(iz,ir,0) = zero
             if (iz.eq.0.and.isym.eq.0) crhoo(iz,ir,0) = zero
c
c            the  subroutine dglnq(f,x1,x2,o) computes the integral
c               int_{x1}^{x2} f(t) dt
c            with gauss legendre of order o.
c            the function f has to be provided as subroutine.
c
             if (rad.lt.1.d-2) then
              qt1 = zero
             else
              radpow = 8.0
              call qgaus(fo1, zero, rad, qt1)
              qt1 = qt1/rad7
             endif
             radpow = one
             call qgaus(fo1, rad , binup, qt2)
             coulo(iz,ir,0) = fo0*(qt1+qt2)*y30p
c   mind periodic boundaries:
             if (iz.eq.nz) coulo(iz,ir,0) = zero
             if (iz.eq.0.and.isym.eq.0) coulo(iz,ir,0) = zero
          enddo
          enddo
          call obs(crhoo,comono,rmso,codip,coquad,coct,cohex,qsqco
     &             ,zero)
          write(6,'(a/3x,7(1pg11.3))')
     &   ' from "crhoo": comono,rmso,codip,coquad,coct,cohex,qsqco=',
     &    comono,rmso,codip,coquad,coct,cohex,qsqco
        if(itest.ge.1) then
          call coldif(coulo,crhoo,wcolt,' octupole background')
c         call coldf2(coulo,crhoo,wcolt,' octupole background')
        endif
        else
          do ir=0,nr
          do iz=0,nz
            coulo(iz,ir,0) = 0.0
            crhoo(iz,ir,0) = 0.0
          enddo
          enddo
        endif
      endif
c
      return
      end
c-----fm1--------------------------------------------------------------
c
      function fm1(x)
      include 'paramr.inc'
c
      common /colinp/ am2,rad02,aq2,rad04,ah2,radpow,          ! for sym=1
     &                ad2,rad0,ao2,rad03                       ! for sym=0
c
      expa  = x*x*am2
      if (expa.le.eighty) then
         expa = exp(-expa)
      else
         expa = zero
      endif
      if(x.ne.zero) then
       fm1 =  x**radpow * expa
      elseif(radpow.eq.zero) then
       fm1 = expa
      else
       fm1 = zero
      endif
c
      return
      end
c-----fq1--------------------------------------------------------------
c
      function fq1(x)
      include 'paramr.inc'
c
      common /colinp/ am2,rad02,aq2,rad04,ah2,radpow,          ! for sym=1
     &                ad2,rad0,ao2,rad03                       ! for sym=0
c
      x2    = x*x
      expa  = x2*aq2
      if (expa.le.eighty) then
         expa = exp(-expa)
      else
         expa = zero
      endif
      if(radpow.eq.zero) then
       fq1 =  expa/(rad02+x2)
      elseif(x.ne.zero) then
       fq1 =  x**radpow*expa/(rad02+x2)
      else
       fq1 = zero
      endif
c
      return
      end
c-----fh1--------------------------------------------------------------
c
      function fh1(x)
      include 'paramr.inc'
c
      common /colinp/ am2,rad02,aq2,rad04,ah2,radpow,          ! for sym=1
     &                ad2,rad0,ao2,rad03                       ! for sym=0
c
      x2    = x*x
      expa  = x2*ah2
      if (expa.le.eighty) then
         expa = exp(-expa)
      else
         expa = zero
      endif
      if(radpow.eq.zero) then
       fh1   =    expa/(rad04+x2*x2)
      elseif(x.ne.zero) then
       fh1   =    x**radpow*expa/(rad04+x2*x2)
      else
       fh1   = zero
      endif
c
       return
       end
c-----fd1--------------------------------------------------------------
c
      function fd1(x)
      include 'paramr.inc'
c
      common /colinp/ am2,rad02,aq2,rad04,ah2,radpow,          ! for sym=1
     &                ad2,rad0,ao2,rad03                       ! for sym=0
c
c     x2    = x*x
      x0    = 5.
      x2    = (x-x0)*(x-x0)
      expa  = x2*ad2
      if (expa.le.eighty) then
         expa = exp(-expa)
      else
         expa = zero
      endif
      if(radpow.eq.zero) then
       fd1   =    expa/(rad0+x)
      elseif(x.ne.zero) then
       fd1   =    x**radpow*expa/(rad0+x)
      else
       fd1   = zero
      endif
c
       return
       end
c-----fo1--------------------------------------------------------------
c
      function fo1(x)
      include 'paramr.inc'
c
      common /colinp/ am2,rad02,aq2,rad04,ah2,radpow,          ! for sym=1
     &                ad2,rad0,ao2,rad03                       ! for sym=0
c
      x2    = x*x
      x3    = x2 * x
      expa  = x2*ao2
      if (expa.le.eighty) then
         expa = exp(-expa)
      else
         expa = zero
      endif
      if(radpow.eq.zero) then
       fo1   =    expa/(rad03 + x3)
      elseif(x.ne.zero) then
       fo1   =    x**radpow*expa/(rad03 + x3)
      else
       fo1   = zero
      endif
c
       return
       end
c----qgaus-------------------------------------------------------------
c
c     numerical recipes routine qgaus (20 point)
c
      subroutine qgaus(func,a,b,ss)
      include 'paramr.inc'
c
      dimension x(10),w(10)
c
      external func
c
      data x/.076526521133497, .227785851141645, .373706088715419,
     &       .510867001950827, .636053680726515, .746331906460150,
     &       .839116971822218, .912234428251325, .963971927277913,
     &       .993128599185094/
      data w/.152753387130725, .149172986472603, .142096109318382,
     &       .131688638449176, .118194531961518, .101930119817240,
     &       .083276741576704, .062672048334109, .040601429800386,
     &       .017614007139152/
c
      xm = 0.5d0*(b+a)
      xr = 0.5d0*(b-a)
      sum= 0.0d0
      do j=1,10
        dx = xr*x(j)
        sum= sum+w(j) * (func(xm+dx) + func(xm-dx))
      enddo
      ss = xr*sum
      return
      end

c-----colin3----------------------------------------------------------
c
      subroutine colin3(rhoini,wcol,mangmo)
c
c     initialises the monpole and quadrupole background fields for
c     the subroutine colrun, here using solution by radial integration
c     of the green's function
c
c      'mangmo' is the angular momentum about the z-axis
c               mangmo=0 is standard for ground state calculations
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
      dimension  wcol(0:nzd,0:nrd),rhoini(0:nzd,0:nrd)
c
c     define super-array which contains all reference fields
c     from /coulom/ adressable by pointers:
c       l=0 density    <--> lpoint = 0
c       l=0 potential  <--> lpoint = 1
c       l=1 density    <--> lpoint = 2
c       .
c       .
c       .
c       l=4 potential  <--> lpoint = 9
c     make sure that the fields in common /coulom/ are arranged in
c     the proper order!
c
      dimension colref(0:nzd,0:nrd,0:KANGMO,0:9)
      equivalence (colref(0,0,0,0),crhom(0,0,0))
c
c     parameters of the reference densities for the long range part,
c     provided her up l=4:
c       aexp   = width parameter for the Gaussian
c       bdiv   = weight parameter for the denominator
c       radpow = power of prefactor for integration, transfer to function call
c       lpow   = power of actual 'l', transfer to function call
c
      common /colinp3/ aexp(0:4),bdiv(0:4),rshif(0:4),radpow,lpow
c
c      rs  = Wigner-Seitz radius for nuclei in units [fm]
c
      parameter (rs=1.16)
c
c     switches to required cases  l=0..4
c
      logical  tcase(0:4)
c
      external refrho            ! function computing the reference density
c
c----------------------------------------------------------------------
c
c     initialisation of wcol
c
      do ir=0,nr
      do iz=0,nz
        wcol(iz,ir) = rhoini(iz,ir)
        wcol(iz,ir) = 0.0
      enddo
      enddo
c
c     general parameters for background fields
c
      zmax    = ((nz)+sz)*dz
      rmax    = ((nr)+sr)*dr
      zrmin   = min(zmax,rmax)
      zrmax   = max(zmax,rmax)
      binup   = zrmax+5.0d0
      fzrm    = zrmin/5.5d0
      rad0    = rs*anumb**edr
c
      agen    = sqrt( -log( oh ) )/fzrm
      agen    = one/min(rad0*oh,one/agen)
      agen    = agen*agen
      do l=0,4
        aexp(l) = agen
       if(l.ge.1) then
          bdiv(l)  = zrmax**l
          rshif(l) = rad0
        else
          rshif(0) = zero
        endif
      enddo
c
      tcase(0) = .true.
      tcase(1) = ndrun.eq.1 .and. isym.ne.1
      tcase(2) = nqrun.eq.1
      tcase(3) = isym.eq.0  .and. norun.eq.1
      tcase(4) = nhrun.eq.1
c
c
c     initialize the multipole fields
c
      do 20 l=0,4
      lpow  = l
      lpoint = l+l
      if(tcase(l)) then
c                                      normalization
        radpow = 2*l+2
        call qgaus(refrho, zero, binup, xnorm)
        rnorm   = one/xnorm
        vnorm   = rnorm*4.0*pi*chg/(l+l+one)
c
        do ir=0,nr
          r = (ir+sr)*dr
          do iz=0,nz
            z     = (iz+sz)*dz
            rad2  = z*z+r*r
            rad   = sqrt(rad2)
            ylmp  = ylm(z,r,l,mangmo)
c                                      density
            radpow = zero
            colref(iz,ir,mangmo,lpoint) = rnorm*refrho(rad)*ylmp
c                                      potential
            if (rad.lt.1.d-2) then
             t1 = zero
            else
             radpow = 2*l+2
             call qgaus(refrho, zero, rad, t1)
             t1 = t1/rad**(l+l+one)
            endif
            radpow = one
            call qgaus(refrho, rad , binup, t2)
            colref(iz,ir,mangmo,lpoint+1) = vnorm*(t1+t2)*ylmp
          enddo
        enddo
c
          call obs3(colref(0,0,mangmo,lpoint),cqmono,rmsq,cqdip,
     &              cqquad,cqoct,cqhex,betacq,mangmo)
          write(6,'(a,i2,a/5x,8(1pg11.3))')
     &    ' from "l=',l,
     &    '": cqmono,cqdip,cqquad,cqoct,cqhex,rmsq,betacq=',
     &        cqmono,cqdip,cqquad,cqoct,cqhex,rmsq,betacq
          if(l.eq.0) then
            cmono = cqmono
          elseif(l.eq.1) then
            cdip  = cqdip
          elseif(l.eq.2) then
            cquad = cqquad
          elseif(l.eq.3) then
            coct  = cqoct
          elseif(l.eq.4) then
            chex  = cqhex
          endif
c          if(itest.ge.1)
c     &      call coldif(colref(0,0,mangmo,lpoint+1),
c     &                  colref(0,0,mangmo,lpoint),
c     &                  wcolt,mangmo,' multipole back.')
          if(itest.ge.1) then
            call testpr(colref(0,0,mangmo,lpoint+1),0,nz,0,9,' Coul.')
            call testpr(colref(0,0,mangmo,lpoint),0,nz,0,9,' Dens.')
          endif
      else
        do ir=0,nr
        do iz=0,nz
          colref(iz,ir,mangmo,lpoint+1) = 0.0
          colref(iz,ir,mangmo,lpoint)   = 0.0
        enddo
        enddo
      endif
   20 continue
c
      return
      end

c-----refrho-----------------------------------------------------------
c
      function refrho(x)
c
      include 'paramr.inc'
c
c     reference density at 'x'
c     for given angular momentum 'lpow' from common /colinp/
c     and with prefactor 'radpow' from common /colinp/
c
      common /colinp3/ aexp(0:4),bdiv(0:4),rshif(0:4),radpow,lpow
c
      x2    = x-rshif(lpow)
      x2    = x2*x2
      arg   = x2*aexp(lpow)
      if (arg.le.eighty) then
         expa = exp(-arg)
      else
         expa = zero
      endif
      if(lpow.gt.0) then
        if(radpow.eq.zero) then
          refrho =  expa/(bdiv(lpow)+x**lpow)
        elseif(x.ne.zero) then
          refrho =  x**radpow*expa/(bdiv(lpow)+x**lpow)
        else
          refrho = zero
        endif
      else
        expa   = expa*(one+arg+oh*arg*arg*(one+edr*arg))
        if(x.ne.zero) then
          refrho =  x**radpow * expa
        elseif(radpow.eq.zero) then
          refrho = expa
        else
          refrho = zero
        endif
      endif
c
      return
      end
c-----ylm-------------------------------------------------------------
c
      function ylm(z,r,l,m)
c
      include 'paramr.inc'
c
c     spherical harmonics: R**l*Y_lm(r,z), R=sqrt(r*r+z*z)
c
c
c     prefactors for spherical harmonics
c
      data y00fac/0.282094791/     !     1/(2*SQRT(PI))
      data y10fac/0.488602512/     !   -SQRT(3)/(2*SQRT(PI))
      data y11fac/0.345494149/     !    SQRT(6)/(4*SQRT(PI))
      data y20fac/0.315391565/     !    SQRT(5)/(4*SQRT(PI))
      data y21fac/0.772548404/     !   -SQRT(30)/(4*SQRT(PI))
      data y22fac/0.386274200/     !
      data y30fac/0.746352665/     !   -2*SQRT(7)/(4*SQRT(PI))
      data y31fac/0.323180184/     !    SQRT(21)/(8*SQRT(PI))
      data y32fac/1.0        /     !    
      data y33fac/1.0        /     !    
      data y40fac/0.846284375/     !    24/(16*SQRT(PI))
      data y41fac/1.892349392/     !   -12*SQRT(5)/(8*SQRT(PI))
      data y42fac/1.0        /     !    
      data y43fac/1.0        /     !    
      data y44fac/1.0        /     !    

      double precision r,z
      integer          l,m
      
c
c----------------------------------------------------------------------
c
c     check forbidden cases
c
      if (l.lt.0) then
        stop ' l<0 is nonsense in FUNCTION YLM '
      endif
c
      if (m.gt.l) then
        stop ' m>l is nonsense in FUNCTION YLM '
      endif
c
      if (l.gt.4) then
        stop ' l>4 not yet provided in FUNCTION YLM '
      endif
c
c     determin phasefactor for negative m
c     
      if (m.lt.0) then
         phasefac = (-1)**m
         m        = -m
      else
         phasefac = +1
      endif
c
c     spherical harmonics in separate cases
c
      if (l.eq.0) then
            ylm = y00fac
c
      elseif (l.eq.1) then
         if (m.eq.0) then
            ylm = y10fac*z
         else
            ylm = y11fac*r
         endif
c     
      elseif (l.eq.2) then
         if (m.eq.0) then
            ylm = y20fac*(2.0D0*(z*z)-(r*r))
         elseif (m.eq.1) then
            ylm = y21fac*z*r
         else
            ylm = y22fac*r*r
         endif
c     
      elseif(l.eq.3) then
         if (m.eq.0) then
            ylm = y30fac*z*((z*z)-1.5D0*(r*r))
         elseif (m.eq.1) then
            ylm = y31fac*r*(4.0D0*(z*z)-(r*r))
         elseif (m.eq.2) then
            ylm = y32fac*(r*r)*z
         else
            ylm = y33fac*(r*r)*r
         endif
c     
      else
         if (m.eq.0) then
            ylm = y40fac*((z*z)*(z*z)
     &                    -3.0D0*(z*z)*(r*r)+0.375D0*(r*r)*(r*r))
         elseif (m.eq.1) then
            ylm = y41fac*z*r*((z*z)-0.75D0*(r*r))
         elseif (m.eq.2) then
            ylm = y42fac*1.0D0 ! <-----------------
         elseif (m.eq.3) then
            ylm = y43fac*1.0D0 ! <-----------------
         else
            ylm = y44fac*1.0D0 ! <-----------------
         endif
c     
      endif

      ylm = ylm * phasefac
c
      return
      end

c
c-----obs--------------------------------------------------------------
c
      subroutine obs3(rhoini,atnum,rms,d,q,o,h,qsq,mangmo)
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
      dimension rhoini(0:nzd,0:nrd)
c
      data three,seven,eight,twfour /3.0d0,7.0,8.0,24.0/
c
c----------------------------------------------------------------------
c
      if(mangmo.ne.0)
     &   stop 'OBS3 not yet prepared for MANGMO.NE.0'
c
c     summing up the wavefunction
c
      drezw  = one/12.0
      zcm    = 0.0D0            ! set default for Coulomb solver
c
cbm   at1    = oh*rhoini(0,0)
      at1    = rhoini(0,0)
      if (isym.eq.1) at1 = at1 * oh
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
         d3    = rhot*z   + d3
         q3    = rhot*qt  + q3
         qsq3  = rhot*q2t  + qsq3
         o3    = rhot*ot  + o3
         h3    = rhot*ht  + h3
       enddo
       rhot  = rhoini(nz,ir)*dir
       d3    = d3 - rhot*z
       o3    =  o3 - rhot*ot
      enddo
c
c     the different parts are summed up:
c
      atnum   =  at1+at2+at3
      rms     =  rms1 +rms2 +rms3
      d       =  d1 +d2 +d3
      q       =  q1 +q2 +q3
      qsq     =  qsq1 +qsq2 +qsq3
      o       =  o1 +o2 +o3
      h       =  h1 +h2 +h3
c
c     the final result
c
      rms   = sqrt(abs(rms/max(atnum,smal)))
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


ctstc-------coltst---------------------------------------------------------
ctstc
ctst      subroutine coltst(rhoin,wcol,idef)
ctstc
ctstc     test routine for coulomb solver
ctstc
ctst      include 'paramr.inc'
ctst      include 'levels.inc'
ctst      include 'fields.inc'
ctst      include 'dynpar.inc'
ctstc
ctst      dimension  wcol(0:nzd,0:nrd),rhoin(0:nzd,0:nrd)
ctstc
ctstc     pfacd = (3/(2*pi))**(1.5), rp: ca. rms radius of protons
ctstc
ctstcold      data pfacd,opf,tpf /0.32992261,1.5,2.5/
ctstc
ctstc
ctstc     monopole field and density
ctstc
ctst      qdef    = 40.0
ctst      pro     = (npart(1))
ctst      zmax    = ((nz-icm)+sz)*dz
ctst      rmax    = ((nr)+sr)*dr
ctst      zrmin   = min(zmax,rmax)
ctst      zrminh  = zrmin/6.5
ctstc
ctst      dis     = zrmin*0.06
ctst      xnq     = qdef/(two*two*dis*dis)
ctst      xaq     = sqrt(-log(oh))/zrminh
ctst      xaq2    = xaq  * xaq
ctst      xaq3    = xaq2 * xaq
ctst      d0q     = xaq3*xnq/pi**1.5
ctstcold      pfacq   = two*pi/xaq2
ctstc
ctst      am      = sqrt(-log(oh))/zrminh
ctst      am2     = am*am
ctst      am3     = am2*am
ctst      d0m     = am3*pro/pi**1.5
ctst      pfac2   = two*pi*d0m/am2
ctst      do 20 ir=0,nr
ctst      r = (ir+sr)*dr
ctst      r2= r*r
ctst      do 20 iz=0,nz
ctst         z     = (iz+sz-icm)*dz
ctst         z2    = z*z
ctst         rad2  = z2+r2
ctst         rad   = sqrt(rad2)
ctstc
ctst         expa  = rad2*am2
ctst         if (expa.le.eighty) then
ctst            expf = exp(-expa)
ctst         else
ctst            expf = zero
ctst         endif
ctst         dichm        = d0m*expf
ctst         radam        = rad*am
ctstc
ctstc        the function subroutine errf2 provides the integral
ctstc          errf2(x) = 2/gamma(3/2) * int_{0}^{x} exp(-t*t) t**2 dt
ctstc
ctst         if (rad.lt.1.d-2) then
ctst          qt1 = zero
ctst         else
ctst          qt1 = errf2(radam)
ctst          qt1 = qt1*pro/rad
ctst         endif
ctst         qt2 = pfac2 * expf
ctst         feldm        = chg*(qt1+qt2)
ctstc
ctstc        quadrupol:
ctstc
ctst       zpd    = z + dis
ctst       zpd2   = zpd*zpd
ctst       zmd    = z - dis
ctst       zmd2   = zmd*zmd
ctst       radpd2 = r2+zpd2
ctst       radmd2 = r2+zmd2
ctstc       radpd  = sqrt(radpd2)
ctstc       radmd  = sqrt(radmd2)
ctst       expaq  = xaq2*rad2
ctst       expapd = xaq2*radpd2
ctst       expamd = xaq2*radmd2
ctst       if (expaq.le.eighty) then
ctst         expfq  = exp(-expaq)
ctst       else
ctst         expfq  = zero
ctst       endif
ctst       if (expapd.le.eighty) then
ctst         expfpd = exp(-expapd)
ctst       else
ctst         expfpd = zero
ctst       endif
ctst       if (expamd.le.eighty) then
ctst         expfmd = exp(-expamd)
ctst       else
ctst         expfmd = zero
ctst       endif
ctst       dens   = d0q*(expfmd-2.0*expfq+expfpd)
ctstc
ctst       dichq        = dens
ctst       if (idef.eq.0) then
ctst        rhoin(iz,ir) = dichm
ctst        wcol(iz,ir)= feldm
ctst       else if (idef.eq.1) then
ctst        rhoin(iz,ir) = dichm+dichq
ctst        wcol(iz,ir)= feldm+dichq
ctst       endif
ctst20    continue
ctstc     do 10 ir=0,nr
ctstc     do 10 iz=0,nz
ctstc       wcol(iz,ir) = coulm(iz,ir)
ctstc10   continue
ctstc
ctstc     call testpr(coulm,0,9,0,9,' coulm')
ctstc     call testpr(crhom,0,9,0,9,' crhom')
ctstc
ctst      return
ctst      end
ctstc-----errf2 -----------------------------------------------------------
ctstc
ctst      function errf2(x)
ctstc
ctstc        the function subroutine errf2 provides the integral
ctstc          errf2(x) = 2/gamma(3/2) * int_{0}^{x} exp(-t*t) t**2 dt;
ctstc        uses reduction to error function and tschebyscheff approximation.
ctstc        (for the latter see 'numerical recipes', chapter 6.2).
ctstc
ctstc
ctst      include 'paramr.inc'
ctstc
ctstc     compute error function for x with fractional error less 1.2e-7.
ctstc
ctst      z     = abs(x)
ctst      t     = 2.0/(2.0+z)
ctst      erfc  = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
ctst     &       t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
ctst     &       t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
ctst      if(x.lt.0) erfc = 2.0-erfc
ctst      erf   = 1.0-erfc
ctstc
ctstc     now compose the final function
ctstc
ctst      t     = 1.772453851*exp(x*x)
ctst      errf2 = -(x+x-t*erf)/t
ctstc
ctst      return
ctst      end


c=======================================================================
c
c     ini_coulomb(mangmol,mangmom)
c
c     initialize the coulomb solver (background fields and densities,
c     damping mass.
c
c       mangmol :  maximum angular momentum l of background fields
c       mangmom :  maximum z-component m of angular momentum
c
c=======================================================================
      subroutine ini_coulomb(mangmol,mangmom)
      
      include 'paramr.inc'
      include 'constr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'

      integer     mangmol
      integer     mangmom
      integer     ir, iz

      data ddw/0.4/          ! width of damping
      data ddmass/1.0/       ! factor on damping mass
      data d2/2.0/           ! thickness of surface zone
      data dqr/0.0/          ! to make this variable static
c
c     ...
c-----------------------------------------------------------------------
      zmax    = ((nz-icm)+sz)*dz
      rmax    = ((nr)+sr)*dr
      drmin = min(zmax,rmax) - d2
      drmin2 = drmin * drmin
      dzadd = zmax - drmin - d2
      dradd = rmax - drmin - d2
      dtermr = drmin2 - dradd*dradd
      dtermz = drmin2 - dzadd*dzadd
c
c     initialize damping mass
c-----------------------------------------------------------------------
      do iz=0,nz
        z  = (iz+sz-icm)*dz
        z2 = z * z
c
        if( dzadd.ge.dradd ) then
          if (z.ge.dzadd) then
            dqr = sqrt(dtermz + d2*dzadd*z)
          else
            dqr = sqrt(z2 + drmin2)
          endif
        endif
c
        do ir=0,nr
          r  = (ir+sr)*dr
          r2 = r * r
          drr = dsqrt(z2 + r2)
          if( dradd.gt.dzadd ) then
            if (r.ge.dradd) then
              dqr = sqrt(dtermr + d2*dradd*r)
            else
              dqr = sqrt(r2 + drmin2)
            endif
          endif
          sqrm(iz,ir) = ddmass / (one + exp((dqr-drr)/ddw))
        enddo
      enddo
c
      if(itest.gt.0) call testpr(sqrm,0,nz,0,nr,'sqrm')     
c
c     initialize background fields
c-----------------------------------------------------------------------

c
c     initialize background densities
c-----------------------------------------------------------------------



      
      return
      end
