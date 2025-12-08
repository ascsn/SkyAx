c-----yukrun------------------------------------------------------------
c
      subroutine yukrun(rhoin,wyuk,iter,width)
c
c     this subroutine solves the Laplace equation for mass 'ymass'.
c       rhoin    = density which is to be folded
c       wyuk     = Yukawa folded density
c       iter     = number of outer iteration
c       width    = spatial width of Yukawa folding
c
      include 'paramr.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'observ.inc'
c
      dimension wyuk(0:nzd,0:nrd),rhoin(0:nzd,0:nrd),ywork(0:nzd,0:nrd)
c
      data ephimx/0.0/     ! to make it static
      LOGICAL :: testyuk=.FALSE.,tfirst=.TRUE.
c
      if (iter.le.10) then
        irepmx = 10*irepc
      else
        irepmx = 4*irepc            ! use the same iteration number as Coulomb
      endif
      dfac    = 0.1D0 !one/cmssqr
      ymass2 = 6D0/width**2  
      IF(testyuk) write(6,'(a,4g11.3)')
     &     'yukrun: cmssqr,dfac=',cmssqr,dfac,ymass2,width
c
      ephimx=0D0
      IF(tfirst) wyuk(0:nz,0:nr) = rhoin(0:nz,0:nr)/ymass2
!      tfirst = .FALSE.
      do 40 irep=1,irepmx
!      do 40 irep=1,0
c     if (mod(irep,20).eq.0) then
c     write(6,'(a,i3)') ' iteration nr=',irep-1
c     call testpr(wcol,0,9,0,9,' wcol  ')
c     call testpr(wcol,nz-9,nz,nr-9,nr,' wcol  ')
c     endif
c
      call laplace_scalar(wyuk,ywork,1,0,1,-one)  
c
      do ir=0,nr
      do iz=0,nz
         ywork(iz,ir) = ywork(iz,ir)+ymass2*wyuk(iz,ir) - rhoin(iz,ir)
      enddo
      enddo
c
      call damp_meson(ywork,ywork,one)
cPK   call mulz(zdmes(0,0,1),ywork,wrk1,one,0,-1)
cPK   call mulr(rdmes(0,0,1),wrk1,ywork,one,0,-1)
c
      ephimx = -one
      do ir=0,nr
      do iz=0,nz
       ephi    = ywork(iz,ir)
       if (abs(ephi).gt.ephimx) then
         ephimx = abs(ephi)
       endif
       wyuk(iz,ir) = wyuk(iz,ir) - dfac*ephi
      enddo
      enddo

c     if (mod(irep,1).eq.0) then
      IF(testyuk) write(6,'(/a,g15.5,3(a,i3)/)') 
     &    ' max error=',ephimx,' rep=',irep
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
c
      anormin=f1int(rhoin)
      anormout=f1int(wyuk)
      wyuk(0:nz,0:nr) = wyuk(0:nz,0:nr)*(anormin/anormout)
      anormout=f1int(wyuk)
      if (abs(ephimx).ge.colerr .OR. testyuk) then
c       write(6,*)colerr,abs(ephimx),irepmx
       IF(testyuk) THEN
         write(6,'(/2(a,g15.5),i4,a/)') ' yukrun: max error=',
     &                   colerr,' actual error=',abs(ephimx),
     &                   irep,' iterations'
        call obs(rho(0,0,1),apnum,rmsp,dp,qp,op,hp,qsqp,zcm)
        call obs(wyuk(0,0),aynum,rmsy,dp,qp,op,hp,qsqp,zcm)
        WRITE(6,'(a,i5,5(1pg13.5))') 'Yukawa: iter,a,r=',
     &      iter,apnum,aynum,rmsp,rmsy,anormout
       END IF
      endif
c      if(itest.ge.1)
c     &  call coldif(wcol,rhoin,ywork,' full solution')
c
      return
      end
