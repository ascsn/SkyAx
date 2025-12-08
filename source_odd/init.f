c=======================================================================
c
c     block #2: initialization of the derivative matrices
c
c=======================================================================
c-----inider------------------------------------------------------------
c
      subroutine inider
c
c     the different derivative matrices in the ft basis are
c     constructed.
c     the different cases are explained in the include member
c     'deriv'.
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'wavefr.inc'
c
      call procos(rpro,nrd,nr,icutr) !!  needs parity ??
      call procos(zpro,nzd,nz,icutz) !!  needs parity ??
c
      if (idrtyp.eq.0) then
         call ini_grid(-1,-1)
      elseif (idrtyp.eq.1) then
         call ini_grid(9,-1)
      elseif (idrtyp.eq.2) then
         call ini_grid(-1,9)
      elseif (idrtyp.eq.3) then
         call ini_grid(9,9)
      endif
c
      return
      end
c-----sder  -----------------------------------------------------------
c
      subroutine sder(der,ipar,nmax,n,dx,icut)
c
c     1. derivative in r-direction using fourier-definition
c
      include 'paramr.inc'
      dimension der(0:nmax,0:nmax)
c
      dn    = (n)
      afac  = pi/dn
      pfac  = 2.0*afac/(dn*dx)
cbm
      if (ipar.eq.0) then
       iodd = mod(n+1,2)
       if (iodd.ne.0) then
        write(7,'(a)') ' warning: nz has to be odd !'
        stop
       endif
       icn   = (n+1)/2
       dn    = (icn)
       afac  = pi / dn
       pfac  = (one/(dn))*(afac/dx)
c
c
       do l=0,n
         dl  = (l)
         do i=0,n
           di     = (i)
           sum    = zero
           do j=1,icn
             dj  = (j)
c            dj2  = dj*dj
             sum  = sum  - dj*sin(afac*dj*(dl-di))
           enddo
cpgr        sum = sum - oh * dfloat((-1)**(l-i)) * dn**2
           der(i,l) = - pfac * sum
         enddo
       enddo
cbm
c
      elseif (ipar.eq.1) then
c
       do l=0,n
         dl  = (l)
         do i=0,n
           di     = (i)
           sum    = zero
           do j=1,n-icut
             dj  = (j)
             sum  = sum  + dj*sin(afac*dj*di)*cos(afac*dj*dl)
           enddo
           if (l.eq.0.or.l.eq.n) then
             der(i,l) = -0.5*pfac * sum
           else
             der(i,l) = -      pfac * sum
           endif
         enddo
       enddo
c
      else if (ipar.eq.2) then
c
       do l=0,n
         dl  = (l)
         do i=0,n
           di     =  (i)
           sum    =  zero
           do j=0,n
             dj  = (j)
             sum = sum + dj*cos(afac*dj*di)*sin(afac*dj*dl)
           enddo
           der(i,l) = pfac * sum
         enddo
       enddo
c
      endif
c
      return
      end
c-----sder2-------------------------------------------------------------
c
      subroutine sder2(der,ipar,nmax,n,dx,icut)
c
c     2. derivative in r-or-z-direction using fourier-definition
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'observ.inc'
      dimension der(0:nmax,0:nmax)
c
cbm
      if (ipar.eq.0) then
       iodd = mod(n+1,2)
       if (iodd.ne.0) then
        write(7,'(a)') ' warning: nz has to be odd !'
        stop
       endif
       icn   = (n+1)/2
       dn    = (icn)
       afac  = pi/dn
       pfac  = (one/(dn))*(afac/dx)**2
c
c
       do l=0,n
        dl  = (l)
        do i=0,n
         di     = (i)
         sum    = zero
         do j=1,icn
          dj  = (j)
          dj2  = dj*dj
          sum  = sum  + dj2*cos(afac*dj*(dl-di))
         enddo
        sum = sum - oh * real((-1)**(l-i)) * dn**2
        der(i,l) = - pfac * sum
       enddo
      enddo
cbm
c
      else if (ipar.eq.1) then
c
      dn    = (n)
       afac  = pi/dn
      pfac  = (2.0/dn)*(afac/dx)**2
c
      do l=0,n
      dl  = (l)
      do i=0,n
       di     = (i)
       sum    = zero
       do j=1,n-icut
        dj  = (j)
        dj2  = dj*dj
        if(j.eq.n) then
c         dj2  = (4.0*(n-1)**2-(n-2)**2)/3.0
         sum  = sum  + 0.5*dj2*cos(afac*dj*di)*cos(afac*dj*dl)
        else
         sum  = sum  + dj2*cos(afac*dj*di)*cos(afac*dj*dl)
        endif
       enddo
       if (l.eq.0.or.l.eq.n) then
        der(i,l) = -0.5*pfac * sum
       else
        der(i,l) = -      pfac * sum
       endif
      enddo
      enddo
c
      else if (ipar.eq.2) then
c
      dn    = (n)
       afac  = pi/dn
      pfac  = (2.0/dn)*(afac/dx)**2
      do l=0,n
      dl  = (l)
      do i=0,n
       di     =  (i)
       sum    =  zero
       do j=0,n
        dj  = (j)
        dj2  = -dj*dj
        sum = sum + dj2*sin(afac*dj*dl)*sin(afac*dj*di)
       enddo
       der(i,l) = pfac * sum
      enddo
      enddo
c
      endif
c
cbm
c   test spectral properties
cbm
      if (ipar.ne.0) then
        nk   = n
      else
        nk   = (n+1)/2
      endif
c
      if(itest.ge.1) then
       if (ipar.eq.0.or.ipar.eq.1) then
       do 991 k=0,nk
         ptest  = afac*k
         sumnor = zero
         sumder = zero
         do 992 i=0,n
           cosi   = cos(i*ptest)
           sumnor = cosi*cosi + sumnor
           sumt   = zero
           do 993 j=0,n
             sumder  = sumder  + cos(j*ptest)*der(j,i)*cosi
             sumt    = sumt    + cos(j*ptest)*der(i,j)
  993      continue
           if(abs(cosi).ge.1.e-10) then
            wrk2(i,k) = sumt/cosi
           else
            wrk2(i,k) = zero
           endif
  992    continue
         wrk1(k,0) = sumder/sumnor
  991  continue
c
c     print spectral properties
c
        write(6,'(a,10(/1x,5g12.4))')
     &  ' test spectrum of 2. derivative:',
     &  (wrk1(k,0),k=0,nk)
        write(6,'(a,10(/1x,5g12.4))')
     &  ' reproduction of 1. state:',
     &  (wrk2(i,1),i=0,n)
        write(6,'(a,10(/1x,5g12.4))')
     &  ' reproduction of last state:',
     &  (wrk2(i,nk),i=0,n)
       endif
c     
       if (ipar.eq.0.or.ipar.eq.2) then
       do 994 k=1,nk
        ptest  = afac*k
        sumnor = zero
        sumder = zero
        do 995 i=0,n
          cosi   = sin(i*ptest)
          sumnor = cosi*cosi + sumnor
          sumt   = zero
          do 996 j=1,n
            sumder  = sumder  + sin(j*ptest)*der(j,i)*cosi
            sumt    = sumt    + sin(j*ptest)*der(i,j)
  996     continue
           if(abs(cosi).ge.1.e-10) then
            wrk2(i,k) = sumt/cosi
           else
            wrk2(i,k) = zero
           endif
  995     continue
        wrk1(k,0) = sumder/sumnor
  994 continue
c
c     print spectral properties
c
        write(6,'(a,10(/1x,5g12.4))')
     &  ' test spectrum of 2. derivative:',
     &  (wrk1(k,0),k=0,nk)
        write(6,'(a,10(/1x,5g12.4))')
     &  ' reproduction of 1. state:',
     &  (wrk2(i,1),i=0,n)
        write(6,'(a,10(/1x,5g12.4))')
     &  ' reproduction of last state:',
     &  (wrk2(i,nk-1),i=0,n)
       endif
c
      endif
c
c
      return
      end
c-----fder2-------------------------------------------------------------
c
      subroutine fder2(der,ipar,nmax,n,dx)
c
c     second derivative with 9 point formula
c
      include 'paramr.inc'
      include 'wavefr.inc'
      dimension der(0:nmax,0:nmax)
      parameter (korder=9, kdim=korder/2)
c5      parameter (korder=5, kdim=korder/2)
      dimension points(-kdim:kdim)
c
c      data zero,pi /0.0,3.141592653589793/
c
c     preset 9 point formula
c
      dx2  = 1.0/(dx*dx)
      points(-4) = -9/5040.*dx2
      points(-3) = 128/5040.*dx2
      points(-2) = -1008/5040.*dx2
      points(-1) = 8064/5040.*dx2
      points( 0) = -14350/5040.*dx2
      points( 1) = 8064/5040.*dx2
      points( 2) = -1008/5040.*dx2
      points( 3) = 128/5040.*dx2
      points( 4) = -9/5040.*dx2
c5      points(-2) = -1/12.*dx2
c5      points(-1) = 16/12.*dx2
c5      points( 0) = -30/12.*dx2
c5      points( 1) = 16/12.*dx2
c5      points( 2) = -1/12.*dx2
c
c     preset zero
c
      do j=0,n
      do i=0,n
      der(i,j) = zero
      enddo
      enddo
c
c     now fill the nonvanishing elements.
c     parity is used near zero
c     and at the outer boundary.
c
      if(ipar.eq.1) then
       parity = 1.0
      else
       parity = -1.0
      endif
c
c     straightforward mid-region
c
      do i=0,n-kdim
cfull      do i=0,n
       do j=i-kdim,i+kdim
       if(j.le.n) then
        if(j.lt.0) then
         der(i,abs(j)) = parity*points(j-i) + der(i,abs(j))
        else
         der(i,j) = points(j-i) + der(i,j)
        endif
       else
        jeff  = j-n
        jeff  = n-jeff
        der(i,jeff) = parity*points(j-i) + der(i,jeff)
       endif
ctest       write(7,'(a,3i5,f9.4)') ' kdim,i,j,der:',kdim,i,j,der(i,j)
       enddo
      enddo
c
c     lower precision at the outer boundary
c
c     5-point precision
c
      der(n-3,n-5) = -1/12.*dx2
      der(n-3,n-4) = 16/12.*dx2
      der(n-3,n-3) = -30/12.*dx2
      der(n-3,n-2) = 16/12.*dx2
      der(n-3,n-1) = -1/12.*dx2
c
      der(n-2,n-4) = -1/12.*dx2
      der(n-2,n-3) = 16/12.*dx2
      der(n-2,n-2) = -30/12.*dx2
      der(n-2,n-1) = 16/12.*dx2
      der(n-2,n)   = -1/12.*dx2
c
c     3-point precision
c
      der(n-1,n-2) =  dx2
      der(n-1,n-1) =  -2*dx2
      der(n-1,n)   =  dx2
      der(n,n-1)   =  dx2
      der(n,n)     =  -2*dx2
c
c
c     test spectral properties
c
ctest      afac  = pi/n
ctest      do k=0,n
ctest        ptest  = afac*k
ctest        sumnor = zero
ctest        sumder = zero
ctest        do i=0,n
ctest          cosi   = cos(i*ptest)
ctest          sumnor = cosi*cosi + sumnor
ctest          sumt   = zero
ctest          do j=0,n
ctest            sumder  = sumder  + cos(j*ptest)*der(j,i)*cosi
ctest            sumt    = sumt    + cos(j*ptest)*der(i,j)
ctest          enddo
ctest          if(cosi.ne.zero) wrk2(i,k) = sumt/cosi
ctest        enddo
ctest        wrk1(k,0) = sumder/sumnor
ctest      enddo
ctest      write(7,'(a,10(/1x,5g12.4))')
ctest     &  ' test spectrum od 2. derivative:',
ctest     &  (wrk1(k,0),k=0,n)
ctest      write(7,'(a,10(/1x,5g12.4))')
ctest     &  ' reproduction of 1. state:',
ctest     &  (wrk2(i,1),i=0,n)
ctest      write(7,'(a,10(/1x,5g12.4))')
ctest     &  ' reproduction of last state:',
ctest     &  (wrk2(i,n),i=0,n)
ctest      write(7,'(a)')
ctest     &  ' print der(i,j):'
ctest      do j=0,n
ctest        write(7,'(10f8.4)') (der(i,j),i=0,n)
ctest      enddo
c
c
      return
      end
c-----fder-------------------------------------------------------------
c
      subroutine fder(der,ipar,nmax,n,dx)
c
c     first derivative with 7-point formula
c
      include 'paramr.inc'
      dimension der(0:nmax,0:nmax)
c
      parameter (korder=11, kdim=korder/2)
      dimension points(-kdim:kdim)
c      data for 7-point formula 'korder=7'
c      data points/-0.01666666666667,0.15,-0.75,0.0,0.75,-0.15,
c     &            0.01666666666667/
c
c
c     set the field 'points' for the 11-point formula
c
      points(  -5) =  -1/(1260*dx)
      points(  -4) =  5/(504*dx)
      points(  -3) =  -5/(84*dx)
      points(  -2) =  5/(21*dx)
      points(  -1) =  -5/(6*dx)
      points(   0) =  0.0
      points(   1) =  5/(6*dx)
      points(   2) =  -5/(21*dx)
      points(   3) =  5/(84*dx)
      points(   4) =  -5/(504*dx)
      points(   5) =  1/(1260*dx)
c
c     preset zero
c
      do j=0,n
      do i=0,n
      der(i,j) = zero
      enddo
      enddo
c
c     now fill the nonvanishing elements.
c     parity is used near zero.
c
      if(ipar.eq.1) then
       parity = 1.0
      else
       parity = -1.0
      endif
c
c     straightforward mid-region
c
      do i=0,n-kdim
       do j=i-kdim,i+kdim
       if(j.le.n) then
        if(j.lt.0) then
         der(i,abs(j)) = parity*points(j-i) + der(i,abs(j))
        else
         der(i,j) = points(j-i) + der(i,j)
        endif
       endif
       enddo
      enddo
c
c     lower precision at the outer boundary
c
c     9-point precision
c
      der(n-4,n-8) =  1/(280*dx)
      der(n-4,n-7) =  -4/(105*dx)
      der(n-4,n-6) =  1/(5*dx)
      der(n-4,n-5) =  -4/(5*dx)
      der(n-4,n-4) =  0.
      der(n-4,n-3) =  4/(5*dx)
      der(n-4,n-2) =  -1/(5*dx)
      der(n-4,n-1) =  4/(105*dx)
      der(n-4,n)   =  -1/(280*dx)
c
c     7-point precision
c
      der(n-3,n-6) =  -1/(60*dx)
      der(n-3,n-5) =  3/(20*dx)
      der(n-3,n-4) =  -3/(4*dx)
      der(n-3,n-3) =  0.0
      der(n-3,n-2) =  3/(4*dx)
      der(n-3,n-1) =  -3/(20*dx)
      der(n-3,n)   =  1/(60*dx)
c
c     5-point precision
c
      der(n-2,n-4) =  1/(12*dx)
      der(n-2,n-3) =  -2/(3*dx)
      der(n-2,n-2) =  0.0
      der(n-2,n-1) =  2/(3*dx)
      der(n-2,n)   =  -1/(12*dx)
c
c     3-point precision
c
      der(n-1,n-2) =  -1/2/dx
      der(n-1,n-1) =  0.0
      der(n-1,n)   =  1/2/dx
      der(n,n-1)   =  0.0
      der(n,n)     =  1/2/dx
c
      return
      end
c-----procos-----------------------------------------------------------
c
      subroutine procos(p,nmax,n,icut)
c
c     construction of projection operator. icut: number of cut off
c     components in momentum space. icut=0 => unity matrix.
c
      include 'paramr.inc'
      dimension p(0:nmax,0:nmax)
c
c
      dn    = (n)
      afac  = pi/dn
      pfac  = 2.0/dn
      do 50 j=0,n
      dj  = (j)
      do 50 m=0,n
       dm     =  (m)
       sum    =  zero
cold       do 40 i=0,n-icut
       do 40 i=n-icut+1,n
        di  = (i)
        if (i.eq.0.or.i.eq.n) then
         sum = sum + 0.5*cos(afac*dj*di)*cos(afac*dm*di)
        else
         sum = sum +       cos(afac*dj*di)*cos(afac*dm*di)
        endif
 40    continue
       if(m.eq.j) then
        p(j,j) = 1.0
       else
        p(j,m) = 0.0
       endif
       if (m.eq.0.or.m.eq.n) then
        p(j,m) = p(j,m)-0.5*pfac * sum
       else
        p(j,m) = p(j,m)    -pfac * sum
       endif
 50   continue
c
c      write(7,'(a)') ' in procos:'
c      do 160 k=0,n
c      write(7,'(a,i3/100(5f12.5/))')
c     &    'k=',k,(p(k,j),j=0,n)
c 160  continue
c
      return
      end

c-----prifou-----------------------------------------------------------
c
      subroutine prifou(field,text)
c
c     print fourier components of 'field' along z- and -r-axis
c
      include 'paramr.inc'
c
      include 'dynpar.inc'
      include 'wavefr.inc'
c
      character*(*) text
      dimension field(0:nzd,0:nrd)
c
c      data zero,one /0.0,1.0/
c      data pi /3.141592653589793/
c
c     r-direction
c
      dn    = (nr)
      afac  = pi/dn
c      pfac  = 2.0/dn
      do kr=0,nr
        cfac = kr*afac
        sum  = 0.0
        do ir=0,nr
          cosfac = cos(cfac*ir)
          if(ir.eq.0 .or. ir.eq.nr) then
            cosfac = 0.5*cosfac
          endif
          sum  = 0.5*cosfac*field(0,ir)+sum
          do iz=1,nz-1
           sum  = cosfac*field(iz,ir)+sum
          enddo
          sum  = 0.5*cosfac*field(nz,ir)+sum
        enddo
        wrk1(0,kr) = sum
      enddo
      write(7,'(3a,10(/1x,5g12.4))')
     &  'cosine transform of ',text,' in r-direction:',
     &  (wrk1(0,ir),ir=0,nr)
c
c     z-direction
c
      dn    = (nz)
      afac  = pi/dn
c      pfac  = 2.0/dn
      do kz=0,nz
        cfac = kz*afac
        sum  = 0.0
        do iz=0,nz
          cosfac = cos(cfac*iz)
          if(iz.eq.0 .or. iz.eq.nz) then
            cosfac = 0.5*cosfac
          endif
          sum  = 0.5*cosfac*field(iz,0)+sum
          do ir=1,nr-1
           sum  = cosfac*field(iz,ir)+sum
          enddo
          sum  = 0.5*cosfac*field(iz,nr)+sum
        enddo
        wrk1(kz,0) = sum
      enddo
      write(7,'(3a,10(/1x,5g12.4))')
     &  'cosine transform of ',text,' in z-direction:',
     &  (wrk1(iz,0),iz=0,nz)
c
      return
      end

c=======================================================================
c
c     block #3: several other initializations
c
c=======================================================================
c-----const -----------------------------------------------------------
c
      subroutine const
c
c     sets some constants for common /const/ in 'paramr.inc'
c     and in 'dynpar.inc' .
c     this is done dynamically here rather than employing
c     a block data subroutine
c
      include 'paramr.inc'
      include 'dynpar.inc'
c
c     begin of the data definitions
c
c      pi     = 3.141592653589793
      hc     = 197.3285879
      e2     = 0.00729735*hc
c      zero   = 0.0
c      edr    = 1.0/3.0
c      oh     = 0.5
c      one    = 1.0
c      two    = 2.0
c      eighty = 80.0
c      smal   = 1.0e-30
      r0     = 1.400
      b_strength = 0D0
c
c     initial file number for rho- and psi plot
c
      iou = 18
c
      return
      end
c=======================================================================
c-----param -----------------------------------------------------------
c
      subroutine param(set)
*
************************************************************************
*                                                                      *
*     in order to make input easier the used parameter sets are        *
*     listed and read from here                                        *
*                                                                      *
************************************************************************
*
      include 'paramr.inc'
c
      parameter (ksets=8)
c
      character*20 set,names(ksets)
      include 'dynpar.inc'
c
      dimension parin(12,ksets)
      data names/'skyrme M ',
     &           'skyrme M , no ls',
     &            'ext.osc. ',
     &           ' m*/m = 1.0 ',
     &           ' m*/m = 0.90 ',
     &           ' m*/m = 0.81 ',
     &           ' m*/m = 0.74 ',
     &           ' m*/m = 0.68 '/
      data parin/-2645.0,385.0,-120.0,15595.0,130.0,
     &           0.09,0.0,0.0,0.0,0.0,0.16666666667,20.75,
     &           -2645.0,385.0,-120.0,15595.0,0.0,
     &           0.09,0.0,0.0,0.0,0.0,0.16666666667,20.75,
     &  0.0,0.0,0.0,0.0,0.0,    0.0,0.0,0.0,0.0,0.0,0.0,20.75,
     &   -1600.735,   205.894,   -13.894, 11875.792,   112.351,
     &   -0.350069, -0.266870,  9.864514, -0.746603,       0.0,
     &   0.3860479,20.7525,
     &   -1663.007,   267.154,    28.846, 11388.487,   115.369,
     &   -0.337644, -0.297870, -6.116109, -0.772582,       0.0,
     &   0.3536788,20.7525,
     &   -1798.334,   330.052,    69.948, 11262.011,   118.871,
     &   -0.333179, -0.318210, -3.073331, -0.777484,       0.0,
     &   0.3000000,20.7525,
     &   -1846.054,   389.115,   110.885, 10809.466,   121.553,
     &   -0.323949, -0.332954, -2.258584, -0.799724,       0.0,
     &   0.2760962,20.7525,
     &   -2010.772,   450.169,   149.831, 10974.273,   124.609,
     &   -0.325938, -0.344503, -1.901577, -0.790546,       0.0,
     &   0.2271105,20.7525/
c
c     further parametrisations for use through free input.
c     Variation of "C3P":
c C3P=0 :
c     T0    ,    T1    ,     T2   ,     T3   ,     T4   ,
c  -1798.334,   330.052,    69.948, 11262.011,   118.871,
c     X0    ,    X1    ,     X2   ,     X3   ,    B4P,  POWER  ,
c  -0.333179, -0.318210, -3.073331, -0.777484,    0.0,0.3000000,
c C3P=-2000 :
c     T0    ,    T1    ,     T2   ,     T3   ,     T4   ,
c  -1800.013,   338.518,    61.482, 11256.033,   119.006,
c     X0    ,    X1    ,     X2   ,     X3   ,    B4P,  POWER  ,
c   0.222217, -0.322757, -3.427706,  0.147806,    0.0,0.3000000,
c C3P=-4000 :
c     T0    ,    T1    ,     T2   ,     T3   ,     T4   ,
c  -1800.934,   345.313,    54.687, 11246.145,   119.073,
c     X0    ,    X1    ,     X2   ,     X3   ,    B4P,  POWER  ,
c   0.777115, -0.326244, -3.791438,  1.074626,    0.0,0.3000000,
c C3P=-6000 :
c     T0    ,    T1    ,     T2   ,     T3   ,     T4   ,
c  -1800.290,   350.147,    49.853, 11223.874,   118.995,
c     X0    ,    X1    ,     X2   ,     X3   ,    B4P,  POWER  ,
c   1.333038, -0.328643, -4.110627,  2.005839,    0.0,0.3000000,
c SCANLARG.RES 2.set:
c     T0    ,    T1    ,     T2   ,     T3   ,     T4   ,
c  -1719.445,   442.017,  2229.465,  9506.987,   122.352,
c     X0    ,    X1    ,     X2   ,     X3   ,    B4P,  POWER  ,
c  -0.308600, -4.399423, -1.291047, -0.544043,    0.0, 0.30000,
c All parametrizations require CHGX=1 .
c
c
************************************************************************
*                                                                      *
*     retrieve the parameters from data fields.                        *
*                                                                      *
************************************************************************
      if(iparam.gt.8)
     &        stop ' this paramter set not yet implemented'
c
      if(iparam.gt.0) then
       set    = names(iparam)
       t0     = parin(1,iparam)
       t1     = parin(2,iparam)
       t2     = parin(3,iparam)
       t3     = parin(4,iparam)
       t4     = parin(5,iparam)
       x0     = parin(6,iparam)
       x1     = parin(7,iparam)
       x2     = parin(8,iparam)
       x3     = parin(9,iparam)
       b4p    = parin(10,iparam)
       power  = parin(11,iparam)
       h2m(1) = parin(12,iparam)
       h2m(2) = h2m(1)
       ifzpe = 1
      endif
c
c
      end
c
c-----inicw -----------------------------------------------------------
c
      subroutine inicw(nmass,nchg)
c
c     initialises integration weights, constraint terms
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'dynpar.inc'
      include 'inipar.inc'
      include 'observ.inc'
      include 'levels.inc'
c
      chg   = chg * e2        ! ?? changing input data is dangerous !
c
c     axial symmetrical volume element
c
      nrm1  = nr-1
      nzm1  = nz-1
      if(itest.ge.2) then
       nro   = nr
       nzo   = nz
      else
       nro   = 9
       nzo   = 9
      endif
      drsqr = dr*dr
      dzsqr = dz*dz
      dvol  = drsqr*dz*two*pi
      isym  = int(sym)
      icm   = (1-isym)*(nz+1)/2
      zcm   = zero
      sz    = zero
      sr    = zero
c
c     initialisation of the constraint
c
      if(itcons.gt.0) qlamb = 0.0
c
c     storing particle numbers
c
      npart(1) = nchg
      npart(2) = nmass-nchg
      anumb  = npart(1)+npart(2)
c
c     compute maximum oscillator shell
c
      eocc(1) = 1.5*npart(1)**edr + deocc(1)
      eocc(2) = 1.5*npart(2)**edr + deocc(2)
c
      return
      end
CR-11-21-98
c-----inipai-----------------------------------------------------------
c
      subroutine inipai(fermip,fermin,ggp,ggn,nmass)
c
c     initialises pairing parameters.
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'dynpar.inc'
      include 'inipar.inc'
      include 'observ.inc'
      include 'levels.inc'
      include 'resden.inc'
c
      fermx(1)=fermip
      fermx(2)=fermin
      if(ipair.eq.0) then
       if(ggp.le.smal) ggp = 11.2 / sqrt(one*(nmass))
       if(ggn.le.smal) ggn = 11.2 / sqrt(one*(nmass))
      elseif(ipair.eq.1) then
c       if(ggp.le.smal) ggp = 23.0/nmass
       if(ggp.le.smal) ggp = 17.0/nmass
       if(ggn.le.smal) ggn = 23.0/nmass
      elseif(ipair.eq.2 .or. ipair.eq.3) then
       if(ipair.eq.3) then
         rspln(1)  = 0.0
         rspln(2)  = 0.0
         suv2rs(1) = 1.0
         suv2rs(2) = 1.0
         ilnrsp = 1
         if(ilnrsp.eq.1) then
           write(7,'(a)') ' INCLUDING linear response from L.N.'
         else
           write(7,'(a)') ' WITHOUT linear response from L.N.'
         endif
       endif
       call nucmat(rhonm,eovera,compr,effmas,asymme)
       write(6,'(a,g12.4)') ' in NUCMAT: m*/m=',effmas
       write(7,'(a,g12.4)') ' in NUCMAT: m*/m=',effmas
       ibcs = 3-ipair
       numsta(1)  = 0                            ! compute number of states
       do ist=1,nst(1)
         numsta(1) = watep(ist) + numsta(1)
       enddo
       rhoscp = zero
       ggp  = gavrg(ibcs,npart(1),npart(2),1,rhoscp,1,numsta(1),
     &                   effmas,0.58*0.960d0)              ! estimate surfp=0.96
       numsta(2)  = 0                            ! compute number of states
       do ist=nst(1)+1,nst(1)+nst(2)
         numsta(2) = watep(ist) + numsta(2)
       enddo
       rhoscn = zero
       ggn  = gavrg(ibcs,npart(1),npart(2),0,rhoscn,1,numsta(2),
     &                   effmas,0.58*0.960d0)              ! estimate surfp=0.96
       write(7,'(a,2i5,4g12.4)')
     &  ' GAVRG: numsta(1,2),rhoscp,rhoscn,ggp,ggn=',
     &    numsta,rhoscp,rhoscn,ggp,ggn
       write(7,'(a,4i5)') ' NST,NUMSTA=',nst,numsta
      elseif(ipair.eq.4) then
       if(ggp.eq.zero .or. ggn.eq.zero) then
	write(6,'(a,2f8.4,a)')
     &  ' Temperature(=GGP/N) =',ggp,ggn,' will not work. Stop!'
	stop ' zero temperature!'
       endif
      endif
*
***********************************************************************
*                                                                     *
*        initialize variables for the pairing part of the             *
*        calculation. if ipair=0 this block   block is.               *
*                   ipair -   index for pairing force as follows      *
*                             ipair=-1  no pairing                    *
*                             ipair=0  constant gap model             *
*                             ipair=1  constant force model           *
*                             ipair=2  constant force model, BCS      *
*                             ipair=3  Lipkin Nogami                  *
*                             ipair=4  thermal CHF                    *
*                   fer(q)-   1 st guess for the fermi energy for     *
*                             q=1 (protons) and q=2 (neutrons)        *
*                   gap(q)-   pairing gap (mev) for                   *
*                             q=1 (protons), and q=2 (neutrons)       *
*                   ggg(q)-    pairing strength (mev) for             *
*                             q=1 (protons), and q=2 (neutrons)       *
*                   vsqf  -   bcs occupation numbers                  *
*                                                                     *
***********************************************************************
*
      do i=1,2
        fer(i)=  0.0                           !?? guess -5
        feradd(i) = 0.0
        gap(i)=11.2 / sqrt(one*(nmass))
      enddo
*
      ggg(1)=ggp
      ggg(2)=ggn
      if(ipair.eq.0 .or. ipair.eq.4) then
         gap(1)=ggg(1)
         gap(2)=ggg(2)
         pairfc(1) = 0.0D0
         pairfc(2) = 0.0D0
      else
         gap(1) = 11.2 / sqrt(one*(nmass))
         gap(2) = gap(1)
      end if
      if(ipair.eq.5 .or. ipair.eq.7) then
        pairfc(1) = ggp
        pairfc(2) = ggn
      endif
      if(ipair.eq.9 .or. ipair.eq.11) then
        pairfc(1) = ggp
        pairfc(2) = ggn
        rspln(1)  = 0.0
        rspln(2)  = 0.0
        suv2rs(1) = 1.0
        suv2rs(2) = 1.0
        ilnrsp = 1
c        ilnrsp = 0
        if(ilnrsp.eq.1) then
           write(7,'(a)') ' INCLUDING linear response from L.N.'
        else
           write(7,'(a)') ' WITHOUT linear response from L.N.'
        endif
      endif
      if(ipair.ge.5) then
        do ist=1,nst(1)+nst(2)
         wtpair(ist) = 1.0
        enddo
      endif
c
      fer2(1) = ggp/4.0
      fer2(2) = ggn/4.0
c
      return
      end
CR-11-21-98
c
c-----inicon-----------------------------------------------------------
c
      subroutine inicon
c
c     initialises on 'qfac' the quadrupole constraint with
c     Saxon-Wood damping at the outer surface.
c     the damping factor as such is stored on 'fda'.
c
c     option "iflin=0": quadrupole as such is damped
c     option "iflin=1": linear asymptotic and softened damping
c
c     The "linear asymtotic" constraint is:
c
c                    4*z^2           r^2
c       Q = c*sqrt(1+-----)-c - ---------------
c                      c        sqrt(1+4*z^2/c)
c
c     where
c                2
c       c = 16  z
c                switch
c
c     and the "z_switch" indicates the z-value where the quadratic
c     behaviour switches to linear behaviour. "z_switch" is determined
c     as "swifac*1.2*A**1/3" and "swifac" comes from input.
c
c     the 'z' is taken relative to the center-of-mass 'zcm'. 
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
c
      data qr/0.0/                 ! to make it static
      data y40fac/0.846284375/     !    24/(16*SQRT(PI))
c
c-----------------------------------------------------------------------
c
      octfac = oh * sqrt(7.0/pi)
      a      = (npart(1)+npart(2))
      rcd    = 1.2*a**(edr)
      besoll = sqrt(5.0*pi)*qsoll/(3.0*rcd*rcd*a)
      rmin   = rcd+dalp
      rmin2  = rmin*rmin
      radd1  = besoll*rcd*sqrt(5.0/pi)/(2.0)
c      zswitc = swifac*rcd
c      cswitc = 16.0*zswitc*zswitc
      if (qsoll.gt.zero) then
        radd2  = sqrt( (qsoll           )/(0.5*a) )*0.5
        if (radd1.le.radd2) then
         radd = radd1
        else
         radd = radd2
        endif
      else if (qsoll.eq.zero) then
        radd  = zero
      else
        radd  =-0.5*radd1
      endif
      dterm = rmin2-radd*radd
      do iz=0,nz
      z      = ((iz)+sz-icm)*dz - zcm
      zsqr   =  z*z
      if (qsoll.ge.zero) then
        if (abs(z).ge.radd) then
         qr         = sqrt(dterm+2.0*radd*abs(z))
        else
         qr         = sqrt(zsqr+rmin2)
        endif
      endif
      do ir=0,nr
      r      = ((ir)+sr)*dr
      rsqr   =  r*r
      rr     = zsqr+rsqr
      sqrtrr = sqrt(rr)
c
c     the bare constraining force on save field 'qcfac'
c
      if(iflin.eq.1) then
c        sqswit = sqrt(one+4*zsqr/cswitc)
c        qcfac(iz,ir) = (cswitc*sqswit-cswitc-rsqr/sqswit**3)   ! asymp. linear
        qcfac(iz,ir) = z
      else
        if(lqcon.eq.2) then
          qcfac(iz,ir) = (zsqr + zsqr - rsqr)                ! cart. quadrupole
        elseif(lqcon.eq.3) then
          qcfac(iz,ir) = octfac*z*(zsqr - 1.5*rsqr)          ! spher. octupole
        elseif(lqcon.eq.4) then
          qcfac(iz,ir) =y40fac*(zsqr*zsqr-3.0*zsqr*rsqr+0.375*rsqr*rsqr)
        else
          stop ' this LQCON not provided '
        endif
      endif
!     additional for quadr. constr.
      IF(qlambquad.NE.0D0) qfacquad(iz,ir) = (zsqr+zsqr-rsqr)

c
c     the damping functional
c
      if (qsoll.lt.zero) then
        if (r.ge.radd) then
         qr         = sqrt(dterm+2.0*radd*r)
        else
         qr         = sqrt(rmin2+rsqr)
        endif
      endif
      dampfac = one/(one+exp((sqrtrr-qr)/dgam))
      qfac(iz,ir)  =  dampfac * qcfac(iz,ir)

!     additional for quadr. constr.
      IF(qlambquad.NE.0D0) 
     &   qfacquad(iz,ir) = qlambquad*dampfac * (zsqr+zsqr-rsqr)

c
      enddo
      enddo
      if(itest.gt.0) then
        call testpr(qcfac,0,nz,0,9,' qcfac')
        call testpr(qfac,0,nz,0,9,' qfac')
      endif
c
      return
      end

c
c-----inicon-----------------------------------------------------------
c
      subroutine inipotdamp
c
c     initialises on 'potdamp' the additional potential
c     to ensure value and derivative equal to zero in
c     the neighbourhood of the boundaries
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
      include 'fields.inc'
c
c-----------------------------------------------------------------------
c
c
      integer          ir,iz
      double precision r,z
      double precision sph_rsq,sph_Rmax,sph_Rdamp,dlt
      double precision rmax,zmax1,zmax2,zmin,zmax
      double precision damp_rad,damp_hgt
c
c     set parameters of damping field
c 
      damp_rad = 5.0D0     ! width of damping field in fm
      damp_hgt = 100.0D0   ! strength of damping field in MeV
c
c     ------------------------------------------------------------------
c     initialize damping field                    
c     ------------------------------------------------------------------
      rmax=(nr+sr)*dr
      zmax1=(nz+sz-icm)*dz
      zmax2=-(sz-icm)*dz

      if (isym.eq.0) then
         if (zmax1.ge.zmax2) then     ! find maximum z-value
            zmax=zmax1
            zmin=zmax2
         else
            zmax=zmax2
            zmin=zmax1
         endif
      else
         zmin = zmax1
      endif

      if (zmin.le.rmax) then          ! find maximum radius
         sph_Rmax = zmin
      else
         sph_Rmax = rmax
      endif
      
      sph_Rdamp = sph_Rmax - damp_rad

      do ir=0,nr
        r=(ir+sr)*dr
        do iz=0,nz
           z=(iz+sz-icm)*dz
           
           sph_rsq   = sqrt(r**2+z**2)
           dlt       = sph_rsq-sph_Rdamp

           if (sph_rsq.le.sph_Rdamp) then
              potdamp(iz,ir) = 0.0D0
           else if (sph_rsq.ge.sph_Rmax) then
              potdamp(iz,ir) = 1.0D0
           else                         
              potdamp(iz,ir) = 1.0D0-exp(dlt**2/(dlt**2-damp_rad**2))
           endif
           
        enddo
      enddo
c
c     weight damping field
c
      do ir=0,nr
         do iz=0,nz
            potdamp(iz,ir) = damp_hgt*potdamp(iz,ir)
         enddo
      enddo
c
c      call wrt_2D_scalar(potdamp,'potdamp___',1,'tst')
c
      return
      end

c-----idanuc------------------------------------------------------------
c
      subroutine idanuc
c
c     in this subroutine the matrices needed for the damping
c     are initialised using the derivative matrices of the subroutine
c     inider.
c     the damping operator d is given by:
c
c                              1              1
c      d = x0inv/e0inv -------------    ------------
c                         1+t_r/e0inv   1+t_z/e0inv
c
c     where t_r is the kinetic energy in r_direction and t_z in
c     z-direction in axial symmetry.
c
c     output via
c     common derivs:  zdnuc  -->  damping matrix, x-dir
c                     rdnuc  -->  damping matrix, y-dir,
c                                angular momentum term is included
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'derivs.inc'
      include 'wavefr.inc'
      dimension zaux(0:nzd,0:nzd),raux(0:nrd,0:nrd)
!      equivalence (zaux,psi),(raux,psi)
      equivalence (zaux,wrk1),(raux,wrk2)
c
c     set up unit matrix for 'e0inv=0'
c
      if(e0inv.eq.zero) then
       do ipi=1,4
         do i=0,nz
         do j=0,nz
           zdnuc(i,j,ipi) = zero
           if(i.eq.j) zdnuc(i,j,ipi) = x0inv
         enddo
         enddo
         do i=0,nr
         do j=0,nr
           rdnuc(i,j,ipi) = zero
           if(i.eq.j) rdnuc(i,j,ipi) = x0inv
         enddo
         enddo
       enddo
      else
c
c     now the standard case for damped iteration
c
      befac  = h2mav/(two*e0inv)
!      write(7,'(a,g12.4,a)') ' test idanuc: befac=',befac,' , derz2:'
c
      do ipi=1,4
c
c     nucleon damping in z-direction
c
        do i=0,nz
          do j=0,nz
            zdnuc(j,i,ipi) = -befac * derz2f(j,i,ipi)
          enddo
          zdnuc(i,i,ipi) = one + zdnuc(i,i,ipi)
        enddo
c
c     inversion yields final result for zdnuc
c
        call minv(zdnuc(0,0,ipi),zaux,nzd+1,nz+1)
c
        if (icutz.gt.1) then
          call matrml(zdnuc(0,0,ipi),zpro,zaux,nzd,nz)
          call matrml(zpro,zaux,zdnuc(0,0,ipi),nzd,nz)
        endif
        if (ipi.eq.2.or.ipi.eq.4) then
          do i=0,nz
            zdnuc(0,i,ipi) = 0.0D0
            zdnuc(i,0,ipi) = 0.0D0  ! addition PK
          enddo
        endif
c
      enddo
c
c     test damping
c
      if(itest.eq.1) then
       do ir=0,nr
       do iz=0,nz
c        wrk1(iz,ir) = exp(-iz*iz*one/20.0)
        wrk1(iz,ir) = max(1-abs(iz-ir),0)
        wrk2(iz,ir) = wrk1(iz,ir)
       enddo
       enddo
       call mulz(derz2f(0,0,1),wrk1,wrk2,-befac,1,-1)
       call mulz(zdnuc(0,0,1),wrk2,wrk3,one,0,-1)
       sum    = 0.0
       errmax = 0.0
       do ir=0,nr
       do iz=0,nz
         sum    = (wrk3(iz,ir)-wrk1(iz,ir))**2 + sum
         errmax = max(abs(wrk3(iz,ir)-wrk1(iz,ir)),errmax)
       enddo
       enddo
       write(7,'(a,3(/a,10(1pg20.12)))')
     &   ' test z-damping:',
     &   ' input:',(wrk1(iz,1),iz=0,9),
     &   ' kin. :',(wrk2(iz,1),iz=0,9),
     &   ' test :',(wrk3(iz,1),iz=0,9)
       write(7,'(a,2(1pg13.5))') 
     &   ' error: average,maximum=',sqrt(sum/((nr+1)*(nz+1))),errmax
      endif
c
c     write(7,'(a)') ' z-matrix in inidam:'
c     do  i=0,nz
c      write(7,'(2(a,i3))') ' row=',i,' ipi=',ipi
c      write(7,'(30(5g15.4/))') (zdnuc(i,j,ipi),j=0,nz)
c     enddo
c
c     the nucleon damping in r-direction
c
       do impsi=0,kmpsi
         do j=0,nr
         do i=0,nr
           rdnuc(i,j,impsi) = -befac * derr2af(i,j,impsi,1)
         enddo
         enddo
         do i=0,nr
           rdnuc(i,i,impsi) =  one + rdnuc(i,i,impsi)
         enddo
c
         call minv(rdnuc(0,0,impsi),raux,nrd+1,nr+1)
c
c
         if (icutr.gt.1) then
           call matrml(rdnuc(0,0,impsi),rpro,raux,nrd,nr)
           call matrml(rpro,raux,rdnuc(0,0,impsi),nrd,nr)
         endif
         if (impsi.gt.0) then
ctest       if (impsi.lt.0) then
           do i=0,nr
             rdnuc(0,i,impsi) = zero
             rdnuc(i,0,impsi) = zero ! addition PK
           enddo
         endif
c
c       test damping
c
        if(itest.eq.1) then
         do ir=0,nr
         do iz=0,nz
c          wrk1(iz,ir) = exp(-iz*iz*one/20.0)
          wrk1(iz,ir) = max(1-abs(iz-ir),0)
          wrk2(iz,ir) = wrk1(iz,ir)
         enddo
         enddo

         call mulr(derr2af(0,0,impsi,1),wrk1,wrk2,-befac,1,-1)
         call mulr(rdnuc(0,0,impsi),wrk2,wrk3,one,0,-1)
         sum    = 0.0
         errmax = 0.0
         do ir=0,nr
         do iz=1,nz
           sum    = (wrk3(iz,ir)-wrk1(iz,ir))**2 + sum
           errmax = max(abs(wrk3(iz,ir)-wrk1(iz,ir)),errmax)
         enddo
         enddo
         write(7,'(a,i3,a,3(/a,10(1pg20.12)))')
     &   ' test r-damping for m=',impsi,':',
     &   ' input:',(wrk1(iz,1),iz=0,9),
     &   ' kin. :',(wrk2(iz,1),iz=0,9),
     &   ' test :',(wrk3(iz,1),iz=0,9)
         write(7,'(a,2(1pg13.5))') 
     &   ' error: average,maximum=',sqrt(sum/((nr+1)*(nz+1))),errmax
        endif
c
      enddo
c
      endif
c
      return
      end
c
c-----idames------------------------------------------------------------
c
      subroutine idames
c
c     this subroutine initialises the damping matrices for the
c     meson damping. the meson damping is given by:
c
c               1                   1                    1
c            -------   *    ----------------- *  -----------------  ,
c            dams**2        1-delta_r/dams**2    1-delta_z/dams**2
c
c     where dams is a numerical damping mass and delta_r/z is the
c     is the delta operator in r/z direction in axial symmetry.
c
c     output via common:
c             rdmes --> damping in r-direction
c             zdmes --> damping in z-direction
c
c
      include 'paramr.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'wavefr.inc'
      dimension zaux(0:nzd,0:nzd),raux(0:nrd,0:nrd)
      equivalence (zaux,psi),(raux,psi)
c
c     workspaces , equivalenced to save space.
c
c      dimension zaux(0:nzd,0:nzd),raux(0:nrd,0:nrd)
c      equivalence (zaux,psi),(raux,psi)
c
      imes=1                    ! only coulomb is initialized now
c
      tms  = -dmsfac/cmssqr
c
      do i=0,nr
      do j=0,nr
        rdmes(j,i,imes) =  tms  * derr2af(j,i,0,1)
      enddo
      enddo
c
      do i=0,nr
        rdmes(i,i,imes)  =  rdmes(i,i,imes) + one
      enddo
c
      call minv(rdmes(0,0,imes),raux,nrd+1,nr+1)
      if(itest.ge.2) then
c       call testpr(rdmes,0,9,0,9,' rdmes')
      endif
c      call matrml(rpro,rdmes(0,0,imes),rdmes(0,0,imes),nrd,nr)
c
      do i=0,nz
      do j=0,nz
        zdmes(j,i,imes) = tms * derz2f(j,i,1)
      enddo
      enddo
c
      do i=0,nz
        zdmes(i,i,imes) = zdmes(i,i,imes) + one
      enddo
c
      call minv(zdmes(0,0,imes),zaux,nzd+1,nz+1)
      if(itest.ge.2) then
c       call testpr(zdmes,0,9,0,9,' zdmes')
      endif
c      call matrml(zpro,zdmes(0,0,imes),zdmes(0,0,imes),nzd,nz)
c
 50   continue                                            ! obsolete ?
c
      return
      end
c
c-----ssetup-----------------------------------------------------------
c
      subroutine ssetup
*
***********************************************************************
*                                                                     *
*     initialises the 'psi' with harmonic oscillator wave functions   *
*     in the upper components; lowers are set to zero.                *
*     initialises the bookkeeping for the buffering algorithm.        *
*     initialises pairing fields.                                     *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'inipar.inc'
      include 'observ.inc'
      include 'constr.inc'
c
      dimension  xkappa(2),xmue(2),spnc(2),eocct(2)
      dimension  phi(nnzxd1,0:nzd)
      logical trestart
      data trestart/.false./
c
c----------------------------------------------------------------------
c
      write(7,'(a)') ' SSETUP entered'
*
***********************************************************************
*                                                                     *
*        setup the deformed harmonic oscillator parameters in terms   *
*        of the mass and initial deformation parameter betain         *
*        compute the oscillator sizes and frequencies in the r and z  *
*        directions and scale eocc by the frequency                   *
*                                                                     *
***********************************************************************
*
      xpi=sqrt(4.0/(5.0*pi))
      amas=npart(1)+npart(2)
      hboma=41.0/amas**(edr)
      hbomz=hboma*exp(-xpi*betain)
      hbomr=hboma*exp(xpi/two*betain)
      brosc=sqrt(two*h2mav/hbomr)
      bzosc=sqrt(two*h2mav/hbomz)
      eocct(1)=eocc(1)*hboma
      eocct(2)=eocc(2)*hboma
      if(isym.eq.0) then
        octstp = octini*0.5*sqrt(7.0/pi)/(1.7*amas*amas)     ! 1.7*amas = R^3
      else
        octstp = zero
      endif
      write(6,'(a,2g12.4)') ' octini,octstp=',octini,octstp
*
***********************************************************************
*                                                                     *
*        parameters for l*s and l**2 force in nilsson model           *
*                                                                     *
***********************************************************************
*
       if(npart(1).le.50) then
           xkappa(1)=0.08
           xmue(1)  =0.0
       elseif(npart(1).le.82) then
           xkappa(1)=0.0637
           xmue(1)  =0.0
       else
           xkappa(1)=0.0577
           xmue(1)  =0.0
       endif
*
*
*
       if(npart(2).le.50) then
           xkappa(2)=0.08
           xmue(2)  =0.0
       elseif(npart(2).le.82) then
           xkappa(2)=0.0637
           xmue(2)  =0.0
       elseif(npart(2).le.126) then
           xkappa(2)=0.0637
           xmue(2)  =0.0
       else
           xkappa(2)=0.0635
           xmue(2)  =0.0
       endif
*
***********************************************************************
*                                                                     *
*     zero wavefn                                                     *
*                                                                     *
***********************************************************************
*
      do ist=1,npsid
      do is=1,nsd
      do ir=0,nrd
      do iz=0,nzd
        psi(iz,ir,is,ist) = zero
      enddo
      enddo
      enddo
      enddo
*
***********************************************************************
*                                                                     *
*        initialize the harmonic oscillator wavefunctions             *
*                                                                     *
***********************************************************************
*
      if(isym.eq.1) then 
        nzho=nz
      else
        nzho=(nz/2)+1
      endif
      do j=0,nzho
        call zwf(j,bzosc,phi)
      enddo
*
***********************************************************************
*                                                                     *
*        set loop to allocate the single particle states in terms     *
*        of axially symmetric harmonic oscillators                    *
*        new variables are as follows:                                *
*                      lst     -  total no. of states                 *
*                      spn(*)  -  normalization of wavefunctions      *
*                      iofk    -  pointer to the 1 st index of a      *
*                                 state of charge q                   *
*                      lbuf    -  buffer index for direct i/o         *
*                                                                     *
***********************************************************************
*
      do i=1,npsid
        spn(i)=0.0
      enddo
      lst=0
      iofk=1
      kold=0
      iqold=1
      nom=1
      mhi(1)=1
      mlo(1)=1
      miq(1)=iqold
      mim(1)=kold
      jst=0
      mst=0
      lno=locc(1)
      locou=1
*
***********************************************************************
*                                                                     *
*        main loop over charge states begins here                     *
*        for each q mmax(q) is the maximum number of allowed  m       *
*        values                                                       *
*        note: m is k-1/2, k being the jz quantum number              *
*                                                                     *
***********************************************************************


* isTtmp = 1 only half the particles
* isTtmp = 0  full number of particles

      if(isTSym.ge.1)then
      isTtmp = 0
      else
      isTtmp = 1
      endif

      do 1000 iq=1,2
      mmaxq=mmax(iq)-1
      mminq= -mmax(iq)*(1-isTtmp)
      
***************
* Test: Alternative intialisierungsreihenfolge
*      if(isTSym.eq.0) then
*       mmaxreal = 2*mmaxq +1
*      else
*       mmaxreal = mmaxq
*      endif
*      do 500 k2=0,mmaxreal
*      if(isTSym.eq.1) then
*        k1 = k2
*       else
*         if(mod(k2,2).eq.0) then
*           k1 = 0.5*k2
*         else
*           k1 = -0.5*(k2+1)
*         endif
*       endif
*
* Ende Test 
      do 500 k1=mminq,mmaxq
c       do 500 k1=mmaxq,mminq,-1
c      k1 = -k2
c       k1=-k2-1
c      do 500 k1=1,mmaxq
*      k1=mm-1
*
*
***********************************************************************
*                                                                     *
*        nnz and nnr are the principle quantum numbers of the         *
*        cylindrical oscillators- nnzxd and nnrxd are their maximum   *
*        allowed values                                               *
*                                                                     *
***********************************************************************
*
      do 400 nnz=1,nnzxd
cold      nnzm1=nnz-1
      do 350 nnr=1,nnrxd
cold      nnrm1=nnr-1
*
***********************************************************************
*                                                                     *
*        ispgen counts the two possible spin orientations per         *
*        given mm. ispgen=1 occupies dominantly the upper             *
*        component whereas ispgen=2 occupies the lower one.           *
*        small lower, or upper, components are added by per-          *
*        turbative treatment of the l*s force.                        *
*                                                                     *
***********************************************************************
*
      do 300 ispgen=1,2
      spin=1.5-(ispgen)
      mm1=abs(k1)+ispgen-1
      if(k1.lt.zero) then 
        mm1 = mm1-1
      endif     
      mmr=mm1+1
*
***********************************************************************
*                                                                     *
*         ecur is the single particle energy of the deformed          *
*         oscillator state plus l*s and l**2 force.                   *
*         this must be less than eocc for each charge state           *
*         to be included in the basis                                 *
*         new variables as follows:                                   *
*                       mpsi(*)  -  jz quantum number minus 1/2       *
*                                    note that upper comp. has m=mpsi *
*                                    but lower comp. has m=mpsi+1     *
*                       npar(*)  -  z parity of the upper comp.       *
*                       wate(*)  -  occupation numbers in the         *
*                                   absence of pairing                *
*                       icharg(*)   isospin quantum number            *
*                                                                     *
***********************************************************************
*
      ecur= hbomz*(nnz- 0.5)+ hbomr*(2*nnr+ mmr-2)
     x      -(spin+spin)*hboma*mm1*xkappa(iq)
     x      -hboma*mm1*mm1*xkappa(iq)*xmue(iq)
*
***********************************************************************
*                                                                     *
*        if block for selecting the initial occupied states           *
*                                                                     *
***********************************************************************
*
      isigis  = 3-2*ispgen
      if (ecur.le.eocct(iq).and.mst+1.eq.lno) then
       mst   = mst+1
       locou = locou+1
       lno   = locc(locou)
       goto 300
      else if(ecur.le.eocct(iq)) then
         mst         = mst+1
         lst         = lst+1
         if(lst.gt.npsid) then
            write(7,'(a,i4,a,i4)') ' error: initial states',
     &      lst,' are too much, max. allowed:',npsid
            stop 'ssetup: initial states too much'
         end if
         jst         = jst+1
         spenrg(lst) = (ecur-70.0)
         mpsi(lst)   = k1
         nqr(lst)    = nnr
         nqz(lst)    = nnz
c         sigis       = 3.0-2.0*(ispgen)
         isigis      = 3-2*ispgen
         if (k1+1.gt.kmpsi) then
           write(7,'(a,i3)')
     &         ' kmpsi (dim of debe for ang mom) too small:',kmpsi
           stop
         endif
c        nrpar(lst)  = isigis*(1-mod(k1,2)*2)
         npar(lst)   = isigis*(mod(nnz,2)*2-1)*sign(1,k1)
* djd
* change to nparz(1,lst)=nparz(2,lst)=1 for sym=0, and as is otherwise
*
         if(isym.eq.1)then
            nparz(1,lst)= (3-npar(lst))/2
            nparz(2,lst)= 3-nparz(1,lst)
         elseif(isym.eq.0)then                   ! djd addition
           nparz(1,lst)= 1                       ! djd addition
           nparz(2,lst)= 1                       ! djd addition
         end if

c axtest npar(lst)   = isigis*(mod(1,2)*2-1)
         nrpar(lst)=   (1-mod(abs(k1),2)*2)
c        npar(lst)   = (mod(nnz,2)*2-1)
         icharg(lst) = iq
codd     wate(lst)   = two
         wate(lst)   = one*(1+isTtmp)
         watep(lst)  = one
c         xkim(lst)   = zero
*
***********************************************************************
*                                                                     *
*        construction of the wave function, dominant spinor part      *
*                                                                     *
***********************************************************************
*
         phinor=one/sqrt(pi+pi)
         do ir=0,nr
           r      = ((ir)+sr)*dr
           r2     = r*r
           call radwf(ir,br	osc,nnr,mmr,psirad)
           do iz=0,nz
             z      = ((iz)+sz-icm)*dz
             izp=iz-icm     !+1          cbm
             izm=icm-iz     !+2-ispgen   cbm
             octadd = octstp*z*(z*z-1.5*r2)
             if(isym.eq.0) then
               if(iz.lt.icm) then
c                 psi(iz,ir,ispgen,jst) = (one-octstp*qfac(iz,ir))
                 psi(iz,ir,ispgen,jst) = (one-octadd)
     x                 *float(mod(nnz,2)*2-1)*psirad*phi(nnz,izm)
               else
                 psi(iz,ir,ispgen,jst) =
     x                 (one-octadd)*psirad*phi(nnz,izp)
c     x                 (one-octstp*qfac(iz,ir))*psirad*phi(nnz,izp)
               endif
             else
                psi(iz,ir,ispgen,jst) = psirad*phi(nnz,iz)
             end if
             psi(iz,ir,ispgen,jst)=phinor*psi(iz,ir,ispgen,jst)
           enddo
         enddo
*
*
         if(itest.ge.1) then
           write(7,*) ' 1. time in ssetup: lst,jst,k1,iq=',
     x               lst,jst,k1,iq
           if (ispgen.eq.1) then
           call testpr(psi(0,0,1,jst),0,nzo,0,nro,' psi 1.comp')
           else if (ispgen.eq.2) then
           call testpr(psi(0,0,2,jst),0,nzo,0,nro,' psi 2.comp')
           endif
         endif
         ielse = 3 - ispgen
         do ir=0,nr
         do iz=0,nz
           psi(iz,ir,ielse,jst)=zero
         enddo
         enddo
*
***********************************************************************
*                                                                     *
*        calculate the normalization of the single particle state     *
*        and renormalize the wavefunctions                            *
*                                                                     *
***********************************************************************
*
c        projection of the last component in momentum space
c
c        goto 417
codd
         if(mpsi(jst).ge.zero) then
         muplo  =  mpsi(jst)+ispgen-1
         npuplo =  npar(jst)*isigis
         else
         muplo = abs(mpsi(jst)+1)+ispgen -1
         npuplo =  (-1)*npar(jst)*isigis
         endif
codd
c        if (npuplo.gt.0) then
         if ((npuplo.gt.0 .or. icutz.gt.1)
     &       .and.(mod(muplo,2).eq.0 .or. icutr.gt.1)) then
          do ir=0,nr
          do iz=0,nz
            wrk2(iz,ir) = psi(iz,ir,ispgen,jst)
          enddo
          enddo
          call mulz(zpro,wrk2,wrk1,one,0,-1)
          call mulr(rpro,wrk1,wrk2,one,0,-1)
c        goto 317
         else if (npuplo.gt.0 .or. icutz.gt.1)  then
          do ir=0,nr
          do iz=0,nz
            wrk1(iz,ir) = psi(iz,ir,ispgen,jst)
          enddo
          enddo
          call mulz(zpro,wrk1,wrk2,one,0,-1)
c        else
         else if (mod(muplo,2).eq.0 .or. icutr.gt.1) then
          do ir=0,nr
          do iz=0,nz
            wrk1(iz,ir) = psi(iz,ir,ispgen,jst)
          enddo
          enddo
          call mulr(rpro,wrk1,wrk2,one,0,-1)
         endif
c317     continue
         if ((npuplo.gt.0 .or. icutz.gt.1).or.
     &       (mod(muplo,2).eq.0 .or. icutr.gt.1)) then
         do ir=0,nr
         do iz=0,nz
           if (muplo.gt.0.and.ir.eq.0) then
            psi(iz,ir,ispgen,jst) = 0.0
           else if (npuplo.lt.0.and.iz.eq.0) then
            psi(iz,ir,ispgen,jst) = 0.0
           else
            psi(iz,ir,ispgen,jst) = wrk2(iz,ir)
           endif
          enddo
          enddo
         if (mod(muplo,2).eq.1) then
          do iz=0,nz
            psi(iz,nr,ispgen,jst) = 0.0
          enddo
         endif
         endif
c417      continue
c
******************************
* Falls m kleiner Null, wende Zeispiegel-Operator an
******************************
         if(k1.lt.0) then
         do iz=0,nz
          do ir=0,nr 
           copy = psi(iz,ir,1,lst) 
           psi(iz,ir,1,lst)=(-1)*psi(iz,ir,2,lst)
            psi(iz,ir,2,lst)=copy
          enddo
         enddo
         endif

codd

         spn(lst)=zero
         do 241 is=1,2
          spnc(is) = wfovpc(psi(0,0,is,lst))
c         write(7,*) ' lst=',lst,' is=',is,' spnc=',spnc(is)
241      continue
         spn(lst)=spnc(1)+spnc(2)
         cmpn(1,lst)=spnc(1)/spn(lst)
         cmpn(2,lst)=spnc(2)/spn(lst)
         spnc(1)=1.0/spn(lst)
         spnc(2)=spnc(1)
         if(itest.ge.2) then
          write(7,'(a,i3/3(a,f10.4/))') ' lst=',lst,' norm,1.comp',
     &    cmpn(1,lst),' norm,2.comp',cmpn(2,lst),' tot norm',spn(lst)
         endif
         do is=1,2
         do ir=0,nr
         do iz=0,nz
           psi(iz,ir,is,jst)=psi(iz,ir,is,jst)/sqrt(spn(lst))
         enddo
         enddo
         enddo
         spn(lst)=sqrt(spn(lst))
c        xrwt(0)=zero
         if(itest.ge.2) then
           write(7,*) ' after proj. in ssetup: lst,jst,k1,iq=',
     x               lst,jst,k1,iq
           call testpr(psi(0,0,1,jst),0,nz,0,nr,' psi 1.comp')
           call testpr(psi(0,0,2,jst),0,nz,0,nr,' psi 2.comp')
c        plot of the wavefunctions
c          call psiplt(psi(0,0,1,jst),lst,jst,k1,iq,iou)
         endif
codd     testprint of the fields
*           write(7,*) mpsi(jst), ' ' , jst
*           call testpr(psi(0,0,1,jst),0,nz,0,nr,' psi 1.comp')
*           call testpr(psi(0,0,2,jst),0,nz,0,nr,' psi 2.comp')

codd
*
***********************************************************************
*                                                                     *
*        the indices for a wavefunctins having a definite value of    *
*        mm must be stored. each buffer is comprized of complete      *
*        blocks of wavefunctions for a particular mm as follows       *
*        mhi(nom) - array storing the upper index for the             *
*                        wavefunctions in the buffer lbuf             *
*        mlo(nom) - array storing the upper index for the             *
*                        wavefunction in the buffer lbuf              *
*        miq(nom) - array storing the iq value for the wave-          *
*                        function in the buffer lbuf                  *
*        mim(nom) - array storing the m value for the wave-           *
*                        function in the buffer lbuf                  *
*        nom      - number of mm blocks in the buffer lbuf            *
*                                                                     *
*        write out the wavefunctions to disk if the buffer is full    *
*                                                                     *
***********************************************************************
*
         if (k1.ne.kold.or.iq.ne.iqold) then
            mhi(nom)=jst-1
            nom=nom+1
            mlo(nom)=jst
            kold=k1
            iqold=iq
            miq(nom)=iqold
            mim(nom)=kold
         end if
         nbuf(jst) = nom

*
***********************************************************************
*                                                                     *
*        end if block for allowed quantum numbers                     *
*        end nnz and nnr loops                                        *
*                                                                     *
***********************************************************************
*
      end if
300   continue
350   continue
400   continue
*
***********************************************************************
*                                                                     *
*        end loop over m quantum number                               *
*        set the number of single particle states in nst              *
*                                                                     *
***********************************************************************
*
  500 continue
      nst(iq)= lst
      iokk(iq)=iofk
      iofk=iofk+lst
*
***********************************************************************
*                                                                     *
*         end loop over iq                                            *
*                                                                     *
***********************************************************************
*
 1000 continue
      nst(2)=nst(2)-nst(1)
      mhi(nom)=jst
      len=mhi(nom)
      jst=0
*
***********************************************************************
*                                                                     *
*        calculate the occupation numbers and the ordering of the     *
*        levels. note:                                                *
*                     iofk is the pointer to the 1 st index           *
*                     of state q                                      *
*                                                                     *
***********************************************************************
*
      do iq=1,2
        iofk=iokk(iq)
        call vsqtr(vsqf(iofk),wate(iofk),nst(iq),
     x             spenrg(iofk),npart(iq),iorder(iofk))
      enddo
      ifocsv = 0                                   ! reset flag
      write(6,'(a)') ' reset IFOCSV=0'
c
c     print levels
c
      write(7,'(/a,2f8.3//18x,a,12x,a/a,a)')
     &  ' The initial levels, chosen up to EOCC(1/2)=',EOCC,
     &  'parities','components',
     &  '   Nr Iso  j_z total  r    z   v^2  1.   2.',
     &  ' osc.quant'
      do ist=1,len
        write(7,'(i5,i4,f6.1,3i5,f6.3,2f5.1,f7.3)')
     &   ist,icharg(ist),mpsi(ist)+0.5,npar(ist),nrpar(ist),
     &   nparz(1,ist),
     &   vsqf(ist),cmpn(1,ist),cmpn(2,ist),(spenrg(ist)+70.0)/hboma
      enddo
      write(7,'(a)') ' '
c
*
***********************************************************************
*                                                                     *
*        initial pairing gap for constant strength models             *
*                                                                     *
***********************************************************************
*
*
c
c     test particle number
c
      na = 0
      do iq=1,2
        sumprt = zero
        sumnon = zero
        sumtot = zero
        do n=1,nst(iq)
          na  = 1 + na
          v2  = vsqf(na)
          wgt = wate(na)
          sumtot = wgt + sumtot
          sumprt = v2*wgt + sumprt
          sumnon = (one-v2)*wgt + sumnon
        enddo
        write(6,'(a,i2,a,g12.4/14x,2(a,g12.4))')
     &  ' after ssetup: for type=',iq,' nr. of particles=',sumprt,
     &  ' nr. of empty places=',sumnon,
     &  ' nr. of possible particles=',sumtot
        if(sumtot.lt.npart(iq)-err)
     &   stop ' not enough states to reach particle number'
cold        if((ipair.eq.5 .or. ipair.eq.7) .and. sumtot.lt.npart(iq)*two)
        if(ipair.ge.5 .and. sumtot.lt.npart(iq)*two)
     &   stop ' IPAIR=5 or 7 requires larger space, enhance DEOCC!'
      enddo
c
c     branch to restart: is activated if 'restart2D.data' exists
c     
      IF(KRESTA) THEN
        INQUIRE(FILE='restart2D.data',EXIST=trestart)
        IF(trestart) then
          call read2D()
          return
        endif
      ENDIF
*
*
      WRITE(6,'(a)') 'buffering scheme:'
      WRITE(6,'(a)') ' im  mhi  mlo  miq  mim'
      DO im=1,nom
        WRITE(6,'(5i5)') im,mhi(im),mlo(im),miq(im),mim(im)
      END DO
      WRITE(6,'(a,5000i4)') 'nbuf,iso=1:',nbuf(1:nst(1))
      WRITE(6,'(a,5000i4)') 'nbuf,iso=2:',nbuf(nst(1)+1:nst(1)+nst(2))

      write(7,'(a)') ' SSETUP finished'
      return
      end
c
c-----zwf   -----------------------------------------------------------
c
      subroutine zwf(i,bz,phi)
*
***********************************************************************
*                                                                     *
*         zwf (double precision) calculates the normalized            *
*         one-dimensional wave functions phi(n,i).                    *
*         n denotes the number of nodes along the z-axis.             *
*         i labels the z-coordinate.                                  *
*                                                                     *
*         the wavefunctions are normalized such that the integral     *
*         from minus infinity to infinity of  dz phi**2 = 1.          *
*                                                                     *
*         bz - harmonic oscillator size constant                      *
*         zi - z coordinate                                           *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
*
      dimension  phi(nnzxd1,0:nzd)
*
      root2=sqrt(two)
      rootpi=sqrt(pi)
cpgr      z=(i+sz-icm)*dz
      z = (i) * dz
      x=z/bz
      earg = x*x/two
      if (earg.le.60.0*one) then
       expf = exp(-earg)
      else
       expf = zero
      endif
      phi(1,i)=expf/sqrt(rootpi*bz)
      phi(2,i)=root2*x*phi(1,i)
      if(nnzxd.lt.2)   go to 20
      do  nnz=2,nnzxd
        nzt=nnz-1
        en=nzt
        phi(nnz+1,i)= sqrt(two/(en+one))*x*phi(nnz,i)
     &               -sqrt(en/(en+one))*phi(nzt,i)
      enddo
   20 continue
*
      return
      end
c
c-----radwf -----------------------------------------------------------
c
      subroutine radwf(i,br,nnr,mmr,psirad)
*
***********************************************************************
*                                                                     *
*         radwf (double precision) calculates the normalized          *
*         radial wave functions psi(n,l,i).                           *
*         n denotes the number of radial nodes.                       *
*         l denotes the magnitude of the z-component of angular       *
*         momentum.
*
*         i labels the cylindrical coordinate r.                      *
*
*         the wave functions are normalized such that the integral    *
*         from zero to infinity  dr r psi**2  is equal to one.        *
*                                                                     *
*         br - harmonic oscillator size constant                      *
*         ri - cylindrical coordinate                                 *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
*
      dimension  psi(nnrxd1), gamma(mxd1)
c
      data ifirst/-1367/
      data gamma/mxd1*0.0/                  ! to make this a static variable
c
c     initialize the factorials
c
      if(ifirst.eq.-1367) then
       if(mxd.gt.33) write(7,'(a,i5,a)')
     &  ' mxd=',mxd,' probably to large for single precision'
       gamma(1)=one
       do n=1,mxd
         xxx=(n)
         gamma(n+1)=xxx*gamma(n)
       enddo
       ifirst = 0
      endif
*
      r=((i)+sr)*dr
      rooty=r/br
      if (rooty.lt.one*1.e-4) then
       y = zero
      else
       y = rooty**2
      endif
      earg=y*oh
      if (earg.le.60.0*one) then
       expf = exp(-earg)
      else
       expf = zero
      endif
      do 20 ll=1,mmr
      l=ll-1
      el=l
      if (l.eq.0) then
       ryl = one
      else if (rooty.lt.1.d-4) then
       ryl = zero
      else
       ryl = rooty**l
      endif
      psi(1) = sqrt(two/gamma(l+1))/br*ryl*expf
      psi(2) = (el+one-y)/sqrt(el+one)*psi(1)
      if(nnr.lt.2)   go to 15
      do 10 nnrr=2,nnr
      nrt=nnrr-1
      enr=nrt
      psi(nnrr+1)=((two*enr+el+one-y)*psi(nnrr)
     &      -sqrt((enr+el)*enr)*psi(nrt))/sqrt((enr+el+one)*(enr+one))
   10 continue
   15 continue
   20 continue
      psirad = psi(nnr)
*
      return
      end
