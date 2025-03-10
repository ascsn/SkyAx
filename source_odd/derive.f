c=======================================================================
c
c     initialization of the derivative matrices 
c     
c       drv    :  array to store the derivative matrix
c
c       in     :  order of derivative
c       ipnt   :  number of points into each direction from central point
c                 or -1 for fourier derivatives
c       iparmin:  parity at lower boundary or 0 for reflecting
c       iparmax:  parity at upper boundary or 0 for reflecting
c
c       dx  :  spacing of gridpoints
c       inx :  number of gridpoints
c       inxd:  total dimension of array used to store derivatives
c              (array is only used up to index inx)
c
c
c     to do:  - fourier derivatives only tested up to order 2
c             - fourier derivatives with negative parity at 
c               upper end of box not yet supported 
c               (it should be possible to obtain them by mirroring)
c
c=======================================================================
c
      subroutine ini_drv_weights(der,
     &                           in,ipnt,
     &                           iparmin,iparmax,
     &                           dx,inx,inxd)
c
      parameter (imaxdrvrng=20)           ! maximum range of derivatives
      
      dimension der(0:inxd,0:inxd),
     &          A(0:2*imaxdrvrng,0:2*imaxdrvrng), 
     &          b(0:2*imaxdrvrng),
     &          jeff(0:4*inxd),
     &          isign(0:4*inxd)

c     ------------------------------------------------------------------
c     check consistency of input parameters
c     ------------------------------------------------------------------
      if (in.le.0) then   
         write(6,'(a)') 
     &        'order of derivative must be >0'
         stop 'error in ini_drv_weights'
      endif
c
      if ((ipnt.ne.-1).and.(ipnt.lt.(in+1)/2)) then   
         write(6,'(i2 a i2 a)') 
     &        in,'. order derivative needs at least ',ipnt,' points'
         stop 'error in ini_drv_weights'
      endif
c     
      if ((ipnt.lt.-1).or.(ipnt.gt.imaxdrvrng).or.(ipnt.eq.0)) then
         write(6,'(ai2)') 
     &        'number of points must be -1,1,2,...,',imaxdrvrng
         stop 'error in ini_drv_weights'
      endif
c
      if ((iparmin.lt.-1).or.(iparmin.gt.1)) then
         write(6,'(a)') 
     &        'parity at lower boundary must be -1,0,+1'
         stop 'error in ini_drv_weights'
      endif
c
      if ((iparmax.lt.-1).or.(iparmax.gt.1)) then
         write(6,'(a)') 
     &        'parity at upper boundary must be -1,0,+1'
         stop 'error in ini_drv_weights'
      endif    
c
      if (dx.le.0.0) then
         write(6,'(a)') 
     &        'spacing of gridpoints must be >0'
         stop 'error in ini_drv_weights'
      endif            
c
      if (inx.gt.inxd) then
         write(6,'(a)') 
     &        'size of matrix exceeded maximum dimension'
         write(6,'(a i3 a i3 a)')
     &        '(size: ',inx,', allowed: ', inxd,')' 
         stop 'error in ini_drv_weights'
      endif 

      if (ipnt.eq.-1) then   
c     ------------------------------------------------------------------
c     determin weights for the n-th fourier derivative 
c     ------------------------------------------------------------------
      pi=3.14159265358979
c     ------------------------------------------------------------------
c     the odd fourier derivatives
c     ------------------------------------------------------------------
      if (mod(in,2).eq.1) then

      if ((iparmin.eq.+1).and.(iparmax.ne.-1)) then
         dn    = (inx)
         afac  = pi/dn
         pfac  = (2.0/dn)*(afac/dx)**in
         sfac  = (-1)**((in+1)/2)
         do l=0,inx
            dl  = (l)
            do i=0,inx
               di     = (i)
               sum    = 0.0
               do j=0,inx
                  dj  = (j)
                  sum  = sum  + dj**in*sin(afac*dj*di)*cos(afac*dj*dl)
               enddo
               if (l.eq.0.or.l.eq.inx) then
                  der(i,l) =  0.5* sfac*pfac * sum
               else
                  der(i,l) =       sfac*pfac * sum
               endif
            enddo
         enddo
      elseif ((iparmin.eq.-1).and.(iparmax.ne.-1)) then
         dn    = (inx)
         afac  = pi/dn
         pfac  = (2.0/dn)*(afac/dx)**in
         sfac  = (-1)**(in/2)
         do l=0,inx
            dl  = (l)
            do i=0,inx
               di     =  (i)
               sum    =  0.0
               do j=0,inx
                  dj  = (j)
                  sum = sum + dj**in*cos(afac*dj*di)*sin(afac*dj*dl)
               enddo
               if (l.eq.0.or.l.eq.inx) then
                  der(i,l) =  0.5* sfac*pfac * sum
               else
                  der(i,l) =       sfac*pfac * sum
               endif
            enddo
         enddo
      elseif ((iparmin.eq.0).and.(iparmax.eq.0)) then
         if (mod(inx+1,2).ne.0) then
            write(6,'(a)') ' warning: nz has to be odd !'
            stop
         endif
         icn   = (inx+1)/2
         dn    = (icn)
         afac  = pi/dn
         pfac  = (1.0/dn)*(afac/dx)**in
         sfac  = (-1)**((in+3)/2)
         do l=0,inx
            dl  = (l)
            do i=0,inx
               di     = (i)
               sum    = 0.0
               do j=1,icn
                  dj  = (j)
                  sum  = sum  + dj**in*sin(afac*dj*(dl-di))
               enddo
               if (l.eq.0.or.l.eq.inx) then
                  der(i,l) =  0.5* sfac*pfac * sum
               else
                  der(i,l) =       sfac*pfac * sum
               endif
            enddo
         enddo
      else
         write(6,'(a,i2,a,i2,a)') 
     &        'parity ',iparmin,' at lower and ',iparmax,
     &        'at upper boundary not supported'
         stop 'error in ini_drv_weights'
      endif

c     ------------------------------------------------------------------
c     the even fourier derivatives
c     ------------------------------------------------------------------
      elseif (mod(in,2).eq.0) then

      if ((iparmin.eq.+1).and.(iparmax.ne.-1)) then
         dn    = (inx)
         afac  = pi/dn
         pfac  = (2.0/dn)*(afac/dx)**in
         sfac  = (-1)**(in/2)
         do l=0,inx
            dl  = (l)
            do i=0,inx
               di     = (i)
               sum    = 0.0
               do j=0,inx
                  dj  = (j)
                  sum  = sum  + dj**in*cos(afac*dj*di)*cos(afac*dj*dl)
               enddo
               if (l.eq.0.or.l.eq.inx) then
                  der(i,l) =  0.5* sfac*pfac * sum
               else
                  der(i,l) =       sfac*pfac * sum
               endif
            enddo
         enddo
      elseif ((iparmin.eq.-1).and.(iparmax.ne.-1)) then
         dn    = (inx)
         afac  = pi/dn
         pfac  = (2.0/dn)*(afac/dx)**in
         sfac  = (-1)**(in/2)
         do l=0,inx
            dl  = (l)
            do i=0,inx
               di     =  (i)
               sum    =  0.0
               do j=0,inx
                  dj  = (j)
                  sum = sum + dj**in*sin(afac*dj*di)*sin(afac*dj*dl)
               enddo
               if (l.eq.0.or.l.eq.inx) then
                  der(i,l) =  0.5* sfac*pfac * sum
               else
                  der(i,l) =       sfac*pfac * sum
               endif
            enddo
         enddo
      elseif ((iparmin.eq. 0).and.(iparmax.eq. 0)) then
         if (mod(inx+1,2).ne.0) then
            write(6,'(a)') ' warning: nz has to be odd !'
            stop
         endif
         icn   = (inx+1)/2
         dn    = (icn)
         afac  = pi/dn
         pfac  = (1.0/dn)*(afac/dx)**in
         sfac  = (-1)**((in+1)/2)
         do l=0,inx
            dl  = (l)
            do i=0,inx
               di     = (i)
               sum    = 0.0
               do j=1,icn
                  dj  = (j)
                  sum  = sum  + dj**in*cos(afac*dj*(dl-di))
               enddo
               if (l.eq.0.or.l.eq.inx) then
                  der(i,l) =  0.5* sfac*pfac * sum
               else
                  der(i,l) =       sfac*pfac * sum
               endif
            enddo
         enddo
      else
         write(6,'(a,i2,a,i2,a)') 
     &        'parity ',iparmin,' at lower and ',iparmax,
     &        'at upper boundary not supported'
         stop 'error in ini_drv_weights' 
      endif

      endif  ! end of calculation of the even fourier derivatives
         
      else             
c     ------------------------------------------------------------------
c     determin weights for the n-th finite derivative using p points 
c     into each direction from the central point ((2*p+1)-formula).
c     ------------------------------------------------------------------
      
c     setup the linear equation system which defines the weights
      do j=0,2*ipnt
         do i=-ipnt,ipnt
            A(j,i+ipnt) = i**j
         enddo
         if (j.eq.in) then
            b(j)=1.0 
         else 
            b(j)=0.0
         endif
      enddo

c     forward gauss elimination
      do i=0,2*ipnt
         do j=i+1,2*ipnt
            buf = A(j,i) / A(i,i)
            do k=i+1,2*ipnt
               A(j,k) = A(j,k) - A(i,k) * buf
            enddo
            b(j) = b(j) - b(i) * buf
         enddo
      enddo

c     backward substitution
      do i=2*ipnt,0,-1
         buf = 0.0
         do j=i+1,2*ipnt
            buf = buf + A(i,j) * b(j)
         enddo
         b(i) = (b(i) - buf) / A(i,i)
      enddo

c     calculate faculty of order n
      buf=1.0
      do i=1,in 
         buf = 1.0/dx * buf * i
      enddo

c     compose the weights for the finite differences
      do i=0,2*ipnt
         b(i) = buf * b(i)
      enddo
         
c     setup finite derivative matrix
      do i=0,inx
         do j=0,inx
            der(i,j)=0.0
         enddo
      enddo
      
      do i=0,inx
        do j=-ipnt,ipnt
           ! straight forward mid region
           if (((i+j).ge.0).and.((i+j).le.inx)) then
              der(i,i+j)    = der(i,i+j)         +         b(j+ipnt)
           ! take respect of parity at lower boundary
           elseif ((i+j).lt.0) then
              der(i,-(i+j)) = der(i,-(i+j))      + iparmin*b(j+ipnt)
           elseif ((i+j).gt.inx) then
              der(i,2*inx-(i+j)) 
     &                      = der(i,2*inx-(i+j)) + iparmax*b(j+ipnt)
           endif
        enddo
      enddo

      endif   ! end of finite definition of derivatives

c     finally output the derivative matrix 
c      do i=0,inx
c         write(6,'(31f7.3)') (der(i,j),j=0,inx)
c      enddo

      return
      end
