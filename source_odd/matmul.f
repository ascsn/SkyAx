c=======================================================================
c
c     mulz(a,x,y,factor,iaccu,iband)
c
c     matrix multiplication for z-direction.
c     computes
c
c          factor*(a*x)+ y --> y    for   iaccu=1
c          factor*(a*x)    --> y    for   iaccu=0
c
c     a is assumed to be a bandmatrix with A_{i,i+j}=0 for |j|>iband
c     or a full matrix if iband==-1
c
c=======================================================================
      subroutine mulz(a,x,y,factor,iaccu,iband)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc' ! nrd,nzd
      include 'dynpar.inc' ! nr,nz
c
c     ------------------------------------------------------------------
c
      double precision a(0:nzd,0:nzd),x(0:nzd,0:nrd),y(0:nzd,0:nrd)
      double precision acol(0:nzd)
      double precision sum,sum1,sum2,sum3,sum4
      double precision factor
      integer          iaccu,iband
      integer          iystart
      integer          ir,ix,iy,iz,ixmin,ixmax
c
c     ------------------------------------------------------------------
c
c     set target field to zero in case of no accumulation
c     ------------------------------------------------------------------
      if (iaccu.eq.0) then
         do ir=0,nr
            do iz=0,nz
               y(iz,ir) = zero
            enddo
         enddo
      endif
c     
c     full matrix multiplication
c     ------------------------------------------------------------------
      if (iband.eq.-1) then
         if (abs(factor).gt.smal) then
cV1            do ir=0,nr
cV1               do iy=0,nz
cV1                  sum = zero 
cV1                  do ix=0,nz
cV1                     sum     = a(iy,ix) * x(ix,ir) + sum
cV1                  end do
cV1                  y(iy,ir) = factor * sum + y(iy,ir)
cV1               end do
cV1            end do
!!            do iy=0,nz
!!               do ix=0,nz
!!                  acol(ix) = a(iy,ix)
!!               enddo
!!               do ir=0,nr
!!                  sum = zero 
!!                  do ix=0,nz
!!                     sum     = acol(ix) * x(ix,ir) + sum
!!                  end do
!!                  y(iy,ir) = factor * sum + y(iy,ir)
!!               end do
!!            end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            iystart=mod(nz+1,4)
c
            if (iystart.ne.0) then
               do iy=0,iystart-1
                  do ir=0,nr
                     sum1=zero
                     do ix=0,nz
                        sum1 = sum1 + a(iy,ix) * x(ix,ir)
                     end do
                     y(iy,ir) = y(iy,ir) + sum1 * factor
                  end do
               enddo
            endif
c            
            do iy=iystart,nz,4
               do ir=0,nr
                  sum1=zero
                  sum2=zero
                  sum3=zero
                  sum4=zero
                  do ix=0,nz
                     sum1 = sum1 + a(iy,ix) * x(ix,ir)
                     sum2 = sum2 + a(iy+1,ix) * x(ix,ir)
                     sum3 = sum3 + a(iy+2,ix) * x(ix,ir)
                     sum4 = sum4 + a(iy+3,ix) * x(ix,ir)
                  end do
                  y(iy,ir) = y(iy,ir) + sum1 * factor
                  y(iy+1,ir) = y(iy+1,ir) + sum2 * factor
                  y(iy+2,ir) = y(iy+2,ir) + sum3 * factor
                  y(iy+3,ir) = y(iy+3,ir) + sum4 * factor
               end do
            enddo
c
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!schlecht!            call DGEMM('N','N',nzd+1,nrd+1,nzd+1,factor,a(0,0),nzd+1,
!     &                  x(0,0),nzd+1,1.d0,y(0,0),nzd+1) 
         endif
!!c
c     band matrix multiplication
c     ------------------------------------------------------------------
      else
         if (abs(factor).gt.smal) then
cV1            do ir=0,nr
cV1               do iy=0,nz
cV1                  sum = zero
cV1                  ixmin=max(0,iy-iband)
cV1                  ixmax=min(nz,iy+iband)
cV1                  do ix=ixmin,ixmax
cV1                     sum = a(iy,ix) * x(ix,ir) + sum
cV1                  enddo
cV1                  y(iy,ir) = factor * sum + y(iy,ir)
cV1               end do
cV1            end do
            do iy=0,nz
               ixmin=max(0,iy-iband)
               ixmax=min(nz,iy+iband)
               do ix=ixmin,ixmax
                  acol(ix) = a(iy,ix)
               enddo
               do ir=0,nr
                  sum = zero
                  do ix=ixmin,ixmax
                     sum = acol(ix) * x(ix,ir) + sum
                  enddo
                  y(iy,ir) = factor * sum + y(iy,ir)
               end do
            end do
         endif
      endif
c
c-----------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     mulr(a,x,y,factor,iaccu,iband)
c
c     matrix multiplication for r-direction.
c     computes
c
c          factor*(a*x)+ y --> y    for   iaccu=1
c          factor*(a*x)    --> y    for   iaccu=0
c
c     a is assumed to be a bandmatrix with A_{i,i+j}=0 for |j|>iband
c     or a full matrix if iband==-1
c
c=======================================================================
      subroutine mulr(a,x,y,factor,iaccu,iband)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc' ! nrd,nzd
      include 'dynpar.inc' ! nr,nz
c
c     ------------------------------------------------------------------
c
      double precision a(0:nrd,0:nrd),x(0:nzd,0:nrd),y(0:nzd,0:nrd)
      double precision s(0:nzd)
      double precision factor,tmp
      integer          iaccu,iband,irstart
      integer          ir,iy,iz,ixmin,ixmax
c
c     ------------------------------------------------------------------
c
      do iz=0,nz
         s(iz) = 0.0D0
      end do
c
c     set target field to zero in case of no accumulation
c     ------------------------------------------------------------------
      if (iaccu.eq.0) then
        do ir=0,nr
           do iz=0,nz
              y(iz,ir) = zero
           enddo
        enddo
      endif
c
c     full matrix multiplication
c     ------------------------------------------------------------------
      if (iband.eq.-1) then
         if (abs(factor).gt.smal) then
!!            do iy=0,nr
!!               do ix=0,nr
!!                  do iz=0,nz
!!                     s(iz) = a(iy,ix) * x(iz,ix) + s(iz)
!!                  end do
!!               enddo
!!               do iz=0,nz
!!                  y(iz,iy) = factor * s(iz) + y(iz,iy)
!!                  s(iz)    = zero
!!               end do
!!            end do
            irstart=mod(nr+1,4)
            if (irstart.ne.0) then
               do iy=0,irstart-1
                  do ix=0,nr
                     do iz=0,nz
                        y(iz,iy) = y(iz,iy) + factor*x(iz,ix)*a(iy,ix)
                     end do
                  enddo
               end do
            endif
c
            do iy=irstart,nr,4
               do ix=0,nr
                  do iz=0,nz
                     tmp=factor*x(iz,ix)
                     y(iz,iy)   = y(iz,iy)   + tmp * a(iy,ix)
                     y(iz,iy+1) = y(iz,iy+1) + tmp * a(iy+1,ix) 
                     y(iz,iy+2) = y(iz,iy+2) + tmp * a(iy+2,ix) 
                     y(iz,iy+3) = y(iz,iy+3) + tmp * a(iy+3,ix) 
                  end do
               enddo
            end do
c
         endif
c
c     band matrix multiplication
c     ------------------------------------------------------------------
      else
         if (abs(factor).gt.smal) then
            do iy=0,nr
               ixmin=max(0,iy-iband)
               ixmax=min(nr,iy+iband)
               do ix=ixmin,ixmax
                  do iz=0,nz
                     s(iz) = a(iy,ix) * x(iz,ix) + s(iz)
                  enddo
               enddo
               do iz=0,nz
                  y(iz,iy) = factor * s(iz) + y(iz,iy)
                  s(iz)    = zero
               end do
            end do
         endif
      endif
c
c-----------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     mulr0(a,x,y,factor,iaccu,iband)
c
c     matrix multiplication for r-direction at r=0.
c     computes
c
c          factor*(a*x)+ y --> y    for   iaccu=1
c          factor*(a*x)    --> y    for   iaccu=0
c
c     a is assumed to be a bandmatrix with A_{i,i+j}=0 for |j|>iband
c     or a full matrix if iband==-1
c
c=======================================================================
      subroutine mulr0(a,x,y,factor,iaccu,iband)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc' ! nrd,nzd
      include 'dynpar.inc' ! nr,nz
c
c     ------------------------------------------------------------------
c
      double precision a(0:nrd,0:nrd),x(0:nzd,0:nrd),y(0:nzd,0:nrd)
      double precision s(0:nzd)
      double precision factor
      integer          iaccu,iband
      integer          ir,iy,iz,ixmin,ixmax
c
c     ------------------------------------------------------------------
c
c     set target field to zero in case of no accumulation
c     ------------------------------------------------------------------
      if (iaccu.eq.0) then
         do iz=0,nz
            y(iz,0) = zero
         enddo
      endif
c
c     full matrix multiplication
c     ------------------------------------------------------------------
      if (iband.eq.-1) then
         if (abs(factor).gt.smal) then
            do ix=0,nr
               do iz=0,nz
                  s(iz) = a(0,ix) * x(iz,ix) + s(iz)
               end do
            enddo
            do iz=0,nz
               y(iz,0) = factor * s(iz) + y(iz,0)
               s(iz)    = zero
            end do
         endif
c
c     band matrix multiplication
c     ------------------------------------------------------------------
      else
         if (abs(factor).gt.smal) then
            ixmax=min(nr,iband)
            do ix=0,ixmax
               do iz=0,nz
                  s(iz) = a(0,ix) * x(iz,ix) + s(iz)
               enddo
            enddo
            do iz=0,nz
               y(iz,0) = factor * s(iz) + y(iz,0)
               s(iz)    = zero
            end do
         endif
      endif
c
c-----------------------------------------------------------------------
c
      return
      end
