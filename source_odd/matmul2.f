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

      IF(iband.NE.-1) STOP 'MULZ: finite bandwidth not implemented'

      beta=1D0
      IF (iaccu == 0) beta = 0D0

      IF (ABS(factor) > smal) 
     &  CALL DGEMM('N','N',nz+1,nr+1,nz+1,factor,a,nzd+1,x,nzd+1,beta,y,
     &             nzd+1) 

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

      IF(iband.NE.-1) STOP 'MULR: finite bandwidth not implemented'
      beta=1D0
      IF (iaccu == 0) beta = 0D0

      IF (ABS(factor) > smal) CALL DGEMM('N','T',nz+1,nr+1,nr+1,
     &                       factor,x,nzd+1,a,nrd+1,beta,y,nzd+1) 

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
      IF(iband.NE.-1) STOP 'MULR0: finite bandwidth not implemented'
      beta=1D0
      IF (iaccu == 0) beta = 0D0

      IF (ABS(factor) > smal) CALL DGEMV('N',nz+1,nr+1,factor,
     &       x,nzd+1,a,nzd+1,beta,y,1)


c
c-----------------------------------------------------------------------
c
      return
      end
