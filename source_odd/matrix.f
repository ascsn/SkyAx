c=======================================================================
c
c     matrml(a,b,c,imax,nmax)
c
c     full matrix multiplication a*b=c. If b=c the matrix b is 
c     overwritten.
c
c       a,b,c  :  matrices
c       nmax   :  static dimension of matrix field
c       imax   :  dynamic dimesnion of matrix field  
c
c=======================================================================
      subroutine matrml(a,b,c,imax,nmax)
c
c     ------------------------------------------------------------------
c
      double precision a(0:nmax,0:nmax)
      double precision b(0:nmax,0:nmax)
      double precision c(0:nmax,0:nmax) 
      double precision work(0:nmax)
      double precision sum
      integer          i,j,l
      integer          imax,nmax
c
c     ------------------------------------------------------------------
c
c     check dimensions
c     ------------------------------------------------------------------
      if(imax.gt.nmax) stop ' dimension > "nmax" in MATRML'
c
c     matrix multiplication
c     ------------------------------------------------------------------
      do j=0,imax
         do i=0,imax
            sum = 0.0
            do l=0,imax
               sum = sum + a(l,i)*b(j,l)
            enddo
            work(i) = sum
         enddo
         do i=0,imax
            c(j,i) = work(i)
         enddo
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end

c=======================================================================
c
c     matrml2(a,b,c,imax,nmax)
c
c     full matrix multiplication a*b=c. If b=c the matrix b is 
c     overwritten.
c
c       a,b,c  :  matrices
c       nmax   :  static dimension of matrix field
c       imax   :  dynamic dimesnion of matrix field  
c
c     same subroutine as matrml but here starts the sum with 1
c
c=======================================================================
      subroutine matrml2(a,b,c,imax,nmax)
c
c     ------------------------------------------------------------------
c
      double precision a(1:nmax,1:nmax)
      double precision b(1:nmax,1:nmax)
      double precision c(1:nmax,1:nmax) 
      double precision work(1:nmax)
      double precision sum
      integer          i,j,l
      integer          imax,nmax
c
c     ------------------------------------------------------------------
c
c     check dimensions
c     ------------------------------------------------------------------
      if(imax.gt.nmax) stop ' dimension > "nmax" in MATRML'
c
c     matrix multiplication
c     ------------------------------------------------------------------
      do j=1,imax
         do i=1,imax
            sum = 0.0
            do l=1,imax
               sum = sum + a(l,i)*b(j,l)
            enddo
            work(i) = sum
         enddo
         do i=1,imax
            c(j,i) = work(i)
         enddo
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end

c=======================================================================
c
c     minv(a,imax,nmax)
c
c     matrix inversion of the matrix a
c
c       a      :  matrix to invert 
c       work   :  workspace
c       nmax   :  static dimension of matrix field
c       imax   :  dynamic dimension of matrix field  
c
c=======================================================================
      subroutine minv(a,work,nzrd,nzr)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
c
c     ------------------------------------------------------------------
c
      double precision  a(nzrd,nzrd),work(nzrd,nzrd)
c      dimension a(nzrd,nzrd),work(nzrd,nzrd)
      integer   index(nmx)
c      integer   nzr,nzrd,index(nzr)
c
c     ------------------------------------------------------------------
c     
c     initialize unity matrix
c     ------------------------------------------------------------------
      do i=1,nzr
        do j=1,nzr
          work(i,j) = 0.0
        enddo
        work(i,i) = 1.0
      enddo
c
c     lower-upper decomposition
c     ------------------------------------------------------------------
      call ludcmp(a,nzr,nzrd,index,sign)
c
c     backward substitution
c     ------------------------------------------------------------------
      do i=1,nzr
        call lubksb(a,nzr,nzrd,index,work(1,i))
      enddo
c
c     copy inverse matrix back to original field
c     ------------------------------------------------------------------
      do i=1,nzr
        do j=1,nzr
          a(i,j) = work(i,j)
        enddo
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end




c=======================================================================
c
c     lubksb(...)
c
c     backward substitution (needed for matrix inversion)
c
c=======================================================================
      SUBROUTINE lubksb(a,n,np,indx,b)
c
c     ------------------------------------------------------------------
c
      INTEGER n,np,indx(n)
      DOUBLE PRECISION a(np,np),b(n)
      INTEGER i,ii,j,ll
      DOUBLE PRECISION sum
c
c     ------------------------------------------------------------------
c
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.d0) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
c
c     ------------------------------------------------------------------
c
      return
      END


c=======================================================================
c
c     ludcmp(...)
c
c     backward substitution (needed for matrix inversion)
c
c=======================================================================c
      SUBROUTINE ludcmp(a,n,np,indx,d)
c
c     ------------------------------------------------------------------
c
      INTEGER n,np,indx(n),NMAX
      DOUBLE PRECISION d,a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0d-20)
      INTEGER i,imax,j,k
      DOUBLE PRECISION aamax,dum,sum,vv(NMAX)
      data imax/0/
c
c     ------------------------------------------------------------------
c

 



      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.d0) pause 'singular matrix in LUDCMP'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.d0)a(j,j)=TINY
        if(j.ne.n)then
          dum=1.d0/a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
c
c     ------------------------------------------------------------------
c      
      return
      END






c=======================================================================
c
c     choldc(...)
c
c     needed for matrix inversion and calculation of a determinant
c
c     works only for symmetric and pos. definit matrices
c
c=======================================================================c
      SUBROUTINE choldc(a,n,np,p)
c
c     ------------------------------------------------------------------
c


      INTEGER n,np
      DOUBLE PRECISION p(n),a(np,np)
      INTEGER i,j,k
      DOUBLE PRECISION sum

c
c     ------------------------------------------------------------------
c




      do 13 i=1,n
      do 12 j=i,n
         sum=a(i,j)
      do 11 k=i-1,1,-1
         sum=sum-a(i,k)*a(j,k)
 11   continue
      if(i.eq.j)then
         if(sum.le.0.)pause 'choldc failed' !a, with rounding errors, is not
         p(i)=sqrt(sum)                     !positive definite.
      else
         a(j,i)=sum/p(i)
      endif

 12   continue
 13   continue



c
c     ------------------------------------------------------------------
c      
      return
      END




c=======================================================================
c
c     cholsl(...)
c
c     backward substitution
c
c     needed for matrix inversion and calculation of a determinant
c
c     works only for symmetric and pos. definit matrices
c
c=======================================================================c
      SUBROUTINE cholsl(a,n,np,p,b,x)
c
c     ------------------------------------------------------------------
c


      INTEGER n,np
      DOUBLE PRECISION p(n),a(np,np),x(n),b(n)
      INTEGER i,k
      DOUBLE PRECISION sum

c
c     ------------------------------------------------------------------
c

      do 12 i=1,n 
         sum=b(i)
         do 11 k=i-1,1,-1
            sum=sum-a(i,k)*x(k)

 11      continue 
         x(i)=sum/p(i)
 12   continue 
      do 14 i=n,1,-1 
         sum=x(i)
         do 13 k=i+1,n
            sum=sum-a(k,i)*x(k)
 13      continue 
         x(i)=sum/p(i)
 14   continue







c
c     ------------------------------------------------------------------
c      
      return
      END





c=======================================================================
c
c     minv2(a,imax,nmax)
c
c     matrix inversion of the matrix a
c
c       a      :  matrix to invert 
c       work   :  workspace
c       nmax   :  static dimension of matrix field
c       imax   :  dynamic dimension of matrix field  
c
c=======================================================================
      subroutine minv2(a,work,nzrd,nzr)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
c
c     ------------------------------------------------------------------
c
      double precision a(nzrd,nzrd),work(nzrd,nzrd),p(nzr)
      integer   index(nmx)
c
c     ------------------------------------------------------------------
c     
c     initialize unity matrix
c     ------------------------------------------------------------------
      do i=1,nzr
        do j=1,nzr
          work(i,j) = 0.0
        enddo
        work(i,i) = 1.0
      enddo
c
c     lower-upper decomposition
c     ------------------------------------------------------------------

      call choldc(a,nzr,nzrd,p)



c      call ludcmp(a,nzr,nzrd,index,sign)
c
c     backward substitution
c     ------------------------------------------------------------------
      do i=1,nzr

        call cholsl(a,nzr,nzrd,p,work(1,i),work(1,i))

c        call lubksb(a,nzr,nzrd,index,work(1,i))
      enddo
c
c     copy inverse matrix back to original field
c     ------------------------------------------------------------------
      do i=1,nzr
        do j=1,nzr
          a(i,j) = work(i,j)
        enddo
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end
