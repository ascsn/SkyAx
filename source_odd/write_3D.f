c-----writ3D-----------------------------------------------------------
c
      subroutine writ3D()
*
***********************************************************************
*                                                                     *
*     print of densities if iter matches print modulus.               *
*     evaluate energies, mass number , deformation, radii, etc.       *
*     prints those quantities if iter matches print modulus.          *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
c      PARAMETER (NXYZ3D)
c      PARAMETER (NX3=NXYZ3D,NY3=KWRITX,NZ3=KWRITX)    ! define 3D grid
*
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
c      include 'derivs.inc'
*
      real*8     rval(4),wfval(4)                 ! for interpolation
      real*8     cmfield(3)
      real*8     rfield3D(NX3,NY3,NZ3)
      complex*16 cfield3D(NX3,NY3,NZ3,2)
      complex*16 circle
      integer    npmin3(2),npsi3(2)                 ! 3D book-keeping
      integer    npmin2(2),npsi2(2)                 ! 2D pointers
      equivalence (rfield3D,cfield3D)
c
c-------------------------------------------------------------------------
c
c     generate output unit on standard name for 3D code
c
      open(55,status='unknown',form='unformatted',
     &        file='restart.data')
c
c     3D code uses reverse ordering of isospin - reshuffle
c
      npmin3(1) = 1
      npsi3(1)  = (nst(2))
      npmin3(2) = 1+npsi3(1)
      npsi3(2)  = nst(1)+nst(2)
      npmin2(1) = nst(1)+1
      npsi2(1)  = nst(1)+nst(2)
      npmin2(2) = 1
      npsi2(2)  = nst(1)
      cmfield(1) = 0.0D0
      cmfield(2) = 0.0D0
      cmfield(3) = icm*dz
      write(6,'(a,5(2i4,2x))')
     & ' nst,npmin3,npsi3,npmin2,npsi2=',nst,npmin3,npsi3,npmin2,npsi2
c
c     write single-particle information , remind doubling of states
c
      DO iq = 1,2  
        write(6,'(a,3i4)') ' iq,npmin,npsi=',iq,npmin3(iq),npsi3(iq)
        write(6,'(a,3i4)') 
     &      ' nst(3-iq)=',nst(3-iq),npsi2(iq)-npmin2(iq)+1,
     &         npsi3(iq)-npmin3(iq)+1
        IF(iq.eq.1) THEN 
          WRITE(55) 1,2*nst(2),cmfield
        ELSE
          WRITE(55) 2*nst(2)+1,2*(nst(1)+nst(2)),cmfield
        ENDIF
        WRITE(55) (1,ist=1,2*nst(3-iq))
        WRITE(55) ((vsqf(ist),is=1,2),ist=npmin2(iq),npsi2(iq))
        WRITE(55) (0.0D0,ist=1,2*nst(3-iq))
        WRITE(55) ((spenrg(ist),is=1,2),ist=npmin2(iq),npsi2(iq))
        WRITE(55) ((1.0D0*npar(ist),is=1,2),ist=npmin2(iq),npsi2(iq))
        WRITE(55) ((spn(ist),is=1,2),ist=npmin2(iq),npsi2(iq))
        WRITE(55) ((0.0D0,is=1,2),ist=npmin2(iq),npsi2(iq))
        WRITE(55) ((spevar(ist),is=1,2),ist=npmin2(iq),npsi2(iq))
      END DO
c
c      write three times a dummy array: 
c      2 times for Yukawa, once for Coulomb
c
       WRITE(55) rfield3D
       WRITE(55) rfield3D
       WRITE(55) rfield3D
c
c      interpolation of wavefunctions on 3D grid 
c      (map spin-up and spin-down partner)
c
      DO iq = 1,2  
      DO ist=npmin2(iq),npsi2(iq)
        DO is=1,2
        mmact = mpsi(ist)+is-1
        DO iz=1,NZ3
        DO iy=1,NY3
        DO ix=1,NX3
          x1 = ix*dr-(NX3+1)*dr*0.5D0
          y1 = iy*dr-(NY3+1)*dr*0.5D0
          r1 = sqrt(x1*x1+y1*y1)
          angle = ATAN(y1/x1)
          IF(y1.LT.0.0D0) angle = angle+PI
          ir  = r1/dr
          DO ival=1,4
            iract = ir+ival-2
            rval(ival) = iract*dr+sr
            IF(iract.GE.0) THEN
              wfval(ival) = psi(iz,iract,is,ist)
            ELSE
              wfval(ival) = (2*mod(is,2)-1)*nrpar(ist)
     &                     *psi(iz,1,is,ist)
            ENDIF
          END DO          
          call polint(rval,wfval,4,r1,wfact,error)
          r1 = max(r1,smal)
          circle = CMPLX(x1/r1,y1/r1)
          cfield3D(ix,iy,iz,is) =  wfact*circle**mmact
c     &         wfact*cexp(CMPLX(0.0D0,mmact*angle))
c     &         wfact*cexp(CMPLX(0.0D0,-mmact*angle))
        END DO
        END DO
        END DO
        END DO
        WRITE(55) (((cfield3D(ix,iy,iz,1),ix=1,NX3),
     &                iy=1,NY3),iz=1,NZ3)
        WRITE(55) (((cfield3D(ix,iy,iz,2),ix=1,NX3),
     &                iy=1,NY3),iz=1,NZ3)
        WRITE(55) (((CONJG(cfield3D(ix,iy,iz,2)),ix=1,NX3),
     &                iy=1,NY3),iz=1,NZ3)
        WRITE(55) (((CONJG(cfield3D(ix,iy,iz,1)),ix=1,NX3),
     &                iy=1,NY3),iz=1,NZ3)
      END DO
      END DO
      CLOSE(55)
      write(6,'(a)') ' 3D restart file written'
c
c     all done
c
      RETURN
      END       
c-----------------------------------------------------------------------
c
      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=abs(x-xa(1))
      do 11 i=1,n
        dift=abs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.d0)pause 'failure in polint'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
