c======================================================================
c
c     tst_grid()
c
c=======================================================================
      subroutine tst_grid()
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc' ! nzd,nrd,kmpsi
      include 'derivs.inc' ! derivative matrices, iband_derz,iband_derr
      include 'dynpar.inc' ! isym,dr,nr,sr,nz,dz
c
c     ------------------------------------------------------------------
c
      double precision fac, sumnor, sumexp, sumvar, weight
      double precision cosi,cosj,sini,sinj
      integer          ikz,nkz,iz,jz,jz2
c
c     ------------------------------------------------------------------
c     set dimensions of grid in fourier space     
c     ------------------------------------------------------------------
c
      if (isym.eq.0) then
        nkz  = (nz+1)/2
        nkr  =  nr
        ip1  = 1
        ip2  = 1
      else
        nkz  =  nz
        nkr  =  nr
        ip1  = 1
        ip2  = 2
      endif
c
c     ------------------------------------------------------------------
c     open the test file                         
c     ------------------------------------------------------------------
c
      open(85,status='unknown',form='formatted',file='test.grid')
c
c     ------------------------------------------------------------------
c     test 1st derivatives in z
c     ------------------------------------------------------------------
c
      write(85,'(a)') '# test of 1st derivatives (A) in z (even par.)'
      write(85,'(a)') '#      k       -i<A>      <A>-A_exact     Var(A)'
c     
      do ikz=0,nkz
c
         fac    = pi/nz*ikz    ! k*dz
         sumnor = 0.0D0
         sumexp = 0.0D0
         sumvar = 0.0D0
c
         do iz=0,nz
c           ! summation of the norm
            cosi   = cos(iz*fac)
            if (iz.eq.0 .or. iz.eq.nz) then
               sumnor = 0.5D0*cosi*cosi + sumnor
            else 
               sumnor = cosi*cosi + sumnor
            endif
c
            do jz=0,nz
c              ! summation of the expectation
               cosj = cos(jz*fac)
               sinj = sin(jz*fac)
               if (jz.eq.0 .or. jz.eq.nz) then
                  sumexp = 0.5D0*sinj*derz1(jz,iz,ip1)*cosi + sumexp
               else 
                  sumexp = sinj*derz1(jz,iz,ip1)*cosi + sumexp
               endif
c
               do jz2=0,nz
c                 ! summation of the variance
                  if (jz.eq.0 .or. jz.eq.nz) then
                     sumvar = 0.5D0*cosj*derz1(jz,jz2,ip2)*
     &                              derz1(jz2,iz,ip1)*cosi + sumvar
                  else 
                     sumvar = cosj*derz1(jz,jz2,ip2)*
     &                              derz1(jz2,iz,ip1)*cosi + sumvar
                  endif
c
               enddo
            enddo
         enddo
c
         write(85,'(g12.4,g12.4,g12.4,g12.4)') 
     &            fac/dz, sumexp/sumnor, sumexp/sumnor+(fac/dz),
     &            sumvar/sumnor+(sumexp/sumnor)**2
      enddo
c
      write(85,'(a)') ''
      write(85,'(a)') '# test of 1st derivatives (A) in z (odd par.)'
      write(85,'(a)') '#      k       -i<A>      <A>-A_exact     Var(A)'
c     
      do ikz=0,nkz
c
         fac    = pi/nz*ikz    ! k*dz
         sumnor = 0.0D0
         sumexp = 0.0D0
         sumvar = 0.0D0
c
         do iz=0,nz
c           ! summation of the norm
            sini   = sin(iz*fac)
            if (iz.eq.0 .or. iz.eq.nz) then
               sumnor = 0.5D0*sini*sini + sumnor
            else 
               sumnor = sini*sini + sumnor
            endif
c
            do jz=0,nz
c              ! summation of the expectation
               cosj = cos(jz*fac)
               sinj = sin(jz*fac)
               if (jz.eq.0 .or. jz.eq.nz) then
                  sumexp = 0.5D0*cosj*derz1(jz,iz,ip2)*sini + sumexp
               else 
                  sumexp = cosj*derz1(jz,iz,ip2)*sini + sumexp
               endif
c
               do jz2=0,nz
c                 ! summation of the variance
                  if (jz.eq.0 .or. jz.eq.nz) then
                     sumvar = 0.5D0*sinj*derz1(jz,jz2,ip1)*
     &                              derz1(jz2,iz,ip2)*sini + sumvar
                  else 
                     sumvar = sinj*derz1(jz,jz2,ip1)*
     &                              derz1(jz2,iz,ip2)*sini + sumvar
                  endif
c
               enddo
            enddo
         enddo
c
         if (fac/dz.ne.0.0D0) then
            write(85,'(g12.4,g12.4,g12.4,g12.4)') 
     &           fac/dz, sumexp/sumnor, sumexp/sumnor-(fac/dz),
     &           sumvar/sumnor+(sumexp/sumnor)**2
         endif
      enddo
c
c     ------------------------------------------------------------------
c     test 1st derivatives in r
c     ------------------------------------------------------------------
c
      write(85,'(a)') ''
      write(85,'(a)') '# test of 1st derivatives (A) in r (even par.)'
      write(85,'(a)') '#      k       -i<A>      <A>-A_exact     Var(A)'
c     
      do ikr=0,nkr
c
         fac    = pi/nr*ikr    ! k*dr
         sumnor = 0.0D0
         sumexp = 0.0D0
         sumvar = 0.0D0
c
         do ir=0,nr
c           ! summation of the norm
            cosi   = cos(ir*fac)
            if (ir.eq.0 .or. ir.eq.nr) then
               sumnor = 0.5D0*cosi*cosi + sumnor
            else 
               sumnor = cosi*cosi + sumnor
            endif
c
            do jr=0,nr
c              ! summation of the expectation
               cosj = cos(jr*fac)
               sinj = sin(jr*fac)
               if (jr.eq.0 .or. jr.eq.nr) then
                  sumexp = 0.5D0*sinj*derr1(jr,ir,1)*cosi + sumexp
               else 
                  sumexp = sinj*derr1(jr,ir,1)*cosi + sumexp
               endif
c
               do jr2=0,nr
c                 ! summation of the variance
                  if (jr.eq.0 .or. jr.eq.nr) then
                     sumvar = 0.5D0*cosj*derr1(jr,jr2,2)*
     &                              derr1(jr2,ir,1)*cosi + sumvar
                  else 
                     sumvar = cosj*derr1(jr,jr2,2)*
     &                              derr1(jr2,ir,1)*cosi + sumvar
                  endif
c
               enddo
            enddo
         enddo
c
         write(85,'(g12.4,g12.4,g12.4,g12.4)') 
     &            fac/dr, sumexp/sumnor, sumexp/sumnor+(fac/dr),
     &            sumvar/sumnor+(sumexp/sumnor)**2
      enddo
c
      write(85,'(a)') ''
      write(85,'(a)') '# test of 1st derivatives (A) in r (odd par.)'
      write(85,'(a)') '#      k       -i<A>      <A>-A_exact     Var(A)'
c     
      do ikr=0,nkr
c
         fac    = pi/nr*ikr    ! k*dr
         sumnor = 0.0D0
         sumexp = 0.0D0
         sumvar = 0.0D0
c
         do ir=0,nr
c           ! summation of the norm
            sini   = sin(ir*fac)
            if (ir.eq.0 .or. ir.eq.nr) then
               sumnor = 0.5D0*sini*sini + sumnor
            else 
               sumnor = sini*sini + sumnor
            endif
c
            do jr=0,nr
c              ! summation of the expectation
               cosj = cos(jr*fac)
               sinj = sin(jr*fac)
               if (jr.eq.0 .or. jr.eq.nr) then
                  sumexp = 0.5D0*cosj*derr1(jr,ir,2)*sini + sumexp
               else 
                  sumexp = cosj*derr1(jr,ir,2)*sini + sumexp
               endif
c
               do jr2=0,nr
c                 ! summation of the variance
                  if (jr.eq.0 .or. jr.eq.nr) then
                     sumvar = 0.5D0*sinj*derr1(jr,jr2,1)*
     &                              derr1(jr2,ir,2)*sini + sumvar
                  else 
                     sumvar = sinj*derr1(jr,jr2,1)*
     &                              derr1(jr2,ir,2)*sini + sumvar
                  endif
c
               enddo
            enddo
         enddo
c
         if (fac/dr.ne.0.0D0) then
            write(85,'(g12.4,g12.4,g12.4,g12.4)') 
     &           fac/dr, sumexp/sumnor, sumexp/sumnor-(fac/dr),
     &           sumvar/sumnor+(sumexp/sumnor)**2
         endif
      enddo
c
c     ------------------------------------------------------------------
c     test of 2nd derivative in z-direction
c     ------------------------------------------------------------------
c
      write(85,'(a)') ''
      write(85,'(a)') '# test of 2nd derivatives (A) in z (even par.)'
      write(85,'(a)') '#      k         <A>      <A>-A_exact     Var(A)'
c     
      do ikz=0,nkz
c
         fac    = pi/nz*ikz    ! k*dz
         sumnor = 0.0D0
         sumexp = 0.0D0
         sumvar = 0.0D0
c
         do iz=0,nz
c           ! summation of the norm
            cosi   = cos(iz*fac)
            if (iz.eq.0 .or. iz.eq.nz) then
               sumnor = 0.5D0*cosi*cosi + sumnor
            else 
               sumnor = cosi*cosi + sumnor
            endif
c
            do jz=0,nz
c              ! summation of the expectation
               cosj = cos(jz*fac)
               if (jz.eq.0 .or. jz.eq.nz) then
                  sumexp = 0.5D0*cosj*derz2(jz,iz,ip1)*cosi + sumexp
               else 
                  sumexp = cosj*derz2(jz,iz,ip1)*cosi + sumexp
               endif
c
               do jz2=0,nz
c                 ! summation of the variance
                  if (jz.eq.0 .or. jz.eq.nz) then
                     sumvar = 0.5D0*cosj*derz2(jz,jz2,ip1)*
     &                              derz2(jz2,iz,ip1)*cosi + sumvar
                  else 
                     sumvar = cosj*derz2(jz,jz2,ip1)*
     &                              derz2(jz2,iz,ip1)*cosi + sumvar
                  endif
c
               enddo
            enddo
         enddo
c
         write(85,'(g12.4,g12.4,g12.4,g12.4)') 
     &            fac/dz, sumexp/sumnor, sumexp/sumnor+(fac/dz)**2,
     &            sumvar/sumnor-(sumexp/sumnor)**2
      enddo
c
      write(85,'(a)') ''
      write(85,'(a)') '# test of 2nd derivatives (A) in z (odd par.)'
      write(85,'(a)') '#      k         <A>      <A>-A_exact     Var(A)'
c     
      do ikz=0,nkz
c
         fac    = pi/nz*ikz    ! k*dz
         sumnor = 0.0D0
         sumexp = 0.0D0
         sumvar = 0.0D0
c
         do iz=0,nz
c           ! summation of the norm
            sini   = sin(iz*fac)
            if (iz.eq.0 .or. iz.eq.nz) then
               sumnor = 0.5D0*sini*sini + sumnor
            else 
               sumnor = sini*sini + sumnor
            endif
c
            do jz=0,nz
c              ! summation of the expectation
               sinj = sin(jz*fac)
               if (jz.eq.0 .or. jz.eq.nz) then
                  sumexp = 0.5D0*sinj*derz2(jz,iz,ip2)*sini + sumexp
               else 
                  sumexp = sinj*derz2(jz,iz,ip2)*sini + sumexp
               endif
c
               do jz2=0,nz
c                 ! summation of the variance
                  if (jz.eq.0 .or. jz.eq.nz) then
                     sumvar = 0.5D0*sinj*derz2(jz,jz2,ip2)*
     &                              derz2(jz2,iz,ip2)*sini + sumvar
                  else 
                     sumvar = sinj*derz2(jz,jz2,ip2)*
     &                              derz2(jz2,iz,ip2)*sini + sumvar
                  endif
c
               enddo
            enddo
         enddo
c
         if (fac/dz.ne.0.0D0) then
            write(85,'(g12.4,g12.4,g12.4,g12.4)') 
     &           fac/dz, sumexp/sumnor, sumexp/sumnor+(fac/dz)**2,
     &           sumvar/sumnor-(sumexp/sumnor)**2
         endif
      enddo
c
c     ------------------------------------------------------------------
c     test spectral properties of radial laplace
c     ------------------------------------------------------------------
c
      write(85,'(a)') ''
      write(85,'(a)') '# test of radial laplace (even par.)'
      write(85,'(a)') '#      k       <A>      <A>-A_exact     Var(A)'
c     
      write(85,'(a)') ' WARNING: test routines for radial laplace'//
     &                ' not yet implemented'
c
      close(85)
c
c     ------------------------------------------------------------------
c
      return
      end
