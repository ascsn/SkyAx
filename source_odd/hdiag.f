c-----collmatr-----------------------------------------------------------
c
      subroutine collmatr(hmatr,hpsi,nlo,jst,ndimact)
c
c     collects Hamiltonian matrix for block
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'crankin.inc'
c
      dimension hmatr(ndimact,ndimact)
c      dimension psi(0:nzd,0:nrd,nsd,npsid)
      dimension hpsi(0:nzd,0:nrd,nsd)
c
c----------------------------------------------------------------------
c
      do nhst=1,ndimact
        nact = jst-nlo+1
        kst  = nlo+nhst-1
c        write(6,'(a,6i5)')
c     &   ' collmatr: jst,nlo,ndimact,nhst,nact,kst=',
c     &               jst,nlo,ndimact,nhst,nact,kst
c        write(7,'(a,6i5)')
c     &   ' collmatr: jst,nlo,ndimact,nhst,nact,kst=',
c     &               jst,nlo,ndimact,nhst,nact,kst
        if (nparz(1,jst).eq.nparz(1,kst)) then
           hmatr(nhst,nact) = wfovlp(dpsir(0,0,1,kst),hpsi(0,0,1))
        else
           hmatr(nhst,nact) = 0.0D0
        endif
        if(nhst.eq.nact)
     &   hmatr(nhst,nact)= hmatr(nhst,nact)+spenrg(jst) 
c        hmatr(nhst,nact) = wfovlp(psi(0,0,1,nst),psi(0,0,1,nst))
      enddo
c
      return
      end

c-----diagmatr-----------------------------------------------------------
c
      subroutine diagmatr(hmatr,eigvec,nlo,ndimact)
c
c     collects Hamiltonian matrix for block
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'crankin.inc'
c
      dimension hmatr(ndimact,ndimact)
      dimension eigvec(ndimact,ndimact)
      dimension workv1(ndimact,ndimact)
      dimension workv2(ndimact,ndimact)
      dimension root(ndimact)
c
      parameter (itdiagmax=10)  ! 200)
      parameter (diagstep=0.01)
      parameter (diagprecis=1D-6)
      LOGICAL, PARAMETER :: tprint=.FALSE.
      LOGICAL, PARAMETER :: tgivens=.FALSE.

c
c----------------------------------------------------------------------
c
      if(ndimact.eq.1) return
c
      IF(tgivens) THEN
        CALL diag_real_symm(hmatr,root,eigvec,ndimact,ndimact,
     &                      (ndimact*(ndimact+1))/2)
        IF(tprint) WRITE(6,'(a,i3,200(1pg13.5))') 
     &         'hdiag: ndim,roots=',ndimact,root
      ELSE
      if(tprint) then
        write(6,'(a)') 'H matrix:'
        do na=1,ndimact
          write(6,'(10(1pg13.5))')
     &        (hmatr(na,nb),nb=1,ndimact)
        enddo
      endif
c
c     prepare stepping matrix and eigenvectors
c
      hmin = 1D99
      do na=1,ndimact
        hmin = min(hmatr(na,na),hmin)
      enddo
c      
      do na=1,ndimact
      do nb=1,ndimact
        hmatr(na,nb) = -diagstep*hmatr(na,nb)
        eigvec(na,nb) = 0.0D0
        if(na.eq.nb) then
          hmatr(na,nb)=1.0D0+hmatr(na,nb)+diagstep*hmin
          eigvec(na,nb)=1.0D0
        endif
      enddo
      enddo
c
c     iterate
c
      do itdiag=1,itdiagmax
c
c       analysis of variance
c
        do na=1,ndimact
        do nb=1,ndimact
c
          sum = 0.0D0
          do nc=1,ndimact
            sum = hmatr(na,nc)*eigvec(nc,nb) + sum
          enddo
          workv1(na,nb) = sum
c
        enddo
        enddo
c
        do na=1,ndimact
        do nb=1,ndimact
c
          sum = 0.0D0
          do nc=1,ndimact
            sum = hmatr(na,nc)*workv1(nc,nb) + sum
          enddo
          workv2(na,nb) = sum
c
        enddo
        enddo
c
        sumvar = 0.0D0
        do na=1,ndimact
          sumeig = 0.0D0
          sumeig2 = 0.0D0
          do nc=1,ndimact
            sumeig = eigvec(nc,na)*workv1(nc,na) 
     &              + sumeig
            sumeig2 = eigvec(nc,na)*workv2(nc,na)
     &               + sumeig2
          enddo
          variance = sumeig2-sumeig**2
          sum = variance + sum
        enddo
        variance = sqrt(max(abs(variance),smal))
        if(tprint) 
     &    write(6,'(a,2i5,1pg12.4)')
     &    'diag: block,iter,variance:',ndimact,itdiag,
     &    variance
        FLUSH(6)
        if(variance.lt.diagprecis) goto 19
c
c       stepping
c
        do na=1,ndimact
        do nb=1,ndimact
          eigvec(na,nb) = workv1(na,nb)
        enddo
        enddo
        call orthmatr(eigvec,ndimact)
c
c
      enddo
 19     continue
      ENDIF

c
c     transform wavefunctions
c
      do na=1,ndimact
      ja = nlo+na-1
        do is=1,nsd
        do ir=0,nr
        do iz=0,nz
          dpsir(iz,ir,is,ja) = 0.0D0
        enddo
        enddo
        enddo
        do nb=1,ndimact
        jb = nlo+nb-1
c
          do is=1,nsd
          do ir=0,nr
          do iz=0,nz
            dpsir(iz,ir,is,ja) = dpsir(iz,ir,is,ja)
     &         + psi(iz,ir,is,jb)*eigvec(nb,na)
          enddo
          enddo
          enddo
c
        enddo
      enddo
c
ctest      return
c
      do na=1,ndimact
      ja = nlo+na-1
        do is=1,nsd
        do ir=0,nr
        do iz=0,nz
          psi(iz,ir,is,ja) = dpsir(iz,ir,is,ja) 
        enddo
        enddo
        enddo
      enddo
c
      return
      end

c-----orthmatr-----------------------------------------------------------
c
      subroutine orthmatr(eigvec,ndimact)
c
c     ortho-normalizes matrix 'eigvec'
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
c
      dimension eigvec(ndimact,ndimact)
c
      logical ttest
      parameter (ttest=.false.)
c
c----------------------------------------------------------------------
c
c             
      do na=1,ndimact
        do nb=1,na-1
          ov = 0.0D0
          do nc=1,ndimact
            ov = ov + eigvec(nc,nb)*eigvec(nc,na)
          enddo
          do nc=1,ndimact
            eigvec(nc,na) = eigvec(nc,na)-ov*eigvec(nc,nb)
          enddo
        enddo
        ov = 0.0D0
        do nc=1,ndimact
          ov = ov + eigvec(nc,na)*eigvec(nc,na)
        enddo
        ov = 1.0D0/sqrt(ov)
        do nc=1,ndimact
          eigvec(nc,na) = ov*eigvec(nc,na)
        enddo
      enddo
c
c     test
c
      if(.not.ttest) return
c
      write(6,'(a)') ' test orthonormalization:'
      do na=1,ndimact
      do nb=1,ndimact
        ov = 0.0D0
        do nc=1,ndimact
          ov = ov + eigvec(nc,nb)*eigvec(nc,na)
        enddo
        write(6,'(2i5,1pg13.5)') na,nb,ov
      enddo
      enddo
     
c
      return
      end


c-----testmatr-----------------------------------------------------------
c
      subroutine testmatr(hmatr,psit,nlo,ndimact)
c
c     collects Hamiltonian matrix for block
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'crankin.inc'
c
      dimension hmatr(ndimact,ndimact)
      dimension psit(0:nzd,0:nrd,nsd,npsid)
c
c----------------------------------------------------------------------
c
      do na=1,ndimact
        jst = nlo+na-1
        call hpsi(psit(0,0,1,jst),psin,jst)
        do nb=1,ndimact
          kst  = nlo+nb-1
          hmatr(nb,na) = wfovlp(psit(0,0,1,kst),psin(0,0,1))
          if(na.eq.nb)
     &      hmatr(nb,na)= hmatr(nb,na)+spenrg(jst) 
        enddo
      enddo      

      write(6,'(a)') 'H matrix:'
      do na=1,ndimact
        write(6,'(10(1pg13.5))')
     &        (hmatr(na,nb),nb=1,ndimact)
      enddo
c
      return
      end
c-----testdiag-----------------------------------------------------------
c
      subroutine testdiag(hmatr,nlo,ndimact)
c
c     Tests diagonalization of mean-field Hamiltonian
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'crankin.inc'
c
      dimension hmatr(ndimact,ndimact)
      dimension eigvec(ndimact,ndimact)
c
c----------------------------------------------------------------------
c
c     copy wf to workspace
c
      nhi = nlo+ndimact-1
      do jst=nlo,nhi
      do is=1,nsd
      do ir=0,nr
      do iz=0,nz
        dpsir(iz,ir,is,jst) = psi(iz,ir,is,jst) 
      enddo
      enddo
      enddo
      enddo
c
c     compute matrix and diagonalize
c
      write(6,'(a)') 'before diagonalization'
      call  testmatr(hmatr,dpsir,nlo,ndimact)
      call diagmatr(hmatr,eigvec,nlo,ndimact)
c      write(6,'(a)') ' eigvec:'
c      do na=1,ndimact
c        write(6,'(100(1pg13.5))')
c     &   (eigvec(nb,na),nb=1,ndimact)
c      enddo
c
c     test matrix
c
      write(6,'(a)') 'after diagonalization'
c      call  testoverlp(hmatr,dpsir,nlo,ndimact)
c      call  testoverlp(hmatr,psi,nlo,ndimact)
      call  testmatr(hmatr,dpsir,nlo,ndimact)
c
      return
      end

c-----testmatr-----------------------------------------------------------
c
      subroutine testoverlp(omatr,psit,nlo,ndimact)
c
c     collects Hamiltonian matrix for block
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
c      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'crankin.inc'
c
      dimension omatr(ndimact,ndimact)
      dimension psit(0:nzd,0:nrd,nsd,npsid)
c
c----------------------------------------------------------------------
c
      do na=1,ndimact
        jst = nlo+na-1
        do nb=1,ndimact
          kst  = nlo+nb-1
          omatr(nb,na) = wfovlp(psit(0,0,1,kst),psit(0,0,1,jst))
        enddo
      enddo
      
      write(6,'(a)') 'overlap matrix:'
      do na=1,ndimact
        write(6,'(10(1pg13.5))')
     &        (omatr(na,nb),nb=1,ndimact)
      enddo
c
      return
      end
c-----diag_real_symm------------------------------------------------
c
      subroutine diag_real_symm(rmatrix,root,vect,
     &                          nmax,ndim,ndimlin)
c
c     Diagonalaizes a real symmetric matrix.     
c     List parameters are:
c      rmatrix   = input matrix in quadratic storage
c      root      = linear array of eigenvalues
c      vect      = rectangular array of eigenvectors
c      nmax      = actual dimension of the matrix
c                  and numbers of roots
c      ndim      = general dimension of the matrix
c      ndimlin   = nmax*(nmax+1)/2, used for allocation of
c                  workspace
c
      implicit none
c
      integer nmax,ndim,nvect,ndimlin
      double precision rmatrix(ndim,ndim)
      double precision vect(ndim,ndim)
      double precision root(ndim)
c
      integer i,j,ipl
      double precision work(ndimlin)
c
c--------------------------------------------------------------------
c

c    
c     reformat input matrix into upper triangular storage
c
      ipl = 0
      do i=1,nmax
      do j=1,i
        ipl = 1+ipl
        work(ipl) = rmatrix(i,j)
      enddo
      enddo

c
c     diagonalize by call to other routine
c
      call givens(work,root,vect,nmax,nmax,ndim)

      return
      end



c-----givens------------------------------------------------------------------
c
         subroutine givens(a,root,vect,nx,nrootx,njx)
c
c 62.3  givens  -eigenvalues and eigenvectors by the givens method.
c      by franklin prosser, indiana university.
c      september, 1967
c      calculates eigenvalues and eigenvectors of real symmetric matrix
c      stored in packed upper triangular form.
c
c      thanks are due to f. e. harris (stanford university) and h. h.
c      michels (united aircraft research laboratories) for excellent
c      work on numerical difficulties with earlier versions of this
c      program.
c
c      the parameters for the routine are...
c          nx     order of matrix
c          nrootx number of roots wanted.  the nrootx smallest (most
c                  negative) roots will be calculated.  if no vectors
c                  are wanted, make this number negative.
c          njx    row dimension of vect array.  see 'vect' below.
c                  njx must be not less than nx.
c          a      matrix stored by columns in packed upper triangular
c                 form, i.e. occupying nx*(nx+1)/2 consecutive
c                 locations.
c          b      scratch array used by givens.  must be at least
c                  nx*5 cells.
c          root   array to hold the eigenvalues.  must be at least
c                 nrootx cells long.  the nrootx smallest roots are
c                  ordered largest first in this array.
c          vect   eigenvector array.  each column will hold an
c                  eigenvector for the corresponding root.  must be
c                  dimensioned with 'njx' rows and at least 'nrootx'
c                  columns, unless no vectors
c                  are requested (negative nrootx).  in this latter
c                  case, the argument vect is just a dummy, and the
c                  storage is not used.
c
c      the arrays a and b are destroyed by the computation.  the results
c      appear in root and vect.
c
c      for proper functioning of this routine, the result of a floating
c      point underflow should be a zero.
c
c      to convert this routine to double precision (e.g. on ibm 360
c      machines), be sure that all real variables and function
c      references are properly made double precision.
c
c      the original reference to the givens technique is in oak ridge
c      report number ornl 1574 (physics), by wallace givens.
c      the method as presented in this program consists of four steps,
c      all modifications of the original method...
c      first, the input matrix is reduced to tridiagonal form by the
c      householder technique (j. h. wilkinson, comp. j. 3, 23 (1960)).
c      the roots are then located by the sturm sequence method (j. m.
c      ortega (see reference below).  the vectors of the tridiagonal
c      form are then evaluated (j. h. wilkinson, comp. j. 1, 90 (1958)),
c      and last the tridiagonal vectors are rotated to vectors of the
c      original array (first reference).
c      vectors for degenerate (or near-degenerate) roots are forced
c      to be orthogonal, using a method suggested by b. garbow, argonne
c      national ldabs (private communication, 1964).  the gram-schmidt
c      process is used for the orthogonalization.
c
c      an excellent presentation of the givens technique is found in
c      j. m. ortega's article in 'mathematics for digital computers,'
c      volume 2, ed. by ralston and wilf, wiley (1967), page 94.
c
c
c ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c **   users please note...
c **   the following two parameters, eta and theta, should be adjusted
c **   by the user for his particular machine.
c **   eta is an indication of the precision of the floating point
c **   representation on the computer being used (roughly 10**(-m),
c **   where m is the number of decimals of precision ).
c **   theta is an indication of the range of numbers that can be
c **   expressed in the floating point representation (roughly the
c **   largest number or the reciprocal of the smallest number,
c **   whichever is smaller).
c **   some recommended values follow.
c **   for ibm 7094, univac 1108, etc. (27-bit binary fraction, 8-bit
c **   binary exponent), eta=1.e-8, theta=1.e37.
c **   for control data 3600 (36-bit binary fraction, 11-bit binary
c **   exponent), eta=1.e-11, theta=1.e307.
c **   for control data 6600 (48-bit binary fraction, 11-bit binary
c **   exponent), eta=1.e-14, theta=1.e295.
c **   for ibm 360/50 and 360/65 double precision (56-bit hexadecimal
c **   fraction, 7-bit hexadecimal exponent), eta=1.e-16, theta=1.e75.
c **
      implicit double precision (a-h,o-z)
c
!      include 'params.inc'
!      parameter (n1=krpmax)
      parameter (n1=500)
c
      dimension b(n1,5)
c
      dimension a(nx*(nx+1)/2),root(nrootx),vect(njx,nrootx)
c
      data  eta,theta/.1d-14,1.d88/
csing      data  eta,theta/5.0e-7,1.0e33/
c
c ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     check size
c
      if(nx.gt.n1) stop ' givens was called with nx too large '
c
       del1 = eta/100.
       delta = eta**2*100.
       small = eta**2/100.
       delbig = theta*delta/1000.
       theta1 = 1000./theta
c      toler  is a factor used to determine if two roots are close
c      enough to be considered degenerate for purposes of orthogonali-
c      zing their vectors.  for the matrix normed to unity, if the
c      difference between two roots is less than toler, then
c      orthogonalization will occur.
       toler = eta*100.
c
c      initial value for pseudorandom number generator... (2**23)-3
       rpower = 8388608.
       rpow1 = rpower/2.
       rand1 = rpower - 3.
c
      n = nx
      nroot = iabs(nrootx)
      if (nroot.eq.0) go to 1001
      if (n-1) 1001,1003,105
1003  root(1) = a(1)
      if (nrootx.gt.0) vect(1,1) = 1.0
      go to 1001
105   continue
c     nsize    number of elements in the packed array
      nsize = (n*(n+1))/2
      nm1 = n-1
      nm2 = n-2
c
c     scale matrix to euclidean norm of 1.  scale factor is anorm.
      factor = 0.
      do 70 i=1,nsize
70    factor = max(factor,abs(a(i)))
c     if (factor.ne.0.) go to 72
      if(abs(factor).gt.0.) go to 72
c     null matrix.  fix up roots and vectors, then exit.
      do 78 i=1,nroot
      if (nrootx.lt.0) go to 78
      do 77 j=1,n
77    vect(j,i) = 0.
      vect(i,i) = 1.0
78    root(i) = 0.
      go to 1001
c
72    anorm = 0.
      j = 1
      k = 1
86    do 80 i=1,nsize
      if (i.ne.j) go to 81
      anorm = anorm + (a(i)/factor)**2/2.
      k = k+1
      j = j+k
      go to 80
81    anorm = anorm + (a(i)/factor)**2
80    continue
83    anorm = sqrt(anorm*2.)*factor
      do 91 i=1,nsize
91    a(i) = a(i)/anorm
      alimit = 1.0
c
c      tridia section.
c      tridiagonalization of symmetric matrix
       id = 0
       ia = 1
       if (nm2.eq.0) go to 201
       do 200  j=1,nm2
c      j       counts row  of a-matrix to be diagonalized
c      ia      start of non-codiagonal elements in the row
c      id      index of codiagonal element on row being codiagonalized.
       ia = ia+j+2
       id = id + j + 1
       jp2 = j+2
c      sum squares of non-codiagonal elements in row j
       ii = ia
       sum = 0.0
       do 100 i=jp2,n
       sum = sum + a(ii)**2
100    ii = ii + i
       temp = a(id)
       if (sum.gt.small) go to 110
c      no transformation necessary if all the non-codiagonal
c      elements are tiny.
120    b(j,1) = temp
       a(id) = 0.
       go to 200
c      now complete the sum of off-diagonal squares
110    sum = sqrt(sum + temp**2)
c      new codiagonal element
       b(j,1) = -sign(sum,temp)
c      first non-zero element of this w-vector
       b(j+1,2) = sqrt((1.0 + abs(temp)/sum)/2.0)
c      form rest of the w-vector elements
       temp = sign(0.5/(b(j+1,2)*sum),temp)
       ii = ia
       do 130 i=jp2,n
       b(i,2) = a(ii)*temp
130    ii = ii + i
c      form p-vector and scalar.  p-vector = a-matrix*w-vector.
c      scalar = w-vector*p-vector.
       ak = 0.0
c      ic      location of next diagonal element
       ic = id + 1
       j1 = j + 1
       do 190  i=j1,n
       jj = ic
       temp = 0.
       do 180  ii=j1,n
c      i       runs over the non-zero p-elements
c      ii      runs over elements of w-vector
       temp = temp + b(ii,2)*a(jj)
c      change incrementing mode at the diagonal elements.
       if (ii.lt.i) go to 210
140    jj = jj + ii
       go to 180
210    jj = jj + 1
180    continue
c      build up the k-scalar (ak)
       ak = ak + temp*b(i,2)
       b(i,1) = temp
c      move ic to top of next a-matrix 'row'
190    ic = ic + i
c      form the q-ve+tor
       do 150  i=j1,n
150    b(i,1) = b(i,1) - ak*b(i,2)
c      transform the rest of the a-matrix
c      jj      start-1 of the rest of the a-matrix
       jj = id
c      move w-vector into the old a-matrix locations to save space
c      i       runs over the significant elements of the w-vector
       do 160  i=j1,n
       a(jj) = b(i,2)
       do 170  ii=j1,i
       jj = jj + 1
170    a(jj) = a(jj) - 2.0*(b(i,1)*b(ii,2) + b(i,2)*b(ii,1))
160    jj = jj + j
200    continue
c      move last codiagonal element out into its proper place
201    continue
       b(nm1,1) = a(nsize-1)
       a(nsize-1) = 0.
       b(n,1)=0.
c
c     sturm section.
c     sturm sequence iteration to obtain roots of tridiagonal form.
c     move diagonal elements into second n elements of b-vector.
c     this is a more convenient indexing position.
c     also, put square of codiagonal elements in third n elements.
      jump=1
      do 320 j=1,n
      b(j,2)=a(jump)
      b(j,3) = b(j,1)**2
320   jump = jump+j+1
      do 310 i=1,nroot
310   root(i) = +alimit
      rootl = -alimit
c     isolate the roots.  the nroot lowest roots are found, lowest first
      do 330 i=1,nroot
c     find current 'best' upper bound
      rootx = +alimit
      do 335 j=i,nroot
335   rootx = min(rootx,root(j))
      root(i) = rootx
c     get improved trial root
500   trial = (rootl + root(i))*0.5
c     if (trial.eq.rootl.or.trial.eq.root(i)) go to 330
      if(.not.abs(trial-rootl).gt.0..or..not.abs(trial-root(i)).gt.0.)
     1 goto 330
c     form sturm sequence ratios, using ortega's algorithm (modified).
c     nomtch is the number of roots less than the trial value.
350   continue
      nomtch = n
      j = 1
360   f0 = b(j,2) - trial
370   continue
      if (abs(f0).lt.theta1) go to 380
      if (f0.ge.0.) nomtch = nomtch - 1
      j = j + 1
      if (j.gt.n) go to 390
c     since matrix is normed to unity, magnitude of b(j,3) is less than
c     one, so overflow is not possible at the division step, since
c     f0 is greater than theta1.
c     f0 = b(j,2) - trial - b(j-1,3)/f0
      bb = b(j-1,3)/f0
c     if (f0.eq.0) bb=0.
      if(.not.abs(f0).gt.0.) bb=0.
      f0=b(j,2)-trial-bb
      go to 370
380   j = j + 2
      nomtch = nomtch - 1
      if (j.le.n) go to 360
390   continue
c     fix new bounds on roots
      if (nomtch.ge.i) go to 540
      rootl = trial
      go to 500
540   root(i) = trial
      nom = min0(nroot,nomtch)
      root(nom) = trial
      go to 500
330   continue
c     reverse the order of the eigenvalues, since custom dictates
c     'largest first'.  this section may be removed if desired without
c     affecting the remainder of the routine.
c     nrt = nroot/2
c     do 10 i=1,nrt
c     save = root(i)
c     nmip1 = nroot - i + 1
c     root(i) = root(nmip1)
c 10  root(nmip1) = save
c
c     trivec section.
c     eigenvectors of codiagonal form
807   continue
c     quit now if no vectors were requested.
      if (nrootx.lt.0) go to 1002
c     initialize vector array.
      do 15 i=1,n
      do 15 j=1,nroot
15    vect(i,j) = 1.0
      do 700 i=1,nroot
      aroot = root(i)
c     orthogonalize if roots are close.
      if (i.eq.1) go to 710
c     the absolute value in the next test is to assure that the trivec
c     section is independent of the order of the eigenvalues.
715   if (abs(root(i-1)-aroot).lt.toler) go to 720
710   ia = -1
720   ia = ia + 1
      elim1 = a(1) - aroot
      elim2 = b(1,1)
      jump = 1
      do 750  j=1,nm1
      jump = jump+j+1
c     get the correct pivot equation for this step.
      if (abs(elim1).le.abs(b(j,1))) go to 760
c     first (elim1) equation is the pivot this time.  case 1.
      b(j,2) = elim1
      b(j,3) = elim2
      b(j,4) = 0.
      temp = b(j,1)/elim1
      elim1 = a(jump) - aroot - temp*elim2
      elim2 = b(j+1,1)
      go to 755
c     second equation is the pivot this time.  case 2.
760   b(j,2) = b(j,1)
      b(j,3) = a(jump) - aroot
      b(j,4) = b(j+1,1)
      temp = 1.0
      if (abs(b(j,1)).gt.theta1) temp = elim1/b(j,1)
      elim1 = elim2 - temp*b(j,3)
      elim2 = -temp*b(j+1,1)
c     save factor for the second iteration.
755   b(j,5) = temp
750   continue
      b(n,2) = elim1
      b(n,3) = 0.
      b(n,4) = 0.
      b(nm1,4) = 0.
      iter = 1
      if (ia.ne.0) go to 801
c     back substitute to get this vector.
790   l = n + 1
      do 780 j=1,n
      l = l - 1
786   continue
      if (l-nm1) 787,788,789
789   elim1=vect(l,i)
      go to 785
788   elim1=vect(l,i)-vect(l+1,i)*b(l,3)
      go to 785
787   elim1=vect(l,i)-vect(l+1,i)*b(l,3)-vect(l+2,i)*b(l,4)
785   continue
c     if overflow is conceivable, scale the vector down.
c     this approach is used to avoid machine-dependent and system-
c     dependent calls to overflow routines.
      if (abs(elim1).gt.delbig) go to 782
      temp = b(l,2)
      if (abs(b(l,2)).lt.delta) temp = delta
      vect(l,i) = elim1/temp
      go to 780
c     vector is too big.  scale it down.
 782  do 784 k=1,n
784   vect(k,i) = vect(k,i)/delbig
      go to 786
780   continue
      go to (820,800), iter
c     second iteration.  (both iterations for repeated-root vectors).
820   iter = iter + 1
890   elim1 = vect(1,i)
      do 830 j=1,nm1
c     if (b(j,2).eq.b(j,1)) go to 840
      if(.not.abs(b(j,2)-b(j,1)).gt.0.) go to 840
c     case one.
850   vect(j,i) = elim1
      elim1 = vect(j+1,i) - elim1*b(j,5)
      go to 830
c     case two.
840   vect(j,i) = vect(j+1,i)
      elim1 = elim1 - vect(j+1,i)*temp
830   continue
      vect(n,i) = elim1
      go to 790
c     produce a random vector
801   continue
      do 802 j=1,n
c     generate pseudorandom numbers with uniform distribution in (-1,1).
c     this random number scheme is of the form...
c     rand1 = amod((2**12+3)*rand1,2**23)
c     it has a period of 2**21 numbers.
      rand1 = mod(4099.*rand1,rpower)
802   vect(j,i) = rand1/rpow1 - 1.
      go to 790
c
c     orthogonalize this repeated-root vector to others with this root.
800   if (ia.eq.0) go to 885
      do 860 j1=1,ia
      k = i - j1
      temp = 0.
      do 870 j=1,n
870   temp = temp + vect(j,i)*vect(j,k)
      do 880 j=1,n
880   vect(j,i) = vect(j,i) - temp*vect(j,k)
860   continue
885   go to (890,900), iter
c     normalize the vector
900   elim1 = 0.
      do 904 j=1,n
904   elim1 = max(abs(vect(j,i)),elim1)
      temp=0.
      do 910 j=1,n
      elim2=vect(j,i)/elim1
910   temp = temp + elim2**2
      temp=1.0/(sqrt(temp)*elim1)
      do 920 j=1,n
      vect(j,i) = vect(j,i)*temp
      if (abs(vect(j,i)).lt.del1) vect(j,i) = 0.
920   continue
700   continue
c
c      simvec section.
c      rotate codiagonal vectors into vectors of original array
c      loop over all the transformation vectors
       if (nm2.eq.0) go to 1002
       jump = nsize - (n+1)
       im = nm1
       do 950  i=1,nm2
       j1 = jump
c      move a transformation vector out into better indexing position.
       do 955  j=im,n
       b(j,2) = a(j1)
955    j1 = j1 + j
c      modify all requested vectors.
       do 960  k=1,nroot
       temp = 0.
c      form scalar product of transformation vector with eigenvector
       do 970  j=im,n
970    temp = temp + b(j,2)*vect(j,k)
       temp = temp + temp
       do 980  j=im,n
980    vect(j,k) = vect(j,k) - temp*b(j,2)
960    continue
       jump = jump - im
950    im = im - 1
1002   continue
c      restore roots to their proper size.
       do 95 i=1,nroot
95     root(i) = root(i)*anorm
1001   return
       end
