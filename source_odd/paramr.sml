c include block: paramr
c
c  to make the program run double precision, just enable the 'implicit'
c  statement
c
       implicit real*8 (a-h,o-z)
c                                                                               
***********************************************************************
*        parameter statments for setting dimensions:                  *         
*                                                                     *         
*                  nzd    -   no. of mesh points in z direction.      *         
*                  nrd    -   no. of mesh points in r direction.      *         
*                  nsd    -   no. of spinor components.               *         
*                  kqm    -   no. of isospins  (2)                    *
*                  npsid  -   maximum no. of wavefunctions            *         
*                  nhwd   -   maximum no. of oscillator quanta        *         
*                  kmpsi  -   maximal number of mpsi for damping      *         
*                             in r-direction                          *         
*                  npa    -   parity +/- one                          *         
*                                                                     *         
*        derived parameters:                                          *         
*                                                                     *         
*                  nmx    -   max(nrd,nzd)                            *         
*                  nzd1   -   nzd-1                                   *         
*                  nrd1   -   nrd-1                                   *         
*                  nrzd   -   (nrd+1)*(nzd+1)                         *         
*                  mxd    -   maximum no. of allowed m values.        *         
*                             =nhwd+1                                 *         
*                  mxd1   -   mxd+1                                   *         
*                  nbufd  -   maximum no. of buffer blocks            *         
*                             nbufd=npsid                             *         
*                                                                     *         
*                  nnrxd  -   (nhwd+1)/2 -                            *         
*                             maximum principle quantum number for    *         
*                             the cylindrical oscillator functions    *         
*                             (r-direction)                           *         
*                  nnzxd  -   nhwd -                                  *         
*                             maximum principle quantum number for    *         
*                             the cylindrical oscillator functions    *         
*                             (z-direction)                           *         
*                  nnzxd1     nnzxd+1                                 *         
*                  nnrxd1     nnrxd+1                                 *         
*                                                                     *         
*                  lcurre -   length of common block curren in words  *         
*                  ldiffc -   length of common block diffcf in words  *         
*                  lwavef -   length of common block wavefn in words  *         
*                  lyukww -   leng-3nrdp1 of c.block yukwww in words  *         
*                  linvco -   length of common block invcof in words  *         
*                  lgrdwg -   length of common block gridwg in words  *         
*                  ldmpar -   length-3 of common block dampar         *         
*                                                                     *         
***********************************************************************         
      parameter  (nzd= 27, nrd= 19, nmx= 27, nsd=2, npsid=30,
     x            nhwd=30, npa=  2, kqm=  2, kmpsi=11,
     x            nzdp1=nzd+1, nrdp1=nrd+1, mxd=nhwd+1, mxd2=2*mxd,             
     x            mxd1=mxd+1, nzd1= nzd-1, nrd1=nrd-1, nbufd=npsid,             
     x            nrzd=(nrd+1)*(nzd+1), n4rzd=4*nrzd, nnzxd=nhwd,               
     x            nnrxd=(nhwd+1)/2, nnrxd1=nnrxd+1, nnzxd1=nnzxd+1,             
     x            lcurre=n4rzd, linvco=10*nrzd, lmeson=n4rzd,                   
     x        ldiffc=2*nsd*(nrzd+1+nrdp1), lyukww=nrzd+3*(nzdp1+nrdp1),         
ctestx            ldiffc=6*nsd*nrzd, lyukww=nrzd+3*(nzdp1+nrdp1),               
ctestx            ldiffc=10*nsd*nrzd, lyukww=nrzd+3*(nzdp1+nrdp1),              
     x            ldmpar=0,                                                     
     x            lgrdwg=15*(nrd+1)+2*(nzd+1)*nsd,                              
     x            lwavef=npsid*nsd*nrzd, n2zd=2*nzd             )               
c
c
c
      common  /constv/  zero,pi,hc,e2,r0,
     &                  edr,oh,one,two,eighty,smal,
     &                  iou
*                                                                               
***********************************************************************         
*                                                                     *         
*        physical and mathematical constants                          *         
*          pi     = 3.1415926                                         *         
*          h2m    = hbar**2/baryonmass                                *         
*          hc     = hbar*c                                            *         
*          e2     = electron charge squared                           *         
*          r0     = pion range (1.40fm)                               *         
*          some double precision numbers:                             *
*          edr    = 1.0/3.0                                           *
*          dzero  = 0.0d0                                             *         
*          oh     = 0.5d0                                             *         
*          one    = 1.0d0                                             *         
*          two    = 2.0d0                                             *         
*          eighty = 80.0d0                                            *         
*          smal   = 1.0e-30                                           *         
*          iou    = unit number for plot of fields                    *
***********************************************************************         
c end include block                                                             
