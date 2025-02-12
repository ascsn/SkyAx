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
*                  kqm    -   no. of isospins  (2)
*
*                  npsid  -   maximum no. of wavefunctions            *        
*                  nhwd   -   maximum no. of oscillator quanta        *        
*                  kmpsi  -   maximal number of mpsi for damping      *        
*                             in r-direction                          *        
*                  npa    -   parity +/- one                          *        
*                  kcrank -   switch for cranking part of code        *
*                  kintpl -   switch to write s.p. properties for     *
*                             later interpolation of cranking mass    *
*                  klocal -   switch to print of localization         *
*                  kwrit3 -   switch to writing wf's on 3D grid       *
*                  kangmo -   number of a.m. for Coulomb solver       *
*                             =0 suffices for ground state            *
*                             =2 is recommended for linear response   *
*                  ktestm -   switches to test rectangular storage    *
*                             using m=0 mode                          *
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
c
      logical KTESTM,KRESTA
c        
c      parameter  (nzd= 99, nrd= 49, nmx= 99, nsd=2, npsid=590,
c       parameter  (nzd= 75, nrd= 50, nmx= 75, nsd=2, npsid=490,
c       parameter  (nzd= 75, nrd= 50, nmx= 75, nsd=2, npsid=590,
c      parameter  (nzd= 150, nrd= 70, nmx= 219, nsd=2, npsid=590,
       parameter  (nzd= 61,  nrd= 61, nmx= 67, nsd=2, npsid=890,
c       parameter  (nzd= 26,  nrd= 26, nmx= 67, nsd=2, npsid=890,
c       parameter  (nzd= 31,  nrd= 23, nmx= 31, nsd=2, npsid=990,
     x            nhwd=120, npa=  2, kqm=  2, kmpsi=19,  
     x            kcrank=1,kintpl=0,
     x            klocal=0,KANGMO=0,KTESTM=.FALSE.,KRESTA=.FALSE.,
     x            NX3=1,NY3=1,NZ3=1,               ! >1 activates 3D prints
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
     x            lwavef=npsid*nsd*nrzd, n2zd=2*nzd,                           
     x            nbliq=2, nblint=90)                        ! integers for blocking
c
c
c
      parameter (zero=0.0D0)
      parameter (one=1.0D0)
      parameter (two=2.0D0)
      parameter (pi=3.141592653589793238462643383279502884197169399D0)
      parameter (edr=0.333333333333333)
      parameter (oh=0.5D0)
      parameter (eighty=80.0D0)
      parameter (smal=1.0D-30)
      parameter (ifoccs=0)                       ! switch to print occ.numbers
      LOGICAL,PARAMETER :: ttrack=.TRUE.
!      REAL(8),PARAMETER :: g_magnetic(2)=(/5.584694712D0,-3.82608544D0/)
      REAL(8),PARAMETER :: g_magnetic(2)=(/0D0,1D0/)
      REAL(8),PARAMETER :: v2lowlim=1D-20        ! lower limit for v**2
c
      common  /constv/  hc,e2,r0,
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


	common /timesym /  isTSym


c end include block                                                           







