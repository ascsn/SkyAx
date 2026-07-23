      parameter (nzdc=kcrank*nzd,   nrdc=kcrank*nrd)
      parameter (ktrian=(npsid*npsid+npsid)/2) ! space for triangular storage
      parameter (ktriac=kcrank*ktrian+1)
c
c     fields for calculation the SC cranking mass
c      dpsir   = storage for d/dr derivative of wave functions
cc      dpsiz   = storage for d/dz derivative of wave functions
c      curqr   = respsonse current in r-direction
c      curqz   = respsonse current in z-direction
c      curqph  = respsonse current in phi-direction (only m.neq.0)
c      spinr   = respsonse spin-densiity in r-direction (only m.neq.0)
c      spinz   = respsonse spin-densiity in z-direction (only m.neq.0)
c      spinph  = respsonse spin-densiity in phi-direction
c      dynpoz  = response potential 'A' in z-direction
c      dynpor  = response potential 'A' in r-direction
c      dynpop  = response potential 'A' in phi-direction
c      spipor  = spin response potential 'S' in r-direction
c      spipoz  = spin response potential 'S' in z-direction
c      spipop  = spin response potential 'S' in phi-direction
c      rhop    = P-response density
c      taup    = P-response kinetic density
c      drhlsp  = P-response spin-orbit density
c      paircp  = P-response pair current
c      e2quas  = two quasi-particle energies
c      iposmx  = maximum number of matrix elements 
c
      common  /crankw/  
     &     dpsir(0:nzdc,0:nrdc,nsd,npsid),
c     &     dpsiz(0:nzdc,0:nrdc,nsd,npsid),
     &     curqr(0:nzdc,0:nrdc,kqm),
     &     curqz(0:nzdc,0:nrdc,kqm),
     &     curqph(0:nzdc,0:nrdc,kqm),
     &     spinr(0:nzdc,0:nrdc,kqm),
     &     spinz(0:nzdc,0:nrdc,kqm),
     &     spinph(0:nzdc,0:nrdc,kqm),
     &     dynpor(0:nzdc,0:nrdc,kqm),
     &     dynpoz(0:nzdc,0:nrdc,kqm),
     &     dynpop(0:nzdc,0:nrdc,kqm),
     &     spipor(0:nzdc,0:nrdc,kqm),
     &     spipoz(0:nzdc,0:nrdc,kqm),
     &     spipop(0:nzdc,0:nrdc,kqm),
     &     rhop(0:nzdc,0:nrdc,kqm),
     &     taup(0:nzdc,0:nrdc,kqm),
     &     drhlsp(0:nzdc,0:nrdc,kqm),
     &     paircp(0:nzdc,0:nrdc,kqm),
     &	   qresp(ktriac),presp(ktriac),ptest(ktriac),	
     &     qab01(ktriac),qab10(ktriac),qab11(ktriac),
     &     e2quas(ktriac),
     &     iposmx

      


