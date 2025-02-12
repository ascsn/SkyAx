c include block: observ                                                       
c
      logical tmass
      common /observ/ atnum,apnum,annum,rms,rmsp,rmsn,q,qp,qn,h,hp,hn,       
     &                beta1,beta1p,beta1n,beta2,beta2p,beta2n,
     &                beta3,beta3p,beta3n,beta4,beta4p,beta4n,
     &                zcm,rchrg2,rdmsc,surfc,rmsc,r4mom,
     &                crmaso,crmasp,crmasq,xinert,bmassc,bmassr,scmass,
     &                pnmass,pnmas0,zpetri,omgtri,
     &                pwid,pwido,ammass,amwid,ampnms,
     &                zpecmx,widcmx,zpecmz,widcmz,widdxy,widdz,
     &                energyFy,
     &                symcondi(2),symcondf(2),
     &                tmass  
c                                                                             
c     the observables                                                         
c
       common  /cntrlv/  econ,eold,energy,energr,ekin,epot,ecoul,etpair,
     &                  ecolex,
     &                  eshell,erearr,equad,devh2,err,qerr,qold,qcon,
     &                  devmf2,devoc2,devgap,varqsq,qq,qdh,symmiso(2),
     &                  itmax,ivs,itest,iprint,iplot,ifcran
*                                                                              
***********************************************************************        
*                                                                     *        
*        control variables for the iteration,printing,and plotting    *        
*        for the energies,densities,etc.                              *
*      energr = total energy from shell- and rearrangement-energy     *
*      energy = total energy directly from densities                  *
*      etpair = total pairing energy (or -T*S in case of thermal)     *
*      egrand = grand canonical potential                             *
*                                                                     * 
***********************************************************************         
c end include block                                                             






