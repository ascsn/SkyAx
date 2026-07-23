c include block: derivs                                                     
c                                                                            
      common /derivs/ 
c
     &                derz1(0:nzd,0:nzd,4),
     &                derz2(0:nzd,0:nzd,4),
     &                derr1(0:nrd,0:nrd,4),
     &                derr2a(0:nrd,0:nrd,0:kmpsi,2),
c
     &                derz2f(0:nzd,0:nzd,4),
     &                derr1f(0:nrd,0:nrd,4),
     &                derr2f(0:nrd,0:nrd,4),
     &                derr2af(0:nrd,0:nrd,0:kmpsi,2),
c
     &                r1inv(0:nrd,0:nrd),
     &                r2inv(0:nrd,0:nrd),
c
     &                deri(0:nrd,0:nrd),xipo(0:nrd),
     &                zdnuc(0:nzd,0:nzd,4),                                
     &                rdnuc(0:nrd,0:nrd,0:kmpsi),
     &                zdmes(0:nzd,0:nzd,1),
     &                rdmes(0:nrd,0:nrd,1),
     &                rpro(0:nrd,0:nrd),
     &                zpro(0:nzd,0:nzd),
     &                fldbuf(0:nzd,0:nrd),
c
     &                iband_derr, 
     &                iband_derz
                       
c*********************************************************************** 
c   the derivative matrices:                                           *
c                                                                      *
c   derz1  --> first derivative in z-direction                         *       
c   derz2  --> second derivative in z-direction                        *       
c   derr1  --> first derivative in r-direction                         *       
c   derr2a --> bessel operator  in r-direction                         *
c                                                                      * 
c   derz1f --> first derivative in z-direction (fourier definition)    *       
c   derz2f --> second derivative in z-direction (fourier definition)   *      
c   derr1f --> first derivative in r-direction (fourier definition)    *       
c   derr2af--> bessel operator  in r-direction (fourier definition)    * 
c                                                                      *
c   r1inv  --> 1/r                                                     *
c   r2inv  --> 1/(r*r)                                                 *
c                                                                      *
c   zdnuc  --> damping operator for the nucleons in z-direction        *       
c   rdnuc  --> damping operator for the nucleons in r-direction        *       
c   zdmes  --> damping operator for the photon in z-direction          *
c   rdmes  --> damping operator for the photon in r-direction          *
c   npa = 2 are the two different parities.                            *       
c***********************************************************************       
c end include block                                                            
