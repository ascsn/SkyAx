c include block: wavefr                                                       
c                                                                             
      common  /wavefr/  psi(0:nzd,0:nrd,nsd,npsid),
     &                  dpsiz(0:nzd,0:nrd,nsd,npsid),
     x                  psip(0:nzd,0:nrd,nsd),psin(0:nzd,0:nrd,nsd),
     x                  psi1(0:nzd,0:nrd,nsd),psi2(0:nzd,0:nrd,nsd),
     x                  psi3(0:nzd,0:nrd,nsd),psi4(0:nzd,0:nrd,nsd),
     x                  wrk1(0:nzd,0:nrd),wrk2(0:nzd,0:nrd),
     x                  wrk3(0:nzd,0:nrd),wrk4(0:nzd,0:nrd)
*
***********************************************************************       
*                                                                     *        
*        the wavefunctions:                                           *
*         psi      = the main storage for all s.p. states             *
*         dpsiz    = storage for gradient step                        *
*                    also used for z-derivative in crankin.f          *
*         psip     = auxiliary storage for one s.p. state
*         psin     = auxiliary storage for one s.p. state
*         psi1..4  = auxiliary fields (with spin)
*         wrk1..4  = auxiliary fields (without spin)
*                                                                     *        
***********************************************************************        
c end include block                                                            
