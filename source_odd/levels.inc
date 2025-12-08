c include block: levels                                                        
c
       REAL(8) :: j2matrix,j2squar                                        
        common  /levels/j2matrix(npsid,npsid),j2squar(npsid),
     &                  wate(npsid)  ,spn(npsid)     ,
     x                  spenrg(npsid),sphc(npsid)    ,
     x                  spevar(npsid),speold(npsid),
     x                  watep(npsid) ,aux(npsid)     ,oc(npsid)     ,          
     x                  vsqf(npsid)  ,dvsqf(npsid)   ,
     x                  vsqfol(npsid),cmpn(nsd,npsid),vsqobs(npsid),
     x                  deltaf(npsid),wtpair(npsid),
     x                  qspav(npsid) ,tcdiff(npsid),rtpartn(0:npsid),
     &                  angular(npsid),
     &                  akblock(1:nblint,2),
     &                  amblock(1:nblint,2),
     &                  aisoblo(1:nblint,2),
     &                  aistblo(1:nblint,2),testblock(nbliq),
     &                  edifblo(1:nblint,2),ediffblock(2),
     &                  testblock2(nbliq),updo(1:nblint,2),
     &			spinblo(nbliq),istcblock(nbliq),
     x                  mpsi(npsid)  ,npar(npsid)    ,nrpar(npsid) ,
     x                  nparz(2,npsid) ,
     x                  nqr(npsid)   ,nqz(npsid)     ,                         
     x                  icharg(npsid),iorder(npsid),diff(npsid),
     &             kblock(nbliq),mblock(nbliq),isoblo,istblo(nbliq),
     x                  ifocsv,maxblist(2),        ! flag for saving on 'oc'
     &			ktcblock(nbliq),
     &                  itpartn(0:npsid)






*                                                                     *        
***********************************************************************        
*                                                                     *        
*        single particle energies,occupation probabilities,           *        
*        and quantum numbers in the global ordering, i.e.             *        
*        each state has a label 'ist' which uniquely acceses          *        
*        the above quantum numbers. the fields are, in detail:        *        
*          mpsi   = k-0.5, where k is the z-component of total        *        
*                   angular momentum. note that mpsi is also the      *        
*                   z-component of the orbital a.m. for the upper     *        
*                   component with spin up.                           *        
*          npar   = z-parity for the upper component spin up          *        
*                   (component is=1). the subsequent components       *        
*                   have then z-parities                              *        
*                     -npar for upper,spin down (is=2)                *        
*                     -npar for lower,spin up   (is=3)                *        
*                     +npar for lower,spin down (is=4)                *        
*          nrpar  = r-parity for the upper component spin up          *        
*                   (component is=1). the subsequent components       *        
*                   have then r-parities                              *        
*                     -nrpar for upper,spin down (is=2)               *        
*                     +nrpar for lower,spin up   (is=3)               *        
*                     -nrpar for lower,spin down (is=4)               *        
*          icharg = isospin ( 1 for protons, 2 for neutrons)          *        
*          iorder = pointer field for reordering according to         *        
*                   the actual energy                                 *        
*          wate   = multiplicity due to spin degeneracy, in the       *        
*                   problem here always beeing 2.0 (kramers           *        
*                   degeneracy                                        *        
*          spn    = norm of the s.p. wave function, used as a         *        
*                   countercheck (it should be 1.0000)                *        
*          spenrg = single particle energy                            *        
*          spevar = variance of the single particle energy
*                                                                     *
*          vsqf   = occupation probability. in case of pure           *     
*                   slater determinants 1.0 throughout. in case       *        
*                   of frozen occupation the valence states are       *        
*                   occupied with 0.0<vsqf<1.0 . in case of           *        
*                   pairing the vsqf is computed dynamically          *        
*                   according to the actual spenrg in the sub-        *        
*                   routine vsqtr.                                    *        
*          dvsqf  = change in occupation v**2                         *
*          vsqfol = storage for old occupation probability from the   *
*              previous iteration.                                    *
*          vsqobs = occupation probablities for observables, from     *
*                   Lipkin-Nogami                                     *
*          deltaf = state dependent gap for dynamic BCS               *
*          wtpair = dynamic weight for BCS phase space                *
*          qspav  = diagonal element of constraint                    *
*          cmpn   = norm of the four components (is=1,..,nsd)         *        
*                   separately                                        *        
*          mpmax  = maximum of mpsi                                   *        
*       kblock   =  2*K values of blocked s.p. state                  *
*       mblock   =  M values of blocked s.p. state                    *
*       isoblo   =  isospin of blocked s.p. state                     *
*       istblo   =  number of blocked s.p. state in level sequence    *
*                                                                     *        
*	block.f							      *
*								      *
*       akblock  =  array for different kblock                        *
*       amblock  =  array for different mblock                        *
*	aisoblo  =  array for different isoblo                        *
*	aistblo  =  array for different istblo                        *
*	testblock  						      *
*								      *
***********************************************************************        
c
      common  /nucsys/  anumb,
     &                  nst(2),npart(2),iokk(2),numsta(2),lst
*                                                                     *        
***********************************************************************        
*                                                                     *        
*        arrays defining the nucleus (running over iso=1,2):          *        
*          nst    = total no. of states                               *        
*          npart  = particle number                                   *        
*          iokk   = pointer to the 1 st index of state with iso       *
*          numsta = number of available states for particle sort      *
*          lst    = total no. of states                               *        
*                                                                     *        
***********************************************************************        
c
      common  /bufcom/  len,nom,mhi(mxd),                                     
     x                  mlo(mxd),miq(mxd),mim(mxd),nbuf(npsid)
*                                                                             
***********************************************************************       
*                                                                     *        
*        variables defining the i/o blocking for buffering the        *        
*        wavefunctions to and from disk                               *        
*        mlo = lower state in block                                   *
*        mhi = upper state in block                                   *
*        miq = isospin in block                                       *
*        mim = m-value in block                                       *
*        nbuf = nr. of buffer for given state                         *
*                                                                     *        
***********************************************************************        
c end include block                                                            
