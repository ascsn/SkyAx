c include block: fields
c
      common  /fields/  pot(0:nzd,0:nrd,kqm),
     x                  coulpt(0:nzd,0:nrd),
     x                  bmass(0:nzd,0:nrd,kqm),
     x                  potls(0:nzd,0:nrd,kqm,0:2),
     x                  drbmas(0:nzd,0:nrd,kqm),dzbmas(0:nzd,0:nrd,kqm),
     x                  potadd(0:nzd,0:nrd,kqm),rhosur(0:nzd,0:nrd,kqm),
     x                  pairpt(0:nzd,0:nrd,kqm),potdamp(0:nzd,0:nrd),
     x                  potS(0:nzd,0:nrd,kqm,2),
     x                  potC(0:nzd,0:nrd,kqm,2),
     x                  potA(0:nzd,0:nrd,kqm)
*                                                                     *
***********************************************************************
*                                                                     *
*        the potential fields:                                        *
*           pot    -->  local hf potential                            *
*           bmass  -->  inverse mass                                  *
*           potls  -->  spin orbit potential                          *
*           coulpt -->  coulomb potential                             *
*           drbmas -->  d/dr bmass                                    *
*           dzbmas -->  d/dz bmass                                    *
*           potadd -->  additional potential to stabilize positive    *
*                       energy states                                 *
*           pairpt -->  gap potential for dynamic BCS                 *
*                                                                     *
***********************************************************************
c
      common  /curren/  rho(0:nzd,0:nrd,kqm),
     x                  crho(0:nzd,0:nrd),
     x                  tau(0:nzd,0:nrd,kqm),
     x                  drhols(0:nzd,0:nrd,kqm),
     x                  rholsr(0:nzd,0:nrd,kqm),        ! r-comp. l*s current
     x                  rholsz(0:nzd,0:nrd,kqm),        ! z-comp. l*s current
     x                  d2rho(0:nzd,0:nrd,kqm),
     x                  rhoti(0:nzd,0:nrd,kqm),
     x                  paircr(0:nzd,0:nrd,kqm),
     x                  rotjr(0:nzd,0:nrd,kqm),         ! r-comp rot j
     x                  rotjz(0:nzd,0:nrd,kqm),         ! z-comp rot j
     x                  rotspden(0:nzd,0:nrd,kqm),	! rot spin-density
     x                  spdenr(0:nzd,0:nrd,kqm),	! r-comp. spin-desity
     x                  spdenz(0:nzd,0:nrd,kqm),	! z-comp. spin-density
     x                  d2spdenr(0:nzd,0:nrd,kqm),	! r-comp. Laplace spin-desity
     x                  d2spdenz(0:nzd,0:nrd,kqm),	! z-comp. Laplace spin-density
     x                  curphi(0:nzd,0:nrd,kqm),	! phi-comp. current
     x                  spdenkr(0:nzd,0:nrd,kqm),       !r-comp.kin spin-den
     x                  spdenkz(0:nzd,0:nrd,kqm)        !z-comp.kin spin-den


*
***********************************************************************
*                                                                     *
*        densities and currents:                                      *
*          rho --> density                                            *
*          tau --> kinetic energy density                             *
*          rhols --> spin orbit current                               *
*          crho --> coulomb  density                                  *
*          d2rho --> laplacian on 'rho'                               *
*          rhoti --> the density for the constraint calculation       *
*          paircr --> pairing density                                 *
*                                                                     *
***********************************************************************
c
      common  /coulom/  crhom(0:nzd,0:nrd,0:KANGMO),
     &                  coulm(0:nzd,0:nrd,0:KANGMO),
     &                  crhoq(0:nzd,0:nrd,0:KANGMO),
     &                  coulq(0:nzd,0:nrd,0:KANGMO),
     &                  crhoh(0:nzd,0:nrd,0:KANGMO),  
     &                  coulh(0:nzd,0:nrd,0:KANGMO),
     &                  crhod(0:nzd,0:nrd,0:KANGMO),
     &                  could(0:nzd,0:nrd,0:KANGMO),
     &                  crhoo(0:nzd,0:nrd,0:KANGMO),
     &                  coulo(0:nzd,0:nrd,0:KANGMO),
     &                  wcolt(0:nzd,0:nrd),sqrm(0:nzd,0:nrd),
     &                  nqrun,nhrun,ndrun,norun
*                                                                     *
***********************************************************************
c
c     fields for the coulomb solver:
c     monopole and quadrupole densities and fields
c
***********************************************************************
*
c end include block
