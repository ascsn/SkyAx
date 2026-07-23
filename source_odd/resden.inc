* 
* $RCSfile: resden.inc,v $
* 
* $Revision: 1.2 $
* $Date: 1995/04/12 13:22:54 $
* $Author: bender $
* $Locker: bender $
* $State: Exp $
*
**********************************************************
* 
* $Log: resden.inc,v $
* Revision 1.2  1995/04/12  13:22:54  bender
* Response-Beitraege in common-Block fuer Ausgabe in respsurf.plt
*
* Revision 1.1  1995/03/31  13:03:41  bender
* Initial revision
*
* 
**********************************************************
c include block: response densities                                             
c                                                                               
      common  /resden/  rrho(0:nzd,0:nrd,kqm),
     x                  rcrho(0:nzd,0:nrd),
     x                  rtau(0:nzd,0:nrd,kqm),
     x                  rrhols(0:nzd,0:nrd,kqm),
     x                  rd2rho(0:nzd,0:nrd,kqm),
     x                  rcoupt(0:nzd,0:nrd),
     x                  wcoltt(0:nzd,0:nrd),
     x                  wgtres(npsid)
*                                                                               
***********************************************************************         
*                                                                     *         
*         response-densities and -currents: (M.B.21.02.1995)          *
*         rrho   --> density                                          *
*         rtau   --> kinetic energy density                           *
*         rrhols --> spin orbit current                               *
*         rcrho  --> coulomb  density                                 *
*         rd2rho --> laplacian on 'rho'                               *
*                                                                     *         
***********************************************************************         
*
      common /linolr/  rspln(1:2),suv2rs(1:2),
     x                 xpotp,xkinp,xpd2p,xspop,xcol,xcex,
     x                 xpotn,xkinn,xpd2n,xspon,
     x                 ilnrsp,ilnini
* 
c end include block
