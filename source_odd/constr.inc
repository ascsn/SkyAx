c include block: constr
c
*kr******************
      logical lrho(0:nzd,0:nrd)
      REAL(8),PARAMETER :: c0_suv=0.2D0
      common /constr/
     &                qcfac(0:nzd,0:nrd), distsf(0:nzd,0:nrd),
     &                qfac(0:nzd,0:nrd),
     &                qfacquad(0:nzd,0:nrd),   ! for additional quadr. constr.
     &                qcst,qsoll,qtil,qtilc,qtild,qsq,qsqd,qlamb,qcorr,
     &                c0,d0,qeps,dalp,dgam,swicon,swifac,
     &                qmin,qmax,qw,
     &                tempmn,tempmx,tempw,
     &                qlambquad,               ! for additional quadr. constr.
     &                suv_wanted(2),suv_mue(2),suv_lambda(2), ! for sum(uv) constraint
     &                itcons,iflin,idampc,lqcon,iflmb0,
     &                lrho
*kr******************
*
************************************************************************
*
*    contains parameters used for the quadrupol-constraint calculation
c      qcfac   = raw constraining force, e.g. quadrupole
c      qfac    = damped constraing force, actually used
c      distsf  = field of distances to the chosen surface line
*
************************************************************************
*
c end include block
