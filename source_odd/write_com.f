c--------write_com---------------------------------------------------
c

      subroutine write_com(outnam)



      include 'paramr.inc'
*
!      character*80 title
      character*20 set
      character*3  outnam
c
      include 'constr.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'inipar.inc'
      include 'observ.inc'
      include 'regime.inc'






      open(101,status='unknown',form='formatted',
     & file='reconstr.'//outnam)
      write(101,*)    
     &             qcfac(0:nzd,0:nrd), distsf(0:nzd,0:nrd),
     &             qfac(0:nzd,0:nrd),
     &             qcst,qsoll,qtil,qtilc,qtild,qsq,qsqd,qlamb,qcorr,
     &             c0,d0,qeps,dalp,dgam,swicon,swifac,
     &             qmin,qmax,qw,
     &             itcons,iflin,idampc,lqcon,iflmb0,
     &             lrho(0:nzd,0:nrd)
      
c      open(102,status='unknown',form='formatted',
c     & file='recoulomb.'//outnam)
c      write(102,*)
c     &                  cbgrho(0:nzd,0:nrd,0:KANGMOL,0:KANGMOM),
c     &                  cbgfld(0:nzd,0:nrd,0:KANGMOL,0:KANGMOM),
c     &                  cbgmom(0:KANGMOL),
c     &                  cfldbuf(0:nzd,0:nrd),
c     &                  cdmpmas(0:nzd,0:nrd)
      open(103,status='unknown',form='formatted',
     & file='rederivs.'//outnam)
      write(103,*)
     &                derz1(0:nzd,0:nzd,1:4),
     &                derz2(0:nzd,0:nzd,1:4),
     &                derr1(0:nrd,0:nrd,1:4),
     &                derr2a(0:nrd,0:nrd,0:kmpsi,1:2),
     &                derz2f(0:nzd,0:nzd,1:4),
     &                derr1f(0:nrd,0:nrd,1:4),
     &                derr2f(0:nrd,0:nrd,1:4),
     &                derr2af(0:nrd,0:nrd,0:kmpsi,1:2),
     &                r1inv(0:nrd,0:nrd),
     &                r2inv(0:nrd,0:nrd),
     &                deri(0:nrd,0:nrd),xipo(0:nrd),
     &                zdnuc(0:nzd,0:nzd,1:4),                                
     &                rdnuc(0:nrd,0:nrd,0:kmpsi),
     &                zdmes(0:nzd,0:nzd,1),
     &                rdmes(0:nrd,0:nrd,1),
     &                rpro(0:nrd,0:nrd),
     &                zpro(0:nzd,0:nzd),
     &                fldbuf(0:nzd,0:nrd),
     &                iband_derr, 
     &                iband_derz
      
      open(104,status='unknown',form='formatted',
     & file='reinterpr.'//outnam)
      write(104,*)
     &                  cmono,cquad,chex,cdip,coct,
     &                  cmssqr,colerr,snl,snlhc,yukerr,dmsfac,
     &                  errenr,errvar,
     &                  irepm,irepc,itpfix,
     &                  itert,swdiag     
      
      open(105,status='unknown',form='formatted',
     & file='regridpr.'//outnam)
      write(105,*)    
     &                  dr,dz,sr,sz,drsqr,dzsqr,dvol,sym,
     &                  nr,nz,nrm1,nzm1,nro,nzo,isym,icm,
     &                  icutr,icutz,idrtyp
      open(106,status='unknown',form='formatted',
     & file='redampar.'//outnam)
      write(106,*)
     &           x0inv,e0inv,addrho,addpot,addocc,x0fac,adyrho,
     &           adypot,mpmax,itsub

      open(107,status='unknown',form='formatted',
     & file='reforces.'//outnam)
      write(107,*)
     &        x0inv,e0inv,addrho,addpot,addocc,x0fac,adyrho,
     &        adypot,mpmax,itsub
      open(108,status='unknown',form='formatted',
     & file='repairpr.'//outnam)
      write(108,*) 
     &   ferp(1:2),ggg(1:2),gap(1:2),fer(1:2),fermx(1:2),ewin(1:2),
     &   epair(1:2),fer2(1:2),pairfc(1:2),feradd(1:2),pairstab,
     &   rho0pr,addpar,epspar,wfac,errf,erra,
     &   stwid,stslop,rhocut,ifstab,
     &   npiter,ipair,itpini
      
      open(109,status='unknown',form='formatted',
     & file='refields.'//outnam)
      write(109,*)
     &         pot(0:nzd,0:nrd,1:kqm),
     &         coulpt(0:nzd,1:nrd),
     &         bmass(0:nzd,0:nrd,1:kqm),
     &         potls(0:nzd,0:nrd,1:kqm,1:2),
     &         drbmas(0:nzd,0:nrd,1:kqm),dzbmas(0:nzd,0:nrd,1:kqm),
     &         potadd(0:nzd,0:nrd,1:kqm),rhosur(0:nzd,0:nrd,1:kqm),
     &         pairpt(0:nzd,0:nrd,1:kqm),potdamp(0:nzd,0:nrd)
      
      open(110,status='unknown',form='formatted',
     & file='recurren.'//outnam)
      write(110,*)   
     &                  rho(0:nzd,0:nrd,1:kqm),
     &                  crho(0:nzd,0:nrd),
     &                  tau(0:nzd,0:nrd,1:kqm),
     &                  drhols(0:nzd,0:nrd,1:kqm),
     &                  rholsr(0:nzd,0:nrd,1:kqm),    
     &                  rholsz(0:nzd,0:nrd,1:kqm),    
     &                  d2rho(0:nzd,0:nrd,1:kqm),
     &                  rhoti(0:nzd,0:nrd,1:kqm),
     &                  paircr(0:nzd,0:nrd,1:kqm)
      
      open(111,status='unknown',form='formatted',
     & file='recoulom.'//outnam)
      write(111,*)    
     &                  crhom(0:nzd,0:nrd,0:KANGMO),
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
      
      open(112,status='unknown',form='formatted',
     & file='reinipar.'//outnam)
      write(112,*)
     &                  hbomz,hbomr,betain,octini,bzosc,brosc,
     &                  eocc(1:2),deocc(1:2),mmax(1:2)
c      ,locc(1:10)

      open(113,status='unknown',form='formatted',
     & file='relevels.'//outnam)
      write(113,*)
     &         wate(1:npsid)  ,spn(1:npsid)     ,
     &         spenrg(1:npsid),sphc(1:npsid)    ,
     &         spevar(1:npsid),speold(1:npsid),
     &        watep(1:npsid) ,aux(1:npsid)     ,oc(1:npsid),          
     &        vsqf(1:npsid)  ,dvsqf(1:npsid) ,
     &        vsqfol(1:npsid),cmpn(nsd,1:npsid),vsqobs(1:npsid),
     &        deltaf(1:npsid),wtpair(1:npsid),
     &        qspav(1:npsid) ,
     &        mpsi(1:npsid)  ,npar(1:npsid)    ,nrpar(1:npsid) ,
     &        nparz(1:2,1:npsid) ,
     &        nqr(1:npsid)   ,nqz(1:npsid)     ,                         
     &        icharg(1:npsid),iorder(1:npsid),diff(1:npsid),
c     &        kblock(1:nbliq),mblock(1:nbliq),isoblo,istblo(1:nbliq),
     &        ifocsv,maxblist(1:2),        ! flag for saving on 'oc'
     &        akblock(1:nblint,1:2),
     &        amblock(1:nblint,1:2),
     &        aisoblo(1:nblint,1:2),
     &        aistblo(1:nblint,1:2)


      
      open(114,status='unknown',form='formatted',
     & file='renucsys.'//outnam)
      write(114,*)
     &                   anumb,
     &                  nst(1:2),npart(1:2),iokk(1:2),numsta(1:2),lst
      
      open(115,status='unknown',form='formatted',
     & file='rebufcom.'//outnam)
      write(115,*)
     &     len,nom,mhi(1:mxd),                                     
     &     mlo(1:mxd),miq(1:mxd),mim(1:mxd)
      
      open(116,status='unknown',form='formatted',
     & file='reobserv.'//outnam)
      write(116,*)   
     &                atnum,apnum,annum,rms,rmsp,rmsn,q,qp,qn,h,hp,hn,       
     &                beta1,beta1p,beta1n,beta2,beta2p,beta2n,
     &                beta3,beta3p,beta3n,beta4,beta4p,beta4n,
     &                zcm,rchrg2,rdmsc,surfc,
     &                crmaso,crmasp,crmasq,xinert,bmassc,bmassr,scmass,
     &                pnmass,pnmas0,zpetri,omgtri,
     &                pwid,pwido,ammass,amwid,ampnms,
     &                zpecmx,widcmx,zpecmz,widcmz,widdxy,widdz,
     &                tmass  
      
      open(117,status='unknown',form='formatted',
     & file='recntrlv.'//outnam)
      write(117,*)    
     &                  econ,eold,energy,energr,ekin,epot,ecoul,etpair,
     &                  eshell,erearr,equad,devh2,err,qerr,qold,qcon,
     &                  devmf2,devoc2,devgap,varqsq,qq,qdh,
     &                  itmax,ivs,itest,iprint,iplot,ifcran
      
      open(118,status='unknown',form='formatted',
     & file='reord.'//outnam)
      write(118,*)
     &             iord2(1:npsid),spen2(1:npsid),
     &             ireg(1:2)

      open(119,status='unknown',form='formatted',
     & file='rewavefr.'//outnam)
      write(119,*)
     &                  psi(0:nzd,0:nrd,1:nsd,1:npsid),
     &                  dpsiz(0:nzd,0:nrd,1:nsd,1:npsid),
     &                  psip(0:nzd,0:nrd,1:nsd),psin(0:nzd,0:nrd,1:nsd),
     &                  psi1(0:nzd,0:nrd,1:nsd),psi2(0:nzd,0:nrd,1:nsd),
     &                  psi3(0:nzd,0:nrd,1:nsd),psi4(0:nzd,0:nrd,1:nsd),
     &                  wrk1(0:nzd,0:nrd),wrk2(0:nzd,0:nrd),
     &                  wrk3(0:nzd,0:nrd),wrk4(0:nzd,0:nrd)

      close(101)
c      close(102)
      close(103)
      close(104)
      close(105)
      close(106)
      close(107)
      close(108)
      close(109)
      close(110)
      close(111)
      close(112)
      close(113)
      close(114)
      close(115)
      close(116)
      close(117)
      close(118)
      close(119)

      return

      end 
      
     
