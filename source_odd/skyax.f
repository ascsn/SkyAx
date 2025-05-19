c-----skyax -----------------------------------------------------------
c
c
      program skyax
c
***********************************************************************
*                                                                     *
*                                                                     *
*   skyax  -   pgm to solve the static mean field equations           *
*              in cylindrical coordinates using Skyrme forces.        *
*              The numerical representation uses a full matrix        *
*              for the first and second derivatives. These            *
*              matrices are built from a Fourier definition.          *
*                                                                     *
*              pgm uses ansi fortran 77 standard coding.              *
*                                                                     *
*              this program has been developed from a                 *
*              relativistic mean-field code which, in turn,           *
*              had is roots in a Skyrme mean-field code using         *
*              finite differences.                                    *
*              maintainance of the program should carry a             *
*              version number, the history of the changes,            *
*              and the corresponding authors.                         *
*                                                                     *
*              History:                                               *
*               ~1982        --> axial Skyrm-Hf with 3-point F.D. .   *
*                                S.Umar, M.R. Strayer, S.J.Lee,       *
*                                P.-G.Reinhard                        *
*               ~1985        --> change to Relativistic MF with       *
*                                3-point F.D. and shifted grids.      *
*                                S.Umar, M.R. Strayer, S.J.Lee        *
*                                J. Maruhn, J. Fink, P.-G.Reinhard    *
*               ~1991        --> change of the Relativistic code      *
*                                to Fourier representations of        *
*                                the derivative matrices.             *
*                                J. Maruhn, J. Fink, P.-G.Reinhard    *
*               March 1993   --> change to Skyrme HF code using       *
*                                Fourier representations of           *
*                                the derivative matrices.             *
*                                P.-G.Reinhard                        *
*               1993-1993    --> testing and extension                *
*                                (e.g. LN pairing)                    *
*                                P.-G.Reinhard, M.Bender, J.Maruhn    *
*               Oct/Nov1996  --> z-asymmetry, delta-force BCS         *
*                                D.Dean, P.-G.Reinhard, M.R.Strayer   *
*                                                                     *
***********************************************************************
*                                                                     *
*      Version     Date     Status                                    *
*       1.0.0    03/13/93   first setup of program.                   *
*                           Coulomb, spin-orbit, and pairing run.     *
*                           Agreement with the spherical code is      *
*                           tested to be better than 0.03% for        *
*                           doubly magic nuclei (i.e. spherical       *
*                           and without pairing).                     *
*       1.0.1    03/17/93   Input of initial maximum oscillator shell *
*                           is no done with 'deocc', relative to a    *
*                           first guess of the appropriate valence    *
*                           shell for given 'npart'. It has also been *
*                           added a check whether sufficient states   *
*                           are supplied in the 'ssetup'.             *
*       1.0.2    03/22/93   More forces in PARAM and flexible input   *
*                           of arbitrary parameters;                  *
*                           clean-up of plot files;                   *
*                           new form of printing s.p.properties;      *
*                           new pairing routine;                      *
*                           two different convergence criteria;       *
*                           fixing some minor bugs.                   *
*       1.0.3    05/28/93   damping of propagation of occupation      *
*                           numbers with 'addocc',                    *
*                           dynamical mixture of densities with       *
*                           'adddyn' replacing 'addnew'.              *
*                           minor changes in output formats.          *
*       1.0.4    09/15/93   direct evaluation of energy,              *
*                           distinction between 'addrho' and 'addpot'.*
*       1.0.5    09/16/94   new damping of constraint,                *
*                           Lipkin-Nogami (preliminary),              *
*                           various test options                      *
*       1.0.6    10/17/94   stabilizing wall (STAPOT),                *
*                           alternatve distance (DELRMN) for damping  *
*                           of constraint,                            *
*                           debugging such that variables can be      *
*                           supposed as 'dynamic' in subroutines      *
*       2.0.1.   11/16/96   Extension to asymmetry along z-axis,      *
*                           octupole constraint,                      *
*                           initial octupole "kick",                  *
*                           BCS-pairing with delta- or surface-force, *
*                           revised structure of input file.          *
*                                                                     *
***********************************************************************
*                                                                     *
*      Here follows a list of topics for the near future development: *
*                     --------------         -----------              *
*         -   Deformation properties and pairing should still be      *
*             tested, although both seems superficially o.k.          *
*             In particular, watch the convergence of the various     *
*             deformations 'qtil', 'qold', and 'q' towards each       *
*             other; there may still be a problem.                    *
*             This problem exists, in particular, if the previous     *
*             solution is reused in a deformation or temperature      *
*             loop.                                                   *
*         -   The iteration of the quadrupole constraint may also     *
*             be speeded up if that turns out to be necessary.        *
*         -   Change input of desired (constrained) deformation to    *
*             the dimensionless quadrupole momentum 'beta'.           *
*         -   Improve the treatment of the boundary 'r=0' for the     *
*             spin-orbit current.                                     *
*         -   The handling of the c.m. zero-point energy may be       *
*             improved by computing <P^2> microscopically. A similar  *
*             zero-point energy for rotation may be added.            *
*                                                                     *
*         !   yet to code: s.c. resp[onse for Lipkin Nogami           *
*                          use of 'vsqobs' for spatial observables    *
*                                                                     *
*                                                                     *
***********************************************************************
*                                                                     *
*      sample input file:                                             *
*
* 'net: SkI1, 20Ne PES'             title (first 3 digits determine output name)
* 10,10,                            nucleus: nprot,nneut
* -1913.619, 439.809, 2697.594 , 10592.267, ,     t0,t1,t2,t3
* -0.954536, -5.782388, -1.287379, -1.561421,     x0,x1,x2,x3
* 62.130,   62.130, 0.25  , 20.7525 ,             b4=t4/2,b4p,alpha,mass
* 1,0,1,1.0,                        force option:ex,ls,cm,CHG
* 5,  285.209,  291.384, -0.9       pairing:ipair,ggp,ggn,pairstab,rho0pr
* 0, 23, 49, 0.95, 0.95, 0, 0, 0,   grid: SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP
* 0.2, 0.0, 0.7,0.7, 11,11,2         init: BETAIN, OCTINI,  DEOCC(2),MMAX(2),isTSym
* 1,  0.0,  80.0, 20.0, 3,1,        CHF: swicon,bmin,bmax,bw,lqcon,ifresh
* 1.0,1.0,1.0,                      thermal: tempmn,tempmx,tempw
* 200,2000,200,0,1000000,           print opt.: IPRINT,IPLOT,ITMAX,ITEST,ITMXT
* 50.0, 0.4, 0.3, 0.9, 0.5,         step iter:E0INV,X0INV,ADDrho,adpot,ADDOCC
*  1.E-6,  1.0E-2,                  limit iter:ERRENR,ERhttp://www.golem.de/RVAR,
* 9999, 1,1.5,                      ext.step: IVS,ITSUB,X0FAC
* 1.E-6,1.E-6,3000,3000,100,11,     pair iter:ERRF,ERRA,FERMIP,FERMIN,NPITER,ITPINI
* 0.5,1.0, -2.1,                    pair iter:ADDPAR,EPSPAR,WFAC
* 1,1,1,1,                          Coul.BC:  ndrun,NQRUN,norun,NHRUN
* 1.D-5,20,6.0D-1,                  Coul.Iter:COLERR,IREPC,CMSsQR,
* 1.D-7 ,-50.0,   20,   1.0,        Yukawa Iter:YUKERR,SNL,IREPM,DMSFAC
* 0.1, 0.8, 0.0,  0.4,  999, 0,     CHF iter:QERR,C0,D0,QEPS,ITCONS(999),idampc
* 0,  10.0,  15.0,  0.03,           stab.CHF:ifstab,stwid,stslop,rhocut,
* 1.5,0.30,    0,   0.05,           damp.CHF:DALP,DGAM,IFLIN,SWIFAC,
*     0,                            switch type force IPARAM: SkM,SkM(NoLs),ext.osc, m*/m
*  0,0,0,0,0,0,0,0,0,               ireg(2),kblock(2),mblock(2),isoblo,istblo(2)
*   -1,,,,,,,,,,,,,,,               occup.selection:  LOCC(15)  (inactive)
*
*
*  The input parameters can be grouped in two categories,
*
*  I. parameters which must be specified by the user:
*  --------------------------------------------------
*
*     title      - a character string carrying the title of the run
*                  the first three digits are reserved for naming the
*                  output files. the name 'plt' is given if these are blank.
*
*     nprot      - proton number
*     nneut      - neutron number
*
*     t0..t3,x0..x3 ->  Skyrme parameters (only active for iparam=0)
*                       iparam=-1 switches to Woods-Saxon whose
*                       parameters are read in right after the 'iparam'
*     b4         - =t4/2 spin-orbit parameter
*     b4p        - new parameter for isotopic spin-orbit force
*                  old convention is recovered for b4p=0.0
*                  vanishing b4p is to be entered as a small b4p>1.E-20
*     alpha      - power of density dependence
*     mass       - hbar^2/2m
*
*     ifex       - switch for coulomb exchange in slater approximation
*     ifls       - switch for contribution from t1,t2 to spin-orbit force
*     ifcm       - switch for c.m. correction
*                   0 --> m_eff = m-m/A
*                   1 --> E_cm=<P^2>/2M approx. as 17.2MeV/A^0.2
*     chg        -  effective charge for direct field
*                   (use 1.0 for full Coulomb,
*                    0.0 for No Coulomb)
*
c     ipair      - type of pairing:
C                  -1 = filling approximation
C                   0 = constant Delta approximation
C                   1 = constant Force approximation (BCS)
C                   2 = constant Force approximation with
c                       s.p. correction (BCS)
C                   3= Lipkin-Nogami
C                   4 = thermal HF, temperature handled by "DELTP"
C                   5 = BCS with delta force
C                   6 = HFB with delta force
C                   7 = BCS with dens.dep. delta-force
C                   8 = HFB with dens.dep. delta-force
C                   9 = BCS-LN with delta force
C                  11 = BCS-LN with dens.dep. delta-force
C                 note: ipair=5,7,9 and 11 require available nr. of states
C                       high enough to cover 2*nr.of.particles
*     ggp,ggn    - gap or strength parameters for the BCS pairing
*                   0 --> defaults are used for ipair=0,1,2,3
*     pairstab   - stabilizing parameter for SBCS (0.0=BCS, has to be <0)
*                  pairstab= - E_cutp^2 (paper)  
*     b_strength - strength of external B-field 
*     sym        - sym=1.0 for calculations which are
*                  reflection symmetric about z=0.0
*                  otherwise sym=0.0
*     nz         - no. of z direction mesh points used in the
*                  calculation. note nz < nzd.
*                  'nz=0' calls for default setting according to 'nmass'
*     nr         - no. of r direction mesh points used in the
*                  calculation. note nr < nrd.
*                  'nr=0' calls for default setting according to 'nmass'
*     dz         - mesh spacing in the z direction.
*     dr         - mesh spacing in the r direction.
*     icutr      - momentum comp. cutoff in r-direction, ! set 0 !
*     icutz      - momentum comp. cutoff in z-direction  ! set 0 !
*     idrtyp     - switch for type of derivative
*                  0 -> fourier definition in z&r direction ;
*                  1 -> fourier definition z direction, but
*                       finite diff. of 11. order in r direction;
*                  2 -> fourier definition r direction, but
*                       finite diff. of 11. order in z direction;
*                  3 -> finite diff. 11. order in r&z direction.
*     isTSym     - changing the functional
*                  0 -> half spectra, only for even nuclei
*                  1 -> full spectra, only time-even densities
*                  2 -> full spectra, minimal galilaen-invariant func
*                  3 -> full spectra, minimal + spin-densities
*                  4 -> full spectra, complete functional
*
*     betain     - deformaion parameter for the initial axial
*                  harmonic oscillator wavefunctions,
*                  for 'betain=0' a default 'betain' is deduced
*                  from the actual deformation
*     octini     - initial kick towards octupole deformation
*                  deforms with 'octini'*Q_30/<r>^3;
*                  use small values to avoid curious nods in the w.f.
*                  for 'octini=0' a default 'octini' is deduced
*                  from the actual deformation
*     deocc(1/2) - addition to the first guess of 'eocc' with
*                  'eocc = 1.4*a**1/3 + deocc'. where 'eocc(1/2)'
*                  is the cutoff for the number of initial
*                  oscillator shells. 'deocc' can be near 0 for
*                  magic nucleon number. somewhat larger values
*                  are recommended far away from magic shells.
*                  A value of 0.8-1.2 is recommended for ipair=5,7,9 and 11.
*     mmax(1/2)  - maximum value for angular momentum 'm' in the
*                  setup phase. 'mm(1)' stands for protons and
*                  'mm(2)' for neutrons.
*
*     swicon      - switch for the constraint quadrupole
*                   1 --> constraint,  0 --> no constraint
*     betmin      - minimum deformation
*                   the given deformations are also used for initialization
*                   if 'betain=0' or 'octini=0' indicate this default choice,
*                   note that 'betmin' allows to transfer initialization
*                   even in case of unconstrained iteration.
*     betmax      - maximum deformation
*     betaw       - step in deformation
*                   these entries use dimensionless deformations
*                   of quadrupole or octupole type depending on 'lqcon'
*     lqcon       - type of constraint
*                    2 --> quadrupole
*                    3 --> octupole (only for ISYM=0)
*                    4 --> hexadacopole
*     ifresh      - switch for refresh of inital w.f. at each new deformation
*                    0 --> continue next deformation from old w.f.
*                    1 --> start next deformation with new w.f. from setup
*
*     tempmn      - minimum temperatur in T-loop for ipair=4
*     tempmn      - maximum temperatur in T-loop for ipair=4
*     tempw       - step in temperatur in T-loop for ipair=4
*
*  II. parameters which should only be changed by the experienced user:
*  --------------------------------------------------------------------
*
*     iprint     - modulus for printing output
*     itmax      - maximum no. of static iterations for one given
*                  quadrupole momentum
*     iplot      - modulus for printer plot of the density contours
*     itest      - test level = level test prints
*     itpfix     - iteration after which pairing is fixed
*     itmxt      - maximum number of all iterations (over all
*                  deformations)
*
*     e0inv      - energy cutoff for damping operator
*     x0inv      - expansion parameter for dpa step
*     addrho     - weight for adding new density
*     addpot     - weight for adding new potential
*                    'addrho' and 'addpot' have about the same effect,
*                    slowing down the propagation of the mean field.
*                    'addpot' is to be preferred if one aims top trace
*                    the energy convergence using 'subroutine energ' as
*                    this gives always the actual energy.
*     addocc     - weight for adding new occupation numbers from
*                  the pairing iteration.
*
*     errenr     - first convergence criterion:
*                  change in energy must come below 'errenr'
*     errvar     - second convergence criterion:
*                  the variance must come below 'errvar'
*
*     ivs        - level are reordered (and reoccupied) according
*                  to their s.p. energies as long as
*                      iter < ivs  .
*                  no reordering is done for
*                      iter > ivs.
*     itsub      - number of sub-iterations without orthogonalization  !mod 5/93
*                  presently active are:
*                   1 --> standard iteration with one 'hpsi' per iteration
*                   2 --> two damped steps per iteration, second step
*                         is weighted with 'x0fac'
*                  -2 --> second step is undamped and has stepsize x0fac
*     x0fac      - factor on 'x0inv' at even subiterations             !mod 5/93
*
*     errf        - allowed error in the fermi energies
*     erra        - allowed error in the particle no.
*     fermip      - upper limit for the proton fermi
*                   energy
*     fermin      - upper limit for the neutron fermi
*                   energy
*     npiter      - limit on the no. of iterations for
*                   self-consistent solutions to the
*                   gap equations
*     itpini      - nr. of iterations BEFORE pairing is switched on
*
*     addpar      -  mixing of new and old occ. numbers
*     epspar      -  damping of the paring step
*     wfac        -  pairing window: wfac*41/a**(1/3)
*
*     ndrun      - switch for dipole field in coulomb solver
*     nqrun      - switch for quadrupole field in coulomb solver
*     norun      - switch for octupole field in coulomb solver
*     nhrun      - switch for hexadecoupole field in coulomb solver
*
*     colerr     - level of precision for coulomb solver
*     irepc      - maximum number of iterations for coulomb solver
*     cmssqr     - damping mass (inverse length) for coulomb solver
*
*     yukerr     - level of precision for yukawa folding
*     snl        - damping mass (inverse length) for yukawa folding
*     irepm      - maximum number of iterations for yukawa folding
*
*     qerr        - accuracy for quadrupol deformation
*     c0,d0       - convergence-parameters for constraint
*     qeps        - effective stepsize in quadrupole iteration
*     itcons      - maximum number of constraint steps within iteration.
*                   the iteration returns to the standard step with
*                   the now fixed constraining force 'qlambd*qfac(iz,ir)'.
*                   a consistent quadrupole-constrained iteration is
*                   only obtained for 'itcons>itmax'.
*                   The 'itcons>0' serves also as a switch to iteration with
*                   fixed constraining force 'qlam': for 'itcons=0', the
*                   'qmin,qmax,qw' are deduced from 'betmin,betmax,betaw'
*                   and are used as the values for 'qlam'.
*                   The 'itcons<0' is used to switch to unconstrained 
*                   iteration after |itcons| inital constrained steps.
*                   This serves to drice pure HF iteration near a 
*                   desired deformation minimum.
*                   Thus use a very large 'itcons' to be on the safe side!
*                   ?? note: this feature may be removed some time later.
*     idampc      - option to switch the damping of the cnstraint
*                   0 --> old form with Woods-Saxon hat
*                   1 --> new dynamical damping in 'damcst'
*
*     ifstab      - switch for stabilizing external potentianl,
*                   recommended for high temeperatures and/or
*                   exotic nuclei
*     stwid       - parameter for the energy-width of the states included
*                   in the generating density
*     stslop      - slope of the stabilizing potential
*     rhocut      - 'rhocut*rho0' is the density from which on the
*                   asymptotic evaluation of the additional potential
*                   starts to count
*
*     dalp        - radius of the woods-saxon which damps the quadrupole
*                   operator at large distances
*     dgam        - width  of the woods-saxon which damps the quadrupole
*                   operator at large distances
*     iflin       - switch to linear asymptotic constraint
*     swifac      - factor on switch in z-direction,
*                   switch from quadratic to linear shape comes at
*                     swifac*1.2*A**1/3
*
*     iparam - number of the parametriztion used
*              0 -> direct input 't0,...,power' from for005
*              1,...,8 -> select predefined forces in subroutine 'param'
*              -1 -> deformed WS 
*
*     ireg(1..2)  - switch to printing of levels in Nilsson scheme
*       kblock   =  2*K values of blocked s.p. state                  *
*       mblock   =  M values of blocked s.p. state                    *
*       isoblo   =  isospin of blocked s.p. state                     *
*       istblo   =  number of blocked s.p. state in level sequence    *
*     locc        - explicit preoccupation
*                   --> presently disabled, use -1,,,,
*
*
*
***********************************************************************
*
*
c
      include 'paramr.inc'
*
      character*20 set
      character*3  outnam

c      integer blist   ! blocking

c      integer nneut, nprot               !implicit
c      integer iq
      
c      logical mothordaug
      
      integer decayjter,decayiter ! for beta decay


      logical blnprot       ! for blocking 
      logical blnneut

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

      include 'blocking.inc'
      LOGICAL, PARAMETER :: tglan=.FALSE.
      INTEGER, PARAMETER :: knumstateblock=200
      INTEGER :: istateblocking(2,knumstateblock)

c
cMPI      include 'mpir.inc'
c
      data iallop/0/                    ! protocol of opening files
      data itest/0/
c      logical tafter                    ! switch to "afterburn"  !also defined in iteration-loop -> not necessary
c      data tafter/.true./
!     set defaults
      tpairspace=.FALSE.         ! if TRUE then v**2 space, else uv space
      suv_wanted=0D0
      suv_mue=0D0

*
*
***********************************************************************
*                                                                     *
*     the prologue: print version number, set constants,              *
*                   and check parameters.
*                                                                     *
***********************************************************************
*
c      ifcran  = 1                    ! switch for cranking masses
C      swdiag  = 0                    ! switch for pre-diagonalisation
c
c      open(5,status='old',form='formatted',file='for005')
c      open(06,status='unknown',form='formatted',file='for006')
c
c
cMPI      call MPI_INIT(mpierror)
cMPI      call MPI_COMM_SIZE(MPI_COMM_WORLD,mpisize,mpierror)
cMPI      call MPI_COMM_RANK(MPI_COMM_WORLD,mpirank,mpierror)
c
      if(nzd.gt.nmx) stop ' wrong "nmx" in paramr.inc, check "nzd"!'
      if(nrd.gt.nmx) stop ' wrong "nmx" in paramr.inc, check "nrd"!'
*
*
***********************************************************************
*                                                                     *
*     outer loop: start each time a new run                           *
*                                                                     *
***********************************************************************
*
      do 1009 iread=1,99
      read(5,*,err=1010,end=1010) itmxt
      if(itmxt.le.0) goto 1010                       ! terminate ITMXT<0
c
      call const                                     ! set constants
*
*
***********************************************************************
*                                                                     *
*        read the parameters.                                         *
*                                                                     *
***********************************************************************
*
      read(5,*,err=991,end=991)  title
c
      if(title(1:3).eq.'   ') then
        outnam = 'plt'
      else
        outnam = title(1:3)
      endif
c
      open(7,status='unknown',form='formatted',file='for007.'//outnam)
      write(7,'(10x,a,i5)') ' no. of total iterations=',itmxt
      write(7,'(5x,2a)')
     &    ' Skyax  - the axial Skyrme-Hartree-Fock Code ',
     &    '          >>> version 2.1.0-odd <<<  '
      if(ifoccs.eq.1) then
        write(7,'(a)') ' WARNING: running IFOCCS=1 is a test version !'
        write(6,'(a)') ' WARNING: running IFOCCS=1 is a test version !'
      endif
c
      read(5,*,err=992,end=992)  nprot,nneut
      nmass  = nprot+nneut
      nchg   = nprot
      read(5,*,err=993,end=993)  t0,t1,t2,t3
      read(5,*,err=993,end=993)  x0,x1,x2,x3
      h2m(2) = 0D0 !!!!
      dglan = 0D0
      dglanp = 0D0
      IF(tglan) THEN
        read(5,*,err=993,end=993)  b4,b4p,power,h2m,dglan,dglanp
      ELSE
        read(5,*,err=993,end=993)  b4,b4p,power,h2m
      END IF
      WRITE(6,'(a,6(1pg13.5))') 'b4,b4p,power,h2m,dglan,dglanp=',
     &        b4,b4p,power,h2m,dglan,dglanp
      if(h2m(2).eq.0D0) h2m(2) = h2m(1) !!!!
      ioptsk = 0
      icoulfold = 0
      ifrso = -1      ! set the model for l*s contrib. to charge density
      read(5,*,err=993,end=993)  chgx,ifls,ifzpe,CHG,b_strength,ioptsk,
     &   icoulfold,ifrso
      WRITE(6,'(a,f4.1,2i3,2f10.5,2i4)') 'chgx...=',
     &   chgx,ifls,ifzpe,CHG,b_strength,ioptsk,ifrso
!      WRITE(*,*) 'chgx,ifls,ifzpe,CHG,b_strength,ioptsk=',
!     &     chgx,ifls,ifzpe,CHG,b_strength,ioptsk
      cutfac = 1.65D0                ! default for Skyrme forces
      IF(ioptsk==13) THEN            ! Fayans functional
        cutfac = -15.0D0                ! default for Fayans functionals
        tpairFy=.TRUE.
        READ(5,*) avp,h1vp,h2vp     
        WRITE(6,*) avp,h1vp,h2vp     
        READ(5,*) avm,h1vm,h2vm
        WRITE(6,*) avm,h1vm,h2vm
        READ(5,*) sigmavol,rho0fy
        WRITE(6,*) sigmavol,rho0fy
        READ(5,*) asp,hsp,hgradp 
        WRITE(6,*) asp,hsp,hgradp 
        READ(5,*) asm,hsm,hgradm 
        WRITE(6,*) asm,hsm,hgradm 
        READ(5,*) hcoul,sigmacoul
        WRITE(6,*) hcoul,sigmacoul
        READ(5,*) fxi,h1xi,hgradxi,gammaxi,rho0xi
        WRITE(6,*) fxi,h1xi,hgradxi,gammaxi,rho0xi
        READ(5,*) delfxi,ecstab,cutfac !,isovectormodel
        WRITE(6,*) delfxi,ecstab,cutfac,isovectormodel
!         delfxi=0D0
        IF(isovectormodel==0) THEN
          fxip=fxi
          h1xip=h1xi
          hgradxip=hgradxi
          fxin=fxi
          h1xin=h1xi
          hgradxin=hgradxi
        ELSE IF(isovectormodel==2) THEN
          fxip=fxi+delfxi
          h1xip=h1xi
          hgradxip=hgradxi
          fxin=fxi-delfxi
          h1xin=h1xi
          hgradxin=hgradxi
        ELSE IF(isovectormodel==1) THEN
          xipfac=(fxi+delfxi)/fxi
          fxip=fxi*xipfac
          h1xip=h1xi*xipfac
          hgradxip=hgradxi*xipfac
          xinfac=(fxi-delfxi)/fxi
          fxin=fxi*xinfac
          h1xin=h1xi*xinfac
          hgradxin=hgradxi*xinfac
        END IF
        READ(5,*) c0tau,c1tau 
        WRITE(6,*) c0tau,c1tau 
        READ(5,*) c0nabj,c1nabj,c0tens,c1tens
        WRITE(6,*) c0nabj,c1nabj,c0tens,c1tens
        CALL ini_fyfunc()
      END IF

      read(5,*,err=993,end=993)  ipairin,ggp,ggn,pairstab,rho0pr,cutfac
      WRITE(6,'(i5,2(1pg14.6),0pf5.2,1pg14.6,0pf5.1)') 
     &          ipairin,ggp,ggn,pairstab,rho0pr,cutfac
      IF(ioptsk==13) pairstab=-ecstab
      IF(ipairin==6 .OR. ipairin==8) THEN
        thfb = .TRUE.
        ipair = ipairin-1
        WRITE(7,*) ' now HFB:',ipair,thfb
      ELSE
        thfb = .FALSE.
        ipair = ipairin
      END IF
      read(5,*,err=997,end=997)  SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP
      WRITE(6,*)  'GRID:',SYM,NR,NZ,DR,DZ,ICUTR,ICUTZ,IDRTYP
      isym = sym+0.0000001
      read(5,*,err=994,end=994)  BETAIN,OCTINI,DEOCC,MMAX,isTSym
      WRITE(6,*)  'INIT:', BETAIN,OCTINI,DEOCC,MMAX,isTSym
c      read(5,*,err=994,end=994)  distin,npextra,nnextra
      qlambquad = 0D0
      read(5,*,err=995,end=995)  
     &       swicon,betmin,betmax,betaw,lqcon,ifresh,ifcran
     &       ,qlambquad       ! for additional quadr. constr.
      WRITE(6,*) 'CONSTR:',
     &       swicon,betmin,betmax,betaw,lqcon,ifresh,ifcran
     &       ,qlambquad       ! for additional quadr. constr.
      read(5,*,err=995,end=995)  tempmn,tempmx,tempw,suv_mue
      WRITE(6,*)   'TEMP:',tempmn,tempmx,tempw,suv_mue
      read(5,*,err=996,end=996)  IPRINT,IPLOT,ITMAX,ITPFIX
      WRITE(6,*)   'PLOT:',IPRINT,IPLOT,ITMAX,ITPFIX
      read(5,*,err=998,end=998)  E0INV,X0INV,ADDrho,addpot,ADDOCC
      read(5,*,err=998,end=998)  ERRENR,ERRVAR
      read(5,*,err=998,end=998)  IVS,ITSUB,X0FAC
      read(5,*,err=998,end=998)  ERRF,ERRA,FERMIP,FERMIN,NPITER,ITPINI
      read(5,*,err=998,end=998)  ADDPAR,EPSPAR,WFAC
      read(5,*,err=998,end=998)  ndrun,NQRUN,norun,NHRUN
      read(5,*,err=998,end=998)  COLERR,IREPC,CMSsQR
      read(5,*,err=998,end=998)  YUKERR,SNL,IREPM,DMSFAC
      read(5,*,err=998,end=998)  QERR,C0,D0,QEPS,ITCONS,idampc
      read(5,*,err=999,end=999)  ifstab,stwid,stslop,rhocut
      read(5,*,err=999,end=999)  DALP,DGAM,IFLIN,SWIFAC
      read(5,*,err=999,end=999)  iparam
      if(iparam.EQ.-1) then
        do iq=1,2
          read(5,*,err=999,end=999) wsrad0(iq),ws0(iq),wskap(iq),
     &       wsalfa,wsb0(iq),wsb2(iq),wsb4(iq)
          wswid(iq) = one/wsalfa
        enddo  
      endif
      read(5,*,err=999,end=999)  ireg,kblock,mblock,isoblo,istblo !,spinblo
      read(5,*,err=999,end=999)  LOCC
      WRITE(*,*) 'LOCC=',locc
      IF(SUM(kblock).NE.0) THEN
        IF(locc(1)==-1) THEN
          numstateblock=NINT((betmax-betmin)/betaw)+1
          DO i=1,numstateblock
            IF(i>knumstateblock) STOP 'too many deformations'
            istateblocking(1,i)=kblock(isoblo)
            istateblocking(2,i)=istblo(isoblo)
          END DO
          WRITE(*,'(a,200i4)') 'istateblocking - kblock:',
     &        istateblocking(1,1:numstateblock)
          WRITE(*,'(a,200i4)') 'istateblocking - istblo:',
     &        istateblocking(2,1:numstateblock)
        ELSE
          DO i=1,knumstateblock
            READ(5,*,ERR=19,END=19) betact,istateblocking(:,i)
            WRITE(*,*) 'istateblocking=',i,istateblocking(:,i)
            numstateblock=i
          END DO
 19       continue
          WRITE(*,*) 'numstateblock:',
     &           numstateblock,NINT((betmax-betmin)/betaw)+1
          IF( NINT((betmax-betmin)/betaw)+1 .NE. numstateblock) 
     &      STOP "supply of blocked states does not match deformations"
        END IF
      END IF
      if(max(mmax(1),mmax(2)).gt.kmpsi+1)
     & stop ' MMAX too large. enhance KPSI in PARAMR.INC!'
      if(ifcran.eq.2 .and. swicon.gt.oh) 
     &  stop ' Separable RPA based on ground state calculation '
      if(itpfix.lt.10) itpfix=itmax+1
      if(iparam.EQ.-1 .AND. swicon.GT.oh) 
     &  stop ' WS potential (IPARAM=-1) excludes CHF'

      istblototal = SUM(istblo)                                ! to memorize initial setting

      write(7,*) 'ireg etc:',
     &   ireg,kblock,mblock,isoblo,istblo,iparam,SWIFAC

c      read(5,'(1x)')
c
c      if(ifls.eq.1)
c     &  stop ' spin-orbit contribution from T1,T2 not yet active'
      t4    =  2.0*b4
      call param(set)                                ! force parameters
      
      if(b4p.eq.zero) b4p = t4*oh


!      istblo = 0
      istcblock = 0 ! time reversed state of the blocked state

c TEST increase pairing strength
c      ggp = ggp*1.2 
c TEST ENDE

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
c     Loop for calculation of the matrix element of a beta-decay              c 
c     or to test the overlap of the GOA                                       c
c                                                                             c
c                                                                             c 
c     to make an iteration for the initial and the final core                 c
c     and to calculate the overlap of the single particle wave functions      c
c     afterwards                                                              c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


c
c     convert deformation, using the relations:
c       Q_cartesian = 1.08 * A**5/3  * beta_2
c       Q_40        = 0.495* A**7/3  * beta_4
c       Q_30        = 0.41 * A**2    * beta_3
c       Q_20        = 0.34 * A**5/3  * beta_2
c
      amass = nmass
      if(lqcon.eq.2) then
        cfact  = 1.09 * amass**1.666666667
      elseif(lqcon.eq.3) then
        cfact  = 0.41 * amass**2
      else
        cfact  = 0.495 * amass**2.3333333333
      endif
      qmin   = betmin*cfact
      qmax   = betmax*cfact
      qw     = betaw *cfact
c
c     automatic setting of grid size
c
      if(nz.eq.0) then
        if(nr.eq.0) then
          reff  = 1.2*amass**0.3333333
          if(lqcon.eq.2) then
            b2z   = max(betmin,betmax)
            b2r   = min(betmin,betmax)
            b3max = abs(octini)
          elseif(lqcon.eq.3) then
            b2z   = max(betain,zero)
            b2r   = max(-betain,zero)
            b3max = max(abs(betmin),abs(betmax))
          else
            b2z   = max(betain,zero)
            b2r   = max(-betain,zero)
            b3max = abs(octini)
          endif
        else
          stop ' NZ, NR must be both 0 or both nonzero'
        endif
        reffz  = reff*(1.0+b2z+0.75*b3max)
        reffr  = reff*(1.0-b2r/3.0+0.65*b3max)
!        nz     = (2*reffz+6)/dz
        nz     = INT((reffz+6.4D0)/dz*0.8736)!+1     ! factor 0.8736 for 2D vs. 1D
        if(sym.eq.zero) nz = 2*nz+1
!        nr     = (2*reffr+6)/dr
        nr     = INT((reffr+6.4)/dr*0.8736)!+1
        write(6,'(a/6f8.2,2i5)')
     &  ' b2z,b2r,b3max,reff,reffz,reffr,nz,nr=',
     &    b2z,b2r,b3max,reff,reffz,reffr,nz,nr
      else
        if(nr.eq.0)
     &     stop ' NZ, NR must be both 0 or both nonzero'
      endif
c
c     check consistency of input
c
      if(octini.ne.0.0D0) then
        if(swicon.gt.oh .and. lqcon.eq.3) then
          octini = 0.0
          write(6,'(a)')
     &    ' constraint calc. --> initial octupole kick deactivated'
c        else
c          lqcon = 3
        endif
      endif
      if(lqcon.ne.2 .and. lqcon.ne.3 .and. lqcon.ne.4)
     &  stop ' constraint must be L=2,3 or 4 --> change LQCON'
      if(lqcon.eq.3 .and. isym.eq.1)
     &  stop ' octupole constraint requires ISYM=0'
      if((lqcon.eq.3.or.lqcon.eq.4) .and. iflin.eq.1)
     &  stop ' linear asympt. constraint not enabled with L=3 or 4'
      if(icutr.gt.0)
     &  stop '  icutr>0  presently disabled. ---> stop'
      if(icutz.gt.0)
     &  stop '  icutz>0  presently disabled. ---> stop'

*
***********************************************************************
*                                                                     *
*        write the parameters and options                             *
*                                                                     *
***********************************************************************
*
      write(7,5001) title,outnam,nmass,nchg,
     &              nz,nzd1,nr,nrd1,dz,dr,mmax,icutr,icutz
5001  format('1',//,10x,a,///,
     &               10x,' appendix of output files=',a,/,
     &               10x,' mass number=',i5,/,
     &               10x,'total charge=',i5,//,
     x               10x,'mesh dimensions',//,
     x               10x,'nz =',i5,'  maximum allowed =',i5,/,
     x               10x,'nr =',i5,'  maximum allowed =',i5,//,
     x               10x,'mesh spacings  ',//,
     x               10x,'dz =',f5.3,'  fm',/,
     x               10x,'dr =',f5.3,'  fm',
     x               10x,'mmax =',2i5,/,
     x               10x,'icutr,icutz=',2i5)
        if(idrtyp.eq.0) then
         write(7,'(10x,a,2i5)')
     &      ' d/dr & d/dz from fourier.'
        elseif(idrtyp.eq.1) then
         write(7,'(10x,a,2i5)')
     &      ' d/dr from finite difference & d/dz from fourier.'
        elseif(idrtyp.eq.2) then
         write(7,'(10x,a,2i5)')
     &      ' d/dr from fourier & d/dz from finite difference.'
        elseif(idrtyp.eq.3) then
         write(7,'(10x,a,2i5)')
     &      ' d/dr & d/dz from finite difference.'
        else
         stop ' invalid "idrtyp"'
        endif
*
      IF(ioptsk.NE.13) THEN
        write(7,5002) set,
     &   t0,t1,t2,t3,t4,x0,x1,x2,x3,b4p,power,e2,h2m,ifzpe,ifls,ifrso
5002    format(//,10x,a,//,
     x           10x,'t0     =',f10.3,'  mev*fm^3',6x,/,
     x           10x,'t1     =',f10.3,'  mev*fm^5',6x,/,
     x           10x,'t2     =',f10.3,'  mev*fm^5',6x,/,
     x           10x,'t3     =',f10.3,'  mev*fm^3+3*power',6x,/,
     x           10x,'t4     =',f10.3,'  mev*fm^4',6x,/,
     x           10x,'x0     =',f10.6,/,
     x           10x,'x1     =',f10.6,/,
     x           10x,'x2     =',f10.6,/,
     x           10x,'x3     =',f10.6,/,
     x           10x,'b4p    =',f10.6,/,
     x           10x,'power  =',f10.6,/,
     x           10x,'e2     =',f10.3,/,
     x           10x,'h2m    =',2f10.3,/,            !!!!
     x           10x,'ifzpe  =',i5,/,
     x           10x,'ifls   =',i5,/,
     x           10x,'ifrso  =',i5)
        IF(b_strength.NE.0D0) WRITE(7,'(10x,a,f7.3)') 
     &    'B_strength=',b_strength
      ELSE
        WRITE(7,'(5(a,3(1pg13.5)/),5(a,2(1pg13.5)/),a,i3)')
     &     'avp,h1vp,h2vp=   ',avp,h1vp,h2vp,
     &     'avm,h1vm,h2vm=   ',avm,h1vm,h2vm,
     &     'asp,hsp,hgrad=   ',asp,hsp,hgradp ,
     &     'asm,hsm,hgrad=   ',asm,hsm,hgradm ,
     &     'fxi,h1xi,hgradxi=',fxi,h1xi,hgradxi,
     &     'gammaxi,rho0xi=  ',gammaxi,rho0xi,
     &     'sigmavol,rho0fy= ',sigmavol,rho0fy,
     &     'c0nabj,c1nabj=   ',c0nabj,c1nabj,
     &     'hcoul,sigma=     ',hcoul,sigmacoul,
     &     'delfxi,ecstab=   ',delfxi,ecstab,
     &     'isovectormodel=  ',isovectormodel
        WRITE(7,'(a,2(1pg13.5))') 'fxip,fxin=        ',fxip,fxin
        WRITE(7,'(a,2(1pg13.5))') 'h1xip,h1xin=      ',h1xip,h1xin
        WRITE(7,'(a,2(1pg13.5))') 'hgradxip,hgradxin=',hgradxip,hgradxin
        WRITE(7,'(a,2f10.3,3(/a,i5))')
     x     'h2m=             ',h2m,
     x     'ifzpe=           ',ifzpe,
     &     'ifls=            ',ifls,
     &     'ifrso=           ',ifrso
      END IF
      write(7,5006) yukerr,snl,irepm
5006  format(//,10x,'meson equation parameters',//,
     x               10x,'yukerr =',g10.2,'  error for yuk iter',/,
     x               10x,'snl    =',g10.2,'  damp mass for nonlin.',/,
     x               10x,'irepm  =',i10,'  no. of iterations   ')
*
      write(7,5017) colerr,irepc,ndrun,nqrun,norun,nhrun,cmssqr
5017  format(//,10x,'laplace equation parameters',//,
     x   10x,'colerr = ',g10.2,'  error for col iter',/,
     x   10x,'irepc  = ',i10,'  no. of iterations  ',/,
     x   10x,'ndrun  = ',i10,'  dipole bound (y=1/n=0)  ',/,
     x   10x,'nqrun  = ',i10,'  quad bound (y=1/n=0)  ',/,
     x   10x,'norun  = ',i10,'  oct bound (y=1/n=0)  ',/,
     x   10x,'nhrun  = ',i10,'  hex bound (y=1/n=0)  ',/,
     x   10x,'cmssqr = ',g10.2,'  mass for col damping')
      if(sym.eq.one) then
       ndrun = 0
       norun = 0
      else
       ndrun = 1
       norun = 1
      endif
*
      write(7,5007) chg,chgx,e0inv,x0inv,addrho,addpot,addocc,ivs
5007  format(/,
     x10x,'chg    =',f5.1,'  effective charge for direct field',/,
     x10x,'chgx   =',f5.1,'  effective charge for exchange field',///,
     x10x,'pade damping parameters',//,
     x10x,'e0inv   =',f10.3,'  damping energy scale (mev)',/,
     x10x,'x0inv   =',g10.3,'  expansion parameter',/,
     x10x,'addrho  =',g10.3,'  weight for accumulating density',/,
     x10x,'addpot  =',g10.3,'  weight for accumulating potential',/,
     x10x,'addocc  =',g10.3,'  weight for accumulating occ.numbers',/,
     x10x,'ivs     =',i10  ,'  up to iter=ivs: level ordering ' )
*
      write(7,5004) ipairin,npiter,itpini,ggp,ggn,pairstab,rho0pr,
     x              errf,erra,fermip,fermin,wfac,cutfac,tpairspace
 5004 format('1',//,10x,'bcs pairing model parameters',//,
     x 10x,'ipair   =',i10,'  ipair=-1 nopairing',/,
     x 29x,                '  ipair=0 constant gap model',/,
     x 29x,                '  ipair=1 constant force model',/,
     x 29x,                '  ipair=2 constant force BCS',/,
     x 29x,                '  ipair=3 Lipkin-Nogami',/,
     x 29x,                '  ipair=4 temperature CHF',/,
     x 29x,                '  ipair=5 BCS delta-force',/,
     x 29x,                '  ipair=7 BCS dens.dep. delta-force',//,
     x 29x,                '  ipair=9 BCS-LN delta-force',/,
     x 29x,                '  ipair=11 BCS-LN dens.dep. delta-force',//,
     x 10x,'npiter  =',i10,'  no. of pairing iterations',//,
     x 10x,'itpini  =',i10,'  no. iterations before pairing',//,
     x 10x,'ggp     =',f10.3,'  gap or strength for p',' (mev)',/,
     x 10x,'ggn     =',f10.3,'  gap or strength for n',' (mev)',/,
     x 10x,'           if ggn,ggp eq. 0 => standard model',/,
     x 10x,'pairstab=',f10.3,'  BCS stabilizing parameter',//,
     x 10x,'rho0pr  =',f10.3,'  density parameter in DDDI',//,
     x 10x,'errf    =',f10.3,'  error in fermi energy',/,
     x 10x,'erra    =',f10.3,'  error in mass',//,
     x 10x,'fermip  =',f10.3,'  approximate fermi energy p',' (mev)',/,
     x 10x,'fermin  =',f10.3,'  approximate fermi energy n',' (mev)',/,
     x 10x,'wfac    =',f10.3,' pairing window: wfac*41/a**(1/3)  ',/,
     x 10x,'cutfac  =',f10.3,' pairing cutoff   ',/,
     x 10x,'v2 space=',l10,' v**2 pairing space (else uv space)')
       write(7,5014)    addpar,epspar
 5014  format(
     x 10x,'addpar  =',f10.3,' mixing of new and old occ. numbers',/,
     x 10x,'epspar  =',f10.3,' damping of the paring step ')
c
      if(nr.gt.nrd1) then
         write(7,6001) nr,nrd
6001     format(' insufficient r dimension - nr =',i5,' nrd = ',i5)
         stop 'dimension in setup'
      end if
c
      if(nz.gt.nzd1) then
         write(7,6002) nz,nzd
6002     format(' insufficient z dimension - nz =',i5,' nzd = ',i5)
         stop 'dimension in setup'
      end if
c
      if((ipair.eq.5.or.ipair.eq.7)) then
         if (pairstab.gt.0.0D0) then
            write(7,6003) pairstab
 6003    format(' BCS stab. parameter must not be positive - pairstab=',
     &          f10.3)
         stop 'wrong setup for stabilized BCS'
         endif
      else if (pairstab.ne.0.0D0) then
         write(7,6004)
 6004    format(' stabilized BCS only allowed for ipair=5,7')
         stop 'wrong setup for stabilized BCS'
      endif
c
      if(itcons.ge.0) then
        iflmb0 = 0
      else
        iflmb0 = 1
      endif
      itcons = abs(itcons)
*
      write(7,5008) nmass,nchg,deocc,betain
5008  format('1',///,10x,'static input parameters',//,
     x10x,'nmass   =',i10,'  atomic mass number',/,
     x10x,'nchg    =',i10,'  atomic charge number',/,
     x10x,'deocc(1)=',f10.3,'  additional no. of quanta for proton ',
     x                      '  shell ',/,
     x10x,'deocc(2)=',f10.3,'  additional no. of quanta for neutron',
     x                      '  shell ',/,
     x10x,'betain  =',f10.3,'  deformation parameter ')
      if (locc(1).gt.0) then
      write(7,5088) (locc(i),i=1,10)
5088  format(/,10x,'locc(10)=',10(i3),' cut off')
      endif
      write(7,5009) err
5009  format(/,
     x10x,'err     =',g10.3,'  energy accuracy')
      write(7,5010) itmax,iprint,iplot,itest
5010  format(/,
     x10x,'itmax   =',i10,'  maximum no. of static iterations',/,
     x10x,'iprint  =',i10,'  iteration modulus for printing',/,
     x10x,'iplot   =',i10,'  iteration modulus for density plot'/,
     x10x,'itest   =',i10,'  level for testprints')
      write(7,'(/10x,a,f10.3,a)')
     x 'swicon =',swicon,'  switch for the constraint'
      if (swicon.gt.oh) then
      write(7,'(/10x,a)') 'constraint calculation:'
      write(7,5011)  qerr,c0,d0,qeps,itcons
c     write(7,5011)  qlamb
c011  format(/,
c    x10x,'qlamb   =',f12.8,'  parameter qlamb     ')
c    x10x,'qsoll   =',f10.5,'  constraint cartesian quadrupolmoment',/,
5011  format(/,
     x10x,'qerr    =',g10.3,'  accuracy for quadrupol deformation',/,
     x10x,'c0      =',f10.5,'  stepsize quadr.-correction',/,
     x10x,'d0      =',f10.5,'  parameter d0 for constraint calc  '/,
     x10x,'qeps    =',f10.5,'  stepsize qlambdacorrection ',/,
     x10x,'itcons  =',  i10,'  max. constraint iterations')
      write(7,5020)  betmin,betmax,betaw,lqcon
5020  format(/,
     x10x,'betmin =',f14.5,' minimal beta_soll     ',/,
     x10x,'betmax =',f14.5,' maximal beta_soll     ',/,
     x10x,'betaw  =',f14.5,' step of beta_soll',/,
     x10x,'lqcon  =',i14,  ' type of constraint')
      write(7,5012) dalp,dgam
5012  format(/,
     x10x,'dalp   =',f10.4,'  parameter det. the beg. of damping.',/,
     x10x,'dgam   =',f10.4,'  parameter det. the width of damp.  ',/)
      else
      write(7,'(/10x,a)') 'calculation without constraint'
      endif
*
***********************************************************************
*                                                                     *
*     open the necessary files                                        *
*                                                                     *
***********************************************************************
*
      if(iallop.eq.0) then
       open(18,status='unknown',form='formatted',file='fields.'//outnam)
       open(71,status='unknown',form='formatted',file='conver.'//outnam)
       write(71,'(a,a)') '  Plot convergence properties for:',title
       open(72,status='unknown',form='formatted',file='quadru.'//outnam)
       write(72,'(a,a)') ' Plot deformations for:',title
       if(itest.ge.1) then
       open(73,status='unknown',form='formatted',file='constc.'//outnam)
       open(74,status='unknown',form='formatted',file='constr.'//outnam)
       endif

	open(17,status='unknown',form='formatted',file='block.'//outnam)

c       if(klocal.eq.1) then
c       open(78,status='unknown',form='formatted',file='localz.'//outnam)
c       write(78,'(a,2i5)') '# plot localization for nprot,nneut=',npart
c       endif
       iallop = 1
       if(ifcran.eq.1) then
         open(75,status='unknown',
     &           form='formatted',file='crankn.'//outnam)
         write(75,'(a/a,2i4/a,i1,3a)')
     &   title(1:78),
     &   ' A,Z=',nmass,nchg,
     &  ' beta_',lqcon,'   SC_mass    Ing.mass   In.ms.-od  ',
     &  ' a.m. mass    a.m. wid  beta_wid    energy    cm_corr',
     &  '  PNmass     PNMas0    ZPE_triax   PNmas_a.m.'
         close(75)
c
         open(76,status='unknown',
     &           form='formatted',file='zpecor.'//outnam)
         write(76,'(a/a,2i4/a,i1,3a)')
     &   title(1:78),
     &   ' A,Z=',nmass,nchg,
     &  ' beta_',lqcon,'   ZPE_SC    ZPE_beta   ZPE_beta-od    ZPE_rot',
     &  '    a.m. wid   E_ZPE_cm   ZPE_triax  omega_tri  r**2_charge',
     &  ' beta_pr  rdms_C  surf_C  cm_mod'
         close(76)
c
c         if (lqcon.eq.2) then              !!!!
c            call ini_files(title,outnam)   !!!!
c         endif                             !!!!
c
         if(kintpl.eq.1) then
           open(77,status='unknown',
     &             form='unformatted',file='intpol.'//outnam)
         endif
       endif
      endif

*
***********************************************************************
*                                                                     *
*     initialisation of integration weights, pairing etc.             *
*                                                                     *
***********************************************************************
*
      if(ifzpe.eq.0) then
       h2meff(1)  = h2m(1)-h2m(1)/nmass  !!!!                    ! optionally ZPE by m*
       h2meff(2)  = h2m(2)-h2m(2)/nmass  !!!!
       write(7,'(a,f10.3)')
     &   ' ZPE with m-m/A, h2m_eff now=',h2meff
      else
       h2meff(1)  = h2m(1) !!!!
       h2meff(2)  = h2m(2) !!!!
       write(7,'(a)')
     &  ' ZPE with <P^2>/2M, only approximately as 17.3MeV/A^0.2 '
      endif
      h2mav = 0.5D0*(h2meff(1)+h2meff(2)) !!!!
!     scale Landau parameter g0
      call nucmat(rhonm,eovera,compr,effmas,asymme)
      write(6,'(a,5(1pg13.5))') 'NUCMAT:',
     &   rhonm,eovera,compr,effmas,asymme
      write(7,'(a,5(1pg13.5))') 'NUCMAT:',
     &   rhonm,eovera,compr,effmas,asymme
!      fermp2 = (3.*pi*pi/2.*rhonm)**(0.66666667)
!      hqqd2m = 0.5*(h2m(1)+h2m(2))
!      gfac   = sqrt(fermp2)*effmas/(pi*pi*hqqd2m)
      gfac   = effmas/150D0   !  first version: 75D0
      dglan0=dglan
      dglan = dglan0/gfac
      dglanp0=dglanp
      dglanp = dglanp0/gfac
      WRITE(7,'(a,6(1pg13.5))') 
     &  ' dG0,dG0before,dG0p,dG0pbgefore,gfac=',
     &   dglan,dglan0,dglanp,dglanp0,gfac


      call inicw(nmass,nchg)                       ! grid etc
      if(2*eocc(1).gt.mmax(1)+3.0)
     & stop ' MMAX too small for given DEOCC for protons'
      if(2*eocc(2).gt.mmax(2)+3.0)
     & stop ' MMAX too small for given DEOCC for neutrons'
      call inider               ! derivative matrices
      call idanuc               ! nucleon damping
      call idames               ! coulomb damping
*kr**********************
ckr      call inimul
*kr**********************
c
      if (swicon.lt.oh) then
c         qmin = zero
         qmax = qmin
      end if
c
      if ((isym.eq.0 .and. octini.ne.zero) ) then
       call inicon
      endif
      init77 = 0                        ! counter for initialization of 77



c
c     deciding if even-even or even-odd (odd-even)
c


      blnneut = (mod(nneut,2).eq.1)   
      blnprot = (mod(nprot,2).eq.1)

      
      if(blnneut.or.blnprot)then  ! for even-even only the positiv mpsi
c      isTSym = 0                 ! isTsym = 1,2,3,4 depending on functional
      else
!      isTSym = 0 
      endif   
      WRITE(7,*) 'istsym=',istsym



     


!      if(ifresh.eq.1) then
        call ssetup                     ! the wavefunctions
        call inipai(fermip,fermin,ggp,ggn,nmass)
        if(kintpl.eq.1 .and. init77.eq.0) then
           write(77)    title(1:78),nmass,nchg,nst
           write(7,*)   title(1:78),nmass,nchg,nst
           write(6,*)   title(1:78),nmass,nchg,nst
           init77 = 1
        endif
!      endif





!      if (lqcon.eq.2) then                  !!!!
            call ini_files(title,outnam)
!      endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      Part to calculate the blocked state for an odd nucleus with         c
c      the minimal total energy                                            c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      eold   = energy
      x0sav  = x0inv
      
      

cJE      mothordaug = (mod(decayiter,2).eq.1)         ! odd -> mother , even -> daughter


      
      if((blnneut.or.blnprot).and.(ipair.ne.-1))then   ! for HF -> filling approx.
c         swicon = 0.0                ! searching the optimal blocked without constraint   
         swdiag = 0.0                ! swdiag run into trouble due to the contribution of sigma (shell closure) 
c        swdiag = 1.0
      endif



!      if((blnneut.or.blnprot).and.(ipair.ne.-1) .AND.
!     &    SUM(istblo)==0 .AND. swicon < 0.5D0)then
      if((blnneut.or.blnprot).and.(ipair.ne.-1) .AND.
     &    SUM(istblo)==0 .AND. swicon<0.5D0
     &    .AND. SUM(suv_mue)==0D0)then
        qsoll = qmin
        WRITE(*,*) 'run BLOCK before Q loop'
        call block(blnprot,blnneut,outnam,x0sav)     ! main subroutine to calculate the optimal blocked state
        WRITE(*,*) 'returned from 1. BLOCK in SKYAX'
        IF(ttrack) THEN
          WRITE(7,*) 'returned from 1. BLOCK in SKYAX'
          WRITE(17,*) 'returned from 1. BLOCK in SKYAX'
        END IF
        STOP
      endif

      WRITE(*,*) ' after BLOCK'

      write(17,*)'SKYAX: modneut,modprot',mod(nneut,2),mod(nprot,2)
      IF(ttrack) WRITE(17,*) ' kblock etc=',kblock,mblock,istblo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
c     Here ends the part to calculate the optimal blocking state       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      if((blnneut.or.blnprot).and.(ipair.ne.-1))then 
        if(betmin.ne.betmax)then 
          swicon = 1.0      ! to calculate also PES of odd nuclei 
!          ifresh = 1          ! for a smooth PES
          qmin   = betmin*cfact
          qmax   = betmax*cfact
          qw     = betaw *cfact
          WRITE(*,*) 'qmin,qmax set'
        endif
      endif 



ccccccccccc



      itert  = 0                         ! set counter for total iterations
      if (swicon.eq.1D0) then              ! first constraint iteration 
         swdiag=0                        ! to stabilize MF at the
      else                               ! right deformation

        if(isTSym.eq.0 .AND. swicon<1D-1)then   
           swdiag=1
        else   
           swdiag=0                        ! ping pong 
        endif   

      endif
                                         ! without pre-diagonalisation

c      call wrt_spprop(outnam,1)

*
***********************************************************************
*                                                                     *
*     loop over deformations                                          *
*                                                                     *
***********************************************************************
*

      WRITE(*,*) ' befor Q-loop. SWICON=',swicon,qmin,qmax,qw
!      do 1100 qact=qmin,qmax,qw
!      DO qact=qmin,qmax,qw
      iqdmax=NINT((qmax-qmin)/qw)
      DO iqact=0,iqdmax
        qact = iqact*qw+qmin

        WRITE(*,*) 'Q-loop: qact=',qact,iqact,iqdmax

!       restore initial setting of blocking
        IF(istblototal==0) THEN
          kblock=0
          mblock=0
          istateblocking=0  
        END IF        

!       recover initial setting if non zero blocking
        IF(SUM(kblock).NE.0) THEN
          iact=NINT((qact-qmin)/qw)+1
          WRITE(*,*) 'iact etc:',
     &       iact,kblock(isoblo),istateblocking(1,iact)
          IF(kblock(isoblo).NE.istateblocking(1,iact)) 
     &      STOP "KBLOCK does not match blocking input"
          istblo(isoblo)=istateblocking(2,iact)
          WRITE(7,*) 'now running with istblo=',istblo(isoblo)
          WRITE(*,*) 'now running with istblo=',istblo(isoblo)
        END IF

      if(ifresh.eq.1 .or. qact.eq.qmin) then
ctest        if(lqcon.eq.2 .and. betain.eq.zero) betain = qact/cfact
        if(lqcon.eq.3 .and. octini.eq.zero) octini = qact/cfact


        IF(ttrack) THEN
          WRITE(17,*) '2. kblock etc=',kblock,mblock,istblo
          WRITE(*,*) '2. kblock etc=',kblock,mblock,istblo
        END IF
        if((blnneut.or.blnprot) .AND. SUM(istblo)==0 
     &     .AND. swicon<0.5D0)then
           call read_com(iter,outnam) ! for odd nuclei initialize with first calculation (filling aproximation)
          WRITE(*,*) ' previous calculation read'
        goto 988
        endif
       
        IF(ifresh==1)  THEN
          call ssetup       ! re-initialize wavefunctions
          call inipai(fermip,fermin,ggp,ggn,nmass)
        END IF

 988    continue
        IF(ttrack) THEN
          WRITE(17,*) '3. kblock etc=',kblock,mblock,istblo
          WRITE(*,*) '3. kblock etc=',kblock,mblock,istblo
        END IF


        call inipai(fermip,fermin,ggp,ggn,nmass)
        if(kintpl.eq.1 .and. init77.eq.0) then
           write(77)    title(1:78),nmass,nchg,nst
           init77 = 1
        endif
 

      endif
c
      if(itcons.eq.0) then
       qlamb = qact
       qsoll = betain * 1.09*(one*npart(1)+npart(2))**1.6666667  ! for 'inicms'
      else
       qsoll = qact
      endif
      betsol = qsoll/cfact
      if(itest.ge.1) then
       write(73,'(a,a)')
     & '   qtilc      qcst     varqsc      qq       ratio    qdh/qq  ',
     & '  qlcstp     ratio'
       write(74,'(a,a)')
     & '    qtil        q      varqsq      qq       ratio    qdh/qq  ',
     & '  qlstep     ratio'
      endif
c
c     initialisation of damping mass for coulomb solver and
c     constraint potential
c
      call inicms
      if (swicon.gt.oh) then
       call inicon
       do iwu=6,7
         write(iwu,'(/a,2g12.4,i4/)') 
     &   ' now running qsoll,beta,lqcon=',qsoll,betsol,lqcon
       enddo
      endif
c
c     the loop over temperatures in case of IPAIR=4.
c     take the loop parameters from TEMPMX,TEMPMN and TEMPW
c
      if(ipair.ne.4 .AND. suv_mue(2)==0D0) then    
       tempmn = zero
       tempmx = zero
       tempw  = one
      endif
      do 1210 tempac=tempmn,tempmx,tempw
      if(ipair.eq.4) then
        ggp = tempac
        ggn = tempac
        gap(1) = tempac
        gap(2) = tempac
        write(6,'(/a,g12.4/)') ' now running T=',tempac
        write(7,'(/a,g12.4/)') ' now running T=',tempac
      else
        suv_wanted(2) = tempac
        write(6,'(/a,g12.4/)') ' now running suv_wanted(2)=',tempac
        write(7,'(/a,g12.4/)') ' now running suv_wanted(2)=',tempac
      endif

*
***********************************************************************
*                                                                     *
*     the static iteration                                            *
*                                                                     *
***********************************************************************
*
      eold   = energy
      x0sav  = x0inv
      

      WRITE(*,*) 'test1: swicon,istblo:',  swicon,istblo
      WRITE(7,*) 'test1: swicon,istblo:',  swicon,istblo
      if(ABS(swicon-1D0)<1D-9 .AND. SUM(istblo)==0)then    ! for calc. of PES only fix mpsi
       do iq=1,2
        istblo(iq) = 0
        istcblock(iq) = 0
       enddo
      else
       call findconjpart(1)
       istcblock(1) = itpartn(istblo(1))
       call findconjpart(2)
       istcblock(2) = itpartn(istblo(2))
       write(6,*)'blocking(p/n):',istblo,istcblock  
       write(7,*)'blocking(p/n):',istblo,istcblock  
       write(17,*)'SKYAX LOOP: blocking(p/n)=',istblo(1),istblo(2)
       IF(ttrack) THEN
          WRITE(17,*) '4. kblock etc=',kblock,mblock,istblo
          WRITE(17,*) '4. kblock etc=',kblock,mblock,istblo
       END IF
      endif   


   
      WRITE(*,*) 'test2: swicon,istblo:',  swicon,istblo,istcblock
      WRITE(7,*) 'test2: swicon,istblo:',  swicon,istblo,istcblock
      IF(istblototal==0) THEN
        WRITE(*,*) 'before BLOCK'
        call block(blnprot,blnneut,outnam,x0sav)     ! main subroutine to calculate the optimal blocked state
        WRITE(*,*) 'returned from 2. BLOCK in SKYAX'
        IF(ttrack) THEN
          WRITE(7,*) 'returned from 2. BLOCK in SKYAX'
          WRITE(17,*) 'returned from 2. BLOCK in SKYAX'
        END IF
      ELSE
        call iteration_loop(outnam,x0sav,qact)
      END IF
c      call pairinfo 

      WRITE(*,*) 'one loop terminated qact,tempac,tempmn,tempmx,tempw='
     &         ,qact,tempac,tempmn,tempmx,tempw

c
 1210 continue
c
      WRITE(*,*) 'T loop terminated qact=',qact,qmin,qmax,qw
!      if (swicon.eq.1) swdiag=1          ! turn on pre-diagonalisation for
      WRITE(*,*) 'T loop terminated qact=',qact,iqact,iqdmax
                                         ! the next constraint iteration


c
! 1100 continue
      END DO
      if(Min(NX3,NY3,NZ3).gt.1) call writ3D()
      if(KRESTA) call writ2D()
*
***********************************************************************
*                                                                     *
*     end of loop over deformations                                   *
*                                                                     *
***********************************************************************
*
 917  continue
*
*     end of outer loop
*




 1009 continue
 1010 continue
      if(ireg(1).GT.0) call psiout(outnam)




  
c
cMPI      call MPI_Finalize(mpierror)
cMPIc
cMPI      if(mpirank.eq.0) then
cMPI         stop 'normal termination '
cMPI      else
cMPI         stop 
cMPI      endif
c
      stop 'normal termination '
c
c     error exits
c

  990 stop ' block structure of mother nad daughter is not the same'  
      
  991 stop ' error or EoF while reading TITLE'
  992 stop ' error or EoF while reading nucleus'
  993 stop ' error or EoF while reading force'
  997 stop ' error or EoF while reading grid'
  994 stop ' error or EoF while reading initialization'
  995 stop ' error or EoF while reading Q- and T-loops'
  996 stop ' error or EoF while reading print options'
  998 stop ' error or EoF while reading print iteration parameters'
  999 stop ' error or EoF while reading print CHF extras'
      end
c
c      include 'pairing.f'
c      include 'crankin.f'
c      include 'seprpa.f'
c      include 'qtest.f'
ctest      include 'xtime.for'


