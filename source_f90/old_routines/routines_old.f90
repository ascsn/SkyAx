!**********************************************************************
!                                                                     *
!                                                                     *
!   skyax  -   pgm to solve the static mean field equations           *
!              in cylindrical coordinates using Skyrme forces.        *
!              The numerical representation uses a full matrix        *
!              for the first and second derivatives. These            *
!              matrices are built from a Fourier definition.          *
!                                                                     *
!              This program has been developed from a                 *
!              relativistic mean-field code which, in turn,           *
!              had is roots in a Skyrme mean-field code using         *
!              finite differences.                                    *
!              maintainance of the program should carry a             *
!              version number, the history of the changes,            *
!              and the corresponding authors.                         *
!                                                                     *
!              History:                                               *
!               ~1982        --> axial Skyrm-Hf with 3-point F.D. .   *
!                                S.Umar, M.R. Strayer, S.J.Lee,       *
!                                P.-G.Reinhard                        *
!               ~1985        --> change to Relativistic MF with       *
!                                3-point F.D. and shifted grids.      *
!                                S.Umar, M.R. Strayer, S.J.Lee        *
!                                J. Maruhn, J. Fink, P.-G.Reinhard    *
!               ~1991        --> change of the Relativistic code      *
!                                to Fourier representations of        *
!                                the derivative matrices.             *
!                                J. Maruhn, J. Fink, P.-G.Reinhard    *
!               March 1993   --> change to Skyrme HF code using       *
!                                Fourier representations of           *
!                                the derivative matrices.             *
!                                P.-G.Reinhard                        *
!               1993-1993    --> testing and extension                *
!                                (e.g. LN pairing)                    *
!                                P.-G.Reinhard, M.Bender, J.Maruhn    *
!               Oct/Nov1996  --> z-asymmetry, delta-force BCS         *
!                                D.Dean, P.-G.Reinhard, M.R.Strayer   *
!                                                                     *
!**********************************************************************
!                                                                     *
!      Version     Date     Status                                    *
!       1.0.0    03/13/93   first setup of program.                   *
!                           Coulomb, spin-orbit, and pairing run.     *
!                           Agreement with the spherical code is      *
!                           tested to be better than 0.03% for        *
!                           doubly magic nuclei (i.e. spherical       *
!                           and without pairing).                     *
!       1.0.1    03/17/93   Input of initial maximum oscillator shell *
!                           is no done with 'deocc', relative to a    *
!                           first guess of the appropriate valence    *
!                           shell for given 'npart'. It has also been *
!                           added a check whether sufficient states   *
!                           are supplied in the 'ssetup'.             *
!       1.0.2    03/22/93   More forces in PARAM and flexible input   *
!                           of arbitrary parameters;                  *
!                           clean-up of plot files;                   *
!                           new form of printing s.p.properties;      *
!                           new pairing routine;                      *
!                           two different convergence criteria;       *
!                           fixing some minor bugs.                   *
!       1.0.3    05/28/93   damping of propagation of occupation      *
!                           numbers with 'addocc',                    *
!                           dynamical mixture of densities with       *
!                           'adddyn' replacing 'addnew'.              *
!                           minor changes in output formats.          *
!       1.0.4    09/15/93   direct evaluation of energy,              *
!                           distinction between 'addrho' and 'addpot'.*
!       1.0.5    09/16/94   new damping of constraint,                *
!                           Lipkin-Nogami (preliminary),              *
!                           various test options                      *
!       1.0.6    10/17/94   stabilizing wall (STAPOT),                *
!                           alternatve distance (DELRMN) for damping  *
!                           of constraint,                            *
!                           debugging such that variables can be      *
!                           supposed as 'dynamic' in subroutines      *
!       2.0.1.   11/16/96   Extension to asymmetry along z-axis,      *
!                           octupole constraint,                      *
!                           initial octupole "kick",                  *
!                           BCS-pairing with delta- or surface-force, *
!                           revised structure of input file.          *
!                                                                     *
!**********************************************************************
!                                                                     *
!      Here follows a list of topics for the near future development: *
!                     --------------         -----------              *
!         -   Deformation properties and pairing should still be      *
!             tested, although both seems superficially o.k.          *
!             In particular, watch the convergence of the various     *
!             deformations 'qtil', 'qold', and 'q' towards each       *
!             other; there may still be a problem.                    *
!             This problem exists, in particular, if the previous     *
!             solution is reused in a deformation or temperature      *
!             loop.                                                   *
!         -   The iteration of the quadrupole constraint may also     *
!             be speeded up if that turns out to be necessary.        *
!         -   Change input of desired (constrained) deformation to    *
!             the dimensionless quadrupole momentum 'beta'.           *
!         -   Improve the treatment of the boundary 'r=0' for the     *
!             spin-orbit current.                                     *
!         -   The handling of the c.m. zero-point energy may be       *
!             improved by computing <P^2> microscopically. A similar  *
!             zero-point energy for rotation may be added.            *
!                                                                     *
!         !   yet to code: s.c. resp[onse for Lipkin Nogami           *
!                          use of 'vsqobs' for spatial observables    *
!                                                                     *
!                                                                     *
!**********************************************************************

!  The input parameters can be grouped in two categories,

!  I. parameters which must be specified by the user:
!  --------------------------------------------------

!     title      - a character string carrying the title of the run
!                  the first three digits are reserved for naming the
!                  output files. the name 'plt' is given if these are blank.

!     nprot      - proton number
!     nneut      - neutron number

!     t0..t3,x0..x3 ->  Skyrme parameters (only active for iparam=0)
!                       iparam=-1 switches to Woods-Saxon whose
!                       parameters are read in right after the 'iparam'
!     b4         - =t4/2 spin-orbit parameter
!     b4p        - new parameter for isotopic spin-orbit force
!                  old convention is recovered for b4p=0.0
!                  vanishing b4p is to be entered as a small b4p>1.E-20
!     alpha      - power of density dependence
!     mass       - hbar^2/2m

!     ifex       - switch for coulomb exchange in slater approximation
!     ifls       - switch for contribution from t1,t2 to spin-orbit force
!     ifcm       - switch for c.m. correction
!                   0 --> m_eff = m-m/A
!                   1 --> E_cm=<P^2>/2M approx. as 17.2MeV/A^0.2
!     chg        -  effective charge for direct field
!                   (use 1.0 for full Coulomb,
!                    0.0 for No Coulomb)

!     ipair      - type of pairing:
!                  -1 = filling approximation
!                   0 = constant Delta approximation
!                   1 = constant Force approximation (BCS)
!                   2 = constant Force approximation with
!                       s.p. correction (BCS)
!                   3= Lipkin-Nogami
!                   4 = thermal HF, temperature handled by "DELTP"
!                   5 = BCS with delta force
!                   7 = BCS with dens.dep. delta-force
!                   9 = BCS-LN with delta force
!                  11 = BCS-LN with dens.dep. delta-force
!                 note: ipair=5,7,9 and 11 require available nr. of states
!                       high enough to cover 2*nr.of.particles
!     ggp,ggn    - gap or strength parameters for the BCS pairing
!                   0 --> defaults are used for ipair=0,1,2,3
!     pairstab   - stabilizing parameter for SBCS (0.0=BCS, has to be <0)

!     sym        - sym=1.0 for calculations which are
!                  reflection symmetric about z=0.0
!                  otherwise sym=0.0
!     nz         - no. of z direction mesh points used in the
!                  calculation. note nz < nzd.
!                  'nz=0' calls for default setting according to 'nmass'
!     nr         - no. of r direction mesh points used in the
!                  calculation. note nr < nrd.
!                  'nr=0' calls for default setting according to 'nmass'
!     dz         - mesh spacing in the z direction.
!     dr         - mesh spacing in the r direction.
!     icutr      - momentum comp. cutoff in r-direction, ! set 0 !
!     icutz      - momentum comp. cutoff in z-direction  ! set 0 !
!     idrtyp     - switch for type of derivative
!                  0 -> fourier definition in z&r direction ;
!                  1 -> fourier definition z direction, but
!                       finite diff. of 11. order in r direction;
!                  2 -> fourier definition r direction, but
!                       finite diff. of 11. order in z direction;
!                  3 -> finite diff. 11. order in r&z direction.

!     betain     - deformaion parameter for the initial axial
!                  harmonic oscillator wavefunctions,
!                  for 'betain=0' a default 'betain' is deduced
!                  from the actual deformation
!     octini     - initial kick towards octupole deformation
!                  deforms with 'octini'*Q_30/<r>^3;
!                  use small values to avoid curious nods in the w.f.
!                  for 'octini=0' a default 'octini' is deduced
!                  from the actual deformation
!     deocc(1/2) - addition to the first guess of 'eocc' with
!                  'eocc = 1.4*a**1/3 + deocc'. where 'eocc(1/2)'
!                  is the cutoff for the number of initial
!                  oscillator shells. 'deocc' can be near 0 for
!                  magic nucleon number. somewhat larger values
!                  are recommended far away from magic shells.
!                  A value of 0.8-1.2 is recommended for ipair=5,7,9 and 11.
!     mmax(1/2)  - maximum value for angular momentum 'm' in the
!                  setup phase. 'mm(1)' stands for protons and
!                  'mm(2)' for neutrons.
!     dipin      - strength of external isovector dipole potential

!     swicon      - switch for the constraint quadrupole
!                   1 --> constraint,  0 --> no constraint
!     betmin      - minimum deformation
!                   the given deformations are also used for initialization
!                   if 'betain=0' or 'octini=0' indicate this default choice,
!                   note that 'betmin' allows to transfer initialization
!                   even in case of unconstrained iteration.
!     betmax      - maximum deformation
!     betaw       - step in deformation
!                   these entries use dimensionless deformations
!                   of quadrupole or octupole type depending on 'lqcon'
!     lqcon       - type of constraint
!                    2 --> quadrupole
!                    3 --> octupole (only for ISYM=0)
!                    4 --> hexadacopole
!     ifresh      - switch for refresh of inital w.f. at each new deformation
!                    0 --> continue next deformation from old w.f.
!                    1 --> start next deformation with new w.f. from setup

!     tempmn      - minimum temperatur in T-loop for ipair=4
!     tempmn      - maximum temperatur in T-loop for ipair=4
!     tempw       - step in temperatur in T-loop for ipair=4

!  II. parameters which should only be changed by the experienced user:
!  --------------------------------------------------------------------

!     iprint     - modulus for printing output
!     itmax      - maximum no. of static iterations for one given
!                  quadrupole momentum
!     iplot      - modulus for printer plot of the density contours
!     itest      - test level = level test prints
!     itpfix     - iteration after which pairing is fixed
!     itmxt      - maximum number of all iterations (over all
!                  deformations)

!     e0inv      - energy cutoff for damping operator
!     x0inv      - expansion parameter for dpa step
!     addrho     - weight for adding new density
!     addpot     - weight for adding new potential
!                    'addrho' and 'addpot' have about the same effect,
!                    slowing down the propagation of the mean field.
!                    'addpot' is to be preferred if one aims top trace
!                    the energy convergence using 'subroutine energ' as
!                    this gives always the actual energy.
!     addocc     - weight for adding new occupation numbers from
!                  the pairing iteration.

!     errenr     - first convergence criterion:
!                  change in energy must come below 'errenr'
!     errvar     - second convergence criterion:
!                  the variance must come below 'errvar'

!     ivs        - level are reordered (and reoccupied) according
!                  to their s.p. energies as long as
!                      iter < ivs  .
!                  no reordering is done for
!                      iter > ivs.
!     itsub      - number of sub-iterations without orthogonalization  !mod 5/93
!                  presently active are:
!                   1 --> standard iteration with one 'hpsi' per iteration
!                   2 --> two damped steps per iteration, second step
!                         is weighted with 'x0fac'
!                  -2 --> second step is undamped and has stepsize x0fac
!     x0fac      - factor on 'x0inv' at even subiterations             !mod 5/93

!     errf        - allowed error in the fermi energies
!     erra        - allowed error in the particle no.
!     fermip      - upper limit for the proton fermi
!                   energy
!     fermin      - upper limit for the neutron fermi
!                   energy
!     npiter      - limit on the no. of iterations for
!                   self-consistent solutions to the
!                   gap equations
!     itpini      - nr. of iterations BEFORE pairing is switched on

!     addocc      -  mixing of new and old occ. numbers

!     ndrun      - switch for dipole field in coulomb solver
!     nqrun      - switch for quadrupole field in coulomb solver
!     norun      - switch for octupole field in coulomb solver
!     nhrun      - switch for hexadecoupole field in coulomb solver

!     colerr     - level of precision for coulomb solver
!     irepc      - maximum number of iterations for coulomb solver
!     cmssqr     - damping mass (inverse length) for coulomb solver

!     yukerr     - level of precision for yukawa folding
!     snl        - damping mass (inverse length) for yukawa folding
!     irepm      - maximum number of iterations for yukawa folding

!     qerr        - accuracy for quadrupol deformation
!     c0,d0       - convergence-parameters for constraint
!     qeps        - effective stepsize in quadrupole iteration
!     itcons      - maximum number of constraint steps within iteration.
!                   the iteration returns to the standard step with
!                   the now fixed constraining force 'qlambd*qfac(iz,ir)'.
!                   a consistent quadrupole-constrained iteration is
!                   only obtained for 'itcons>itmax'.
!                   The 'itcons>0' serves also as a switch to iteration with
!                   fixed constraining force 'qlam': for 'itcons=0', the
!                   'qmin,qmax,qw' are deduced from 'betmin,betmax,betaw'
!                   and are used as the values for 'qlam'.
!                   The 'itcons<0' is used to switch to unconstrained
!                   iteration after |itcons| inital constrained steps.
!                   This serves to drice pure HF iteration near a
!                   desired deformation minimum.
!                   Thus use a very large 'itcons' to be on the safe side!
!                   ?? note: this feature may be removed some time later.
!     idampc      - option to switch the damping of the cnstraint
!                   0 --> old form with Woods-Saxon hat
!                   1 --> new dynamical damping in 'damp_constraint'

!     ifstab      - switch for stabilizing external potentianl,
!                   recommended for high temeperatures and/or
!                   exotic nuclei
!     stwid       - parameter for the energy-width of the states included
!                   in the generating density
!     stslop      - slope of the stabilizing potential
!     rhocut      - 'rhocut*rho0' is the density from which on the
!                   asymptotic evaluation of the additional potential
!                   starts to count

!     dalp        - radius of the woods-saxon which damps the quadrupole
!                   operator at large distances
!     dgam        - width  of the woods-saxon which damps the quadrupole
!                   operator at large distances
!     iflin       - switch to linear asymptotic constraint
!     swifac      - factor on switch in z-direction,
!                   switch from quadratic to linear shape comes at
!                     swifac*1.2*A**1/3

!     iparam - number of the parametriztion used
!              0 -> direct input 't0,...,power' from for005
!              1,...,8 -> select predefined forces in subroutine 'param'
!              -1 -> deformed WS

!     ireg(1..2)  - switch to printing of levels in Nilsson scheme
!       kblock   =  2*K values of blocked s.p. state                  *
!       mblock   =  M values of blocked s.p. state                    *
!       isoblo   =  isospin of blocked s.p. state                     *
!       istblo   =  number of blocked s.p. state in level sequence    *
!     locc        - explicit preoccupation
!                   --> presently disabled, use -1,,,,

SUBROUTINE testls

!     test of l*s current by checking vector field against its divergence

!USE paramr
USE grid
USE wavefr, ONLY: wrk3
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'



!-----------------------------------------------------------

CALL div(wrk3,rholsr(0,0,1),rholsz(0,0,1))
WRITE(6,'(/a/a)') ' test spin-orbit fields for protons:',  &
    '  ir   iz   div-direct    div-deduced '
DO ir=1,nr
  DO iz=0,nz
    WRITE(6,'(2i5,3(1pg13.4))')  &
        ir,iz,drhols(iz,ir,1),wrk3(iz,ir),drhols(iz,ir,1)-wrk3(iz,ir)
  END DO
END DO

CALL div(wrk3,rholsr(0,0,2),rholsz(0,0,2))
WRITE(6,'(/a/a)') ' test spin-orbit fields for neutrons:',  &
    '  ir   iz   div-direct    div-deduced '
DO ir=1,nr
  DO iz=0,nz
    WRITE(6,'(2i5,3(1pg13.4))')  &
        ir,iz,drhols(iz,ir,2),wrk3(iz,ir),drhols(iz,ir,2)-wrk3(iz,ir)
  END DO
END DO

RETURN
END SUBROUTINE testls

!-----pripot-----------------------------------------------------------

SUBROUTINE pripot(title,outnam)

!     print various potentials on file 'potentials.<name>'

!USE paramr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'
!INCLUDE 'wavefr.inc'


CHARACTER (LEN=80), INTENT(OUT)          :: title
CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam



!---------------------------------------------------------------------------

OPEN(31,STATUS='unknown',FORM='formatted', FILE='potentials.'//outnam)
DO ispool=1,69999
  READ(31,'(1x)',END=9999)
END DO
STOP " file count exhausted for print of potentials "
9999 CONTINUE

WRITE(31,'(a/2a/a,2i5/a,2f6.1,2i5/a,3f7.3)')  &
    ' potentials:',' title=',title, ' A,Z=',npart(1)+npart(2),npart(1),  &
    ' dz,dr,nz,nr=',dz,dr,nz,nr, ' beta2,beta2p,beta2n=',beta2,beta2p,beta2n
WRITE(31,'(2a)') '  iz ir    z     r     V_proton      V_neutron   ',  &
    '    V_Coulomb   V_prot-V_coul '
DO ir=1,nr
  DO iz=0,nz
    r  = (ir+sr)*dr
    z  = (iz+sz-icm)*dz
    WRITE(31,'(1x,2i3,2f6.1,4g15.7)')  &
        iz,ir,z,r,pot(iz,ir,1),pot(iz,ir,2),coulpt(iz,ir),  &
        pot(iz,ir,1)-coulpt(iz,ir)
  END DO
END DO

CLOSE(31)
WRITE(6,'(a)') ' potentials written on file.'

RETURN
END SUBROUTINE pripot
!testc-----yuktest-----------------------------------------------------------
!testc
!test      subroutine yuktest
!testc
!testc**********************************************************************
!testc                                                                     *
!testc     tests the yukawa solver                                         *
!testc                                                                     *
!testc**********************************************************************
!testc
!test      include 'paramr.inc'
!test      include 'fields.inc'
!test      include 'dynpar.inc'
!testc
!test      dimension rhs(0:nzd,0:nrd)
!testc
!testc
!test      a      = 4.0
!test      b      = sqrt(h2meff*a**(one/3)/41.0)
!testc
!test      do 5 ir=0,nr
!test      rsqr          = ((ir)+sr)**2 * drsqr
!test      do 5 iz=0,nz
!test      zsqr          = ((iz)+sz)**2 * dzsqr
!test      rhosqr        = rsqr + zsqr
!testcc    exarg1        = (rho-rad)/sur
!testcc    exarg2        =-(rho+rad)/sur
!testc
!testcc    if (earg1.ge.-eighty) then
!testcc      exfac1      = exp(exarg1)
!testcc    else
!testcc      exfac1      = 0.0
!testcc    endif
!testcc    if (exarg2.ge.-eighty) then
!testcc      exfac2      = exp(exarg2)
!testcc    else
!testcc      exfac1      = 0.0
!testcc    endif
!testc
!testcc    opex1         = (one+exfac1)
!testcc    opex2         = (one+exfac2)
!testcc    fac1          =-exfac1/(sur*opex1**2*opex2)
!testcc    fac2          = exfac2/(sur*opex2**2*opex1)
!testcc    fac3          = 2.0*exfac1**2/(sur*opex1**3*opex2)
!testcc    fac4          =-2.0*erd/(sur*opex1**2*opex2**2)
!testcc    fac5          = 2.0*exfac2**2/(sur*opex1*opex2**3)
!testc
!testcc    td2           = (fac1-fac2+fac3+fac4+fac5)/sur
!testcc    d1            =  fac1+fac2
!testcc    td1           = 2.0*d1/rho
!testcc    f             = one/(opex1*opex2)
!testcc    tmass         = xmass**2*f
!testc
!testcc    rhs(iz,ir)    = (-td2-td1+tmass)*v0
!testcc    smeson(iz,ir) = f*v0
!test      pfac    =  4.0*vcoupl/(pi**1.5*b**3)
!test      earg    =  (-rhosqr/b**2)
!test      if (abs(earg).le.50.0*one) then
!test         rhs(iz,ir) = pfac*exp(earg)
!test      else
!test         rhs(iz,ir) = 0.0
!test      endif
!test 5    continue
!test      call mesini(vmeson,rhs,vmshc)
!test      call mesrun(vmeson,rhs,vmshc,2,1)
!testc     call mesini(vmeson,rhs,vaccel,vmass)
!testc     call mesrun(vmeson,rhs,vaccel,vmass)
!testc
!testcc    errmx = -one
!testcc    do 6 ir=0,nr
!testcc    do 6 iz=0,nz
!testcc     err= abs(vmeson(iz,ir)-smeson(iz,ir))
!testcc     if (err.gt.errmx) then
!testcc       errmx = err
!testcc       irmx  = ir
!testcc       izmx  = iz
!testcc     endif
!testc6    continue
!testc     tsum = 0.0
!testc     do 6 ir=0,nr
!testc     r  = (sr+ir)
!testc     do 6 iz=0,nz
!testc       tsum = tsum+vmeson(iz,ir)*r
!testc6    continue
!testc     tsum = tsum*dz*dr
!testc     write(7,'(a,g15.5)') ' testsum=',tsum
!test      call test_print(rhs,0,9,0,9,' density')
!testcc    write(7,'(/a,g15.5,2(a,i3)/)') ' error(ana-num)=',errmx,' ir=',
!testcc   &                             irmx,' iz=',izmx
!testc
!testc      call test_print(vmeson,0,9,0,9,' vmeson')
!testc     call test_print(smeson,0,9,0,9,' ana result')
!testc
!testc     call test_print(vmeson,17,25,26,34,' num result')
!testc     call test_print(smeson,17,25,26,34,' ana result')
!testc
!testc     call test_print(vmeson,0,9,11,20,' num result')
!testc     call test_print(smeson,0,9,11,20,' ana result')
!testc
!testc     call test_print(vmeson,11,20,11,20,' num result')
!testc     call test_print(smeson,11,20,11,20,' ana result')
!testc
!test      stop
!test      end

SUBROUTINE prifou(field,text)

!     print fourier components of 'field' along z- and -r-axis

!USE paramr
!IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'dynpar.inc'
!INCLUDE 'wavefr.inc' ! wrk1


REAL(8), INTENT(IN)                  :: field(0:nzd,0:nrd)
CHARACTER (LEN=*), INTENT(IN)        :: text


ALLOCATE(wrk1(0:nz,0:nr))

!      data zero,one /0.0,1.0/
!      data pi /3.141592653589793/

!     r-direction

dn    = (nr)
afac  = pi/dn
!      pfac  = 2.0/dn
DO kr=0,nr
  cfac = kr*afac
  sum  = 0.0
  DO ir=0,nr
    cosfac = COS(cfac*ir)
    IF(ir == 0 .OR. ir == nr) THEN
      cosfac = 0.5*cosfac
    END IF
    sum  = 0.5*cosfac*field(0,ir)+sum
    DO iz=1,nz-1
      sum  = cosfac*field(iz,ir)+sum
    END DO
    sum  = 0.5*cosfac*field(nz,ir)+sum
  END DO
  wrk1(0,kr) = sum
END DO
WRITE(7,'(3a,10(/1x,5g12.4))')  &
    'cosine transform of ',text,' in r-direction:', (wrk1(0,ir),ir=0,nr)

!     z-direction

dn    = (nz)
afac  = pi/dn
!      pfac  = 2.0/dn
DO kz=0,nz
  cfac = kz*afac
  sum  = 0.0
  DO iz=0,nz
    cosfac = COS(cfac*iz)
    IF(iz == 0 .OR. iz == nz) THEN
      cosfac = 0.5*cosfac
    END IF
    sum  = 0.5*cosfac*field(iz,0)+sum
    DO ir=1,nr-1
      sum  = cosfac*field(iz,ir)+sum
    END DO
    sum  = 0.5*cosfac*field(iz,nr)+sum
  END DO
  wrk1(kz,0) = sum
END DO
WRITE(7,'(3a,10(/1x,5g12.4))')  &
    'cosine transform of ',text,' in z-direction:', (wrk1(iz,0),iz=0,nz)
DEALLOCATE(wrk1)

RETURN
END SUBROUTINE prifou

!tstc-------coltst---------------------------------------------------------
!tstc
!tst      subroutine coltst(rhoin,wcol,idef)
!tstc
!tstc     test routine for coulomb solver
!tstc
!tst      include 'paramr.inc'
!tst      include 'levels.inc'
!tst      include 'fields.inc'
!tst      include 'dynpar.inc'
!tstc
!tst      dimension  wcol(0:nzd,0:nrd),rhoin(0:nzd,0:nrd)
!tstc
!tstc     pfacd = (3/(2*pi))**(1.5), rp: ca. rms radius of protons
!tstc
!tstcold      data pfacd,opf,tpf /0.32992261,1.5,2.5/
!tstc
!tstc
!tstc     monopole field and density
!tstc
!tst      qdef    = 40.0
!tst      pro     = (npart(1))
!tst      zmax    = ((nz-icm)+sz)*dz
!tst      rmax    = ((nr)+sr)*dr
!tst      zrmin   = min(zmax,rmax)
!tst      zrminh  = zrmin/6.5
!tstc
!tst      dis     = zrmin*0.06
!tst      xnq     = qdef/(two*two*dis*dis)
!tst      xaq     = sqrt(-log(oh))/zrminh
!tst      xaq2    = xaq  * xaq
!tst      xaq3    = xaq2 * xaq
!tst      d0q     = xaq3*xnq/pi**1.5
!tstcold      pfacq   = two*pi/xaq2
!tstc
!tst      am      = sqrt(-log(oh))/zrminh
!tst      am2     = am*am
!tst      am3     = am2*am
!tst      d0m     = am3*pro/pi**1.5
!tst      pfac2   = two*pi*d0m/am2
!tst      do 20 ir=0,nr
!tst      r = (ir+sr)*dr
!tst      r2= r*r
!tst      do 20 iz=0,nz
!tst         z     = (iz+sz-icm)*dz
!tst         z2    = z*z
!tst         rad2  = z2+r2
!tst         rad   = sqrt(rad2)
!tstc
!tst         expa  = rad2*am2
!tst         if (expa.le.eighty) then
!tst            expf = exp(-expa)
!tst         else
!tst            expf = zero
!tst         endif
!tst         dichm        = d0m*expf
!tst         radam        = rad*am
!tstc
!tstc        the function subroutine errf2 provides the integral
!tstc          errf2(x) = 2/gamma(3/2) * int_{0}^{x} exp(-t*t) t**2 dt
!tstc
!tst         if (rad.lt.1.d-2) then
!tst          qt1 = zero
!tst         else
!tst          qt1 = errf2(radam)
!tst          qt1 = qt1*pro/rad
!tst         endif
!tst         qt2 = pfac2 * expf
!tst         feldm        = chg*(qt1+qt2)
!tstc
!tstc        quadrupol:
!tstc
!tst       zpd    = z + dis
!tst       zpd2   = zpd*zpd
!tst       zmd    = z - dis
!tst       zmd2   = zmd*zmd
!tst       radpd2 = r2+zpd2
!tst       radmd2 = r2+zmd2
!tstc       radpd  = sqrt(radpd2)
!tstc       radmd  = sqrt(radmd2)
!tst       expaq  = xaq2*rad2
!tst       expapd = xaq2*radpd2
!tst       expamd = xaq2*radmd2
!tst       if (expaq.le.eighty) then
!tst         expfq  = exp(-expaq)
!tst       else
!tst         expfq  = zero
!tst       endif
!tst       if (expapd.le.eighty) then
!tst         expfpd = exp(-expapd)
!tst       else
!tst         expfpd = zero
!tst       endif
!tst       if (expamd.le.eighty) then
!tst         expfmd = exp(-expamd)
!tst       else
!tst         expfmd = zero
!tst       endif
!tst       dens   = d0q*(expfmd-2.0*expfq+expfpd)
!tstc
!tst       dichq        = dens
!tst       if (idef.eq.0) then
!tst        rhoin(iz,ir) = dichm
!tst        wcol(iz,ir)= feldm
!tst       else if (idef.eq.1) then
!tst        rhoin(iz,ir) = dichm+dichq
!tst        wcol(iz,ir)= feldm+dichq
!tst       endif
!tst20    continue
!tstc     do 10 ir=0,nr
!tstc     do 10 iz=0,nz
!tstc       wcol(iz,ir) = coulm(iz,ir)
!tstc10   continue
!tstc
!tstc     call test_print(coulm,0,9,0,9,' coulm')
!tstc     call test_print(crhom,0,9,0,9,' crhom')
!tstc
!tst      return
!tst      end
!tstc-----errf2 -----------------------------------------------------------
!tstc
!tst      function errf2(x)
!tstc
!tstc        the function subroutine errf2 provides the integral
!tstc          errf2(x) = 2/gamma(3/2) * int_{0}^{x} exp(-t*t) t**2 dt;
!tstc        uses reduction to error function and tschebyscheff approximation.
!tstc        (for the latter see 'numerical recipes', chapter 6.2).
!tstc
!tstc
!tst      include 'paramr.inc'
!tstc
!tstc     compute error function for x with fractional error less 1.2e-7.
!tstc
!tst      z     = abs(x)
!tst      t     = 2.0/(2.0+z)
!tst      erfc  = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+
!tst     &       t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+
!tst     &       t*(1.48851587+t*(-0.82215223+t*0.17087277)))))))))
!tst      if(x.lt.0) erfc = 2.0-erfc
!tst      erf   = 1.0-erfc
!tstc
!tstc     now compose the final function
!tstc
!tst      t     = 1.772453851*exp(x*x)
!tst      errf2 = -(x+x-t*erf)/t
!tstc
!tst      return
!tst      end



!=======================================================================

!     ini_coulomb(mangmol,mangmom)

!     initialize the coulomb solver (background fields and densities,
!     damping mass.

!       mangmol :  maximum angular momentum l of background fields
!       mangmom :  maximum z-component m of angular momentum

!=======================================================================

SUBROUTINE ini_coulomb(mangmol,mangmom)

USE paramr
USE constr
USE fields
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)


INTEGER, INTENT(IN OUT)                  :: mangmol
INTEGER, INTENT(IN OUT)                  :: mangmom


INTEGER :: ir, iz

REAL(8),PARAMETER :: ddw=0.4D0          ! width of damping
REAL(8),PARAMETER :: ddmass=1D0       ! factor on damping mass
REAL(8),PARAMETER :: d2=2D0           ! thickness of surface zone
REAL(8),SAVE :: dqr=0D0          ! to make this variable static

!     ...
!-----------------------------------------------------------------------
zmax    = ((nz-icm)+sz)*dz
rmax    = ((nr)+sr)*dr
drmin = MIN(zmax,rmax) - d2
drmin2 = drmin * drmin
dzadd = zmax - drmin - d2
dradd = rmax - drmin - d2
dtermr = drmin2 - dradd*dradd
dtermz = drmin2 - dzadd*dzadd

!     initialize damping mass
!-----------------------------------------------------------------------
DO iz=0,nz
  z  = (iz+sz-icm)*dz
  z2 = z * z
  
  IF( dzadd >= dradd ) THEN
    IF (z >= dzadd) THEN
      dqr = SQRT(dtermz + d2*dzadd*z)
    ELSE
      dqr = SQRT(z2 + drmin2)
    END IF
  END IF
  
  DO ir=0,nr
    r  = (ir+sr)*dr
    r2 = r * r
    drr = DSQRT(z2 + r2)
    IF( dradd > dzadd ) THEN
      IF (r >= dradd) THEN
        dqr = SQRT(dtermr + d2*dradd*r)
      ELSE
        dqr = SQRT(r2 + drmin2)
      END IF
    END IF
    sqrm(iz,ir) = ddmass / (one + EXP((dqr-drr)/ddw))
  END DO
END DO

IF(itest > 0) CALL test_print(sqrm,0,nz,0,nr,'sqrm')

RETURN
END SUBROUTINE ini_coulomb


!======================================================================

!     wrt_spprop()

!     print sinlge-particle-properties

!=======================================================================


SUBROUTINE wrt_spprop(outnam,isinit)

!----------------------------------------------------------------------

USE paramr
USE dynpar
USE observ
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'

!----------------------------------------------------------------------


CHARACTER (LEN=3), INTENT(IN)            :: outnam
INTEGER, INTENT(IN)                      :: isinit
INTEGER :: ist

CHARACTER (LEN=5) :: outnum


!----------------------------------------------------------------------

!     if initialization should be done (isinit=1)
!----------------------------------------------------------------------
IF (isinit == 1) THEN
  DO ist=1,nst(1)+nst(2)
    
!     generate a string containing the index of the file
!----------------------------------------------------------------------
    IF (ist <= nst(1)) THEN
      WRITE(UNIT=outnum, FMT='(f5.3)') 1.0D0+ist/1000D0
    ELSE
      WRITE(UNIT=outnum, FMT='(f5.3)') 2.0D0+(ist-nst(1))/1000.d0
    END IF
    
!     open file and write header
!----------------------------------------------------------------------
    OPEN(77,FILE='levels.'//outnum//'.'//outnam)
    WRITE(77,'(//a)') '# beta2, epsilon_i, Var_i, w_i, Delta_i, wtpair_i'
    CLOSE(77)
    
  END DO
  
ELSE
  
!     standard call (isinit=0)
!----------------------------------------------------------------------
  DO ist=1,nst(1)+nst(2)
    
!     generate a string containing the index of the file
!----------------------------------------------------------------------
    IF (ist <= nst(1)) THEN
      WRITE(UNIT=outnum, FMT='(f5.3)') 1.0D0+ist/1000D0
    ELSE
      WRITE(UNIT=outnum, FMT='(f5.3)') 2.0D0+(ist-nst(1))/1000.d0
    END IF
    
!     open file and write header
!----------------------------------------------------------------------
    OPEN(77,FILE='levels.'//outnum//'.'//outnam, POSITION='append')
    WRITE(77,'(f7.4,1(1pg14.7),1(1pg14.5),3(1pg14.5))') beta2,  &
        spenrg(ist),spevar(ist), vsqf(ist),deltaf(ist),wtpair(ist)
    CLOSE(77)
    
  END DO
END IF

!----------------------------------------------------------------------

RETURN
END SUBROUTINE wrt_spprop

!=======================================================================

!     restore_boundary(psia,imu,izpar,ist)

!     restore boundary values of wavefunction

!=======================================================================

SUBROUTINE restore_boundary(psia,imu,izpar,ist)

!------------------------------------------------------------------------

USE paramr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'dynpar.inc'


REAL(8), INTENT(OUT)                        :: psia(0:nzd,0:nrd,nsd)
INTEGER, INTENT(IN)                      :: imu
INTEGER, INTENT(IN)                      :: izpar
INTEGER, INTENT(IN OUT)                  :: ist



INTEGER :: ir,iz

!------------------------------------------------------------------------

RETURN                           !!!!!!!!!!!!!!!!!!!!!!!!!

IF (mpsi(ist)+imu == 0) THEN
  DO iz=0,nz
    psia(iz, 0,2)=0D0
    psia(iz,nr,1)=0D0
    psia(iz,nr,2)=0D0
  END DO
ELSE IF (mpsi(ist)+imu == -1) THEN
  DO iz=0,nz
    psia(iz, 0,1)=0D0
    psia(iz,nr,1)=0D0
    psia(iz,nr,2)=0D0
  END DO
ELSE
  DO iz=0,nz
    psia(iz, 0,1)=0D0
    psia(iz,nr,1)=0D0
    psia(iz,nr,2)=0D0
  END DO
END IF

IF (isym == 0) THEN
  DO ir=0,nr
    psia(0, ir,1)=0D0
    psia(nz,ir,1)=0D0
    psia(0, ir,2)=0D0
    psia(nz,ir,2)=0D0
  END DO
ELSE
  IF (npar(ist)*izpar == -1) THEN
    DO ir=0,nr
      psia(0 ,ir,1)=0D0
      psia(nz,ir,1)=0D0
      psia(nz,ir,2)=0D0
    END DO
  ELSE
    DO ir=0,nr
      psia(0, ir,2)=0D0
      psia(nz,ir,1)=0D0
      psia(nz,ir,2)=0D0
    END DO
  END IF
END IF

!------------------------------------------------------------------------

RETURN
END SUBROUTINE restore_boundary

!=======================================================================

!     upd_static_MF(psia,drpsia,dzpsia,dppsia,ddpsia,zpar,imu,bpar)

!     calculates the static mean-field densities and
!     mean-field potentials

!=======================================================================

SUBROUTINE upd_static_mf(phi,drphi,dzphi,dpphi,ddphi,occphi, mfp_u,  &
    mfp_b,mfp_drb,mfp_dzb, mfp_wz,mfp_wr)

USE paramr
USE levels
USE grid
USE wavefr ! ???
USE dynpar !??
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'dynpar.inc'
!INCLUDE 'wavefr.inc'   


REAL(8), INTENT(IN OUT)         :: phi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN)             :: occphi(npsid)
REAL(8), INTENT(OUT)            :: mfp_u(0:nzd,0:nrd,2)
REAL(8), INTENT(OUT)            :: mfp_b(0:nzd,0:nrd,2)
REAL(8), INTENT(OUT)            :: mfp_drb(0:nzd,0:nrd,2)
REAL(8), INTENT(OUT)            :: mfp_dzb(0:nzd,0:nrd,2)
REAL(8), INTENT(OUT)            :: mfp_wz(0:nzd,0:nrd,2)
REAL(8), INTENT(OUT)            :: mfp_wr(0:nzd,0:nrd,2)


INTEGER :: istart(2),iend(2)
INTEGER :: iz,ir,iiso
REAL(8) :: occ
REAL(8) :: b0,b0p,b1,b1p,b2,b2p,b3,b3p,b4,b5,b5p
REAL(8) :: wgt_b0 ,wgt_b1 ,wgt_b2 ,wgt_b3 ,wgt_b4 ,wgt_b5
REAL(8) :: wgt_b0c,wgt_b1c,wgt_b2c,wgt_b3c,wgt_b4c,wgt_b5p
REAL(8) :: fac1,fac2,fac3
REAL(8) :: buff1,buff2,buff3

REAL(8) :: mfd_rho(0:nzd,0:nrd,2)
REAL(8) :: mfd_tau(0:nzd,0:nrd,2)
REAL(8) :: mfd_jz(0:nzd,0:nrd,2)
REAL(8) :: mfd_jr(0:nzd,0:nrd,2)
REAL(8) :: mfd_grhor(0:nzd,0:nrd,2)
REAL(8) :: mfd_grhoz(0:nzd,0:nrd,2)
REAL(8) :: mfd_drho(0:nzd,0:nrd,2)
REAL(8) :: mfd_divj(0:nzd,0:nrd,2)
REAL(8) :: mfp_unl(0:nzd,0:nrd,2)

!     ------------------------------------------------------------------

istart(1)=1
iend(1)  =nst(1)
istart(2)=nst(1)+1
iend(2)  =nst(1)+nst(2)

DO iiso=1,2
  
!        set the mean-field densities to zero
!        ---------------------------------------------------------------
  DO ir=0,nr
    DO iz=0,nz
      mfd_rho(iz,ir,iiso)   = 0D0
      mfd_tau(iz,ir,iiso)   = 0D0
      mfd_jr(iz,ir,iiso)    = 0D0
      mfd_jz(iz,ir,iiso)    = 0D0
    END DO
  END DO
  
!        accumulate the mean-field densities
!        ------------------------------------------------------------
  DO ist=istart(iiso),iend(iiso)
    occ=2.0D0*occphi(ist)
    DO ir=0,nr
      DO iz=0,nz
        mfd_rho(iz,ir,iiso) = mfd_rho(iz,ir,iiso) + occ*(  &
            phi(iz,ir,1,ist)*phi(iz,ir,1,ist)  &
            + phi(iz,ir,2,ist)*phi(iz,ir,2,ist) )
        
        mfd_tau(iz,ir,iiso) = mfd_tau(iz,ir,iiso) + occ*(  &
            drphi(iz,ir,1,ist)*drphi(iz,ir,1,ist)  &
            + drphi(iz,ir,2,ist)*drphi(iz,ir,2,ist)  &
            + dzphi(iz,ir,1,ist)*dzphi(iz,ir,1,ist)  &
            + dzphi(iz,ir,2,ist)*dzphi(iz,ir,2,ist)  &
            + dpphi(iz,ir,1,ist)*dpphi(iz,ir,1,ist)  &
            + dpphi(iz,ir,2,ist)*dpphi(iz,ir,2,ist) )
        
        mfd_jr(iz,ir,iiso)  = mfd_jr(iz,ir,iiso)  + occ*(  &
            + dpphi(iz,ir,1,ist)*phi(iz,ir,1,ist)  &
            - dpphi(iz,ir,2,ist)*phi(iz,ir,2,ist)  &
            + dzphi(iz,ir,2,ist)*phi(iz,ir,1,ist)  &
            - dzphi(iz,ir,1,ist)*phi(iz,ir,2,ist) )
        
        mfd_jz(iz,ir,iiso)  = mfd_jz(iz,ir,iiso)  + occ*(  &
            - dpphi(iz,ir,1,ist)*phi(iz,ir,2,ist)  &
            - dpphi(iz,ir,2,ist)*phi(iz,ir,1,ist)  &
            - drphi(iz,ir,2,ist)*phi(iz,ir,1,ist)  &
            + drphi(iz,ir,1,ist)*phi(iz,ir,2,ist) )
        
      END DO
    END DO
  END DO
  
!        set the differential MF-densities to zero
!        (or to 2*tau -- in case of the laplacian of rho)
!        ---------------------------------------------------------------
  DO ir=0,nr
    DO iz=0,nz
      mfd_grhor(iz,ir,iiso) = 0D0
      mfd_grhoz(iz,ir,iiso) = 0D0
      mfd_drho(iz,ir,iiso)  = 2.0D0*mfd_tau(iz,ir,iiso)
      mfd_divj(iz,ir,iiso)  = 0D0
    END DO
  END DO
  
!        accumulate the differential mean-field densities
!        ---------------------------------------------------------------
  DO ist=istart(iiso),iend(iiso)
    occ = 4.0D0*occphi(ist)
    DO ir=0,nr
      DO iz=0,nz
        
        mfd_divj(iz,ir,iiso)  = mfd_divj(iz,ir,iiso)  + occ*(  &
            + drphi(iz,ir,1,ist)*dzphi(iz,ir,2,ist)  &
            - drphi(iz,ir,2,ist)*dzphi(iz,ir,1,ist)  &
            + drphi(iz,ir,1,ist)*dpphi(iz,ir,1,ist)  &
            - drphi(iz,ir,2,ist)*dpphi(iz,ir,2,ist)  &
            - dpphi(iz,ir,1,ist)*dzphi(iz,ir,2,ist)  &
            - dpphi(iz,ir,2,ist)*dzphi(iz,ir,1,ist)     )
        
        mfd_drho(iz,ir,iiso)  = mfd_drho(iz,ir,iiso)  + occ*(  &
            ddphi(iz,ir,1,ist)*phi(iz,ir,1,ist)  &
            + ddphi(iz,ir,2,ist)*phi(iz,ir,2,ist) )
        
        mfd_grhor(iz,ir,iiso) = mfd_grhor(iz,ir,iiso) + occ*(  &
            drphi(iz,ir,1,ist)*phi(iz,ir,1,ist)  &
            + drphi(iz,ir,2,ist)*phi(iz,ir,2,ist) )
        
        mfd_grhoz(iz,ir,iiso) = mfd_grhoz(iz,ir,iiso) + occ*(  &
            dzphi(iz,ir,1,ist)*phi(iz,ir,1,ist)  &
            + dzphi(iz,ir,2,ist)*phi(iz,ir,2,ist) )
        
      END DO
    END DO
  END DO
  
END DO  ! END OF ISOSPIN LOOP

!     recouple the Skyrme-Parameters
!     ------------------------------------------------------------------
b0  =                t0 + 0.5D0*t0*x0
b0p =          0.5D0*t0 +       t0*x0
b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
b2  = 0.125D0*(3.0D0*t1 + 1.5D0*t1*x1 -       t2 - 0.5D0*t2*x2)
b2p = 0.125D0*(1.5D0*t1 + 3.0D0*t1*x1 +   0.5*t2 +       t2*x2)
b3  = 0.250D0*(      t3 + 0.5D0*t3*x3)
b3p = 0.250D0*(0.5D0*t3 +       t3*x3)
b4  = 0.500D0*       t4
b5  = 0.125D0*(   x1*t1 +       x2*t2) * riftm
b5p = 0.125D0*(      t1 -          t2) * riftm

wgt_b0  =  (b0-b0p)
wgt_b0c =   b0
wgt_b1  =  (b1-b1p)
wgt_b1c =   b1
wgt_b2  = -(b2-b2p)
wgt_b2c = - b2
IF (b4p == 0D0) THEN    ! no isotopic spin-orbit force
  wgt_b4  =   2.0D0*b4
  wgt_b4c =         b4
ELSE
  wgt_b4  =  (b4+b4p)
  wgt_b4c =   b4
END IF
IF(riftm /= 0D0) THEN
!      if(ifls.eq.1) then
  wgt_b5  = -(b5+b5p)
  wgt_b5c = - b5
ELSE
  wgt_b5  = 0D0
  wgt_b5c = 0D0
END IF

!     set the MF-potentials to zero
!     ------------------------------------------------------------------
DO iiso=1,2
  DO ir=0,nr
    DO iz=0,nz
      mfp_u(iz,ir,iiso)   = 0D0
      mfp_b(iz,ir,iiso)   = h2mav
      mfp_dzb(iz,ir,iiso) = 0D0
      mfp_drb(iz,ir,iiso) = 0D0
      mfp_wz(iz,ir,iiso)  = 0D0
      mfp_wr(iz,ir,iiso)  = 0D0
    END DO
  END DO
END DO

!     accumulate the non linear contributions of the potential U
!     ---------------------------------------------------------------
fac1   =  b3  * (power+2.0D0)/3.0D0 -b3p *  power       /3.0D0
fac2   = -b3p *        2.0D0 /3.0D0
fac3   = +b3p *  power*2.0D0 /3.0D0
DO ir=0,nr
  DO iz=0,nz
    buff1 = mfd_rho(iz,ir,1)
    buff2 = mfd_rho(iz,ir,2)
    buff3 = mfd_rho(iz,ir,1) + mfd_rho(iz,ir,2)
    IF (buff3 /= 0D0) THEN
      buff1 = buff1/buff3
      buff2 = buff2/buff3
      buff3 = buff3**(power+1.0D0)
      mfp_unl(iz,ir,1) = buff3*( fac1 + fac2*buff1 + fac3*buff1*buff2 )
      mfp_unl(iz,ir,2) = buff3*( fac1 + fac2*buff2 + fac3*buff1*buff2 )
    ELSE
      mfp_unl(iz,ir,1) = 0D0
      mfp_unl(iz,ir,2) = 0D0
    END IF
  END DO
END DO

!     the coulomb-exchange in Slater approximation
!     ------------------------------------------------------------------
IF(ifex == 1) THEN
  fac1   =   -1.418033D0*chgx
  fac2   =    1.0D0/3.0D0
  DO ir=0,nr
    DO iz=0,nz
      IF (mfd_rho(iz,ir,1) /= 0D0) THEN
        mfp_unl(iz,ir,1) = mfp_unl(iz,ir,1) + fac1*mfd_rho(iz,ir,1)**(fac2)
      END IF
    END DO
  END DO
END IF

!     accumulate effective mass B and its gradients components
!     ------------------------------------------------------------------
DO iiso=1,2
  DO ir=0,nr
    DO iz=0,nz
      mfp_b(iz,ir,iiso)   = mfp_b(iz,ir,iiso)           +  &
          wgt_b1 *mfd_rho(iz,ir,  iiso) + wgt_b1c*mfd_rho(iz,ir,3-iiso)
      
      mfp_drb(iz,ir,iiso) = mfp_drb(iz,ir,iiso)         +  &
          wgt_b1 *mfd_grhor(iz,ir,  iiso) + wgt_b1c*mfd_grhor(iz,ir,3-iiso)
      
      mfp_dzb(iz,ir,iiso) = mfp_dzb(iz,ir,iiso)         +  &
          wgt_b1 *mfd_grhoz(iz,ir,  iiso) + wgt_b1c*mfd_grhoz(iz,ir,3-iiso)
      
      mfp_u(iz,ir,iiso)   = mfp_u(iz,ir,iiso)  &
          + wgt_b0 *mfd_rho(iz,ir,   iiso) + wgt_b0c*mfd_rho(iz,ir, 3-iiso)  &
          + wgt_b1 *mfd_tau(iz,ir,   iiso) + wgt_b1c*mfd_tau(iz,ir, 3-iiso)  &
          + wgt_b2 *mfd_drho(iz,ir,  iiso) + wgt_b2c*mfd_drho(iz,ir,3-iiso)  &
          - wgt_b4 *mfd_divj(iz,ir,  iiso) - wgt_b4c*mfd_divj(iz,ir,3-iiso)  &
          + mfp_unl(iz,ir,iiso)
      
      mfp_wz(iz,ir,iiso)  = mfp_wz(iz,ir,iiso)  &
          + wgt_b4 *mfd_grhoz(iz,ir,   iiso) + wgt_b4c*mfd_grhoz(iz,ir, 3-iiso)  &
          + wgt_b5 *mfd_jz(iz,ir,      iiso) + wgt_b5c*mfd_jz(iz,ir,    3-iiso)
      
      mfp_wr(iz,ir,iiso)  = mfp_wr(iz,ir,iiso)  &
          + wgt_b4 *mfd_grhor(iz,ir,   iiso) + wgt_b4c*mfd_grhor(iz,ir, 3-iiso)  &
          + wgt_b5 *mfd_jr(iz,ir,      iiso) + wgt_b5c*mfd_jr(iz,ir,    3-iiso)
    END DO
  END DO
END DO  ! END OF ISOSPIN-LOOP

!     solve the coulomb problem (or copy it from the main routine)
!     ------------------------------------------------------------------
IF (chg /= 0D0) THEN
  DO ir=0,nr
    DO iz=0,nz
      mfp_u(iz,ir,1) = mfp_u(iz,ir,1) !+ CPcoupot(iz,ir)
    END DO
  END DO
END IF

RETURN
END SUBROUTINE upd_static_mf

!-----div   -----------------------------------------------------------

SUBROUTINE div(divout,fldrin,fldzin)

!     divergence of vector field:
!        fldrin  =  r-component of vector field, input
!        fldzin  =  z-component of vector field, input
!        divout  =  divergence, output

USE paramr
USE grid
USE wavefr, ONLY: wrk1,wrk2
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'wavefr.inc'


REAL(8), INTENT(OUT)                        :: divout(0:nzd,0:nrd)
REAL(8), INTENT(IN)                         :: fldrin(0:nzd,0:nrd)
REAL(8), INTENT(IN OUT)                     :: fldzin(0:nzd,0:nrd)




!-----------------------------------------------------------


!---- calculate derivatives of fields

CALL nablaz_scalar(fldzin(0,0),wrk1(0,0), 1,1,one)
CALL nablar_scalar(fldrin(0,0),wrk2(0,0), 0,1,one)

!---- value at 'r=0'

DO iz=0,nz
  divout(iz,0) = 2.0*wrk2(iz,0)+wrk1(iz,0)
END DO

!---- all other 'r'

DO ir=1,nr
  r      = (sr+ir)*dr
  rinv = one/r
  DO iz=0,nz
    divout(iz,ir) = wrk2(iz,ir)+rinv*fldrin(iz,ir)+wrk1(iz,ir)
  END DO
END DO

RETURN
END SUBROUTINE div

SUBROUTINE sort_states()

! energy sorting of wavefunctions and correspdonding occupations
! in the 'm' blocks.

USE paramr
USE wavefr
USE levels
USE fields
USE cranking, ONLY: dpsir


IMPLICIT NONE

LOGICAL,PARAMETER :: ttest=.true.

REAL(8), ALLOCATABLE :: esav(:)
INTEGER, ALLOCATABLE :: isort(:),isav(:)

REAL(8) :: emin,xe,fpair
INTEGER :: nstates,i,jst,ist
INTEGER :: im,nhi,nlo,nadd,kpsic,iq

nstates=nst(1)+nst(2)
ALLOCATE(esav(1:nstates),isort(1:nstates),isav(1:nstates))

esav(1:nstates)=spenrg(1:nstates)
DO i=1,nstates; isort(i)=i; END DO

! the outer loop over 'm'&'iq' blocks
DO  im=1,nom         ! loop over all blocks
  nhi=mhi(im)        ! uppermost state for given 'im'
  nlo=mlo(im)        ! lowest state for given 'im''
  nadd=nhi-nlo+1     ! number of states in 'im' block
  kpsic=mim(im)      ! actual 'm' value for 'im,'
  iq=miq(im)         ! actual isospin 'iq' for given 'im'

  ! evalute pointer field 'isort'
  DO  ist=nlo,nhi
    emin = esav(ist)
    DO jst=ist,nhi
      IF(esav(jst)<emin) THEN
        emin=esav(jst)
        esav(jst)=esav(ist)
        esav(ist)=emin
        i=isort(jst)
        isort(jst)=isort(ist)
        isort(ist)=i
      END IF
    END DO
  END DO
  IF(ttest) THEN
    WRITE(6,'(a,5i5)') ' for block: im,iq,kpsic,nlo,nhi=',&
          im,iq,kpsic,nlo,nhi
    WRITE(6,'(a,50(1pg13.5))') 'spenerg:',spenrg(nlo:nhi)
    WRITE(6,'(a,50(1pg13.5))') 'sorted: ',esav(nlo:nhi)
    WRITE(6,'(a,50i5)') 'isort :',isort(nlo:nhi)
  END IF

  ! map the arrays to new order
  esav(nlo:nhi)=spenrg(nlo:nhi)
  DO ist=nlo,nhi; spenrg(ist)=esav(isort(ist)); END DO
  IF(ttest) WRITE(6,'(a,50(1pg13.5))') 'spenerg:',spenrg(nlo:nhi)
  esav(nlo:nhi)=wate(nlo:nhi)
  DO ist=nlo,nhi; wate(ist)=esav(isort(ist)); END DO
  esav(nlo:nhi)=spn(nlo:nhi)
  DO ist=nlo,nhi; spn(ist)=esav(isort(ist)); END DO
  esav(nlo:nhi)=spevar(nlo:nhi)
  DO ist=nlo,nhi; spevar(ist)=esav(isort(ist)); END DO
  esav(nlo:nhi)=vsqf(nlo:nhi)
  DO ist=nlo,nhi; vsqf(ist)=esav(isort(ist)); END DO
  esav(nlo:nhi)=deltaf(nlo:nhi)
  DO ist=nlo,nhi; deltaf(ist)=esav(isort(ist)); END DO
  esav(nlo:nhi)=wtpair(nlo:nhi)
  DO ist=nlo,nhi; wtpair(ist)=esav(isort(ist)); END DO
  isav(nlo:nhi)=iorder(nlo:nhi)
  DO ist=nlo,nhi; iorder(isort(ist))=isav(ist); END DO
  isav(nlo:nhi)=nparz(1,nlo:nhi)
  DO ist=nlo,nhi; nparz(1,ist)=isav(isort(ist)); END DO
  isav(nlo:nhi)=nparz(2,nlo:nhi)
  DO ist=nlo,nhi; nparz(2,ist)=isav(isort(ist)); END DO

  dpsir(0:nz,0:nr,1:nsd,nlo:nhi) = psi(0:nz,0:nr,1:nsd,nlo:nhi)
  DO ist=nlo,nhi
    psi(0:nz,0:nr,1:nsd,ist) = dpsir(0:nz,0:nr,1:nsd,isort(ist))
  END DO
  
!end the outer loop  
END DO

! refresh 'iorder'
esav=spenrg
DO i=1,nstates; iorder(i)=i; END DO

DO  ist=1,nst(1)
  emin = esav(ist)
  DO jst=ist,nst(1)
    IF(esav(jst)<emin) THEN
      emin=esav(jst)
      esav(jst)=esav(ist)
      esav(ist)=emin
      i=iorder(jst)
      iorder(jst)=iorder(ist)
      iorder(ist)=i
    END IF
  END DO
END DO
iorder(1:nst(1))=iorder(1:nst(1))-iokk(1)

DO  ist=nst(1)+1,nst(2)
  emin = esav(ist)
  DO jst=ist,nst(2)
    IF(esav(jst)<emin) THEN
      emin=esav(jst)
      esav(jst)=esav(ist)
      esav(ist)=emin
      i=iorder(jst)
      iorder(jst)=iorder(ist)
      iorder(ist)=i
    END IF
  END DO
END DO
iorder(nst(1)+1:nst(2))=iorder(nst(1)+1:nst(2))-iokk(2)


! new densities
rho=0D0;tau=0D0;drhols=0D0;rholsr=0D0;rholsz=0D0
DO jst=1,nstates
  fpair  = vsqf(jst)*wate(jst)
  IF(jst.LE.nst(1)) THEN
    iq=1
  ELSE
    iq=2
  END IF
  CALL densit(psi(0,0,1,jst),vsqf(jst),jst,  &
              rho(0,0,iq),tau(0,0,iq),drhols(0,0,iq),  &
              rholsr(0,0,iq),rholsz(0,0,iq) )
END DO

DEALLOCATE(esav,isort,isav)


END SUBROUTINE sort_states


!-----pri_2D -----------------------------------------------------------

SUBROUTINE pri_2d()

!**********************************************************************
!                                                                     *
!     print 2D info of densities                                      *
!                                                                     *
!**********************************************************************

USE paramr
USE dynpar
USE fields
USE observ
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'


!CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam
CHARACTER (LEN=3) :: outnumber

INTEGER,SAVE :: icallno=0

!----------------------------------------------------------------------

icallno = 1+icallno
OPEN(77,STATUS='scratch')
WRITE(77,'(i4)') icallno+1000
REWIND(77)
READ(77,'(1x,a3)') outnumber
CLOSE(77)
OPEN(77,FILE='fields_2D.'//outnumber//'.'//outnam)
WRITE(77,'(a,2i5,a,f6.0//a,3f7.3//a,1pg14.6,1pg12.4)')  &
    '# Nr. proton/neutron=',npart, ' total Nr. Nucleons=',anumb,  &
    '# deformations beta2,beta3,beta4=',beta2,beta3,beta4,  &
    '# energy, rms-radius=',energy,rms

DO ir=1,nr
  r   = (sr+ir)*dr
  DO iz=0,nz
    z     = (sz+iz-icm)*dz
    WRITE(77,'(2i5,2f7.2,3(1pg12.4))') ir,iz,r,z,  &
        rho(iz,ir,1),rho(iz,ir,2),rho(iz,ir,1)+rho(iz,ir,2)
  END DO
  WRITE(77,'(1x)')
END DO

CLOSE(77)

RETURN
END SUBROUTINE pri_2d


!======================================================================

!     wrt_2D_spin(fldin,filename,label,incridx,outnam)

!     print 2D info of a spinor field (wavefunction, or
!     isospin-dependent fields)

!=======================================================================

SUBROUTINE wrt_2d_spin(fldin,filename,incridx)

!----------------------------------------------------------------------

USE paramr
USE dynpar
USE fields
USE observ
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!----------------------------------------------------------------------


DOUBLE PRECISION, INTENT(IN OUT)         :: fldin(0:nzd,0:nrd,2)
CHARACTER (LEN=10), INTENT(IN)           :: filename
INTEGER, INTENT(IN)                      :: incridx
!CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam

CHARACTER (LEN=3) :: outnumber

CHARACTER (LEN=80) :: label

INTEGER :: nfld,ir,iz
DOUBLE PRECISION :: r,z
DATA             icallno/0/

!----------------------------------------------------------------------

!     generate a string containing the index of the file
!----------------------------------------------------------------------
IF (incridx == 1) THEN    ! increase index if incridx==1
  icallno = 1 + icallno
END IF

OPEN(77,STATUS='scratch')
WRITE(77,'(i4)') icallno + 1000
REWIND(77)
READ(77,'(1x,a3)') outnumber
CLOSE(77)

!     open file and write header
!----------------------------------------------------------------------
OPEN(77,FILE=filename//'.'//outnumber//'.'//outnam)
WRITE(77,'(a,2i5,a,f6.0/a,3f7.3/a,1pg14.6,1pg12.4)')  &
    '# Nr. proton/neutron=',npart, ' total Nr. Nucleons=',anumb,  &
    '# deformations beta2,beta3,beta4=',beta2,beta3,beta4,  &
    '# energy, rms-radius=',energy,rms

!     write the field (could either be a spinor or a isospinor)
!----------------------------------------------------------------------
DO ir=0,nr
  r   = (sr+ir)*dr
  DO iz=0,nz
    z     = (sz+iz-icm)*dz
    WRITE(77,'(2i5,2f7.2,2(1pg16.5))') ir,iz,r,z,  &
        fldin(iz,ir,1),fldin(iz,ir,2)
  END DO
  WRITE(77,'(1x)')
END DO

CLOSE(77)

!----------------------------------------------------------------------

RETURN
END SUBROUTINE wrt_2d_spin

!======================================================================

!     wrt_2D_scalar(fldin,filename,label,incridx,outnam)

!     print 2D info of a scalar field (like coulomb density)

!=======================================================================

SUBROUTINE wrt_2d_scalar(fldin,filename,incridx)

!----------------------------------------------------------------------

USE paramr
USE dynpar
USE fields
USE observ
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!----------------------------------------------------------------------


DOUBLE PRECISION, INTENT(IN OUT)         :: fldin(0:nzd,0:nrd)
CHARACTER (LEN=10), INTENT(IN)           :: filename
INTEGER, INTENT(IN)                      :: incridx
!CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam

CHARACTER (LEN=3) :: outnumber

CHARACTER (LEN=80) :: label

INTEGER :: nfld,ir,iz
DOUBLE PRECISION :: r,z
DATA             icallno/0/

!----------------------------------------------------------------------

!     generate a string containing the index of the file
!----------------------------------------------------------------------
IF (incridx == 1) THEN    ! increase index if incridx==1
  icallno = 1 + icallno
END IF

OPEN(77,STATUS='scratch')
WRITE(77,'(i4)') icallno + 1000
REWIND(77)
READ(77,'(1x,a3)') outnumber
CLOSE(77)

!     open file and write header
!----------------------------------------------------------------------
OPEN(77,FILE=filename//'.'//outnumber//'.'//outnam)
WRITE(77,'(a,2i5,a,f6.0/a,3f7.3/a,1pg14.6,1pg12.4)')  &
    '# Nr. proton/neutron=',npart, ' total Nr. Nucleons=',anumb,  &
    '# deformations beta2,beta3,beta4=',beta2,beta3,beta4,  &
    '# energy, rms-radius=',energy,rms

!     write the scalar field
!----------------------------------------------------------------------
DO ir=0,nr
  r   = (sr+ir)*dr
  DO iz=0,nz
    z     = (sz+iz-icm)*dz
    WRITE(77,'(2i5,2f7.2,1pg16.5)') ir,iz,r,z,  &
        fldin(iz,ir)
  END DO
  WRITE(77,'(1x)')
END DO

CLOSE(77)

!----------------------------------------------------------------------

RETURN
END SUBROUTINE wrt_2d_scalar


!-----ondiag-----------------------------------------------------------

SUBROUTINE ondiag(psih,jst,nparps,istadd,nhi)

!**********************************************************************
!                                                                     *
!        ondiag - subroutine to force the diagonalization amongst     *
!                 the occupied states.                                *
!                 Input h*psi for state nr. 'istadd' on 'psih'.       *
!                 Computes all h-overlaps up to 'nhi' and performs    *
!                 a perturbative correction towards diagonality.      *
!                 'jst' is the nr. of the state associated to 'psih'  *
!                 in absolute addres mode.                            *
!                                                                     *
!**********************************************************************

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'

REAL(8), INTENT(IN)                         :: psih(0:nzd,0:nrd,nsd)
INTEGER, INTENT(IN)                      :: jst
INTEGER, INTENT(IN OUT)                  :: nparps
INTEGER, INTENT(IN)                      :: istadd
INTEGER, INTENT(IN)                      :: nhi

REAL(8),PARAMETER :: stepmx=0.1D0      ! maximum admixture of other states

stepac = 0.5*stepmx/pi

!**********************************************************************
!                                                                     *
!     first loop (do 2000) runs over the states to be projected.      *
!                                                                     *
!**********************************************************************

ist    = istadd
DO  in=1,nhi
  ist     = ist+1
  npari   = npar(ist)
!        test if parity is a good quantum number select those states  *
  IF((nparps == npari .OR. sym == zero).AND. ist /= jst) THEN
    
!     compute overlap
    
    ovlp   = wfovlp(psi(0,0,1,ist),psih(0,0,1))
    corfac = stepac*ATAN(ovlp/((spenrg(jst)-spenrg(ist))*stepac))
    
!     extract psi(in)-component in psih
    
    DO is=1,nsd
      DO ir=0,nr
        DO iz=0,nz
          psi(iz,ir,is,jst)=psi(iz,ir,is,jst) +corfac*psi(iz,ir,is,ist)
        END DO
      END DO
    END DO
    
  END IF
  
END DO

RETURN
END SUBROUTINE ondiag

SUBROUTINE corzcm(zcm)

!     shifts the wavefunctions such that the center-of-mass tends
!     to zero.
!     zcm  = c.m. at input, corrected c.m. at output

USE paramr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'



REAL(8), INTENT(IN OUT)                    :: zcm

REAL(8),PARAMETER :: addzcm=-0.5D0
LOGICAL,PARAMETER :: ttest=.false.


!----------------------------------------------------------------------

WRITE(6,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
WRITE(7,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
stepfc = zcm*addzcm
DO iq=1,2
  DO ir=0,nr
    DO iz=0,nz
      rhoti(iz,ir,iq) = zero
    END DO
  END DO
END DO

!     loop over wavefunctions in blocks of same quantum numbers

DO im=1,nom
  
  nhi=mhi(im)
  nlo=mlo(im)
  nadd=nhi-nlo+1
  kpsic=mim(im)
  iq=miq(im)
  
!       shift by 'stepfc'
  
  IF(ttest) THEN
    DO jst=nlo,nhi
      CALL nablaz_spin(psi(0,0,1,jst),psi1,nparz(1,jst),1,one)
      DO is=1,2
!PK          call mulz(derz1(0,0,nparz(is,jst)),psi(0,0,is,jst),
!PK     &                wrk4,one,0)
        DO ir=0,nr
          DO iz=0,nz
            psi(iz,ir,is,jst) = psi(iz,ir,is,jst) - stepfc*psi1(iz,ir,is)
          END DO
        END DO
      END DO
    END DO
    
!      Schmidt orthonormmalization of the states
    
    CALL schmit(psi(0,0,1,nlo),nlo-1,nadd)
  END IF
  
!       new density
  
  DO jst=nlo,nhi
    CALL sdenti(psi(0,0,1,jst),rhoti(0,0,iq),vsqf(jst),jst)
  END DO
  
END DO

zcm = zcmobs(rhoti)
WRITE(6,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
WRITE(7,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm

RETURN
END SUBROUTINE corzcm


!---------------------------------------------------------------------------  
! DESCRIPTION: nablap_scalar
!> @brief
!! Calculate the \f$\phi\f$-component of the gradient of a scalar function.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] imfld
!> INTEGER, takes the m-component of the field.
!> @param[in] ippfld
!> INTEGER, takes the p-parity of the field (cos(m*phi)=+1, sin(m*phi)=-1).
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 
SUBROUTINE nablap_scalar(fldin,fldout,imfld,ippfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr)
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ippfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: irpidx,ibpidx

irpidx = 1+MOD(ABS(imfld),2)
ibpidx = irpidx + 2*ibpfld - 2

CALL mulr(r1inv(0:nr,0:nr), fldin(0:nz,0:nr),fldout(0:nz,0:nr),-factor*imfld*ippfld,0)
IF (imfld == 1.OR.imfld == -1) &
  CALL mulr0(derr1(0:nr,0:nr,ibpidx), fldin(0:nz,0:nr),fldout(0:nz,0:nr),-factor*imfld*ippfld,0)


RETURN
END SUBROUTINE nablap_scalar

!---------------------------------------------------------------------------  
! DESCRIPTION: colin3
!> @brief
!!Initialises the monpole and quadrupole background fields for the 
!!subroutine colrun, here using solution by radial integration of the 
!!Green's function.
!>
!> @param[out] wcol
!> REAL(8), array, returns the Coulomb potential.
!> @param[in] mangmo
!> REAL(8), array, takes the z-projection of the angular momentum.
!---------------------------------------------------------------------------
SUBROUTINE colin3(wcol,mangmo)

USE paramr
USE fields
USE levels
USE grid
USE dynpar
IMPLICIT NONE

REAL(8), INTENT(OUT)              :: wcol(0:nz,0:nr)
INTEGER, INTENT(IN)               :: mangmo

REAL(8), PARAMETER :: rs=1.16D0 ! Wigner-Seitz radius for nuclei in units [fm]

LOGICAL :: tcase(0:4)           ! Switches to required cases  l=0..4

INTEGER :: lpoint,l,ir,iz

REAL(8) :: agen,betacq,binup,cqdip,cqhex,cqmono,cqoct,cqquad,fzrm,rad,rad2,&
           rmax,rmsq,rnorm,vnorm,xnorm,ylmp,zmax,zrmax,zrmin

! Initialisation of wcol
wcol(0:nz,0:nr) = 0.0

! General parameters for background fields
zmax    = ((nz)+sz)*dz
rmax    = ((nr)+sr)*dr
zrmin   = MIN(zmax,rmax)
zrmax   = MAX(zmax,rmax)
binup   = zrmax+5.0D0
fzrm    = zrmin/5.5D0
rad0    = rs*anumb**edr

agen    = SQRT( -LOG( oh ) )/fzrm
agen    = one/MIN(rad0*oh,one/agen)
agen    = agen*agen
DO l=0,4
  aexp(l) = agen
  IF(l >= 1) THEN
    bdiv(l)  = zrmax**l
    rshif(l) = rad0
  ELSE
    rshif(0) = zero
  END IF
END DO

tcase(0) = .true.
tcase(1) = ndrun == 1 .AND. isym /= 1
tcase(2) = nqrun == 1
tcase(3) = isym == 0  .AND. norun == 1
tcase(4) = nhrun == 1

!     initialize the multipole fields
DO  l=0,4
  lpow  = l
  lpoint = l+l
  IF(tcase(l)) THEN
!                                      normalization
    radpow = 2*l+2
    CALL qgaus(refrho, zero, binup, xnorm)
    rnorm   = one/xnorm
    vnorm   = rnorm*4.0*pi*chg/(l+l+one)
    
    DO ir=0,nr
      DO iz=0,nz
        rad2  = (z(iz)+icm*dz)**2+r(ir)**2
        rad   = SQRT(rad2)
        ylmp  = ylm((z(iz)+icm*dz),r(ir),l,mangmo)
!                                      density
        radpow = zero
        colref(iz,ir,mangmo,lpoint) = rnorm*refrho(rad)*ylmp
!                                      potential
        IF (rad < 1.d-2) THEN
          t1 = zero
        ELSE
          radpow = 2*l+2
          CALL qgaus(refrho, zero, rad, t1)
          t1 = t1/rad**(l+l+one)
        END IF
        radpow = one
        CALL qgaus(refrho, rad , binup, t2)
        colref(iz,ir,mangmo,lpoint+1) = vnorm*(t1+t2)*ylmp
      END DO
    END DO
    
    CALL obs3(colref(0:nz,0:nr,mangmo,lpoint),cqmono,rmsq,cqdip,  &
        cqquad,cqoct,cqhex,betacq,mangmo)
    WRITE(6,'(a,i2,a/5x,8(1pg11.3))') ' from "l=',l,  &
        '": cqmono,cqdip,cqquad,cqoct,cqhex,rmsq,betacq=',  &
        cqmono,cqdip,cqquad,cqoct,cqhex,rmsq,betacq
    IF(l == 0) THEN
      cmono = cqmono
    ELSE IF(l == 1) THEN
      cdip  = cqdip
    ELSE IF(l == 2) THEN
      cquad = cqquad
    ELSE IF(l == 3) THEN
      coct  = cqoct
    ELSE IF(l == 4) THEN
      chex  = cqhex
    END IF

    IF(itest >= 1) THEN
      CALL test_print(colref(0:nz,0:nr,mangmo,lpoint+1),0,nz,0,9,' Coul.')
      CALL test_print(colref(0:nz,0:nr,mangmo,lpoint),0,nz,0,9,' Dens.')
    END IF
  ELSE
    colref(0:nz,0:nr,mangmo,lpoint:lpoint+1) = zero
  END IF
END DO

RETURN
END SUBROUTINE colin3


!---------------------------------------------------------------------------  
! DESCRIPTION: obs3
!> @brief
!!This routine computes observables from an input density such as particle 
!!number, radius and moments. 
!> @detail
!!The properties are evaluated as integrals.
!>
!> @param[in] rhoini
!> REAL(8), array, takes the input density.
!> @param[out] atnum
!> REAL(8), returns the particle number.
!> @param[out] rms
!> REAL(8), returns the rms radius.
!> @param[out] d
!> REAL(8), returns the spherical dipole moment.
!> @param[out] q
!> REAL(8), returns the cartesian quadrupole moment.
!> @param[out] o
!> REAL(8), returns the spherical octupole moment.
!> @param[out] h
!> REAL(8), returns the spherical hexadecupole moment.
!> @param[out] qsq
!> REAL(8), returns the squared quadrupole moment as an estimate of variance.
!> @param[in] mangmo
!> REAL(8), takes the z-projection of the angular momentum.
!--------------------------------------------------------------------------- 
SUBROUTINE obs3(rhoini,atnum,rms,d,q,o,h,qsq,mangmo) !TODO Can we just use obs instead???

USE paramr
USE fields
USE grid, ONLY: nr,nz,dr,dvol,dz,icm,isym,sr,sym,sz

IMPLICIT NONE

REAL(8), INTENT(IN)       :: rhoini(0:nz,0:nr)
REAL(8), INTENT(OUT)      :: atnum
REAL(8), INTENT(OUT)      :: rms
REAL(8), INTENT(OUT)      :: d
REAL(8), INTENT(OUT)      :: q
REAL(8), INTENT(OUT)      :: o
REAL(8), INTENT(OUT)      :: h
REAL(8), INTENT(OUT)      :: qsq
INTEGER, INTENT(IN)       :: mangmo


REAL(8),PARAMETER :: three=3D0,seven=7D0,eight=8D0,twfour=24D0

INTEGER ir,iz

REAL(8) :: at1,at2,at3,d1,d2,d3,dir,drezw,h1,h2,h3,ht,o1,o2,o3,ot,&
           q1,q2,q2t,q3,qsq1,qsq2,qsq3,qt,r,r2,r2t,rhot,rms1,rms2,&
           rms3,vol,z,z2,zcm

IF(mangmo /= 0) STOP 'OBS3 not yet prepared for MANGMO.NE.0'

!     summing up the wavefunction

drezw  = one/12.0
zcm    = 0.0D0            ! set default for Coulomb solver

at1    = rhoini(0,0)
IF (isym == 1) at1 = at1 * oh
rms1   = zero
d1     = zero
q1     = zero
qsq1   = zero
o1     = zero
h1     = zero

DO iz=1,nz
  z     = ((iz)+sz-icm)*dz-zcm
  z2    = z*z
  r2t   = z2
  qt    = z2 + z2
  q2t   = qt*qt
  ot    = two * z * z2
  ht    = eight*z2*z2
  rhot  = rhoini(iz,0)
  at1   = rhot     + at1
  rms1  = rhot*r2t + rms1
  d1    = rhot*z   + d1
  q1    = rhot*qt  + q1
  qsq1  = rhot*q2t  + qsq1
  o1    = rhot*ot  + o1
  h1    = rhot*ht  + h1
END DO
d1     = d1 - rhot*z
o1    =  o1 - rhot*ot
at1    = drezw * at1
rms1   = drezw * rms1
d1     = drezw * d1
q1     = drezw * q1
qsq1   = drezw * qsq1
o1     = drezw * o1
h1     = drezw * h1

at2    = zero
rms2   = zero
d2     = zero
q2     = zero
qsq2   = zero
o2     = zero
h2     = zero
z      = (sz) * dz - zcm
DO ir=1,nr
  dir   = (ir)+sr
  z2    = z*z
  r     = dir*dr
  r2    = r*r
  r2t   = r2 + z2
  qt    = z2 + z2 - r2
  q2t   = qt*qt
  ot    = two*z*z2 - three*z*r2
  ht    = eight*z2*z2-twfour*z2*r2+three*r2*r2
  rhot  = rhoini(0,ir)*dir
  at2   = rhot     + at2
  rms2  = rhot*r2t + rms2
  d2    = rhot*z   + d2
  q2    = rhot*qt  + q2
  qsq2  = rhot*q2t  + qsq2
  o2    = rhot*ot  + o2
  h2    = rhot*ht  + h2
END DO
IF (isym == 1) THEN
  at2    = oh * at2
  rms2   = oh * rms2
  d2     = oh * d2
  q2     = oh * q2
  qsq2   = oh * qsq2
  o2     = oh * o2
  h2     = oh * h2
END IF

at3    = zero
rms3   = zero
d3     = zero
q3     = zero
qsq3   = zero
o3     = zero
h3     = zero
DO ir=1,nr
  dir   = (ir)+sr
  r     = dir*dr
  r2    = r*r
  DO iz=1,nz
    z     = ((iz)+sz-icm)*dz-zcm
    z2    = z*z
    r2t   = z2+r2
    qt    = z2+z2-r2
    q2t   = qt*qt
    ot    = two*z*z2 - three*z*r2
    ht    = eight*z2*z2-twfour*z2*r2+three*r2*r2
    rhot  = rhoini(iz,ir)*dir
    at3   = rhot     + at3
    rms3  = rhot*r2t + rms3
    d3    = rhot*z   + d3
    q3    = rhot*qt  + q3
    qsq3  = rhot*q2t  + qsq3
    o3    = rhot*ot  + o3
    h3    = rhot*ht  + h3
  END DO
  rhot  = rhoini(nz,ir)*dir
  d3    = d3 - rhot*z
  o3    =  o3 - rhot*ot
END DO

!     the different parts are summed up:

atnum   =  at1+at2+at3
rms     =  rms1 +rms2 +rms3
d       =  d1 +d2 +d3
q       =  q1 +q2 +q3
qsq     =  qsq1 +qsq2 +qsq3
o       =  o1 +o2 +o3
h       =  h1 +h2 +h3

!     the final result

rms   = SQRT(ABS(rms/MAX(atnum,smal)))

vol   = (sym+one)*dvol

atnum = atnum*vol

q     = q  * vol     !!!!!!!!!!!!!
qsq   = qsq* vol
h     = h  * vol * three/(two*eight*SQRT(pi))

IF (isym == 1) THEN
  d = zero
  o = zero
ELSE
  d = d * vol * oh * SQRT(three/pi)
  o = o * vol * oh * oh * SQRT(seven/pi)
END IF
WRITE(*,*) atnum
WRITE(*,*) rms
WRITE(*,*) d
WRITE(*,*) q
WRITE(*,*) o
WRITE(*,*) h
WRITE(*,*) qsq
STOP 'obs3'
RETURN
END SUBROUTINE obs3

!---------------------------------------------------------------------------  
! DESCRIPTION: damint
!> @brief
!!Alternative routine to evluate distance to surface. Defines surface line 
!!by linear interpolation and searches distance to this line.
!> @details
!!This routine is still in testing stage.
!>
!> @param[in] rhoin
!> REAL(8), array, takes the total density.
!> @param[in] rhosw
!> REAL(8), takes the rho-value defining surface.
!> @param[in] delr
!> REAL(8), returns the distance to surface.
!--------------------------------------------------------------------------- 
SUBROUTINE damint(rhoin, rhosw, delr)

USE paramr
USE dynpar
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)   :: rhoin(0:nz,0:nr)   ! total density
REAL(8), INTENT(IN)   :: rhosw         ! rho-value defining surface
REAL(8), INTENT(OUT)  :: delr(0:nz,0:nr)    ! distance to surface

INTEGER, PARAMETER :: maxp=nzd*nrd

REAL(8) :: surf(maxp,2)                    ! make ALLOCATABLE
REAL(8) :: frac, diag, maxdist, dist
INTEGER :: iz, ir, izz, irr, point, dir, p

diag = SQRT(2.0)
maxdist = 1.0E10
point = 1
DO iz=0,nz
  DO ir=0,nr
    DO  dir=1,4
      IF(dir == 1 .AND. ir /= 0 .AND. iz /= 0) THEN
        izz=iz-1
        irr=ir-1
      ELSE IF(dir == 2 .AND. iz /= 0) THEN
        izz=iz-1
        irr=ir
      ELSE IF(dir == 3 .AND. ir /= nr .AND. iz /= 0) THEN
        izz=iz-1
        irr=ir+1
      ELSE IF(dir == 4 .AND. ir /= nr) THEN
        izz=iz
        irr=ir+1
      ELSE
        CYCLE
      END IF
      IF((rhoin(iz,ir) < rhosw .AND. rhoin(izz,irr) < rhosw) .OR.  &
          (rhoin(iz,ir) >= rhosw .AND. rhoin(izz,irr) >= rhosw)) CYCLE
      frac=(rhosw-rhoin(iz,ir))/(rhoin(izz,irr)-rhoin(iz,ir))
      IF(iz == izz) THEN
        surf(point,1) = (iz+sz-icm) * dz
      ELSE IF(ir == irr) THEN
        surf(point,1) = (iz+sz-frac-icm) * dz
      ELSE
        surf(point,1) = (iz+sz-frac/diag-icm) * dz
      END IF
      IF(ir == irr) THEN
        surf(point,2) = (ir+sr) * dr
      ELSE IF(iz == izz) THEN
        surf(point,2) = (ir+sr+frac) * dr
      ELSE IF(irr < ir) THEN
        surf(point,2) = (ir+sr-frac/diag) * dr
      ELSE
        surf(point,2) = (ir+sr+frac/diag) * dr
      END IF
      point = point + 1
    END DO
  END DO
END DO

DO iz=0,nz
  DO ir=0,nr
    delr(iz,ir) = maxdist
    DO p=1,point-1
      dist = (z(iz)-surf(p,1))**2 + (r(ir)-surf(p,2))**2
      IF(delr(iz,ir) > dist) delr(iz,ir) = dist
    END DO
    delr(iz,ir) = SQRT(delr(iz,ir))
    IF(rhoin(iz,ir) > rhosw) THEN
      delr(iz,ir) = -delr(iz,ir)
    END IF
  END DO
END DO

RETURN
END SUBROUTINE damint

!---------------------------------------------------------------------------  
! DESCRIPTION: delrmn
!> @brief
!!To compute the field of minimum distance to the surface line defined by 
!!the "switch-density":  \c rho(r,z) = \c rhosw.
!> @details
!!The definition is:
!!\f[ \Delta_R^2(r,z)={\rm Min} [(r-r^\prime)^2+(z-z^\prime)^2] \f]
!!with \f$ \rho(r^\prime,r^\prime)= \texttt{rhosw} \f$
!!The minimum is searched in two steps:
!!  -# The surface line r_sw= r_sw(z) is established,
!!  -# For each point on the grid, the minimum distance to
!!     the points on r_sw is looked for. This search proceeds
!!     in strips of constant z while restarting from the minimum
!!     point found previously.
!!This routine is in testing stage and it works presently only for 
!!the case SYM=1
!!The program assumes presently the same grid size in r and z
!>
!> @param[in] rhoin
!> REAL(8), array, takes the input density.
!> @param[in] rhosw
!> REAL(8), takes the switch density.
!> @param[out] delr
!> REAL(8), array, returns the output field of minimal distances.
!--------------------------------------------------------------------------- 
SUBROUTINE delrmn(rhoin,rhosw,delr)

USE paramr
USE dynpar
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)     :: rhoin(0:nz,0:nr) ! input density
REAL(8), INTENT(IN)     :: rhosw            ! switch density
REAL(8), INTENT(OUT)    :: delr(0:nz,0:nr)  ! output field of min. distances

LOGICAL :: tifzmn                ! flag for "cap" at lower z
REAL(8) :: rsw(-1:nz)           ! r values for surface line
REAL(8) :: zsw(-1:nz)           ! z values for surface line

INTEGER :: mbase,izmin,izmax,izm,izini,izend,iz1,iz0,iz,ira,ir,idir

REAL(8) :: x,sqr1,sqr2,drsqp1,drsqmn,drsqm1,drsq0,addz


IF(dr /= dz) STOP ' subroutine "DELRMN" requires DR=DZ'
IF(itest >= 2) THEN
  WRITE(7,'(a,g12.4)') ' DELRMN called with rhosw=',rhosw
  CALL test_print(rhoin(0:nz,0:nr),0,nz,0,nr,'rhoin')
END IF

! Search surface line and store on 'rsw' with linear interpolation:
DO iz=0,nz
  DO ir=nr,0,-1
    ira  = ir
    IF(rhoin(iz,ir) > rhosw) THEN
      ira = -1
      EXIT
    END IF
  END DO
  IF(ira == -1) THEN
    rsw(iz) = -nr            ! signals that 'rhoin < rhosw' always
  ELSE
    rsw(iz) = ira +(rhoin(iz,ira)-rhosw)/(rhoin(iz,ira)-rhoin(iz,ira+1))
    zsw(iz) = iz
  END IF
END DO

! Determine minimum and maximum relevant 'iz'
DO iz=0,nz
  izmin = iz
  IF(rsw(iz) > zero) EXIT
END DO

tifzmn = izmin > 0
DO iz=nz,0,-1
  izmax = iz
  IF(rsw(iz) > zero) EXIT
END DO

! determine 'z' values of the end points, reset 'rsw'
IF(izmin > 0) THEN
  sqr1 = SQRT(rsw(izmin))
  sqr2 = SQRT(rsw(izmin+1))
  izmin = izmin-1
  addz  = MIN(one,sqr2/(sqr1+sqr2))
  zsw(izmin) = izmin+addz
  rsw(izmin) = zero
  zsw(izmin-1) = zsw(izmin+1)      ! symm. point for cap-parabola
  rsw(izmin-1) = -rsw(izmin+1)
ELSE
  izmin = izmin-1
  zsw(-1) = izmin
  IF(isym == 1)rsw(-1) = rsw(1)                 ! requires SYM=1
  IF(isym == 0)rsw(-1) = rsw(nz)                ! requires SYM=0
END IF
IF(izmax < nz) THEN
  sqr1 = SQRT(rsw(izmax))
  sqr2 = SQRT(rsw(izmax-1))
  izmax = izmax+1
  addz  = MIN(one,sqr1/(sqr2-sqr1))
  zsw(izmax) = izmax-1+addz
  rsw(izmax) = zero
  zsw(izmax+1) = zsw(izmax-1)      ! symm. point for cap-parabola
  rsw(izmax+1) = -rsw(izmax-1)
END IF
IF(itest >= 2) THEN
  WRITE(7,'(a,2i5)') ' izmin,izmax=',izmin,izmax
  WRITE(7,'(a,20(/1x,5g12.4))') ' rsw:',(rsw(iz),iz=-1,nz)
  WRITE(7,'(a,20(/1x,5g12.4))') ' zsw:',(zsw(iz),iz=-1,nz)
END IF

! Scan all grid points and determine minimum distance to 'rsw'
IF(itest >= 2) WRITE(7,'(a)') '  iz  irmbasidirizinizen iz0   drsqmn   delr'
DO iz=0,nz
  mbase = 0              ! relative base point in 'iz' for search
  DO ir=0,nr
    iz0     = MIN(MAX(izmin+1,iz+mbase),izmax-1)
    izm     = iz0-1
!                             first triplet determines direction of search
    drsq0   = (rsw(iz0)-ir)**2+(zsw(iz0)-iz)**2
    drsqm1  = (rsw(izm)-ir)**2+(zsw(izm)-iz)**2
    drsqp1  = (rsw(iz0+1)-ir)**2+(zsw(iz0+1)-iz)**2
    IF(drsqp1 > drsq0) THEN
      IF(drsqm1 < drsq0) THEN
        izini  = MIN(MAX(iz0-2,izmin),izmax)
        IF(izini >= 0) THEN
          idir   = -1
          izend  = izmin
          x      = drsqm1          ! order the "drsq.." for interpolation
          drsqm1 = drsqp1
          drsqp1 = x
        ELSE
          idir   = 0
        END IF
      ELSE
        idir   = 0             ! signals end of search
      END IF
    ELSE
      idir   = +1
      izini  = MIN(MAX(iz0+2,izmin),izmax)
      izend  = izmax
    END IF
!                             search in 'idir' until distance rises
    IF(idir /= 0) THEN
      DO iz1=izini,izend,idir
        iz0    = iz1
        drsqm1 = drsq0
        drsq0  = drsqp1
        drsqp1 = (rsw(iz0)-ir)**2+(zsw(iz0)-iz)**2
        IF(drsqp1 > drsq0) EXIT
      END DO
    END IF
    
! check outer cap-discontinuity
    IF(idir > 0 .AND. iz0 == izmax) THEN
      IF(drsqp1 < drsqm1) THEN  ! complement cap
        drsqm1 = drsq0
        drsq0  = drsqp1
        drsqp1 = (rsw(iz0+1)-ir)**2+(zsw(iz0+1)-iz)**2
      ELSE
        drsqp1 = drsq0
        drsq0  = drsqm1
        drsqm1 = (rsw(iz0-3)-ir)**2+(zsw(iz0-3)-iz)**2
      END IF
    END IF
    
! check inner cap-discontinuity
    IF(tifzmn .AND. idir < 0 .AND. iz0 == izmin) THEN
      IF(drsqp1 < drsqm1) THEN  ! complement cap
        drsqm1 = drsq0
        drsq0  = drsqp1
        drsqp1 = (rsw(iz0-1)-ir)**2+(zsw(iz0-1)-iz0)**2
      ELSE
        drsqp1 = drsq0
        drsq0  = drsqm1
        drsqm1 = (rsw(iz0+3)-ir)**2+(zsw(iz0+3)-iz0)**2
      END IF
    END IF
! quadratic interpolation to minimum distance
    x      = drsqp1-drsqm1
    drsqmn = drsq0-x*x/(8.0*(drsqp1-drsq0-drsq0+drsqm1))
! save and privide proper sign
    delr(iz,ir) = SIGN(SQRT(ABS(drsqmn)),rhosw-rhoin(iz,ir))
    IF(itest >= 2) WRITE(7,'(7i4,2(1pg12.4))')  &
        iz,ir,mbase,idir,izini,izend,iz0,drsqmn,delr(iz,ir)
    mbase  = iz-iz0     ! save minimum position as base for next search
    mbase  = MAX(0,mbase)
  END DO
END DO

RETURN
END SUBROUTINE delrmn


!---------------------------------------------------------------------------  
! DESCRIPTION: damp_constraint
!> @brief
!!Computes \c distsf=distance to \c drho0l surface and applies subsequently 
!!the damping to the constraint. The \c distsf for a point in the region 
!!\c rho < \c drho0l is the minimum distance to a point in 
!!\c rho > \c drho0l, and vice versa.
!>
!> @param[in] iiter
!> INTEGER, takes the current iteration number.
!> @param[out] drho0l
!> REAL(8), returns the density at which the surface line is defined.
!--------------------------------------------------------------------------- 
SUBROUTINE damp_constraint(iiter,drho0l)

USE paramr
USE fields
USE grid, ONLY: nr,nz
USE dynpar
IMPLICIT NONE

INTEGER, INTENT(IN)          :: iiter
REAL(8), INTENT(OUT)         :: drho0l
LOGICAL :: lch

INTEGER :: iz,ir,iz2,ir2

REAL(8) :: drho0,dzmax,dzmin,dzmin0,z
  
! largest distance in whole grid
dzmin0 = dz*dz * REAL((nz+1)*(nz+1)) + dr*dr * REAL((nr+1)*(nr+1))

! maximum density
drho0 = MAX(MAXVAL(rho(0:nz,0:nr,1)+rho(0:nz,0:nr,2)),d0)
drho0l = drho0/10.0D0

! check whether surface has moved since last update
lch = .false.
DO ir = 0, nr
  DO iz = 0, nz
    IF( (rho(iz,ir,2)+rho(iz,ir,2) >= drho0l) .NEQV. lrho(iz,ir) ) THEN
      lrho(iz,ir) = .NOT. lrho(iz,ir)
      lch = .true.
    END IF
  END DO
END DO
IF (iiter > -1) lch = .true.

!  determine distance to surface for each point on the grid
IF( lch ) THEN
  DO ir = 0, nr
    DO iz = 0, nz
      dzmin = dzmin0
      dzmax = dzmin0
      DO ir2 = 0, nr
        DO iz2 = 0, nz
          z = dz*dz*(iz2-iz)*(iz2-iz) + dr*dr*(ir2-ir)*(ir2-ir)
          IF( lrho(iz2,ir2) ) THEN
            dzmin = MIN(dzmin, z)
          ELSE
            dzmax = MIN(dzmax, z)
          END IF
        END DO
      END DO
      IF( lrho(iz,ir) ) THEN
        distsf(iz,ir) = -SQRT(dzmax)
      ELSE
        distsf(iz,ir) = SQRT(dzmin)
      END IF
    END DO
  END DO
END IF

! apply damping 
distsf(0:nz,0:nr) = 1D0/(1D0+EXP((distsf(0:nz,0:nr)-dalp)/dgam))

IF(tswicon0) qfac0(0:nz,0:nr)=qcfac0(0:nz,0:nr)*distsf(0:nz,0:nr)
IF(tswicon2) qfac2(0:nz,0:nr)=qcfac2(0:nz,0:nr)*distsf(0:nz,0:nr)
IF(tswicon3) qfac3(0:nz,0:nr)=qcfac3(0:nz,0:nr)*distsf(0:nz,0:nr)
IF(tswicon4) qfac4(0:nz,0:nr)=qcfac4(0:nz,0:nr)*distsf(0:nz,0:nr)

RETURN
END SUBROUTINE damp_constraint

!---------------------------------------------------------------------------  
! DESCRIPTION: pronuc
!> @brief
!!Projection of high Fourier components on the wave function.
!> @details
!!WARNING: This routine might not be working well!
!>
!> @param[in] psiin
!> REAL(8), array, takes input wave function (destroyed after execution).
!> @param[in] psiout
!> REAL(8), array, returns the projection.
!> @param[in] ist
!> INTEGER, takes the idex of the input wave function.
!--------------------------------------------------------------------------- 
SUBROUTINE pronuc(psiin,psiout,ist)

IMPLICIT NONE

REAL(8), INTENT(IN OUT)         :: psiin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)            :: psiout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)             :: ist

INTEGER  :: iang,is,mm,nparr

WRITE(6,'(a)') ' warning: the routine PRONUC may be dangerous'
WRITE(7,'(a)') ' warning: the routine PRONUC may be dangerous'

DO  is=1,2
  mm     = mpsi(ist)+is-1
  nparr  = 1-2*MOD(mm,2)
  IF (is <= 2) THEN
    iang      = mpsi(ist)+is-1
  ELSE
    iang      = mpsi(ist)+is-3
  END IF
  IF (nparz(is,ist) == 1.AND.nparr == 1) THEN
    CALL mulr(rpro(0:nr,0:nr),psiin(0:nz,0:nr,is), psiout(0:nz,0:nr,is),one,0)

    psiin(0:nz,0:nr,is) = psiout(0:nz,0:nr,is)

    CALL mulz(zpro(0:nz,0:nz),psiin(0:nz,0:nr,is),psiout(0:nz,0:nr,is),one,0)
  ELSE IF (nparz(is,ist) == 1) THEN
    CALL mulz(zpro(0:nz,0:nz),psiin(0:nz,0:nr,is), psiout(0:nz,0:nr,is),one,0)
  ELSE IF (nparr == 1) THEN
    CALL mulr(rpro(0:nr,0:nr),psiin(0:nz,0:nr,is), psiout(0:nz,0:nr,is),one,0)
  END IF
  
  IF (nparz(is,ist) == 2) THEN
    psiout(0,0:nr,is)  =  zero
  END IF
  IF (iang > 0) THEN
    psiout(0:nz,0,is)  =  zero
  END IF
END DO

RETURN
END SUBROUTINE pronuc


!---------------------------------------------------------------------------  
! DESCRIPTION: procos
!> @brief
!!Construction of projection operator to project out higher Fourier
!!components. Level of projection determined by 'icut'.
!>
!> @param[out] p
!> REAL(8), array, returns the projection operator.
!> @param[out] n
!> INTEGER, takes the grid size.
!> @param[in] icut
!> INTEGER, takes the number components below the cutoff.
!--------------------------------------------------------------------------- 
SUBROUTINE procos(p,n,icut)

IMPLICIT NONE

REAL(8), INTENT(OUT)         :: p(0:n,0:n)
INTEGER, INTENT(IN)          :: n
INTEGER, INTENT(IN)          :: icut

REAL(8)                      :: dm,dn,dj,di,sum,afac,pfac
INTEGER                      :: i,j,m

dn    = (n)
afac  = pi/dn
pfac  = 2.0/dn
DO  j=0,n
  dj  = (j)
  DO  m=0,n
    dm     =  (m)
    sum    =  zero
    DO  i=n-icut+1,n
      di  = (i)
      IF (i == 0.OR.i == n) THEN
        sum = sum + 0.5*COS(afac*dj*di)*COS(afac*dm*di)
      ELSE
        sum = sum +     COS(afac*dj*di)*COS(afac*dm*di)
      END IF
    END DO
    IF(m == j) THEN
      p(j,j) = 1.0
    ELSE
      p(j,m) = 0.0
    END IF
    IF (m == 0.OR.m == n) THEN
      p(j,m) = p(j,m)-0.5*pfac * sum
    ELSE
      p(j,m) = p(j,m)    -pfac * sum
    END IF
  END DO
END DO
WRITE(*,*)'*****************p:'
DO i=0,n
  WRITE(*,'(999f4.1)')p(:,i)
END DO
RETURN
END SUBROUTINE procos
