!------------------------------------------------------------------------------
! MODULE: dynpar
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the variables for the force definitions, termination 
!!criteria, parameters for iterations and dynamic variables for pairing. The 
!!routine \c read_force reads forces from the input file or from forces.data.
!------------------------------------------------------------------------------
MODULE dynpar

USE paramr, ONLY : xxxl,ixxxl,oh,one
!USE FyFunc

IMPLICIT NONE

! numerical parameters for iteration etc
REAL(8) :: x0inv=0.3D0     !<Overall factor in front of the damping operator.
REAL(8) :: e0inv=120D0     !<Witdh of the damping operator.
REAL(8) :: addrho=0.2D0    !<Fraction of new density to be added to the old, used to compute 'adyrho'.
REAL(8) :: addpot=1D0      !<Fraction of new potential to be added to the old.
REAL(8) :: addocc=0.3D0    !<Fraction of new occupation numbers 'VSQF' to be added to the old.
REAL(8) :: errenr=1D-9     !<Termination criterion for energy change.
REAL(8) :: errvar=1D-1     !<Termination criterion for energy variance.
REAL(8) :: adyrho          !<Dynamical fraction of new density.
REAL(8) :: adypot          !<Dynamical  fraction of new potential.
INTEGER :: mpmax           !<Max no. of z component of angular momentum.
INTEGER :: itert           !<Overall total number of iterations.
INTEGER :: iterbcs=300     !<iteration nr. from which on HFB is active
INTEGER :: iterhfb=600     !<nr. of iterations with enhanced pairing
INTEGER :: iternat=0       !<nr. of iterations in configuration space (not yet installed)
INTEGER :: modnat=9999999  !<modulus of switxhing to iterations in configuration space
REAL(8) :: pairenhance=2D0 !<initial enhancement of pairing strength
LOGICAL :: tnobcs=.TRUE.   !<suppresses initial BCS phase after first run along a deformation path
CHARACTER (LEN=3) :: outnam = "tst"  !< 3-digit qualifier for output files
CHARACTER (LEN=80) :: title=" "  !< Title for the actual run, printed in the header of output files.


! print switches
INTEGER :: ipri_fields=0  !< switches printing of fields on 'fields.*'
LOGICAL :: tpri_occup=.FALSE.   !< switches printing of l-sorted occupation numbers on 'occup_lz_*.*'
LOGICAL :: tpri_conver=.FALSE.  !< switches printing of convergence protocol on 'conver.*'
LOGICAL :: tpri_bubble=.FALSE.  !< switches printing of central depresion on 'bubble.*'
LOGICAL :: tpri_2d=.FALSE.      !< switches printing of 2D plots in 'plot\_2D.***' for converged calculations.
LOGICAL :: tpri_radii=.TRUE.      !< switches printing of charge radii and surface thickness
LOGICAL :: tpri_formf=.FALSE.      !< switches printing of formfactors and charge distribution



CHARACTER(8) :: name=''   
!>@name Skyrme parameters. 
!>@{
REAL(8),PUBLIC :: t0=xxxl,t1=xxxl,t2=xxxl,t3=xxxl,t4=xxxl
REAL(8),PUBLIC :: power=xxxl                !<Power of density dependence in 't3'.
!>@}
!>@name Skyrme exchange parameters.
!>@{
REAL(8),PUBLIC :: x0=xxxl,x1=xxxl,x2=xxxl,x3=xxxl,b4p=xxxl
!>@}
!>@name Derived Skyrme parameters.
!>@{
REAL(8),PUBLIC :: b0,b0p,b3,b3p,t3a,t3b,b1,b1p,b2,b2p,b4,excfac,&
                  t0a,t0b,t3as,t3bs,t1p2,t1m2,t13p2,t13m2,t1m20,t1p2x,&
                  t7b,t7a,b7a,b7pa,b7pb
!>@}

!>@name Fayans parameters, obligatory input preset to 'xxxl'.
!>@{
LOGICAL :: tpairFy   !???
LOGICAL :: tFayans=.FALSE. !< switch to Fayans functional
REAL(8) :: avp=xxxl    !< Fayans parameter
REAL(8) :: h1vp=xxxl   !< Fayans parameter
REAL(8) :: h2vp=xxxl   !< Fayans parameter
REAL(8) :: avm=xxxl   !< Fayans parameter
REAL(8) :: h1vm=xxxl   !< Fayans parameter
REAL(8) :: h2vm=xxxl   !< Fayans parameter
REAL(8) :: asp=xxxl   !< Fayans parameter
REAL(8) :: hsp=xxxl   !< Fayans parameter
REAL(8) :: hgradp=xxxl   !< Fayans parameter
REAL(8) :: asm=0D0   !< Fayans parameter
REAL(8) :: hsm=-2D3   !< Fayans parameter
REAL(8) :: hgradm=0D0   !< Fayans parameter
REAL(8) :: fxi=xxxl   !< Fayans parameter
REAL(8) :: h1xi=xxxl   !< Fayans parameter
REAL(8) :: hgradxi=xxxl   !< Fayans parameter
REAL(8) :: delfxi=xxxl    !< Fayans parameter
REAL(8) :: c0tau=0D0   !< Fayans parameter
REAL(8) :: c1tau=0D0   !< Fayans parameter
REAL(8) :: c0nabj=xxxl   !< Fayans parameter
REAL(8) :: c1nabj=xxxl   !< Fayans parameter
REAL(8) :: c0tens=0D0   !< Fayans parameter
REAL(8) :: c1tens=0D0   !< Fayans parameter
REAL(8) :: gammaxi=0.66666666667   !< Fayans parameter
REAL(8) :: hcoul=0D0   !< Fayans parameter
INTEGER :: isovectormodel=ixxxl !< Fayansparameter
REAL(8),PARAMETER :: rho0fy=0.08D0   !< Fayans parameter
REAL(8),PARAMETER :: sigmavol=0.333333333333D0   !< frozen Fayans parameter
REAL(8),PARAMETER :: rho0xi=0.16D0   !< frozen Fayans parameter
REAL(8),PARAMETER :: sigmacoul=0.333333333333D0   !< frozen Fayans parameter
!>@}


REAL(8),PUBLIC :: h2m(2)=(/xxxl,xxxl/)      !<\f$ \hbar^2/m_B \f$.
REAL(8),PUBLIC :: h2meff(2)                 !<Effective mass (for 1-1/A correction).
REAL(8),PUBLIC :: alpinv=0D0                !<Density dependent force parameters.
REAL(8),PUBLIC :: betinv=0D0                !<Density dependent force parameters.
!>@name extended Skyrme  parameters
!>@{
REAL(8),PUBLIC :: b7=0D0,b7p=0D0,pow7=1D0
!>@}
REAL(8),PUBLIC :: chg=1D0                   !<Effective charge in units \f$ e^2 \f$ , usually 1.0.
REAL(8),PUBLIC :: chgx=1D0                  !<Effective charge for exchange term in units \f$ e^2 \f$ , usually 1.0.
REAL(8),PUBLIC :: riftm=0D0                 !<Parameter for contribution of t1 and t2 term for spin-orbit.
REAL(8),PRIVATE:: v0neut=xxxl               !<The strength of pairing for protons in MeV. (Only for input)
REAL(8),PRIVATE:: v0prot=xxxl               !<The strength of pairing for neutrons in MeV. (Only for input)
INTEGER,PRIVATE:: ex=ixxxl                  !<Some forces are fitted excluding the Coulombexchange term. (Only for input)
INTEGER,PRIVATE:: zpe=ixxxl                 !<If 0, m-m/A recipe for c.m. correction, if 1 no c.m. correction. (Only for input)

REAL(8) :: h2mav                            !<Average mass paramter
REAL(8) :: riftm1                           !<Parameter for contribution of t1 and t2 term for spin-orbit.
INTEGER :: ifls                             !<Switch for contribution from t1,t2 to spin-orbit force.
INTEGER :: ifzpe                            !<If 0, m-m/A recipe for c.m. correction, if 1 no c.m. correction.
INTEGER :: ifex                             !<Some forces are fitted excluding the Coulombexchange term. 
!!For <tt> ex=1 </tt> it is included (this is the normal case), for <tt> ex=0 </tt> not.


LOGICAL :: tpair                            !<Switch for density-dependent pairing.

REAL(8) :: pairfc(2)=(/xxxl,xxxl/)          !<The strength of pairing for protons and neutrons in MeV.
REAL(8) :: rho0pr=xxxl                      !<Density parameter for the density-dependent delta pairing.
REAL(8) :: pairstab=0D0                     !<Parameter for pairing stabilization
REAL(8) :: pair_cutoff=1.65D0               !<Parameter for soft cutoff of pairing space (<0 = WS profile, >0 = nr. of states)
REAL(8) :: erra=1D-9                        !<Termination criterion for error in particle number.
LOGICAL :: tcutoff_uvspace=.TRUE.           !<Switch to pairing cutoff in u*v weight, else in v*v weight.
LOGICAL :: tcutoff_WS                       !<Switch to pairing cutoff by Woods-Saxon profile
LOGICAL, PARAMETER :: thfbvariat=.FALSE.    !<Switch to variational cutoff -- placeholder, not yet fully coded

INTEGER :: ipair=ixxxl                      !<Defines the pairing model.
INTEGER :: npiter=100                       !<Maximum number of iterations for the solution of the gap equation.
INTEGER :: itpini=5                         !<Number of iterations before pairing is switched on.
INTEGER :: itpfix=0                         !<Number of iterations after which pairing occupations are kept fixed (set 0 to disable).

!REAL(8) :: ggg(2)                           !<Pairing force constant.
REAL(8) :: gap(2)                           !<Pairing gap.
REAL(8) :: fer(2)                           !<Fermi energy.
REAL(8) :: epair(2)                         !<Pairing energy.
REAL(8) :: fer2(2)                          !<Phase space fermi energy.
REAL(8) :: feradd(2)                        !<Difference between phase space Fermi energy and Fermi energy.

REAL(8) :: stwid=10D0                       !<Width for stabilized boundary conditions for constraint calculations.
REAL(8) :: stslop=15D0                      !<Slope for stabilized boundary conditions for constraint calculations.
REAL(8) :: rhocut=0.03D0                    !<Density from which stabilizing potential is added.
LOGICAL :: tstab=.FALSE.                    !<Switch for stabilizing boundaries.
LOGICAL :: thfb=.FALSE.                     !<Indicates thet HFB is active
LOGICAL :: tregpair=.FALSE.                 !<Switches regularized pairing

  !> This contains the parameters for pairing:
  TYPE Pairing
     REAL(8) :: v0prot!<Strength of pairing for protons in MeV.
     REAL(8) :: v0neut!<Strength of pairing for neutrons in MeV.
     REAL(8) :: rho0pr!<Density parameter for the density-dependent delta pairing.
  END TYPE Pairing
  
  !> This contains the parameters for the Skyrme force:
  TYPE Force
     CHARACTER(8) :: name         !<Name of the force used to identify it.
     INTEGER :: ex                !<Some forces are fitted excluding the Coulomb
     !!exchange term. For <tt> ex=1 </tt> it is included (this is the normal
     !!case), for <tt> ex=0 </tt> not.
     INTEGER :: zpe               !<Index for the treatment of the center-of-mass correction.
     REAL(8) :: h2m(2)            !<Value of \f$ \frac{\hbar^2}{2m} \f$, separately for neutrons and protons.
     !>@name Skyrme parameters. 
     !>@{
     REAL(8) :: t0,t1,t2,t3,t4
     !>@}
     !>@name Skyrme exchange parameters.
     !>@{
     REAL(8) :: x0,x1,x2,x3,b4p
     !>@}
     REAL(8) :: power             !<Exponent in the nonlinear (originally three-body) term.
     TYPE(Pairing) :: vdi         !<Parameter set for the volume-delta pairing case. 
     TYPE(Pairing) :: dddi        !<Parameter set for the density-dependent delta pairing case.
  END TYPE Force
  !> This contains the parameters for the Skyrme force:
  TYPE ForceFy
     CHARACTER(8) :: name         !<Name of the force used to identify it.
     INTEGER :: ex                !<Some forces are fitted excluding the Coulomb
     !!exchange term. For <tt> ex=1 </tt> it is included (this is the normal
     !!case), for <tt> ex=0 </tt> not.
     INTEGER :: zpe               !<Index for the treatment of the center-of-mass correction.
     REAL(8) :: h2m(2)            !<Value of \f$ \frac{\hbar^2}{2m} \f$, separately for neutrons and protons.
     !>@name frequently used Fayans parameters 
     !>@{
     REAL(8) :: avp,h1vp,h2vp,avm,h1vm,h2vm,asp,hgradp,c0nabj,c1nabj,fxi,h1xi,hgradxi,delfxi,pair_cutoff,pairstab
     !>@}
     !>@name rarely used Fayans parameters 
     !>@{
     REAL(8) :: hsp,asm,hsm,hgradm,c0tau,c1tau,hcoul
     !>@}
     !>@name set the isovector model
     !>@{
     INTEGER :: isovectormodel
     !>@}
     !>@name determine type of cutoff space (uv weight or v**2 weight)
     !>@{
     LOGICAL :: tcutoff_uvspace
     !>@}
  END TYPE ForceFy
  ! include predefined forces
  INCLUDE 'forces.data'
  !
  TYPE(Force) :: forcact         !<This contains parameters for
  !!the Skyrme force actually used in the present calculation, packed
  !into the derived-type \c Force.
  TYPE(Pairing) :: pairact       !<Pairing parameters used in the
  !!present Skyrme calculation. This is separate from the force itself: the
  !!force definition usually contains suggestions for the associated
  !!pairing, but this often overridden, e.g, by turning off
  !!pairing.
  TYPE(ForceFy) :: forcactFy     !<This contains parameters for
  !!the Fayans force actually used in the present calculation, packed
  !into the derived-type \c Force.


CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: read_force
!> @brief
!!The purpose of the subroutine is to read the force and pairing
!!definitions. The NAMELIST \c force contains all the
!!defining values for a Skyrme force (but now as individual variables,
!!not in a derived type) plus a selection of pairing type and strengths
!!(see input description). In addition there is a logical variable
!!\c turnoff_zpe which allows turning off the center-of-mass correction.
!>
!> @details
!!Some quantities are first set large to see whether required input
!!is missing. Then a predefined Skyrme force with the given name is
!!sought; if it is found, it is simply copied into \c f. If no force
!!of this name is found, a new one is composed from the numbers given in
!!the input.
!!
!!If the input varable \c turnoff_zpe is true, the indicator <tt>
!!f\%zpe </tt> is set to 1, which implies not doing anything about the
!!center-of-mass correction. 
!!
!!Now the "b" coefficients are calculated straightforwardly. 
!!
!!The variable \c pairing from the namelist then determines the
!!pairing. If it is set to \c 'NONE', no pairing is
!!included. Otherwise the strength parameters are taken from the input
!!or from the predefined force. If this process does not find a
!!reasonable pairing combination, stop with an error message.
!--------------------------------------------------------------------------- 
  SUBROUTINE read_force
  
    IMPLICIT NONE
    
    REAL(8) :: dummy
    INTEGER :: i
    LOGICAL :: predefined

    NAMELIST /force/ name, &
         chg,ex,zpe,h2m,t0,t1,t2,t3,t4,x0,x1,x2,x3,b4p,power, &
         ipair,v0prot,v0neut,rho0pr,pairstab,pair_cutoff,riftm, &
         tFayans, &
         avp,h1vp,h2vp,avm,h1vm,h2vm,asp,hsp,hgradp,asm,hsm,hgradm, &
         fxi,h1xi,hgradxi,delfxi,c0tau,c1tau,c0nabj,c1nabj,hcoul, &
         isovectormodel,tcutoff_uvspace,tregpair

    READ(5,force)
    predefined=.FALSE.
    DO i=1,nforce
       IF(TRIM(name)==TRIM(pforces(i)%name)) THEN
          predefined=.TRUE.
          forcact=pforces(i)
          tFayans=.FALSE.
       ENDIF
    END DO
    DO i=1,nforceFy
       IF(TRIM(name)==TRIM(pforcesFy(i)%name)) THEN
          predefined=.TRUE.
          forcactFy=pforcesFy(i)
          tFayans=.TRUE.
       ENDIF
    END DO

    IF(tFayans) THEN
      CALL prepare_Fayans()
    ELSE
      CALL prepare_Skyrme()
    END IF

    IF(pair_cutoff < 0D0) THEN
      tcutoff_WS = .TRUE.
      pair_cutoff=ABS(pair_cutoff)
    ELSE
      tcutoff_WS = .FALSE.
    END IF

CONTAINS

  SUBROUTINE prepare_Skyrme()
    ! check whether model parameters are set explicitly
    IF(.NOT.predefined) THEN
      IF(t0==xxxl .AND. t1==xxxl .AND. t2==xxxl .AND. t3==xxxl .AND.  &
           x0==xxxl .AND. x1==xxxl .AND. x2==xxxl .AND. x3==xxxl .AND. & 
           t4==xxxl .AND. b4p==xxxl .AND. &
           power==xxxl .AND. h2m(1)==xxxl .AND. h2m(2)==xxxl .AND. &
           ex==ixxxl .AND.  zpe==ixxxl) THEN
         name='SV-min'
         forcact=pforces(1)
         predefined=.TRUE.
      ELSE IF(t0==xxxl .OR. t1==xxxl .OR. t2==xxxl .OR. t3==xxxl .OR.  &
           x0==xxxl .OR. x1==xxxl .OR. x2==xxxl .OR. x3==xxxl .OR. & 
           t4==xxxl .OR. b4p==xxxl .OR. &
           power==xxxl .OR. h2m(1)==xxxl .OR. h2m(2)==xxxl .OR. &
           ex==ixxxl .OR.  zpe==ixxxl) THEN
         STOP 'definition of Skyrme force incomplete'
      ELSE
        WRITE(*,*) 'Using new Skyrme force defined explicitly'
      END IF
    END IF
    IF(predefined) THEN
       WRITE(*,*) 'Using predefined force ',name
!       forcactname=name
       ex=forcact%ex; zpe=forcact%zpe
       h2m=forcact%h2m
       t0=forcact%t0; t1=forcact%t1; t2=forcact%t2; t3=forcact%t3; t4=forcact%t4
       x0=forcact%x0; x1=forcact%x1; x2=forcact%x2; x3=forcact%x3; b4p=forcact%b4p;
       power=forcact%power; ex=forcact%ex; zpe=forcact%zpe
       IF(ABS(b4p)<1.0d-9) b4p=forcact%t4/2.0D0
       dummy=h2m(1); h2m(1)=h2m(2); h2m(2)=dummy  ! reorder proton<-->neutron
    ENDIF
    ifzpe=zpe
    ifex = ex
    !     
    ! now set up pairing 
    IF(predefined .AND. (v0prot==xxxl .OR. v0neut==xxxl)) THEN
       IF(ipair==5 .OR. ipair==9) THEN 
         pairact=forcact%vdi
         v0prot=pairact%v0prot; v0neut=pairact%v0neut
       ELSE IF(ipair==7 .OR. ipair==8 .OR. ipair==11) THEN
         pairact=forcact%dddi
         v0prot=pairact%v0prot; v0neut=pairact%v0neut; rho0pr=pairact%rho0pr
       END IF
       IF(v0prot<0D0 .OR. v0neut<0D0) STOP 'Pairing not defined properly'
    ENDIF
    ! stop if this has not delivered meaningful parameters
    IF(v0prot==xxxl .OR. v0neut==xxxl) STOP 'Pairing not defined properly'
    IF((ipair==7 .OR. ipair==11) .AND. rho0pr==xxxl) &
        STOP 'DDDI-Pairing: undefined RHO0PR'
    pairfc(1)=v0prot; pairfc(2)=v0neut
    
    b0     = t0*(one+oh*x0)
    b0p    = t0*(oh+x0)
    b3     = t3*(one+oh*x3)/4.0
    b3p    = t3*(oh+x3)/4.0
    t3a    = t3*(one+oh*x3)/6.0
    t3b    = t3*(oh+x3)/6.0
    b1     = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
    b1p    = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
    b2     = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
    b2p    = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
    b4     = oh*t4
    excfac = -1.418033*chgx               ! coulomb exchange in Slater appr.
    
    t0a    = t0*(one+oh*x0)
    t0b    = t0*(oh+x0)
    t3as   = oh*power*t3a
    t3bs   = oh*power*t3b
    t1p2   = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
    t1m2   = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
    t13p2  = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
    t13m2  = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)

    IF(riftm /= 0D0) THEN
      t1m20  = (t1-t2)/8.0 * riftm1
      t1p2x  = (x1*t1+x2*t2)/8.0 * riftm
    ELSE
      t1m20  = 0.0D0
      t1p2x  = 0.0D0
    END IF
    
    t7b  = b7p/6D0
    t7a  = b7/6D0+0.5D0*t7b
    b7a  = t7a*(pow7+2.0D0)*0.5D0
    b7pa = t7b
    b7pb = t7b*pow7*0.5D0

  END SUBROUTINE prepare_Skyrme

  SUBROUTINE prepare_Fayans()

    IF(.NOT.predefined) THEN
      IF(ex==ixxxl .OR. zpe==ixxxl .OR. isovectormodel==ixxxl .OR. &
       h2m(1)==xxxl .OR. h2m(2)==xxxl  .OR. &
       avp==xxxl .OR. h1vp==xxxl .OR. h2vp==xxxl .OR. &
       avm==xxxl .OR. h1vm==xxxl .OR. h2vm==xxxl .OR. asp==xxxl .OR. &
       hsp==xxxl .OR. hgradp==xxxl .OR. &
       fxi==xxxl .OR. h1xi==xxxl .OR. hgradxi==xxxl .OR. delfxi==xxxl .OR. &
       c0nabj==xxxl .OR. c1nabj==xxxl .OR. isovectormodel==ixxxl) THEN
!       asm==xxxl .OR. &
!       hsm==xxxl .OR. &
!       hgradm==xxxl .OR. &
!       c0tau==xxxl .OR. &
!       c1tau==xxxl .OR. &
!       hcoul==xxxl) THEN
         STOP 'definition of Fayans force incomplete'
      ELSE
        WRITE(*,*) 'Using new Fayans force defined explicitly'
      END IF
    END IF
    IF(predefined) THEN
       WRITE(*,*) 'Using predefined force ',name
!       forcactname=name
       ex=forcactFy%ex; zpe=forcactFy%zpe
       h2m=forcactFy%h2m
       avp=forcactFy%avp
       h1vp=forcactFy%h1vp
       h2vp=forcactFy%h2vp
       avm=forcactFy%avm
       h1vm=forcactFy%h1vm
       h2vm=forcactFy%h2vm
       asp=forcactFy%asp
       hsp=forcactFy%hsp
       hgradp=forcactFy%hgradp
       asm=forcactFy%asm
       hsm=forcactFy%hsm
       hgradm=forcactFy%hgradm
       fxi=forcactFy%fxi
       h1xi=forcactFy%h1xi
       hgradxi=forcactFy%hgradxi
       delfxi=forcactFy%delfxi
       pair_cutoff=forcactFy%pair_cutoff
       pairstab=forcactFy%pairstab
       c0tau=forcactFy%c0tau
       c1tau=forcactFy%c1tau
       c0nabj=forcactFy%c0nabj
       c1nabj=forcactFy%c1nabj
       hcoul=forcactFy%hcoul
       isovectormodel=forcactFy%isovectormodel
       tcutoff_uvspace=forcactFy%tcutoff_uvspace
       dummy=h2m(1); h2m(1)=h2m(2); h2m(2)=dummy  ! reorder proton<-->neutron
    ENDIF
    ifzpe=zpe
    ifex = ex
     WRITE(7,force)
!    CALL ini_fyfunc()

  END SUBROUTINE prepare_Fayans


  END SUBROUTINE read_force


END MODULE dynpar
