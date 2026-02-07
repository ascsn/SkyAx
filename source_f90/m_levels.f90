!------------------------------------------------------------------------------
! MODULE: levels
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains all variables defining the single-particle levels.
!------------------------------------------------------------------------------
MODULE levels

USE paramr, ONLY : mxd,nsd

IMPLICIT NONE

REAL(8),ALLOCATABLE :: wate(:)    !<Multiplicity due to spin degeneracy, in the problem 
                                  !!here always 2.0 (kramers degeneracy).
REAL(8),ALLOCATABLE :: spn(:)     !<Norm of the s.p. wave function, used as a countercheck 
                                  !!(it should be 1.0)
REAL(8),ALLOCATABLE :: spenrg(:)  !<Single particle energy.
REAL(8),ALLOCATABLE :: spevar(:)  !<Variance of the single particle energy.
REAL(8),ALLOCATABLE :: spevhh(:)  !<storage for h*h variance
REAL(8),ALLOCATABLE :: speold(:)  !<Single particle energy from previous iteration.
REAL(8),ALLOCATABLE :: vsqf(:)    !<Occupation probability. In case of pure slater determinants 
                                  !!1.0 throughout. in case of frozen occupation the valence 
                                  !!states are occupied with 0.0<vsqf<1.0 . in case of pairing 
                                  !!the \c vsqf is computed dynamically according to the actual 
                                  !!\c spenrg in the subroutine \c vsqtr. 
REAL(8),ALLOCATABLE :: dvsqf(:)   !<Change in occupation \f$ v^2 \f$.
REAL(8),ALLOCATABLE :: vsqfol(:)  !<Storage for old occupation probability from the previous iteration.
REAL(8),ALLOCATABLE :: oc(:)      !<Intermediate storage of occupation probabilities.
REAL(8),ALLOCATABLE :: cmpn(:,:)  !<Norm of the single components (\c is=1,.., \c nsd) separately.
REAL(8),ALLOCATABLE :: vsqobs(:)  !<Occupation probablities for observables, from Lipkin-Nogami.
REAL(8),ALLOCATABLE :: deltaf(:)  !<State dependent gap for dynamic BCS.
!REAL(8),ALLOCATABLE :: wtpair(:)  !<Dynamic weight for BCS phase space.
REAL(8),ALLOCATABLE :: wtpair_uv(:)  !<Dynamic weight for BCS phase space, u*v wieght
REAL(8),ALLOCATABLE :: wtpair_v2(:)  !<Dynamic weight for BCS phase space, v**2 weight
REAL(8),ALLOCATABLE :: qspav(:)   !<Diagonal element of constraint.
INTEGER,ALLOCATABLE :: mpsi(:)    !<k-0.5, where k is the z-component of total angular momentum. 
                                  !!Note that mpsi is also the z-component of the orbital a.m. 
                                  !!for the upper component with spin up.
INTEGER,ALLOCATABLE :: npar(:)    !<z-parity for the upper component spin up (component \c is=1). 
                                  !!the subsequent components have then z-parities
                                  !!  - -npar for spin down (is=2)
INTEGER,ALLOCATABLE:: nrpar(:)    !<r-parity for the upper component spin up (component \c is=1). The 
                                  !!subsequent components have then r-parities
                                  !!  - -nrpar for spin down (is=2)
INTEGER,ALLOCATABLE :: nparz(:,:) !< z-parity for upper and lower spin component.
INTEGER,ALLOCATABLE :: nqr(:)     !< Radial quantum numbers.
INTEGER,ALLOCATABLE :: nqz(:)     !< Longitudinal quantum numbers.
INTEGER ,ALLOCATABLE:: icharg(:)  !<Isospin (1 for protons, 2 for neutrons).
INTEGER,ALLOCATABLE :: iorder(:)  !<Pointer field for reordering according to the actual energy.

!NEW ODD
REAL(8),ALLOCATABLE :: rtpartn(:) !<Overlap os state and conjugate partner state (auxiliary)
INTEGER,ALLOCATABLE :: itpartn(:) !<Pointer field for conjugate partner states
INTEGER :: kblock(2)=(/0,0/)      !<2*K values of blocked s.p. state. (Default 0)
INTEGER :: ktcblock(2)=(/0,0/)    !<2*K values of conjugate to blocked s.p. state. (Default 0)
INTEGER :: mblock(2)=(/0,0/)      !<M values of blocked s.p. state. 
INTEGER :: istblo(2)=(/0,0/)      !<Number of blocked s.p. state in level sequence.
INTEGER :: istcblo(2)=(/0,0/)     !<Number of conjugate to blocked s.p. state in level sequence.
INTEGER :: isoblo=0               !<Isospin of blocked s.p. state. To be set at compile time.
!END NEW ODD
INTEGER :: ifocsv           !<Switch to store occupation numbers intermediately (field 'oc').

INTEGER :: nprot=0          !<Proton number.
INTEGER :: nneut=0          !<Neutron number.

REAL(8) :: anumb            !<Total particle number.
INTEGER :: nst(2)           !<Total no. of states.
INTEGER :: npsi(2)          !<Upper index of last wf with isospin i.
INTEGER :: npart(2)         !<Particle numbers.
INTEGER :: iokk(2)          !<Pointer to the 1 st index of state with \c iso.
INTEGER :: numsta(2)        !<Number of available states for particle sort.
INTEGER :: ls              !<Total no. of states.
 
INTEGER :: nom              !<Number of m blocks (maximum m quantum number +1).
INTEGER :: mhi(mxd)         !<array storing the upper index for the wave functions 
INTEGER :: mlo(mxd)         !<array storing the lower index for the wave functions
INTEGER :: miq(mxd)         !<array storing the iq value for the wave functions
INTEGER :: mim(mxd)         !<array storing the m value for the wave functions
INTEGER,ALLOCATABLE :: nbuf(:)        !<array storing the m value for the wave functions
 
REAL(8) :: betain=0.1D0               !< Deformation of initial oscillator wavefunctions.
REAL(8) :: octini=0D0                 !< Initial octupole deformation (by oczupole squeezing).
REAL(8) :: deocc(2)=(/1.9D0,1.9D0/)   !< Number of additional oscillator shells (above Fermi state) for initialization.
INTEGER :: mmax(2)=(/29,29/)          !< Maximum m quantum number of protons and neutrons.


REAL(8) :: eocc(2)                   !< Number of oscillator shells for in initialization, derived from basic number (according to N or Z) plus 'deocc'.

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_levels
!> @brief
!!Allocates all arrays describing the energy levels.
!>
!> @param[in] npsi
!> INTEGER, takes the total number of wave function.
!--------------------------------------------------------------------------- 
SUBROUTINE alloc_levels(npsi)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: npsi
  ALLOCATE(&
           wate(npsi),  spn(npsi),   spenrg(npsi),spevar(npsi),  speold(npsi),vsqf(npsi),&
           dvsqf(npsi), vsqfol(npsi),oc(npsi),    cmpn(nsd,npsi),vsqobs(npsi),deltaf(npsi),&
           wtpair_uv(npsi),wtpair_v2(npsi),qspav(npsi),mpsi(npsi),npar(npsi),nrpar(npsi),&
            nparz(2,npsi),nqr(npsi),   nqz(npsi),   icharg(npsi),iorder(npsi),  spevhh(npsi), &
           nbuf(npsi)       )
  wate=0.0D0  ;spn=0.0D0   ;spenrg=0.0D0;spevar=0.0D0;speold=0.0D0;vsqf=0.0D0
  dvsqf=0.0D0 ;vsqfol=0.0D0;oc=0.0D0    ;cmpn=0.0D0  ;vsqobs=0.0D0;deltaf=0.0D0
  wtpair_uv=0.0D0;wtpair_v2=0.0D0;qspav=0.0D0 ;mpsi=0  ;npar=0  ;nrpar=0 ;nparz=0
  nqr=0   ;nqz=0   ;icharg=0;iorder=0;spevhh=0.0D0; nbuf=0
!NEW ODD
  ALLOCATE(itpartn(npsi),rtpartn(npsi))
  itpartn=0;rtpartn=0D0
!END NEW ODD
  RETURN
END SUBROUTINE alloc_levels
!---------------------------------------------------------------------------  
! DESCRIPTION: dealloc_levels
!> @brief
!!Deallocates all arrays describing the energy levels.
!--------------------------------------------------------------------------- 
SUBROUTINE dealloc_levels
  IMPLICIT NONE
  DEALLOCATE(&
           wate,  spn,   spenrg,spevar,  speold,vsqf,&
           dvsqf, vsqfol,oc,    cmpn,    vsqobs,deltaf,&
           wtpair_uv,wtpair_v2,qspav, mpsi,  npar,    nrpar, nparz,&
           nqr,  nqz,   icharg, iorder,  spevhh, nbuf &
           )
  DEALLOCATE(itpartn,rtpartn)    !NEW ODD
  RETURN
END SUBROUTINE dealloc_levels



END MODULE levels
