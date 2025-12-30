!------------------------------------------------------------------------------
! MODULE: fields
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!Contains variables for densities and potentials.
!------------------------------------------------------------------------------
MODULE fields

USE paramr, ONLY : kqm
USE grid,   ONLY : nz,nr
USE dynpar, ONLY : tregpair

IMPLICIT NONE

REAL(8),ALLOCATABLE :: pot(:,:,:)      !<Local HF potential.
REAL(8),ALLOCATABLE :: coulpt(:,:)     !<Coulomb potential.
REAL(8),ALLOCATABLE :: bmass(:,:,:)    !<Inverse mass.
REAL(8),ALLOCATABLE :: potls(:,:,:,:)  !<Spin-orbit potential.
REAL(8),ALLOCATABLE :: drbmas(:,:,:)   !<d/dr \c bmass.
REAL(8),ALLOCATABLE :: dzbmas(:,:,:)   !<d/dz \c bmass.
REAL(8),ALLOCATABLE :: potadd(:,:,:)   !<Additional potential to stabilize positive energy states.
REAL(8),ALLOCATABLE :: rhosur(:,:,:)   !<Auxiliary field for potential wall.
REAL(8),ALLOCATABLE :: pairpt(:,:,:)   !<Gap potential for dynamic BCS.
REAL(8),ALLOCATABLE :: potdamp(:,:)    !<Damping potential for potential wall.
!NEW ODD
REAL(8),ALLOCATABLE :: potS(:,:,:,:)   !<(time-odd) spin potential
REAL(8),ALLOCATABLE :: potC(:,:,:,:)   !<(time odd) spin-orbit potential.
REAL(8),ALLOCATABLE :: potA(:,:,:,:)   !<(time odd) A potential
!END NEW ODD

REAL(8),ALLOCATABLE :: rho(:,:,:)      !<Density.
REAL(8),ALLOCATABLE :: crho(:,:)       !<Coulomb density.
REAL(8),ALLOCATABLE :: tau(:,:,:)      !<Kinetic energy density.
REAL(8),ALLOCATABLE :: drhols(:,:,:)   !<Gradient of spin-orbit current.
REAL(8),ALLOCATABLE :: rholsr(:,:,:)   !<Spin-orbit current in r-direction.
REAL(8),ALLOCATABLE :: rholsz(:,:,:)   !<Spin-orbit current in z-direction.
REAL(8),ALLOCATABLE :: d2rho(:,:,:)    !<Laplaceian on \c rho.
REAL(8),ALLOCATABLE :: rhoti(:,:,:)    !<Density for constraint calculation.
REAL(8),ALLOCATABLE :: paircr(:,:,:)   !<Pairing density.
REAL(8),ALLOCATABLE :: regpairadd(:,:,:)!<Rearranegment potential in regularized pairing
!NEW ODD
REAL(8),ALLOCATABLE :: rotjr(:,:,:)         !< r-comp rot j
REAL(8),ALLOCATABLE :: rotjz(:,:,:)         !< z-comp rot j
REAL(8),ALLOCATABLE :: rotspden(:,:,:)	!< rot spin-density
REAL(8),ALLOCATABLE :: spdenr(:,:,:)	!< r-comp. spin-desity
REAL(8),ALLOCATABLE :: spdenz(:,:,:)	!< z-comp. spin-density
REAL(8),ALLOCATABLE :: d2spdenr(:,:,:)	!< r-comp. Laplace spin-desity
REAL(8),ALLOCATABLE :: d2spdenz(:,:,:)	!< z-comp. Laplace spin-density
REAL(8),ALLOCATABLE :: curphi(:,:,:)	!< phi-comp. current
REAL(8),ALLOCATABLE :: spdenkr(:,:,:)       !<r-comp. kinetic spin-den
REAL(8),ALLOCATABLE :: spdenkz(:,:,:)       !<z-comp. kinetic spin-den
!END NEW ODD

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: alloc_fields
!> @brief
!!Allocates all fields needed in static HF+BCS calculations.
!--------------------------------------------------------------------------- 
SUBROUTINE alloc_fields()
ALLOCATE(&
         pot   (0:nz,0:nr,kqm),coulpt (0:nz,0:nr)    ,bmass (0:nz,0:nr,kqm),potls (0:nz,0:nr,kqm,2),&
         drbmas(0:nz,0:nr,kqm),dzbmas (0:nz,0:nr,kqm),potadd(0:nz,0:nr,kqm),rhosur(0:nz,0:nr,kqm)  ,&
         pairpt(0:nz,0:nr,kqm),potdamp(0:nz,0:nr)    ,rho   (0:nz,0:nr,kqm),crho  (0:nz,0:nr)      ,&
         tau   (0:nz,0:nr,kqm),drhols (0:nz,0:nr,kqm),rholsr(0:nz,0:nr,kqm),rholsz(0:nz,0:nr,kqm)  ,&
         d2rho (0:nz,0:nr,kqm),rhoti  (0:nz,0:nr,kqm),paircr(0:nz,0:nr,kqm)                         &
        )
IF(tregpair) ALLOCATE(regpairadd(0:nz,0:nr,kqm))

         pot   =0.0d0; coulpt =0.0d0 ;bmass =0.0d0; potls =0.0d0
         drbmas=0.0d0; dzbmas =0.0d0 ;potadd=0.0d0; rhosur=0.0d0 
         pairpt=0.0d0; potdamp=0.0d0 ;rho   =0.0d0; crho  =0.0d0
         tau   =0.0d0; drhols =0.0d0 ;rholsr=0.0d0; rholsz=0.0d0
         d2rho =0.0d0; rhoti  =0.0d0 ;paircr=0.0d0 

!NEW ODD
IF(isTsym>0) THEN
  ALLOCATE(&
    rotjr(0:nzd,0:nrd,kqm),  rotjz(0:nzd,0:nrd,kqm),  rotspden(0:nzd,0:nrd,kqm),&
    spdenr(0:nzd,0:nrd,kqm),  spdenz(0:nzd,0:nrd,kqm),  d2spdenr(0:nzd,0:nrd,kqm), &
    d2spdenz(0:nzd,0:nrd,kqm), curphi(0:nzd,0:nrd,kqm),  spdenkr(0:nzd,0:nrd,kqm), spdenkz(0:nzd,0:nrd,kqm)&
  )
  ALLOCATE(&
    potS(0:nzd,0:nrd,kqm,2),  potC(0:nzd,0:nrd,kqm,2),  potA(0:nzd,0:nrd,kqm)&
  )
  potS=0D0;  potC=0D0;  potA=0D0
  rotjr=0D0;  rotjz=0D0;  rotspden=0D0;  spdenr=0D0;  spdenz=0D0
  d2spdenr=0D0;  d2spdenz=0D0;  curphi=0D0;  spdenkr=0D0;  spdenkz=0D0
END IF
!END NEW ODD

IF(tregpair) regpairadd=0.0d0
RETURN
END SUBROUTINE alloc_fields

END MODULE fields
