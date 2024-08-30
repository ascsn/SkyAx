
!------------------------------------------------------------------------------
! MODULE: grid
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the grid variables as well as the matrices and routines 
!!for computing derivatives and integrals.
!> @details
!!Derivatives can be calculated either using a standart n-point formula or the 
!!Fourier definition. In both cases matrices are constructed and functions have 
!!to be multiplied with those to obtain the derivatives. Special attention has 
!!to be given to symmetries. The symmetry in z-direction is controlled with the 
!!variable \c tsym, \c isym, and \c sym). 
!!The logical varaibel \c tsym indicates symmetry or not in the $z$-direction.
!! \c sym and \c isym are double and integer variables which are 1 for symmetry and 0 for asymmetry. These are used to make some calculations more elegant.
!!Also the z-parity and m-component is taken into 
!!account.
!------------------------------------------------------------------------------
!
!                                                   tsym true   tsym flase
!       z-parity of spinor function:               +1    -1     none
!       z-parity label 'is' of spinor-function:     1     2       1
!       z-parity label 'is' of upper component:     1     2       1
!       z-parity label 'is' of lower component:     2     1       1

!       K-value of spinor-function:                 m+1/2
!       'm'-label of spinor function:                 m
!       m-value of upper component:                   m
!       m-value of lower component:                  m+1
!
!       The last label of z-derivative arrays (derrz1,...) stands for
!       z-parity in combination with Fourier basis (sin or cos). The
!       assignment is, for the example of derz1:
!            matrix        z-parity-label  basis
!          derz1(:,:,1)          1          sin
!          derz1(:,:,2)          2          sin
!          derz1(:,:,3)          1          cos
!          derz1(:,:,4)          2          cos
!       The entries 3 and 4 are provided for linear response
!       and play no role for the present static part of the code.
!
MODULE grid

USE paramr

IMPLICIT NONE

REAL(8) :: dr=0.7D0    !<Grid spacing in r-direction.
REAL(8) :: dz=0.7D0    !<Grid spacing in z-direction.
REAL(8) :: sym=1D0   !<double precision variable for symmetry
INTEGER :: nr=0      !<Grid size in r-direction.
INTEGER :: nz=0      !<Grid size in z-direction
INTEGER :: idrtyp=0  !<Switch for type of derivative.

REAL(8) :: sr=0D0    !<First point in r-direction.
REAL(8) :: sz        !<First point in z-direction.
REAL(8) :: dvol      !<Volume element for space integration.
INTEGER :: nr_testp  !<Output length in r-direction for printing (test switch).
INTEGER :: nz_testp  !<Output length in z-direction for printing (test switch).
LOGICAL  :: tsym=.TRUE.  ! controls symmetry
INTEGER :: isym=1    !<Integer Variable for symmetry.
INTEGER :: icm       !<Point on z-grid for Center of nucleus. If tsym=.FALSE.: 1, if tsym=.TRUE.: nz/2+1.

REAL(8),ALLOCATABLE :: r(:)                !<r-coordinates
REAL(8),ALLOCATABLE :: z(:)                !<z-coordinates

REAL(8),ALLOCATABLE :: derz1(:,:,:)        !<First derivative in z-direction.
REAL(8),ALLOCATABLE :: derz2(:,:,:)        !<Second derivative in z-direction.
REAL(8),ALLOCATABLE :: derr1(:,:,:)        !<First derivative in r-direction.
REAL(8),ALLOCATABLE :: derr2a(:,:,:,:)     !<Bessel operator  in r-direction.
REAL(8),ALLOCATABLE :: derz2f(:,:,:)       !<Second derivative in z-direction (fourier definition).
REAL(8),ALLOCATABLE :: derr1f(:,:,:)       !<First derivative in r-direction (fourier definition).
REAL(8),ALLOCATABLE :: derr2af(:,:,:,:)    !<Bessel operator  in r-direction (fourier definition).
REAL(8),ALLOCATABLE :: gaufolz(:,:,:)      !<Gaussian folding in z-direction.
REAL(8),ALLOCATABLE :: gaufolr(:,:,:,:)    !<Gaussian folding in r-direction.
REAL(8),ALLOCATABLE :: r1inv(:,:)          !<\f$ 1/r\f$.
REAL(8),ALLOCATABLE :: r2inv(:,:)          !<\f$ 1/r^2\f$.
REAL(8),ALLOCATABLE :: zdnuc(:,:,:)        !<Damping operator for the nucleons in z-direction.
REAL(8),ALLOCATABLE :: rdnuc(:,:,:)        !<Damping operator for the nucleons in r-direction.
REAL(8),ALLOCATABLE :: zdmes(:,:,:)        !<Damping operator for the photon in z-direction.
REAL(8),ALLOCATABLE :: rdmes(:,:,:)        !<Damping operator for the photon in r-direction.

INTEGER :: iband_derr                      !<Band width of derivatives in r-direction.
INTEGER :: iband_derz                      !<Band width of derivatives in z-direction.
REAL(8), ALLOCATABLE  :: int_weight(:,:)   !<Integration weight for scalar products.

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: init_grid
!> @brief
!!Initialization of the grid operations routines.
!> @details
!!This routine determines the band widths for derivatives and calls the 
!!appropriate routines for initializations for the derivative matrices.
!!  - For idrtyp = 0: Fourier derivatives.
!!  - For idrtyp = 1: 9-point in r-direction, Fourier in z-direction.
!!  - For idrtyp = 2: Fourier in r-direction, 9-point in z-direction.
!!  - For idrtyp = 3: 9-point derivatives.
!--------------------------------------------------------------------------- 

SUBROUTINE init_grid()

IMPLICIT NONE

INTEGER     :: iderpntr, iderpntz !-1 for Fourier and k for k-point derivative
INTEGER     :: ir,iz

ALLOCATE(derz1(0:nz,0:nz,4), derz2(0:nz,0:nz,4),         derr1(0:nr,0:nr,4),&
         derz2f(0:nz,0:nz,4),derr2a(0:nr,0:nr,0:kmpsi,2),derr1f(0:nr,0:nr,4),&
         derr2af(0:nr,0:nr,0:kmpsi,2), r1inv(0:nr,0:nr), r2inv(0:nr,0:nr),&
         zdnuc(0:nz,0:nz,4), rdnuc(0:nr,0:nr,0:kmpsi),   zdmes(0:nz,0:nz,1),&
         rdmes(0:nr,0:nr,1),r(0:nr),                     z(0:nz))
         
derz1=0D0;   derz2=0D0;  derr1=0D0
derz2f=0D0;  derr2a=0D0; derr1f=0D0
derr2af=0D0; r1inv=0D0;  r2inv=0D0
zdnuc=0D0;   rdnuc=0D0;  zdmes=0D0
rdmes=0D0;   r=0D0;       z=0D0

IF (idrtyp == 0) THEN
  iderpntr=-1
  iderpntz=-1
ELSE IF (idrtyp == 1) THEN
  iderpntr=9
  iderpntz=-1
ELSE IF (idrtyp == 2) THEN
  iderpntr=-1
  iderpntz=9
ELSE IF (idrtyp == 3) THEN
  iderpntr=9
  iderpntz=9
ELSE
  STOP 'Derivative mode is not correctly defined.'
END IF

IF(itest >= 2) THEN !set dimensions for test print
  nr_testp   = nr
  nz_testp   = nz
ELSE
  nr_testp   = 9
  nz_testp   = 9
END IF

dvol  = dr*dr*dz*two*pi
!isym  = NINT(sym)
IF(tsym) then
   icm=0
ELSE
   icm   = (nz+1)/2
END IF
sz    = zero
sr    = zero

!     Check input for derivatives and printout
!     ------------------------------------------------------------------
IF (iderpntz == -1) THEN
  iband_derz = -1
  IF(itest>0) WRITE(7,'(a)') '    z - derivatives : FOURIER'
ELSE IF ((MOD(iderpntz,2) == 1).AND.(iderpntz >= 2)) THEN
  iband_derz = (iderpntz-1)/2
  IF(itest>0) WRITE(7,'(a,i2,a)') '    z - derivatives : FINITE (',iderpntz,' points)'
ELSE
  WRITE(7,'(a)') ' type of derivatives must be -1 for fourier or'//  &
      ' 2n+1 for finite differences'
  WRITE(6,'(a)') ' type of derivatives must be -1 for fourier or'//  &
      ' 2n+1 for finite differences'
  STOP 'unsupported derivative type in INI_GRID'
END IF

IF (iderpntr == -1) THEN
  iband_derr = -1
  IF(itest>0) WRITE(7,'(a)') '    r - derivatives : FOURIER'
ELSE IF ((MOD(iderpntr,2) == 1).AND.(iderpntr >= 2)) THEN
  iband_derr = (iderpntr-1)/2
  IF(itest>0) WRITE(7,'(a,i2,a)') '    r - derivatives : FINITE (',iderpntr,' points)'
ELSE
  WRITE(7,'(a)') ' type of derivatives must be -1 for fourier or'//  &
      ' 2n+1 for finite differences'
  WRITE(6,'(a)') ' type of derivatives must be -1 for fourier or'//  &
      ' 2n+1 for finite differences'
  STOP 'unsupported derivative type in INI_GRID'
END IF

!     initialize derivatives in z-direction
!     ------------------------------------------------------------------
IF(itest>1) WRITE(*,*) ' before DERZ1'
CALL init_derive(derz1  ,1,iband_derz,tsym,dz,nz,nz)
IF(itest>1) WRITE(*,*) ' before DERZ2'
CALL init_derive(derz2  ,2,iband_derz,tsym,dz,nz,nz)
IF(itest>1) WRITE(*,*) ' before DERZ2F'
CALL init_derive(derz2f ,2,-1        ,tsym,dz,nz,nz)  ! Fourier definition used for damping operators

!     initialize radial derivatives
!     ------------------------------------------------------------------
IF(itest>1) WRITE(*,*) ' before DERR1'
CALL init_derive(derr1  ,1,iband_derr,   .TRUE.,dr,nr,nr)
IF(itest>1) WRITE(*,*) ' before DERR1F'
CALL init_derive(derr1f ,1,-1        ,   .TRUE.,dr,nr,nr)
IF(itest>1) WRITE(*,*) ' before DERR2A'
CALL init_derive(derr2a ,0,iband_derr,   .TRUE.,dr,nr,nr)
IF(itest>1) WRITE(*,*) ' before DERR2AF'
CALL init_derive(derr2af,0,-1        ,   .TRUE.,dr,nr,nr)

!     initialize r, z, 1/r and 1/r^2
!     ------------------------------------------------------------------
r1inv=0.0D0
r2inv=0.0D0
DO ir=0,nr
  r(ir)=REAL(ir+sr)*dr
END DO
DO iz=0,nz
  z(iz)=REAL(iz+sz-icm)*dz
END DO
DO ir=1,nr
  r1inv(ir,ir)=1.0D0/((ir+sr)*dr)
  r2inv(ir,ir)=r1inv(ir,ir)*r1inv(ir,ir)
END DO

RETURN
END SUBROUTINE init_grid


!---------------------------------------------------------------------------  
! DESCRIPTION: nablaz_scalar
!> @brief
!! Calculate the z-component of the gradient of a scalar function.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] izpfld
!> INTEGER, takes the z-parity 1 for positive, 2 for negative and 1 and \c
!! isym = 0 for no good parity. 
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 
SUBROUTINE nablaz_scalar(fldin,fldout,izpfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)         :: fldin(0:nz,0:nr)
REAL(8), INTENT(OUT)        :: fldout(0:nz,0:nr)
INTEGER, INTENT(IN)         :: izpfld
INTEGER, INTENT(IN)         :: ibpfld
REAL(8), INTENT(IN)         :: factor

INTEGER :: izpidx,ibpidx

izpidx = izpfld*isym + (1-isym)
ibpidx = izpidx + 2*ibpfld - 2

CALL mulz(derz1(:,:,ibpidx), fldin,fldout,factor,0)

RETURN
END SUBROUTINE nablaz_scalar

!---------------------------------------------------------------------------  
! DESCRIPTION: nablar_scalar
!> @brief
!! Calculate the r-component of the gradient of a scalar function.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] imfld
!> INTEGER, takes the m-component of the field.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 
SUBROUTINE nablar_scalar(fldin,fldout,imfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr)
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: irpidx,ibpidx

irpidx = 1+MOD(ABS(imfld),2)
ibpidx = irpidx + 2*ibpfld - 2

CALL mulr(derr1(:,:,ibpidx), fldin,fldout,factor,0)

RETURN
END SUBROUTINE nablar_scalar

!---------------------------------------------------------------------------  
! DESCRIPTION: nablaz_spin
!> @brief
!! Calculate the z-component of the gradient of a spinor.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] izpfld
!> INTEGER, takes the z-parity of upper spinor (1 for positive, 2 for 
!!negative, 1 and isym=0 for no good z-parity).
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE nablaz_spin(fldin,fldout,izpfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)        :: izpfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: is,izpact,izpidx,ibpidx

DO is=1,nsd
  izpact = 3*(izpfld-1) - (2*izpfld-3)*is
  izpidx = izpact*isym  + (1-isym)
  ibpidx = izpidx + 2*ibpfld - 2
  
  CALL mulz(derz1(:,:,ibpidx),fldin(:,:,is),fldout(:,:,is),factor,0)
END DO

RETURN
END SUBROUTINE nablaz_spin

!---------------------------------------------------------------------------  
! DESCRIPTION: nablar_spin
!> @brief
!! Calculate the r-component of the gradient of a spinor.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] imfld
!> INTEGER, takes the m-component of the upper component.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE nablar_spin(fldin,fldout,imfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: imact,irpidx,ibpidx,is

DO is=1,nsd
  imact     = imfld+is-1
  irpidx    = 1+ABS(MOD(imact,2))
  ibpidx    = irpidx + 2*ibpfld - 2
  
  CALL mulr(derr1(:,:,ibpidx),fldin(:,:,is),fldout(:,:,is),factor,0)
END DO

RETURN
END SUBROUTINE nablar_spin

!---------------------------------------------------------------------------  
! DESCRIPTION: nablap_spin
!> @brief
!! Calculate the \f$\phi\f$-component of the gradient of a spinor.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] imfld
!> INTEGER, takes the m-component of the upper component.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE nablap_spin(fldin,fldout,imfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: imact,irpidx,ibpidx,is

DO is=1,nsd
  imact     = imfld+is-1
  irpidx    = 1+ABS(MOD(imact,2))
  ibpidx    = irpidx + 2*ibpfld - 2
  
  CALL mulr(r1inv(:,:),fldin(:,:,is),fldout(:,:,is),factor*imact,0)
  IF (imact == 1.OR.imact == -1) THEN
    CALL mulr0(derr1(:,:,ibpidx), fldin(:,:,is),fldout(:,:,is),factor*imact,0)
  END IF
END DO

RETURN
END SUBROUTINE nablap_spin

!---------------------------------------------------------------------------  
! DESCRIPTION: laplace_scalar
!> @brief
!! Calculate the action of the laplace operator.
!!\f[ f*[(d/dr)^2 + 1/r (d/dr) - m^2/r^2 + (d/dz)^2] \f]
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] izpfld
!> INTEGER, takes the z-parity of field (1 for positive, 2 for negative,
!!1 and isym=0 for no good z-parity)
!> @param[in] imfld
!> INTEGER, takes the m-component of the field.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE laplace_scalar(fldin,fldout,izpfld,imfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr)
INTEGER, INTENT(IN)        :: izpfld
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: izpidx, ibpidx

izpidx = izpfld*isym + (1-isym)
ibpidx = izpidx + 2*ibpfld - 2

CALL mulr(derr2a(:,:,imfld,ibpfld),fldin,fldout,factor,0)
CALL mulz(derz2(:,:,ibpidx),       fldin,fldout,factor,1)

RETURN
END SUBROUTINE laplace_scalar

!---------------------------------------------------------------------------  
! DESCRIPTION: laplace_spin
!> @brief
!! Calculate the action of the laplace operator on a spinor.
!!\f[ f*[(d/dr)^2 + 1/r (d/dr) - m^2/r^2 + (d/dz)^2] \f]
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] izpfld
!> INTEGER, takes the z-parity of upper component (1 for positive, 2 for negative,
!!1 and isym=0 for no good z-parity)
!> @param[in] imfld
!> INTEGER, takes the z-component of angular momentum m of the upper component.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE laplace_spin(fldin,fldout,izpfld,imfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)        :: izpfld
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

INTEGER :: imact,izpact,ibpidx,is,izpidx

DO is=1,nsd
  imact  = ABS(imfld+is-1)
  izpact = 3*(izpfld-1) - (2*izpfld-3)*is
  izpidx = izpact*isym  + (1-isym)
  ibpidx = izpidx + 2*ibpfld - 2
  
  CALL mulz(derz2(:,:,ibpidx),fldin(:,:,is),fldout(:,:,is),factor,0)
  CALL mulr(derr2a(:,:,imact,ibpfld),fldin(:,:,is),fldout(:,:,is),factor,1)

END DO

RETURN
END SUBROUTINE laplace_spin

!---------------------------------------------------------------------------  
! DESCRIPTION: laplacer_scalar
!> @brief
!! Calculate the action of the radial part of the laplace operator.
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] imfld
!> INTEGER, takes the m-component of the field.
!> @param[in] ibpfld
!> INTEGER, takes boundary parity. 1 for negative (zero-value), 2 for positive 
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE laplacer_scalar(fldin,fldout,imfld,ibpfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)        :: fldin(0:nz,0:nr)
REAL(8), INTENT(OUT)       :: fldout(0:nz,0:nr)
INTEGER, INTENT(IN)        :: imfld
INTEGER, INTENT(IN)        :: ibpfld
REAL(8), INTENT(IN)        :: factor

CALL mulr(derr2a(:,:,imfld,ibpfld),fldin,fldout,factor,0)

RETURN
END SUBROUTINE laplacer_scalar

!---------------------------------------------------------------------------  
! DESCRIPTION: damp_meson
!> @brief
!! Calculate action of the meson-damping operator on a scalar function.
!!\f[ \frac{1}{\texttt{dams}^2} \cdot \frac{1}{1-\delta_r/\texttt{dams}^2} 
!!\cdot \frac{1}{1-\delta_z\texttt{dams}^2} \f]
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 
SUBROUTINE damp_meson(fldin,fldout,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)         :: fldin(0:nz,0:nr)
REAL(8), INTENT(OUT)        :: fldout(0:nz,0:nr)
REAL(8), INTENT(IN)         :: factor

REAL(8)                     :: fldbuf(0:nz,0:nr)

CALL mulz(zdmes(:,:,1), fldin,fldbuf ,one,0)
CALL mulr(rdmes(:,:,1), fldbuf,fldout,factor,0)

RETURN
END SUBROUTINE damp_meson

!---------------------------------------------------------------------------  
! DESCRIPTION: damp_nucleon
!> @brief
!! Calculate action of the nucleon-damping operator on a spinor.
!!\f[ \frac{\texttt{x0inv}}{\texttt{e0inv}} \cdot \frac{1}{1+t_r/\texttt{e0inv}} 
!!\cdot \frac{1}{1+t_z\texttt{e0inv}} \f]
!>
!> @param[in] fldin
!> REAL(8), array, takes the input field.
!> @param[out] fldout
!> REAL(8), array, returns the output field.
!> @param[in] izpfld
!> INTEGER, takes the z-parity of upper spinor.
!> @param[in] imfld
!> INTEGER, takes the m-component of the field.
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE damp_nucleon(dampr,dampz,fldin,fldout,izpfld,imfld,factor)

IMPLICIT NONE

REAL(8), INTENT(IN)      :: dampr(0:nr,0:nr,0:kmpsi)
REAL(8), INTENT(IN)      :: dampz(0:nz,0:nz,4)
REAL(8), INTENT(IN)      :: fldin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)     :: fldout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)      :: izpfld
INTEGER, INTENT(IN)      :: imfld
REAL(8), INTENT(IN)      :: factor

REAL(8)                     :: fldbuf(0:nz,0:nr)

INTEGER :: imact
INTEGER :: izpact
INTEGER :: is
INTEGER :: isigma

DO is=1,nsd
  
! lower spin component has m-component imact+1
  imact  = ABS(imfld+is-1)
  
! lower spin component has opposite z-parity if isym=1
! if isym=0 izpact should be always +1
  isigma = 3-2*is
  izpact = (3*(1-isigma)/2+izpfld*isigma)*isym + (1-isym)
  
  CALL mulr(dampr(:,:,imact),  fldin(:,:,is),fldbuf        ,one   ,0)
  CALL mulz(dampz(:,:,izpact), fldbuf,       fldout(:,:,is),factor,0)
END DO

RETURN
END SUBROUTINE damp_nucleon


!---------------------------------------------------------------------------  
! DESCRIPTION: init_derive
!> @brief
!! Initialization of derivative matrices
!>
!> @param[out] d
!> REAL(8), array, returns derivative matrix.
!> @param[in] itype
!> INTEGER, type of derivative: 0: Bessel operator, 1: 1st derivative, 2: 
!!2nd derivative.
!> @param[in] ipnt
!> INTEGER, takes the band width of derivative.
!> @param[in] tsyma
!> INTEGER, switch for good parity for reflection at lower boundary.
!> @param[in] dx
!> REAL(8), takes spacing of grid points.
!> @param[in] ngrid
!> INTEGER, number of grid points.
!> @param[in] mgrid
!> INTEGER, total dimension of array to store derivatives.
!--------------------------------------------------------------------------- 
SUBROUTINE init_derive(d,itype,ipnt,tsyma,dx,ngrid,mgrid)

IMPLICIT NONE

REAL(8), INTENT(OUT)        :: d(0:mgrid,0:mgrid,*)
INTEGER, INTENT(IN)         :: itype
INTEGER, INTENT(IN)         :: ipnt
LOGICAL, INTENT(IN)         :: tsyma
REAL(8), INTENT(IN)         :: dx
INTEGER, INTENT(IN)         :: ngrid
INTEGER, INTENT(IN)         :: mgrid


INTEGER, PARAMETER :: imaxdrvrng=21

REAL(8),ALLOCATABLE:: der1(:,:,:)
REAL(8),ALLOCATABLE:: der2(:,:,:)
REAL(8)            :: d1(0:mgrid+1,0:mgrid+1,4)
REAL(8)            :: d2(0:mgrid+1,0:mgrid+1,4)
REAL(8),ALLOCATABLE:: db(:,:,:,:)

REAL(8),ALLOCATABLE:: sbase(:,:,:)
REAL(8),ALLOCATABLE:: cbase(:,:,:)
REAL(8),ALLOCATABLE:: volx(:), kval(:)

REAL(8)            :: a(0:2*imaxdrvrng,0:2*imaxdrvrng)
REAL(8)            :: b(0:2*imaxdrvrng)

INTEGER            :: iparmax, ibuf,in
INTEGER            :: inx,k,nx,i,im,ix,j,jx
REAL(8)            :: buf,buf1,buf2

IF(itest>1) WRITE(*,*) ' ini_deriv entered'
in=2

IF (.NOT.tsyma) THEN
  inx=ngrid+1
ELSE
  inx=2*ngrid+2
END IF

ALLOCATE(db(0:mgrid+1,0:mgrid+1,0:kmpsi,2),der1(0:inx,0:inx,2),der2(0:inx,0:inx,2),&
         sbase(0:2,0:inx,0:inx),cbase(0:2,0:inx,0:inx),volx(0:inx),kval(0:inx))
db=0D0;   der1=0D0; der2=0D0
sbase=0D0;cbase=0D0;volx=0D0;kval=0D0

!     ------------------------------------------------------------------
!     check consistency of input parameters
!     ------------------------------------------------------------------
IF ((ipnt /= -1).AND.(ipnt < (in+1)/2)) THEN
  WRITE(6,'(i2,a,i2,a)')  &
      in,'. order derivative needs at least ',(in+1)+1,' points'
  STOP 'error in init_derive'
END IF

IF ((2*ipnt+1) > inx-2) THEN
  WRITE(6,'(i2,a,i2,a)') (2*ipnt+1),' formula needs at least ',  &
      (2*ipnt+3),' gridpoints'
  STOP 'error in init_derive'
END IF

IF ((ipnt < -1).OR.(ipnt > imaxdrvrng).OR.(ipnt == 0)) THEN
  WRITE(6,'(a,i2)') 'number of points must be -1,1,2,...,',imaxdrvrng
  STOP 'error in init_derive'
END IF

IF (dx <= 0.0) THEN
  WRITE(6,'(a)') 'spacing of gridpoints must be >0'
  STOP 'error in init_derive'
END IF

!     ------------------------------------------------------------------
!
!     FOURIER derivatives
!
!     ------------------------------------------------------------------
IF (ipnt == -1) THEN
!        ---------------------------------------------------------------
!        initialize cos and sin-base
!        ---------------------------------------------------------------
  DO k=0,inx
    
    kval(k) = k*pi/(inx*dx)
    
    buf1 = SQRT(2.0D0/(inx*dx))
    IF (k == 0 .OR. k == inx) THEN
      buf2 = SQRT(1.0D0/(inx*dx))
      buf1 = 0.0D0
    ELSE
      buf2 = buf1
    END IF
    
    DO i=1,inx-1
      sbase(0,k,i) =  buf1 *  SIN(k*i*pi/inx)
      cbase(0,k,i) =  buf2 *  COS(k*i*pi/inx)
      sbase(1,k,i) =  kval(k)    * cbase(0,k,i)
      cbase(1,k,i) = -kval(k)    * sbase(0,k,i)
      sbase(2,k,i) = -kval(k)**2 * sbase(0,k,i)
      cbase(2,k,i) = -kval(k)**2 * cbase(0,k,i)
    END DO
    
    sbase(0,k,0)    =  0.0D0
    cbase(0,k,0)    =  buf2
    sbase(1,k,0)    =  kval(k)    * buf2
    cbase(1,k,0)    =  0.0D0
    sbase(2,k,0)    =  0.0D0
    cbase(2,k,0)    = -kval(k)**2 * buf2
    sbase(0,k,inx)  =  0.0D0
    cbase(0,k,inx)  =  buf2 *(-1)**k
    sbase(1,k,inx)  =  kval(k)    * buf2 *(-1)**k
    cbase(1,k,inx)  =  0.0D0
    sbase(2,k,inx)  =  0.0D0
    cbase(2,k,inx)  = -kval(k)**2 * buf2 *(-1)**k
    
  END DO
  IF(itest>1) WRITE(*,*) ' after cos basis'
  
!        ---------------------------------------------------------------
!        initialize volume element (cartesian)
!        ---------------------------------------------------------------

  volx(0)       =  0.5D0*dx
  volx(1:inx-1) =  1.0D0*dx
  volx(inx)     =  0.5D0*dx
  
!        ---------------------------------------------------------------
!        initialize 1st and 2nd fourier-derivative
!        ---------------------------------------------------------------
  der1(0:inx,0:inx,:)=0.0D0
  der2(0:inx,0:inx,:)=0.0D0
  IF(itest>1) WRITE(*,*) ' after zeroing ders'
  
  DO i=0,inx
    DO j=0,inx
      der1(j,i,1) = der1(j,i,1) + sum(sbase(1,0:inx,i)*sbase(0,0:inx,j)*volx(j))
      der2(j,i,1) = der2(j,i,1) + sum(sbase(2,0:inx,i)*sbase(0,0:inx,j)*volx(j))
      der1(j,i,2) = der1(j,i,2) + sum(cbase(1,0:inx,i)*cbase(0,0:inx,j)*volx(j))
      der2(j,i,2) = der2(j,i,2) + sum(cbase(2,0:inx,i)*cbase(0,0:inx,j)*volx(j))
    END DO
  END DO
  IF(itest>1) WRITE(*,*) ' after ders'
  
!     ------------------------------------------------------------------
!
!     FINITE derivatives
!
!     ------------------------------------------------------------------
ELSE    
  DO iparmax=-1,+1,2
    DO in=1,2
!        ---------------------------------------------------------------
!        determin weights for the n-th finite derivative using p points
!        into each direction from the central point ((2*p+1)-formula).
!        ---------------------------------------------------------------         
!           setup the linear equation system which defines the weights
      DO j=0,2*ipnt
        DO i=-ipnt,ipnt
          a(j,i+ipnt) = i**j
        END DO
        IF (j == in) THEN
          b(j)=1.0
        ELSE
          b(j)=0.0
        END IF
      END DO
      
!           forward gauss elimination
      DO i=0,2*ipnt
        DO j=i+1,2*ipnt
          buf = a(j,i) / a(i,i)
          DO k=i+1,2*ipnt
            a(j,k) = a(j,k) - a(i,k) * buf
          END DO
          b(j) = b(j) - b(i) * buf
        END DO
      END DO
      
!           backward substitution
      DO i=2*ipnt,0,-1
        buf = 0.0
        DO j=i+1,2*ipnt
          buf = buf + a(i,j) * b(j)
        END DO
        b(i) = (b(i) - buf) / a(i,i)
      END DO
      
!           calculate faculty of order n
      buf=1.0
      DO i=1,in
        buf = 1.0/dx * buf * i
      END DO
      
!           compose the weights for the finite differences
      DO i=0,2*ipnt
        b(i) = buf * b(i)
      END DO
      
!           setup 1st finite derivative matrix
      IF (in == 1) THEN
        
        IF (iparmax == -1) THEN
          k=1
        ELSE
          k=2
        END IF
        
        der1(0:inx,0:inx,k)=0.0

        IF (.NOT.tsyma) THEN
          DO i=1,inx
            DO j=-ipnt,ipnt
! straight forward mid region
              IF (((i+j) > 0).AND.((i+j) <= inx)) THEN
                der1(i+j,i,k) = der1(i+j,i,k)   +         b(j+ipnt)
! take respect of parity at lower boundary
              ELSE IF ((i+j) < 0) THEN
                der1(-(i+j),i,k) = der1(-(i+j),i,k)+ iparmax*b(j+ipnt)
! take respect of parity at upper boundary
              ELSE IF ((i+j) > inx) THEN
                der1(2*inx-(i+j),i,k)  &
                    = der1(2*inx-(i+j),i,k) + iparmax*b(j+ipnt)
              END IF
            END DO
          END DO
        ELSE
          DO i=0,inx
            DO j=-ipnt,ipnt
! straight forward mid region
              IF (((i+j) >= 0).AND.((i+j) <= inx)) THEN
                der1(i+j,i,k) = der1(i+j,i,k)   +         b(j+ipnt)
! take respect of parity at lower boundary
              ELSE IF ((i+j) < 0) THEN
                der1(-(i+j),i,k) = der1(-(i+j),i,k)+ iparmax*b(j+ipnt)
! take respect of parity at upper boundary
              ELSE IF ((i+j) > inx) THEN
                der1(2*inx-(i+j),i,k)  &
                    = der1(2*inx-(i+j),i,k) + iparmax*b(j+ipnt)
              END IF
            END DO
          END DO
        END IF
        
!           setup 2nd finite derivative matrix
      ELSE
        
        IF (iparmax == -1) THEN
          k=1
        ELSE
          k=2
        END IF
        
        der2(0:inx,0:inx,k)=0.0
        
        IF (.NOT.tsyma) THEN
          DO i=1,inx
            DO j=-ipnt,ipnt
! straight forward mid region
              IF (((i+j) > 0).AND.((i+j) <= inx)) THEN
                der2(i+j,i,k) = der2(i+j,i,k)   +         b(j+ipnt)
! take respect of parity at lower boundary
              ELSE IF ((i+j) < 0) THEN
                der2(-(i+j),i,k) = der2(-(i+j),i,k)+ iparmax*b(j+ipnt)
! take respect of parity at upper boundary
              ELSE IF ((i+j) > inx) THEN
                der2(2*inx-(i+j),i,k)  &
                    = der2(2*inx-(i+j),i,k) + iparmax*b(j+ipnt)
              END IF
            END DO
          END DO
        ELSE
          DO i=0,inx
            DO j=-ipnt,ipnt
! straight forward mid region
              IF (((i+j) >= 0).AND.((i+j) <= inx)) THEN
                der2(i+j,i,k) = der2(i+j,i,k)   +         b(j+ipnt)
! take respect of parity at lower boundary
              ELSE IF ((i+j) < 0) THEN
                der2(-(i+j),i,k) = der2(-(i+j),i,k)+ iparmax*b(j+ipnt)
! take respect of parity at upper boundary
              ELSE IF ((i+j) > inx) THEN
                der2(2*inx-(i+j),i,k)  &
                    = der2(2*inx-(i+j),i,k) + iparmax*b(j+ipnt)
              END IF
            END DO
          END DO
        END IF
        
      END IF
      
    END DO
  END DO
  
  IF(itest>1) WRITE(*,*) ' after finite differences'
END IF  ! end of finite definition of derivatives

!     ===============================================================
!     reduce to grid with/without symmetry at x=0
!     ===============================================================
IF (tsyma) THEN
  DO i=inx/2,inx
    DO j=inx/2,inx
      d1(i-inx/2,j-inx/2,1)=der1(i,j,1) + der1(inx-i,j,1)
      d1(i-inx/2,j-inx/2,2)=der1(i,j,1) - der1(inx-i,j,1)
      d1(i-inx/2,j-inx/2,3)=der1(i,j,2) + der1(inx-i,j,2)
      d1(i-inx/2,j-inx/2,4)=der1(i,j,2) - der1(inx-i,j,2)
      
      d2(i-inx/2,j-inx/2,1)=der2(i,j,1) + der2(inx-i,j,1)
      d2(i-inx/2,j-inx/2,2)=der2(i,j,1) - der2(inx-i,j,1)
      d2(i-inx/2,j-inx/2,3)=der2(i,j,2) + der2(inx-i,j,2)
      d2(i-inx/2,j-inx/2,4)=der2(i,j,2) - der2(inx-i,j,2)
    END DO
  END DO
  d1(0,0:inx/2,:)=0.5D0*d1(0,0:inx/2,:)
  d2(0,0:inx/2,:)=0.5D0*d2(0,0:inx/2,:)
  
  d1(0:inx/2,0,1)=0.0D0
  d1(0:inx/2,0,3)=0.0D0
  d2(0:inx/2,0,2)=0.0D0
  d2(0:inx/2,0,4)=0.0D0
  
ELSE  
  d1(0:inx,0:inx,1)=der1(0:inx,0:inx,1)
  d1(0:inx,0:inx,3)=der1(0:inx,0:inx,2)
      
  d2(0:inx,0:inx,1)=der2(0:inx,0:inx,1)
  d2(0:inx,0:inx,3)=der2(0:inx,0:inx,2)
END IF

!     ===============================================================
!     initialize bessel operator
!     ===============================================================
DO im=0,kmpsi
  
  DO i=0,inx/2
    IF(i /= 0) THEN
      kval(i)=1.0D0/(i*dx)
    ELSE
      kval(i)=0.d0
    END IF
    db(0:inx/2,i,im,1) = 0.0D0
    db(0:inx/2,i,im,2) = 0.0D0
  END DO
  
  IF (MOD(im,2) == 0) THEN
    
    DO i=0,inx/2
      DO j=0,inx/2
        db(j,i,im,1) = db(j,i,im,1) + kval(i)*d1(j,i,1) + d2(j,i,1)
        db(j,i,im,2) = db(j,i,im,2) + kval(i)*d1(j,i,3) + d2(j,i,3)
      END DO
    END DO
    
    IF (im == 0) THEN
      db(0:inx/2,0,im,1) = 2.0D0*d2(0:inx/2,0,1)
      db(0:inx/2,0,im,2) = 2.0D0*d2(0:inx/2,0,3)
    ELSE
      db(0:inx/2,0,im,1) = 0.0D0
      db(0:inx/2,0,im,2) = 0.0D0
      
      DO j=0,inx/2
        db(j,j,im,1) = db(j,j,im,1) - im*im*kval(j)*kval(j)
        db(j,j,im,2) = db(j,j,im,2) - im*im*kval(j)*kval(j)
      END DO
    END IF
    
  ELSE
    
    DO i=0,inx/2
      DO j=0,inx/2
        db(j,i,im,1) = db(j,i,im,1) + kval(i)*d1(j,i,2) + d2(j,i,2)
        db(j,i,im,2) = db(j,i,im,2) + kval(i)*d1(j,i,4) + d2(j,i,4)
      END DO
    END DO
    
    db(0:inx/2,0,im,1) = 0.0D0
    db(0:inx/2,0,im,2) = 0.0D0
    
    DO j=0,inx/2
      db(j,j,im,1) = db(j,j,im,1) - im*im*kval(j)*kval(j)
      db(j,j,im,2) = db(j,j,im,2) - im*im*kval(j)*kval(j)
    END DO
    
  END IF
  
END DO

!     ---------------------------------------------------------------
!     return derivative matrix
!     ---------------------------------------------------------------
IF (.NOT.tsyma) THEN
  nx = inx/2
ELSE
  nx = inx
END IF

IF (itype == 1) THEN
  DO jx=0,ngrid
    DO ix=0,ngrid
      d(ix,jx,1:4)=d1(jx,ix,1:4)
    END DO
  END DO
  
ELSE IF (itype == 2) THEN
  DO jx=0,ngrid
    DO ix=0,ngrid
      d(ix,jx,1:4)=d2(jx,ix,1:4)
    END DO
  END DO
  
ELSE IF (itype == 0) THEN
  DO im=0,kmpsi
    ibuf = im+1
    DO jx=0,ngrid
      DO ix=0,ngrid
        d(ix,jx,ibuf        )=db(jx,ix,im,1)
        d(ix,jx,ibuf+kmpsi+1)=db(jx,ix,im,2)
      END DO
    END DO
  END DO

END IF

DEALLOCATE(db,der1,der2,sbase,cbase,volx,kval)
RETURN
END SUBROUTINE init_derive

!---------------------------------------------------------------------------  
! DESCRIPTION: f2int
!> @brief
!!Computes the integral over \c fua * \c fub.
!>
!> @param[in] fua
!> REAL(8), array, takes first input function.
!> @param[in] fub
!> REAL(8), array, takes second input function.
!--------------------------------------------------------------------------- 

REAL(8) FUNCTION f2int(fua,fub)

IMPLICIT NONE

REAL(8), EXTERNAL :: DDOT
REAL(8), INTENT(IN)                         :: fua(0:nz,0:nr)
REAL(8), INTENT(IN)                         :: fub(0:nz,0:nr)

  IF (.NOT.ALLOCATED(int_weight)) CALL init_int_weight
  !f2int=SUM(int_weight*fua*fub)  
  f2int=DDOT((nz+1)*(nr+1),fua*fub,1,int_weight,1)
RETURN
END FUNCTION f2int

!---------------------------------------------------------------------------  
! DESCRIPTION: f1int
!> @brief
!!Computes the integral over \c fua.
!>
!> @param[in] fua
!> REAL(8), array, takes input function.
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION f1int(fua)

IMPLICIT NONE

REAL(8), EXTERNAL :: DDOT
REAL(8), INTENT(IN)          :: fua(0:nz,0:nr)

  IF (.NOT.ALLOCATED(int_weight)) CALL init_int_weight
!  f1int = sum(fua * int_weight)
  f1int=DDOT((nz+1)*(nr+1),fua,1,int_weight,1)
RETURN
END FUNCTION f1int

!---------------------------------------------------------------------------  
! DESCRIPTION: matr_mul
!> @brief
!!Full matrix multiplication a*b=c.
!>
!> @param[in] a
!> REAL(8), array, takes the matrix a.
!> @param[in] a
!> REAL(8), array, takes the matrix b.
!> @param[out] c
!> REAL(8), array, returns the matrix c.
!> @param[in] imax
!> INTEGER, takes the actual dimension of matrices a, b and c.
!> @param[in] nmax
!> INTEGER, takes the maximal dimension of the matrices.
!--------------------------------------------------------------------------- 
SUBROUTINE matr_mul(a,b,c,imax,nmax)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nmax,0:nmax)
REAL(8), INTENT(IN)             :: b(0:nmax,0:nmax)
REAL(8), INTENT(OUT)            :: c(0:nmax,0:nmax)
INTEGER, INTENT(IN)             :: imax
INTEGER, INTENT(IN)             :: nmax

IF(imax > nmax) STOP ' dimension > "nmax" in MATRML'
STOP 'matrml'
CALL DGEMM('N','N',imax+1,imax+1,imax+1,one,a,nmax+1,b,nmax+1,zero,c,nmax+1)

RETURN
END SUBROUTINE matr_mul

!---------------------------------------------------------------------------  
! DESCRIPTION: matr_inv
!> @brief
!!Matrix inversion of matrix, implemented with LAPACK routines.
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix a and returns the inverse.
!> @param[in] imax
!> INTEGER, takes the actual dimension of matrix a.
!> @param[in] nmax
!> INTEGER, takes the maximal dimension of the matrices.
!--------------------------------------------------------------------------- 

SUBROUTINE matr_inv(a,nmax,imax)

IMPLICIT NONE

REAL(8), INTENT(INOUT)          :: a(nmax,nmax)
INTEGER, INTENT(IN)             :: nmax
INTEGER, INTENT(IN)             :: imax

INTEGER                         :: info,ipiv(imax)
REAL(8)                         :: work(imax)

!     lower-upper decomposition
!     ------------------------------------------------------------------
CALL DGETRF(imax,imax,a,nmax,ipiv,info)

!     backward substitution
!     ------------------------------------------------------------------
CALL DGETRI(imax,a,nmax,ipiv,work,imax,info)

RETURN
END SUBROUTINE matr_inv

!---------------------------------------------------------------------------  
! DESCRIPTION: mulz
!> @brief
!!Matrix-matrix multiplication for operations in z-direction. 
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix \c a.
!> @param[inout] x
!> REAL(8), array, takes the matrix \c x.
!> @param[inout] y
!> REAL(8), array, takes the matrix \c y.
!> @param[in] factor
!> REAL(8), takes the value \c factor.
!> @param[in] iaccu
!> INTEGER, 0 for no accumulation, 1 for accumulation of \c y.

!--------------------------------------------------------------------------- 

SUBROUTINE mulz(a,x,y,factor,iaccu)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nz,0:nz)
REAL(8), INTENT(IN)             :: x(0:nz,0:nr)
REAL(8), INTENT(INOUT)          :: y(0:nz,0:nr)
REAL(8), INTENT(IN)             :: factor
INTEGER, INTENT(IN)             :: iaccu

REAL(8) :: beta

beta=one
IF (iaccu == 0) beta = zero

IF (ABS(factor) > smal) CALL DGEMM('N','N',nz+1,nr+1,nz+1,factor,a,nz+1,x,nz+1,beta,y,nz+1) 

RETURN
END SUBROUTINE mulz

!---------------------------------------------------------------------------  
! DESCRIPTION: mulr
!> @brief
!!Matrix-matrix multiplication for operations in r-direction. 
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix \c a.
!> @param[inout] x
!> REAL(8), array, takes the matrix \c x.
!> @param[inout] y
!> REAL(8), array, takes the matrix \c y.
!> @param[in] factor
!> REAL(8), takes the value \c factor.
!> @param[in] iaccu
!> INTEGER, 0 for no accumulation, 1 for accumulation of \c y.
!--------------------------------------------------------------------------- 

SUBROUTINE mulr(a,x,y,factor,iaccu)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nr,0:nr)
REAL(8), INTENT(IN)             :: x(0:nz,0:nr)
REAL(8), INTENT(INOUT)          :: y(0:nz,0:nr)
REAL(8), INTENT(IN)             :: factor
INTEGER, INTENT(IN)             :: iaccu

REAL(8) :: beta

beta=one
IF (iaccu == 0) beta = zero

IF (ABS(factor) > smal) CALL DGEMM('N','T',nz+1,nr+1,nr+1,factor,x,nz+1,a,nr+1,beta,y,nz+1) 

RETURN
END SUBROUTINE mulr

!---------------------------------------------------------------------------  
! DESCRIPTION: mulr0
!> @brief
!!Matrix-matrix multiplication for operations in r-direction (only r=0). 
!>
!> @param[inout] a
!> REAL(8), array, takes the matrix \c a.
!> @param[inout] x
!> REAL(8), array, takes the matrix \c x.
!> @param[inout] y
!> REAL(8), array, takes the matrix \c y.
!> @param[in] factor
!> REAL(8), takes the value \c factor.
!> @param[in] iaccu
!> INTEGER, 0 for no accumulation, 1 for accumulation of \c y.
!--------------------------------------------------------------------------- 

SUBROUTINE mulr0(a,x,y,factor,iaccu)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: a(0:nr,0:nr)
REAL(8), INTENT(IN)             :: x(0:nz,0:nr)
REAL(8), INTENT(INOUT)          :: y(0:nz,0:nr)
REAL(8), INTENT(IN)             :: factor
INTEGER, INTENT(IN)             :: iaccu

REAL(8) :: beta

beta=one
IF (iaccu == 0) beta = zero

IF (ABS(factor) > smal) CALL DGEMV('N',nz+1,nr+1,factor,x,nz+1,a,nz+1,beta,y,1)

RETURN
END SUBROUTINE mulr0

!---------------------------------------------------------------------------  
! DESCRIPTION: init_int_weight
!> @brief
!!Initializes integration weight for wave function overlaps and scalar integrals.
!--------------------------------------------------------------------------- 
SUBROUTINE init_int_weight

  IMPLICIT NONE

  REAL(8)               :: r
  REAL(8), PARAMETER    :: drezw=1.0D0/12.0D0
  INTEGER               :: ir

  ALLOCATE (int_weight(0:nz,0:nr))
  int_weight = one
  int_weight(:,0) = int_weight(:,0) * drezw
  int_weight(0,:) = int_weight(0,:) * oh

  r = sr
  DO ir=1,nr
    r = one + r
    int_weight(0:nz,ir) = int_weight(0:nz,ir) * r
  END DO
  int_weight=int_weight*(isym+1)*dvol
END SUBROUTINE

!---------------------------------------------------------------------------  
! DESCRIPTION: test_print
!> @brief
!!Testprint of real 'field' in the range nzmin<iz<nzmax and nrmin<ir<nrmax. 
!!Maximally 10 ir are printed. Print is named by 'header'. Maximum and 
!!minimum value of \c field is evaluated and printed (in whole range). 
!>
!> @param[in] field
!> REAL(8), array, takes the field to print.
!> @param[in] nzmin
!> INTEGER, takes minimum z index for printing.
!> @param[in] nzmax
!> INTEGER, takes maximum z index for printing.
!> @param[in] nrmin
!> INTEGER, takes minimum r index for printing.
!> @param[in] nrmax
!> INTEGER, takes maximum r index for printing.
!> @param[in] header
!> CHARACTER, array, takes maximum r index for printing.
!--------------------------------------------------------------------------- 

SUBROUTINE test_print(field,nzmin,nzmax,nrmin,nrmax,header)

IMPLICIT NONE

REAL(8), INTENT(IN)                  :: field(0:nz,0:nr)
INTEGER, INTENT(IN)                  :: nzmin
INTEGER, INTENT(IN)                  :: nzmax
INTEGER, INTENT(IN)                  :: nrmin
INTEGER, INTENT(IN)                  :: nrmax
CHARACTER (LEN=*), INTENT(IN)        :: header

REAL(8) :: fldmax,fldmin
INTEGER :: ir,irmin,iz,nrpmax

WRITE(7,'(a,a,a,1pg13.5)') ' testprint of field:',header, &
   ',  total norm=',f1int(field(0:nz,0:nr))

DO  irmin=nrmin,nrmax,10
  nrpmax = MIN0(nrmax,irmin+9)
  WRITE(7,'(a,3x,a,i2,9(7x,a,i2))') ' nz',('ir=',ir,ir=irmin,nrpmax)
  DO  iz=nzmin,nzmax
    WRITE(7,'(i3,10g12.4)') iz,(field(iz,ir),ir=irmin,nrpmax)
  END DO
END DO

fldmin = +1.e30
fldmax = -1.e30
DO  ir=0,nr
  DO  iz=0,nz
    fldmin = MIN(ABS(field(iz,ir)),fldmin)
    fldmax = MAX(ABS(field(iz,ir)),fldmax)
  END DO
END DO
WRITE(7,'(/a,g12.4,a,g12.4)') ' abs maximum value = ',fldmax,  &
    ' abs minimum value = ',fldmin

RETURN
END SUBROUTINE test_print

END MODULE grid
