SUBROUTINE user_interface_1()
! user defined corrections to energy, radii etc  (called in 'skyax.f90')

! here: compute charge density by 2D Fourier transformation

USE paramr, ONLY: pi
USE grid, ONLY: nr,nz,z,r,icm,dr,dz,f1int
USE fields, ONLY: rho
USE dynpar, ONLY: outnam
USE observ, ONLY: cmwid,zcmobs

IMPLICIT NONE

REAL(8), PARAMETER    :: drezw=1.0D0/12.0D0
LOGICAL,PARAMETER :: tprintff=.FALSE.

REAL(8),ALLOCATABLE :: qr(:)          ! q-values for radial grid
REAL(8),ALLOCATABLE :: qz(:)          ! q-values in z direction
REAL(8),ALLOCATABLE :: form2D(:,:,:)    ! formfactor in 2D
REAL(8),ALLOCATABLE :: rhoC2D(:,:)    ! charge density in 2D
REAL(8),ALLOCATABLE :: rhoform(:,:)   ! workspace
REAL(8),ALLOCATABLE :: besstrans(:,:) ! Bessel functions for transformation
REAL(8),ALLOCATABLE :: sintrans(:,:)  ! sinus basis for transformation

REAL(8) :: qmax=20D0        ! maximum momentum in Fourier space

INTEGER :: nqrmax,nqzmax,n,iz,ir,iq,iqr,iqz
REAL(8) :: radmax,zmax,delqr,delqz,qact,facnorm
REAL(8) :: r2p,r2c,r4p,r4c,zcm

!  determine set of q-values, prepare Besselfunctions
  qmax = pi/dr
  radmax = r(nr)
  delqr = pi/radmax
  nqrmax=nr-1  !INT(qmax/delqr)-1
  zmax = z(nz)
  delqz = pi/(zmax)
  nqzmax=nz-1 !INT(qmax/delqz)-1
  ALLOCATE(qr(0:nqrmax),form2D(0:nqzmax,0:nqrmax,1:2))
  ALLOCATE(qz(0:nqzmax))
  ALLOCATE(besstrans(0:nr,0:nqrmax))
  ALLOCATE(sintrans(0:nz,0:nqzmax))
  ALLOCATE(rhoform(0:nqzmax,0:nr))
  ALLOCATE(rhoC2D(0:nz,0:nr))

  DO n=0,nqrmax
    qr(n) = besszero((n+0.5D0)*pi,(n+1.5D0)*pi)/radmax
    DO ir=0,nr
      besstrans(ir,n)=bessel_j0(qr(n)*r(ir))
    END DO
    facnorm =   SUM(besstrans(:,n)*besstrans(:,n)*r)+&
       dr*drezw*besstrans(0,n)*besstrans(0,n)
    besstrans(:,n)= besstrans(:,n)/SQRT(facnorm)
  END DO
  CALL orthon_1Dr()

  DO iq=0,nqzmax
    qact = (iq+0.5D0)*delqz
    qz(iq) = qact
    DO iz=0,nz
      sintrans(iz,iq)=cos(qact*z(iz))
    END DO
  END DO

! Fourier transformation, use rhoC
  DO iq=1,2
    DO ir=0,nr
      DO iqz=0,nqzmax  
        rhoform(iqz,ir) = SUM(sintrans(:,iqz)*rho(:,ir,1)) &
                -0.5D0*sintrans(0,iqz)*rho(0,ir,1)
      END DO
    END DO
    DO iqz=0,nqzmax
      DO iqr=0,nqrmax  
        form2D(iqz,iqr,iq) = SUM(besstrans(:,iqr)*rhoform(iqz,:)*r(:))&
         +dr*drezw*besstrans(0,iqr)*rhoform(iqz,0)
      END DO
    END DO
  END DO


! apply nucleon formfactors for 'form2D' in place

  DO iqz=0,nqzmax
    DO iqr=0,nqrmax  
      form2D(iqz,iqr,1) = &
         charge_formfactor(qr(iqr)**2+qz(iqz)**2, &
                           cmwid,form2D(iqz,iqr,1),form2D(iqz,iqr,2))
    END DO
  END DO



! Fourier backtransformation

  DO iqz=0,nqzmax
    DO ir=0,nr  
      rhoform(iqz,ir) = SUM(besstrans(ir,:)*form2D(iqz,:,1)) !*qr(:))
    END DO
  END DO

  DO ir=0,nr
    DO iz=0,nz  
      rhoC2D(iz,ir) = SUM(sintrans(iz,:)*rhoform(:,ir)) 
    END DO
  END DO

! normalize


  facnorm = f1int(rho(:,:,1))/f1int(rhoC2D)
  rhoC2D = facnorm*rhoC2D
  DEALLOCATE(rhoform)
  ALLOCATE(rhoform(0:nz,0:nr))
  zcm = zcmobs(rho(0:nz,0:nr,1))
  facnorm = f1int(rho(:,:,1))
  DO ir=0,nr; DO iz=0,nz
    rhoform(iz,ir) = rho(iz,ir,1)*(r(ir)**2+(z(iz)-zcm)**2)
  END DO; END DO
  r2p=SQRT(f1int(rhoform)/facnorm)
  DO ir=0,nr; DO iz=0,nz
    rhoform(iz,ir) = rho(iz,ir,1)*(r(ir)**2+(z(iz)-zcm)**2)**2
  END DO; END DO
  r4p=(f1int(rhoform)/facnorm)**0.25D0
  zcm = zcmobs(rhoC2D)
  DO ir=0,nr; DO iz=0,nz
    rhoform(iz,ir) = rhoC2D(iz,ir)*(r(ir)**2+(z(iz)-zcm)**2)
  END DO; END DO
  r2c=SQRT(f1int(rhoform)/facnorm)
  DO ir=0,nr; DO iz=0,nz
    rhoform(iz,ir) = rhoC2D(iz,ir)*(r(ir)**2+(z(iz)-zcm)**2)**2
  END DO; END DO
  r4c=(f1int(rhoform)/facnorm)**0.25D0

  OPEN(307,FILE='rhocharge.'//outnam)
  WRITE(307,'(a,4(1pg13.5))') '# radii: r_2C,r_2p,r_4C,r_4p=',r2c,r2p,r4c,r4p
  WRITE(307,'(a)') '# 2D print of charge density',&
          '#   z     r   rho_charge  rho_proton  difference'
  DO ir=0,nr 
    DO iz=0,nz
      WRITE(307,'(2f6.1,3(1pg12.4))') (z(iz)-zcm),r(ir),rhoC2D(iz,ir),&
             rho(iz,ir,1),rhoC2D(iz,ir)-rho(iz,ir,1)
    END DO
    WRITE(307,'(1x)')
  END DO

  IF(tprintff) THEN
    WRITE(307,'(1x)')
    WRITE(307,'(a,4i5)') '# nz,nqzmax,nr,nqrmax=',nz,nqzmax,nr,nqrmax
    WRITE(307,'(a,20(/5(1pg13.5)))') '# z:',z
    WRITE(307,'(a,20(/5(1pg13.5)))') '# r:',r
    WRITE(307,'(a,20(/5(1pg13.5)))') '# qr:',qr

    WRITE(307,*) '# besstrans'
    DO n=0,nqrmax
      WRITE(307,'(5(1pg13.5))') besstrans(:,n)
    END DO

    WRITE(307,*) '# sintrans'
    DO iq=0,nqzmax
      WRITE(307,'(5(1pg13.5))')  sintrans(:,iq)
    END DO

    WRITE(307,*) '# besstrans orthog'
    DO ir=0,nqrmax
      WRITE(307,'(5(1pg13.5))') &
         (SUM(besstrans(:,ir)*besstrans(:,n)*r)+&
         dr*drezw*besstrans(0,ir)*besstrans(0,n),n=0,nqrmax)
    END DO

    WRITE(307,*) '# sintrans orthog'
    DO ir=0,nqzmax
      WRITE(307,'(5(1pg13.5))') &
         (SUM(sintrans(:,ir)*sintrans(:,n))-&
         0.5D0*sintrans(0,ir)*sintrans(0,n),n=0,nqzmax)
    END DO

  END IF

  CLOSE(307)
  DEALLOCATE(rhoC2D,sintrans,besstrans,rhoform,form2D,qr)

RETURN

CONTAINS 
REAL(8) FUNCTION besszero(xlow,xup)

IMPLICIT REAL(8) (a-h,o-z)

REAL(8),INTENT(IN) :: xlow,xup
INTEGER, PARAMETER :: numiter=50
REAL(8), PARAMETER :: precis=1D-14

INTEGER :: iter
REAL(8) :: f1,f2,f3,x1,x2,x3

f1=bessel_j0(xlow)
f2=bessel_j0(xup)
x1=xlow
x2=xup

DO iter=1,numiter
  x3=0.5D0*(x1+x2)
  f3=bessel_j0(x3)
  IF(f3*f1<0D0) THEN
    x2=x3
    f2=f3
  ELSE IF(f2*f1<0D0) THEN
    x1=x3
    f1=f3
  ELSE IF(f3==0D0) THEN
    EXIT
  ELSE
    STOP 'case not allowed'
  END IF  
  IF(ABS(x2-x1)<precis) EXIT
END DO

besszero=x3
RETURN

END  FUNCTION besszero

REAL(8) FUNCTION overlap_1Dr(func1,func2)

REAL(8), INTENT(IN) :: func1(0:nr),func2(0:nr)

  overlap_1Dr = SUM(func1(:)*func2(:)*r)+dr*drezw*func1(0)*func2(0)

  RETURN

END FUNCTION overlap_1Dr

SUBROUTINE orthon_1Dr()

INTEGER :: n1,n2
REAL(8) :: factor

  DO n1=0,nqrmax  
    DO n2=0,n1-1
      factor = overlap_1Dr(besstrans(:,n1),besstrans(:,n2))
      besstrans(:,n1) = besstrans(:,n1)-factor*besstrans(:,n2)
    END DO
    factor = overlap_1Dr(besstrans(:,n1),besstrans(:,n1))
!    WRITE(307,'(a,i4,1pg13.5)') '# ortho:',n1,factor
    besstrans(:,n1) = besstrans(:,n1)/SQRT(factor)
  END DO

END SUBROUTINE orthon_1Dr

!---------------------------------------------------------------------------  
! DESCRIPTION: charge_formfactor
!> @brief
!!Functions which computes the charge formfactor for given formfactors
!!of proton and neutron distribution together with the nucleon
!!formfactors evaluated within this function routine.

!> @details
!!The present version omits the contribution from the magnetic formfactor 
!!of the nucleon.
!>
!> @param[in] q2
!> REAL(8), takes the squared radial momentum at which the 
!! fourier-transformation is to be evaluated.
!> @param[in] cmwid
!> REAL(8), parameter for center of mass correction.
!> @param[in] fop
!> REAL(8), the formfactor of the proton distribution at 'q1'
!> @param[in] fon
!> REAL(8),  the formfactor of the proton distribution at 'q1'
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION charge_formfactor(q2,cmwid,fop,fon)

USE fields
USE grid
USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)      :: q2
REAL(8), INTENT(IN OUT)  :: cmwid
REAL(8), INTENT(IN)      :: fop
REAL(8), INTENT(IN)      :: fon

REAL(8), PARAMETER ::  half=0.5D0

REAL(8) :: darfac,exphf,qa,schsen,schsep,schses,schsev

!      the parameters in data are for the sachs-formfactors:
!       sa..,sm.. are for the isoscalar electric formfactor
!       va..,vm.. are for the isovector electric formfactor
!       gp..,gm.. are for the magnetic formfactors of proton and neutron

!       isoscalar and isovector formfactors are recoupled to
!       the electric formfactors of proton and neutron.

!     the sachs-formfactors are then recoupled to the nucleon
!     vector formfactor (form1..) and the nucleon tensor
!     formfactor (form2..).

!     parameters for scalar and vector electric formfactor
!     (walther, priv.comm., sept. 86)

REAL(8),PARAMETER :: sa1=2.2907,sa2=-0.6777,sa3=-0.7923,sa4=0.1793,&
                     sm1=15.75, sm2=26.68,  sm3=41.01,  sm4=134.2
REAL(8),PARAMETER :: va1=0.3681,va2=1.2263,va3=-0.6316,va4=0.0372,&
                     vm1=5.00,  vm2=15.02,  vm3=44.08,  vm4=154.2

!     magnetic formfactor , simon et al n.p. a333 (1980) 381

REAL(8),PARAMETER :: gp1=0.694,gp2=0.719,gp3=-0.418,gp4=0.005,&
                     gm1=8.50, gm2=15.02,gm3=44.08, gm4=355.4

!     the anomanlous magnetic moments and the darwin factor 'd4m'

!REAL(8),PARAMETER :: apro0=0.74d0,apro11=0.0d0,apro12=0.0d0,bpro11=0.0d0
!REAL(8),PARAMETER :: aneu0=0.74d0,aneu11=0.0d0,aneu12=0.0d0,bneu11=0.0d0
!REAL(8),PARAMETER :: xmgne=0.04197d0,xmgpr=-.06131d0,d4m=.010989d0


!      isoscalar and isovector electrical formactors on 'schses'
!      and 'schsev', recoupled then to proton and neutron electric
!      formfactors  'schsep' and 'schsen'. proton and neutron
!      magnetic formfactor on 'schsmp' and 'schsmn'.

qa     = q2
schses = sa1/(1D0+qa/sm1)+sa2/(1D0+qa/sm2)+sa3/(1D0+qa/sm3) +sa4/(1D0+qa/sm4)
schsev = va1/(1D0+qa/vm1)+va2/(1D0+qa/vm2)+va3/(1D0+qa/vm3) +va4/(1D0+qa/vm4)
schsep = half*(schses+schsev)
schsen = half*(schses-schsev)
darfac = 1D0                      ! Darwin term -- ignored 
exphf  = EXP(cmwid*qa)            ! c.m. correction

charge_formfactor = (fop*schsep+fon*schsen)*exphf*darfac

RETURN
END FUNCTION charge_formfactor



END SUBROUTINE user_interface_1

SUBROUTINE user_interface_2()
  ! use hdf5
  use fields
  use hdf5_generic
  
  
  implicit none

  real(8), dimension(4,2) :: testArr

  call write_array(testArr)
  ! real(8), dimension(4,2) :: testArr

  ! call write_array_2d(testArr)
  
  ! call write_array(pot)

!   REAL(8),ALLOCATABLE :: pot(:,:,:)      !<Local HF potential.
! REAL(8),ALLOCATABLE :: coulpt(:,:)     !<Coulomb potential.
! REAL(8),ALLOCATABLE :: bmass(:,:,:)    !<Inverse mass.
! REAL(8),ALLOCATABLE :: potls(:,:,:,:)  !<Spin-orbit potential.
! REAL(8),ALLOCATABLE :: drbmas(:,:,:)   !<d/dr \c bmass.
! REAL(8),ALLOCATABLE :: dzbmas(:,:,:)   !<d/dz \c bmass.
! REAL(8),ALLOCATABLE :: potadd(:,:,:)   !<Additional potential to stabilize positive energy states.
! REAL(8),ALLOCATABLE :: rhosur(:,:,:)   !<Auxiliary field for potential wall.
! REAL(8),ALLOCATABLE :: pairpt(:,:,:)   !<Gap potential for dynamic BCS.
! REAL(8),ALLOCATABLE :: potdamp(:,:)    !<Damping potential for potential wall.

! REAL(8),ALLOCATABLE :: rho(:,:,:)      !<Density.
! REAL(8),ALLOCATABLE :: crho(:,:)       !<Coulomb density.
! REAL(8),ALLOCATABLE :: tau(:,:,:)      !<Kinetic energy density.
! REAL(8),ALLOCATABLE :: drhols(:,:,:)   !<Gradient of spin-orbit current.
! REAL(8),ALLOCATABLE :: rholsr(:,:,:)   !<Spin-orbit current in r-direction.
! REAL(8),ALLOCATABLE :: rholsz(:,:,:)   !<Spin-orbit current in z-direction.
! REAL(8),ALLOCATABLE :: d2rho(:,:,:)    !<Laplaceian on \c rho.
! REAL(8),ALLOCATABLE :: rhoti(:,:,:)    !<Density for constraint calculation.
! REAL(8),ALLOCATABLE :: paircr(:,:,:)   !<Pairing density.
  write(6,*) "Have I found the right file?"
  write(6,*) shape(rho)

  RETURN
END SUBROUTINE user_interface_2


SUBROUTINE user_interface_3()
!
! This routine initializes user defined pairing
!
  RETURN
END SUBROUTINE user_interface_3



SUBROUTINE user_interface_4()
!
! This routine points to user defined pairing
!
  RETURN
END SUBROUTINE user_interface_4


SUBROUTINE user_interface_5()
! user defined modifications of profile of constraint
  RETURN
END SUBROUTINE user_interface_5

SUBROUTINE user_interface_6()
! user defined alternative to mean-field Hamiltonian
RETURN
END SUBROUTINE user_interface_6

