!------------------------------------------------------------------------------
! MODULE: steputil
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routines for in and output of wave functions and helper 
!!routines for printout.
!------------------------------------------------------------------------------
MODULE util

IMPLICIT NONE

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: write_skyax_header
!> @brief
!! Writes header for main output.
!--------------------------------------------------------------------------- 

SUBROUTINE write_skyax_header()
IMPLICIT NONE

WRITE(7,'(10x,a)')'============================================='
WRITE(7,'(10x,a)')'|| SSSSSSS K     K Y     Y    A    X     X ||'
WRITE(7,'(10x,a)')'|| S       K   K    Y   Y    A A    X   X  ||'
WRITE(7,'(10x,a)')'|| S       K K       Y Y    A   A    X X   ||'
WRITE(7,'(10x,a)')'|| SSSSSSS K          Y    A     A    X    ||'
WRITE(7,'(10x,a)')'||       S K K        Y    AAAAAAA   X X   ||'
WRITE(7,'(10x,a)')'||       S K   K      Y    A     A  X   X  ||'
WRITE(7,'(10x,a)')'|| SSSSSSS K     K    Y    A     A X     X ||'
WRITE(7,'(10x,a)')'============================================='
WRITE(7,'(10x,a)')'     the axial Skyrme-Hartree-Fock Code'
WRITE(7,'(10x,a)')'            >>> version x.x.x <<<       '

WRITE(6,'(10x,a)')'============================================='
WRITE(6,'(10x,a)')'|| SSSSSSS K     K Y     Y    A    X     X ||'
WRITE(6,'(10x,a)')'|| S       K   K    Y   Y    A A    X   X  ||'
WRITE(6,'(10x,a)')'|| S       K K       Y Y    A   A    X X   ||'
WRITE(6,'(10x,a)')'|| SSSSSSS K          Y    A     A    X    ||'
WRITE(6,'(10x,a)')'||       S K K        Y    AAAAAAA   X X   ||'
WRITE(6,'(10x,a)')'||       S K   K      Y    A     A  X   X  ||'
WRITE(6,'(10x,a)')'|| SSSSSSS K     K    Y    A     A X     X ||'
WRITE(6,'(10x,a)')'============================================='
WRITE(6,'(10x,a)')'     the axial Skyrme-Hartree-Fock Code'
WRITE(6,'(10x,a)')'            >>> version x.x.x <<<       '
END SUBROUTINE

!---------------------------------------------------------------------------  
! DESCRIPTION: plot_2d
!> @brief
!! Plots densities in the file 'plot_2D.???'
!--------------------------------------------------------------------------- 

SUBROUTINE plot_2d()

USE paramr
USE dynpar
USE fields
USE observ
USE levels, ONLY : anumb,npart,nst
USE grid

IMPLICIT NONE

LOGICAL :: tfirst=.true.

INTEGER :: ir,iz
INTEGER,SAVE :: cit=1

IF(tfirst) THEN
  OPEN(77,FILE='plot_2D.'//outnam)
  WRITE(77,'(a,a)')'# ',title(1:78)
  WRITE(77,'(a,2i5,a,f10.5)') '# Nr. proton/neutron=',npart, ' total Nr. Nucleons=',anumb
  WRITE(77,'(a,2i5)') '# Nr. of states=',nst
  WRITE(77,'(a,2(i5,X),2(f6.3,X),L5)') '# nr,nz,dr,dz,tsym= ',nr,nz,dr,dz,tsym
  tfirst = .false.
ELSE
  OPEN(77,FILE='plot_2D.'//outnam,POSITION='APPEND')
  WRITE(77,'(1x)')
END IF
WRITE(77,'(a,i5)')'# number ',cit
WRITE(77,'(a,f9.5,2f7.3)')'# deformations beta2,beta3,beta4=',beta2,beta3,beta4
WRITE(77,'(a,2f13.5)') '# energy, rms-radius= ',energy,rms
DO ir=0,nr
  DO iz=0,nz
    WRITE(77,'(2f6.1,9(1pg12.4))') r(ir),z(iz)-icm*dz,  &
        rho(iz,ir,1),rho(iz,ir,2),tau(iz,ir,1),tau(iz,ir,2),&
        pot(iz,ir,1),pot(iz,ir,2),paircr(iz,ir,1),paircr(iz,ir,2),coulpt(iz,ir)
  END DO
  WRITE(77,'(1x)')
END DO

CLOSE(77)
cit=cit+1
RETURN
END SUBROUTINE plot_2d

!---------------------------------------------------------------------------  
! DESCRIPTION: polint
!> @brief
!!Polynomial interpolation.
!>
!> @param[in] xa
!> REAL(8), array, takes x values.
!> @param[in] ya
!> REAL(8), array, takes corresponding y values.
!> @param[in] n
!> INTEGER, takes order of polynomial for interpolation.
!> @param[in] x
!> REAL(8), takes x value at which interpolated function is wanted.
!> @param[out] y
!> REAL(8), returns interpolated function value.
!> @param[out] dy
!> REAL(8), returns error estimate for resulting y.
!--------------------------------------------------------------------------- 

SUBROUTINE polint(xa,ya,n,x,y,dy)

IMPLICIT NONE

REAL(8), INTENT(IN)             :: xa(n)
REAL(8), INTENT(IN)             :: ya(n)
INTEGER, INTENT(IN)             :: n
REAL(8), INTENT(IN)             :: x
REAL(8), INTENT(OUT)            :: y
REAL(8), INTENT(OUT)            :: dy


INTEGER, PARAMETER :: nmax=10

INTEGER :: i,m,ns
REAL(8) :: den,dif,dift,ho,hp,w,c(nmax),d(nmax)

ns=1
dif=ABS(x-xa(1))
DO  i=1,n
  dift=ABS(x-xa(i))
  IF (dift < dif) THEN
    ns=i
    dif=dift
  END IF
  c(i)=ya(i)
  d(i)=ya(i)
END DO
y=ya(ns)
ns=ns-1
DO  m=1,n-1
  DO  i=1,n-m
    ho=xa(i)-x
    hp=xa(i+m)-x
    w=c(i+1)-d(i)
    den=ho-hp
    IF(den == 0.d0)  STOP 'failure in polint'
    den=w/den
    d(i)=hp*den
    c(i)=ho*den
  END DO
  IF (2*ns < n-m)THEN
    dy=c(ns+1)
  ELSE
    dy=d(ns)
    ns=ns-1
  END IF
  y=y+dy
END DO
RETURN
END SUBROUTINE polint

!---------------------------------------------------------------------------  
! DESCRIPTION: protfold
!> @brief
!!Folding of given density \c rhoin with Gaussian of width \c rcut1. Result is 
!!returned on \c rhofold. 
!> @details
!!Uses Taylor expansion of Gaussian in momentum space.
!>
!> @param[in] rhoin
!> REAL(8), array, takes the input density.
!> @param[out] rhofold
!> REAL(8), array, returns the folded density.
!--------------------------------------------------------------------------- 

SUBROUTINE protfold(rhoin,rhofold)

USE paramr
USE grid

IMPLICIT NONE

REAL(8), INTENT(IN)                         :: rhoin(0:nz,0:nr)
REAL(8), INTENT(OUT)                        :: rhofold(0:nz,0:nr)

REAL(8), PARAMETER :: rcut1=0.8D0  !0.55D0
REAL(8), PARAMETER :: rcut2=rcut1*rcut1/6D0

REAL(8),ALLOCATABLE :: wrk1(:,:)
LOGICAL,PARAMETER :: ttestnegative=.TRUE.

ALLOCATE(wrk1(0:nz,0:nr))
wrk1=0D0

!     first folding term, stored temporarily on 'rhofold'
CALL laplace_scalar(rhoin(0:nz,0:nr),rhofold(0:nz,0:nr), 1, 0, 1, rcut2)

!     second folding term, stored on 'wrk1'
CALL laplace_scalar(rhofold(0:nz,0:nr), wrk1(0:nz,0:nr), 1, 0, 1, rcut2/2D0)

!     compose
rhofold(0:nz,0:nr) = rhoin(0:nz,0:nr)+rhofold(0:nz,0:nr)+wrk1(0:nz,0:nr)
DEALLOCATE(wrk1)

IF(ttestnegative) WRITE(6,*) 'PROTFOLD negative density:',MINVAL(rhofold(0:nz,0:nr)),MAXVAL(rhofold(0:nz,0:nr))

RETURN
END SUBROUTINE protfold

SUBROUTINE stimer(iaction,text)

!     Initializes timer and retrieves CPU times                                                                                                                      
!       iaction  = 1  for initial call                                                                                                                               
!                  2  for subsequent call                                                                                                                            
!     The relative time is taken between step 2 and 1.                                                                                                               

!USE dynpar

INTEGER, INTENT(IN) :: iaction
CHARACTER(*), INTENT(IN) :: text

INTEGER, SAVE :: itimeold=0,itimenew=0,itimeinit

!-----------------------------------------------------------------                                                                                                   


IF(iaction == 1) THEN
  CALL system_clock(itimeold)
  itimeinit = itimeold
ELSE IF(iaction == 2) THEN
  CALL system_clock(itimenew)
  WRITE(6,'(a,2(1pg13.5))') text//' time=',&
   (itimenew-itimeold)*1D-4,(itimenew-itimeinit)*1D-4
  itimeold = itimenew
ELSE
  STOP 'SYS.TIMER: this IACTION is not valid'
END IF

RETURN
END SUBROUTINE stimer


 




END MODULE

