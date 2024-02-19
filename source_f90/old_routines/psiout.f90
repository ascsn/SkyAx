
!---------------------------------------------------------------------------  
! DESCRIPTION: psiout
!> @brief
!!Subroutine \c psiout assigns Nillson quantum numbers to the s-p levels, 
!!computes the complete space parity \c ntpar, writes single-particle wave
!!functions, energies, quantum numbers and occupation numbers to the output 
!!file SPOUT and writes the local s-p potential to the output file 'sppoten'.
!!Used as interface to the Dubna SRPA code.
!>
!> INTEGER, array, 
!--------------------------------------------------------------------------- 

SUBROUTINE psiout()

USE paramr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT NONE

INTEGER :: iord2(npsi(2)) !iorder (i) - numbers of the ordered levels 
!(protons from 1 to nst(1) and neutrons from 1 to to nst(2)). For example,
!if there are 20 proton and 30 neutron levels, then iorder changes in the 
!intervals 1-20 for protons and 1-30 for neutrons. At the same time 'ist' 
!changes from 1 to 50.
REAL(8) :: spen2(npsi(2))

!INTEGER,INTENT(IN) :: ireg(2)

REAL(8) :: ak(npsi(2))

INTEGER, PARAMETER :: nsemax=25!max. number of s-p states of EVEN-parity and 
!of fixed isospin and  magnetic number Lambda
INTEGER, PARAMETER :: nlemax=8!max. number of Lambda's for EVEN-parity states 
!of a fixed isospin
INTEGER, PARAMETER :: nsomax=30!max. number of s-p states of ODD-parity and 
!of fixed isospin and  magnetic number Lambda
INTEGER, PARAMETER :: nlomax=10!max. number of Lambda's for ODD-parity states 
!of a fixed isospin

INTEGER :: niprev(nsemax,nlemax),niprod(nsomax,nlomax)
INTEGER :: niobev(nsemax,nlemax),niobod(nsomax,nlomax)
INTEGER :: icount(2,10),nqn(npsi(2))!Nillson quantum numbers of the state 'ist'

INTEGER :: i,i1,i2,i3,ifin,init,ipar,ir,is,ist,iz,lam,np
REAL(8) :: rdist,zdist

DATA niprev/0,220,211,200,440,431,420,411,400,   &
    660,651,640,631,620,611,600, 880,871,860,851,840,831,820,811,800,  &
    211,202,431,422,411,402,                        &
    651,642,631,622,611,602, 871,862,851,842,831,822,811,802,  &
    0,0,0,0,0, 202,422,413,402,642,633,622,613,602,           &
    862,853,842,833,822,813,802, 0,0,0,0,0,0,0,0,0,  &
    413,404,633,624,613,604,                     &
    853,844,833,824,813,804, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0, 404,624,615,604,844,835,824,815,804,      &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,  &
    615,606,835,826,815,806,                       &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,  &
    606,826,817,806,                               &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0, 817,808,                                    &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0/

DATA niobev/0,200,211,220,400,411,420,431,440,   &
    600,611,620,631,640,651,660, 800,811,820,831,840,851,860,871,880,  &
    202,211,402,411,422,431,                       &
    602,611,622,631,642,651, 802,811,822,831,842,851,862,871,  &
    0,0,0,0,0, 202,402,413,422,602,613,622,633,642,            &
    802,813,822,833,842,853,862, 0,0,0,0,0,0,0,0,0,  &
    404,413,604,613,624,633,                        &
    804,813,824,833,844,853, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0, 404,604,615,624,804,815,824,835,844,            &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,  &
    606,615,806,815,826,835,                        &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,  &
    606,806,817,826,                                &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0, 808,817,                                       &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0/

!     prolate, odd parity:
!    line 1     K=1/2-
!    lines 2/3  K=3/2-

DATA niprod/110,101,330,321,310,301,           &
    550,541,530,521,510,501, 770,761,750,741,730,721,710,701,  &
    990,981,970,961,950,941,930,921,910,901,  &
    101,321,312,301,541,532,521,512,501,          &
    761,752,741,732,721,712,701, 981,972,961,952,941,932,921,912,901,  &
    0,0,0,0,0, 312,303,532,523,512,503,                      &
    752,743,732,723,712,703, 972,963,952,943,932,923,912,903,  &
    0,0,0,0,0,0,0,0,0,0, 303,523,514,503,743,734,723,714,703,          &
    963,954,943,934,923,914,903, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0, 514,505,734,725,714,705,                      &
    954,945,934,925,914,905, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0,0,0, 505,725,716,705,945,936,925,916,905,           &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0, 716,707,936,927,916,907,                       &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0, 707,927,918,907,                               &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0, 918,909,                                       &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0,0,0, 909,                                           &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0,0,0,0/

!     oblate, odd parity:
DATA niobod/101,110,301,310,321,330,           &
    501,510,521,530,541,550, 701,710,721,730,741,750,761,770,  &
    901,910,921,930,941,950,961,970,981,990,  &
    101,301,312,321,501,512,521,532,541,          &
    701,712,721,732,741,752,761, 901,912,921,932,941,952,961,972,981,  &
    0,0,0,0,0, 303,312,503,512,523,532,                      &
    703,712,723,732,743,752, 903,912,923,932,943,952,963,972,  &
    0,0,0,0,0,0,0,0,0,0, 303,503,514,523,703,714,723,734,743,          &
    903,914,923,934,943,954,963, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0, 505,514,705,714,725,734,                      &
    905,914,925,934,945,954, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0,0,0, 505,705,716,725,905,916,925,936,945,           &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0, 707,716,907,916,927,936,                       &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0, 707,907,918,927,                               &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0, 909,918,                                       &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0,0,0, 909,                                           &
    0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,  &
    0,0,0,0,0,0,0,0,0/
!------write single-particle potential:
OPEN(81,STATUS='unknown',FORM='formatted', FILE='sppoten.'//outnam)
DO iz=0,nz
  zdist=dz*iz
  DO ir=0,nr
    rdist=dr*ir
    WRITE(81,'(3e13.5)') zdist,rdist,pot(iz,ir,1)
  END DO
END DO
!------end write single-particle potential
OPEN(80,STATUS='unknown',FORM='unformatted', FILE='nillson.'//outnam)
!IF (ireg(1) > 0) THEN
  WRITE(80) atnum,apnum,annum,beta2,beta2p,beta2n,  &
      beta4,beta4p,beta4n,nst(1),nst(2),nzd,nrd,dz,dr
!     &  gap(1),gap(2),fer(1),fer(2)
!END IF
DO  np=1,2
!------ protons (np=1) and neutrons (np=2)
  IF(np == 1) THEN
    WRITE(7,'(10x,a)') 'PROTONS:'
  ELSE
    WRITE(7,'(10x,a)') 'NEUTRONS:'
  END IF
  WRITE(7,'(10x,a)') 'Reinhards quasiordering'
  init=1+nst(1)*(np-1)
  ifin=nst(1)+nst(2)*(np-1)
  DO i=init,ifin
    iord2(i)=iorder(i)
    ist=iord2(i)+nst(1)*(np-1)
    spen2(i)=spenrg(ist)
    WRITE(7,'(2x,2i3,2x,e11.5)') i,iord2(i),spen2(i)
  END DO
!----------------------------------------
!        Ordering of Reinhard  is approximate.
!        So, before the asignment, we repeat the ordering.
!----------------------------------------
!********************
  CALL order(init,ifin,iord2,spen2)
!********************
  WRITE(7,'(2x,a)') '  i  ist [Nn_zl]  k   occ    energy'
!--------assignment of asypt. Nillson quantum numbers:

! 1st index means parity but not isospin
! 2nd index means number of K for a given parity.
icount(1:2,1:10)=0
  DO  i=init,ifin
    ist=iord2(i)+nst(1)*(np-1)
    lam=mpsi(ist)
    ipar=npar(ist)*nrpar(ist)
    IF(ipar == 1) THEN
      i1=1
      i2=lam+1
!   i2=mpsi+1=(mpsi+1/2)+1/2=K+1/2
!   K   i2
!   1/2  1
!   3/2  2
!   ...  ...
      icount(i1,i2)=icount(i1,i2)+1
      i3=icount(i1,i2)
!   i3 number of [Nn_z|Lambda] for a given K^{\pi}
      IF(beta2 >= 0.0) THEN
        nqn(ist)=niprev(i3,i2)
      ELSE
        nqn(ist)=niobev(i3,i2)
      END IF
    ELSE
      i1=2
      i2=lam+1
      icount(i1,i2)=icount(i1,i2)+1
      i3=icount(i1,i2)
      IF(beta2 >= 0.0) THEN
        nqn(ist)=niprod(i3,i2)
      ELSE
        nqn(ist)=niobod(i3,i2)
      END IF
    END IF
    ak(ist)=mpsi(ist)+0.5
    WRITE(7,'(2x,3(i3,2x),f5.1,f8.5,2x,e11.5)') i,  &
        ist,nqn(ist),ak(ist),vsqf(ist),spenrg(ist)
    IF(.NOT.(np == 1.AND.vsqf(ist) > 0D0)) CYCLE
!----write w.f.
    WRITE(80) (((psi(iz,ir,is,ist), iz=0,nz),ir=0,nr),is=1,nsd)
    WRITE(80) spenrg(ist),mpsi(ist),npar(ist),nrpar(ist),  &
        ipar,icharg(ist),vsqf(ist),iorder(ist), nqn(ist),spn(ist)
!----end write w.f.
  END DO
END DO
CLOSE(80)
CLOSE(81)
RETURN
END SUBROUTINE psiout

!---------------------------------------------------------------------------  
! DESCRIPTION: order
!> @brief
!!Orders \c iord2 and \c spen2 with respect to \c spen2 in the range of 
!!\c imin to \c imax.
!>
!> INTEGER, array, 
!--------------------------------------------------------------------------- 

SUBROUTINE order(imin,imax,iord2,spen2)

USE paramr
USE levels, ONLY: npsi

IMPLICIT NONE

INTEGER,INTENT(IN OUT) :: iord2(npsi(2)),imin,imax
REAL(8),INTENT(IN OUT) :: spen2(npsi(2))

INTEGER :: i,ip,j,j1,k,max1,min1
REAL(8) :: x1


IF(imin >= imax) RETURN
MAX1=imax-1
MIN1=imin

DO WHILE(MAX1 > MIN1)
  ip=1
  j1=ip
  j=j1
  DO  i=MIN1,MAX1
    IF(spen2(i) <= spen2(i+1)) CYCLE 
    !swap spen2
    x1=spen2(i)
    spen2(i)=spen2(i+1)
    spen2(i+1)=x1
    !swap ird2
    k=iord2(i)
    iord2(i)=iord2(i+1)
    iord2(i+1)=k
    IF(ip == 1) THEN
      j=i-1
      ip=0
    END IF
    j1=i
  END DO
  MAX1=j1
  MIN1=MAX(j,imin)
END DO

RETURN
END SUBROUTINE order
