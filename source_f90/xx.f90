!---------------------------------------------------------------------------  
! DESCRIPTION: Qvariance
!> @brief
!! Variance of quadrupole momentum, presently only for istsym>0
!--------------------------------------------------------------------------- 
SUBROUTINE Qvariance(DelQ2)

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid

IMPLICIT NONE

REAL(8), INTENT(OUT) :: DelQ2(-2:2)

INTEGER,PARAMETER :: intest=0,itstpr=0
REAL(8),PARAMETER :: fac20=SQRT(5D0/(16D0*pi))
REAL(8),PARAMETER :: fac21=SQRT(5D0/(16D0*pi))
REAL(8),PARAMETER :: fac22=SQRT(5D0/(16D0*pi))

REAL(8),ALLOCATABLE :: psi1(:,:,:),psin(:,:,:),psip(:,:,:)

REAL(8) :: deltai,deltaj,edifi,edifj,eferm,equasi,equasj,pij,&
          sumdir,tstovl,ui,uj,vi,vj,r2
INTEGER :: iz,ir,is,ist,jst,jstmax,jstmin,mist,mm,nparr,npist


IF(istsym.LE.0) THEN
  DelQ2=0D0
END IF

ALLOCATE(psi1(0:nz,0:nr,nsd),psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
psi1=0D0;psip=0D0;psin=0D0

!     1. double loop to accumulate Q_20 width
IF(intest == 1) WRITE(7,'(a)') ' Q_20:'
sumdir = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
! apply Q_20
  DO ir=0,nr
    r2     = r(ir)*+2 !(dr*(ir+sr))**2
    DO iz=0,nz
!      z     = ((iz)+sz-icm)*dz
      psip(iz,ir,:) = fac20*(2*z(iz)**2-r2)*psi(iz,ir,:,ist)
    END DO
  END DO
  
!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  sumdir = 4.0*vi*vi*wfovlp(psin,psin)*wtpair_v2(ist) +sumdir
  IF(intest == 1) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <sp|Q**2|sp>=',ist,tstovl
  DO jst=jstmin,jstmax
    IF(  npist == nparz(1,jst) .AND. mist == mpsi(jst) ) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      IF(isTsym.LE.0 .AND. mist == 0 .AND. mpsi(jst) == 0) THEN
        pij    = wfovoc(psi(0:nz,0:nr,1,jst),psin(0:nz,0:nr,2))  &
            + wfovoc(psi(0:nz,0:nr,2,jst),psin(0:nz,0:nr,1))
      ELSE 
        pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      END IF
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
      IF(intest == 1) WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij, pij*(vj*ui-uj*vi)
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
!               pij    = pij*(vj*ui-uj*vi)
!               pij    = pij*pij
!               sumexy = pij/(equasi+equasj) + sumexy
!               sumwxy  = pij + sumwxy
    END IF
  END DO
END DO
IF(isTsym>0) sumdir=sumdir/2D0
DelQ2(0) = sumdir/2.0

!     2. loop to accumulate Q_21-variance
IF(intest == 1) WRITE(7,'(a)') ' Q_21:'
sumdir = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
!  action of Q_21 operator
  DO ir=0,nr
!    r     = (dr*(ir+sr))
    DO iz=0,nz
!      z     = ((iz)+sz-icm)*dz
      psip(iz,ir,:) = fac21*z(iz)*r(ir)*psi(iz,ir,:,ist)
    END DO
  END DO

!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  IF(intest /= 0) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <Delta_zr**2>=',ist,tstovl
  sumdir = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir
  DO jst=jstmin,jstmax
    IF(npist == nparz(2,jst) .AND. mist+1 == mpsi(jst)) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      IF(intest == 1)   WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij*pij,  &
          (pij*(vj*ui-uj*vi))**2
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
    END IF
  END DO
END DO
IF(isTsym>0) sumdir=sumdir/2D0
DelQ2(1) = sumdir/2.0
DelQ2(-1) = sumdir/2.0



!     3. loop to accumulate Q_22-variance
IF(intest == 1) WRITE(7,'(a)') ' Q_22:'
sumdir = zero
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
  END IF
  
!  action of Q_21 operator
  DO ir=0,nr
    r2     = r(ir)**2 !(dr*(ir+sr))**2
    DO iz=0,nz
!      z     = ((iz)+sz-icm)*dz
      psip(iz,ir,:) = fac22*r2*psi(iz,ir,:,ist)
    END DO
  END DO

!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  IF(mist<0) ui=-ui
  deltai = deltaf(ist)*wtpair_uv(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
  IF(intest /= 0) WRITE(7,'(a,i3,1pg12.4)')  &
      ' test: ist, <Delta_zr**2>=',ist,tstovl
  sumdir = 2.0*vi*vi*tstovl*wtpair_v2(ist) +sumdir
  DO jst=jstmin,jstmax
    IF(npist == nparz(2,jst) .AND. mist+2 == mpsi(jst)) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      IF(mpsi(jst)<0) uj=-uj
      deltaj = deltaf(jst)*wtpair_uv(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
      IF(intest == 1)   WRITE(7,'(a,3(2i2,1x),2(1pg12.4))')  &
          ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',  &
          ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij*pij,  &
          (pij*(vj*ui-uj*vi))**2
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
      sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)*wtpair_v2(ist)*wtpair_v2(jst)
    END IF
  END DO
END DO
IF(isTsym>0) sumdir=sumdir/2D0
DelQ2(2) = sumdir/2.0
DelQ2(-2) = sumdir/2.0



DEALLOCATE(psi1,psin,psip)

RETURN
END SUBROUTINE Qvariance
