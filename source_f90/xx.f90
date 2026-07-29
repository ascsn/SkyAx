
! Code converted using TO_F90 by Alan Miller
! Date: 2026-07-16  Time: 09:56:38
 
!-----crnkam-----------------------------------------------------------

SUBROUTINE crnkam

!     and computes Inglis cranking mass for angular momentum

INCLUDE 'paramr.inc'
INCLUDE 'constr.inc'
INCLUDE 'observ.inc'
INCLUDE 'dynpar.inc'
INCLUDE 'fields.inc'
INCLUDE 'levels.inc'
INCLUDE 'wavefr.inc'
INCLUDE 'derivs.inc'

DATA itstpr/0/            ! switch to testprints of Inglis Cranking

!     buffered SQRT as internal statement function

sqrtb(x) = SQRT(MAX(x,smal))

!-----------------------------------------------------------------------


!     double loop to accumulate cranking mass for angular momentum

sumj2  = zero
sumj2e = zero
sumtst = 0.0
ampnms = 0.0
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    jstmin = 1
    eferm  = fer(1)
!R            efermo = ferolp
  ELSE
    jstmax = nst(1)+nst(2)
    jstmin = nst(1)+1
    eferm  = fer(2)
!R            efermo = feroln
  END IF
  
!     ang.mom. operator J+ on 'psip' and J- on 'psin'
  
!                                     first, the L+ and L-
  CALL nablaz_spin(psi(0,0,1,ist),psi1,nparz(1,ist),1,one)
  CALL nablar_spin(psi(0,0,1,ist),psi2,   mpsi(ist),1,one)
  
  DO is=1,nsd
    mm        = mist+is-1
    npist  = nparz(is,ist)
    nparr     = 1+MOD(mm,2)
!PK         call mulr(derr1(0,0,nparr),psi(0,0,is,ist),psi2(0,0,is),one,0)
!PK         call mulz(derz1(0,0,npist),psi(0,0,is,ist),psi1(0,0,is),one,0)
    ir = 0
    IF(mm == 1) THEN
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
        psip(iz,ir,is) = 0.0
!              psin(iz,ir,is) = (z+z)*psi2(iz,ir,is)
        psin(iz,ir,is) = -(z+z)*psi2(iz,ir,is)
      END DO
    ELSE
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
!              psip(iz,ir,is) = zero
        psip(iz,ir,is) = z*psi2(iz,ir,is)
!              psin(iz,ir,is) = psip(iz,ir,is)
        psin(iz,ir,is) = -psip(iz,ir,is)
      END DO
    END IF
    DO ir=1,nr
      r    = (ir+sr)*dr
      rinv = mm/r
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
        psip(iz,ir,is) = z*(psi2(iz,ir,is)-rinv *psi(iz,ir,is,ist))  &
            -r*psi1(iz,ir,is)
!              psin(iz,ir,is) = z*(wrk2(iz,ir)+rinv*psi(iz,ir,is,ist))
!     &                        -r*wrk4(iz,ir)
        psin(iz,ir,is) = -z*(psi2(iz,ir,is)+rinv *psi(iz,ir,is,ist))  &
            +r*psi1(iz,ir,is)
      END DO
    END DO
  END DO
!                                          the S+ and S-
  DO ir=0,nr
    DO iz=0,nz
      psip(iz,ir,1) = psip(iz,ir,1)+psi(iz,ir,2,ist)
      psin(iz,ir,2) = psin(iz,ir,2)+psi(iz,ir,1,ist)
    END DO
  END DO
!                                   test       the S+ and S-
!t        do ir=0,nr
!t        do iz=0,nz
!t          psip(iz,ir,1) = -psi(iz,ir,2,ist)
!t          psip(iz,ir,2) = 0.0
!t          psin(iz,ir,1) = 0.0
!t          psin(iz,ir,2) = psi(iz,ir,1,ist)
!t        enddo
!t        enddo
!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  deltai = deltaf(ist)*wtpair(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip,psip)
!test         write(7,'(a,i3,1pg12.4)') ' test: ist, <J**2>=',ist,tstovl
  sumtst = tstovl*vi*vi + sumtst
  DO jst=jstmin,jstmax
    IF((npist == nparz(2,jst) .AND.  &
          (mist+1 == mpsi(jst) .OR. mist-1 == mpsi(jst)) ) .OR.  &
          (npist == nparz(1,jst) .AND.  &
          (mist == 0 .AND. mpsi(jst) == 0))            ) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      deltaj = deltaf(jst)*wtpair(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      IF(mist == 0 .AND. mpsi(jst) == 0) THEN
!test            if((ist.eq.2 .and. jst.eq.5) .or.
!test     &         (ist.eq.5 .and. jst.eq.2))
!test     &         write(7,'(a,2(2x,2(1pg12.4)))')
!test     &        ' 2&5: upper,lower=',
!test     &        wfovoc(psi(0,0,1,jst),psin(0,0,2)),
!test     &        wfovoc(psi(0,0,2,jst),psin(0,0,1)),
!test     &        wfovoc(psi(0,0,1,jst),psi(0,0,1,ist)),
!test     &        wfovoc(psi(0,0,2,jst),psi(0,0,2,ist))
!            do ir=0,nr
!            do iz=0,nz
!              wrk2(iz,ir)   = psin(iz,ir,1)
!              psin(iz,ir,1) = psin(iz,ir,2)
!              psin(iz,ir,2) = wrk2(iz,ir)
!              psin(iz,ir,2) = -psav
!            enddo
!            enddo
        pij    = +wfovoc(psi(0,0,2,jst),psin(0,0,1))  &
            -wfovoc(psi(0,0,1,jst),psin(0,0,2))
        pij    = pij*(vj*ui-uj*vi)
!            write(7,'(a,3(2i2,1x),2(1pg12.4))')
!     &     ' 0:ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
!     &       ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
!     &       pij*(vj*ui-uj*vi)
      ELSE IF(mist+1 == mpsi(jst)) THEN
        pij    = wfovlp(psi(0,0,1,jst),psip)
        pij    = pij*(vj*ui-uj*vi)
!            write(7,'(a,3(2i2,1x),2(1pg12.4))')
!     &     ' +:ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
!     &       ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
!     &       pij*(vj*ui-uj*vi)
      ELSE IF(mist-1 == mpsi(jst)) THEN
        pij    = wfovlp(psi(0,0,1,jst),psin)
        pij    = pij*(vj*ui-uj*vi)
!            write(7,'(a,3(2i2,1x),2(1pg12.4))')
!     &     ' -:ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
!     &       ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
!     &       pij*(vj*ui-uj*vi)
      END IF
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
!PK               pij    = pij*(vj*ui-uj*vi)
      pij    = pij*pij
      sumj2e = pij/(equasi+equasj) + sumj2e
      sumj2  = pij + sumj2
      ampnms = pij/(equasi+equasj)**2 *(ui*ui-vi*vi+uj*uj-vj*vj)  &
          + ampnms
    END IF
  END DO
END DO
ammass = sumj2e/2.0
amwid  = sumj2/2.0
ampnms = ampnms/2.0/ammass**2
WRITE(7,'(a,1pg12.4)') ' test <J**2>=',sumtst

RETURN
END SUBROUTINE crnkam


