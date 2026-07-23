!-----lagpar
!
      subroutine lagpar(hfbpsi,symcond,istmin,istmax,iter,
     &                  devmf2,devh2,devgap,qq,qh,qhdig,tqq)

!
!     Subroutine to compute the action of the HFB Hamiltonian.
!     Computes the matrix of HFB Lagrange parameters and
!     with that the action of the HFB Hamiltonian on the 
!     s.p. wavfunctions for both types of nucleons separately.
!
!      hfbpsi  = set of s.p. wavefunctions after action of
!                the HFB Hamiltonian
!      symcond = accumulator for symmetry condition
!      istmin  = minimal state nr.
!      istmax  = maximal state nr.
!      iter    = actual number of iteration
!      devmf2,devgap,qq,qh,qhdig  = accumulators for variances
!
      include 'paramr.inc'
      include 'constr.inc'
!      include 'observ.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'resden.inc'
      include 'crankin.inc'

      REAL(8),INTENT(OUT) :: hfbpsi(0:nzd,0:nrd,nsd,npsid)
      REAL(8),INTENT(IN OUT) :: symcond  
      REAL(8),INTENT(IN OUT) :: devmf2,devh2,devgap,qq,qh,qhdig
      INTEGER,INTENT(IN) :: iter,istmin,istmax 
      LOGICAL,INTENT(IN) :: tqq

      REAL(8), ALLOCATABLE :: hlagr(:,:)
      REAL(8), ALLOCATABLE :: symatrix(:,:)
      REAL(8),PARAMETER :: xsmall=1D-20

      LOGICAL,PARAMETER :: tbcs=.FALSE.   ! testswitch
      LOGICAL,PARAMETER :: ttest=.FALSE.   ! testswitch
      LOGICAL,PARAMETER :: tvarh=.FALSE.   ! full variance (not yet proper)
!      LOGICAL,PARAMETER :: tHFBvariat=.FALSE.   ! variational HFB Hamiltonian
!
!      DO iq=1,2
!        IF(iq==1) THEN
!          istmin = 1
!          istmax = nst(1)
!        ELSE
!          istmin = nst(1)+1
!          istmax = nst(1)+nst(2)
!        END IF
      IF(istmin.LE.nst(1)) THEN
       iq = 1
      ELSE
        iq = 2
      END IF

!     compute and store action of HFB Hamiltonian
      DO jst=istmin,istmax
        call hpsi(psi(0,0,1,jst),hfbpsi(0,0,1,jst),jst)
!          hfbpsi(0:nz,0:nr,1:2,jst) = hfbpsi(0:nz,0:nr,1:2,jst)
!     &                               +psi(0:nz,0:nr,1:2,jst)*spenrg(jst)
        IF(jst==istcblock(iq) .OR. jst==istblo(iq)) THEN    !corrected
          v2a=1D0
          uva=0D0
!          WRITE(7,'(a,4i5,4f8.4)')
!     &        'jump=',iq,jst,istblo(iq),istcblock(iq),v2a,uva,
!     &         wfovlp(psi(0,0,1,jst),hfbpsi(0,0,1,jst)),spenrg(jst)
        ELSE
          v2a = vsqf(jst)
          uva = SQRT(MAX(v2a-v2a*v2a,xsmall))
        END IF
        IF(tvarh) THEN
          psip(0:nz,0:nr,1) = 
     &                  +pairpt(0:nz,0:nr,iq)*psi(0:nz,0:nr,1,jst)
          psip(0:nz,0:nr,2) = 
     &                  +pairpt(0:nz,0:nr,iq)*psi(0:nz,0:nr,2,jst)
          accsub=0D0
          vvj=vsqf(jst)
          uvj=SQRT(vvj-vvj*vvj)
          uuj=1D0-vvj
          DO kst=istmin,istmax
            IF(kst==jst) CYCLE
            vvk=vsqf(kst)
            uvk=SQRT(vvk-vvk*vvk)
            uuk=1D0-vvk
            if (nparz(1,jst).eq.nparz(1,kst)) then
              hovlp=wfovlp(psi(0,0,1,kst),hfbpsi(0,0,1,jst))
              accsub = accsub 
     &           +( hovlp*(vvj*vvk-uvj*uvk)
     &             +wfovlp(psi(0,0,1,kst),psip)**2*(uuj*uuk-vvj*vvk))**2
     &          -hovlp**2*vvj
            endif
          END DO
          accsub=wate(jst)*accsub
          devh2  =vsqf(jst)*wate(jst)*spevar(jst)**2 + accsub 
     &            + devh2
        END IF
        IF(tbcs) THEN
          v2a=1D0
          uva=0D0
          IF(ttest) 
     &      WRITE(7,'(a,i4,1pg13.5)') 'ist,spenrg=',jst,spenrg(jst)
          hfbpsi(0:nz,0:nr,is,jst) = hfbpsi(0:nz,0:nr,is,jst)
     &        -psi(0:nz,0:nr,is,jst)*spenrg(jst)
          IF(ttest) THEN
            WRITE(7,'(a,i4)') 'ortho-normality before: jst=',jst
            WRITE(7,'(100(1pg13.5))') 
     &        (wfovlp(psi(0,0,1,ist),psi(0,0,1,jst)),ist=istmin,istmax)
          END IF
        ELSE
          IF(.NOT.tpairspace) uva=uva*wtpair(jst)
          IF(tHFBvariat) THEN
            fackincut = -(1D0-wtpair(jst))*deltaf(jst)/(-0.1D0*cutfac)
            v2avar = v2a + uva*fackincut
          ELSE
            v2avar = v2a
          END IF
          DO is=1,2
            hfbpsi(0:nz,0:nr,is,jst) = hfbpsi(0:nz,0:nr,is,jst)*v2avar
     &                  +pairpt(0:nz,0:nr,iq)*psi(0:nz,0:nr,is,jst)*uva
          END DO
!          IF(jst==istcblock(iq) .OR. jst==istblo(iq)) 
!     &      WRITE(7,'(a,i4,4f8.4)') 'after jump=',jst,
!     &        wfovlp(psi(0,0,1,jst),hfbpsi(0,0,1,jst)),spenrg(jst)
        END IF
      END DO


!
!        accumulation of variances, <QQ>, <Qh>, devmf, devgap
!
      DO jst=istmin,istmax
         IF(tbcs) devmf2  = vsqf(jst)*wate(jst)*spevar(jst)**2 + devmf2
         vv2     = vsqf(jst)
         IF(tpairspace) THEN
           devgap  = wate(jst)*
     &            ((spenrg(jst)-fer(iq))*2.0*sqrt(vv2-vv2*vv2)
     &             +deltaf(jst)*(one-vv2-vv2) )
     &            + devgap
         ELSE
           devgap  = wate(jst)*
     &            ((spenrg(jst)-fer(iq))*2.0*sqrt(vv2-vv2*vv2)
     &             +deltaf(jst)*wtpair(jst)*(one-vv2-vv2) )
     &            + devgap
         END IF
         if(tqq)  then
           do is=1,nsd
             psip(0:nz,0:nr,is) = qfac(0:nz,0:nr)*psi(0:nz,0:nr,is,jst)
           enddo
c            write(6,'(a,i5)') ' CALL PROJEC at jst=',jst
c            call projec(psip,npar(jst),nlo-1,nhi-nlo+1,offnor)
           qq = vsqf(jst)*wate(jst)*wfovlp(psip,psip) + qq
           qh = vsqf(jst)*wate(jst)*wfovlp(psip,psin) + qh
           va = vsqf(jst)
           ua = sqrt(one-va)
           va = sqrt(va)
           qhdig = (ua*ua-va*va)*2.0*ua*va*deltaf(jst) + qhdig
           do jstp=istmin,istmax  !nlo,nhi
             IF(jst.LE.0 .OR. jstp.LE.0) THEN
!               WRITE(6,'(a,6i10)') 'nlo,nhi,jst,jstp=',nlo,nhi,jst,jstp
               CYCLE
             END IF
             if(nparz(1,jst).eq.nparz(1,jstp)) then
               hab    = wfovlp(psi(0,0,1,jstp),psin)
               qab    = wfovlp(psi(0,0,1,jstp),psip)
               vb     = vsqf(jstp)
               ub     = sqrt(one-vb)
               vb     = sqrt(vb)
               offwgt = va*vb*(ua*ub-va*vb)*wate(jst)
               qh     = offwgt*hab*qab + qh
               qq     = offwgt*qab*qab + qq
             endif
           enddo
         endif
      END DO

!     compute matrix of Lagrange parameters
      ALLOCATE(hlagr(istmin:istmax,istmin:istmax))
      ALLOCATE(symatrix(istmin:istmax,istmin:istmax))
      hlagr = 0D0
      symatrix = 0D0
      DO ist=istmin,istmax
      mist = mpsi(ist)
      npzist = nparz(1,ist)
      DO jst=istmin,ist
        IF(mpsi(jst)==mist .AND. nparz(1,jst)==npzist) THEN
          hlagr(ist,jst) = wfovlp(psi(0,0,1,ist),hfbpsi(0,0,1,jst))
          symatrix(ist,jst) = hlagr(ist,jst)
          IF(ist.NE.jst) THEN
            IF(ist==istblo(iq) .OR. jst==istblo(iq) .OR.                    !corrected
     &         ist==istcblock(iq) .OR. jst==istcblock(iq)) THEN
              hlagr(jst,ist)=0D0
              hlagr(ist,jst)=0D0
              symatrix(ist,jst)=0D0
              symatrix(jst,ist)=0D0
            ELSE
              hlagr(jst,ist) = wfovlp(psi(0,0,1,jst),hfbpsi(0,0,1,ist))
              symatrix(ist,jst) = hlagr(ist,jst)-hlagr(jst,ist)
              symatrix(jst,ist) = symatrix(ist,jst)
              hlagr(ist,jst) = 0.5D0*(hlagr(ist,jst)+hlagr(jst,ist))
              hlagr(jst,ist) = hlagr(ist,jst) 
            END IF
          ELSE
!            IF(jst==istcblock(iq)) WRITE(7,'(a,i4,2(1pg13.5))') 
!     &       'blocked s.p.energy=',jst,hlagr(jst,jst),spenrg(jst)
          END IF
        ELSE
          hlagr(ist,jst) = 0D0
        END IF
      END DO
      END DO

!     accumulate symmetry condition
      DO ist=istmin,istmax
      DO jst=istmin,ist-1
        symcond = symcond + symatrix(ist,jst)**2
      END DO
      END DO

!     subtract Lagrange matrix from 'hfbpsi'
      IF(tbcs) hlagr=0D0
!      IF(istmin==9)  hlagr(9:16,9:16)=0D0
!      IF(istmin==1) OPEN(933,file='hlagr')
!      WRITE(933,'(a)') '# ist  m(ist)  jst  m(jst)  hlagr'
      DO ist=istmin,istmax
!        IF(ist==istcblock(iq)) WRITE(7,'(a,50(1pg13.5))') 'HLAGR:',
!     &    (hlagr(jst,ist),jst=9,16)
        mist = mpsi(ist)
        npzist = nparz(1,ist)
        DO jst=istmin,istmax
        IF(mpsi(jst)==mist .AND. nparz(1,jst)==npzist)
     &    hfbpsi(0:nz,0:nr,1:2,ist) = hfbpsi(0:nz,0:nr,1:2,ist)
     &                      -psi(0:nz,0:nr,1:2,jst)*hlagr(jst,ist)
!     &      /MAX(vsqf(jst),xsmall)
!        IF(mpsi(jst)==mist .AND. nparz(1,jst)==npzist)
!     &  WRITE(933,'(4i6,1pg13.5)') ist,mist,jst,mpsi(jst),hlagr(jst,ist)
        END DO
        IF(.NOT.tbcs) devmf2  = devmf2 + wate(ist)*
     &      wfovlp(hfbpsi(:,:,:,ist),hfbpsi(:,:,:,ist))
      END DO
!      IF(istmax==80) CLOSE(933)
      DEALLOCATE(hlagr,symatrix)

!     maximum gap on the grid
      pptmax = MAXVAL(ABS(pairpt(0:nz,0:nr,iq)))
!      pptmax = 0D0
!      DO ir=0,nr
!      DO iz=0,nz
!        pptmax = MAX(pptmax,ABS(pairpt(iz,ir,iq)))
!      END DO
!      END DO
      IF(ttest) WRITE(7,'(a,i3,1pg13.5)') 'iq,pptmax=',iq,pptmax

!     kinetic energy damping 
      DO ist=istmin,istmax
        IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
          WRITE(7,'(a,i4,3(1pg13.5))') 'HFB overlaps before:',
     &     ist,wfovlp(psi(0,0,1,ist),hfbpsi(0,0,1,ist)),
     &     wfovlp(hfbpsi(0,0,1,ist),hfbpsi(0,0,1,ist)),
     &     wfovlp(psi(0,0,1,ist),psi(0,0,1,ist))
        END IF
        CALL damp_HFB(hfbpsi(:,:,:,ist),vsqf(ist),pptmax,ist)
        IF(tbcs .AND. iq==1 .AND. mpsi(ist)==0) THEN
          WRITE(7,'(a,i4,3(1pg13.5))') 'HFB overlaps after:',
     &     ist,wfovlp(psi(0,0,1,ist),hfbpsi(0,0,1,ist)),
     &     wfovlp(hfbpsi(0,0,1,ist),hfbpsi(0,0,1,ist)),
     &     wfovlp(psi(0,0,1,ist),psi(0,0,1,ist))
        END IF
        psi(:,:,:,ist) = psi(:,:,:,ist)-hfbpsi(:,:,:,ist)
        IF(tbcs) THEN
          WRITE(7,'(a,i4)') 'ortho-normality after: ist=',ist
          WRITE(7,'(100(1pg13.5))') 
     &      (wfovlp(psi(0,0,1,jst),psi(0,0,1,ist)),jst=istmin,istmax)
        END IF
      END DO

!      END DO   ! end of outermost loop over isospin

      RETURN
      END
!-----damp_HFB
!
      SUBROUTINE damp_HFB(psiact,v2act,pptmax,istact)
!
!     Nucleon damping for HFB iteration.
!     The damping matrices are computed for each case anew.
!     Presently, the standard damping strategy is recycled.
!     This should be replaced later by direct inversion.
!
!     Input parameters are:
!      psiact    = s.p. wavefunction to be damped
!      v2act     = occupation of the actual state
!      pptmax    = maximum gap potential
!      istact    = number of actual state
!
      include 'paramr.inc'
!      include 'constr.inc'
!      include 'observ.inc'
!      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'dynpar.inc'
!      include 'resden.inc'
!      include 'crankin.inc'

      REAL(8),INTENT(IN OUT) :: psiact(0:nzd,0:nrd,nsd)
      REAL(8),INTENT(IN) :: v2act,pptmax
      INTEGER,INTENT(IN) :: istact

!      REAL(8),PARAMETER :: v2limh=1D-6,delmix=0.8D0  ! from 'hfbparams' via 'dynpar.inc'
!      REAL(8),PARAMETER :: v2limh=1D0,delmix=0D0   ! setting to switch off
      REAL(8),PARAMETER :: xsmall=1D-20

      LOGICAL,PARAMETER :: ttest=.FALSE.,tdamphfb=.TRUE.

      e0sav = e0inv
      x0sav = x0inv
      iq = icharg(istact)
      IF(istact==istblo(iq) .OR. istact==istcblock(iq)) THEN
        e0inv = e0sav
        x0inv = x0sav        
      ELSE
        uvact = SQRT(MAX(v2act-v2act*v2act,xsmall))
        e0inv = e0inv+pptmax*uvact/MAX(v2act,v2limh)*delmix
        x0inv = x0inv/MAX(v2act,v2limh)
      END IF

      IF(ttest) WRITE(7,'(a,3(1pg13.5))') 
     &    'IDANUC with v2act,e0inv,x0inv=',v2act,e0inv,x0inv
      IF(tdamphfb) THEN
        CALL damnuc_HFB(psiact,psip,istact,e0inv,x0inv)
!     &                  e0inv+pptmax*uvact/MAX(v2act,v2limh)*delmix,
!     &                  x0inv/MAX(v2act,v2limh))
      ELSE
        CALL idanuc
!        step   = x0inv/e0inv
        CALL damnuc(psiact,psip,istact)
      END IF
!      IF(ttest .AND. istact==1) WRITE(7,'(a,2(1pg13.5))')
      IF(ttest) THEN
         WRITE(7,'(a,i5)')  'after DAMNUC:',istact
         WRITE(7,'(a,2(1pg13.5))')
     &   'after DAMNUC:',wfovlp(psiact,psip)
     &    ,wfovlp(psi(0,0,1,istact),psip)
      END IF
      psiact = psip
!      call damp_nucleon(rdnuc,zdnuc,
!     &                  psiact,psiact,nparz(1,ist),mpsi(ist),step)

!
!       set boundary zeroes if parity requires it.                              
!
      do is=1,2
        mm = mpsi(istact)+is-1
        if (ABS(mm).gt.0) psiact(0:nz,0,is)  =  0D0
        if (nparz(is,istact).eq.2) psiact(0,0:nr,is)  =  zero
      enddo


      e0inv = e0sav 
      x0inv = x0sav 
      IF(.NOT.tdamphfb) CALL idanuc


      RETURN
      END


      SUBROUTINE damnuc_HFB(psiin,psiout,istact,e0act,x0act)

!      IMPLICIT NONE
      include "paramr.inc"
      include "derivs.inc"
      include "levels.inc"
      include "dynpar.inc"

      REAL(8), INTENT(IN)  :: e0act,x0act
      REAL(8), INTENT(IN)  :: psiin(0:nzd,0:nrd,nsd)
      REAL(8), INTENT(OUT) :: psiout(0:nzd,0:nrd,nsd)
      INTEGER, INTENT(IN)  :: istact

      REAL(8),ALLOCATABLE :: zdnuc_dyn(:,:,:)        !<Damping operator for the nucleons in z-direction.
      REAL(8),ALLOCATABLE :: rdnuc_dyn(:,:,:)        !<Damping operator for the nucleons in r-direction.
      REAL(8),ALLOCATABLE :: rauxz(:,:)               !<workspace>
      REAL(8),ALLOCATABLE :: rauxr(:,:)               !<workspace>

      REAL(8) :: step,befac
      INTEGER :: mm, is, i, impsi, ipi

      ALLOCATE(zdnuc_dyn(0:nzd,0:nzd,2), rdnuc_dyn(0:nrd,0:nrd,0:kmpsi),
     &         rauxz(0:nzd,0:nzd),rauxr(0:nrd,0:nrd))

! first, the construction of the damping matrices

        befac  = h2mav/(2D0*e0act)
!        WRITE(6,'(a,i5,f8.4)') 'DAMNUC_HFB: 1, istact=',istact,befac
!     nucleon damping in z-direction
        zdnuc_dyn(0:nz,0:nz,1:2) = -befac * derz2f(0:nz,0:nz,1:2)
        DO i=0,nz
          zdnuc_dyn(i,i,:) = one + zdnuc_dyn(i,i,:)
        END DO
!     inversion 
        DO ipi=1,2   
!          CALL matr_inv(zdnuc_dyn(0:nz,0:nz,ipi),nz+1,nz+1)
          CALL minv(zdnuc_dyn(0,0,ipi),rauxz,nzd+1,nz+1)
          IF (ipi == 2.OR.ipi == 4) THEN
            zdnuc_dyn(0,0:nz,ipi) = 0.0D0
            zdnuc_dyn(0:nz,0,ipi) = 0.0D0  ! addition PK
          END IF
        END DO
!        WRITE(6,'(a,2(1pg13.5))') 'DAMNUC_HFB: 2',
!     &    SUM(zdnuc_dyn(0:nz,0:nz,:)),
!     &    SUM(derz2f(0:nz,0:nz,:))

!     the nucleon damping in r-direction
        rdnuc_dyn(0:nr,0:nr,0:kmpsi) = 
     &     -befac * derr2af(0:nr,0:nr,0:kmpsi,1)
        DO i=0,nr
          rdnuc_dyn(i,i,0:kmpsi) =  one + rdnuc_dyn(i,i,0:kmpsi)
        END DO

        mm = mpsi(istact)
        IF(mm.GE.0) THEN
          mmlow = mm
          mmup = mm+1
        ELSE
          mmlow = mm !-1
          mmup = mm+1
        END IF
!        WRITE(7,'(a,4i5)') 'DAMNUC_HFB: mm etc:',mm,istact
!        IF(mm<0) STOP 'DAMNUC_HFB not yet prepared for negative MPSI'
        DO impsi=mmlow,mmup
!          CALL matr_inv(rdnuc_dyn(0:nr,0:nr,impsi),nr+1,nr+1)
          CALL minv(rdnuc_dyn(0,0,ABS(impsi)),rauxr,nrd+1,nr+1)
          IF (ABS(impsi) > 0) THEN
            rdnuc_dyn(0,0:nr,ABS(impsi)) = 0D0
            rdnuc_dyn(0:nz,0,ABS(impsi)) = 0D0 ! addition PK
          END IF
!          WRITE(6,'(a,i4,2(1pg13.5))') 'DAMNUC_HFB: 3',impsi,
!     &      SUM(rdnuc_dyn(0:nr,0:nr,impsi)),
!     &      SUM(derr2af(0:nr,0:nr,impsi,1))
        END DO
!        WRITE(6,'(a,3i5)') 'DAMNUC_HFB: 3, mpsi=',
!     &       istact,mpsi(istact),nparz(1,istact)


!     the straightforward imaginary-time step in case of "e0act=0"

      IF(e0act == 0D0) THEN
        psiout(0:nz,0:nr,1:2) = x0act*psiin(0:nz,0:nr,1:2)
      ELSE
  
!     now the full kinetic energy damping
        istsav=istact
        CALL damp_nucleon(rdnuc_dyn(0,0,0),zdnuc_dyn(0,0,1),
     &                    psiin(0,0,1),psiout(0,0,1),
     &                    nparz(1,istact),mpsi(istact),x0act/e0act)
!     &                    nr,nz,nsd,kmpsi,isym)

!       set boundary zeroes if parity requires it.
!        WRITE(6,'(a,3i5)') 'DAMNUC_HFB: 4, istact=',istact,istsav
!        istact=istsav
        DO is=1,2
          mm = mpsi(istsav)+is-1
!          WRITE(6,'(a,3i5)') 'DAMNUC_HFB: 5, istact..=',istact,is,mm
          IF (mm > 0) THEN
            psiout(0:nz,0,is)  =  0D0
          END IF
          IF (nparz(is,istsav) == 2) THEN
            psiout(0,0:nr,is)  =  0D0
          END IF
        END DO
  
      END IF
!        WRITE(6,'(a)') 'DAMNUC_HFB: 6'

      DEALLOCATE(zdnuc_dyn, rdnuc_dyn, rauxz, rauxr)

      RETURN
      END SUBROUTINE damnuc_HFB
!      SUBROUTINE damp_nucleon(dampr,dampz,fldin,fldout,izpfld,imfld,
!     &                        factor,nr,nz,nsd,kmpsi,isym)
!
!      IMPLICIT NONE
!
!      REAL(8), INTENT(IN)      :: dampr(0:nr,0:nr,0:kmpsi)
!      REAL(8), INTENT(IN)      :: dampz(0:nz,0:nz,4)
!      REAL(8), INTENT(IN)      :: fldin(0:nz,0:nr,nsd)
!      REAL(8), INTENT(OUT)     :: fldout(0:nz,0:nr,nsd)
!      INTEGER, INTENT(IN)      :: izpfld
!      INTEGER, INTENT(IN)      :: imfld
!      INTEGER, INTENT(IN)      :: nr,nz,nsd,kmpsi,isym
!      REAL(8), INTENT(IN)      :: factor
!
!      REAL(8)                     :: fldbuf(0:nz,0:nr)
!
!      INTEGER :: imact
!      INTEGER :: izpact
!      INTEGER :: is
!      INTEGER :: isigma
!
!      DO is=1,nsd
!  
!! lower spin component has m-component imact+1
!        imact  = ABS(imfld+is-1)
!  
!!       lower spin component has opposite z-parity if isym=1
!!       if isym=0 izpact should be always +1
!        isigma = 3-2*is
!        izpact = (3*(1-isigma)/2+izpfld*isigma)*isym + (1-isym)
!  
!        CALL mulr(dampr(:,:,imact),fldin(:,:,is),fldbuf      ,1D0   ,0)
!        CALL mulz(dampz(:,:,izpact),fldbuf,    fldout(:,:,is),factor,0)
!      END DO
!
!      RETURN
!      END SUBROUTINE damp_nucleon
