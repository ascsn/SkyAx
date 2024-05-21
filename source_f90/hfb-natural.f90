!-----hfb_natural---------------------------------------------------------

! Code converted using TO_F90 by Alan Miller
! Date: 2024-05-21  Time: 13:57:40

SUBROUTINE hfb_natural(iso,efermiso,partnumin)

USE paramr

!     Solves the HFB equations for natural orbitals as matrix equations
!     in the space of given s.p. states for nucleon type 'iso'.
!     Version for the 2D code.

USE wavefr
USE grid
USE dynpar
USE levels
USE fields
USE observ

INTEGER, INTENT(IN)      :: iso
REAL(8), INTENT(IN OUT)  :: efermiso
REAL(8), INTENT(IN)      :: partnumin


INTEGER :: istbloact
REAL(8) :: partnum
REAL(8), ALLOCATABLE :: cm(:,:),cmnew(:,:),symcond(:,:)
REAL(8), ALLOCATABLE :: hcm(:,:)
REAL(8), ALLOCATABLE :: hmatrix(:,:),gapmatrix(:,:)
REAL(8), ALLOCATABLE :: phasfac(:),paircut_uv(:),paircut_v2(:),spenat(:),gapnat(:)
REAL(8), ALLOCATABLE :: gwa(:),v2a(:),uva(:),lambda(:,:),dummy(:)
INTEGER, ALLOCATABLE :: jssav(:)
REAL(8), ALLOCATABLE :: psip(:,:,:),psin(:,:,:)

LOGICAL, PARAMETER :: tprint=.false.,tprint0=.false.
LOGICAL, PARAMETER :: tprint1=.false.
LOGICAL, PARAMETER :: toffdiag=.true.,tnewv2=.false.
LOGICAL, PARAMETER :: tbuffer=.true.    ! use buffering for matrix products
INTEGER, PARAMETER :: ifcut=1
!      INTEGER,PARAMETER :: iternat=0!100
!      REAL(8),PARAMETER :: delstepbas=1.0D0,e0bas=0D0
REAL(8),PARAMETER :: delstepbas=1.0D0,e0bas=20D0

!--------------------------------------------------------------------

IF(iternat<0) RETURN
partnum=partnumin

!?      isoadd = (iso-1)*ngrid
!?      isnpad = (iso-1)*kjmx
IF(iso==1) THEN
  nmaxa=nst(1)
  isonpl=0
  efadd=feradd(1)
  efermiso=fer(1)
  istbloact=istblo
ELSE IF(iso==2) THEN
  nmaxa=nst(2)
  isonpl=nst(1)
  efadd=feradd(2)
  efermiso=fer(2)
  istbloact=istblo-nst(1)
  WRITE(6,'(i4)') npart(2)
ELSE
  STOP 'HFB: ISO must be 1 or 2'
END IF
IF(tprint0) WRITE(6,'(a,2i4,4(1pg13.5))')  &
    'HFB natural for iso,nmaxa,eferm,partnum=',iso,nmaxa,efermiso, partnum
ALLOCATE(lambda(nmaxa,nmaxa),cmnew(nmaxa,nmaxa),gapnat(nmaxa))
ALLOCATE(hmatrix(nmaxa,nmaxa),gapmatrix(nmaxa,nmaxa))
ALLOCATE(gwa(nmaxa),hcm(nmaxa,nmaxa))
ALLOCATE(paircut_uv(nmaxa),paircut_v2(nmaxa),phasfac(nmaxa),spenat(nmaxa))
ALLOCATE(jssav(nmaxa),v2a(nmaxa),uva(nmaxa),cm(nmaxa,nmaxa))
ALLOCATE(dummy(nmaxa),symcond(nmaxa,nmaxa))

!     copy old weights to new workspace

DO nq=1,nmaxa
  nqi=nq+isonpl
  v2a(nq)=vsqf(nqi)
  paircut_uv(nq)=wtpair_uv(nqi)
  paircut_v2(nq)=wtpair_v2(nqi)
END DO
gwa=v2a

!     set up Hamiltonian matrix, gap matrix, and initial value for
!     unitary transformation to natural orbits

!      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)

IF(tprint)  WRITE(6,'(a,i4,1pg12.5)')  &
    'HFB natural at entry: iso=',iso,sum(pairpt(0:nz,0:nr,iso)),  &
    ' n1  ist m par  spe   delta  spenew   gap  variance'
ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))
hmatrix=0D0
gapmatrix=0D0
cm=0D0
DO n1=1,nmaxa
  ist=n1+isonpl
  m1=mpsi(ist)
  np1=nparz(1,ist)
  CALL hpsi(psi(0,0,1,ist),psin,ist)
  DO is=1,nsd
    psip(0:nz,0:nr,is)=psi(0:nz,0:nr,is,ist)*pairpt(0:nz,0:nr,iso)
    psip(0:nz,0:nr,is)=psi(0:nz,0:nr,is,ist)*pairpt(0:nz,0:nr,iso)
  END DO
  DO n2=1,nmaxa
    jst=n2+isonpl
    IF((toffdiag .OR. n2==n1) .AND.  &
          (mpsi(jst)==m1 .AND. nparz(1,jst)==np1) ) THEN
      hmatrix(n2,n1) = wfovlp(psi(0:nz,0:nr,1:2,jst),psin)
      gapmatrix(n2,n1) = wfovlp(psi(0:nz,0:nr,1:2,jst),psip)
    END IF
  END DO
  IF(tprint) THEN
    WRITE(6,'(4i4,6(1pg13.5))') n1,ist,mpsi(ist),nbuf(ist),  &
        spenrg(ist),deltaf(ist),hmatrix(n1,n1)+spenrg(ist),  &
        gapmatrix(n1,n1),SQRT(wfovlp(psin,psin)-hmatrix(n1,n1)**2),  &
        SQRT(sum(hmatrix(1:nmaxa,n1)**2)),hmatrix(n1,n1)
  END IF
  hmatrix(n1,n1)= hmatrix(n1,n1)+spenrg(ist)
  cm(n1,n1) = 1D0
  spenat(n1) = hmatrix(n1,n1)  !+efermiso
  gapnat(n1) = gapmatrix(n1,n1)
END DO
DEALLOCATE(psip,psin)

CALL sortenerg(nmaxa,spenat,iorder(1+isonpl))     ! sorting to energy
IF(tprint1) WRITE(6,'(a,5000(2i4,2x))') 'IORDER:',  &
    (na+isonpl,iorder(na+isonpl)+isonpl,na=1,nmaxa)
IF(tnewv2) THEN
  IF(tprint0 .AND. iso==2) THEN
    WRITE(6,'(a,2(1pg13.5))') 'before PAIRDN:'
    WRITE(6,'(a,200(1pg13.5))') 'v2a:   ',v2a(1:nmaxa)
    WRITE(6,'(a,200(1pg13.5))') 'gapnat:',gapnat(1:nmaxa)
  END IF
  IF(pairstab /= 0)  &
      CALL pair_stab(spenat,gapnat,v2a,wate(isonpl+1),paircut_uv,epair,  &
      pairstab,nmaxa,nmaxa,iab,tcutoff_uvspace)
  CALL pairdn(spenat,gapnat,v2a,wate(isonpl+1),paircut_uv,paircut_v2,partnum,  &
      efermiso,epair(iso),gap(iso),efadd,  &
      nst(iso),nst(iso),npiter,ipair,istbloact,erra,itest-1)  
  IF(tprint0) WRITE(6,'(a,6(1pg13.5))') 'after PAIRDN:',  &
      efermiso,epair(iso),gap(iso),efadd
  IF(tprint0 .AND. iso==2) THEN
    WRITE(6,'(a,200(1pg13.5))') 'v2a:   ',v2a(1:nmaxa)
    WRITE(6,'(a,200(1pg13.5))') 'gapnat:',gapnat(1:nmaxa)
  END IF
END IF
IF(tprint) THEN
  hermh=0D0
  hermg=0D0
  hermhm=0D0
  DO na=1,nmaxa; DO nb=1,nmaxa
    hermh = hermh + (hmatrix(na,nb)-hmatrix(nb,na))**2
    hermg = hermg + (gapmatrix(na,nb)-gapmatrix(nb,na))**2
    IF((hmatrix(na,nb)-hmatrix(nb,na))**2> hermhm) THEN
      hermhm=(hmatrix(na,nb)-hmatrix(nb,na))**2
      namax=na
      nbmax=nb
    END IF
  END DO; END DO
  hermhm= SQRT(hermhm)
  WRITE(6,'(a,3(1pg13.5),2i4)') 'test hermiticity:',  &
      SQRT(hermh/nmaxa**2),SQRT(hermg/nmaxa**2),hermhm,namax+isonpl, nbmax+isonpl
END IF
IF(tcutoff_uvspace) THEN                             ! check that!
  uva(:)=SQRT((1D0-v2a(:))*v2a(:))*paircut_uv(:)      
ELSE
  uva(:)=SQRT((1D0-v2a(:))*v2a(:))
END IF

avsymcond=0D0
symcondmax=0D0
IF(tbuffer) THEN
  DO na=1,nmaxa
    ist=na+isonpl
    istbuf=nbuf(ist)
    np1=nparz(1,ist)
    nblo=mlo(istbuf)-isonpl
    nbhi=mhi(istbuf)-isonpl
    IF(iso==2) WRITE(6,'(a,3i4)') 'istbuf,nblo,nbhi=',nblo,nbhi,istbuf
    DO nb=nblo,nbhi
      IF(iso==2) WRITE(6,'(a,2i4)') 'NA,NB=',na,nb
      IF(nparz(1,nb+isonpl) /= np1) cycle
      symcond(na,nb) = +(v2a(nb)*hmatrix(na,nb)-v2a(na)*hmatrix(nb,na))  &
          +(uva(nb)*gapmatrix(na,nb)-uva(na)*gapmatrix(nb,na))
      symcondmax=MAX(symcondmax,ABS(symcond(na,nb)))
      avsymcond=avsymcond+symcond(na,nb)**2
    END DO
  END DO
ELSE
  DO na=1,nmaxa; DO nb=1,nmaxa
    symcond(na,nb) = +(v2a(nb)*hmatrix(na,nb)-v2a(na)*hmatrix(nb,na))  &
        +(uva(nb)*gapmatrix(na,nb)-uva(na)*gapmatrix(nb,na))
    symcondmax=MAX(symcondmax,ABS(symcond(na,nb)))
    avsymcond=avsymcond+symcond(na,nb)**2
  END DO; END DO
END IF
symcondi(iso)=SQRT(avsymcond/partnum)
WRITE(6,'(a,2i5,1pg13.5)') 'HFB natural: iso,iter,symm=',  &
    iso,0,SQRT(avsymcond/partnum)
!      WRITE(7,'(a,2i5,1pg13.5)') 'HFB natural: iso,iter,symm=',
!     &         iso,0,SQRT(avsymcond/partnum)
IF(tprint0) THEN
  WRITE(6,'(a,i5)') 'HFB natural: iter=',0
  WRITE(6,'(a,7(1pg13.5))') ' partnum,efermi,epair,cutfac=',  &
      partnum,sum(v2a*phasfac)*2D0,efermiso,epair,cutfac
  IF(tprint1) THEN
    WRITE(6,'(a)')'# n  js  spe  v2  uv  cut   phasfac'
    DO na=1,nmaxa
      WRITE(6,'(2i3,5(1pg13.5))') na,jssav(na),spenat(na),v2a(na),  &
          gapnat(na),paircut_uv(na),phasfac(na)
      nb=na+isonpl
      WRITE(6,'(6x,5(1pg13.5))') spenrg(nb),vsqf(nb),  &
          deltaf(nb),paircut_uv(nb),wate(nb)
    END DO
  END IF
!          WRITE(6,'(a)') 'asymmetry matrix:'
!          DO na=1,nmaxa
!            WRITE(6,'(i3,50(1pg13.5))') na,symcond(:,na)
!          END DO
  WRITE(6,'(a,i2,4(1pg13.5))')  &
      'HFB-natural init.symm. condition: iso,max,average=',  &
      iso,symcondmax,symcondi(iso)
  IF(tprint1) THEN
    WRITE(6,'(a)') 'h matrix:'
    DO na=1,nmaxa
      WRITE(6,'(i3,a,500(1pg13.5))') na,': ',hmatrix(:,na)
    END DO
    WRITE(6,'(a)') 'gap matrix:'
    DO na=1,nmaxa
      WRITE(6,'(i3,a,500(1pg13.5))') na,': ',gapmatrix(:,na)
    END DO
  END IF
END IF

IF(iternat <= 0) RETURN

!     iterate

iab=-2
IF(tprint1) iab=-1
DO iter=1,iternat
  
  
!     h->psi in configuration space
  IF(tbuffer) THEN
    DO na=1,nmaxa
      IF(thfbvariat) THEN       ! check for paircut_v2!
        fackincut = -(1D0-paircut_uv(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
        v2avar = v2a(na) + uva(na)*fackincut
      ELSE
        v2avar = v2a(na)
      END IF
      ist=na+isonpl
      istbuf=nbuf(ist)
      np1=nparz(1,ist)
      nblo=mlo(istbuf)-isonpl
      nbhi=mhi(istbuf)-isonpl
      DO n1=nblo,nbhi
        IF(nparz(1,n1+isonpl) /= np1) cycle
        hcm(n1,na) = v2avar*sum(hmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))+  &
            uva(na)*sum(gapmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))
      END DO
    END DO
  ELSE
    DO na=1,nmaxa
      IF(thfbvariat) THEN          ! check for paircut_v2!
        fackincut = -(1D0-paircut_uv(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
        v2avar = v2a(na) + uva(na)*fackincut
      ELSE
        v2avar = v2a(na)
      END IF
      DO n1=1,nmaxa
        hcm(n1,na) = v2avar*sum(hmatrix(n1,:)*cm(:,na))+  &
            uva(na)*sum(gapmatrix(n1,:)*cm(:,na))
      END DO
    END DO
  END IF
!     matrix of Lagrangian parameters and symmetry condition
!        lambda=0D0
!        uva(:)=uva(:)*paircut(:)
  symcond=0D0
  lambda=0D0
  IF(tbuffer) THEN
    DO na=1,nmaxa
      ist=na+isonpl
      istbuf=nbuf(ist)
      np1=nparz(1,ist)
      nblo=mlo(istbuf)-isonpl
      nbhi=mhi(istbuf)-isonpl
      DO nb=nblo,nbhi
        IF(nparz(1,nb+isonpl) /= np1) cycle
        cbhca = sum(cm(nblo:nbhi,nb)*hcm(nblo:nbhi,na))
        cahcb = sum(cm(nblo:nbhi,na)*hcm(nblo:nbhi,nb))
        lambda(na,nb)= 0.5D0*(cbhca+cahcb)
        lambda(nb,na)=lambda(na,nb)
        symcond(na,nb)= (cbhca-cahcb)
        symcond(nb,na)= -symcond(na,nb)
      END DO
    END DO
  ELSE
    DO na=1,nmaxa; DO nb=1,na
      cbhca = sum(cm(:,nb)*hcm(:,na))
      cahcb = sum(cm(:,na)*hcm(:,nb))
      lambda(na,nb)= 0.5D0*(cbhca+cahcb)
      lambda(nb,na)=lambda(na,nb)
      symcond(na,nb)= (cbhca-cahcb)
      symcond(nb,na)= -symcond(na,nb)
    END DO; END DO
  END IF
  IF(tprint0) WRITE(6,'(a,1pg13.5)') 'before step:',  &
      SQRT(sum(symcond(:,:)**2)/partnum)
!        nitsub=1
!        DO itsub=1,nitsub
!     gradient step
  spemin=minval(spenat,nmaxa)
  spemax=maxval(spenat,nmaxa)
  delstep=delstepbas/(spemax-spemin)
  IF(tbuffer) THEN
    DO na=1,nmaxa
      IF(thfbvariat) THEN     !check for paircut_v2 !
        fackincut = -(1D0-paircut_uv(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
        v2avar = v2a(na) + uva(na)*fackincut
      ELSE
        v2avar = v2a(na)
      END IF
      ist=na+isonpl
      istbuf=nbuf(ist)
      np1=nparz(1,ist)
      nblo=mlo(istbuf)-isonpl
      nbhi=mhi(istbuf)-isonpl
      DO n1=nblo,nbhi
        IF(nparz(1,n1+isonpl) /= np1) cycle
        IF(e0bas==0D0) THEN
          dampstep = delstep
        ELSE
          dampstep = delstepbas/(hmatrix(n1,n1)-spemin+e0bas)
!    &                  /nitsub
        END IF
        cmnew(n1,na) = cm(n1,na) - dampstep*(hcm(n1,na)  &
            -  sum(cm(n1,nblo:nbhi)*lambda(nblo:nbhi,na))     )
!     &          v2avar*SUM(hmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))
!     &         +uva(na)*SUM(gapmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))
!     &       -  SUM(cm(n1,nblo:nbhi)*lambda(nblo:nbhi,na))     )
      END DO
    END DO
  ELSE
    DO na=1,nmaxa;
      IF(thfbvariat) THEN    ! check for paircut_v2
        fackincut = -(1D0-paircut_uv(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
        v2avar = v2a(na) + uva(na)*fackincut
      ELSE
        v2avar = v2a(na)
      END IF
      DO n1=1,nmaxa
        IF(e0bas==0D0) THEN
          dampstep = delstep
        ELSE
          dampstep = delstepbas/(hmatrix(n1,n1)-spemin+e0bas)
!    &                  /nitsub
        END IF
        cmnew(n1,na) = cm(n1,na) - dampstep*(  &
            v2avar*sum(hmatrix(n1,:)*cm(:,na))  &
            +uva(na)*sum(gapmatrix(n1,:)*cm(:,na))  &
            -sum(cm(n1,:)*lambda(:,na))     )
      END DO
    END DO
  END IF
  
!     ortho-normalize
  cm=cmnew
  IF(tbuffer) THEN
    DO nae=1,nmaxa               ! loop runs energy sorted
      ist=iorder(nae+isonpl)+isonpl
      na=ist-isonpl
      istbuf=nbuf(ist)
      np1=nparz(1,ist)
      nblo=mlo(istbuf)-isonpl
      nbhi=mhi(istbuf)-isonpl
      DO nbe=1,nae-1
        jst=iorder(nbe+isonpl)+isonpl
        nb=jst-isonpl
        IF(nbuf(jst) /= istbuf .OR. nparz(1,jst) /= np1) cycle
        ovlab=sum(cm(nblo:nbhi,nb)*cm(nblo:nbhi,na))
        cm(nblo:nbhi,na)=cm(nblo:nbhi,na)-cm(nblo:nbhi,nb)*ovlab
      END DO
      ovlaa=1D0/SQRT(sum(cm(nblo:nbhi,na)*cm(nblo:nbhi,na)))
      cm(nblo:nbhi,na)=cm(nblo:nbhi,na)*ovlaa
    END DO
  ELSE
    DO nae=1,nmaxa               ! loop runs energy sorted
      ist=iorder(nae+isonpl)+isonpl
      na=ist-isonpl
      m1=mpsi(ist)
      np1=nparz(1,ist)
      DO nbe=1,nae-1
        jst=iorder(nbe+isonpl)+isonpl
        nb=jst-isonpl
        IF(mpsi(jst) /= m1 .OR. nparz(1,jst) /= np1)  cycle
        ovlab=sum(cm(:,nb)*cm(:,na))
        cm(:,na)=cm(:,na)-cm(:,nb)*ovlab
      END DO
      ovlaa=1D0/SQRT(sum(cm(:,na)*cm(:,na)))
      cm(:,na)=cm(:,na)*ovlaa
    END DO
  END IF
  
!        END DO          ! end of sub-iterations
!     new energies and occupations
  DO na=1,nmaxa
    spenat(na) = 0D0
    gapnat(na) = 0D0
    DO n1=1,nmaxa
      spenat(na) =spenat(na)+cm(n1,na)*sum(hmatrix(n1,:)*cm(:,na))
      gapnat(na) = gapnat(na)+ cm(n1,na)*sum(gapmatrix(n1,:)*cm(:,na))
    END DO
  END DO
  IF(pairstab /= 0)  &
      CALL pair_stab(spenat,gapnat,v2a,wate(isonpl+1),paircut_uv,epair,  &
      pairstab,nmaxa,nmaxa,iab,tcutoff_uvspace)
  CALL pairdn(spenat,gapnat,v2a,wate(isonpl+1),paircut_uv(isonpl+1),paircut_v2(isonpl+1),partnum,  &
      efermiso,epair(iso),gap(iso),efadd,  &
      nst(iso),nst(iso),npiter,ipair,istbloact,erra,itest-1)
  IF(tcutoff_uvspace) THEN                               ! check
    uva(:)=SQRT((1D0-v2a(:))*v2a(:))*paircut_uv(:)
  ELSE
    uva(:)=SQRT((1D0-v2a(:))*v2a(:))
  END IF
  avsymcond=0D0
  symcondmax=0D0
  DO na=1,nmaxa
    avsymcond=avsymcond+sum(symcond(:,na)**2)
    symcondmax=MAX(symcondmax,maxval(ABS(symcond(:,na))))
  END DO
  IF(tprint0) THEN
    WRITE(6,'(a,i5)') 'HFB natural: iter=',iter
    WRITE(6,'(a,6(1pg13.5))') 'partnum,eferm,epair,delst,symm=',  &
        partnum,efermi,epairact,delstep, symcondmax,SQRT(avsymcond/partnum)
  END IF
  IF(tprint1) THEN
    WRITE(6,'(a)')'# n  js  spe  v2  uv  cut   phasfac'
    DO na=1,nmaxa
      WRITE(6,'(2i3,5(1pg13.5))') na,jssav(na),spenat(na),v2a(na),  &
          uva(na),paircut_uv(na),phasfac(na)
    END DO
    IF(tprint) THEN
      WRITE(6,'(a)') 'symmetry condition in detail:'
      DO na=1,nmaxa
        WRITE(6,'(i3,a,500(1pg13.5))') na,': ',symcond(:,na)
      END DO
    END IF
  ELSE
!          WRITE(7,'(a,2i5,1pg13.5)') 'HFB natural: iso,iter,symm=',
!     &         iso,iter,SQRT(avsymcond/partnum)
  END IF
  
END DO


IF(tprint) THEN
  WRITE(6,'(a,3(1pg13.5))') 'HFB natural exit: n1  spenew  gap'
  DO n1=1,nmaxa
    WRITE(6,'(i4,3(1pg13.5))') n1,hmatrix(n1,n1),gapmatrix(n1,n1)
  END DO
END IF
symcondf(iso)=SQRT(avsymcond/partnum)
WRITE(6,'(a,i2,4(1pg13.5))')  &
    'HFB final symmetry condition: iso,max,average=', iso,symcondmax,symcondf(iso)
!      CALL  hfb_nat_densit(iso,nmaxa,cm,v2a,uva,jssav,phasfac)
CALL  hfb_nat_mapwf(iso,nmaxa,cm,v2a,paircut_uv(isonpl+1),paircut_v2(isonpl+1),jssav,phasfac)

DEALLOCATE(lambda,cmnew,gapnat,hmatrix,gapmatrix)
DEALLOCATE(paircut_uv,paircut_v2,phasfac,spenat,jssav,v2a,uva,cm)
DEALLOCATE(dummy,symcond)

IF(iso==1) THEN
  feradd(1)=efadd
  fer(1)=efermiso
ELSE
  feradd(2)=efadd
  fer(2)=efermiso
END IF

!      IF(iso==2) STOP

RETURN
END SUBROUTINE hfb_natural
!-----hfb_nat_mapwf---------------------------------------------------------

SUBROUTINE hfb_nat_mapwf(iso,nmaxa,cm,v2a,paircut_uv,paircut_v2,jssav,phasfac)

USE paramr

!     Map new wavefunctions from configuration space back to coordinate space.
!     Version for 2D. Misuse array 'dpsiz' temporarily as workspace.

USE wavefr
USE dynpar
USE levels


INTEGER, INTENT(IN)    :: iso
INTEGER, INTENT(IN)    :: nmaxa
REAL(8), INTENT(IN)    :: cm(nmaxa,nmaxa)
REAL(8), INTENT(IN)    :: v2a(nmaxa)
REAL(8), INTENT(IN)    :: paircut_uv(nmaxa)
REAL(8), INTENT(IN)    :: paircut_v2(nmaxa)
INTEGER, INTENT(IN)    :: jssav(nmaxa)          ! jssav
REAL(8), INTENT(IN)    :: phasfac(nmaxa)

LOGICAL, PARAMETER :: tprint0=.false.,tprint=.false.
REAL(8), ALLOCATABLE :: psisav(:,:,:,:)

!-------------------------------------------------------------------------------


IF(iso==1) THEN
  IF(nmaxa /= nst(1)) STOP 'HFB_nat_densit: wrong # of states'
  isonpl=0
ELSE IF(iso==2) THEN
  IF(nmaxa /= nst(2)) STOP 'HFB_nat_densit: wrong # of states'
  isonpl=nst(1)
ELSE
  STOP 'HFB: ISO must be 1 or 2'
END IF

!     new wave functions from HFB coefficients

ALLOCATE(psisav(0:nz,0:nr,nsd,npsid))
acctot=0D0
psisav=0D0
DO n1=1,nmaxa
  ist=n1+isonpl
  DO n2=1,nmaxa
    jst=n2+isonpl
    psisav(0:nz,0:nr,1:2,ist) = psisav(0:nz,0:nr,1:2,ist)  &
        + cm(n2,n1)*psi(0:nz,0:nr,1:2,jst)
  END DO
  vsqf(ist)=v2a(n1)
  wtpair_uv(ist)=paircut_uv(n1)
  wtpair_v2(ist)=paircut_v2(n1)
END DO

!     move to correct place

DO n1=1,nmaxa
  ist=n1+isonpl
  psi(0:nz,0:nr,1:2,ist) = psisav(0:nz,0:nr,1:2,ist)
END DO
DEALLOCATE(psisav)

RETURN
END SUBROUTINE hfb_nat_mapwf
!-----pair_stab------------ part of the hartree-fock package -----------------

!     modify gap according to stabilized pairing

SUBROUTINE pair_stab(e,delta,v2occ,ph,pairgw,epair,  &
    pairstab,nmax,ndim,iab,tphpair)

!     input via list:
!      e       = array(1..nmax) of single particle energies.
!      delta   = array(1..nmax) of gap parameters.
!      v2occ   = array(1..nmax) of occupation probabilities,
!      ph      = array(1..nmax) of the multiplicities of the single
!                particle levels.
!                This multiplicities excludes (!) the spin degeneracy,
!                e.g. a 1s1/2 state has ph = 1 .
!      pairgw  = space for dynamical cut-off factor (uv version)  ==> check!
!      pairstab  = stabilizing energy
!      nmax    = maximum number of single particle states
!                (needs to be smaller than the parameter 'ndim')
!      ndim    = dimensioning of the arrays, particularly of 'vpairf',
!                in the calling routine
!      iab     = print level; significant values are -2, -1, 0 and 1;
!                larger value prints more.
!      tphpair = switch to old pairing recipe
!     output via list:
!      delta   = array(1..nmax) of gap parameters.
!      epair   = pairing energy

USE paramr

REAL(8), INTENT(IN OUT)  :: e(ndim)
REAL(8), INTENT(IN OUT)  :: delta(ndim)
REAL(8), INTENT(IN)      :: v2occ(ndim)
REAL(8), INTENT(IN)      :: ph(ndim)
REAL(8), INTENT(IN)      :: pairgw(ndim)
REAL, INTENT(OUT)        :: epair
REAL, INTENT(IN)         :: pairstab
INTEGER, INTENT(IN)      :: nmax
INTEGER, INTENT(IN OUT)  :: ndim
INTEGER, INTENT(IN)      :: iab
LOGICAL, INTENT(IN)      :: tphpair


!LOGICAL :: tphpair

LOGICAL, PARAMETER :: tprint=.false.

!      IF(pairstab==0D0) RETURN

IF(tprint) THEN
  WRITE(6,'(a,i4,200(1pg13.5))') 'PAIR_STAB entry: ',nmax
  WRITE(6,'(a,200(1pg13.5))') ' delta: ',delta(1:nmax)
  WRITE(6,'(a,200(1pg13.5))') ' v2occ: ',v2occ(1:nmax)
  WRITE(6,'(a,200(1pg13.5))') 'pairgw: ',pairgw(1:nmax)
  WRITE(6,'(a,200(1pg13.5))') '    ph: ',ph(1:nmax)
END IF

!     compute pairing energy

sumduv = 0D0
DO  na=1,nmax
  IF(tphpair) THEN
    deltaa = delta(na)*pairgw(na)   !  old cutoff
  ELSE
    deltaa = delta(na)  ! new cutoff
  END IF
!         deltab = delta(na)
!         edif   = e(na)-eferm
!         equasi = sqrt(edif*edif+deltaa*deltaa)
!         v2     = (0.5D0-0.5D0*edif/equasi)
  v2     = v2occ(na)
  uv     = SQRT(MAX(v2-v2*v2,xsmall))
  IF(.NOT.tphpair) THEN
    v2     = v2*pairgw(na)   ! new cutoff
    uv     = uv*pairgw(na)   ! new cutoff
  END IF
  sumduv = ph(na)*uv*deltaa + sumduv
END DO

epair  = -sumduv
IF(iab > 1) WRITE(6,'(a,3g12.4)') 'PAIR_STAB epair:',epair

IF(pairstab==0D0) RETURN

!     modify gaps

DO na=1,nmax
  delta(na) = delta(na)*(1D0+(pairstab/epair)**2)
  delta(na) = MIN(delta(na),5D0)
END DO


RETURN
END SUBROUTINE pair_stab


