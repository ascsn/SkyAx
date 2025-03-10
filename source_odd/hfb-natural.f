c-----hfb_eqs---------------------------------------------------------
c
      subroutine hfb_natural(iso,efermiso,partnumin)
c
      include 'paramr.inc'
c
c     Solves the HFB equations for natural orbitals as matrix equations 
c     in the space of given s.p. states for nucleon type 'iso'.
c     Version for the 2D code.
c
c
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'levels.inc'
      include 'fields.inc'
      include 'observ.inc'
c
      REAL(8), INTENT(IN OUT) :: efermiso
      REAL(8), INTENT(IN) :: partnumin
      INTEGER, INTENT(IN) :: iso  !,iternat
c
      REAL(8) :: partnum,epairstab,conocc
      REAL(8), ALLOCATABLE :: cm(:,:),cmnew(:,:),symcond(:,:)
      REAL(8), ALLOCATABLE :: hcm(:,:)
      REAL(8), ALLOCATABLE :: hmatrix(:,:),gapmatrix(:,:)
      REAL(8), ALLOCATABLE :: phasfac(:),paircut(:),spenat(:),gapnat(:)
      REAL(8), ALLOCATABLE :: gwa(:),v2a(:),uva(:),lambda(:,:),dummy(:)
      INTEGER, ALLOCATABLE :: jssav(:)
      LOGICAL,PARAMETER :: tprint=.FALSE.,tprint0=.TRUE.
      LOGICAL,PARAMETER :: tprint1=.TRUE.
      LOGICAL,PARAMETER :: toffdiag=.TRUE.,tnewv2=.FALSE.
      LOGICAL,PARAMETER :: tbuffer=.FALSE.    ! use buffering for matrix products
      LOGICAL :: toccmix = .FALSE.   ! switch to mixing occupations  old/new
      INTEGER,PARAMETER :: ifcut=1
!      INTEGER,PARAMETER :: iternat=0!100
      REAL(8),PARAMETER :: delstepbas=0.0D0,e0bas=40D0
!      REAL(8),PARAMETER :: delstepbas=1.0D0,e0bas=20D0
c
c--------------------------------------------------------------------
c
      IF(iternat.GE.0) 
     &   WRITE(6,'(a)') 'HFB-NATURAL presently disabled, return early'
!      iternat=-1
!      IF(iternat<0) RETURN
      partnum=partnumin
      IF(addocc<1D0) THEN
        toccmix=.TRUE.
        conocc = 1D0-addocc
      END IF
c
!?      isoadd = (iso-1)*ngrid
!?      isnpad = (iso-1)*kjmx
      IF(iso==1) THEN
        nmaxa=nst(1)
        isonpl=0
        efadd=feradd(1)
        efermiso=fer(1)
!        WRITE(6,'(a,2000(1pg13.5))') 'iso=1 WATE:',wate(1:nst(1))
      ELSE IF(iso==2) THEN
        nmaxa=nst(2)
        isonpl=nst(1)
        efadd=feradd(2)
        efermiso=fer(2)
!        WRITE(6,'(a,3i4)') 'NPART(2):',npart(2),nst(1),nst(1)+nst(2)
!        WRITE(6,'(a,2000(1pg13.5))') 'iso=2 WATE:',
!     &         wate(nst(1)+1:nst(1)+nst(2))
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
      IF(tprint0) WRITE(6,'(a,2i4,4(1pg13.5))') 
     &  'HFB natural for iso,nmaxa,eferm,partnum=',iso,nmaxa,efermiso,
     &   partnum
      ALLOCATE(lambda(nmaxa,nmaxa),cmnew(nmaxa,nmaxa),gapnat(nmaxa))
      ALLOCATE(hmatrix(nmaxa,nmaxa),gapmatrix(nmaxa,nmaxa))
      ALLOCATE(gwa(nmaxa),hcm(nmaxa,nmaxa))
      ALLOCATE(paircut(nmaxa),phasfac(nmaxa),spenat(nmaxa))
      ALLOCATE(jssav(nmaxa),v2a(nmaxa),uva(nmaxa),cm(nmaxa,nmaxa))
      ALLOCATE(dummy(nmaxa),symcond(nmaxa,nmaxa))
c
c     copy old weights to new workspace 
c      
      DO nq=1,nmaxa
        nqi=nq+isonpl
        v2a(nq)=vsqf(nqi)
        paircut(nq)=wtpair(nqi)
        phasfac(nq)=wate(nqi)
      END DO
      gwa=v2a
c
c     set up Hamiltonian matrix, gap matrix, and initial value for
c     unitary transformation to natural orbits
c
!      call derv15(bmass(isoadd+1),d1b(isoadd+1),one,ngrid,rstep)

       IF(tprint)  WRITE(6,'(a,i4,1pg12.5)') 
     &     'HFB natural at entry: iso=',iso,SUM(pairpt(0:nz,0:nr,iso)),
     &     ' n1  ist m par  spe   delta  spenew   gap  variance'
      hmatrix=0D0
      gapmatrix=0D0
      cm=0D0
      DO n1=1,nmaxa
        ist=n1+isonpl
        m1=mpsi(ist)
        np1=nparz(1,ist)
        call hpsi(psi(0,0,1,ist),psin,ist)
        DO is=1,nsd
          psip(0:nz,0:nr,is)=psi(0:nz,0:nr,is,ist)*pairpt(0:nz,0:nr,iso)
          psip(0:nz,0:nr,is)=psi(0:nz,0:nr,is,ist)*pairpt(0:nz,0:nr,iso)
        END DO
        DO n2=1,nmaxa
          jst=n2+isonpl
          IF((toffdiag .OR. n2==n1) .AND. 
     &       (mpsi(jst)==m1 .AND. nparz(1,jst)==np1) ) THEN
            hmatrix(n2,n1) = wfovlp(psi(0,0,1,jst),psin)
            gapmatrix(n2,n1) = wfovlp(psi(0,0,1,jst),psip)
          END IF
        END DO
        IF(tprint) THEN
            WRITE(6,'(4i4,6(1pg13.5))') n1,ist,mpsi(ist),nbuf(ist),
     &      spenrg(ist),deltaf(ist),hmatrix(n1,n1)+spenrg(ist),
     &      gapmatrix(n1,n1),SQRT(wfovlp(psin,psin)-hmatrix(n1,n1)**2),
     &      SQRT(SUM(hmatrix(1:nmaxa,n1)**2)),hmatrix(n1,n1)
        END IF
        hmatrix(n1,n1)= hmatrix(n1,n1)+spenrg(ist)
        cm(n1,n1) = 1D0
          spenat(n1) = hmatrix(n1,n1)  !+efermiso
          gapnat(n1) = gapmatrix(n1,n1)
      END DO
      call sortenerg(nmaxa,spenat,iorder(1+isonpl))     ! sorting to energy
      DO na=1,nmaxa   ! overrule sorting
        iorder(na+isonpl)=na
      END DO
      IF(tprint1) WRITE(6,'(a,5000(2i4,2x))') 'IORDER:',
     &     (na+isonpl,iorder(na+isonpl)+isonpl,na=1,nmaxa)
      IF(tnewv2) THEN
        IF(tprint0 .AND. iso==2) THEN
          WRITE(6,'(a,2(1pg13.5))') 'before PAIRDN:'
          WRITE(6,'(a,200(1pg13.5))') 'v2a:   ',v2a(1:nmaxa)
          WRITE(6,'(a,200(1pg13.5))') 'gapnat:',gapnat(1:nmaxa)
        END IF
        IF(pairstab.NE.0D0)
     &   CALL pair_stab(spenat,gapnat,v2a,wate(isonpl+1),paircut,
     &              epairstab,pairstab,nmaxa,nmaxa,iab,.NOT.tpairspace)
        IF(toccmix) gwa=v2a
        CALL pairdn(spenat,gapnat,v2a,wate(isonpl+1),paircut,partnum,
     &              efermiso,epair(iso),gap(iso),ggg(iso),efadd,
     &              cutfac,nst(iso),nst(iso),npiter,ipair,0,istblo(iso),
     &              istcblock(iso),erra,itest-1,1,tpairspace)
        IF(toccmix) v2a=addocc*v2a+conocc*gwa
        IF(pairstab.NE.0D0) THEN
          WRITE(6,'(a,i3,2(1pg13.5))') '1. EPAIRSTAB:',
     &        iso,epair(iso),epairstab
          epair(iso) = epairstab
        END IF
        IF(tprint0) WRITE(6,'(a,6(1pg13.5))') 'after PAIRDN:',
     &              efermiso,epair(iso),gap(iso),ggg(iso),efadd,cutfac
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
        WRITE(6,'(a,3(1pg13.5),2i4)') 'test hermiticity:',
     &  SQRT(hermh/nmaxa**2),SQRT(hermg/nmaxa**2),hermhm,namax+isonpl,
     &  nbmax+isonpl
      END IF
      IF(.NOT.tpairspace) THEN
        uva(:)=SQRT((1D0-v2a(:))*v2a(:))*paircut(:)
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
          IF(iso==2) WRITE(6,'(a,3i4)') 
     &           'istbuf,nblo,nbhi=',nblo,nbhi,istbuf
          DO nb=nblo,nbhi
            IF(iso==2) WRITE(6,'(a,2i4)') 'NA,NB=',na,nb
            IF(nparz(1,nb+isonpl).NE.np1) CYCLE
            symcond(na,nb) = 
     &         +(v2a(nb)*hmatrix(na,nb)-v2a(na)*hmatrix(nb,na))
     &         +(uva(nb)*gapmatrix(na,nb)-uva(na)*gapmatrix(nb,na))
            symcondmax=MAX(symcondmax,ABS(symcond(na,nb)))
            avsymcond=avsymcond+symcond(na,nb)**2
          END DO 
        END DO
      ELSE
        DO na=1,nmaxa; DO nb=1,nmaxa
          symcond(na,nb) = 
     &         +(v2a(nb)*hmatrix(na,nb)-v2a(na)*hmatrix(nb,na))
     &         +(uva(nb)*gapmatrix(na,nb)-uva(na)*gapmatrix(nb,na))
          symcondmax=MAX(symcondmax,ABS(symcond(na,nb)))
          avsymcond=avsymcond+symcond(na,nb)**2
        END DO; END DO
      END IF
      symcondi(iso)=SQRT(avsymcond/partnum)
      WRITE(6,'(a,2i5,1pg13.5)') 'HFB natural: iso,iter,symm=',
     &         iso,0,SQRT(avsymcond/partnum)
!      WRITE(7,'(a,2i5,1pg13.5)') 'HFB natural: iso,iter,symm=',
!     &         iso,0,SQRT(avsymcond/partnum)
      IF(tprint0) THEN      
          WRITE(6,'(a,i5)') 'HFB natural: iter=',0
          WRITE(6,'(a,7(1pg13.5))') ' partnum,efermi,epair,cutfac=',
     &          partnum,SUM(v2a*phasfac)*2D0,efermiso,epair,cutfac
          IF(tprint1) THEN
            WRITE(6,'(a)')'# n  js  spe  v2  uv  cut   phasfac'
            DO na=1,nmaxa
            WRITE(6,'(2i3,5(1pg13.5))') na,jssav(na),spenat(na),v2a(na),
     &            gapnat(na),paircut(na),phasfac(na)
              nb=na+isonpl
              WRITE(6,'(6x,5(1pg13.5))') spenrg(nb),vsqf(nb),
     &            deltaf(nb),paircut(nb),wate(nb)
            END DO
          END IF
!          WRITE(6,'(a)') 'asymmetry matrix:'
!          DO na=1,nmaxa
!            WRITE(6,'(i3,50(1pg13.5))') na,symcond(:,na)
!          END DO
          WRITE(6,'(a,i2,4(1pg13.5))')
     &       'HFB-natural init.symm. condition: iso,max,average=',
     &        iso,symcondmax,symcondi(iso)
          IF(tprint) THEN
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

      IF(iternat.LE.0) RETURN
c
c     iterate
c
      IF(tprint0) THEN
        WRITE(6,'(a)') 'initial cm matrix:'
        DO na=1,nmaxa
          WRITE(6,'(i3,a,500(1pg13.5))') na,': ',cm(:,na)
        END DO
      END IF
      iab=-2
      IF(tprint1) iab=-1
      DO iter=1,iternat


c     h->psi in configuration space
        IF(tbuffer) THEN
          DO na=1,nmaxa
            IF(tHFBvariat) THEN
              fackincut = 
     &          -(1D0-paircut(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
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
              IF(nparz(1,n1+isonpl).NE.np1) CYCLE
                hcm(n1,na) =
     &             v2avar*SUM(hmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))+
     &             uva(na)*SUM(gapmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))
            END DO 
          END DO
        ELSE
          DO na=1,nmaxa
            ist=na+isonpl
            np1=nparz(1,ist)
            m1=mpsi(ist)
            IF(tHFBvariat) THEN
              fackincut = 
     &          -(1D0-paircut(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
              v2avar = v2a(na) + uva(na)*fackincut
            ELSE
              v2avar = v2a(na)
            END IF
            DO n1=1,nmaxa
              jst=n1+isonpl
              IF(mpsi(jst) /= m1 .OR. nparz(1,jst) /= np1)  CYCLE
              hcm(n1,na) = v2avar*SUM(hmatrix(n1,:)*cm(:,na))+
     &                     uva(na)*SUM(gapmatrix(n1,:)*cm(:,na))
            END DO
          END DO
        END IF
c     matrix of Lagrangian parameters and symmetry condition
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
              IF(nparz(1,nb+isonpl).NE.np1) CYCLE
                cbhca = SUM(cm(nblo:nbhi,nb)*hcm(nblo:nbhi,na))
                cahcb = SUM(cm(nblo:nbhi,na)*hcm(nblo:nbhi,nb))
                lambda(na,nb)= 0.5D0*(cbhca+cahcb)
                lambda(nb,na)=lambda(na,nb)
                symcond(na,nb)= (cbhca-cahcb)
                symcond(nb,na)= -symcond(na,nb)
            END DO 
          END DO
        ELSE
          DO na=1,nmaxa
            ist=na+isonpl
            np1=nparz(1,ist)
            m1=mpsi(ist)
            DO nb=1,na
              jst=nb+isonpl
              IF(mpsi(jst) /= m1 .OR. nparz(1,jst) /= np1)  CYCLE
              cbhca = SUM(cm(:,nb)*hcm(:,na))
              cahcb = SUM(cm(:,na)*hcm(:,nb))
              lambda(na,nb)= 0.5D0*(cbhca+cahcb)
              lambda(nb,na)=lambda(na,nb)
              symcond(na,nb)= (cbhca-cahcb)
              symcond(nb,na)= -symcond(na,nb)
            END DO
          END DO
        END IF
        IF(tprint0) WRITE(6,'(a,1pg13.5)') 'before step:',
     &    SQRT(SUM(symcond(:,:)**2)/partnum)
!        nitsub=1
!        DO itsub=1,nitsub
c     gradient step
        spemin=MINVAL(spenat,nmaxa)
        spemax=MAXVAL(spenat,nmaxa)
        delstep=delstepbas/(spemax-spemin)
        IF(tbuffer) THEN
          DO na=1,nmaxa
            IF(tHFBvariat) THEN
              fackincut = 
     &          -(1D0-paircut(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
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
              IF(nparz(1,n1+isonpl).NE.np1) CYCLE
              IF(e0bas==0D0) THEN
                dampstep = delstep
              ELSE
                dampstep = delstepbas/(hmatrix(n1,n1)-spemin+e0bas)
 !    &                  /nitsub
              END IF
              cmnew(n1,na) = cm(n1,na) - dampstep*(hcm(n1,na)
     &       -  SUM(cm(n1,nblo:nbhi)*lambda(nblo:nbhi,na))     )
!     &          v2avar*SUM(hmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))
!     &         +uva(na)*SUM(gapmatrix(n1,nblo:nbhi)*cm(nblo:nbhi,na))
!     &       -  SUM(cm(n1,nblo:nbhi)*lambda(nblo:nbhi,na))     )
            END DO 
          END DO
        ELSE
          DO na=1,nmaxa; 
            ist=na+isonpl
            np1=nparz(1,ist)
            m1=mpsi(ist)
            IF(tHFBvariat) THEN
              fackincut = 
     &          -(1D0-paircut(na))*gapmatrix(na,na)/(-0.1D0*cutfac)
              v2avar = v2a(na) + uva(na)*fackincut
            ELSE
              v2avar = v2a(na)
            END IF
            DO n1=1,nmaxa
              jst=n1+isonpl
              IF(mpsi(jst) /= m1 .OR. nparz(1,jst) /= np1)  CYCLE
              IF(e0bas==0D0) THEN
                dampstep = delstep
              ELSE
                dampstep = delstepbas/(hmatrix(n1,n1)-spemin+e0bas)
 !    &                  /nitsub
              END IF
              cmnew(n1,na) = cm(n1,na) - dampstep*(
     &          v2avar*SUM(hmatrix(n1,:)*cm(:,na))
     &         +uva(na)*SUM(gapmatrix(n1,:)*cm(:,na))
     &         -SUM(cm(n1,:)*lambda(:,na))     )
            END DO 
          END DO
        END IF

c     ortho-normalize
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
              IF(nbuf(jst).NE.istbuf .OR. nparz(1,jst).NE.np1) CYCLE
              ovlab=SUM(cm(nblo:nbhi,nb)*cm(nblo:nbhi,na))
              cm(nblo:nbhi,na)=cm(nblo:nbhi,na)-cm(nblo:nbhi,nb)*ovlab
            END DO
            ovlaa=1D0/SQRT(SUM(cm(nblo:nbhi,na)*cm(nblo:nbhi,na)))
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
              IF(mpsi(jst).NE.m1 .OR. nparz(1,jst).NE.np1)  CYCLE
              ovlab=SUM(cm(:,nb)*cm(:,na))
              cm(:,na)=cm(:,na)-cm(:,nb)*ovlab
            END DO
            ovlaa=1D0/SQRT(SUM(cm(:,na)*cm(:,na)))
            cm(:,na)=cm(:,na)*ovlaa
          END DO
        END IF

!        END DO          ! end of sub-iterations
c     new energies and occupations
        DO na=1,nmaxa
          ist=na+isonpl
          np1=nparz(1,ist)
          m1=mpsi(ist)
          spenat(na) = 0D0
          gapnat(na) = 0D0
          DO n1=1,nmaxa
            jst=n1+isonpl
            IF(mpsi(jst) /= m1 .OR. nparz(1,jst) /= np1)  CYCLE
            spenat(na) =spenat(na)+cm(n1,na)*SUM(hmatrix(n1,:)*cm(:,na))
            gapnat(na) = gapnat(na)+
     &                   cm(n1,na)*SUM(gapmatrix(n1,:)*cm(:,na))
          END DO
        END DO
        IF(pairstab.NE.0D0)
     &   CALL pair_stab(spenat,gapnat,v2a,wate(isonpl+1),paircut,
     &              epairstab,pairstab,nmaxa,nmaxa,iab,.NOT.tpairspace)
        IF(toccmix) gwa=v2a
        CALL pairdn(spenat,gapnat,v2a,wate(isonpl+1),paircut,partnum,
     &              efermiso,epair(iso),gap(iso),ggg(iso),efadd,
     &              cutfac,nst(iso),nst(iso),npiter,ipair,0,istblo(iso),
     &              istcblock(iso),erra,itest-1,1,tpairspace)
        IF(toccmix) v2a=addocc*v2a+conocc*gwa
        IF(pairstab.NE.0D0) THEN
!          WRITE(6,'(a,i3,2(1pg13.5))') '2. EPAIRSTAB:',
!     &        iso,epair(iso),epairstab
          epair(iso) = epairstab
        END IF
        IF(.NOT.tpairspace) THEN
          uva(:)=SQRT((1D0-v2a(:))*v2a(:))*paircut(:)
        ELSE
          uva(:)=SQRT((1D0-v2a(:))*v2a(:))
        END IF
        avsymcond=0D0
        symcondmax=0D0
        DO na=1,nmaxa
          avsymcond=avsymcond+SUM(symcond(:,na)**2)
          symcondmax=MAX(symcondmax,MAXVAL(ABS(symcond(:,na))))
        END DO
        IF(tprint0) THEN      
          WRITE(6,'(a,i5)') 'HFB natural: iter=',iter
          WRITE(6,'(a,6(1pg13.5))') 'partnum,eferm,epair,delst,symm=',
     &          partnum,efermi,epairact,delstep,
     &          symcondmax,SQRT(avsymcond/partnum)
         END IF
         IF(tprint1) THEN      
          WRITE(6,'(a)')'# n  js  spe  v2  uv  cut   phasfac'
          DO na=1,nmaxa
            WRITE(6,'(2i3,5(1pg13.5))') na,jssav(na),spenat(na),v2a(na),
     &            uva(na),paircut(na),phasfac(na)
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
      WRITE(6,'(a,3(1pg13.5))')
     &     'HFB natural exit: n1  spenew  gap'
        DO n1=1,nmaxa
          WRITE(6,'(i4,3(1pg13.5))')
     &     n1,hmatrix(n1,n1),gapmatrix(n1,n1)
        END DO
      END IF
      symcondf(iso)=SQRT(avsymcond/partnum)
      WRITE(6,'(a,i2,4(1pg13.5))')
     &       'HFB final symmetry condition: iso,max,average=',
     &        iso,symcondmax,symcondf(iso)
!      CALL  hfb_nat_densit(iso,nmaxa,cm,v2a,uva,jssav,phasfac)
      IF(tprint0) THEN
        WRITE(6,'(a/500(1pg13.5))') 'final occuptions:',v2a(:)
        WRITE(6,'(a)') 'final cm matrix:'
        DO na=1,nmaxa
          WRITE(6,'(i3,a,500(1pg13.5))') na,': ',cm(:,na)
        END DO
      END IF
      CALL  hfb_nat_mapwf(iso,nmaxa,cm,v2a,paircut,jssav,phasfac)

      DEALLOCATE(lambda,cmnew,gapnat,hmatrix,gapmatrix)
      DEALLOCATE(paircut,phasfac,spenat,jssav,v2a,uva,cm,gwa)
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
      END
c-----hfb_nat_mapwf---------------------------------------------------------
c
      subroutine hfb_nat_mapwf(iso,nmaxa,cm,v2a,paircut,jssav,phasfac)
c
      include 'paramr.inc'
c
c     Map new wavefunctions from configuration space back to coordinate space.
c     Version for 2D. Misuse attay 'dpsiz' temporarily as workspace.
c
      include 'wavefr.inc'
      include 'dynpar.inc'
      include 'levels.inc'
c
      INTEGER, INTENT(IN) :: nmaxa,iso,jssav(nmaxa)
      REAL(8), INTENT(IN) :: cm(nmaxa,nmaxa),v2a(nmaxa),paircut(nmaxa)
      REAL(8), INTENT(IN) :: phasfac(nmaxa)

      LOGICAL,PARAMETER :: tprint0=.FALSE.,tprint=.FALSE.
c
c      dimension zwj(21)
c
c-------------------------------------------------------------------------------
c
c
      IF(iso==1) THEN
        IF(nmaxa.NE.nst(1)) STOP 'HFB_nat_densit: wrong # of states'
        isonpl=0
      ELSE IF(iso==2) THEN
        IF(nmaxa.NE.nst(2)) STOP 'HFB_nat_densit: wrong # of states'
        isonpl=nst(1)
      ELSE
        STOP 'HFB: ISO must be 1 or 2'
      END IF
c
c     new wave functions from HFB coefficients
c
      acctot=0D0
      dpsiz=0D0
      DO n1=1,nmaxa
        ist=n1+isonpl
        DO n2=1,nmaxa
          jst=n2+isonpl
          dpsiz(0:nz,0:nr,1:2,ist) = dpsiz(0:nz,0:nr,1:2,ist)
     &                   + cm(n2,n1)*psi(0:nz,0:nr,1:2,jst)
        END DO
        vsqf(ist)=v2a(n1)
        wtpair(ist)=paircut(n1)
      END DO
c
c     move to correct place
c
      DO n1=1,nmaxa
        ist=n1+isonpl
        psi(0:nz,0:nr,1:2,ist) = dpsiz(0:nz,0:nr,1:2,ist)
      END DO
c

      RETURN
      END
c-----pair_stab------------ part of the hartree-fock package -----------------
c  
c     modify gap according to stabilized pairing
c
      subroutine pair_stab(e,delta,v2occ,ph,pairgw,epair,
     &                     pairstab,nmax,ndim,iab,tphpair)
c
c     input via list:
c      e       = array(1..nmax) of single particle energies.
c      delta   = array(1..nmax) of gap parameters.
c      v2occ   = array(1..nmax) of occupation probabilities,
c      ph      = array(1..nmax) of the multiplicities of the single
c                particle levels.
c                This multiplicities excludes (!) the spin degeneracy,
c                e.g. a 1s1/2 state has ph = 1 .
c      pairgw  = space for dynamical cut-off factor
c      pairstab  = stabilizing energy
c      nmax    = maximum number of single particle states
c                (needs to be smaller than the parameter 'ndim')
c      ndim    = dimensioning of the arrays, particularly of 'vpairf',
c                in the calling routine
c      iab     = print level; significant values are -2, -1, 0 and 1;
c                larger value prints more.
c      tphpair = switch to old pairing recipe
c     output via list:
c      delta   = array(1..nmax) of gap parameters.
c      epair   = pairing energy
c
      include 'paramr.inc'
      REAL(8) :: e(ndim),ph(ndim),delta(ndim),pairgw(ndim)
      REAL(8) :: v2occ(ndim)
      LOGICAL :: tphpair
      LOGICAL,PARAMETER :: tprint=.FALSE.
c
!      IF(pairstab==0D0) RETURN
c
      IF(tprint) THEN
        WRITE(6,'(a,i4,200(1pg13.5))') 'PAIR_STAB entry: nmax=',nmax
        WRITE(6,'(a,200(1pg13.5))') ' delta: ',delta(1:nmax)
        WRITE(6,'(a,200(1pg13.5))') ' v2occ: ',v2occ(1:nmax)
        WRITE(6,'(a,200(1pg13.5))') 'pairgw: ',pairgw(1:nmax)
        WRITE(6,'(a,200(1pg13.5))') '    ph: ',ph(1:nmax)
      END IF
c
c     compute pairing energy
c
      sumduv = 0D0
      sumv2 = 0D0
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
         sumv2 = ph(na)*v2 + sumv2
      ENDDO
C
!      epair  = -sumduv/2D0
      epair  = sumduv/2D0
!      if(iab.gt.1) write(6,'(a,3g12.4)') 'PAIR_STAB epair:',epair

      IF(pairstab==0D0) RETURN
      IF(tprint) write(6,'(a,2(3g12.4))') 'PAIR_STAB epair:',epair,sumv2
c
c     modify gaps
c
      DO na=1,nmax
        delta(na) = delta(na)*(1D0+(pairstab/epair)**2)
        delta(na) = SIGN(MIN(ABS(delta(na)),5D0),delta(na))
      ENDDO

c
      RETURN
      END


