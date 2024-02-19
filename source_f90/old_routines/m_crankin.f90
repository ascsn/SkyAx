MODULE cranking

USE paramr
USE steputil
USE linresponse, ONLY : psip,psin,ampsi,wrk1,wrk2,wrk3,wrk4
IMPLICIT REAL(8) (a-h,o-z)

LOGICAL,PARAMETER :: ktestm=.FALSE.   !<Switch to test "+-m" storage using m=0 mode.
LOGICAL,PARAMETER :: tprecpri=.false.     !<High precision output for alignments.
INTEGER,PARAMETER :: kcrank=1            !<Switch for cranking part of code, neede: value=1.
INTEGER,PARAMETER :: kintpl=1         !<Switch to writing file for later interpolation.

INTEGER,PARAMETER :: nzdc=kcrank*nzd,  nrdc=kcrank*nrd
INTEGER,PARAMETER :: ktrian=((npsid*npsid+npsid)/2) ! space for triang. storage
INTEGER,PARAMETER :: ktriac=kcrank*ktrian+1

!     fields for calculation the SC cranking mass
!      dpsir   = storage for d/dr derivative of wave functions
!c      dpsiz   = storage for d/dz derivative of wave functions
!      curqr   = respsonse current in r-direction
!      curqz   = respsonse current in z-direction
!      curqph  = respsonse current in phi-direction (only m.neq.0)
!      spinr   = respsonse spin-densiity in r-direction (only m.neq.0)
!      spinz   = respsonse spin-densiity in z-direction (only m.neq.0)
!      spinph  = respsonse spin-densiity in phi-direction
!      dynpoz  = response potential 'A' in z-direction
!      dynpor  = response potential 'A' in r-direction
!      dynpop  = response potential 'A' in phi-direction
!      spipor  = spin response potential 'S' in r-direction
!      spipoz  = spin response potential 'S' in z-direction
!      spipop  = spin response potential 'S' in phi-direction
!      rhop    = P-response density
!      taup    = P-response kinetic density
!      drhlsp  = P-response spin-orbit density
!      paircp  = P-response pair current
!      e2quas  = two quasi-particle energies
!      iposmx  = maximum number of matrix elements

!COMMON  /crankw/ 
REAL(8) :: dpsir(0:nzdc,0:nrdc,nsd,npsid), &
curqr(0:nzdc,0:nrdc,kqm), curqz(0:nzdc,0:nrdc,kqm),  &
    curqph(0:nzdc,0:nrdc,kqm), spinr(0:nzdc,0:nrdc,kqm),  &
    spinz(0:nzdc,0:nrdc,kqm), spinph(0:nzdc,0:nrdc,kqm),  &
    dynpor(0:nzdc,0:nrdc,kqm), dynpoz(0:nzdc,0:nrdc,kqm),  &
    dynpop(0:nzdc,0:nrdc,kqm), spipor(0:nzdc,0:nrdc,kqm),  &
    spipoz(0:nzdc,0:nrdc,kqm), spipop(0:nzdc,0:nrdc,kqm),  &
    rhop(0:nzdc,0:nrdc,kqm), taup(0:nzdc,0:nrdc,kqm),  &
    drhlsp(0:nzdc,0:nrdc,kqm), paircp(0:nzdc,0:nrdc,kqm),  &
    qresp(ktriac),presp(ktriac),ptest(ktriac),qab11(ktriac), &
    e2quas(ktriac)
REAL(8) :: sphc(npsid)      !<Intermediate storage for s.p. energies.
REAL(8),TARGET :: qab01(ktriac),qab10(ktriac)    ! ??
REAL(8) :: betafullp   ! full beta**2 including variance
INTEGER :: iposmx
LOGICAL :: tcranking=.FALSE.  ! switch to activate cranking
LOGICAL :: tatdhf=.FALSE.     ! switch to ATDHF path a posteriori

REAL(8) :: dpsiz(0:nzd,0:nrd,nsd,npsid)                !< Storage for gradient step also used for 
                                                       !!z-derivative in m_crankin.f90.
                                                       

LOGICAL :: tmass    !<
REAL(8) :: betawp   !< 
REAL(8) :: crmaso   !<
REAL(8) :: crmasp   !<
REAL(8) :: crmasq   !<
REAL(8) :: scmass   !<
REAL(8) :: pnmass   !<
REAL(8) :: pnmas0   !<
REAL(8) :: pwid     !<
REAL(8) :: pwido    !<
REAL(8) :: ammass   !<
REAL(8) :: amwid    !<
REAL(8) :: ampnms   !<

CONTAINS


SUBROUTINE all_cranking()
!---------------------------------------------------------------------------  
! DESCRIPTION: all_cranking
!> @brief
!!Central interface to 'skyax' main program. Collects all calls to
!!cranking and linear response routines.

USE linresponse, ONLY: linresp,linprint,ini_linrespfiles
USE observ
USE dynpar, ONLY: errenr,errvar,outnam,title
USE constr, ONLY: qerr,tswicon,swicon2
USE levels, ONLY: nprot,nneut


IMPLICIT NONE

LOGICAL :: tcrankfirst=.TRUE.

NAMELIST /cranking/ tcranking,tatdhf

! 
! First call: read conditions, check consistency, and prepare files
!
  IF(tcrankfirst) THEN
    READ(5,cranking) 
    tcrankfirst=.FALSE.
    IF(.NOT.tcranking) RETURN

    IF(nx3*ny3*nz3 > 1)  &
        STOP ' 3D writing incompatible with cranking'
    IF(swicon2<0.5D0) STOP 'Cranking only with L=2 constraint'
!    IF(ifcran == 2 .AND. tswicon)  &
!        STOP ' Separable RPA based on ground state calculation '

    OPEN(75,STATUS='unknown', FORM='formatted',FILE='crankn.'//outnam)
    WRITE(75,'(a/a,2i4/a,i1,3a)') title(1:78),  &
        ' A,Z=',nprot+nneut,nprot,  &
        ' beta_',2,'   SC_mass    Ing.mass   In.ms.-od  ',  &
        ' a.m. mass    a.m. wid  beta_wid    energy    cm_corr',  &
        '  PNmass     PNMas0    ZPE_triax   PNmas_a.m.'
    CLOSE(75)
    OPEN(175,STATUS='unknown', FORM='formatted',FILE='crank2.'//outnam)
    WRITE(175,'(a/a,2i4/a,i1,3a)') title(1:78),  &
        ' A,Z=',nprot+nneut,nprot,  &
        ' beta_',2,'   SC_mass    Ing.mass   In.ms.-od  ',  &
        ' a.m. mass    a.m. wid  beta_wid    energy    cm_corr',  &
        '  PNmass     PNMas0    ZPE_triax   PNmas_a.m.'
    CLOSE(175)
    
    OPEN(76,STATUS='unknown', FORM='formatted',FILE='zpecor.'//outnam)
    WRITE(76,'(a/a,2i4/a,i1,3a)') title(1:78),  &
        ' A,Z=',nprot+nneut,nprot,  &
        ' beta_',2,'   ZPE_SC    ZPE_beta   ZPE_beta-od    ZPE_rot',  &
        '    a.m. wid   E_ZPE_cm   ZPE_triax  omega_tri  r**2_charge',  &
        '  beta_p   rdms   surf   full-beta full-beta_p'
    CLOSE(76)
    
    
    OPEN(176,STATUS='unknown', FORM='formatted',FILE='zpeco2.'//outnam)
    WRITE(176,'(a/a,2i4/a,i1,3a)') title(1:78),  &
        ' A,Z=',nprot+nneut,nprot,  &
        ' beta_',2,'   ZPE_SC    ZPE_beta   ZPE_beta-od    ZPE_rot',  &
        '    a.m. wid   E_ZPE_cm   ZPE_triax  omega_tri  r**2_charge',  &
        '  beta_p   rdms   surf   full-beta full-beta_p'
    CLOSE(176)

    OPEN(77,STATUS='unknown', FORM='unformatted',FILE='intpol.'//outnam)

    IF (swicon2>0.5D0) THEN
      CALL ini_linrespfiles()
    END IF
    
  END IF

!
! The cranking parts
!
  IF(.NOT.tcranking) RETURN
  
  tmass =  econ <= 20.0*errenr  &
          .AND. devh2 < 20.0*errvar .AND. qcon <= 10.0*qerr
  tmass=.TRUE.

  IF(tmass) THEN
    CALL linresp()
    CALL crankin()
    CALL crnkam()
    CALL trizpe
    CALL linprint()
!        IF(tmass) THEN
!          OPEN(175,STATUS='unknown',FORM='formatted',&
!                  POSITION='append',FILE='crank2.'//outnam)
!          IF(cpbeta(0).NE.0D0 .AND. b(0,3).NE.0D0) &
!            WRITE(175,'(1x,f9.5,1pg20.12,5(1pg12.4),0pf13.5,0pf7.3,5(1pg12.4))')  &
!           cpbeta(0),0.5D0/b(0,3),0.5D0/b(0,2),0D0,0.5D0/b(1,3), &
!           2D0*pvar(1,3),pvar(0,3),  &
!           cpenergy(0),cpzpemod(0),2D0*bn(0,3),2D0*bn(0,2),zpetri,&
!           bn(1,2),bn(1,3)!/b(1,3)**2
!/b0q(0,0)**2/1.5D0
!          CLOSE(175)
!
!          OPEN(176,STATUS='unknown',FORM='formatted',&
!                  POSITION='append',FILE='zpeco2.'//outnam)
!          IF(cpbeta(0).NE.0D0 .AND. b(0,3).NE.0D0) &
!          WRITE(176,'(1x,f9.5,1pg20.12,2(1pg12.4),1pg20.12,1pg12.4,  &
!           1pg20.12,2(1pg12.4),0pf10.4,5f9.5)') cpbeta(0),b(0,3)*pvar(0,3), &
!           b(0,2)*pvar(0,2),0D0,2D0*b(1,3)*pvar(1,3),  &
!           2D0*pvar(1,3),cpenergy(0)-cpzpemod(0),zpetri,omgtri, &
!           rchrg2,betawp,rdmsc,surfc,betafull,betafullp
!          CLOSE(176)
!        END IF
    CALL crinfo()
    IF(.NOT.tswicon) THEN
      WRITE(21,'(a)') ' LINRESP: after'
      WRITE(22,'(a)') ' LINRESP: after'
      CALL linresp()
    END IF
!      ELSE IF(ifcran == 2) THEN
!        CALL seprpa                     ! presently disabled
  END IF


END SUBROUTINE all_cranking


!-----crankin-----------------------------------------------------------


SUBROUTINE crankin

!     saves wavefunctions
!     and computes Inglis cranking mass for quadrupole motion
!     by building derivative from previously save wavefunction

USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'


LOGICAL :: tderv              ! checks if derivative is available

!      real*4 e2quas(ktriac)     ! storage for 2-quasiparticle energies
REAL(8) :: betain2(ktriac)      !  make ALLOCATABLE

INTEGER,PARAMETER :: iqtest=0         ! switch to test given Q-operators
INTEGER,PARAMETER :: iptest=0         ! switch to test 'hcomP'
!      data itstcm/0/            ! switch to test case c.m. motion
!      data itstpr/1/            ! switch to testprints of Inglis Cranking
!      data isavsp/1/            ! save s.p. properties on file for interpol.

!      data delb/0.0/
!      save beto2,beto2p,beto2n,beto3,ferolp,feroln

!     buffered SQRT as internal statement function

!      sqrtb(x) = sqrt(max(x,smal))

!-----------------------------------------------------------------------

IF(kcrank <= 0) RETURN

CALL inie2q(0,1)                         ! 2QP energies on 'e2quas'

CALL quadms(qresp,0,1)                   ! simple quadrupole mass
CALL quadms(betain2,0,1)

INQUIRE(UNIT=98,OPENED=tderv)         ! checks status
CALL cringl(presp)                    ! prepare 'presp'
IF(.NOT.tderv) RETURN

!     first a few tests

IF(iqtest == 1) THEN
  DO icase=1,6
    
!     with cases
!     icase = 1...6
!     1: Q = 2z^2
!     2: Q = -r^2
!     3: Q = 2z^2-r^2
!     4: Q = z
!     5: Q = sigma_r
    
!     6: Q = p_r sigma_r
    
    b4=0.5*t4
    CALL hcot(icase)
    t1sav = t1
    t2sav = t2
    x1sav = x1
    x2sav = x2
    b4sav = b4
    b4psav = b4p
    CALL oabtst(icase,betain2)
    t1 = 0.0
    t2 = 0.0
    x1 = 0.0
    x2 = 0.0
    b4 = 0.0
    b4p = 0.0
    CALL oabtst(icase,betain2)
    t1 = t1sav
    t2 = t2sav
    x1 = x1sav
    x2 = x2sav
    b4 = b4sav
    b4p = b4psav
  END DO
  STOP 'cranking test-stop'
END IF
WRITE(7,'(a,2(1pg12.4))')  &
    ' CRANKN: P-response initialized  <P**2>=',opovlp(presp,presp)

!     optionally test 'hcommP'

IF(iptest == 1 .AND. tderv) THEN
  CALL hcpt(betain2)
!         stop ' after HCPT'
END IF

!     initialize Q-response with PRESP

!      write(7,'(a)') ' ipostr   ist   jst    equasij  deltai  deltaj',
!     &   '  eferm   spei  spej'
!test      ipostr = 0
!test      do ist=1,nst(1)+nst(2)
!test        mist   = mpsi(ist)
!test        npist  = nparz(1,ist)
!test        if(ist.le.nst(1)) then
!test          jstmax = nst(1)
!test        else
!test          jstmax = nst(1)+nst(2)
!test        endif
!test        do jst=ist,jstmax
!test          if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
!test            ipostr = 1 + ipostr
DO ipostr=1,iposmx
  qresp(ipostr) = presp(ipostr)/e2quas(ipostr)
END DO
!test          endif
!test        enddo
!test      enddo
WRITE(7,'(a,2(1pg12.4))')  &
    ' CRANKN: Q-response initialized  <Q**2>=',opovlp(qresp,qresp)
!      iposmx = ipostr
!      write(7,'(a,i10)') ' second IPOSTR=',ipostr

!     iterate    Q_new = Q - stepQ * ( [H,Q] - P ),
!     the commutator [H,Q] is computed in 'hcommQ' and returned on 'ptest'.
!     The collective mass <[Q,[H,Q]]> is also accumulated in 'hcommQ' and
!     returned on 'scmass'. It is printed as protocol of convergence.

stepq  =   0.05                !  this needs to become an input variable
itmaxq =   150                !  this needs to become an input variable
epsq   = 0.0001
scmaso = 0.0
DO itq=1,itmaxq
  CALL hcommq(qresp,presp,betain2,ptest,0)
  IF(scmass > zero) pnmass = pnmass/scmass**2
  IF(itq == 1) pnmas0 = pnmass                      ! save Inglis level
  DO ipostr=1,iposmx
    qresp(ipostr) = qresp(ipostr) - stepq*(ptest(ipostr)-presp(ipostr))  &
        /e2quas(ipostr)                  ! to damp step
  END DO
  WRITE(7,'(a,i4,a,4(1pg12.4))') ' LRiter=',itq,  &
      '  SCmass,diff.,PNmasses=',scmass,scmass-scmaso,pnmass,pnmas0
  IF(itq > 10 .AND. scmass /= zero) THEN
    IF(ABS(scmaso/scmass-one) < epsq .OR. scmass < zero  &
        .OR. scmass > 1.0E20) GO TO 999
!          continue
  ELSE IF(scmass == zero) THEN
    GO TO 999
  END IF
  scmaso = scmass
END DO
999  CONTINUE
WRITE(7,'(a,i5)') ' mass iteration ended at itq=',itq

RETURN
END SUBROUTINE crankin




!-----crnkam-----------------------------------------------------------

SUBROUTINE crnkam

!     and computes Inglis cranking mass for angular momentum

USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
USE linresponse, ONLY: angmomop,psi1,psi2
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'


INTEGER,PARAMETER :: itstpr=0       ! switch to testprints of Inglis Cranking
LOGICAL,PARAMETER :: tangmomsub=.TRUE.
LOGICAL,PARAMETER :: tprintdetails=.FALSE.

!     buffered SQRT as internal statement function

!sqrtb(x) = SQRT(MAX(x,smal))

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
  IF(tangmomsub) THEN
    CALL angmomop(psi(0,0,1,ist),psip,psin,nparz(1,ist),mpsi(ist))  
  ELSE
!                                     first, the L+ and L-
  CALL nablaz_spin(psi(0:nz,0:nr,1:2,ist),psi1(0:nz,0:nr,1:2),nparz(1,ist),1,1D0)
  CALL nablar_spin(psi(0:nz,0:nr,1:2,ist),psi2(0:nz,0:nr,1:2),   mpsi(ist),1,1D0)
  
  DO is=1,nsd
    mm        = mist+is-1
    npist  = nparz(is,ist)
    nparr     = 1+MOD(mm,2)
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
!  angular momentum operator finished
   END IF

   ampsi(:,:,:,ist,1) = psip
   ampsi(:,:,:,ist,2) = psin
   betampsi=beta2

!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  deltai = deltaf(ist)*wtpair(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  tstovl = wfovlp(psip(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2))
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
        pij=wfovlptr(psi(0:nz,0:nr,1:2,jst),psin(0:nz,0:nr,1:2))
        isw=1
      ELSE IF(mist+1 == mpsi(jst)) THEN
        pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip(0:nz,0:nr,1:2))
        isw=3
      ELSE IF(mist-1 == mpsi(jst)) THEN
        pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psin(0:nz,0:nr,1:2))
        isw=2
      END IF
      IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
          ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
!PK               pij    = pij*(vj*ui-uj*vi)
      IF(tprintdetails .AND. ABS(vj*ui-uj*vi)>1D-3) &
         WRITE(993,'(a,f9.6,5i4,4(1pg12.4))')  &
          'b2,i,j,phase,pij=',beta2,ist,jst,mpsi(ist),mpsi(jst),isw,&
               -(vj*ui-uj*vi),pij,(pij*(vj*ui-uj*vi))**2
      pij    = pij*(vj*ui-uj*vi)
      pij    = pij*pij
      sumj2e = pij/(equasi+equasj) + sumj2e
      sumj2  = pij + sumj2
      ampnms = pij/(equasi+equasj)**2 *(ui*ui-vi*vi+uj*uj-vj*vj)  &
          + ampnms
!      IF((vj*ui-uj*vi)**2>1D-5) WRITE(991,'(a,f9.5,4i4,2(1pg13.5))') &
!           'CRNKAM:',beta2,ist,jst,mist,mpsi(jst),pij,(vj*ui-uj*vi)**2
!      CALL flush(991)
    END IF
  END DO
END DO
ammass = sumj2e/2.0
amwid  = sumj2/2.0
ampnms = ampnms/2.0/ammass**2
WRITE(7,'(a,1pg12.4)') ' test <J**2>=',sumtst
WRITE(993,'(a,2(1pg12.4))') ' test <J**2>=',sumtst,sumj2
IF(tprintdetails) CALL flush(993)

RETURN
END SUBROUTINE crnkam

!-----hcommP-----------------------------------------------------------

SUBROUTINE hcommp(popin,qopout)

!     Commutator [H,P] of Hamiltonian with P-response matrix,
!     matrix elements of P are entered on 'popin',
!     result is returned on 'qopout'.
!     These matrices carry the matrix elements of the operator as
!     two-quasiparticle excitations. They are stored in upper
!     triangular form.

!     The routine uses the array 'dpsir' in 'common /crankin/' as
!     workspace.

!     This version computes the variation of the response Hamiltonian
!     directly from the potential by finite differences.

USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)


REAL(8), INTENT(IN)                         :: popin(ktriac)
REAL(8), INTENT(OUT)                        :: qopout(ktriac)

REAL(8),PARAMETER :: stepq0=0.000001D0   ! determines step size --> 'stepr'

!     internal statement function for extrapolation to r=0
REAL(8),PARAMETER :: trhalf=1.5D0,trfift=0.6D0,oneten=0.1D0

REAL(8),ALLOCATABLE :: potp(:,:,:),coulptp(:,:), &
                       bmassp(:,:,:),potlsp(:,:,:,:)


ALLOCATE(potp(0:nzdc,0:nrdc,kqm), coulptp(0:nzdc,0:nrdc),  &
    bmassp(0:nzdc,0:nrdc,kqm), potlsp(0:nzdc,0:nrdc,kqm,2))



!-----------------------------------------------------------------------

IF(ifls /= 0) WRITE(*,*) ' IFLS=0 needed for linear response'
!      if(ifls.ne.0) stop ' IFLS=0 needed for linear response'
IF(ipair /= 5 .AND. ipair /= 9 .AND. ipair /= 11)  &
    WRITE(*,*) ' IPAIR=5 or 9 needed for linear P response'
!     &    stop ' IPAIR=5 or 9 needed for linear P response'

!     1. step: accumulate response currents

DO iq=1,2
  DO ir=0,nr
    DO iz=0,nz
      rhop(iz,ir,iq) = 0.0
      taup(iz,ir,iq) = 0.0
      drhlsp(iz,ir,iq) = 0.0
      paircp(iz,ir,iq) = 0.0
    END DO
  END DO
END DO

stepr = 0.01/SQRT(opovlp(popin,popin))

ipostr = 0
DO ist=1,nst(1)+nst(2)
  vi     = SQRT(vsqf(ist))
  ui     = SQRT(one-vi*vi)
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    eferm  = fer(1)
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    eferm  = fer(2)
    jstmax = nst(1)+nst(2)
  END IF
  addfac = 2.0D0
  DO jst=ist,jstmax
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      xm2    = mist*mist
      xm2l   = (mist+1)*(mist+1)
      ipostr = 1 + ipostr
      vj     = SQRT(vsqf(jst))
      uj     = SQRT(one-vj*vj)
      factor = addfac*(ui*vj+vi*uj)*popin(ipostr)
      factpa = addfac*(ui*uj-vi*vj)*popin(ipostr)
      DO ir=1,nr
        r      = (sr+ir)*dr
        zentfu = xm2/(r*r)
        zentfl = xm2l/(r*r)
        DO iz=0,nz
          rhoacu =  psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
          rhoacl =  psi(iz,ir,2,ist)*psi(iz,ir,2,jst)
          rhop(iz,ir,iq) = factor*(rhoacu+rhoacl) + rhop(iz,ir,iq)
          paircp(iz,ir,iq) = factpa*(rhoacu+rhoacl) + paircp(iz,ir,iq)
          addtau = factor* (dpsir(iz,ir,1,ist)*dpsir(iz,ir,1,jst)  &
              +dpsir(iz,ir,2,ist)*dpsir(iz,ir,2,jst)  &
              +dpsiz(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)  &
              +dpsiz(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)  &
              +zentfu*rhoacu+zentfl*rhoacl         )
          wrk1(iz,ir) = addtau          ! save that for spin-orbit
          taup(iz,ir,iq) = addtau + taup(iz,ir,iq)
        END DO
      END DO
      IF(t4 /= zero) THEN
        DO ir=1,nr
          rinv  = one/((sr+ir)*dr)
          xmup  = mist*rinv
          xmlow = (mist+1)*rinv
          DO iz=0,nz
            spsiui = dpsiz(iz,ir,1,ist)+dpsir(iz,ir,2,ist)  &
                +xmlow*psi(iz,ir,2,ist)
            spsili = -dpsiz(iz,ir,2,ist)+dpsir(iz,ir,1,ist)  &
                -xmup*psi(iz,ir,1,ist)
            spsiuj = dpsiz(iz,ir,1,jst)+dpsir(iz,ir,2,jst)  &
                +xmlow*psi(iz,ir,2,jst)
            spsilj = -dpsiz(iz,ir,2,jst)+dpsir(iz,ir,1,jst)  &
                -xmup*psi(iz,ir,1,jst)
            addls  = factor*(spsiui*spsiuj+spsili*spsilj) -wrk1(iz,ir)
            drhlsp(iz,ir,iq) = -addls + drhlsp(iz,ir,iq)
          END DO
        END DO
      END IF
    END IF
    addfac = 4.0 ! 2.0
  END DO
END DO
!                                                 extrapolate to ir=0
DO iz=0,nz
  rhop(iz,0,1) = extrap(rhop(iz,1,1),rhop(iz,2,1),rhop(iz,3,1))
  rhop(iz,0,2) = extrap(rhop(iz,1,2),rhop(iz,2,2),rhop(iz,3,2))
  taup(iz,0,1) = extrap(taup(iz,1,1),taup(iz,2,1),taup(iz,3,1))
  taup(iz,0,2) = extrap(taup(iz,1,2),taup(iz,2,2),taup(iz,3,2))
  drhlsp(iz,0,1) = extrap(drhlsp(iz,1,1),drhlsp(iz,2,1), drhlsp(iz,3,1))
  drhlsp(iz,0,2) = extrap(drhlsp(iz,1,2),drhlsp(iz,2,2), drhlsp(iz,3,2))
END DO

!     2. step: shift the main density by the response densities,
!              save the main densities on the response fields as
!              workspace

rhopno = wfovlp(rhop(0:nz,0:nr,1:2),rhop(0:nz,0:nr,1:2))
rhopov = wfovlp(rhop(0:nz,0:nr,1:2),rho(0:nz,0:nr,1:2))
rhono  = wfovlp(rho(0:nz,0:nr,1:2),rho(0:nz,0:nr,1:2))
varrho = rhopno/rhono
!      varrho = wfovlp(rhop,rhop)/wfovlp(rho,rho)
taupno = wfovlp(taup(0:nz,0:nr,1:2),taup(0:nz,0:nr,1:2))
taupov = wfovlp(taup(0:nz,0:nr,1:2),tau(0:nz,0:nr,1:2))
tauno  = wfovlp(tau(0:nz,0:nr,1:2),tau(0:nz,0:nr,1:2))
vartau = taupno/tauno
!      vartau = wfovlp(taup,taup)/wfovlp(tau,tau)
stepr  = SQRT(MIN(stepq0/varrho,stepq0/vartau))
WRITE(7,'(a,3(1pg12.4))') ' HCOMMP: stepr,varrho,vartau=',stepr,varrho,vartau
WRITE(7,'(a/5(1pg12.4))')  &
    ' N(rho)  N(rhop)  Overlap   varrho  difference ',  &
    rhono,rhopno,rhopov,varrho,(rhono+rhopno-2.0*rhopov)/rhono
WRITE(7,'(a/5(1pg12.4))')  &
    ' N(tau)  N(taup)  Overlap   vartau  difference ',  &
    tauno,taupno,taupov,vartau,(tauno+taupno-2.0*taupov)/tauno
steph = 0.5*stepr
DO iq=1,2
  DO ir=0,nr
    DO iz=0,nz
      xsav   = rho(iz,ir,iq)
!        x      = sqrt(xsav)
!        rho(iz,ir,iq) = (x+steph*rhop(iz,ir,iq)/x)**2
      rho(iz,ir,iq) = xsav+stepr*rhop(iz,ir,iq)
      rhop(iz,ir,iq) = xsav
      xsav   = tau(iz,ir,iq)
      tau(iz,ir,iq) = xsav+stepr*taup(iz,ir,iq)
      taup(iz,ir,iq) = xsav
      xsav   = drhols(iz,ir,iq)
      drhols(iz,ir,iq) = xsav+stepr*drhlsp(iz,ir,iq)
      drhlsp(iz,ir,iq) = xsav
    END DO
  END DO
END DO
DO ir=0,nr
  DO iz=0,nz
    crho(iz,ir) = rho(iz,ir,1)
  END DO
END DO

!     test new density

rhono  = wfovlp(rhop(0:nz,0:nr,1:2),rhop(0:nz,0:nr,1:2))
rhopov = wfovlp(rhop(0:nz,0:nr,1:2),rho(0:nz,0:nr,1:2))
rhopno = wfovlp(rho(0:nz,0:nr,1:2),rho(0:nz,0:nr,1:2))
tauno  = wfovlp(taup(0:nz,0:nr,1:2),taup(0:nz,0:nr,1:2))
taupov = wfovlp(taup(0:nz,0:nr,1:2),tau(0:nz,0:nr,1:2))
taupno = wfovlp(tau(0:nz,0:nr,1:2),tau(0:nz,0:nr,1:2))
!      stepr  = sqrt(min(1.0E-4/varrho,1.0E-4/vartau))
WRITE(7,'(a/a/5(1pg12.4))') ' HCOMMP:',  &
    ' N(rho)  N(rhop)  Overlap   angle ',  &
    rhono,rhopno,rhopov,one-rhopov*rhopov/(rhono*rhopno)
WRITE(7,'(a/5(1pg12.4))') ' N(tau)  N(taup)  Overlap   angle ',  &
    tauno,taupno,taupov,one-taupov*taupov/(tauno*taupno)

!     3. step: compute Hamiltonian for shifted densities
!              and deduce response Hamiltonian from that.
!              Response fields reside on potential fields.

DO ist=1,nst(1)+nst(2)
  sphc(ist) = spenrg(ist)                        ! save s.p. energy
  spenrg(ist) = 0.0D0
END DO
DO ir=0,nr
  DO iz=0,nz
    coulptp(iz,ir) = coulpt(iz,ir)
    potp(iz,ir,1) = pot(iz,ir,1)
    potp(iz,ir,2) = pot(iz,ir,2)
    bmassp(iz,ir,1) = bmass(iz,ir,1)
    bmassp(iz,ir,2) = bmass(iz,ir,2)
    potlsp(iz,ir,1,1) = potls(iz,ir,1,1)
    potlsp(iz,ir,2,1) = potls(iz,ir,2,1)
    potlsp(iz,ir,1,2) = potls(iz,ir,1,2)
    potlsp(iz,ir,2,2) = potls(iz,ir,2,2)
  END DO
END DO
CALL mfhmlt(0)
sd = -1.0D0/stepr
DO ir=0,nr
  DO iz=0,nz
    coulpt(iz,ir) = (coulptp(iz,ir)-coulpt(iz,ir))*sd
    pot(iz,ir,1) = (potp(iz,ir,1)-pot(iz,ir,1))*sd
    pot(iz,ir,2) = (potp(iz,ir,2)-pot(iz,ir,2))*sd
    bmass(iz,ir,1) = (bmassp(iz,ir,1)-bmass(iz,ir,1))*sd
    bmass(iz,ir,2) = (bmassp(iz,ir,2)-bmass(iz,ir,2))*sd
    potls(iz,ir,1,1) = (potlsp(iz,ir,1,1)-potls(iz,ir,1,1))*sd
    potls(iz,ir,2,1) = (potlsp(iz,ir,2,1)-potls(iz,ir,2,1))*sd
    potls(iz,ir,1,2) = (potlsp(iz,ir,1,2)-potls(iz,ir,1,2))*sd
    potls(iz,ir,2,2) = (potlsp(iz,ir,2,2)-potls(iz,ir,2,2))*sd
  END DO
END DO

CALL nablaz_scalar(bmass(0:nz,0:nr,1),dzbmas(0:nz,0:nr,1), 1,1,one)
CALL nablaz_scalar(bmass(0:nz,0:nr,2),dzbmas(0:nz,0:nr,2), 1,1,one)
CALL nablar_scalar(bmass(0:nz,0:nr,1),drbmas(0:nz,0:nr,1), 0,1,one)
CALL nablar_scalar(bmass(0:nz,0:nr,2),drbmas(0:nz,0:nr,2), 0,1,one)

CALL test_print(pot(0:nz,0:nr,1),0,9,0,9,'resp pot')
CALL test_print(pot(0:nz,0:nr,2),0,9,0,9,'resp pot')
CALL test_print(bmass(0:nz,0:nr,1),0,9,0,9,'resp pot')
CALL test_print(bmass(0:nz,0:nr,2),0,9,0,9,'resp pot')

!     4. step: compute action of response Hamiltonian and
!              compose new Q operator

ipostr = 0
DO ist=1,nst(1)+nst(2)
  vi     = SQRT(vsqf(ist))
  ui     = SQRT(one-vi*vi)
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    jstmax = nst(1)+nst(2)
  END IF
  CALL hpsi(psi(0:nz,0:nr,1:2,ist),psin,ist)
  pairfa = -pairfc(iq)
  DO is=1,2
    DO ir=0,nr
      DO iz=0,nz
        psip(iz,ir,is) = pairfa*paircp(iz,ir,iq)*psi(iz,ir,is,ist)
      END DO
    END DO
  END DO
  DO jst=ist,jstmax
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj     = SQRT(vsqf(jst))
      uj     = SQRT(one-vj*vj)
      qopout(ipostr) = popin(ipostr)*e2quas(ipostr)  &
          +(ui*vj+vi*uj)*wfovlp(psi(0:nz,0:nr,1:2,jst),psin(0:nz,0:nr,1:2))
!     &              +(ui*uj-vi*vj)*wfovlp(psi(0,0,1,jst),psip)/stepr
    END IF
  END DO
END DO
!      write(6,'(a,i6)') ' HCOMMP: ipostr=',ipostr

!     4. step: restore densities and mean field Hamiltonian

!      write(7,'(a)') ' after 3.step in HCOMMP: ist, spenrg'
DO ist=1,nst(1)+nst(2)
  spenrg(ist) = sphc(ist)                        ! restore s.p. energy
!        write(7,'(1x,i5,1pg12.4)') ist,spenrg(ist)
END DO
DO iq=1,2
  DO ir=0,nr
    DO iz=0,nz
      rho(iz,ir,iq) = rhop(iz,ir,iq)
      tau(iz,ir,iq) = taup(iz,ir,iq)
      drhols(iz,ir,iq) = drhlsp(iz,ir,iq)
    END DO
  END DO
END DO
DO ir=0,nr
  DO iz=0,nz
    crho(iz,ir) = rho(iz,ir,1)
    coulpt(iz,ir) = coulptp(iz,ir)
    pot(iz,ir,1) =  potp(iz,ir,1)
    pot(iz,ir,2) =  potp(iz,ir,2)
    bmass(iz,ir,1) =  bmassp(iz,ir,1)
    bmass(iz,ir,2) =  bmassp(iz,ir,2)
    potls(iz,ir,1,1) =  potlsp(iz,ir,1,1)
    potls(iz,ir,2,1) =  potlsp(iz,ir,2,1)
    potls(iz,ir,1,2) =  potlsp(iz,ir,1,2)
    potls(iz,ir,2,2) =  potlsp(iz,ir,2,2)
  END DO
END DO
!      call mfhmlt(0)
DO ist=1,nst(1)+nst(2)
  spenrg(ist) = sphc(ist)
END DO

DEALLOCATE(potp, coulptp,bmassp, potlsp)

RETURN

END SUBROUTINE hcommp

REAL(8) PURE FUNCTION extrap(f1,f2,f3)
REAL(8),INTENT(IN) :: f1,f2,f3
!DATA trhalf,trfift,oneten/1.5,0.6,0.1/
!extrap(f1,f2,f3) = trhalf*f1-trfift*f2+oneten*f3

extrap = 1.5D0*f1-0.6D0*f2+0.1D0*f3

RETURN
END FUNCTION extrap



!-----opovlp--------------------------------------------------------------

REAL(8) PURE FUNCTION opovlp(op1in,op2in)

!     Overlap of 1ph-operators 'op1in' and 'op2in'.

USE wavefr
USE levels
IMPLICIT REAL(8) (a-h,o-z)



REAL(8), INTENT(IN)                         :: op1in(ktriac)
REAL(8), INTENT(IN)                         :: op2in(ktriac)


!-----------------------------------------------------------------------

!     accumulate matrix elements

sum    = 0.0
ipostr = 0
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  IF(ist <= nst(1)) THEN ! protons
    jstmax = nst(1)
  ELSE                   ! neutrons
    jstmax = nst(1)+nst(2)
  END IF
  DO jst=ist,jstmax
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      IF(ist == jst) THEN
        sum = 2.0 * op1in(ipostr)*op2in(ipostr)   + sum
      ELSE
        sum = 4.0 * op1in(ipostr)*op2in(ipostr)    + sum
      END IF
    END IF
  END DO
END DO
opovlp = sum

RETURN
END FUNCTION opovlp

!-----hcommQ-----------------------------------------------------------

SUBROUTINE hcommq(qopin,popin,betainp,popout,icase)
!SUBROUTINE hcommq(qopin,popin,betainp,popout,icase)

!     Commutator [H,Q] of Hamiltonian with Q-response matrix,
!     matrix elements of Q are entered on 'qopin',
  !     result is returned on 'popout'.
!     The collective mass <[Q,[H,Q]]> is also accumulated and
!     returned on 'scmass'.
!     The particle-number mass is accumulated on
!     'pnmass' = with linear response.
!     The routine requires that the derivatives of the wave functions
!     are already computed and stored in 'dpsir' and 'dpsiz'.

!     icase: for test-operatores, icase=0 normal use, for details on
!            icase see source of Oabtst() and hcOT()

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'


REAL(8), INTENT(IN)                         :: qopin(ktriac)
REAL(8), INTENT(IN)                         :: popin(ktriac)
REAL(8), INTENT(IN)                         :: betainp(ktriac)
REAL(8), INTENT(OUT)                        :: popout(ktriac)
INTEGER, INTENT(IN)                      :: icase



INTEGER,PARAMETER :: ifhcts=0             ! switch for testing
INTEGER,PARAMETER :: ifb4ad=1                 ! switch for b4&b4p terms

!     buffered SQRT as internal statement function
!sqrtb(x) = SQRT(MAX(x,smal))

!-----------------------------------------------------------------------

IF(ktestm) STOP ' old HCOMMQ needs KTESTM=.FALSE.'
IF(ifls /= 0) WRITE(6,*) ' Warning: IFLS=0 needed for linear response'
!      if(ifls.ne.0) stop ' IFLS=0 needed for linear response'

!     prepare derivatives for later use

DO ist=1,nst(1)+nst(2)
  CALL nablaz_spin(psi(0:nz,0:nr,1:2,ist),dpsiz(0:nz,0:nr,1:2,ist), nparz(1,ist),1,one)
  CALL nablar_spin(psi(0:nz,0:nr,1:2,ist),dpsir(0:nz,0:nr,1:2,ist), mpsi(ist),1,one)
END DO
!      write(7,'(a)') ' HCOMMQ: derivatives are prepared'

!     1. step: response currents are collected

DO ir=0,nr
  DO iz=0,nz
    curqr(iz,ir,1)  = 0.0
    curqr(iz,ir,2)  = 0.0
    curqz(iz,ir,1)  = 0.0
    curqz(iz,ir,2)  = 0.0
    spinph(iz,ir,1) = 0.0
    spinph(iz,ir,2) = 0.0
  END DO
END DO

ipostr = 0
DO ist=1,nst(1)+nst(2)
  vi     = sqrtb(vsqf(ist))
  ui     = sqrtb(one-vi*vi)
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    eferm  = fer(1)
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    eferm  = fer(2)
    jstmax = nst(1)+nst(2)
  END IF
  DO jst=ist,jstmax
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj     = sqrtb(vsqf(jst))
      uj     = sqrtb(one-vj*vj)
!R030399 check
!old            factor = (ui*ui*vj*vj+vi*vi*uj*uj)*qopin(ipostr)
      IF(ist /= jst) THEN ! exlude diagonal element ? yes, should be zero.
        factor = (ui*ui*vj*vj-vi*vi*uj*uj)*qopin(ipostr)
        DO ir=0,nr
          DO iz=0,nz
            
! factor 1/2 ?
            
!R030399 check factor ?  !0.5 ! j_Q,q;r 
            curqr(iz,ir,iq)  = -2.0       &
                *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)  &
                + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)  &
                - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)  &
                - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist)) *factor  +  curqr(iz,ir,iq)
          END DO
          DO iz=0,nz
!R030399 check   !0.5 ! j_Q,q;z
            curqz(iz,ir,iq)  = -2.0       &
                *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)  &
                + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)  &
                - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)  &
                - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist)) *factor  +  curqz(iz,ir,iq)
!                        spinph(iz,ir,iq) = 0.5 ! sigma_Q,q;z
!     &                       *(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
!     &                       - psi(iz,ir,1,ist)*psi(iz,ir,2,jst)
!     &                       - psi(iz,ir,2,jst)*psi(iz,ir,1,ist)
!     &                       + psi(iz,ir,1,jst)*psi(iz,ir,2,ist))
!     &                       *factor  +  spinph(iz,ir,iq)
          END DO
          DO iz=0,nz
!R030399 check ! sigma_Q,q;z 
            spinph(iz,ir,iq) = 2.0 &
                *(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)  &
                - psi(iz,ir,1,ist)*psi(iz,ir,2,jst)) *factor  +  spinph(iz,ir,iq)
          END DO
        END DO
      END IF
    END IF
  END DO
END DO
iposmx = ipostr

!rold      iact = min(kcranc,1)
!############################################################################
!R250599
! test output of currents

IF(icase /= 0) THEN
  CALL ptstcur(icase)
END IF
!R250599 end of currents test output
!############################################################################


!     2. step: response potentials are composed

b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
b1c = b1-b1p
IF(ifb4ad /= 1) THEN
  b4  = 0.0   ! t4*0.5
  b4c = 0.0   ! b4+b4p
ELSE
  b4  = t4*0.5
  b4c = b4+b4p
END IF

IF(isym == 1)THEN
  nparzs = 2
ELSE
  nparzs = 1
END IF

CALL nablaz_scalar(spinph(0:nz,0:nr,1),wrk1(0:nz,0:nr),nparzs,1,one)
CALL nablaz_scalar(spinph(0:nz,0:nr,2),wrk2(0:nz,0:nr),nparzs,1,one)
CALL nablar_scalar(spinph(0:nz,0:nr,1),wrk3(0:nz,0:nr),     1,1,one)
CALL nablar_scalar(spinph(0:nz,0:nr,2),wrk4(0:nz,0:nr),     1,1,one)

DO ir=0,nr
  DO iz=0,nz
    dynpor(iz,ir,1)  = -2.0*(b1c*curqr(iz,ir,1) + b1*curqr(iz,ir,2))  &
        +     b4c*wrk1(iz,ir)    + b4*wrk2(iz,ir)
    dynpor(iz,ir,2)  = -2.0*(b1c*curqr(iz,ir,2) + b1*curqr(iz,ir,1))  &
        +     b4c*wrk2(iz,ir)    + b4*wrk1(iz,ir)
    dynpoz(iz,ir,1)  = -2.0*(b1c*curqz(iz,ir,1) + b1*curqz(iz,ir,2))  &
        -     b4c*wrk3(iz,ir)    - b4*wrk4(iz,ir)
    dynpoz(iz,ir,2)  = -2.0*(b1c*curqz(iz,ir,2) + b1*curqz(iz,ir,1))  &
        -     b4c*wrk4(iz,ir)    - b4*wrk3(iz,ir)
  END DO
END DO
!     ! rot(j_Q,q)

CALL nablaz_scalar(curqr(0:nz,0:nr,1),wrk1(0:nz,0:nr), 1,1,one)
CALL nablaz_scalar(curqr(0:nz,0:nr,2),wrk2(0:nz,0:nr), 1,1,one)
CALL nablar_scalar(curqz(0:nz,0:nr,1),wrk3(0:nz,0:nr), 0,1,one)
CALL nablar_scalar(curqz(0:nz,0:nr,2),wrk4(0:nz,0:nr), 0,1,one)

!     ! S_Q,q
DO ir=0,nr
  DO iz=0,nz
    spipop(iz,ir,1)  = -b4c*(wrk1(iz,ir)-wrk3(iz,ir))  &
        - b4*(wrk2(iz,ir)-wrk4(iz,ir))
    spipop(iz,ir,2)  = -b4c*(wrk2(iz,ir)-wrk4(iz,ir))  &
        - b4*(wrk1(iz,ir)-wrk3(iz,ir))
  END DO
END DO

!     3. step: accumulate matrix elements by spatial integration with
!              transition currents

pnmass = 0.0
pqnorm = 0.0
scmass = 0.0
tstmas = 0.0              ! Inglis mass as test
ipostr = 0
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    eferm  = fer(1)
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    eferm  = fer(2)
    jstmax = nst(1)+nst(2)
  END IF
  deltai = deltaf(ist)*wtpair(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
  vi     = sqrtb(vsqf(ist))
  ui     = sqrtb(one-vi*vi)
  
  DO jst=ist,jstmax
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj     = sqrtb(vsqf(jst))
      uj     = sqrtb(one-vj*vj)
      deltaj = deltaf(jst)*wtpair(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
!     ! h_Q,alpha,beta
      DO ir=0,nr
        DO iz=0,nz
! sign? 
          wrk1(iz,ir) = -0.5                                      &
              *(psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)  &
              + psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)  &
              - psi(iz,ir,1,jst)*dpsir(iz,ir,1,ist)  &
              - psi(iz,ir,2,jst)*dpsir(iz,ir,2,ist)) * dynpor(iz,ir,iq) &
           -0.5                                    &
              *(psi(iz,ir,1,ist)*dpsiz(iz,ir,1,jst)  &
              + psi(iz,ir,2,ist)*dpsiz(iz,ir,2,jst)  &
              - psi(iz,ir,1,jst)*dpsiz(iz,ir,1,ist)  &
              - psi(iz,ir,2,jst)*dpsiz(iz,ir,2,ist)) * dynpoz(iz,ir,iq) &
          +(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)  &
              - psi(iz,ir,1,ist)*psi(iz,ir,2,jst)) * spipop(iz,ir,iq)
          
!     &                    +0.5
!     &                    *(psi(iz,ir,2,ist)*psi(iz,ir,1,jst)
!     &                    - psi(iz,ir,1,ist)*psi(iz,ir,2,jst)
!     &                    - psi(iz,ir,2,jst)*psi(iz,ir,1,ist)
!     &                    + psi(iz,ir,1,jst)*psi(iz,ir,2,ist))
!     &                    * spipop(iz,ir,iq)  &
        END DO
      END DO
      popout(ipostr) =  qopin(ipostr)*(equasi+equasj)  &
          + (ui*vj-vi*uj)*f1int(wrk1(0:nz,0:nr))
      IF(ist == jst) THEN
        scmass = 2.0 * popout(ipostr)*qopin(ipostr)   + scmass
        pqnorm = 2.0 * popin(ipostr)*betainp(ipostr)   + pqnorm
!test                  tstmas = popin(ipostr)**2/(equasi+equasj)     + tstmas
        IF(ifhcts == 1)  &
            tstmas = 2.0 * qopin(ipostr)**2*(equasi+equasj)+tstmas
        pnmass = 2.0 * qopin(ipostr)**2  &
            *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
      ELSE
        scmass = 4.0*popout(ipostr)*qopin(ipostr)     + scmass
        pqnorm = 4.0*popin(ipostr)*betainp(ipostr)     + pqnorm
!test                  tstmas = 2.0*popin(ipostr)**2/(equasi+equasj) + tstmas
        IF(ifhcts == 1)  &
            tstmas = 4.0*qopin(ipostr)**2*(equasi+equasj) + tstmas
        pnmass = 4.0 * qopin(ipostr)**2  &
            *(ui*ui-vi*vi+uj*uj-vj*vj)       + pnmass
      END IF
    END IF
  END DO
END DO
!PK      if(ifhcts.eq.1) then
WRITE(7,'(a,4(1pg12.4))') ' hcommQ: scmass,tstmas,crmasp,pqnorm=',  &
    scmass,tstmas,crmasp,pqnorm
!PK   endif
!      write(6,'(a,i6)') ' HCOMMQ: ipostr=',ipostr
IF(ifb4ad /= 1) THEN
  b4  = t4*0.5
  b4c = b4+b4p
END IF

!     all done

RETURN
END SUBROUTINE hcommq


!
! Here follow three subroutines used in the initialization phase.
!


!-----inie2q-----------------------------------------------------------

SUBROUTINE inie2q(mcomp,npcomp)

!     Initializes the two-quasi-particle energies on 'e2quas' in
!     'common /crankin/'.
!     The storage sequence depends on angular momentum 'mcomp'
!     and z-parity of the mode 'npcomp'.
!     This routine is to be called first for each mode seperately.
!     It also defines the storage space 'iposmx'.

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'


INTEGER, INTENT(IN)                      :: mcomp
INTEGER, INTENT(IN)                      :: npcomp

LOGICAL, PARAMETER :: ktest=.false.


!-----------------------------------------------------------------------


IF(mcomp == 0 .AND. .NOT.ktestm) THEN
  ipostr = 0             ! set position for lower triangular storage
  DO ist=1,nst(1)+nst(2)
    mist   = mpsi(ist)
    npist  = nparz(1,ist)
    IF(ist <= nst(1)) THEN
      jstmax = nst(1)
      eferm  = fer(1)
    ELSE
      jstmax = nst(1)+nst(2)
      eferm  = fer(2)
    END IF
    deltai = deltaf(ist)*wtpair(ist)
    edifi  = spenrg(ist)-eferm
    equasi = SQRT(edifi*edifi+deltai*deltai)
    
    DO jst=ist,jstmax
      IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
        deltaj = deltaf(jst)*wtpair(jst)
        edifj  = spenrg(jst)-eferm
        equasj = SQRT(edifj*edifj+deltaj*deltaj)
        ipostr = 1 + ipostr
        IF(ipostr > ktriac) STOP "INIE2Q: tridiagonal storage out of range"
        
        e2quas(ipostr)  = equasi+equasj
        
        IF(ktest)  WRITE(6,'(a,6i5,2f5.1,g12.4)')  &
            'INIE2Q: ist,jst,mpsi(ist),mgoal,npist,npgoal,ki,kj=',  &
            ist,jst,mpsi(ist),mpsi(ist),npist,npist,  &
            (2*mpsi(ist)+1)/2.0,(2*mpsi(ist)+1)/2.0, e2quas(ipostr)
        
      END IF
    END DO
  END DO
  WRITE(6,'(a,i6)')  ' INIE2Q: initialized, IPOSTR=',ipostr
  iposmx = ipostr
  
ELSE IF(mcomp > 0 .OR. ktestm) THEN
  rmcomp = mcomp*1.0D0
  ipostr = 0
  DO ksign=-1,1,2                ! to span all 'k' values
    iup   = (3-ksign)/2            ! effective upper component
    idw   = (3+ksign)/2            ! effective lower component
    DO ist=1,nst(1)+nst(2)
      npist  = nparz(1,ist)*ksign
      kist2  = (2*mpsi(ist)+1)*ksign     ! 2*'k'
      kgoal  = kist2+2*mcomp             ! double 'k' for final state
      IF(kgoal > 0) THEN
        mgoal  = (kgoal-1)/2             ! spatial a.m. for final state
        ksig2  = 1                       ! sign of final 'k'
        npgoal = npist*npcomp            ! parity of final state
        jup    = 1                       ! no time reversal
        jdw    = 2
      ELSE
        mgoal  = (-kgoal-1)/2
        ksig2  = -1
        npgoal = -npist*npcomp
        jup    = 2                       ! time reversal
        jdw    = 1
      END IF
      IF(ist <= nst(1)) THEN ! protons
        iq     = 1
        eferm  = fer(1)
        jstmin = 1                       ! lower bound of loop for final
        jstmax = nst(1)                  ! upper bound of loop for final
      ELSE                   ! neutrons
        iq     = 2
        eferm  = fer(2)
        jstmin = nst(1)+1
        jstmax = nst(1)+nst(2)
      END IF
      deltai = deltaf(ist)*wtpair(ist)
      edifi  = spenrg(ist)-eferm
      equasi = SQRT(edifi*edifi+deltai*deltai)
      DO jst=jstmin,jstmax
        IF(mpsi(jst) == mgoal .AND. nparz(1,jst) == npgoal) THEN
          ipostr = 1 + ipostr
          
          IF(ipostr > ktriac) STOP "INE2Q: tridiagonal storage out of range"
          
          deltaj = deltaf(jst)*wtpair(jst)
          edifj  = spenrg(jst)-eferm
          equasj = SQRT(edifj*edifj+deltaj*deltaj)
          e2quas(ipostr)  = equasi+equasj
          
          IF(ktest)  WRITE(6,'(a,6i5,2f5.1,g12.4)')  &
              'INIE2Q: ist,jst,mpsi(ist),mgoal,npist,npgoal,ki,kj=',  &
              ist,jst,mpsi(ist),mgoal,npist,npgoal,  &
              kist2/2.0,kgoal/2.0,e2quas(ipostr)
          
        END IF
      END DO
    END DO
  END DO
  WRITE(6,'(a,i6)')  ' INIE2Q: initialized, IPOSTR=',ipostr
  iposmx = ipostr
  
ELSE
END IF

RETURN
END SUBROUTINE inie2q

!-----quadms-----------------------------------------------------------

SUBROUTINE quadms(qopout,mcomp,npcomp)

!     To compute the simple mass using qadrupole momenta and equasi**-3
!     results is 'crmasq' stored on common .
!     Q operator is returned on 'qopout'.
!     The z-component of angular momentum is 'mcomp' whereby only
!     'mcomp=0' or 'mcomp=2' make sense.
!     The z-parity of the mode 'npcomp'; it must be +1.

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'



REAL(8), INTENT(OUT)                        :: qopout(ktriac)
INTEGER, INTENT(IN)                      :: mcomp
INTEGER, INTENT(IN)                      :: npcomp


!      data itstpr/0/

!     buffered SQRT as internal statement function

!      sqrtb(x) = sqrt(max(x,smal))

!-----------------------------------------------------------------------

IF(npcomp /= 1) STOP 'QUADMS: parity NPCOMP must not be negative'

! PK  qufac  = sqrt(0.2*pi)/(rms*rms*atnum)
qufac  = SQRT(0.2D0*pi)/(0.6D0*1.2D0**2*atnum**(5.0D0/3.0D0))
WRITE(6,'(a,2f8.3)') ' QUADMS: rms,A=',rms,atnum

!     the component 'mcomp=0'

IF(mcomp == 0 .AND. .NOT.ktestm) THEN
  sumq   = zero
  sumq2  = zero
  sumcrs = zero
  sumnor = zero
!  --- test PK ---
  sumbeta= zero
!  --- test PK ---
  ipostr = 0             ! set position for lower triangular storage
  DO ist=1,nst(1)+nst(2)
    mist   = mpsi(ist)
    npist  = nparz(1,ist)
    IF(ist <= nst(1)) THEN
      jstmax = nst(1)
    ELSE
      jstmax = nst(1)+nst(2)
    END IF
    vi    = sqrtb(vsqf(ist))
    ui    = sqrtb(one-vi*vi)
    DO is=1,nsd
      DO ir=0,nr
        r     = dr*(ir+sr)
        r2    = r*r
        DO iz=0,nz
          z     = ((iz)+sz-icm)*dz
          z2    = z*z
          psip(iz,ir,is) = (z2+z2-r2)*psi(iz,ir,is,ist)
        END DO
      END DO
    END DO
    
    DO jst=ist,jstmax
      IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
        vj    = sqrtb(vsqf(jst))
        uj    = sqrtb(one-vj*vj)
        qij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip(0:nz,0:nr,1:2))*qufac
        ipostr = 1 + ipostr
        IF(ipostr > ktriac) STOP "QUADMS: tridiagonal storage out of range"
! --- test PK ---
        IF (ist == jst) THEN
          sumbeta = 2.0D0*vsqf(ist)*qij + sumbeta
        END IF
! --- test PK ---
        
        qij    = qij*(vj*ui+uj*vi)
        qopout(ipostr)  = qij
        IF(ist == jst) THEN
          qij    = qij*qij
          sumq   = qij*vi*vi + sumq
        ELSE
          qij    = qij*qij*2.0
        END IF
        sumcrs  = qij/e2quas(ipostr)**3 + sumcrs
        sumnor  = qij/e2quas(ipostr) + sumnor
        sumq2   = qij + sumq2
!                sumcrs  = qij/(equasi+equasj)**3 + sumcrs
!                sumnor  = qij/(equasi+equasj) + sumnor
      END IF
    END DO
  END DO
  
!     the component 'mcomp=2'
  
ELSE IF(mcomp == 2 .OR. ktestm) THEN
  rmcomp = mcomp*1.0D0
  ipostr = 0
  DO ksign=-1,1,2                ! to span all 'k' values
    iup   = (3-ksign)/2            ! effective upper component
    idw   = (3+ksign)/2            ! effective lower component
    DO ist=1,nst(1)+nst(2)
      npist  = nparz(1,ist)*ksign
      kist2  = (2*mpsi(ist)+1)*ksign     ! 2*'k'
      kgoal  = kist2+2*mcomp
      IF(kgoal > 0) THEN
        mgoal  = (kgoal-1)/2
        ksig2  = 1
        npgoal = npist*npcomp
        jup    = 1
        jdw    = 2
      ELSE
        mgoal  = (-kgoal-1)/2
        ksig2  = -1
        npgoal = -npist*npcomp
        jup    = 2
        jdw    = 1
      END IF
      IF(ist <= nst(1)) THEN ! protons
        iq     = 1
        eferm  = fer(1)
        jstmin = 1
        jstmax = nst(1)
      ELSE                   ! neutrons
        iq     = 2
        eferm  = fer(2)
        jstmin = nst(1)+1
        jstmax = nst(1)+nst(2)
      END IF
      vi    = sqrtb(vsqf(ist))
      ui    = sqrtb(one-vi*vi)
      deltai = deltaf(ist)*wtpair(ist)
      edifi  = spenrg(ist)-eferm
      equasi = SQRT(edifi*edifi+deltai*deltai)
      DO is=1,nsd
        DO ir=0,nr
          r     = dr*(ir+sr)
          r2    = r*r
          DO iz=0,nz
            z     = ((iz)+sz-icm)*dz
            z2    = z*z
            IF(ktestm .AND. mcomp == 0) THEN
              psip(iz,ir,is) = (z2+z2-r2)*psi(iz,ir,is,ist)
            ELSE
              psip(iz,ir,is) = (r2)*psi(iz,ir,is,ist)
            END IF
          END DO
        END DO
      END DO
      DO jst=jstmin,jstmax
        IF(mpsi(jst) == mgoal .AND. nparz(1,jst) == npgoal) THEN
          ipostr = 1 + ipostr
          vj    = sqrtb(vsqf(jst))
          uj    = sqrtb(one-vj*vj)
          qij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip(0:nz,0:nr,1:2))*qufac
          qij    = qij*(vj*ui+uj*vi)
          qopout(ipostr)  = qij
          qij    = qij*qij
          sumcrs  = qij/e2quas(ipostr)**3 + sumcrs
          sumnor  = qij/e2quas(ipostr) + sumnor
        END IF
      END DO
    END DO
  END DO
  
!     everything else
  
ELSE
  STOP 'QUADMS: m-component other than 0 or 2 not possible'
END IF

crmasq = sumcrs/(sumnor*sumnor)
betafull = SQRT(sumbeta**2+sumq2)    ! total squared quadrupole
betafullp = SQRT(beta2p**2+sumq2)
IF(ipostr /= iposmx) STOP 'in QUADMS: inconsistent nr. of matrix elements'
WRITE(7,'(a,6(1pg12.4),i6)')  &
    'QUADMS: sumcrs,sumnor,crmasq,Q,Q2,beta2,ipostr=',  &
    sumcrs,sumnor,crmasq,sumbeta,sumq2,beta2,ipostr

RETURN
END SUBROUTINE quadms
!-----cringl-----------------------------------------------------------

SUBROUTINE cringl(popout)

!     Prepares collective momenum operator on 'popout' and
!     computes Inglis cranking mass for collectiev (quadrupole) motion
!     by building derivative from previously saved wavefunction.
!     Actual wavefunctions are saved on scratch for subsequent call.

!     This routine contains some test switches e.g. for testing
!     'nabla' operator (use only without z-symmetry).
!     The switch 'KTESTM' invokes rectangular storage instead
!     of triangular one; this is provided for testing other
!     routines as e.g. QCURRE etc.

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'


REAL(8), INTENT(OUT)                        :: popout(ktriac)


INTEGER,PARAMETER :: itstcm=0,itstpr=0
INTEGER,SAVE :: ifirst=0
!DATA ifirst/0/
!DATA itstcm/0/            ! switch to test case c.m. motion
!DATA itstpr/0/            ! switch to testprints of Inglis Cranking
!test      data isavsp/1/            ! save s.p. properties on file for interpol.

!old      data e2quas/ktriac*1.0/
!DATA delb/0.0/
REAL(8),SAVE :: beto2,beto2p,beto2n,beto3,ferolp,feroln,delb

!     buffered SQRT as internal statement function
!sqrtb(x) = SQRT(MAX(x,smal))

!-----------------------------------------------------------------------

!     double loop over states with twofold intention:
!      1. to compute and store the matrix of collective derivative 'popout'
!      2. to accumulate Inglis cranking mass from derivatives

IF(itstcm == 1) THEN
  ifirst = 1
  betao2=0.0
  betao3=0.0
  IF(sym > zero) STOP ' test c.m. requires SYM=0'
END IF

IF(ifirst /= 0) THEN
  IF(lqcon == 2) THEN
    delb   = beta2-beto2
  ELSE IF(lqcon == 3) THEN
    delb   = beta3-beto3
  END IF
  IF(itstpr /= 0) WRITE(6,'(a,3(1pg12.4))')  &
      ' BETAs=',delb,beta2,beto2
  sumcrp = zero
  sumdig = zero
  sumwdp = zero
  sumwdd = zero
  ipostr = 0             ! set position for lower triangular storage
  isigmin = 1
  IF(ktestm) isigmin=-1
  DO isigm=isigmin,1,2
    DO ist=1,nst(1)+nst(2)
      mist   = mpsi(ist)
      npist  = nparz(1,ist)
      IF(ist <= nst(1)) THEN
        IF(ktestm) jstmin = 1
        jstmax = nst(1)
!               eferm  = fer(1)
!R               efermo = ferolp
      ELSE
        IF(ktestm) jstmin = nst(1)+1
        jstmax = nst(1)+nst(2)
!               eferm  = fer(2)
!R               efermo = feroln
      END IF
      IF(.NOT.ktestm) jstmin = ist
      IF(itstcm /= 1) THEN
        READ(98) (((psip(iz,ir,is),iz=0,nz),ir=0,nr),is=1,nsd),  &
            speol,vsqol,delol,wtpol
        IF(itstpr > 0) WRITE(7,'(a,i5,4(1pg12.4))')  &
            ' ist,speol,vsqol,delol,wtpol=',ist,speol,vsqol, delol,wtpol
      END IF
!                                         test overlaps
      IF(itstpr /= 0) THEN
        DO jst=1,jstmax
          IF(mist == mpsi(jst) .AND. npist == nparz(1,jst))  &
              WRITE(6,'(a,2i4,2(1pg12.4))') ' test overlaps: ist,jst,ov=',  &
              ist,jst,wfovlp(psip(0:nz,0:nr,1:2),psi(0:nz,0:nr,1:2,jst)),  &
              wfovlp(psi(0:nz,0:nr,1:2,ist),psi(0:nz,0:nr,1:2,jst))
        END DO
      END IF
!                                         derivative
      IF(itstcm /= 1) THEN
        DO is=1,nsd
          DO ir=0,nr
            DO iz=0,nz
              psip(iz,ir,is) = (psi(iz,ir,is,ist)-psip(iz,ir,is))  &
                  /delb
            END DO
          END DO
        END DO
        
!     gradient operator Nabla_z on 'psip'
        
      ELSE
        DO is=1,nsd
          mm        = mist+is-1
          npis     = nparz(is,ist)
          CALL mulz(derz1(0:nz,0:nz,npis),psi(0:nz,0:nr,is,ist),wrk4(0:nz,0:nr),one,0)
          DO ir=0,nr
            DO iz=0,nz
              psip(iz,ir,is) = wrk4(iz,ir)
            END DO
          END DO
        END DO
      END IF
!                                         auxiliary variables
      vi    = sqrtb(vsqf(ist))
      ui    = sqrtb(one-vi*vi)
!            deltai = deltaf(ist)*wtpair(ist)
!            edifi  = spenrg(ist)-eferm
!            equasi = sqrt(edifi*edifi+deltai*deltai)
!                                         diagonal part of cranking mass
      vold   = sqrtb(vsqol)
      uold   = sqrtb(one-vsqol)
      vnew   = vi
      unew   = ui
      uaver  = 0.5*(unew+uold)
      vaver  = 0.5*(vnew+vold)
      IF(vaver > 0.5*one) THEN
        pii    = -(unew-uold)/(vaver*delb)
      ELSE
        pii    = (vnew-vold)/(uaver*delb)
      END IF
!old            ipostr = 1 + ipostr
      IF(ipostr > ktriac) STOP "CRINGL: tridiagonal storage out of range"
!old            if(kcrank.gt.0) then
!old              popout(ipostr)  = pii
!cold              e2quas(ipostr) = equasi+equasi
!old            endif
!            sumdig = pii*pii/(equasi+equasi) + sumdig
      sumwdd = pii*pii + sumwdd
!test            write(7,'(a,5(1pg12.4))')
!test     &        '     uold,vold,unew,vnew,delb=',uold,vold,unew,vnew,delb
!                                         off-diagonal elements
      
      
      DO jst=jstmin,jstmax
        IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
          ipostr = 1 + ipostr
          vj    = sqrtb(vsqf(jst))
          uj    = sqrtb(one-vj*vj)
!                  deltaj = deltaf(jst)*wtpair(jst)
!                  edifj  = spenrg(jst)-eferm
!                  equasj = sqrt(edifj*edifj+deltaj*deltaj)
          pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip(0:nz,0:nr,1:2))
          IF(itstpr /= 0) WRITE(6,'(a,2i4,2f8.2,4(1pg12.4))')  &
              ' i,j,equasi/j,phase,pij=', ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
          pij    = pij*(vj*ui-uj*vi)
          IF(ist == jst) THEN
            IF(kcrank > 0) popout(ipostr)  = pii
            sumdig = pii*pii/e2quas(ipostr) + sumdig
            pij     = pij*pij
            IF(itstpr == 1) WRITE(7,'(a,2i5,2(1pg12.4))')  &
                ' ist,ipostr,pii,equasi=',ist,ipostr,pii,equasi
          ELSE
            IF(kcrank > 0) popout(ipostr)  = pij
            IF(ktestm) THEN
              pij     = pij*pij
            ELSE
              pij     = 2.0*pij*pij
            END IF
          END IF
          sumcrp = pij/e2quas(ipostr) + sumcrp
          sumwdp = pij + sumwdp
        END IF
      END DO
    END DO
    REWIND 98
  END DO
  IF(ktestm) THEN
    crmaso = sumcrp
    crmasp = (sumcrp+sumdig)
    pwido  = sumwdp*0.5D0
    pwid   = (sumwdp+sumwdd)*0.5D0
  ELSE
    crmaso = sumcrp*2.0
    crmasp = (sumcrp+sumdig)*2.0
    pwido  = sumwdp
    pwid   = sumwdp+sumwdd
  END IF
  
  WRITE(6,'(a,4g12.4,i5)') ' CRINGL: crmaso,crmasp,pwido,pwid,ipostr=',  &
      crmaso,crmasp,pwido,pwid,ipostr
  IF(ipostr /= iposmx) STOP 'CRINGL: inconsistent nr. of matrix elements '
  
END IF
IF(ktestm)  WRITE(7,'(a,i10)') ' CRINGL: IPOSTR=',ipostr

!     save actual wavefunction and associated s.p. properties

IF(ifirst == 0) THEN
  OPEN(98,FORM='unformatted',STATUS='scratch')
  ifirst = 1
ELSE
  REWIND 98
END IF
beto2  = beta2
beto3  = beta3
beto2p = beta2p
beto2n = beta2n
ferolp = fer(1)
feroln = fer(2)
DO ist=1,nst(1)+nst(2)
  WRITE(98) (((psi(iz,ir,is,ist),iz=0,nz),ir=0,nr),is=1,nsd),  &
      spenrg(ist),vsqf(ist),deltaf(ist),wtpair(ist)
  WRITE(6,'(a,i4)') ' scratch file written. record nr.',ist
END DO
WRITE(98)   ifirst        ! dummy record
REWIND 98

!     optionally save s.p. properties for later interpolation

!test      if(isavsp.eq.1) then
!test        if(ifirst.eq.0) then
!test           open(97,form='unformatted',file='testint')
!test           write(97) nst
!test        endif
!test        write(97)
!test     &    beta2,crmaso,pwido,fer,
!test     &    (spenrg(ist)-fer(1),vsqf(ist),deltaf(ist)*wtpair(ist),
!test     &     ist=1,nst(1)),
!test     &    (spenrg(ist)-fer(2),vsqf(ist),deltaf(ist)*wtpair(ist),
!test     &     ist=nst(1)+1,nst(1)+nst(2))
!test      endif


RETURN
END SUBROUTINE cringl

!REAL(8) PURE FUNCTION sqrtb(x)
!USE paramr, ONLY: smal
!REAL(8),INTENT(IN) :: x
!sqrtb = SQRT(MAX(x,smal))
!RETURN
!END


!
!  Here follows a block of subroutines for testing commutators etc
!

!-----hcPT--------------------------------------------------------------

SUBROUTINE hcpt(betain2)

!     Test-routine for Commutator [H,P]
!     using P-operator 'presp' from routine 'crankin'

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)

REAL(8),INTENT(IN) :: betain2(ktriac)   

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------

WRITE(7,'(a,1pg12.4)') ' now test HCOMMP: initial <P**2>=',opovlp(presp,presp)

CALL hcommp(presp,qresp)

!     test

WRITE(7,'(a,3(1pg12.4))') ' hcommP: <P**2>,<[P,[H,P]]>,<Q**2>=',  &
    opovlp(presp,presp),opovlp(presp,qresp),opovlp(qresp,qresp)

CALL hcommq(qresp,presp,betain2,ptest,0)
WRITE(7,'(a,3(1pg12.4))') ' hcommQ: <[Q,[H,Q]]>,<P_new**2>=',  &
    opovlp(ptest,qresp),opovlp(ptest,ptest)

RETURN
END SUBROUTINE hcpt

!-----hcQT--------------------------------------------------------------

SUBROUTINE hcqt(icase)

!     Test-routine for Commutator [H,Q] of Hamiltonian
!     with cases {b1=0,b1 not 0}x{Q1=2z^2-r^2, Q2=2z^2, Q3=-r^2}

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------

!     0. step:


IF(itest == 0) THEN       ! only for testing
  
!     test particle number
  
  DO ir=0,nr
    DO iz=0,nz
!               wrk1(iz,ir)  = 0.0
!               wrk2(iz,ir)  = 0.0
      wrk4(iz,ir) = rho(iz,ir,1) + rho(iz,ir,2)
    END DO
  END DO
  WRITE(17,'(/a,1(1pg12.4))')  &
      '> hcQT: test output : Tr{1}  = #particles =', f1int(wrk4(0:nz,0:nr))
  
  DO i=1,2               ! i=1 -> b1=0, i=2 -> b1 not 0
    
!     setting force parametres
    
    IF (i == 1)  THEN
      b1  = 0.0
      b1p = 0.0
    ELSE
      b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
    END IF
    b1c = b1-b1p
    
    WRITE(17,'(a,2(1pg12.4))') '> hcQT: test case: b1, b1p  = ',b1,b1p
    
    
!     calculating test masses icases: 3 , 1 , 2
    
!     mass1 \       / Q1 = 2 z^2 - r^2 => (Nabla Q1)^2 =  16 z^2 +  4 r^2
!     mass2  | <-> |  Q2 = 2 z^2       => (Nabla Q2)^2 =  16 z^2
!     mass3 /       \ Q3 = - r^2       => (Nabla Q3)^2 =   4 r^2
    
!     making integrant in workspaces
    
    DO ir=0,nr
      DO iz=0,nz
        
!     coordinates
        
        r  = dr*(ir+sr)
        r2 = r*r
        z  = dz*(iz+sz-icm)
        z2 = z*z
        
!    h2meffq ??? factor 2.0 : +2 b1 ???
        
        h2meffp =(h2m(1)+b1c*rho(iz,ir,1)+b1*rho(iz,ir,2))
        h2meffn =(h2m(2)+b1c*rho(iz,ir,2)+b1*rho(iz,ir,1))
        
        wh2meff = h2meffp+h2meffn
        
        wrk1(iz,ir) = wh2meff * ( z2 ) * rho(iz,ir,1)
        wrk2(iz,ir) = wh2meff * ( r2 ) * rho(iz,ir,2)
      END DO
    END DO
    
    
!     final integrating, with factor
    
    dummy1 = f1int(wrk1(0:nz,0:nr))
    dummy2 = f1int(wrk2(0:nz,0:nr))
    xmass1 =  8.0  * (4.0*dummy1+dummy2)
    xmass2 = 32.0  * (dummy1)
    xmass3 =  8.0  * (dummy2)
    
!R250599
    
!            write(7,'(a/,a,3(1pg12.4))')
!     &           '> hcQT: 2z^2-r^2,  2z^2,   -r^2',
!     &           '> hcQT: xmass1,  xmass2, xmass3 (inv) = ',
!     &           xmass1,xmass2,xmass3
    
    IF(icase == 1) THEN
      WRITE(17,'(a/,a,1(1pg12.4))') '> hcQT: 2z^2 ',  &
          '> hcQT: xmass2 (inv) = ', xmass2
    ELSE IF(icase == 2) THEN
      WRITE(17,'(a/,a,1(1pg12.4))') '> hcQT: -r^2 ',  &
          '> hcQT: xmass3 (inv) = ', xmass3
    ELSE
      WRITE(17,'(a/,a,1(1pg12.4))') '> hcQT: 2z^2-r^2 ',  &
          '> hcQT: xmass1 (inv) = ', xmass1
    END IF
    
!     all done
    
  END DO
END IF
RETURN
END SUBROUTINE hcqt

!-----Qabtst-----------------------------------------------------------

SUBROUTINE qabtst(icase,betain2)

!     calculates the matrix-elements Q_\alpha\beta of Q=2z^2-r^2 for
!     testing and stores results on 'qresp' where it can be processed
!     further.
!     Three cases are provided:
!      icase = 1  --> Q=2z^2
!      icase = 2  --> Q=-r^2
!      icase = 3  --> Q=2z^2-r^2

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)
REAL(8),INTENT(IN) :: betain2(ktriac)   
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------


IF(icase == 1) THEN
  WRITE(7,'(/a)') ' test HCOMMQ with Q=2*z**2:'
  WRITE(17,'(/a)') '# Qabtst: test HCOMMQ with Q=2*z^2:'
ELSE IF(icase == 2) THEN
  WRITE(7,'(/a)') ' test HCOMMQ with Q=-r**2:'
  WRITE(17,'(/a)') '# Qabtst: test HCOMMQ with Q=-r^2:'
ELSE
  WRITE(7,'(/a)') ' test HCOMMQ with Q=2*z**2-r**2:'
  WRITE(17,'(/a)') '# Qabtst: test HCOMMQ with Q=2*z^2-r^2:'
END IF

sumq   = 0.0
sumq2  = 0.0
ipostr = 0

DO ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
  
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  vi    = SQRT(vsqf(ist))
  ui    = SQRT(one-vi*vi)
  
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    jstmax = nst(1)+nst(2)
  END IF
  
  DO jst=ist,jstmax      ! <-> beta
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj    = SQRT(vsqf(jst))
      uj    = SQRT(one-vj*vj)
      IF(icase == 1) THEN
        DO ir=0,nr
          DO iz=0,nz
            z  = dz*(iz+sz-icm)
            z2 = z*z
            wrk1(iz,ir) =                               &
                2.0*z2*(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                + psi(iz,ir,2,ist)*psi(iz,ir,2,jst))  ! Q=2z^2  
          END DO
        END DO
      ELSE IF(icase == 2) THEN
        DO ir=0,nr
          r  = dr*(ir+sr)
          r2 = r*r
          DO iz=0,nz
            wrk1(iz,ir) =                                &
                -r2*(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                + psi(iz,ir,2,ist)*psi(iz,ir,2,jst))   ! Q=2z^2
          END DO
        END DO
      ELSE
        DO ir=0,nr
          r  = dr*(ir+sr)
          r2 = r*r
          DO iz=0,nz
            z  = dz*(iz+sz-icm)
            wrk1(iz,ir) =  (2.0*(z*z)-r2)              &
                *(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                + psi(iz,ir,2,ist)*psi(iz,ir,2,jst))   ! quadrupole
          END DO
        END DO
      END IF
      
!     integrating and storing
      
      qmatre = f1int(wrk1(0:nz,0:nr))
      qresp(ipostr) = qmatre*(ui*vj+uj*vi)
!R---------------------------------------------------------------------
!r fehlt faktor 2?
      IF(ist == jst) THEN
        sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
        sumq  = qmatre*vi*vi + sumq
      ELSE
        sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
      END IF
!R---------------------------------------------------------------------
      
    END IF
  END DO
END DO
iposmx = ipostr
WRITE(7,'(a,2f8.3,a,2(1pg12.4))')  &
    ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
WRITE(17,'(a,2f8.3,a,2(1pg12.4))')  &
    '# Qabtst: for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2

!     testing by calling commutator

CALL hcommq(qresp,presp,betain2,ptest,icase)

RETURN
END SUBROUTINE qabtst

!-----Sabtst-----------------------------------------------------------
!R040399

SUBROUTINE sabtst(icase,betain2)

!     calculates the matrix-elements Q_\alpha\beta of Q=sigma_r for
!     testing and stores results on 'qresp' where it can be processed
!     further. (more or less the same as Qabtst)
!      icase = 3  --> Q=2z^2-r^2

!     icase is only used for the call of hcommQ

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)
REAL(8),INTENT(IN) :: betain2(ktriac)   
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------


WRITE(7,'(/a)') ' test HCOMMQ with Q=sigma_r:'
WRITE(17,'(/a)') ' test HCOMMQ with Q=sigma_r:'

sumq   = 0.0
sumq2  = 0.0
ipostr = 0

DO ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
  
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  vi     = SQRT(vsqf(ist))
  ui     = SQRT(one-vi*vi)
  
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    jstmax = nst(1)+nst(2)
  END IF
  
  DO jst=ist,jstmax      ! <-> beta
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj    = SQRT(vsqf(jst))
      uj    = SQRT(one-vj*vj)
      
      DO ir=0,nr
        DO iz=0,nz
!                     z  = dz*(iz+sz-icm)
!                     z2 = z*z
          wrk1(iz,ir) = ( psi(iz,ir,1,ist)*psi(iz,ir,2,jst)  &
              + psi(iz,ir,2,ist)*psi(iz,ir,1,jst))   ! Q=sigma_r
        END DO
      END DO
      
!     integrating and storing
      
      qmatre = f1int(wrk1(0:nz,0:nr))
      qresp(ipostr) = qmatre*(ui*vj+uj*vi)
!r--------------------------------------------------------
!R ??? fehlt da nicht ein Faktor 2 ?
      IF(ist == jst) THEN
        sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
        sumq  = qmatre*vi*vi + sumq
      ELSE
        sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
      END IF
!R
!r--------------------------------------------------------
      
    END IF
  END DO
END DO
iposmx = ipostr
WRITE(7,'(a,2f8.3,a,2(1pg12.4))')  &
    ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
WRITE(17,'(a,2f8.3,a,2(1pg12.4))')  &
    '# Sabtst: for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2

!     testing by calling commutator

CALL hcommq(qresp,presp,betain2,ptest,icase)

RETURN
END SUBROUTINE sabtst

!-----------------------------------------------------------------------
!R210499

SUBROUTINE hcdt

!     Test-routine for Commutator [H,Q] of Hamiltonian
!     with the case Q=z

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------

!     0. step:


IF(itest == 0) THEN       ! only for testing
  
!     test particle number
  
  DO ir=0,nr
    DO iz=0,nz
!               wrk1(iz,ir)  = 0.0
!               wrk2(iz,ir)  = 0.0
      wrk4(iz,ir) = rho(iz,ir,1) + rho(iz,ir,2)
    END DO
  END DO
  WRITE(17,'(/a,1(1pg12.4))')  &
      '> hcDT: test output : Tr{1}  = #particles =', f1int(wrk4(0:nz,0:nr))
  
  DO i=1,2               ! i=1 -> b1=0, i=2 -> b1 not 0
    
!     setting force parametres
    
    IF (i == 1)  THEN
      b1  = 0.0
      b1p = 0.0
      b1c = 0.0
    ELSE
      b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      b1c = b1-b1p
    END IF
    
    WRITE(17,'(a,2(1pg12.4))') '> hcDT: test case: b1, b1p  = ',b1,b1p
    
    
!     calculating test masse
    
!     mass2   <->   Q =  z        => (Nabla Q)^2 =  1
    
!     making integrant in workspaces
    
    DO ir=0,nr
      DO iz=0,nz
        
!     coordinates
        
        r  = dr*(ir+sr)
        z  = dz*(iz+sz-icm)
        
!    h2meffq +2.0 in 2.zeile ????
        
        
        h2meffp = h2m(1)+b1c*rho(iz,ir,1)+b1*rho(iz,ir,2)
        h2meffn = h2m(2)+b1c*rho(iz,ir,2)+b1*rho(iz,ir,1)
        
        wrk1(iz,ir) = h2meffp*rho(iz,ir,1) + h2meffn*rho(iz,ir,2)
        
      END DO
    END DO
    
    
!     final integrating, with factor
    
    xmass = f1int(wrk1(0:nz,0:nr))
    
    WRITE(17,'(a,1(1pg12.4))') '> hcDT: xmass (inv) = ',  &
        xmass
    
!     all done
    
  END DO
END IF
RETURN
END SUBROUTINE hcdt

!-----Dabtst-----------------------------------------------------------
!R200499

SUBROUTINE dabtst(icase,betain2)

!     calculates the matrix-elements Q_\alpha\beta of Q=z for
!     testing and stores results on 'qresp' where it can be processed
!     further.

!     icase is only used for the call of hcommQ

!     Dipole-Test

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)
REAL(8),INTENT(IN) :: betain2(ktriac)   
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------


WRITE(7,'(/a)') ' test HCOMMQ with Q=z:'
WRITE(17,'(/a)') ' test HCOMMQ with Q=z:'

sumq   = 0.0
sumq2  = 0.0
ipostr = 0

DO ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
  
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  vi     = SQRT(vsqf(ist))
  ui     = SQRT(one-vi*vi)
  
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    jstmax = nst(1)+nst(2)
  END IF
  
  DO jst=ist,jstmax      ! <-> beta
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj    = SQRT(vsqf(jst))
      uj    = SQRT(one-vj*vj)
      
      DO ir=0,nr
        DO iz=0,nz
          z  = dz*(iz+sz-icm)
          
          wrk1(iz,ir) = ( psi(iz,ir,1,ist)*z*psi(iz,ir,1,jst)  &
              + psi(iz,ir,2,ist)*z*psi(iz,ir,2,jst))  ! Q=z
        END DO
      END DO
      
!     integrating and storing
      
      qmatre = f1int(wrk1(0:nz,0:nr))
      qresp(ipostr) = qmatre*(ui*vj+uj*vi)
!r----------------------------------------------------------
!R fehlt da nicht ein Faktor 2?
      IF(ist == jst) THEN
        sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
        sumq  = qmatre*vi*vi + sumq
      ELSE
        sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
      END IF
!R----------------------------------------------------------
      
    END IF
  END DO
END DO
iposmx = ipostr
WRITE(7,'(a,2f8.3,a,2(1pg12.4))')  &
    ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
WRITE(17,'(a,2f8.3,a,2(1pg12.4))')  &
    '# Dabtst for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2

!     testing by calling commutator

CALL hcommq(qresp,presp,betain2,ptest,icase)

RETURN
END SUBROUTINE dabtst
!#######################################################################
!R250599
!-----PSabtst-----------------------------------------------------------

SUBROUTINE psabtst(icase,betain2)

!     calculates the matrix-elements Q_\alpha\beta of Q= p_r sigma_r for
!     testing and stores results on 'qresp' where it can be processed
!     further.

!     icase is only used for the call of hcommQ

!     p_r*sigma_r Test

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)
REAL(8),INTENT(IN) :: betain2(ktriac)   
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------


WRITE(7,'(/a)') ' test HCOMMQ with Q= p_r sigma_r:'
WRITE(17,'(/a)') ' test HCOMMQ with Q= p_r sigma_r:'

sumq   = 0.0
sumq2  = 0.0
ipostr = 0

DO ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
  
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  vi     = SQRT(vsqf(ist))
  ui     = SQRT(one-vi*vi)
  
  IF(ist <= nst(1)) THEN ! protons
    iq     = 1
    jstmax = nst(1)
  ELSE                   ! neutrons
    iq     = 2
    jstmax = nst(1)+nst(2)
  END IF
  
  DO jst=ist,jstmax      ! <-> beta
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      vj    = SQRT(vsqf(jst))
      uj    = SQRT(one-vj*vj)
      
      DO ir=0,nr
        DO iz=0,nz
          z  = dz*(iz+sz-icm)
          
! Q= p_r sigma_r! p_r = -i d_r  &
          wrk1(iz,ir) = &
          -0.5*( dpsir(iz,ir,1,ist) * psi(iz,ir,2,jst)  &
              + psi(iz,ir,1,ist)   * dpsir(iz,ir,2,jst)  &
              + dpsir(iz,ir,2,ist) * psi(iz,ir,1,jst)  &
              + psi(iz,ir,2,ist)   * dpsir(iz,ir,1,jst) )
        END DO
      END DO
      
!     integrating and storing
      
      qmatre = f1int(wrk1(0:nz,0:nr))
      qresp(ipostr) = qmatre*(ui*vj+uj*vi)
!r----------------------------------------------------------
!R fehlt da nicht ein Faktor 2?
      IF(ist == jst) THEN
        sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
        sumq  = qmatre*vi*vi + sumq
      ELSE
        sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
      END IF
!R----------------------------------------------------------
      
    END IF
  END DO
END DO
iposmx = ipostr
WRITE(7,'(a,2f8.3,a,2(1pg12.4))')  &
    ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
WRITE(17,'(a,2f8.3,a,2(1pg12.4))')  &
    '# PSabtst for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2

!     testing by calling commutator

CALL hcommq(qresp,presp,betain2,ptest,icase)

RETURN
END SUBROUTINE psabtst
!#######################################################################
!R250599
!-----Oabtst-----------------------------------------------------------

SUBROUTINE oabtst(icase,betain2)

!     Major-Routine to calculate Test-Operator-Matrixelements

!     with cases

!     Qabtst:   icase = 1,2,3

!                     1: Q = 2z^2
!                     2: Q = -r^2
!                     3: Q = 2z^2-r^2

!     Dabtst:   icase = 4: Q = z

!     Sabtst:   icase = 5: Q = sigma_r
!     PSabtest: icase = 6: Q = p_r sigma_r

!USE paramr
!USE cranking, ONLY: ktriac
IMPLICIT REAL(8) (a-h,o-z)
REAL(8),INTENT(IN) :: betain2(ktriac)   

IF((icase >= 1).AND.(icase <= 3)) THEN
  CALL qabtst(icase,betain2)
ELSE IF(icase == 4) THEN
  CALL dabtst(icase,betain2)
ELSE IF(icase == 5) THEN
  CALL sabtst(icase,betain2)
ELSE IF(icase == 6) THEN
  CALL psabtst(icase,betain2)
ELSE IF(icase > 6) THEN
  STOP 'Oabtst: unkown icase'
END IF

RETURN
END SUBROUTINE oabtst
!#######################################################################
!R250599
!-----hcOT-----------------------------------------------------------

SUBROUTINE hcot(icase)

!     Major-Routine to call hcommQ-Test-Routines

!     with cases

!     hcQT: icase = 1,2,3

!                     1: Q = 2z^2
!                     2: Q = -r^2
!                     3: Q = 2z^2-r^2

!     hcDT: icase = 4  : Q = z
!     icase 5: sigma_r
!     icase 6: pr_sigma_r

!USE paramr
!USE cranking, ONLY: ktriac
IMPLICIT REAL(8) (a-h,o-z)
!REAL(8),INTENT(IN) :: betain2(ktriac)   

IF((icase >= 1).AND.(icase <= 3)) THEN
!old  CALL hcqt(icase,betain2)
  CALL hcqt(icase)
ELSE IF(icase == 4) THEN
!old  CALL hcdt(betain2)
  CALL hcdt()
ELSE IF(icase > 6) THEN
  STOP 'hcOT: unkown icase'
END IF

RETURN
END SUBROUTINE hcot
!#######################################################################
!R250599
!-----ptstcur-----------------------------------------------------------

SUBROUTINE ptstcur(icase)

!     subroutine of hcommQ: output of test currents on unit 17

!     with cases

!     ptstcur: icase = 1,2,3

!                     1: Q = 2z^2
!                     2: Q = -r^2
!                     3: Q = 2z^2-r^2

!     ptstcur: icase = 4  Q = z
!     ptstcur: icase = 5  Q = sigma_r
!     ptstcur: icase = 6  Q = p_r sigma_r ! todo


!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!----------------------------------------------------------------------

IF(icase <= 3) THEN       ! Qadrupol
  
  IF(icase == 1) THEN    ! plots of currents, sigma_Q,\phi==0
    WRITE(17,'(//a)') '# ptstcur: plot of currents for Q=2z^2'
  ELSE IF(icase == 2) THEN
    WRITE(17,'(//a)') '# ptstcur: plot of currents for Q=-r^2'
  ELSE
    WRITE(17,'(//a)') '# ptstcur: plot of currents for Q=2z^2-r^2'
  END IF
  
  iz0 = (nz+1)/2+1       ! this is z=0+dz
  
  WRITE(17,'(/a)') '# test curqz for ir=1 :',  &
      '# iz  z  curqz_p  curqz/rho_p  curqz_n  curqz/rho_n '
  DO iz=0,nz
    WRITE(17,'(1x,i3,f8.2,4(1pg12.4))') iz,(iz+sz-icm)*dz,curqz(iz,1,1),  &
        curqz(iz,1,1)/rho(iz,1,1), curqz(iz,1,2),curqz(iz,1,2)/rho(iz,1,2)
  END DO
!rold         write(17,'(/a)')
!rold     &        '# test curqr for iz=0 :',
!rold     &        '# ir  r  curqr  curqr/rho '
!rold         do ir=0,nr
!rold            write(17,'(1x,i3,f8.2,4(1pg12.4))')
!rold     &           ir,(ir+sr)*dr,curqr(0,ir,1),curqr(0,ir,1)/rho(0,ir,1),
!rold     &           curqr(0,ir,2),curqr(0,ir,2)/rho(0,ir,2)
!rold         enddo
  WRITE(17,'(/a)') '# test curqr for z=0+dz :',  &
      '# ir  r  curqr_p  curqr/rho_p  curqr_n  curqr/rho_n '
  DO ir=0,nr
    WRITE(17,'(1x,i3,f8.2,4(1pg12.4))') ir,(ir+sr)*dr,curqr(iz0,ir,1),  &
        curqr(iz0,ir,1)/rho(iz0,ir,1), curqr(0,ir,2),curqr(iz0,ir,2)/rho(iz0,ir,2)
  END DO
  
!R030399 test currents j_Q,q Q=2z^2-r^2
!R210499
!rold         write(17,'(/a/)') '# Q=2z^2-r^2 '
!rold         write(17,'(/a,2g12.6)') '# dr,dz= ',
!rold     &        dr,dz
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,p;r /(rho_p*dr))|iz=0,ir=1 = -2r * ',
!rold     &        curqr(0,iact,1) / (-2.0 * rho(0,iact,1) * dr)
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,p;z /(rho_p*dz))|ir=0,iz=1 =  4z * ',
!rold     &        curqz(iact,0,1) / ( 4.0 * rho(iact,0,1) * dz)
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,n;r /(rho_n*dr))|iz=0,ir=1 = -2r * ',
!rold     &        curqr(0,iact,2) / (-2.0 * rho(0,iact,2) * dr)
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,n;z /(rho_n*dz))|ir=0,iz=1 =  4z * ',
!rold     &        curqz(iact,0,2) / ( 4.0 * rho(iact,0,2) * dz)
!rold
  IF(icase == 1) THEN    ! plots of currents, sigma_Q,\phi==0
    WRITE(17,'(//a)') '# ptstcur: Q=2z^2'
  ELSE IF(icase == 2) THEN
    WRITE(17,'(//a)') '# ptstcur: Q=-r^2'
  ELSE
    WRITE(17,'(//a)') '# ptstcur: Q=2z^2-r^2'
  END IF
  
  WRITE(17,'(/a,2g12.6)') '# dr,dz= ', dr,dz
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 = -2r * ',  &
      curqr(iz0,1,1) / (-2.0 * rho(iz0,1,1) * dr)
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  4z * ',  &
      curqz(iz0,1,1) / ( 4.0 * rho(iz0,1,1) * dz)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 = -2r * ',  &
      curqr(iz0,1,2) / (-2.0 * rho(iz0,1,2) * dr)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  4z * ',  &
      curqz(iz0,1,2) / ( 4.0 * rho(iz0,1,2) * dz)
  
ELSE IF(icase == 4) THEN ! Q=z
  iz0 = (nz+1)/2+1 ! this is z=0+dz
!rold         write(17,'(/a)')
!rold     &        '# ptstcur: plot of currents for Q=z'
!roldc
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,p;r /(rho_p*dr))|ir=1,iz=0 =  0 ',
!rold     &        curqr(0,1,1)
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,p;z /(rho_p*dz))|ir=0,iz=1 =  rho /2 ',
!rold     &        curqz(iact,0,1) / ( rho(iact,0,1) * dz) * 1.0
!rold         write(17,'(a,1g12.6)')
!rold     &        '#(j_Q,n;r /(rho_n*dr))|ir=1,iz=0 =  0 ',
!rold     &        curqr(0,1,2)
!rold         write(17,'(a,1g12.6)')
!rold     &        '# (j_Q,n;z /(rho_n*dz))|ir=0,iz=1 =  rho /2 ',
!rold     &        curqz(iact,0,2) / ( rho(iact,0,2) * dz) * 1.0
  WRITE(17,'(/a)') '# ptstcur: plot of currents for Q=z'
  
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 =  0 ',  &
      curqr(iz0,1,1)
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  rho /2 ',  &
      curqz(iz0,1,1) / ( rho(iz0,1,1) * dz) * 1.0
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 =  0 ',  &
      curqr(iz0,1,2)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  rho /2 ',  &
      curqz(iz0,1,2) / ( rho(iz0,1,2) * dz) * 1.0
  
ELSE IF(icase == 5) THEN
  iz0 = (nz+1)/2+1 ! this is z=0+dz
  WRITE(17,'(/a)') '# ptstcur: plot of currents for Q=sigma_r'
!R040399 test currents sigma_Q,q Q=sigma_r
!     for protons ->wrk1, neutons ->wrk2
  
!     cleaning up
  DO ir=0,nr
    DO iz=0,nz
      wrk1(iz,ir) = 0.0
      wrk2(iz,ir) = 0.0
      wrk3(iz,ir) = 0.0
      wrk4(iz,ir) = 0.0
    END DO
  END DO
  
  DO ist=1,nst(1)+nst(2)
    vi     = SQRT(vsqf(ist))
    ui     = SQRT(one-vi*vi)
    mist   = mpsi(ist)
    npist  = nparz(1,ist)
    IF(ist <= nst(1)) THEN ! protons
      iq     = 1
      eferm  = fer(1)
      jstmax = nst(1)
    ELSE                ! neutrons
      iq     = 2
      eferm  = fer(2)
      jstmax = nst(1)+nst(2)
    END IF
    DO jst=ist,jstmax
      IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
        vj     = SQRT(vsqf(jst))
        uj     = SQRT(one-vj*vj)
        IF(ist /= jst) THEN
          factor = (ui*ui*vj*vj-vi*vi*uj*uj)
          DO ir=0,nr
            DO iz=0,nz
              IF(iq == 1) THEN ! protons
                wrk1(iz,ir) = -2.0 *(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                    - psi(iz,ir,2,ist)*psi(iz,ir,2,jst)) *factor  +  wrk1(iz,ir)
              ELSE ! neutrons
                wrk1(iz,ir) = -2.0 *(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                    - psi(iz,ir,2,ist)*psi(iz,ir,2,jst)) *factor  +  wrk2(iz,ir)
              END IF
            END DO
          END DO
        END IF
      END IF
    END DO
  END DO
  
  
  WRITE(17,'(/a)') '# test spinph for ir=1 :',  &
      '# iz  z  spinph_p  exact_p spinph_n  exact_n '
  DO iz=0,nz
    WRITE(17,'(1x,i3,f8.2,4(1pg12.4))')  &
        iz,(iz+sz-icm)*dz,spinph(iz,1,1),wrk1(iz,1), spinph(iz,1,2),wrk2(iz,1)
  END DO
  WRITE(17,'(a)') '# test curqr for z=0+dz :',  &
      '# ir  r  curqr_p  curqr/rho_p  curqr_n  curqr/rho_n '
  DO ir=0,nr
    WRITE(17,'(1x,i3,f8.2,4(1pg12.4))')  &
        ir,(ir+sr)*dr,spinph(iz0,ir,1),wrk1(iz0,ir), spinph(iz0,ir,2),wrk2(iz0,ir)
  END DO
  
  DO ir=0,nr
    DO iz=0,nz
      wrk3(iz,ir)=(wrk1(iz,ir)-spinph(iz,ir,1))**2
      wrk4(iz,ir)=(wrk2(iz,ir)-spinph(iz,ir,2))**2
    END DO
  END DO
  
  WRITE(17,'(/a)') '# ptstcur: plot of currents for Q=sigma_r'
  
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 =  0 ',  &
      curqr(iz0,1,1)
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  0 ',  &
      curqz(iz0,1,1)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 =  0 ',  &
      curqr(iz0,1,2)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  0 ',  &
      curqz(iz0,1,2)
  
  xst1 = f1int(wrk3(0:nz,0:nr))
  xst2 = f1int(wrk4(0:nz,0:nr))
  
  WRITE(17,'(/a,1g12.6)') '# sigma_Q=sigma_r,p: 0?=xst1 ', xst1
  WRITE(17,'(a,1g12.6)') '# sigma_Q=sigma_r,n: 0?=xst2 ', xst2
  
!     cleaning up   (not important, but more clean)
  DO ir=0,nr
    DO iz=0,nz
      wrk1(iz,ir) = 0.0
      wrk2(iz,ir) = 0.0
      wrk3(iz,ir) = 0.0
      wrk4(iz,ir) = 0.0
    END DO
  END DO
  
ELSE IF(icase == 6) THEN
  iz0 = (nz+1)/2+1 ! this is z=0+dz
  WRITE(17,'(/a)') '# ptstcur: plot of currents for Q=p_r sigma_r'
!R250599 test currents sigma_Q,q Q=p_r sigma_r
!     for protons ->wrk1, neutons ->wrk2
  
!     cleaning up
  DO ir=0,nr
    DO iz=0,nz
      wrk1(iz,ir) = 0.0
      wrk2(iz,ir) = 0.0
      wrk3(iz,ir) = 0.0
      wrk4(iz,ir) = 0.0
    END DO
  END DO
  
  DO ist=1,nst(1)+nst(2)
    vi     = SQRT(vsqf(ist))
    ui     = SQRT(one-vi*vi)
    mist   = mpsi(ist)
    npist  = nparz(1,ist)
    IF(ist <= nst(1)) THEN ! protons
      iq     = 1
      eferm  = fer(1)
      jstmax = nst(1)
    ELSE                ! neutrons
      iq     = 2
      eferm  = fer(2)
      jstmax = nst(1)+nst(2)
    END IF
    DO jst=ist,jstmax
      IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
        vj     = SQRT(vsqf(jst))
        uj     = SQRT(one-vj*vj)
        IF(ist /= jst) THEN
          factor = (ui*ui*vj*vj-vi*vi*uj*uj)
          DO ir=0,nr
            DO iz=0,nz
              
! sigma_Q_q = - 2 i sum phi^+ p_r sigma_z phi
              
              IF(iq == 1) THEN ! protons
                wrk1(iz,ir) = -2.0*0.5  &
                    *(dpsir(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                    - dpsir(iz,ir,2,ist)*psi(iz,ir,2,jst)  &
                    + psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)  &
                    - psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst) )  &
                    *factor  +  wrk1(iz,ir)
              ELSE ! neutrons
                wrk1(iz,ir) = -2.0*0.5  &
                    *(dpsir(iz,ir,1,ist)*psi(iz,ir,1,jst)  &
                    - dpsir(iz,ir,2,ist)*psi(iz,ir,2,jst)  &
                    + psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)  &
                    - psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst) )  &
                    *factor  +  wrk2(iz,ir)
              END IF
            END DO
          END DO
        END IF
      END IF
    END DO
  END DO
  
  
  WRITE(17,'(/a)') '# test spinph for ir=1 :',  &
      '# iz  z  spinph_p  exact_p spinph_n  exact_n '
  DO iz=0,nz
    WRITE(17,'(1x,i3,f8.2,4(1pg12.4))')  &
        iz,(iz+sz-icm)*dz,spinph(iz,1,1),wrk1(iz,1), spinph(iz,1,2),wrk2(iz,1)
  END DO
  WRITE(17,'(a)') '# test curqr for z=0+dz :',  &
      '# ir  r  curqr_p  curqr/rho_p  curqr_n  curqr/rho_n '
  DO ir=0,nr
    WRITE(17,'(1x,i3,f8.2,4(1pg12.4))')  &
        ir,(ir+sr)*dr,spinph(iz0,ir,1),wrk1(iz0,ir), spinph(iz0,ir,2),wrk2(iz0,ir)
  END DO
  
  DO ir=0,nr
    DO iz=0,nz
      wrk3(iz,ir)=(wrk1(iz,ir)-spinph(iz,ir,1))**2
      wrk4(iz,ir)=(wrk2(iz,ir)-spinph(iz,ir,2))**2
    END DO
  END DO
  
  WRITE(17,'(/a)') '# ptstcur: plot of currents for Q=p_r sigma_r'
  
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 =  0 ',  &
      curqr(iz0,1,1)
  WRITE(17,'(a,1g12.6)') '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  0 ',  &
      curqz(iz0,1,1)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 =  0 ',  &
      curqr(iz0,1,2)
  WRITE(17,'(a,1g12.6)') '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  0 ',  &
      curqz(iz0,1,2)
  
  xst1 = f1int(wrk3(0:nz,0:nr))
  xst2 = f1int(wrk4(0:nz,0:nr))
  
  WRITE(17,'(/a,1g12.6)') '# sigma_Q=pr sigma_r,p: 0?=xst1 ',xst1
  WRITE(17,'(a,1g12.6)') '# sigma_Q=pr sigma_r,n: 0?=xst2 ', xst2
  
!     cleaning up   (not important, but more clean)
  DO ir=0,nr
    DO iz=0,nz
      wrk1(iz,ir) = 0.0
      wrk2(iz,ir) = 0.0
      wrk3(iz,ir) = 0.0
      wrk4(iz,ir) = 0.0
    END DO
  END DO
  
END IF
RETURN
END SUBROUTINE ptstcur
!R250599 end of currents test output
!############################################################################


!-----seprpa-----------------------------------------------------------

SUBROUTINE seprpa(betain2)

!     separable RPA
!     -------   presently a testing code  --------------------

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
REAL(8),INTENT(IN) :: betain2(ktriac)   

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

!-----------------------------------------------------------------------

IF(kcrank <= 0) RETURN

CALL inie2q(0,1)                           ! 2QP energies on 'e2quas'

!     c.m. test

IF(isym == 0) THEN
  CALL ini_cmz(presp)
  CALL hcomm0(presp,ptest)
  pp   = opovlp(presp,presp)
  ph0p = opovlp(presp,ptest)
  CALL hcommp(presp,qresp)
  php   = opovlp(presp,qresp)
  WRITE(6,'(a,5(1pg12.4))') ' CM TEST: <P*P>,<PHP>,<PH0P>=',pp,php,ph0p
END IF

!     test of 'hcommP' with sum rule model

CALL quadms(qresp,0,1)                    ! initialize Q

CALL hcommq(qresp,ptest,betain2,presp,0)
qq    = opovlp(qresp,qresp)
qhq   = opovlp(qresp,presp)

CALL hcommp(presp,qresp)
pp    = opovlp(presp,presp)
php   = opovlp(presp,qresp)

CALL hcomm0(presp,ptest)
ph0p  = opovlp(presp,ptest)

CALL quadms(qresp,0,1)                    ! initialize Q

CALL hcomm0(qresp,presp)
qh0q  = opovlp(qresp,presp)

CALL hcomm0(presp,ptest)
qh02q = opovlp(presp,ptest)

WRITE(6,'(a,5(1pg12.4))')  &
    ' sumrule estimate: <Q*Q>,<QHQ>,<P*P>,<PHP>,omega=',  &
    qq,qhq,pp,php,SQRT(php/qhq)
WRITE(6,'(a,5(1pg12.4))') '                   <QH0Q>,<PH0P>,omega0=',  &
    qh0q,ph0p,SQRT(qh02q/qh0q)

RETURN
END SUBROUTINE seprpa
!-----hcomm0-----------------------------------------------------------

SUBROUTINE hcomm0(opin,opout)

!     Commutator [H_0,A] of mean-field Hamiltonian with A-response matrix,
!     matrix elements of A are entered on 'opin',
!     result is returned on 'opout'.

!USE paramr
!USE cranking
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'



REAL(8), INTENT(IN)                         :: opin(ktriac)
REAL(8), INTENT(OUT)                        :: opout(ktriac)


!-----------------------------------------------------------------------

ipostr = 0
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  npist  = nparz(1,ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
  ELSE
    jstmax = nst(1)+nst(2)
  END IF
  DO jst=ist,jstmax
    IF(mist == mpsi(jst) .AND. npist == nparz(1,jst)) THEN
      ipostr = 1 + ipostr
      opout(ipostr) = opin(ipostr)*e2quas(ipostr)
    END IF
  END DO
END DO
WRITE(6,'(a,i6)') ' HCOMM=: ipostr=',ipostr

RETURN
END SUBROUTINE hcomm0



!-----ini_cmz-----------------------------------------------------------

SUBROUTINE ini_cmz(popout)

USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
USE linresponse, ONLY: psi1
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'constr.inc'
!INCLUDE 'observ.inc'
!INCLUDE 'dynpar.inc'
!INCLUDE 'fields.inc'

!INCLUDE 'crankin.inc'

REAL(8), INTENT(OUT)                        :: popout(ktriac)

!      data intest/0/
!      data itstpr/1/

!     buffered SQRT as internal statement function
!sqrtb(x) = SQRT(MAX(x,smal))

!-----------------------------------------------------------------------

IF(isym /= 0) STOP ' c.m. test runs only without z-symmetry'

sumwz  = zero
sumez  = zero
sumdir = zero
ipostr = 0
DO ist=1,nst(1)+nst(2)
  mist   = mpsi(ist)
  IF(ist <= nst(1)) THEN
    jstmax = nst(1)
    eferm  = fer(1)
  ELSE
    jstmax = nst(1)+nst(2)
    eferm  = fer(2)
  END IF
  
!     gradient operator Nabla_z on 'psip'
  
  CALL nablaz_spin(psi(0:nz,0:nr,1:2,ist),psi1(0:nz,0:nr,1:2),nparz(1,ist),1,one)
  DO is=1,nsd
    mm        = mist+is-1
    npist  = nparz(is,ist)
    nparr     = 1+MOD(mm,2)
    DO ir=0,nr
      DO iz=0,nz
        psip(iz,ir,is) = psi1(iz,ir,is)
      END DO
    END DO
  END DO
!                                       auxiliary variables
  npist  = nparz(1,ist)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  DO jst=ist,jstmax
    IF(mist == mpsi(jst)) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      pij    = wfovlp(psi(0:nz,0:nr,1:2,jst),psip(0:nz,0:nr,1:2))
      pij    = pij*(vj*ui-uj*vi)
      ipostr = 1 + ipostr
      popout(ipostr) = pij
    END IF
  END DO
END DO
WRITE(6,'(a,i6)') ' C.M. Z-MODE initialized, IPOSTR=',ipostr

RETURN
END SUBROUTINE ini_cmz

!---------------------------------------------------------------------------  
! DESCRIPTION: crinfo
!> @brief
!!Prints quantities from cranking part. 
!--------------------------------------------------------------------------- 
SUBROUTINE crinfo()
USE constr
USE observ
USE dynpar
USE fields
USE levels
USE paramr
USE grid

IMPLICIT NONE

!CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam

REAL(8) :: betawr !,ezp
INTEGER :: ist
!     formfactor and related quantities
CALL radius()
WRITE(7,'(a,4(1pg13.5))') ' testprint R,sigma=',rdmsc,surfc


!     print cranking results
IF(tswicon2) THEN
!  IF(lqcon == 3) THEN
!    betawr = beta3
!    betawp = beta3p
!  ELSE IF(lqcon == 2) THEN
    betawr = beta2
    betawp = beta2p
!  ELSE
!    betawr = LOG(rms)
!    betawp = LOG(rmsp)
  WRITE(7,'(a,2(1pg13.5))')  &
      ' rigid body momentum of inertia and inverse   :', 0.5*xinert,one/xinert
  WRITE(7,'(a,2(1pg13.5))')  &
      ' inv. scaling mass for quadrupole (cart., b_2):', bmassc,bmassr
  IF(crmasp /= zero) WRITE(7,'(a,2(1pg13.5))')  &
      ' Inglis mass for quadrupole (direct,inverse)  :', crmasp,one/crmasp
  IF(crmaso /= zero) WRITE(7,'(a,2(1pg13.5))')  &
      ' off.d. mass for quadrupole (direct,inverse)  :', crmaso,one/crmaso
  IF(crmasq /= zero) WRITE(7,'(a,2(1pg13.5))')  &
      ' 3.Inglis mass for quadrupole (direct,inverse):', crmasq,one/crmasq
  IF(pwid /= zero) WRITE(7,'(a,2(1pg13.5))')  &
      ' GCM width,  ZPE correction                   :', pwid,pwid/(2.0*crmasp)
  IF(pwido /= zero) WRITE(7,'(a,2(1pg13.5))')  &
      ' off diag. GCM width,  ZPE correction         :', pwido,pwido/(2.0*crmaso)
  IF(ammass /= zero) WRITE(7,'(a,3(1pg13.5))')  &
      ' a.m. mass   a.m. width  a.m. ZPE             :',  &
      ammass,amwid,amwid/(2.0*ammass)
  IF(widcmz /= zero) WRITE(7,'(a,4(1pg13.5))')  &
      ' widcm_z&xy, zpe_cm_z&xy                      :',  &
      widcmz,widcmx,zpecmz,zpecmx
  IF(widdz /= zero .AND. tmass) THEN
    OPEN(75,STATUS='unknown',FORM='formatted', &
          POSITION='append',FILE='crankn.'//outnam)

!    IF(scmass.NE.0D0) THEN
    IF(tprecpri) THEN
      WRITE(75,'(1x,15(1pg20.12))')  &
          betawr,scmass,crmasp,crmaso,ammass,amwid,pwid,  &
          energy,zpemod,pnmass,pnmas0,zpetri,ampnms
    ELSE
      WRITE(75,'(1x,f9.5,1pg20.12,5(1pg12.4),0pf13.5,0pf7.3, 4(1pg12.4))')  &
          betawr,scmass,crmasp,crmaso,ammass,amwid,pwid,  &
          energy,zpemod,pnmass,pnmas0,zpetri,ampnms
    END IF
!    END IF
    CLOSE(75)
  END IF
!  IF(tmass .AND. scmass.NE.0D0) THEN
  IF(tmass) THEN
    OPEN(76,STATUS='unknown',FORM='formatted',POSITION='append',&
            FILE='zpecor.'//outnam)
    rchrg2 = rmsp*rmsp + 0.64  - 2.25/(widdz+widdxy)
    WRITE(7,'(a,3f10.3)') ' rmsp2,rchrg2,c.m. r**2=',  &
        rmsp*rmsp,rchrg2,2.25/(widdz+widdxy)
    IF(tprecpri) THEN
      WRITE(76,'(1x,15(1pg20.12))') betawr,pwid/(2.0*scmass),  &
          pwid/(2.0*crmasp),pwido/(2.0*crmaso),amwid/(2.0*ammass),  &
          amwid,energy-zpemod,zpetri,omgtri,rchrg2,betawp,  &
          rdmsc,surfc,betafull,betafullp
    ELSE
      WRITE(76,'(1x,f9.5,1pg20.12,2(1pg12.4),1pg20.12,1pg12.4,  &
          1pg20.12,2(1pg12.4),0pf10.4,5f9.5)') betawr,pwid/(2.0*scmass),  &
          pwid/(2.0*crmasp),pwido/(2.0*crmaso),amwid/(2.0*ammass),  &
          amwid,energy-zpemod,zpetri,omgtri,rchrg2,betawp,  &
          rdmsc,surfc,betafull,betafullp
    END IF
    CLOSE(76)
    IF(kintpl == 1) THEN
      WRITE(77) betawr,scmass,crmasp,crmaso,ammass,amwid,pwido,  &
          energy,zpemod,fer,  &
          (spenrg(ist)-fer(1),vsqf(ist),deltaf(ist)*wtpair(ist), ist=1,nst(1)),  &
          (spenrg(ist)-fer(2),vsqf(ist),deltaf(ist)*wtpair(ist),  &
          ist=nst(1)+1,nst(1)+nst(2))
    END IF
  END IF
END IF

RETURN
END SUBROUTINE crinfo

!---------------------------------------------------------------------------  
! DESCRIPTION: atdhf_path()
!> @brief
!! Computes ATDHF path from simple gradient at the end of all loops
!--------------------------------------------------------------------------- 
SUBROUTINE atdhf_path()
USE constr
USE observ
USE dynpar
USE fields
USE levels
USE paramr
USE grid
USE init, ONLY: idanuc
USE staticstep, ONLY: sstep
USE util, ONLY: plot_2d

IMPLICIT NONE

INTEGER,PARAMETER :: itmax_atdhf=500  ! max. number of points on path
INTEGER :: iter
LOGICAL :: topen173

  IF(.NOT.tatdhf) RETURN

  IF(itmax > itcons2 .AND. iflmb0 == 1) THEN   !  to be generalized
    WRITE(90,'(a)') '#  start direct iteration'
    WRITE(173,'(a)') '# gradient iteration'
    qlamb0=0D0;qlamb2=0D0;qlamb3=0D0;qlamb4=0D0;
    tswicon0=.FALSE.;tswicon2=.FALSE.;tswicon3=.FALSE.;tswicon4=.FALSE.;
    adyrho = 0.5
    adypot = 0.5
    e0inv=0D0
    x0inv=0.0005D0
    CALL idanuc               ! nucleon damping
    DO iter=1,itmax_atdhf   ! itmaxact+1,itmax
      CALL sstep(iter)
      CALL sinfo(iter)
      WRITE(7,'(a,i5,f9.2,2f8.3,f7.3,2f5.1,1pg12.4)')  &
          ' It,E,b2,b3,D^2E,dE/dIt=', iter,energy,  &
          beta2,beta3, devh2,LOG10(MAX(econ,smal))
      CALL flush(7)
      IF(MOD(iter,iplot)==0) THEN
        IF(isym == 0) THEN
          WRITE(90,'(f7.3,11f6.2,f18.11,2f7.1,2f7.2,i5,3f7.3, 2f9.3)')  &
              beta1,beta1p,beta1n,beta2,beta2p,beta2n,  &
              beta3,beta3p,beta3n,beta4,beta4p,beta4n,  &
              energy,epair,fer,iter,rms,rmsp,rmsn,rsepfrag, fragdist(rho(0:nz,0:nr,:))
        ELSE
          WRITE(90,'(f9.5,5f6.2,f18.11,2f7.1,2f7.2,i5,3f7.3,f11.3,  &
              &4f8.3,1pg13.5)') beta2,beta2p,beta2n,beta4,beta4p,beta4n,  &
              energy,epair,fer,iter,rms,rmsp,rmsn,rsepfrag, fragdist(rho(0:nz,0:nr,:))
        END IF
        CALL plot_2d()
        INQUIRE(UNIT=173,OPENED=topen173)
        IF(topen173) THEN
          WRITE(173,'(2f7.3,f10.3,2(1pg13.5))')  &
              beta2,fragdist(rho(0:nz,0:nr,:)),energy,rho(icm,0,1),rho(icm,0,2)
          CALL flush(173)
        END IF
      END IF
    END DO
  END IF
END SUBROUTINE atdhf_path


!---------------------------------------------------------------------------  
! DESCRIPTION: trizpe
!> @brief
!!Zero point energy for triaxial vibrations, deduced simply from quadrupole 
!!momentum.
!--------------------------------------------------------------------------- 
SUBROUTINE trizpe

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT NONE

INTEGER,PARAMETER :: intest=0

REAL(8) :: deltai,deltaj,edifi,edifj,eferm,eqainv,equasi,equasj,pij,r,&
           r2,sumq,sumqe,sumqi,sumqi2,sumqi3,ui,uj,vi,vj,volq
INTEGER :: ir,is,ist,jst,jstmax,jstmin,mist,mistdw,mistup,npist
!     double loop to accumulate r**2 width with 'm' stepping by +-2

REAL(8),ALLOCATABLE :: psip(:,:,:)

ALLOCATE(psip(0:nz,0:nr,nsd))
IF(intest == 1) WRITE(7,'(a)') ' Nabla_xy:'
sumq   = zero
sumqe  = zero
sumqi  = zero
sumqi2 = zero
sumqi3 = zero
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
  
!     transition operator r**2 on 'psip'
  DO is=1,nsd
    DO ir=0,nr
      r    = (ir+sr)*dr
      r2   = r*r
      psip(0:nz,ir,is) = r2*psi(0:nz,ir,is,ist)
    END DO
  END DO
!                                       auxiliary variables
  npist  = nparz(1,ist)
  mistup = mist+2
  mistdw = ABS(mist-2)
  vi    = sqrtb(vsqf(ist))
  ui    = sqrtb(one-vi*vi)
  deltai = deltaf(ist)*wtpair(ist)
  edifi  = spenrg(ist)-eferm
  equasi = SQRT(edifi*edifi+deltai*deltai)
!                                       the overlaps
  DO jst=jstmin,jstmax
    IF(npist == nparz(1,jst) .AND.  &
          (mistup == mpsi(jst) .OR. mistdw == mpsi(jst)) ) THEN
      vj    = sqrtb(vsqf(jst))
      uj    = sqrtb(one-vj*vj)
      deltaj = deltaf(jst)*wtpair(jst)
      edifj  = spenrg(jst)-eferm
      equasj = SQRT(edifj*edifj+deltaj*deltaj)
      pij    = wfovoc(psi(0:nz,0:nr,1,jst),psip(0:nz,0:nr,1))  &
          + wfovoc(psi(0:nz,0:nr,2,jst),psip(0:nz,0:nr,2))
      volq   = (pij*pij)*(ui*vj+uj*vi)**2
      eqainv  = 1.0/(equasi+equasj)
      sumq   = volq                        + sumq
      sumqe  = volq*(equasi+equasj)        + sumqe
      sumqi  = volq*eqainv                 + sumqi
      sumqi2 = volq*eqainv*eqainv          + sumqi2
      sumqi3 = volq*eqainv*eqainv*eqainv   + sumqi3
    END IF
  END DO
END DO
zpetri = 0.5*sumqi2/sumqi3
omgtri = SQRT(sumqi/sumqi3)
WRITE(7,'(a,5(1pg12.4))') ' trizpe: sumQ,sumQe,sumQi,sumQi2,sumQi3=',  &
    sumq,sumqe,sumqi,sumqi2,sumqi3
WRITE(7,'(a,3(1pg12.4))')  &
    '      ZPE,omega,GR=',zpetri,omgtri,SQRT(sumqe/sumqi)

DEALLOCATE(psip)
RETURN
END SUBROUTINE trizpe


END MODULE cranking
