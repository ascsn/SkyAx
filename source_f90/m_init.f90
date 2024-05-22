!------------------------------------------------------------------------------
! MODULE: init
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains routines to initialize output files, several parameters 
!!and variables, and wave functions. 
!------------------------------------------------------------------------------
MODULE init

CONTAINS

!---------------------------------------------------------------------------  
! DESCRIPTION: ini_files
!> @brief
!! Initialises several output files.
!> @details
!!
!--------------------------------------------------------------------------- 
SUBROUTINE ini_files()

USE paramr
USE grid, ONLY: tsym
USE constr
USE dynpar
USE levels
USE observ, ONLY: tempac,tempmn,tempw,tempmx

IMPLICIT NONE

  OPEN(90,STATUS='unknown',FORM='formatted',FILE='ensurf.'//outnam)
  IF (itert==0 .AND. tempac == tempmn) THEN
    IF(ipair /= 4) THEN
      WRITE(90,'(a)') '#'//title(4:80)
      WRITE(90,1098) '#'//title(1:70),npart(1)+npart(2),npart(1),  &
          dz,dr,nz,nr, addrho,addpot,addocc,c0,qeps,  &
          x0inv,e0inv
      1098 FORMAT('# Plot potential energy surface for:',/,1X,a,/,  &
          '# Nmass=',i5,', Ncharge=',i5,/,  &
          '#  dz,dr: ',2F6.3,' nz,nr:',2I4,/,  &
          '# addrho,addpot,addocc: ',3F6.3,'  c0,qeps:',2F6.3,/,  &
          '# x0inv(fm),e0inv(mev): ',2F7.2)
      IF(istblo > 0) WRITE(90,'(a,4i5)')  &
          '# blocking: kblock,mblock,isoblo,istblo=',  &
          kblock,mblock,isoblo,istblo
      IF(.NOT.tsym) THEN
        WRITE(90,'(3a)')  &
            '# bet1   bet1p  bet1n  bet2   bet2p  bet2n  bet3   bet3p  bet3n', &
            '  bet4   bet4p  bet4n   energy           epairp epairn ferp   fern', &
            '  iter   rmsc    rms    rmsp    rmsn     N_boundp  N_boundn'
      ELSE
        WRITE(90,'(3a)') '# bet2   bet2p  bet2n  ',  &
            'bet4   bet4p  bet4n     energy         epairp epairn ferp   fern  iter',  &
            '  rmsc    rms      rmsp    rmsn    N_boundp  N_boundn '
      END IF
    ELSE
      WRITE(90,1092) title,npart(1)+npart(2),npart(1),  &
          tempmn,tempmx,tempw, dz,dr,nz,nr,  &
          addrho,addpot,addocc,c0,qeps, x0inv,e0inv
      1092 FORMAT('# Plot potential energy surface for:',/,a,/,  &
          '# Nmass=',i5,', Ncharge=',i5/,  &
          '# tempmn,tempmx,tempw=',3F6.3,/,  &
          '#  dz,dr: ',2F6.3,' nz,nr:',2I4,/,  &
          '# addrho,addpot,addocc: ',3F6.3,'  c0,qeps:',2F6.3,/,  &
          '# x0inv(fm),e0inv(mev): ',2F7.2,/,  &
          '#   T    beta2   beta2p    beta4   beta4p   ',  &
          '#   energy     efree         T*S   level-d.')
    END IF
  END IF
  CLOSE(90)
  IF(ipri_fields>0) OPEN(18,STATUS='unknown',FORM='formatted',FILE='fields.'//outnam)
  IF(tpri_conver) THEN
    OPEN(71,STATUS='unknown',FORM='formatted',FILE='conver.'//outnam)
    WRITE(71,'(a,a)') '  Plot convergence properties for: ',title
  END IF
  IF(itest >= 1) THEN
    OPEN(73,STATUS='unknown',FORM='formatted',FILE='constraint_damped_details'//outnam)
    OPEN(74,STATUS='unknown',FORM='formatted',FILE='constraint_details'//outnam)
  END IF
  IF(tpri_radii)  THEN
    IF(itert==0) THEN
      OPEN(76,STATUS='unknown',FORM='formatted',FILE='radii.'//outnam)
      WRITE(76,'(a)') '# Print charge form  parameters:',&
        '#  beta_2    energy   r_rms    r_diffr  sigma   r_4'
      CLOSE(76)
    END IF
  END IF

  IF(tpri_formf)  THEN
    IF(itert==0) THEN
      OPEN(176,STATUS='unknown',FORM='formatted',FILE='formfactor.'//outnam)
      WRITE(176,'(a)') '# Print charge formfactors and charge distribution:'
      !CLOSE(76)
    END IF
  END IF


END SUBROUTINE ini_files

!---------------------------------------------------------------------------  
! DESCRIPTION: inipai
!> @brief
!! Initialises pairing parameters. 
!> @details
!!
!--------------------------------------------------------------------------- 
SUBROUTINE inipai(nmass)   

USE paramr
USE constr
USE dynpar
USE levels
USE grid
USE observ
USE pairin, ONLY : rspln,suv2rs,ilnrsp,wcoltt

IMPLICIT NONE

INTEGER, INTENT(IN) :: nmass

IF((ipair==9.OR.ipair==11) .AND. .NOT. ALLOCATED(rspln)) THEN
  ALLOCATE(rspln(1:2),suv2rs(1:2),wcoltt(0:nz,0:nr))
  rspln=0D0;suv2rs=0D0;wcoltt=0D0
END IF

fer(1:2)=0D0                             ! Fermi energy
feradd(1:2)=0D0                          ! Cutoff "Fermi energy"
gap(1:2)=11.2 / SQRT(one*(nmass))        ! estimate for gap

IF(ipair.GE.0 .AND. ipair.LE.3) THEN
  CALL user_interface_3()          ! initializes user defined pairing 
ELSE IF(ipair==4) THEN
  pairfc=0D0
ELSE IF(ipair == 5 .OR.ipair == 6 .OR. ipair == 9) THEN
  IF(pairfc(1)==0D0) pairfc(1) = 300D0
  IF(pairfc(2)==0D0) pairfc(2) = 300D0
ELSE IF(ipair == 7 .OR. ipair == 8 .OR. ipair == 11) THEN
  IF(pairfc(1)==0D0) pairfc(1) = 600D0
  IF(pairfc(2)==0D0) pairfc(2) = 600D0
ELSE
  WRITE(*,*) ipair
  STOP "invalid option for IPAIR"
END IF
IF(ipair == 9 .OR. ipair == 11) THEN
  rspln(1:2)  = 0.0
  suv2rs(1:2) = 1.0
  IF(ilnrsp == 1) THEN
    WRITE(7,'(a)') ' INCLUDING linear response from L.N.'
  ELSE
    WRITE(7,'(a)') ' WITHOUT linear response from L.N.'
  END IF
END IF



RETURN
END SUBROUTINE inipai

!---------------------------------------------------------------------------  
! DESCRIPTION: inipotdamp
!> @brief
!! Initializes the additional potential to ensure value and derivative equal 
!! to zero in the neighbourhood of the boundaries on \c potdamp.
!> @details
!!
!--------------------------------------------------------------------------- 

SUBROUTINE inipotdamp

USE paramr
USE constr
USE dynpar
USE fields
USE levels
USE grid
USE observ
IMPLICIT NONE

INTEGER :: ir,iz
REAL(8) :: sph_rsq,sph_rmax,sph_rdamp,dlt
REAL(8) :: rmax,zmax1,zmax2,zmin,zmax
REAL(8) :: damp_rad,damp_hgt

!     set parameters of damping field

damp_rad = 5.0D0     ! width of damping field in fm
damp_hgt = 100.0D0   ! strength of damping field in MeV

!------------------------------------------------------------------
!     initialize damping field
!------------------------------------------------------------------
rmax=(nr+sr)*dr
zmax1=(nz+sz-icm)*dz
zmax2=-(sz-icm)*dz

IF (.NOT.tsym) THEN
  IF (zmax1 >= zmax2) THEN     ! find maximum z-value
    zmax=zmax1
    zmin=zmax2
  ELSE
    zmax=zmax2
    zmin=zmax1
  END IF
ELSE
  zmin = zmax1
END IF

IF (zmin <= rmax) THEN          ! find maximum radius
  sph_rmax = zmin
ELSE
  sph_rmax = rmax
END IF

sph_rdamp = sph_rmax - damp_rad

DO ir=0,nr
  DO iz=0,nz
    sph_rsq   = SQRT(r(ir)**2+z(iz)**2)
    dlt       = sph_rsq-sph_rdamp
    
    IF (sph_rsq <= sph_rdamp) THEN
      potdamp(iz,ir) = 0.0D0
    ELSE IF (sph_rsq >= sph_rmax) THEN
      potdamp(iz,ir) = 1.0D0
    ELSE
      potdamp(iz,ir) = 1.0D0-EXP(dlt**2/(dlt**2-damp_rad**2))
    END IF
    
  END DO
END DO

!     weight damping field
potdamp(0:nz,0:nr) = damp_hgt*potdamp(0:nz,0:nr)

RETURN
END SUBROUTINE inipotdamp


!---------------------------------------------------------------------------  
! DESCRIPTION: idanuc
!> @brief
!! In this subroutine the matrices needed for the damping are initialised 
!! using the derivative matrices of the subroutine \c inider.
!> @details
!! The damping operator \c d is given by:
!! \f[ d=\frac{\texttt{x0inv}}{\texttt{e0inv}}\frac{1}{1+t_r/\texttt{e0inv}}
!!\frac{1}{1+t_z/\texttt{e0inv}} \f]
!!where \f$ t_r \f$ is the kinetic energy in \f$r\f$-direction and \f$t_z\f$ in
!!\f$z\f$-direction in axial symmetry.
!!
!!Output via \c zdnuc (damping matrix, \f$z\f$-direction) and 
!!\c rdnuc (damping matrix, \f$r\f$-direction), angular momentum term is 
!!included.
!--------------------------------------------------------------------------- 
SUBROUTINE idanuc

USE paramr
USE dynpar
USE wavefr
USE grid
USE observ
IMPLICIT NONE

REAL(8) :: befac,errmax,sum
INTEGER :: i,impsi,ipi,ir,iz

REAL(8), ALLOCATABLE :: zaux(:,:),raux(:,:)
REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:),wrk3(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr),wrk3(0:nz,0:nr))
ALLOCATE(zaux(0:nz,0:nz),raux(0:nr,0:nr))
wrk1=0D0;wrk2=0D0;wrk3=0D0;zaux=0D0;raux=0D0

!     set up unit matrix for 'e0inv=0'

IF(e0inv == zero) THEN
  zdnuc=zero
  rdnuc=zero
  DO i=0,nz
    zdnuc(i,i,:) = x0inv
  END DO
  DO i=0,nr
    rdnuc(i,i,:) = x0inv
  END DO
ELSE
!     now the standard case for damped iteration
  befac  = h2mav/(two*e0inv)
  IF(itest>1) WRITE(7,'(a,g12.4,a)') ' test idanuc: befac=',befac,' , derz2:'

!     nucleon damping in z-direction
  zdnuc(0:nz,0:nz,:) = -befac * derz2f(0:nz,0:nz,:)
  DO i=0,nz
    zdnuc(i,i,:) = one + zdnuc(i,i,:)
  END DO
    
!     inversion yields final result for zdnuc
  DO ipi=1,4    
    CALL matr_inv(zdnuc(0:nz,0:nz,ipi),nz+1,nz+1)
    IF (ipi == 2.OR.ipi == 4) THEN
      zdnuc(0,0:nz,ipi) = 0.0D0
      zdnuc(0:nz,0,ipi) = 0.0D0  ! addition PK
    END IF
    
  END DO
  
!     test damping
  IF(itest == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        wrk1(iz,ir) = MAX(1-ABS(iz-ir),0)
        wrk2(iz,ir) = wrk1(iz,ir)
      END DO
    END DO
    CALL mulz(derz2f(0:nz,0:nz,1),wrk1(0:nz,0:nr),wrk2(0:nz,0:nr),-befac,1)
    CALL mulz(zdnuc(0:nz,0:nz,1),wrk2(0:nz,0:nr),wrk3(0:nz,0:nr),one,0)
    sum    = 0.0
    errmax = 0.0
    DO ir=0,nr
      DO iz=0,nz
        sum    = (wrk3(iz,ir)-wrk1(iz,ir))**2 + sum
        errmax = MAX(ABS(wrk3(iz,ir)-wrk1(iz,ir)),errmax)
      END DO
    END DO
    WRITE(7,'(a,3(/a,10(1pg20.12)))') ' test z-damping:',  &
        ' input:',(wrk1(iz,1),iz=0,9), ' kin. :',(wrk2(iz,1),iz=0,9),  &
        ' test :',(wrk3(iz,1),iz=0,9)
    WRITE(7,'(a,2(1pg13.5))')  &
        ' error: average,maximum=',SQRT(sum/((nr+1)*(nz+1))),errmax
  END IF

!     the nucleon damping in r-direction
  rdnuc(0:nr,0:nr,0:kmpsi) = -befac * derr2af(0:nr,0:nr,0:kmpsi,1)
  DO i=0,nr
    rdnuc(i,i,0:kmpsi) =  one + rdnuc(i,i,0:kmpsi)
  END DO
  
  DO impsi=0,kmpsi
  
    CALL matr_inv(rdnuc(0:nr,0:nr,impsi),nr+1,nr+1)
    IF (impsi > 0) THEN
      rdnuc(0,0:nr,impsi) = zero
      rdnuc(0:nr,0,impsi) = zero ! addition PK
    END IF
    
!       test damping
    IF(itest == 1) THEN
      DO ir=0,nr
        DO iz=0,nz
          wrk1(iz,ir) = MAX(1-ABS(iz-ir),0)
          wrk2(iz,ir) = wrk1(iz,ir)
        END DO
      END DO
      
      CALL mulr(derr2af(0:nr,0:nr,impsi,1),wrk1(0:nz,0:nr),wrk2(0:nz,0:nr),-befac,1)
      CALL mulr(rdnuc(0:nr,0:nr,impsi),wrk2(0:nz,0:nr),wrk3(0:nz,0:nr),one,0)
      sum    = 0.0
      errmax = 0.0
      DO ir=0,nr
        DO iz=1,nz
          sum    = (wrk3(iz,ir)-wrk1(iz,ir))**2 + sum
          errmax = MAX(ABS(wrk3(iz,ir)-wrk1(iz,ir)),errmax)
        END DO
      END DO
      WRITE(7,'(a,i3,a,3(/a,10(1pg20.12)))')  &
          ' test r-damping for m=',impsi,':', ' input:',(wrk1(iz,1),iz=0,9),  &
          ' kin. :',(wrk2(iz,1),iz=0,9), ' test :',(wrk3(iz,1),iz=0,9)
      WRITE(7,'(a,2(1pg13.5))')  &
          ' error: average,maximum=',SQRT(sum/((nr+1)*(nz+1))),errmax
    END IF
  END DO
END IF

DEALLOCATE(zaux,raux,wrk1,wrk2,wrk3)
RETURN
END SUBROUTINE idanuc

!---------------------------------------------------------------------------  
! DESCRIPTION: idames
!> @brief
!! This subroutine initialises the damping matrices for the meson damping.
!> @details
!! The meson damping is given by:
!! \f[ d=\frac{1}{\texttt{dams}^2}\frac{1}{1-\Delta_r/\texttt{dams}^2}
!!\frac{1}{1-\Delta_z/\texttt{dams}^2} \f]
!!where \c dams is a numerical damping mass and \f$\Delta_{r/z}\f$ is the
!! delta operator in \f$r/z\f$ direction in axial symmetry.
!!
!!Output via \c rdmes (damping matrix, \f$r\f$-direction) and 
!!\c zdmes (damping matrix, \f$z\f$-direction).
!--------------------------------------------------------------------------- 

SUBROUTINE idames

USE paramr
USE dynpar
USE wavefr
USE grid
USE observ
USE coulomb, ONLY: dmsfac,cmssqr
IMPLICIT NONE

REAL(8) :: tms
INTEGER :: i,imes

imes=1                    ! only coulomb is initialized now

tms  = -dmsfac/cmssqr

rdmes(0:nr,0:nr,imes) =  tms  * derr2af(0:nr,0:nr,0,1)
DO i=0,nr
  rdmes(i,i,imes)  =  rdmes(i,i,imes) + one
END DO
CALL matr_inv(rdmes(0:nr,0:nr,imes),nr+1,nr+1)

zdmes(0:nz,0:nz,imes) = tms * derz2f(0:nz,0:nz,1)
DO i=0,nz
  zdmes(i,i,imes) = zdmes(i,i,imes) + one
END DO

CALL matr_inv(zdmes(0:nz,0:nz,imes),nz+1,nz+1)

RETURN
END SUBROUTINE idames

!---------------------------------------------------------------------------  
! DESCRIPTION: ssetup
!> @brief
!! Initialises the \c psi with harmonic oscillator wave functions. 
!> @details
!! The routine is run twice. Once to determine the number of states and the 
!! second time to initialize wave functions and levels. After the first run 
!!the arrays are allocated with the exact number of states.
!--------------------------------------------------------------------------- 

SUBROUTINE ssetup

USE paramr
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
USE observ
USE steputil, ONLY: vsqtr

IMPLICIT NONE

REAL(8) :: amas,brosc,bzosc,ecur,hboma,hbomr,hbomz,octadd,octstp,phinor,&
           psirad,spin,sumnon,sumprt,sumtot,v2,wgt,xpi
INTEGER :: ielse,iofk,iq,iqold,ir,is,isigis,ispgen,ist,iz,izm,izp,j,&
           jst,k1,kold,mm,mm1,mmaxq,mmr,mst,muplo,n,na,nnr,nnz,&
           npuplo,nzho,len,irun

REAL(8) ::  xkappa(2),xmue(2),spnc(2),eocct(2)
REAL(8) ::  phi(nnzxd+1,0:nz)
LOGICAL :: init_run

REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
wrk1=0D0;wrk2=0D0

IF(itest>0) WRITE(7,'(a)') ' SSETUP entered'

!**********************************************************************
!                                                                     *
!        setup the deformed harmonic oscillator parameters in terms   *
!        of the mass and initial deformation parameter betain         *
!        compute the oscillator sizes and frequencies in the r and z  *
!        directions and scale eocc by the frequency                   *
!                                                                     *
!**********************************************************************

xpi=SQRT(4.0/(5.0*pi))
amas=npart(1)+npart(2)
hboma=41.0/amas**(edr)
hbomz=hboma*EXP(-xpi*betain)
hbomr=hboma*EXP(xpi/two*betain)
brosc=SQRT(two*h2mav/hbomr)
bzosc=SQRT(two*h2mav/hbomz)
eocct(1:2)=eocc(1:2)*hboma
IF(tsym) THEN
  octstp = zero
ELSE
  octstp = octini/(amas)     !Estimate to reproduce actual beta_3
END IF
IF(itest>0) WRITE(6,'(a,2g12.4)') ' octini,octstp=',octini,octstp

!**********************************************************************
!                                                                     *
!        parameters for l*s and l**2 force in nilsson model           *
!                                                                     *
!**********************************************************************

IF(npart(1) <= 50) THEN
  xkappa(1)=0.08
  xmue(1)  =0.0
ELSE IF(npart(1) <= 82) THEN
  xkappa(1)=0.0637
  xmue(1)  =0.0
ELSE
  xkappa(1)=0.0577
  xmue(1)  =0.0
END IF

IF(npart(2) <= 50) THEN
  xkappa(2)=0.08
  xmue(2)  =0.0
ELSE IF(npart(2) <= 82) THEN
  xkappa(2)=0.0637
  xmue(2)  =0.0
ELSE IF(npart(2) <= 126) THEN
  xkappa(2)=0.0637
  xmue(2)  =0.0
ELSE
  xkappa(2)=0.0635
  xmue(2)  =0.0
END IF
!**********************************************************************
!                                                                     *
!        initialize the harmonic oscillator wavefunctions             *
!                                                                     *
!**********************************************************************

IF(tsym) THEN
  nzho=nz
ELSE
  nzho=(nz/2)+1
END IF
DO j=0,nzho
  CALL zwf(j,bzosc,phi)
END DO
!**********************************************************************
!                                                                     *
!        set loop to allocate the single particle states in terms     *
!        of axially symmetric harmonic oscillators                    *
!        new variables are as follows:                                *
!                      lst     -  total no. of states                 *
!                      spn(*)  -  normalization of wavefunctions      *
!                      iofk    -  pointer to the 1 st index of a      *
!                                 state of charge q                   *
!                      lbuf    -  buffer index for direct i/o         *
!                                                                     *
!**********************************************************************
IF(ALLOCATED(wate)) THEN
  CALL dealloc_levels
  CALL dealloc_wf
END IF
  
init_run=.FALSE.
DO irun=1,2
IF(irun==2) init_run=.TRUE.
lst=0
iofk=1
kold=0
iqold=1
nom=1
mhi(1)=1
mlo(1)=1
miq(1)=iqold
mim(1)=kold
jst=0
mst=0
IF(init_run) THEN
  psi = zero    ! zero wavefn
  spn(:)=0.0
END IF

!**********************************************************************
!                                                                     *
!        main loop over charge states begins here                     *
!        for each q mmax(q) is the maximum number of allowed  m       *
!        values                                                       *
!        note: m is k-1/2, k being the jz quantum number              *
!                                                                     *
!**********************************************************************

DO  iq=1,2
  mmaxq=mmax(iq)
  DO  mm=1,mmaxq
    k1=mm-1
    
!**********************************************************************
!                                                                     *
!        nnz and nnr are the principle quantum numbers of the         *
!        cylindrical oscillators- nnzxd and nnrxd are their maximum   *
!        allowed values                                               *
!                                                                     *
!**********************************************************************
    
    DO  nnz=1,nnzxd
      DO  nnr=1,nnrxd
        
!**********************************************************************
!                                                                     *
!        ispgen counts the two possible spin orientations per         *
!        given mm. ispgen=1 occupies dominantly the upper             *
!        component whereas ispgen=2 occupies the lower one.           *
!        small lower, or upper, components are added by per-          *
!        turbative treatment of the l*s force.                        *
!                                                                     *
!**********************************************************************
        
        DO  ispgen=1,2
          spin=1.5-(ispgen)
          mm1=k1+ispgen-1
          mmr=mm1+1
          
!**********************************************************************
!                                                                     *
!         ecur is the single particle energy of the deformed          *
!         oscillator state plus l*s and l**2 force.                   *
!         this must be less than eocc for each charge state           *
!         to be included in the basis                                 *
!         new variables as follows:                                   *
!                       mpsi(*)  -  jz quantum number minus 1/2       *
!                                    note that upper comp. has m=mpsi *
!                                    but lower comp. has m=mpsi+1     *
!                       npar(*)  -  z parity of the upper comp.       *
!                       wate(*)  -  occupation numbers in the         *
!                                   absence of pairing                *
!                       icharg(*)   isospin quantum number            *
!                                                                     *
!**********************************************************************
          
          ecur= hbomz*(nnz- 0.5)+ hbomr*(2*nnr+ mmr-2)  &
              -(spin+spin)*hboma*mm1*xkappa(iq) -hboma*mm1*mm1*xkappa(iq)*xmue(iq)
          
!**********************************************************************
!                                                                     *
!        if block for selecting the initial occupied states           *
!                                                                     *
!**********************************************************************
          
          isigis  = 3-2*ispgen
          IF(ecur <= eocct(iq)) THEN
            mst         = mst+1
            lst         = lst+1
            jst         = jst+1
            IF(init_run) THEN
              spenrg(lst) = (ecur-70.0)
              mpsi(lst)   = k1
              nqr(lst)    = nnr
              nqz(lst)    = nnz
              IF (k1+1 > kmpsi) THEN
                WRITE(7,'(a,i3)')  &
                    ' kmpsi (dim of debe for ang mom) too small:',kmpsi
                STOP
              END IF
              npar(lst)   = isigis*(MOD(nnz,2)*2-1)
              IF(tsym)THEN
                nparz(1,lst)= (3-npar(lst))/2
                nparz(2,lst)= 3-nparz(1,lst)
              ELSE                                    ! djd addition
                nparz(1,lst)= 1                       ! djd addition
                nparz(2,lst)= 1                       ! djd addition
              END IF
              nrpar(lst)  = (1-MOD(k1,2)*2)
              icharg(lst) = iq
              wate(lst)   = two
              wtpair_uv(lst)  = 1D0
              wtpair_v2(lst)  = 1D0
            
!**********************************************************************
!                                                                     *
!        construction of the wave function, dominant spinor part      *
!                                                                     *
!**********************************************************************
            
              phinor=one/SQRT(pi+pi)
              DO ir=0,nr
                CALL radwf(ir,brosc,nnr,mmr,psirad)
                DO iz=0,nz
                  izp=iz-icm     !+1          cbm
                  izm=icm-iz     !+2-ispgen   cbm
                  octadd = octstp*z(iz)*(z(iz)**2-1.5*r(ir)**2)
                  IF(.NOT.tsym) THEN
                    IF(iz < icm) THEN
                      psi(iz,ir,ispgen,jst) = (one-octadd)  &
                          *FLOAT(MOD(nnz,2)*2-1)*psirad*phi(nnz,izm)
                    ELSE
                      psi(iz,ir,ispgen,jst) = (one-octadd)*psirad*phi(nnz,izp)
                    END IF
                  ELSE
                    psi(iz,ir,ispgen,jst) = psirad*phi(nnz,iz)
                  END IF
                  psi(iz,ir,ispgen,jst)=phinor*psi(iz,ir,ispgen,jst)
                END DO
              END DO
              
              IF(itest >= 1) THEN
                WRITE(7,*) ' 1. time in ssetup: lst,jst,k1,iq=', lst,jst,k1,iq
                IF (ispgen == 1) THEN
                  CALL test_print(psi(0:nz,0:nr,1,jst),0,nz_testp,0,nr_testp,' psi 1.comp')
                ELSE IF (ispgen == 2) THEN
                  CALL test_print(psi(0:nz,0:nr,2,jst),0,nz_testp,0,nr_testp,' psi 2.comp')
                END IF
              END IF
              ielse = 3 - ispgen
              psi(0:nz,0:nr,ielse,jst)=zero
            
!**********************************************************************
!                                                                     *
!        calculate the normalization of the single particle state     *
!        and renormalize the wavefunctions                            *
!                                                                     *
!**********************************************************************
            
!        projection of the last component in momentum space

              muplo  =  mpsi(jst)+ispgen-1
              npuplo =  npar(jst)*isigis
              IF (npuplo > 0 .OR. MOD(muplo,2) == 0) THEN
                IF (muplo > 0) THEN
                  psi(0:nz,0,ispgen,jst) = 0.0d0
                ELSE IF (npuplo < 0) THEN
                  psi(0,0:nr,ispgen,jst) = 0.0d0
                END IF
                IF (MOD(muplo,2) == 1) THEN
                  psi(0:nz,nr,ispgen,jst) = 0.0d0
                END IF
              END IF
            
              spn(lst)=zero
              DO  is=1,2
                spnc(is) = wfovpc(psi(0:nz,0:nr,is,lst))
              END DO
              spn(lst)=spnc(1)+spnc(2)
              cmpn(1,lst)=spnc(1)/spn(lst)
              cmpn(2,lst)=spnc(2)/spn(lst)
              spnc(1)=1.0/spn(lst)
              spnc(2)=spnc(1)
              IF(itest >= 2) THEN
                WRITE(7,'(a,i3/3(a,f10.4/))') ' lst=',lst,' norm,1.comp',  &
                    cmpn(1,lst),' norm,2.comp',cmpn(2,lst),' tot norm',spn(lst)
              END IF
              psi(0:nz,0:nr,1:2,jst)=psi(0:nz,0:nr,1:2,jst)/SQRT(spn(lst))
              spn(lst)=SQRT(spn(lst))
              IF(itest >= 2) THEN
                WRITE(7,*) ' after proj. in ssetup: lst,jst,k1,iq=',  &
                    lst,jst,k1,iq
                CALL test_print(psi(0:nz,0:nr,1,jst),0,9,0,9,' psi 1.comp')
              END IF
            
!**********************************************************************
!                                                                     *
!        the indices for a wavefunctins having a definite value of    *
!        mm must be stored. each buffer is comprized of complete      *
!        blocks of wavefunctions for a particular mm as follows       *
!        mhi(nom) - array storing the upper index for the             *
!                        wavefunctions in the buffer lbuf             *
!        mlo(nom) - array storing the lower index for the             *
!                        wavefunction in the buffer lbuf              *
!        miq(nom) - array storing the iq value for the wave-          *
!                        function in the buffer lbuf                  *
!        mim(nom) - array storing the m value for the wave-           *
!                        function in the buffer lbuf                  *
!        nom      - number of mm blocks in the buffer lbuf            *
!                                                                     *
!        write out the wavefunctions to disk if the buffer is full    *
!                                                                     *
!**********************************************************************
              IF (k1 /= kold.OR.iq /= iqold) THEN
                mhi(nom)=jst-1
                nom=nom+1
                mlo(nom)=jst
                kold=k1
                iqold=iq
                miq(nom)=iqold
                mim(nom)=kold
              END IF
              nbuf(jst) = nom
            
            
!**********************************************************************
!                                                                     *
!        end if block for allowed quantum numbers                     *
!        end nnz and nnr loops                                        *
!                                                                     *
!**********************************************************************
            END IF !END block which is not executed first time (determine nst)
          END IF
        END DO
      END DO
    END DO
    
!**********************************************************************
!                                                                     *
!        end loop over m quantum number                               *
!        set the number of single particle states in nst              *
!                                                                     *
!**********************************************************************
    
  END DO
  nst(iq)= lst
  iokk(iq)=iofk
  iofk=iofk+lst
  
!**********************************************************************
!                                                                     *
!         end loop over iq                                            *
!                                                                     *
!**********************************************************************
  
END DO
IF(.NOT.init_run) THEN
  CALL alloc_levels(lst)
  CALL alloc_wf(lst)
END IF
END DO
npsi(2)=nst(2)
npsi(1)=nst(1)
nst(2)=nst(2)-nst(1)
mhi(nom)=jst
LEN=mhi(nom)
jst=0


!**********************************************************************
!                                                                     *
!        calculate the occupation numbers and the ordering of the     *
!        levels. note:                                                *
!                     iofk is the pointer to the 1 st index           *
!                     of state q                                      *
!                                                                     *
!**********************************************************************

DO iq=1,2
  iofk=iokk(iq)
  CALL vsqtr(vsqf(iofk),wate(iofk),nst(iq),  &
      spenrg(iofk),npart(iq),iorder(iofk))
END DO
ifocsv = 0                                   ! reset flag
WRITE(6,'(a)') ' reset IFOCSV=0'

!     print levels

WRITE(7,'(/a,2f8.3//18x,a,12x,a/a,a)')  &
    ' The initial levels, chosen up to EOCC(1/2)=',eocc,  &
    'parities','components', '   Nr Iso  j_z total  r    z   v^2  1.   2.',  &
    ' osc.quant'
DO ist=1,LEN
  WRITE(7,'(i5,i4,f6.1,3i5,f6.3,2f5.1,f7.3)')  &
      ist,icharg(ist),mpsi(ist)+0.5,npar(ist),nrpar(ist), nparz(1,ist),  &
      vsqf(ist),cmpn(1,ist),cmpn(2,ist),(spenrg(ist)+70.0)/hboma
END DO
WRITE(7,'(a)') ' '

!**********************************************************************
!                                                                     *
!        initial pairing gap for constant strength models             *
!                                                                     *
!**********************************************************************
!     test particle number

na = 0
DO iq=1,2
  sumprt = zero
  sumnon = zero
  sumtot = zero
  DO n=1,nst(iq)
    na  = 1 + na
    v2  = vsqf(na)
    wgt = wate(na)
    sumtot = wgt + sumtot
    sumprt = v2*wgt + sumprt
    sumnon = (one-v2)*wgt + sumnon
  END DO
  WRITE(6,'(a,i2,a,g12.4/14x,2(a,g12.4))')  &
      ' after ssetup: for type=',iq,' nr. of particles=',sumprt,  &
      ' nr. of empty places=',sumnon, ' nr. of possible particles=',sumtot
!  IF(sumtot < npart(iq)-ERR)  &
  IF(sumtot < npart(iq))  &
      STOP ' not enough states to reach particle number'
!old        if((ipair.eq.5 .or. ipair.eq.7) .and. sumtot.lt.npart(iq)*two)
  IF(ipair >= 5 .AND. sumtot < npart(iq)*two)  &
      STOP ' IPAIR=5 or 7 requires larger space, enhance DEOCC!'
END DO


IF(itest>0) WRITE(7,'(a)') ' SSETUP finished'
DEALLOCATE(wrk1,wrk2)

RETURN
END SUBROUTINE ssetup

!---------------------------------------------------------------------------  
! DESCRIPTION: zwf
!> @brief
!! zwf calculates the normalized one-dimensional wave functions phi(n,i). 
!> @details
!!n denotes the number of nodes along the z-axis. i labels the z-coordinate.
!!The wavefunctions are normalized such that the integral from minus infinity 
!!to infinity of dz phi**2 = 1.
!> @param[in] bz
!> REAL(8), takes the harmonic oscillator size constant.
!> @param[in] i
!> INTEGER, takes the z coordinate.
!> @param[out] phi
!> REAL(8), array, returns the normalized wave function.
!--------------------------------------------------------------------------- 

SUBROUTINE zwf(i,bz,phi)

USE paramr
USE dynpar
USE grid
IMPLICIT NONE

INTEGER, INTENT(IN)        :: i
REAL(8), INTENT(IN)        :: bz
REAL(8), INTENT(OUT)       :: phi(nnzxd+1,0:nz)

REAL(8) :: earg,expf,root2,rootpi,x,zsym
INTEGER :: nnz


root2=SQRT(two)
rootpi=SQRT(pi)

zsym = i*dz
x=zsym/bz
earg = x*x/two
IF (earg <= 60.0d0) THEN
  expf = EXP(-earg)
ELSE
  expf = zero
END IF
phi(1,i)=expf/SQRT(rootpi*bz)
phi(2,i)=root2*x*phi(1,i)
IF(nnzxd < 2) RETURN
DO  nnz=2,nnzxd
  phi(nnz+1,i)= SQRT(two/REAL(nnz))*x*phi(nnz,i) -SQRT(REAL(nnz-1)/REAL(nnz))*phi(nnz-1,i)
END DO
RETURN
END SUBROUTINE zwf

!---------------------------------------------------------------------------  
! DESCRIPTION: radwf
!> @brief
!! radwf calculates the normalized radial wave functions phi(n,l,i). 
!> @details
!!n denotes the number of nodes along the z-axis. i labels the z-coordinate.
!!The wavefunctions are normalized such that the integral from minus infinity 
!!to infinity of dz phi**2 = 1.
!> @param[in] br
!> REAL(8), takes the harmonic oscillator size constant.
!> @param[in] i
!> INTEGER, takes the r coordinate.
!> @param[in] nnr
!> INTEGER, takes the number of radial nodes.
!> @param[in] mmr
!> INTEGER, takes the magnitude of the z-component of angular momentum.
!> @param[in] psirad
!> REAL(8), array, the wave functions are normalized such that 
!> \f$ \int_{-\infty}^\infty {\rm d}r\, r \psi^2 = 1 \f$.
!--------------------------------------------------------------------------- 

SUBROUTINE radwf(i,br,nnr,mmr,psirad)
USE paramr
USE dynpar
USE grid
IMPLICIT NONE

INTEGER, INTENT(IN)               :: i
REAL(8), INTENT(IN)               :: br
INTEGER, INTENT(IN)               :: nnr
INTEGER, INTENT(IN)               :: mmr
REAL(8), INTENT(OUT)              :: psirad

REAL(8)      ::  psi(nnrxd+1)
REAL(8),SAVE :: gamma(mxd+1)

REAL(8) :: earg,el,enr,expf,rooty,ryl,xxx,y
INTEGER :: l,ll,n,nnrr,nrt
INTEGER :: ifirst=-1367

!     initialize the factorials

IF(ifirst == -1367) THEN
  IF(mxd > 33.AND.itest>1) WRITE(7,'(a,i5,a)')  &
      ' mxd=',mxd,' probably too large for single precision'
  gamma(1)=one
  DO n=1,mxd
    xxx=(n)
    gamma(n+1)=xxx*gamma(n)
  END DO
  ifirst = 0
END IF

rooty=r(i)/br
IF (rooty < one*1.e-4) THEN
  y = zero
ELSE
  y = rooty**2
END IF
earg=y*oh
IF (earg <= 60.0*one) THEN
  expf = EXP(-earg)
ELSE
  expf = zero
END IF
DO  ll=1,mmr
  l=ll-1
  el=l
  IF (l == 0) THEN
    ryl = one
  ELSE IF (rooty < 1.d-4) THEN
    ryl = zero
  ELSE
    ryl = rooty**l
  END IF
  psi(1) = SQRT(two/gamma(l+1))/br*ryl*expf
  psi(2) = (el+one-y)/SQRT(el+one)*psi(1)
  IF(nnr >= 2) THEN
    DO  nnrr=2,nnr
      nrt=nnrr-1
      enr=nrt
      psi(nnrr+1)=((two*enr+el+one-y)*psi(nnrr)  &
                 -SQRT((enr+el)*enr)*psi(nrt))/SQRT((enr+el+one)*(enr+one))
    END DO
  END IF
END DO
psirad = psi(nnr)

RETURN
END SUBROUTINE radwf

!---------------------------------------------------------------------------  
! DESCRIPTION: comp_block
!> @brief
!! Computes the blocked state after some pre-iterations. Blocked is a state 
!! with \c kblock, \c mblock and \c isoblock as given in m_levels.f90. The 
!! routine finds out what number this state has in the normal ordering 
!! sequence of the code.
!> @details
!!
!--------------------------------------------------------------------------- 
SUBROUTINE comp_block()
USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT NONE

REAL(8) :: signup
INTEGER :: ist,iunit,mbact,nfin,ninit

IF(kblock > 2*mblock) THEN
  signup = 1.0D0
  mbact  = mblock
ELSE
  signup = -1.0D0
  mbact  = mblock-1
END IF
IF(isoblo == 1) THEN
  ninit = 1
  nfin  = nst(1)
ELSE IF(isoblo == 2) THEN
  ninit = nst(1)+1
  nfin  = nst(1)+nst(2)
ELSE
  STOP ' invalid isospin in blocking routine COMP_BLOCK'
END IF

DO ist=ninit,nfin
  IF((mpsi(ist) == mbact) .AND. &
     (signup*(cmpn(1,ist)-cmpn(2,ist)) >= 0.0D0) .AND. &
     (vsqf(ist) < 0.8D0)) THEN
    istblo = ist
    EXIT
  END IF
END DO

DO iunit=6,7
  WRITE(iunit,'(a,4i5)') ' blocking found: kblock,mblock,isoblo,istblo=',  &
      kblock,mblock,isoblo,istblo
END DO

RETURN
END SUBROUTINE comp_block

END MODULE init
