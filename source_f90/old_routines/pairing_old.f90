! file 'pairing_old.f90': 
!  collects the routines for old-fashioned pairing models
!  (constant gap, constant force, LN thereof).
!
!---------------------------------------------------------------------------  
! DESCRIPTION: pairol
!> @brief
!!Determines the fermi-energy \c eferm and the gap \c delta such that the 
!!"gap-equation" and the "particle number condition" are fulfilled 
!!simultanously. this is achieved by iteration with a gradient scheme 
!!(Newton's tangential formula) where the gradient is evaluated by finite 
!!differences.
!> @detail
!!There are two schemes possible:
!!  -# The "constant force" approach, where the pairing force is given
!!     in 'gp'.
!!  -# The "constant gap" approach, where the gap is kept fixed and is given 
!!     in \c delta.
!!  .
!!There is only \c eferm to be determined in the "constant gap" approach. 
!!The gap-equation is used at the end to determine an according force \c gp.
!>
!> @param[in] e
!> REAL(8), array, takes single particle energies.
!> @param[inout] gw
!> REAL(8), array, takes occupation probabilities, input is used as initial 
!! guess for the iteration. Returns the resulting occupation probabilities.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[in] nz
!> INTEGER, takes the desired particle number.
!> @param[in] nmax
!> INTEGER, takes maximum number of single particle states
!! (needs to be smaller than the parameter 'npsid'
!! and smaller than 46, for dimensioning reasons).
!> @param[inout] gp
!> REAL(8), takes size of pairing-force, is input in case of
!! "constant force" approach, output in case of "constat gap".
!> @param[inout] eferm
!> REAL(8), returns Fermi energy, input is used as initial guess; if eferm=0.0 
!! is input then a suitable initial guess is computed from input \c gw.
!> @param[inout] delta 
!> REAL(8), takes pairing gap; is input in case of "constant gap" approach,
!! output in case of "constant force".
!> @param[inout] partnm
!> REAL(8), returns the particle number, from summing \c gw with degeneracy 
!! \c ph.
!> @param[in] iter
!> REAL(8), takes number of pairing iterations, abs(iter) is the number of 
!! iterations.
!> @param[in] ipair
!> INTEGER, takes pairing type (ipair=1: "constant force" approach, using
!! pairing force 'gp'; ipair=0: "constant gap" approach, using given gap 
!! \c delta; ipair=4: thermal occupation with temperature given by \c delta.
!> @param[in] eps
!> INTEGER, takes end condition for pairing iteration; iteration terminates 
!! if error in particle number is smaller than \c eps and if \c eferm and/or 
!! \c delta are smaller than \c eps.
!> @param[in] iab
!> INTEGER, takes print level; significant values are -2, -1, 0 and 1;
!! larger value prints more.
!--------------------------------------------------------------------------- 

SUBROUTINE pairol(e,gw,ph,nz,nmax,gp,eferm,delta,partnm,iter,ipair,eps,iab)

USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)        :: e(npsid)
REAL(8), INTENT(INOUT)     :: gw(npsid)
REAL(8), INTENT(IN)        :: ph(npsid)
INTEGER, INTENT(IN)        :: nz
INTEGER, INTENT(IN)        :: nmax
REAL(8), INTENT(INOUT)     :: gp
REAL(8), INTENT(INOUT)     :: eferm
REAL(8), INTENT(INOUT)     :: delta
REAL(8), INTENT(OUT)       :: partnm
INTEGER, INTENT(INOUT)     :: iter
INTEGER, INTENT(IN)        :: ipair
REAL(8), INTENT(IN)        :: eps
INTEGER, INTENT(IN)        :: iab


REAL(8) :: parnm

REAL(8) :: dparn(3)=(/0D0,0D0,0D0/) !array(3) of variation in the particle number as
!                                    consequence of 'deferm' and 'ddelta'.
REAL(8) :: dgap(3)=(/0D0,0D0,0D0/)  !array(3) of variation in the gap equation as
!                                    consequence of 'deferm' and 'ddelta'.
REAL(8) :: gr(2,2)=0.0d0            !double array of variation of particle number and of
!                                    gap equation, used for construction of the tangantial step.

!      efmax,dlmax,detmin,defmin,ddlmin  = limitations for the step
REAL(8),PARAMETER :: dampin=0.9D0  !damping factor for the tangential step
REAL(8),PARAMETER :: xmaxlw=1D-30,x10lw=1D-10
REAL(8),SAVE      :: deferm=1D0    !step in the fermi-energy, to evaluate
!                                   numerically the derivatives.
REAL(8),SAVE      :: ddelta=0.1D0  !step in the gap, to evaluate numerically the derivatives;
!                                   used only in the "constant force" approach.
REAL(8),PARAMETER :: efmax=3D0,dlmax=0.6D0,detmin=1D-4
REAL(8),PARAMETER :: defmin=1D-4,ddlmin=1D-4

INTEGER :: k,iwarng,itsum,itmaxp,it,isum,i
REAL(8) :: parn1,parn2,parn3,parnum,grsav,gphalf,gapeq,equasi,emin,emax,elam,&
           efstep,efold,eferm1,eferm2,eferm3,dlstep,det,delti,delold,delact,del2
LOGICAL :: converged

!-------------------------------------------------------------------------------
!      set switch variables
iwarng = 0

IF(iab >= 1) WRITE(6,'(/a,i4,3(a,f9.4)/a,10(/1x,10f8.3))')  &
    ' *** pair entered with: iter=',iter,'  delta=',delta,  &
    '  eferm=',eferm,'  gp=',gp,  &
    '                        single particle energies:', (e(i),i=1,nmax)

!      prepare fermi-energy 'eferm', if not yet initialised.
IF(eferm == zero) THEN
  isum   = 0
  DO  it=1,nmax
    itsum  = it
    isum   = INT(ph(it)+.000001) +isum
    IF(2.0d0*isum >= nz) EXIT
  END DO
  IF(2.0d0*isum < nz) STOP ' not enough states too reach particle number'
  eferm  = e(itsum)+deferm
  IF(iab >= 1) WRITE(06,'(1x,a,i5,(a,g11.4),a,i5)')  &
      'initialistion: isum=',isum,'  eferm=',eferm,'  nz=',nz
END IF

IF(ipair == 1) THEN
!-------------------------------------------------------------------------------
! case of constant pairing-force, given by 'gp' gap 'delta', fermi-energy 
! 'eferm', and occupations are computed.
!-------------------------------------------------------------------------------
  delti  = delta
  IF(ABS(delti) < x10lw) delti = one
  
  DO  it=1,ABS(iter)
    
! gap-equation & particle-number at 'eferm' & 'delti' (for k=1)
! and at 'eferm' varied (k=2) and at 'delta' varied (k=3).
    DO  k=1,3,-1  !TODO: There is sth wrong here!!
! actual fermi-energy and gap
      WRITE(*,*)k
      IF(k == 1) THEN
        elam   = eferm
        delact = delti
      ELSE IF(k == 2) THEN
        elam   = eferm+deferm
        delact = delti
      ELSE
        elam   = eferm
        delact = delti+ddelta
      END IF
! compute actual quasiparticle energies and occupation-weights, accumulate 
! gap-equation and particle number
      parnum = zero
      gapeq  = zero
      DO  i=1,nmax
        equasi = SQRT((e(i)-elam)*(e(i)-elam)+delact*delact)
        gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
        gapeq  = ph(i)/equasi   + gapeq
        parnum = gw(i)*ph(i)    + parnum
      END DO
! store actual part.number and gap-equation
      dgap(k)  = 0.25*gp*gapeq-one          !?? extra 1/2 because of WATE=2
      dparn(k) = parnum-nz
    END DO
    
! construction of the gradient-matrix and its inverse
! ( determinant 'det' is limited in value to avoid singularities)
    gr(1,1) = ( dgap(2)- dgap(1))/deferm
    gr(1,2) = (dparn(2)-dparn(1))/deferm
    gr(2,1) = ( dgap(3)- dgap(1))/ddelta
    gr(2,2) = (dparn(3)-dparn(1))/ddelta
    det     = gr(1,1)*gr(2,2)-gr(1,2)*gr(2,1)
    det     = SIGN(MAX(ABS(det),detmin),det)
    IF(ABS(det) < xmaxlw)  det    = SIGN(xmaxlw,det)
    
    IF(iab >= 1) THEN
      WRITE(6,'(a,i4,3(a,g11.3)/a,3g11.3,a,3g11.3/a,4g11.3)')  &
          ' iteration at it=',it,' :  delti=',delti,'  eferm=',eferm,  &
          '  det=',det, ' dgap=',dgap,' ,  dparn=',dparn,  &
          '   gr=',gr
      WRITE(6,'(a,10f7.4,10(/6x,10f7.4))') '   gw=',gw
    END IF
    
    grsav   =  gr(1,1)
    gr(1,1) =  gr(2,2)/det
    gr(1,2) = -gr(1,2)/det
    gr(2,1) = -gr(2,1)/det
    gr(2,2) = grsav/det
    
!       tangential step  (including some step limiters)
    efold  = eferm
    delold = delti
    efstep = dampin*(gr(1,1)*dgap(1)+gr(2,1)*dparn(1))
    dlstep = dampin*(gr(1,2)*dgap(1)+gr(2,2)*dparn(1))
    eferm  = eferm-SIGN(MIN(ABS(efstep),efmax),efstep)
    delti  = ABS(delti-SIGN(MIN(ABS(dlstep),dlmax),dlstep) )
    deferm = efold-eferm
    ddelta = delold-delti
    IF(ABS(deferm) < defmin) deferm = defmin
    IF(ABS(ddelta) < ddlmin) ddelta = ddlmin
    
    IF(iab == 0) WRITE(6,'(a,i4,2(a,g11.3))')  &
        ' pairing iteration it=',it,'  :  delti=',delti, '  eferm=',eferm
    IF(ABS(dgap(1)) < eps .AND. ABS(dparn(1)) < eps)   EXIT
  END DO
  IF(ABS(dgap(1)) > eps .OR. ABS(dparn(1)) > eps) iwarng = 1

  delta  = delti  
ELSE
  
!-------------------------------------------------------------------------------
!      case of constant paring-gap, given by 'delta'.
!      using now the safer secant step.
!-------------------------------------------------------------------------------
  
! determine upper and lower bound for secant step
  itmaxp = ABS(iter)
  eferm1 = eferm
  eferm2 = eferm
  parn1  = parnm(e,gw,ph,nmax,delta,eferm,ipair)
  parn2  = parn1
  converged=.FALSE.
  IF(parn1-nz > zero) THEN
    emin = +1.0E30
    DO  i=1,nmax
      emin = MIN(e(i),emin)
    END DO
    deferm = (eferm-emin)/itmaxp
    IF(iab >= 1) WRITE(6,'(a,3(1pg12.4))')  &
        ' search brackets: parn2,eferm,emin=',parn2,eferm,emin
    DO  it=1,itmaxp
      eferm1 = eferm1-deferm
      parn1  = parnm(e,gw,ph,nmax,delta,eferm1,ipair)
      IF(iab >= 1) WRITE(6,'(a,i5,3(1pg12.4))')  &
          ' search brackets: it,parn1,eferm1=',it,parn1,eferm1
      IF(parn1-nz < zero) GO TO 319
      eferm2 = eferm1
    END DO
    STOP ' particle number not embraced in PAIR'
  ELSE IF(parn1-nz < zero) THEN
    emax = -1.0E30
    DO  i=1,nmax
      emax = MAX(e(i),emax)
    END DO
    deferm = (emax-eferm)/itmaxp
    IF(iab >= 1) WRITE(6,'(a,3(1pg12.4))')  &
        ' search brackets: parn1,eferm,emax=',parn1,eferm,emax
    DO  it=1,itmaxp
      eferm2 = eferm2+deferm
      parn2  = parnm(e,gw,ph,nmax,delta,eferm2,ipair)
      IF(iab >= 1) WRITE(6,'(a,i5,3(1pg12.4))')  &
          ' search brackets: it,parn2,eferm2=',it,parn2,eferm2
      IF(parn2-nz > zero) GO TO 319
      eferm1 = eferm2
    END DO
    STOP ' particle number not embraced in PAIR'
  ELSE
    IF(iab >= 1) WRITE(6,'(a)') ' solution immediately found '
    converged=.TRUE.            ! solution immediately found. jump to end
  END IF
  319  CONTINUE             !  particle number 'NZ' between 'eferm1 <--> 2'
  IF(iab >= 1) WRITE(6,'(a)') ' brackets found '
  
!      now perform the secant search
  IF(.NOT.converged) THEN
  DO it=1,itmaxp
    eferm3 = eferm1 + (nz-parn1)*(eferm2-eferm1)/(parn2-parn1)
    parn3  = parnm(e,gw,ph,nmax,delta,eferm3,ipair)
    IF(iab >= 2) WRITE(6,'(a,i4,2(a,g11.3))')  &
        ' secant iteration it=',it,'  :  eferm=',eferm3, '  dparn(1)=',parn3-nz
    IF(ABS(parn3-nz) < eps) THEN
      converged=.TRUE.
      EXIT
    ELSE IF(parn3-nz > zero) THEN
      eferm2 = eferm3
      parn2  = parn3
    ELSE
      eferm1 = eferm3
      parn1  = parn3
    END IF
  END DO
  END IF
  
!      the secant step did not succeed. use bisection
  IF(.NOT.converged) THEN
    DO it=1,itmaxp
      eferm3 = 0.5*(eferm2+eferm1)
      parn3  = parnm(e,gw,ph,nmax,delta,eferm3,ipair)
      IF(iab >= 0) WRITE(6,'(a,i4,2(a,g11.3))')  &
          ' bisection it=',it,'  :  eferm=',eferm3, '  dparn(1)=',parn3-nz
      IF(ABS(parn3-nz) < eps) THEN
      converged=.TRUE.
      EXIT
      ELSE IF(parn3-nz < zero) THEN
        eferm1 = eferm3
        parn1  = parn3
      ELSE IF(parn3-nz > zero) THEN
        eferm2 = eferm3
        parn2  = parn3
      END IF
    END DO
  END IF
  IF(.NOT.converged) iwarng = 1  ! still not converged. set a flag
  
  eferm = eferm3
  dparn(1) = parn3 - nz
  dgap(1)  = zero
  partnm   = parn3
  
! determine force 'gp' according to given 'delta' and the equilibrium reached.
  IF(ipair == 0) THEN
    gapeq  = zero
    DO i=1,nmax
      gapeq = ph(i)/SQRT((e(i)-eferm)*(e(i)-eferm)+delta*delta) + gapeq
    END DO
    gp     = 4.0/gapeq                       !?? account for WATE=2.0
  END IF
  
!     end of big switch between cases
END IF

!-------------------------------------------------------------------------------
!     finally the gap-equation and particle-number are checked and printed
!-------------------------------------------------------------------------------

!     in case of incomplete convergence print a warning
IF(iab >= -1 .AND. iwarng == 1) WRITE(6,'(2(a,g11.3),2(a,i3))')  &
    ' ---> no convergence in pair! err(gap-eq)=',dgap(1),  &
    '  err(part.nr)=',dparn(1),'  nz=',nz,'  nmzx=',nmax

!     check gap equation
IF(iab >= 0) THEN
  IF(ipair /= 4) THEN
    gphalf = gp*0.5
    del2   = MAX(xmaxlw,delta*delta)
    gapeq  = zero
    partnm = zero
    DO i=1,nmax
      partnm = gw(i)*ph(i) + partnm
      gapeq  = gphalf*ph(i)/SQRT(del2+(e(i)-eferm)**2) + gapeq
    END DO
  ELSE
    gapeq = zero
  END IF
  
  WRITE(6,'(a,i4,2(a,g13.5))')  &
      ' pair finished with ',it,' iterations: gap-eq.=',gapeq, '  part.nr=',partnm
END IF

RETURN
END SUBROUTINE pairol

!---------------------------------------------------------------------------  
! DESCRIPTION: parnm
!> @brief
!!The function computes the particle number in a pairing or temperature 
!!distribution. 
!>
!> @param[in] e
!> REAL(8), array, takes single particle energies.
!> @param[inout] gw
!> REAL(8), array, takes occupation probabilities, input is used as initial 
!! guess for the iteration. Returns the resulting occupation probabilities.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[in] nmax
!> INTEGER, takes number of states in \c e, \c gw and \c ph.
!> @param[in] delta
!> REAL(8), takes pairing gap or temperature.
!> @param[in] elam
!> REAL(8), takes the Fermi energy.
!> @param[in] ipair
!> INTEGER, switch to pairing (0) or temperature (4).
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION parnm(e,gw,ph,nmax,delta,elam,ipair)

USE paramr

IMPLICIT NONE

REAL(8), INTENT(IN)                     :: e(npsid)
REAL(8), INTENT(OUT)                    :: gw(npsid)
REAL(8), INTENT(IN)                     :: ph(npsid)
INTEGER, INTENT(IN)                     :: nmax
REAL(8), INTENT(IN)                     :: delta
REAL(8), INTENT(IN)                     :: elam
INTEGER, INTENT(IN)                     :: ipair

INTEGER :: i
REAL(8) :: equasi,sum

sum    = zero
IF(ipair == 0) THEN
  DO  i=1,nmax
    equasi = SQRT( (e(i)-elam)*(e(i)-elam) + delta*delta )
    gw(i)  = 0.5-0.5*(e(i)-elam)/equasi
    sum    = gw(i)*ph(i)    + sum
  END DO
ELSE IF(ipair == 4) THEN            ! case of temperature
  DO  i=1,nmax
    equasi = (e(i)-elam)/delta
    IF(equasi > 50.0d0) THEN
      gw(i) = 0.0
    ELSE IF(equasi < -50.0d0) THEN
      gw(i)  = 1.0
    ELSE
      gw(i)  = 1.0/(1.0+EXP(equasi))
    END IF
    sum    = gw(i)*ph(i)    + sum
  END DO
END IF
parnm = sum

RETURN
END FUNCTION parnm
!---------------------------------------------------------------------------  
! DESCRIPTION: bcslip
!> @brief
!!Routine to solve the Lipkin-Nogami or BCS equations.
!> @detail
!!  - gap equation: \f[ 0 = \texttt{f1} = \texttt{delta}-\texttt{gpair}\cdot \sum(uv) \f]
!!  - particle number: \f[ 0 = \texttt{f2} = \sum(2v^2)-n \f]
!!  - particle number fluctuation:  \f[ 0 = \texttt{f3} = \texttt{xlam2}-\texttt{gpair}/4 \cdot 
!!                                      (\sum(u^3\cdot v)\sum(v^3* \cdot u)-\sum(u^4*\cdot v^4))/
!!                                      ((\sum(u^2\cdot v^2))^2 - \sum(u^4 \cdot v^4)) \f]
!! Particle number fluctiation is ommited for BCS. The 3 equations are 
!!solved by Newton's method. The matrix of derivatives is computed by 
!!finite differences. 
!>
!> @param[in] spen
!> REAL(8), array, takes single particle energies.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[in] nspd
!> INTEGER, takes lowest level in pairing.
!> @param[in] nspg
!> INTEGER, takes highest level in pairing.
!> @param[in] npart
!> INTEGER, takes the particle number.
!> @param[in] gpair
!> REAL(8), takes the pairing strength.
!> @param[inout] delta
!> REAL(8), takes first guess of pairing gap and returns final pairing gap.
!> @param[inout] xlam
!> REAL(8), takes first guess for Fermi energy and returns final Fermi energy.
!> @param[inout] xlam2
!> REAL(8), takes the initial value for fluctuation constraint and returns 
!! final value.
!> @param[inout] vv2
!> REAL(8), array, takes guess for internal occupation probabilities and 
!! returns final values.
!> @param[out] vln
!> REAL(8), array, returns occupation probabilities for external observables.
!> @param[out] epair
!> REAL(8), returns naive pairing energy , \f$ - \texttt{delta}^2/\texttt{gpair} \f$.
!> @param[out] ecorr
!> REAL(8), returns consistent pairing energy as correction for the hf energy:
!! \f[-\texttt{delta}^2/\texttt{gpair}-\texttt{gpair}\sum{v^4-v^2}-4\cdot \texttt{xlam2}
!! \sum_(uv)^2= -(4*\texttt{xlam2}-\texttt{gpair})\sum_(uv)^2 \f]
!> @param[out] entot
!> REAL(8), returns total s.q.p. energy (= s.p. sum - g * \c npart/2 for normal state).
!> @param[out] erspln
!> REAL(8), takes Skyrme response energy for 'lambda_2', needed for
!! Skyrme H.F. codes, to be provided in a separte routine. In all other 
!! cases set to 0.0 .
!> @param[out] suv2rs
!> REAL(8), \f$ u^2\cdot v^2 \f$ from Skyrme response routine, for 
!!normalization. In all other cases set to one.
!> @param[in] ilipng
!> INTEGER, switch Lipkin-Nogami (=1) and BCS (=0).
!> @param[in] iftest
!> INTEGER, switch for testprints.
!--------------------------------------------------------------------------- 
SUBROUTINE bcslip(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,  &
    vv2,vln,epair,ecorr,entot,erspln,suv2rs, ilipng,iftest)

USE paramr
IMPLICIT NONE

REAL(8), INTENT(IN)           :: spen(npsid)
REAL(8), INTENT(IN)           :: ph(npsid)
INTEGER, INTENT(IN)           :: nspd
INTEGER, INTENT(IN)           :: nspg
INTEGER, INTENT(IN)           :: npart
REAL(8), INTENT(IN)           :: gpair
REAL(8), INTENT(INOUT)        :: delta
REAL(8), INTENT(INOUT)        :: xlam
REAL(8), INTENT(INOUT)        :: xlam2
REAL(8), INTENT(INOUT)        :: vv2(npsid)
REAL(8), INTENT(OUT)          :: vln(npsid)
REAL(8), INTENT(OUT)          :: epair
REAL(8), INTENT(OUT)          :: ecorr
REAL(8), INTENT(OUT)          :: entot
REAL(8), INTENT(OUT)          :: erspln
REAL(8), INTENT(OUT)          :: suv2rs
INTEGER, INTENT(IN)           :: ilipng
INTEGER, INTENT(IN)           :: iftest

REAL(8),PARAMETER :: epsf1=1D-8    !precision on pairing equation 'f1'.
REAL(8),PARAMETER :: epsf2=1D-8    !precision on pairing equation 'f2'.
REAL(8),PARAMETER :: epsf3=1D-8    !precision on pairing equation 'f3'.
REAL(8),PARAMETER :: epsv2=1D-8    !precision on 'vv2'.
REAL(8),PARAMETER :: epsdv2=1D-8   !limit on change of energy in determining 'vv2'.
REAL(8),PARAMETER :: detmin=1D-1   !limiting minimum value for the determinant in linear eq.
REAL(8),PARAMETER :: ddllim=0.2D0  !limit on stepsize in 'delta'.
REAL(8),PARAMETER :: dxllim=3D0    !limit on stepsize in 'xlam'.
REAL(8),PARAMETER :: dx2lim=1D0    !limit on stepsize in 'xlam2'.
REAL(8),PARAMETER :: dampin=0.8D0  !general reduction of the step, top stabilize.
REAL(8),PARAMETER :: delstp=1D-2   !step in 'delta' for first derivative in finite difference.
REAL(8),PARAMETER :: xlstp=1D-2    !step in 'xlam' for first derivative in finite difference.
REAL(8),PARAMETER :: xl2stp=1D-2   !step in 'xlam2' for first derivative in finite difference.
INTEGER,PARAMETER :: itbcsl=50     !max.nr. iterations in 3 pairing parameters.
INTEGER,PARAMETER :: iterv2=50     !max.nr. iterations in 'vv2'.

INTEGER :: itbcs,k
REAL(8) :: b11,b12,b13,b21,b22,b23,b31,b32,b33,c,ddelac,ddelta,ddf1,ddf2,&
           ddf3,det,detinv,dl2f1,dl2f2,dl2f3,dxlam,dxlam2,efluc,epair1,f1,&
           dlf1,dlf2,dlf3,f1del,f1xl,f1xl2,f2,f2del,f2xl,f2xl2,f3,&
           f3del,f3xl,f3xl2,od,s,sspe,sum2,sum4,sumc,sumc2,suv,suv2,sv2,sv4,&
           uv,v2,xmult
LOGICAL :: converged
!---------------------------------------------------------------------

! preparation and tests
IF(nspg > npsid) STOP ' last level > "KMAXLN" in BCSLIP'
IF(ilipng /= 1) THEN
  xlam2 = zero
  dxlam2 = zero
END IF
converged=.FALSE.
! iteration of the three pairing parameters 'delta', 'xlam', and 'xlam2' .
DO  itbcs=1,itbcsl
  IF(iftest > 0) WRITE(6,'(a,i5,a)')  ' iteration nr.',itbcs,' :'
  
! initial call to pairing equations
  CALL eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,  &
      vv2,vv2,erspln,suv2rs,iftest,ilipng, iterv2, epsv2,epsdv2,dampin,  &
      f1,f2,f3)
  IF(iftest > 0) THEN
    WRITE(6,'(a,3g12.3/a,3g12.3)') ' delta,xlam,xlam2=',delta,xlam,xlam2,  &
        '         f1,f2,f3=',f1,f2,f3
  END IF
  
! test termination
  IF(itbcs > 1 .AND. ABS(f1) < epsf1 .AND. ABS(f2) < epsf2  &
                                     .AND. ABS(f3) < epsf3) THEN
    converged=.TRUE.
    EXIT
  END IF
  
! derivatives by finite differences
  CALL eqdef(spen,ph,nspd,nspg,npart,gpair,delta*(one+delstp),  &
      xlam,xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng, &
      iterv2, epsv2,epsdv2,dampin, f1del,f2del,f3del)  ! vln = aux 
  CALL eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam+xlstp,  &
      xlam2,vv2,vln,erspln,suv2rs,iftest,ilipng,   &
      iterv2, epsv2,epsdv2,dampin, f1xl,f2xl,f3xl)  ! vln = aux 
  IF(ilipng == 1)  &
      CALL eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2+xl2stp,  &
      vv2,vln,erspln,suv2rs,iftest,ilipng,             &
      iterv2, epsv2,epsdv2,dampin, f1xl2,f2xl2,f3xl2)   ! vln = aux 
  ddelac = delta*delstp
  ddf1   = (f1del-f1)/ddelac
  ddf2   = (f2del-f2)/ddelac
  dlf1   = (f1xl-f1)/xlstp
  dlf2   = (f2xl-f2)/xlstp
  IF(ilipng == 1) THEN
    dl2f1  = (f1xl2-f1)/xl2stp
    dl2f2  = (f2xl2-f2)/xl2stp
    ddf3   = (f3del-f3)/ddelac
    dlf3   = (f3xl-f3)/xlstp
    dl2f3  = (f3xl2-f3)/xl2stp
    
!    matrix of inverse is computed explicitely and stored on bij
    det    =(ddf1*dlf2*dl2f3-ddf1*dl2f2*dlf3-ddf2*dlf1*dl2f3  &
        +ddf2*dl2f1*dlf3+ddf3*dlf1*dl2f2-ddf3*dl2f1*dlf2)
    det    = SIGN(MAX(ABS(det),detmin),det)        ! lower limit for 'det'
    detinv = one/det
    b11 = (dlf2*dl2f3-dl2f2*dlf3)*detinv
    b12 = -(dlf1*dl2f3-dl2f1*dlf3)*detinv
    b13 = (dlf1*dl2f2-dl2f1*dlf2)*detinv
    b21 = -(ddf2*dl2f3-dl2f2*ddf3)*detinv
    b22 = (ddf1*dl2f3-dl2f1*ddf3)*detinv
    b23 = -(ddf1*dl2f2-dl2f1*ddf2)*detinv
    b31 = (ddf2*dlf3-dlf2*ddf3)*detinv
    b32 = -(ddf1*dlf3-dlf1*ddf3)*detinv
    b33 = (ddf1*dlf2-dlf1*ddf2)*detinv
    
    ddelta = -b11*f1-b12*f2-b13*f3
    dxlam  = -b21*f1-b22*f2-b23*f3
    dxlam2 = -b31*f1-b32*f2-b33*f3
    IF(iftest > 0) WRITE(6,'(a,4g12.3)')  &
        ' det,ddelta,dxlam,dxlam2=',det,ddelta,dxlam,dxlam2
  ELSE
    
!     case of BCS, step in two quantities
    det    = -(dlf1*ddf2-dlf2*ddf1)
    det    = SIGN(MAX(ABS(det),detmin),det)      ! lower limit for 'det'
    ddelta = (f2*dlf1-f1*dlf2)/det
    dxlam  = (f1*ddf2-f2*ddf1)/det
  END IF
  
!     limit stepsize if necessary
  ddelta = dampin*SIGN(MIN(ABS(ddelta/delta),ddllim),ddelta)
  dxlam  = dampin*SIGN(MIN(ABS(dxlam ),dxllim),dxlam )
  IF(ilipng == 1) dxlam2 = dampin*SIGN(MIN(ABS(dxlam2),dx2lim),dxlam2)
  
!     step forward
  delta  = delta*EXP(ddelta)
  xlam   = xlam  + dxlam
  IF(ilipng == 1) xlam2  = xlam2 + dxlam2
END DO

IF(iftest >= 0.AND..NOT.converged) WRITE(6,'(a)')  &
    ' no convergence of pairing parameters in BCSLIP'

IF(iftest >= 0) WRITE(6,'(a,i4,6(1pg11.3))')  &
    ' BCSLIP: iter,delta,xlam,xlam2,f1,f2,f3=', itbcs,delta,xlam,xlam2,f1,f2,f3

!     compute some observables after convergence

suv    = zero
sv2    = zero
sv4    = zero
suv2   = zero
sspe   = zero
DO  k=nspd,nspg
  v2     = vv2(k)
  uv     = SQRT(v2-v2*v2)
  vln(k) = uv                                     ! auxiliary storage
  xmult  = ph(k)
  suv    = suv    + xmult*uv
  sv2    = sv2    + xmult*v2
  sv4    = sv4    + xmult*v2*v2
  suv2   = suv2   + xmult*(v2 - v2*v2)
  sspe   = sspe   + xmult*spen(k)*v2
END DO
ecorr  = -delta*delta/gpair-(4.0*xlam2-gpair)*suv2
entot  =  sspe - gpair * npart/2
epair1 = -gpair*sv2
efluc  = -4.*xlam2*suv2
epair  = -gpair*suv
epair1 = epair1+efluc+epair

!           corrected coefficients, vln(i)=<a^+_i a_i>_ln
IF(ilipng == 1) THEN
  sum2=0.
  sum4=0.
  sumc=0.
  sumc2=0.
  
  DO  k=nspd,nspg
    v2    = vv2(k)
    xmult = ph(k)
    uv    = vln(k)
    s     = uv+uv
    c     = v2+v2-one
    sum2  = sum2  + xmult*s**2
    sum4  = sum4  + xmult*s**4
    sumc  = sumc  + xmult*s*s*c
    sumc2 = sumc2 + xmult*(s*c)**2
  END DO
  
  od=sum2*0.25/(sumc2*sum2-sumc*sumc+0.25*sum2*(sum2*sum2-sum4))
  
  DO  k=nspd,nspg
    v2    = vv2(k)
    uv    = vln(k)
    s     = uv+uv
    c     = v2+v2-one
    vln(k)= v2+od*s*s*(c*sum2-sumc)
  END DO
ELSE
  DO  k=nspd,nspg
    vln(k)= vv2(k)
  END DO
END IF

RETURN
END SUBROUTINE bcslip

!---------------------------------------------------------------------------  
! DESCRIPTION: eqdef
!> @brief
!!Routine to compute the three pairing equations, as explained in the 
!!header of \c bcslip.
!>
!> @param[in] spen
!> REAL(8), array, takes single particle energies.
!> @param[in] ph
!> REAL(8), array, takes the multiplicities of the single particle levels.
!> @param[in] nspd
!> INTEGER, takes lowest level in pairing.
!> @param[in] nspg
!> INTEGER, takes highest level in pairing.
!> @param[in] npart
!> INTEGER, takes the particle number.
!> @param[in] gpair
!> REAL(8), takes the pairing strength.
!> @param[inout] delta
!> REAL(8), takes first guess of pairing gap and returns final pairing gap.
!> @param[inout] xlam
!> REAL(8), takes first guess for Fermi energy and returns final Fermi energy.
!> @param[inout] xlam2
!> REAL(8), takes the initial value for fluctuation constraint and returns 
!! final value.
!> @param[in] vv2in
!> REAL(8), array, takes the initial internal occupation probabilities.
!> @param[out] vv2out
!> REAL(8), array, returns the new occupation weights.
!> @param[out] erspln
!> REAL(8), takes Skyrme response energy for 'lambda_2', needed for
!! Skyrme H.F. codes, to be provided in a separte routine. In all other 
!! cases set to 0.0 .
!> @param[out] suv2rs
!> REAL(8), \f$ u^2\cdot v^2 \f$ from Skyrme response routine, for 
!! normalization. In all other cases set to one.
!> @param[in] ilipng
!> INTEGER, switch Lipkin-Nogami (=1) and BCS (=0).
!> @param[in] iftest
!> INTEGER, takes switch for testprints.
!> @param[in] iterv2
!> INTEGER, takes max.nr. iterations in \c vv2.
!> @param[in] epsv2
!> REAL(8), takes precision on \c vv2.
!> @param[in] epsdv2
!> REAL(8), takes limit on change of energy in determining \c vv2.
!> @param[in] dampin
!> REAL(8), takes general reduction of the step, top stabilize.
!> @param[out] f1
!> REAL(8), returns gap equation.
!> @param[out] f2
!> REAL(8), returns particle number condition.
!> @param[out] f3
!> REAL(8), returns particle number fluctuation.
!--------------------------------------------------------------------------- 
SUBROUTINE eqdef(spen,ph,nspd,nspg,npart,gpair,delta,xlam,xlam2,  &
    vv2in,vv2out,erspln,suv2rs,iftest,ilipng, iterv2, epsv2,epsdv2,dampin,  &
    f1,f2,f3)

!---------------------------------------------------------------------
USE paramr
IMPLICIT NONE

REAL(8), INTENT(IN)       :: spen(npsid)
REAL(8), INTENT(IN)       :: ph(npsid)
INTEGER, INTENT(IN)       :: nspd
INTEGER, INTENT(IN)       :: nspg
INTEGER, INTENT(IN)       :: npart
REAL(8), INTENT(IN)       :: gpair
REAL(8), INTENT(IN)       :: delta
REAL(8), INTENT(IN)       :: xlam
REAL(8), INTENT(IN)       :: xlam2
REAL(8), INTENT(IN)       :: vv2in(npsid)
REAL(8), INTENT(OUT)      :: vv2out(npsid)
REAL(8), INTENT(IN)       :: erspln
REAL(8), INTENT(IN)       :: suv2rs
INTEGER, INTENT(IN)       :: iftest
INTEGER, INTENT(IN)       :: ilipng
INTEGER, INTENT(IN)       :: iterv2
REAL(8), INTENT(IN)       :: epsv2
REAL(8), INTENT(IN)       :: epsdv2
REAL(8), INTENT(IN)       :: dampin
REAL(8), INTENT(OUT)      :: f1
REAL(8), INTENT(OUT)      :: f2
REAL(8), INTENT(OUT)      :: f3

INTEGER :: it,k
REAL(8) :: alpha,alphd2,d2,denom,derf,dx,e,enume1,enume2,eq,f,fnn,&
           su3v,suv,suv2,suv3,suv4,sv2,u2,uv,v2,x,xmult

!     preparation
fnn    = npart
IF(nspg > npsid) STOP ' last level > "npsid"'
d2     = delta*delta
alpha  = 4.*xlam2-gpair
alphd2 = alpha*0.5d0*d2

!     reset accumulators for combinations of 'u' and 'v'.

sv2   = zero      ! sum v^2
suv   = zero      ! sum u*v
suv4  = zero      ! sum u^4 v^4
suv2  = zero      ! sum u^2 v^2
su3v  = zero      ! sum u*v u^2
suv3  = zero      ! sum u*v v^2

DO  k=nspd,nspg   ! each state separately
  
! iteration loop to compute occupation 'vv2'
  e     = spen(k)-xlam
  v2    = vv2in(k)      ! start from previous value
  x     = e+alpha*v2
  
  DO  it=1,iterv2
    EQ    = SQRT(d2+x*x)
    v2    = 0.5d0-0.5d0*x/EQ
    f     = e+alpha*v2-x
    derf  = -alphd2/(EQ*EQ*EQ)-one
    dx    = f/derf
    IF(ABS(f) < epsv2 .AND. ABS(dx) < epsdv2) EXIT
    x     = x-dx
    IF(it==iterv2) WRITE(6,'(a)') ' maximal number of iterations exceeded in "EQDEF"'
  END DO
  
  IF(iftest > 1) WRITE(6,'(a,i5,a,i5,g12.3)')  &
      ' state nr.',k,' :  iterv,v2',it,v2
  vv2out(k)=v2
  
! accumulate sums with new 'v2'
  xmult = ph(k)
  u2    = one-v2
  uv    = SQRT(u2*v2)
  sv2   = sv2  + xmult*v2
  suv   = suv  + xmult*uv
  IF(ilipng == 1) THEN
    suv4  = suv4 + xmult*u2*u2*v2*v2
    suv2  = suv2 + xmult*u2*v2
    su3v  = su3v + xmult*uv*u2
    suv3  = suv3 + xmult*uv*v2
  END IF
END DO

! final assembling of the three pairing equations
f1     = delta-gpair*suv                    ! gap equation
f2     = sv2+sv2-fnn                        ! particle number

! the numerator from Skyrme response is a bit special.
! it was computed with an older "v2". The corresponding
! "sum V2*U2" is stored in "suv2rs". In order to avoid
! all to frequent recomputation of that expensive piece,
! we rescale it here just by the new "suv2".
IF(ilipng == 1) THEN
  denom  = (suv2*suv2-suv4)                   !  denominator/32
  enume1 = -gpair/4.0*(su3v*suv3-suv4)        !  1. numerator/32
  enume2 = suv2*suv2*erspln/(suv2rs*suv2rs*32.0)   !  2. numerator/32
  f3     = xlam2+(enume1+enume2)/denom        ! p.n. fluctuations
END IF

RETURN
END SUBROUTINE eqdef

!---------------------------------------------------------------------------  
! DESCRIPTION: gavrg
!> @brief
!!Calculates the strength of the monopole pairing interaction using the 
!!average gap method.
!> @detail
!!the average gap model is based on moller and nix, nucl. phys.a536 (1992) 
!!20. the average pair density is computed using the semi-classical 
!!approximation  for   the  harmonic oscillator model modified  by  the  
!!isospin-dependent factor derived based on the wigner-kirkwood expansion
!!for the woods-saxon potential (+ spin-orbit + coulomb).
!>
!> @param[in] ibcs
!> INTEGER, switch for pairing mode. 1 for BCS or Lipkin-Nogami for all 
!! other values.
!> @param[in] iz
!> INTEGER, takes the proton number.
!> @param[in] in
!> INTEGER, takes the neutron number
!> @param[in] iprot
!> INTEGER, switch for Neutrons (0) or protons (1).
!> @param[in] rhosc
!> REAL(8), takes semiclassical level density at Fermi surface, determined
!! by function 'semrho' or by simple estimate if 'rhosc=0D0' is given.
!> @param[in] n1
!> INTEGER, takes number of first HF level used in the BCS (LN) 
!!calculations. (Preferably \c n1 = 1)
!> @param[in] n2
!> INTEGER, takes number of last HF level used in the BCS (LN) 
!!calculations. (Preferably \c n2 = \c npart, i.e., \c iz or \c in)
!> @param[in] xmeff
!> REAL(8), takes effective mass (=1 in the Wood-Saxon model).
!> @param[in] a0
!> REAL(8), takes surface parameter in a Saxon Woods shape, is 0.58 * 
!! surface width from formfactor
!--------------------------------------------------------------------------- 
REAL(8) FUNCTION gavrg(ibcs,iz,in,iprot,rhosc,n1,n2,xmeff,a0)

USE paramr
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: ibcs
INTEGER, INTENT(IN)                      :: iz
INTEGER, INTENT(IN)                      :: in
INTEGER, INTENT(IN)                      :: iprot
REAL(8), INTENT(IN OUT)                  :: rhosc
INTEGER, INTENT(IN)                      :: n1
INTEGER, INTENT(IN)                      :: n2
REAL(8), INTENT(IN)                      :: xmeff
REAL(8), INTENT(IN)                      :: a0

REAL(8) :: alphap(2)=(/5.13,-2.35/),bethap(2)=(/-1.28E-2,-2.38E-4/)

INTEGER :: np,ntot
REAL(8) :: a,alpha,betha,delavr,eta,fac,hom,part3,pieta,r,robar,tz,tz2,xm,xmm,y1,y2

a=iz+in

ntot=iz
IF (iprot == 0) ntot=in

IF(rhosc == zero) THEN
  
! pair density at the fermi level
  alpha=alphap(iprot+1)
  betha=bethap(iprot+1)
  fac=(2.*ntot/a)**(1.d0/3.)
  part3=(3.*ntot)**(2.0/3.)
  hom=41/a**(1.0/3)*fac
  tz=(in-iz)/a
  tz2=tz*tz
  rhosc=part3/hom/(1.-1.0/part3)/2.*(1.+alpha*tz2+betha*a*tz2)
  
! here: correction for the effective mass
  r=1.14*a**(1.d0/3)
  eta=a0/r
  xm=LOG(xmeff)
  pieta=(3.1415926*eta)**2
  xmm=xmeff*(1.-eta*xm)**3*(1.+2*eta*xm*pieta/(1.+pieta))
  rhosc=rhosc*xmm**1.5
  
END IF
robar=rhosc


IF (ibcs == 1) THEN
  
! bcs variant
  IF (iprot == 0) THEN
    delavr=13.3/SQRT(a)
  ELSE IF (iprot == 1) THEN
    delavr=13.9/SQRT(a)
  ELSE
    delavr = 0.0D0
  END IF
  
ELSE
  
! ln variant
  IF (iprot == 0) THEN
    delavr=8.6/SQRT(a)
  ELSE IF (iprot == 1) THEN
    delavr=9.9/SQRT(a)
  ELSE
    delavr = 0.0D0
  END IF
  
END IF

np=ntot/2

y1=(-np+n1-1)/robar
y2=(-np+n2  )/robar

gavrg=2./robar/(LOG(SQRT(y2**2+delavr**2)+y2)- LOG(SQRT(y1**2+delavr**2)+y1))

RETURN
END FUNCTION gavrg
