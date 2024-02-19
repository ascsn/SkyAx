MODULE linresponse
!c#include "define.h"

USE paramr
USE dynpar
USE levels
USE grid
USE observ, ONLY: zpemod,zpetri,omgtri
USE steputil, ONLY: damnuc,schmit
IMPLICIT REAL(8) (a-h,o-z)

!  file        : linresp.inc
!  description : global variables of the linear-response module

LOGICAL,PARAMETER :: tprecpri=.false.     !<High precision output for alignments.
!
! final collective properties along the path:
!
REAL(8),SAVE :: b0q(0:2,0:2)         ! deformation and variance
REAL(8),SAVE :: bvar(0:2,0:2)        ! (second index: isospin)
REAL(8),SAVE :: n0q(0:2)             ! particle number and variance
REAL(8),SAVE :: h0q
REAL(8),SAVE :: nvar(0:2)            ! (second index: isospin)

REAL(8),SAVE :: pvar(0:2,0:5)        ! <P^2>  momentum variance
REAL(8),SAVE :: qvar(0:2,0:5)        ! <Q^2>  Qdyn     variance

REAL(8),SAVE :: c(0:2,0:5)           ! Stiffness C

REAL(8),SAVE :: b(0:2,0:5)           ! Mass  B
REAL(8),SAVE :: bn(0:2,0:5)          ! Part. Number Mass  B_N

REAL(8),SAVE :: nrm0(0:2,0:5)        ! <[Beta,P]>
REAL(8),SAVE :: nrm(0:2,0:5)         ! <[Q,P]>
REAL(8),SAVE :: omgsq(0:2,0:5)       ! <[[H,[H,Q],P]>

!>@name Auxiliary fields (without spin).
!>@{
REAL(8) :: wrk1(0:nzd,0:nrd),wrk2(0:nzd,0:nrd)
REAL(8) :: wrk3(0:nzd,0:nrd),wrk4(0:nzd,0:nrd)
!>@}
!>@name Auxiliary fields (with spin).
!>@{
REAL(8) :: psi1(0:nzd,0:nrd,nsd),psi2(0:nzd,0:nrd,nsd)
REAL(8) :: psi3(0:nzd,0:nrd,nsd),psi4(0:nzd,0:nrd,nsd)
!>@}
!>@name Auxiliary storage for one s.p. state.
!>@{
REAL(8) :: psip(0:nzd,0:nrd,nsd),psin(0:nzd,0:nrd,nsd)
!>@}
! temporary
REAL(8) :: ampsi(0:nzd,0:nrd,nsd,npsid,1:2),betampsi  ! store action of l_+,l_

REAL(8) :: fldbuf(0:nzd,0:nrd)             !<Temporary buffer for operators.

!  static parameters

!     MZD,MRD              grid dimensions in r,z
!     MID,MSD              dimension of (iso)spin hilbert
!     MPSID                maximum number of sp-wavefunctions
!     MTD                  maximum number of transitions


!  common block "LRdens" -- response densities

!     rerho                response density           (rho)
!     rekin                response kinetic density   (tau)
!     recur[r,p,z]         response current           (j)
!     respn[r,p,z]         response spin density      (sigma)
!     resoc[r,p,z]         response spin-orbit dens.  (J)
!     reskd[r,p,z]         response kin. spin dens.   (tau) (vector)
!     repair               response pairing-density.  (chi)


!  common block "LRflds" -- response fields

!     relocpot             response local potential   (U)
!     recoupot             response coulomb potential (U_C)
!     reeffmas             response effective mass    (B)
!     redynpot[r,p,z]      response dynamic potential (A)
!     respnpot[r,p,z]      response spin potential    (S)
!     resocpot[r,p,z]      response sp.-orb. pot.     (W)
!     respnmas[r,p,z]      response spin mass         (C)


!  common block "LRtrans" -- transition arrays

!     TRidx                transition table
!     TRlen                length of transition table
!     TRxi                 occ. factor xi for transition
!     TReta                occ. factor eta for transition

!  common block "LRDamp"  --  damping matrices

!     E0                   damping parameter
!     LRDampr              damping matrix in r direction
!     LRDampz              damping matrix in z direction

!  common block "LRpath" -- storage of the collective path

!     CPrho                density
!     CPkin                kinetic density
!     CPsoc[r,z]           spin-orbit-current

!     CPlocpot             stationary local potential      (U)
!     CPcoupot             stationary coulomb potential    (U_c)
!     CPeffmas             stationary effective mass       (B)
!     CPsocpot[r,z]        stationary spin-orbit potential (W)

!     CPconstr             constraint field                (Q)

!     CPpsi                sp-wavefunctions of path
!     CPvsq                occ-amplitudes v*v
!     CPspeng              sp-energies of sp-states
!     CPlen                number of sp-states

!     CPdeltaf             ????
!     CPwtpair             ????
!     CPpairfc             ????

!     CPenergy             total energy
!     CPlambda             constraint energy (lagrange-multiplier)
!     CPfer                fermi energies of path elements
!     CPbeta               total deformation of path elements
!     CPrms                total rms-radius of path elements
!     CPzcm                center of mass
!     CPatnum              total number of nucleons of path elements

!     CPdrpsi              derivatives of sp-wavefunctions
!     CPdppsi
!     CPdzpsi


!  common block "LRbuff" -- buffers for linear-response module

!     LRwrk_scalar         buffer for scalar function
!     LRwrk1_spinor        buffers for spinor functions
!     LRwrk2_spinor

!  setting switching parameters:
!     FULLSPACE
!     D_LRDEBUG
!     D_OLDCRANK

!***********************************************************************

!=======================================================================

!     include static parameters

!=======================================================================

INTEGER,PARAMETER :: mzd=nzd           ! maximum number of gridpoints in z
INTEGER,PARAMETER :: mrd=nrd           ! maximum number of gridpoints in r
INTEGER,PARAMETER :: mid=2             ! dimension of isospin hilbert
INTEGER,PARAMETER :: msd=nsd           ! dimension of spin hilbert
INTEGER,PARAMETER :: mpsid=npsid       ! number of sp-wavefunctions
INTEGER,PARAMETER :: mtd=1000          ! maximum number of transitions per state
INTEGER,PARAMETER :: mtdtot=100000     ! maximum number of transitions
INTEGER,PARAMETER :: lambdamax=3       ! highest collective Mode in cranking
!PARAMETER(mzd=nzd,   mrd   =nrd, mid=2, msd=nsd, mpsid=npsid,  &
!    mtd=1000,  mtdtot=100000,  lambdamax=3)


!=======================================================================

!     common block "damping matrices"

!=======================================================================
!COMMON /lrdamp/ 
REAL(8) :: lrdampr(0:mrd,0:mrd,0:kmpsi),lrdampz(0:mzd,0:mzd,2),&
           lre0inv,lrx0,lrstep


!=======================================================================

!     common block "response densities"

!=======================================================================
!COMMON /lrdens/ rerho, rekin,  &
!    recurr, recurp, recurz, respnr, respnp, respnz,  &
!    resocr, resocp, resocz, reskdr, reskdp, reskdz, repair

REAL(8) :: rerho(0:mzd,0:mrd,mid)  ! response rho
REAL(8) :: rekin(0:mzd,0:mrd,mid)  ! response tau
REAL(8) :: recurr(0:mzd,0:mrd,mid) ! response j_r
REAL(8) :: recurp(0:mzd,0:mrd,mid) ! response j_phi
REAL(8) :: recurz(0:mzd,0:mrd,mid) ! response j_z
REAL(8) :: respnr(0:mzd,0:mrd,mid) ! response sigma_r
REAL(8) :: respnp(0:mzd,0:mrd,mid) ! response sigma_phi
REAL(8) :: respnz(0:mzd,0:mrd,mid) ! response sigma_z
REAL(8) :: resocr(0:mzd,0:mrd,mid) ! response J_r
REAL(8) :: resocp(0:mzd,0:mrd,mid) ! response J_p
REAL(8) :: resocz(0:mzd,0:mrd,mid) ! response J_z
REAL(8) :: reskdr(0:mzd,0:mrd,mid) ! response tau_r
REAL(8) :: reskdp(0:mzd,0:mrd,mid) ! response tau_p
REAL(8) :: reskdz(0:mzd,0:mrd,mid) ! response tau_z
REAL(8) :: repair(0:mzd,0:mrd,mid) ! response pair-dens.


!=======================================================================

!     common block "response fields"

!=======================================================================
!COMMON /lrflds/ relocpot, recoupot, reeffmas,  &
!    redreffmas,redpeffmas,redzeffmas, redynpotr, redynpotp, redynpotz,  &
!    respnpotr, respnpotp, respnpotz, resocpotr, resocpotp, resocpotz,  &
!    respnmasr, respnmasp, respnmasz, locpotwgt1,locpotwgt2,repairpot

REAL(8) :: relocpot(0:mzd,0:mrd,mid)  ! response U
REAL(8) :: recoupot(0:mzd,0:mrd)      ! response U_C
REAL(8) :: reeffmas(0:mzd,0:mrd,mid)  ! response B
REAL(8) :: redzeffmas(0:mzd,0:mrd,mid)  ! response nabla_r B
REAL(8) :: redpeffmas(0:mzd,0:mrd,mid)  ! response nabla_p B
REAL(8) :: redreffmas(0:mzd,0:mrd,mid)  ! response nabla_z B
REAL(8) :: redynpotr(0:mzd,0:mrd,mid) ! response A_r
REAL(8) :: redynpotp(0:mzd,0:mrd,mid) ! response A_phi
REAL(8) :: redynpotz(0:mzd,0:mrd,mid) ! response A_z
REAL(8) :: respnpotr(0:mzd,0:mrd,mid) ! response S_r
REAL(8) :: respnpotp(0:mzd,0:mrd,mid) ! response S_phi
REAL(8) :: respnpotz(0:mzd,0:mrd,mid) ! response S_z
REAL(8) :: resocpotr(0:mzd,0:mrd,mid) ! response W_r
REAL(8) :: resocpotp(0:mzd,0:mrd,mid) ! response W_phi
REAL(8) :: resocpotz(0:mzd,0:mrd,mid) ! response W_z
REAL(8) :: respnmasr(0:mzd,0:mrd,mid) ! response C_r
REAL(8) :: respnmasp(0:mzd,0:mrd,mid) ! response C_phi
REAL(8) :: respnmasz(0:mzd,0:mrd,mid) ! response C_z
REAL(8) :: repairpot(0:mzd,0:mrd,mid) ! response pair-pot.

REAL(8) :: locpotwgt1(0:mzd,0:mrd,mid)! response weight dU_p/drho_q
REAL(8) :: locpotwgt2(0:mzd,0:mrd,mid)! response weight dU_n/drho_q

!=======================================================================

!     common block "transition table"

!=======================================================================
!COMMON /lrtrans/ tridx,  &
!    trxi, treta, trntr,trmmax,trntrtot, trlen,nsttot,origin,partner,nstadd,  &
!    isresact,ismfresp

!INTEGER :: trlen
INTEGER :: trntrtot                ! number of transitions
INTEGER :: trntr(1:mpsid)          ! dto. per state
INTEGER :: tridx(1:mpsid,1:mtd,8)  ! transition table
REAL(8) :: trxi(1:mpsid,1:mtd,2)   ! occ. factor xi
REAL(8) :: treta(1:mpsid,1:mtd,2)  ! occ. factor eta
INTEGER :: origin(mpsid)      ! original WF in K>0 space
INTEGER :: partner(mpsid)     ! Time-Rev. partner
INTEGER :: trmmax             ! maximum m-component
INTEGER :: nsttot             ! total size of ext. sp-base
INTEGER :: nstadd(2)          ! number of added time-reversed
INTEGER :: isresact,ismfresp  !                        states


!=======================================================================

!     common block "collective path"

!=======================================================================
!COMMON /lrpath/ cprho, cpkin, cpsocr, cpsocz,  &
!    cplocpot,cpeffmas,cpsocpotr,cpsocpotz,cpcoupot,  &
!    cpdreffmas,cpdzeffmas,cpconstr, cppsi,cpdrpsi,cpdppsi,cpdzpsi,cpddpsi,  &
!    p2qatdhf,q2qatdhf,p2qchf,q2qchf,  &
!    cpdeltaf,cpwtpair,cpspenrg,cpqpenrg,cpfer,cpvsqf,  &
!    cpenergy,cplambda,cpbeta,cpzcm, cpatnum,cpradrms,cpraddif,cpradsrf,  &
!    cpenergyc,cplen

REAL(8) :: cprho(0:mzd,0:mrd,mid,-1:1)       ! density
REAL(8) :: cpkin(0:mzd,0:mrd,mid,-1:1)       ! kinetic dens.
REAL(8) :: cpsocr(0:mzd,0:mrd,mid,-1:1)      ! sp.-orb. current
REAL(8) :: cpsocz(0:mzd,0:mrd,mid,-1:1)      !
REAL(8) :: cplocpot(0:mzd,0:mrd,mid,-1:1)    ! local potential
REAL(8) :: cpeffmas(0:mzd,0:mrd,mid,-1:1)    ! effective mass
REAL(8) :: cpsocpotr(0:mzd,0:mrd,mid,-1:1)   ! sp.-orb. potential
REAL(8) :: cpsocpotz(0:mzd,0:mrd,mid,-1:1)   !
REAL(8) :: cpcoupot(0:mzd,0:mrd,-1:1)        ! coulomb potential
REAL(8) :: cpdreffmas(0:mzd,0:mrd,mid)       ! r-deriv. of mass
REAL(8) :: cpdzeffmas(0:mzd,0:mrd,mid)       ! z-deriv. of mass

REAL(8) :: cpconstr(0:mzd,0:mrd,-1:1)        ! constraint field

REAL(8) :: cppsi(0:mzd,0:mrd,msd,mpsid,-1:1) ! sp-waves
REAL(8) :: cpdrpsi(0:mzd,0:mrd,msd,mpsid)    ! r-derivatives
REAL(8) :: cpdppsi(0:mzd,0:mrd,msd,mpsid)    ! z-derivatives
REAL(8) :: cpdzpsi(0:mzd,0:mrd,msd,mpsid)    ! phi-derivatives
REAL(8) :: cpddpsi(0:mzd,0:mrd,msd,mpsid)    ! laplace

REAL(8) :: p2qatdhf(mtdtot,0:lambdamax)
REAL(8) :: q2qatdhf(mtdtot,0:lambdamax)
REAL(8) :: p2qchf(mtdtot,0:lambdamax)
REAL(8) :: q2qchf(mtdtot,0:lambdamax)

REAL(8) :: cpbeta(-1:1)           ! deformations
REAL(8) :: cpzcm(-1:1)            ! center of mass
REAL(8) :: cpatnum(-1:1)          ! numbr. of nucl.
REAL(8) :: cpradrms(-1:1)         ! rms-radius
REAL(8) :: cpraddif(-1:1)         ! diff.-radius
REAL(8) :: cpradsrf(-1:1)         ! sigma
REAL(8) :: cpenergy(-1:1)         ! total energy
REAL(8) :: cpenergyc(-1:1)        ! total constraint energy
REAL(8) :: cplambda(-1:1)         ! constraint energy
REAL(8) :: cpzpemod(-1:1)         ! c.m. energy correction

REAL(8) :: cpdeltaf(mpsid,-1:1) ! ????
REAL(8) :: cpwtpair(mpsid,-1:1) ! ????
REAL(8) :: cppairfc(mpsid,-1:1) ! ????
REAL(8) :: cpspenrg(mpsid,-1:1) ! sp-energies
REAL(8) :: cpqpenrg(mpsid,-1:1) ! qp-energies
REAL(8) :: cpvsqf(mpsid,-1:1)   ! occ-amplitudes v**2
REAL(8) :: cpfer(mid,-1:1)        ! fermi energy
INTEGER :: cplen                  ! number of path elements


!=======================================================================

!     common block "workspaces"

!=======================================================================
!COMMON /lrbuff/ lrres,lrbuffer,lrwrk_scalar,lrwrk2_scalar,  &
!    lrwrk1_spinor,lrwrk2_spinor, lrwrk3_spinor,lrwrk4_spinor

REAL(8) :: lrres(0:mzd,0:mrd,msd,mpsid,2)
REAL(8) :: lrbuff1(0:mzd,0:mrd,msd,2)
REAL(8) :: lrbuff2(0:mzd,0:mrd,msd,2)
REAL(8) :: lrwrk_scalar(0:mzd,0:mrd)      ! scalar buffer
REAL(8) :: lrwrk2_scalar(0:mzd,0:mrd)     ! scalar buffer
REAL(8) :: lrwrk1_spinor(0:mzd,0:mrd,msd) ! spinor buffers
REAL(8) :: lrwrk2_spinor(0:mzd,0:mrd,msd)
REAL(8) :: lrwrk3_spinor(0:mzd,0:mrd,msd)
REAL(8) :: lrwrk4_spinor(0:mzd,0:mrd,msd)



!=======================================================================

!     setting some switches

!=======================================================================

!LOGICAL :: fullspace
LOGICAL,PARAMETER :: d_lrdebug=.false., d_oldcrank=.TRUE.
LOGICAL,PARAMETER :: tangmomsub=.FALSE.    ! a.m. subroutines

!     parameter (FULLSPACE=0)
!PARAMETER (d_lrdebug=.false.)
!PARAMETER (d_oldcrank=.false.)



CONTAINS
!=======================================================================

!     linresp()

!       main routine for the cranking module of SkyAx. The subroutine
!       calculates the inverse masses B of a collective Hamiltonian
!       and the curvature (stiffness) of the potential C

!                                          1   d^2 V
!            H = B P^2 + V(Q)         C =  - * -----
!                                          2   dQ^2

!       with

!            P  = 1/(2B) * [H,Q]       C = 1/2 * [P,[H,P]]
!            B  = 1/2 * [Q,[H,Q]]
!            C  = 1/2 * [P,[H,P]]

!       for different models

!        1) simple mass from multipole moments
!        2) inglis cranking
!        3) self-consistent crankting (P from CHF)
!        3) ATDHF cranking (self-consistant P,Q-operators)



!       the phonon energy of the mode described by P,Q is
!       defined by

!            [H,[H,Q]] = omega^2 Q     E_0 = 1/2 * omega


!       Q_ij are choosen real ==> P_ij == purely complex
!       the variable "x2q",etc. stores

!           q_ij = Q_ij          (if Q time-even)
!           p_ij = -i/hbar P_ij  (if P time-odd)

!       If the response of the MF is included the
!       commutators become

!           [H,Q] = [H_0,Q] + H^{resp}_Q
!           [H,P] = [H_0,P] + H^{resp}_P

!       where H_0 is diagonal in quasiparticle-space

!           H_0 = e_i alpha^+_i alpha_i

!       and therefore the 2q-part becomes

!           [H,Q]^{2q} = (e_i + e_j) Q_ij + (H^{resp}_Q)_ij
!           [H,P]^{2q} = (e_i + e_j) P_ij + (H^{resp}_P)_ij

!       usefull commutators:

!        <Q^2>     = 1/2 * Q_ij Q_ij       <P^2> = 1/2 * P_ij P_ij




!=======================================================================

SUBROUTINE linresp()

!     ------------------------------------------------------------------

USE paramr
USE constr
USE observ
USE wavefr
USE levels
USE dynpar, ONLY: outnam
IMPLICIT REAL(8) (a-h,o-z)

!INCLUDE 'linresp.inc'

!INCLUDE 'observ.inc'
!INCLUDE 'constr.inc'

!     ------------------------------------------------------------------


!CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam


REAL(8) :: b2q(1:mtdtot)
REAL(8) :: p2q(1:mtdtot)
REAL(8) :: q2q(1:mtdtot)

REAL(8) :: h2q(1:mtdtot)
REAL(8) :: n2q(1:mtdtot)


REAL(8) :: norm0(0:2),norm0qp(0:2),norm0res(0:2)

REAL(8) :: d2energy,d2radrms,d2raddif,d2radsrf
REAL(8) :: d2energyc
REAL(8) :: buff(10),buffer

INTEGER :: imu,i,ibuff
REAL(8) :: betapow

INTEGER :: il
INTEGER :: pathlen

!     ------------------------------------------------------------------
!     switches
!     ------------------------------------------------------------------
ismfresp=1          ! switch to turn on/off MF-response
il      =2          ! collective mode (2=quadrupole)

IF (swicon > oh) THEN
  pathlen = 4
ELSE
  pathlen = 1   ! 4
END IF


!     update the collective path
!     ------------------------------------------------------------------
CALL upd_path()

!     static afterburn to get a good sp-basis also for unoccupied states
!     ------------------------------------------------------------------
IF(pathlen > 1) THEN
!  CALL static_afterburn(+1,0.01D0,1000)
  CALL static_afterburn(+1,0.01D0,1000)
ELSE
!  CALL static_afterburn(+1,0.01D0,200)
  CALL static_afterburn(+1,0.01D0,200)
END IF
CALL return_wavefunctions(+1)

!     update deformations off last path element (to keep consistence)
!     ------------------------------------------------------------------
!      call ini_transtab(il,0)
!      call ini_matrix_B(B2q,CPbeta(+1),buffer,il,0,0,+1)

!     output the deformations kept in the path storage
!     ------------------------------------------------------------------
WRITE(21,'(a,3f8.4)')  &
    'LINRESP()    : elements of collective path -- beta= ', (cpbeta(i),i=-1,1)
WRITE(6,'(a,3f8.4)')  &
    'LINRESP()    : elements of collective path -- beta= ', (cpbeta(i),i=-1,1)

!     check if path-length is to small
!     ------------------------------------------------------------------
IF (cplen < pathlen) THEN
!      if (CPlen.ne.6) then
  WRITE(21,'(a)') 'LINRESP()    : SKIPPED linear response  -- '//  &
      'collective path incomplete'
  RETURN
END IF

!     prepare derivatives of central element of path
!     ------------------------------------------------------------------
CALL upd_derivatives()

!     update non-linear response weight dU_q'/drho_q
!     ------------------------------------------------------------------
CALL upd_locpotwgt()

!     calculate quadrupole mass for m=0,2
!     ------------------------------------------------------------------
DO imu=0,il,1
!      do imu=0,0,1
  
!        write header to cranking log
!        ------------------------------------------------------------------
  WRITE(21,'(a)')   ' ------------------------------------------'
  WRITE(21,'(a,i2,a,i2)') ' cranking for lambda = ',il,'  mu = ',imu
  WRITE(21,'(a)')   ' ------------------------------------------'
  WRITE(21,'(a,7(a16))') 'METHOD    : ',  &
      'B','C','<P^2>','<Q^2>-Q^2','Omega^2', '<[Q,P]>','<[alpha,P]>'
  
!        initialize transitions
!        ---------------------------------------------------------------
  CALL ini_transtab(il,imu)
  
!        initialize quasiparticle energies, beta-operator
!        ---------------------------------------------------------------
  CALL ini_matrix_h(h2q,h0q,il,imu)
  CALL ini_matrix_b(b2q,b0q(imu,0),bvar(imu,0),il,imu,0,0)
  IF (imu == 0) THEN
    CALL ini_matrix_n(n2q, n0q(0),nvar(0),il,imu,0)
  ELSE
    CALL ini_matrix_n(n2q,buff(1),buff(2),il,imu,0)
  END IF
  
  IF (cplen == pathlen) THEN
    DO itr=1,trntrtot
      q2qchf(itr,imu)  =b2q(itr)/h2q(itr)**2
      p2qchf(itr,imu)  =b2q(itr)/h2q(itr)
      q2qatdhf(itr,imu)=b2q(itr)/h2q(itr)**2
      p2qatdhf(itr,imu)=b2q(itr)/h2q(itr)
    END DO
  END IF
  
  
!        calculate coll. mass using inglis and self-consistant cranking
!        ---------------------------------------------------------------
  IF (imu == 0.OR.imu == 1) THEN
    
!           initialize P-operator from CHF (SkyAx, constraint loop)
    CALL ini_matrix_p(p2q,pvar(imu,1),il,imu)
    CALL eval_comm_bp(nrm0(imu,1),il,imu)
    pvar(imu,2)=pvar(imu,1)    ! SC-, INGLIS- and static GCM
    pvar(imu,3)=pvar(imu,1)    ! have the same P-Variance and
    nrm0(imu,2)=nrm0(imu,1)    ! <[beta,P]> commutator
    nrm0(imu,3)=nrm0(imu,1)
    
!           the GCM-mass
    CALL crank_gcm(b2q,q2q,p2q,h2q, b(imu,1),c(imu,1),  &
        qvar(imu,1),pvar(imu,1), il,imu)
    CALL eval_com_qnq(n2q,q2q,bn(imu,1),il,imu,0)
    pvar(imu,1)=pvar(imu,2)
    
!           inglis cranking mass
    CALL crank_inglis(b2q,q2q,p2q,h2q, b(imu,2),c(imu,2),  &
        qvar(imu,2),buff(1), il,imu)
    CALL eval_com_qnq(n2q,q2q,bn(imu,2),il,imu,0)
    
!           selfconsistent cranking mass
    IF(pathlen > 1 .OR. imu == 1) THEN
      WRITE(22,'(a,2i5)') ' PathLen,imu=',pathlen,imu
      CALL crank_selfc(b2q,q2q,p2q,h2q, b(imu,3),c(imu,3),  &
          qvar(imu,3),buff(1), il,imu)
      CALL eval_com_qnq(n2q,q2q,bn(imu,3),il,imu,0)
    END IF
    
!           debug subroutine to test the responses
!            call crank_test_resp(B2q,Q2q,P2q,H2q,
!     &                         B(imu,0),C(imu,0),
!     &                         QVar(imu,0),PVar(imu,0),
!     &                         il,imu)
    
!           the scaling mass
!            call crank_scaling(B2q,Q2q,P2q,H2q,
!     &                         B(imu,0),C(imu,0),
!     &                         QVar(imu,0),PVar(imu,0),
!     &                         il,imu)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,0),il,imu,0)
!c
!c           CHF in linear response
!            call crank_CHF(B2q,Q2q,P2q,H2q,
!     &                     B(imu,4),C(imu,4),
!     &                     QVar(imu,4),PVar(imu,4),
!     &                     il,imu)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
!            call crank_ATDHF(B2q,Q2q,P2q,H2q,
!     &                     B(imu,4),C(imu,4),
!     &                     QVar(imu,4),PVar(imu,4),
!     &                     il,imu,1)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
!c
!c           full ATDHF
!            call crank_ATDHF(B2q,Q2q,P2q,H2q,
!     &                     B(imu,5),C(imu,5),
!     &                     QVar(imu,5),PVar(imu,5),
!     &                     il,imu,0)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,1),il,imu,0)
    
  ELSE IF (imu > 1) THEN
    
!           debug subroutine to test the responses
!            call crank_test_resp(B2q,Q2q,P2q,H2q,
!     &                         B(imu,0),C(imu,0),
!     &                         QVar(imu,0),PVar(imu,0),
!     &                         il,imu)
    
!           the scaling mass
!            call crank_scaling(B2q,Q2q,P2q,H2q,
!     &                         B(imu,0),C(imu,0),
!     &                         QVar(imu,0),PVar(imu,0),
!     &                         il,imu)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,0),il,imu,0)
!c
!c           CHF in linear response
!            call crank_CHF(B2q,Q2q,P2q,H2q,
!     &                     B(imu,4),C(imu,4),
!     &                     QVar(imu,4),PVar(imu,4),
!     &                     il,imu)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
!            call crank_ATDHF(B2q,Q2q,P2q,H2q,
!     &                     B(imu,4),C(imu,4),
!     &                     QVar(imu,4),PVar(imu,4),
!     &                     il,imu,1)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
!c
!c           the GCM-mass
!c            call crank_GCM(B2q,Q2q,P2q,H2q,
!c     &                        B(imu,1),C(imu,1),
!c     &                        QVar(imu,1),buff(1),
!c     &                        il,imu)
!c            call eval_com_QNQ(N2q,Q2q,BN(imu,1),il,imu,0)
!c
!c           full ATDHF
!            call crank_ATDHF(B2q,Q2q,P2q,H2q,
!     &                     B(imu,5),C(imu,5),
!     &                     QVar(imu,5),PVar(imu,5),
!     &                     il,imu,0)
!            call eval_com_QNQ(N2q,Q2q,BN(imu,5),il,imu,0)
    
  END IF
  
!        calculate widths (and deformations) for particle sorts seperatly
!        ---------------------------------------------------------------
  CALL ini_matrix_b(b2q,b0q(imu,1),bvar(imu,1),il,imu,1,0) ! protons
  CALL ini_matrix_b(b2q,b0q(imu,2),bvar(imu,2),il,imu,2,0) ! neutrons
  
!         PVar(imu,2)=PVar(imu,1)     ! <P^2> is equal for GCM,SC,INGLIS
!         PVar(imu,3)=PVar(imu,1)     !
!         Nrm0(imu,2)=Nrm0(imu,1)     ! <[Q,P]> is equal for GCM,SC,INGLIS
!         Nrm0(imu,3)=Nrm0(imu,1)
  
!        calculate particle numbers and variances
!        ---------------------------------------------------------------
  IF (imu == 0) THEN
    CALL ini_matrix_n(n2q,n0q(0),nvar(0),il,imu,0)
    CALL ini_matrix_n(n2q,n0q(1),nvar(1),il,imu,1)
    CALL ini_matrix_n(n2q,n0q(2),nvar(2),il,imu,2)
  END IF
  
END DO

!     write collective parameters of Hamiltonian
!     ------------------------------------------------------------------

IF(tprecpri) THEN
  DO i=0,5
    WRITE(23+i,'(f15.12,21(1pg20.12))') cpbeta(0),  &
        b(0,i), c(0,i), pvar(0,i), qvar(0,i), omgsq(0,i),  &
        b(1,i), c(1,i), pvar(1,i), qvar(1,i), omgsq(1,i),  &
        b(2,i), c(2,i), pvar(2,i), qvar(2,i), omgsq(2,i),  &
        bn(0,i),bn(1,i),bn(2,i), nrm0(0,i),nrm0(1,i),nrm0(2,i)
  END DO
ELSE
  DO i=0,5
    WRITE(23+i,'(f8.5,21(1pg12.4))') cpbeta(0),  &
        b(0,i), c(0,i), pvar(0,i), qvar(0,i), omgsq(0,i),  &
        b(1,i), c(1,i), pvar(1,i), qvar(1,i), omgsq(1,i),  &
        b(2,i), c(2,i), pvar(2,i), qvar(2,i), omgsq(2,i),  &
        bn(0,i),bn(1,i),bn(2,i), nrm0(0,i),nrm0(1,i),nrm0(2,i)
  END DO
END IF

!     calculate curvature directly from path
!     ------------------------------------------------------------------
d2energy        = (cpenergy(+1) - 2.0*cpenergy( 0) +cpenergy(-1))  &
    /(0.5*(cpbeta( +1) - cpbeta(-1)))**2/2.0
d2energyc       = (cpenergyc(+1)-2.0*cpenergyc( 0)+cpenergyc(-1))  &
    /(0.5*(cpbeta( +1) - cpbeta(-1)))**2/2.0
d2radrms        = (cpradrms(+1) - 2.0*cpradrms( 0) +cpradrms(-1))  &
    /(0.5*(cpbeta( +1) - cpbeta(-1)))**2/2.0
d2raddif        = (cpraddif(+1) - 2.0*cpraddif( 0) +cpraddif(-1))  &
    /(0.5*(cpbeta( +1) - cpbeta(-1)))**2/2.0
d2radsrf        = (cpradsrf(+1) - 2.0*cpradsrf( 0) +cpradsrf(-1))  &
    /(0.5*(cpbeta( +1) - cpbeta(-1)))**2/2.0

!     write the radii
!     -----------------------------------------------------------------
IF(tprecpri) THEN
  WRITE(29,'(6(1pg20.12))') cpbeta(0),cpradrms(0),cpraddif(0),cpradsrf(0),  &
      d2energy, d2energyc
ELSE
  WRITE(29,'(6(1pg14.5))') cpbeta(0),cpradrms(0),cpraddif(0),cpradsrf(0),  &
      d2energy, d2energyc
END IF

!     write the deformations and deformation-widths (and part. number)
!     -----------------------------------------------------------------
IF(tprecpri) THEN
  WRITE(30,'(20(1pg20.12))') b0q(0,0),b0q(0,1),b0q(0,2),  &
      bvar(0,0),bvar(1,0),bvar(2,0), bvar(0,1),bvar(1,1),bvar(2,1),  &
      bvar(0,2),bvar(1,2),bvar(2,2), n0q(0),n0q(1),n0q(2),  &
      nvar(0),nvar(1),nvar(2), pvar(0,1),pvar(1,1)
ELSE
  WRITE(30,'(20(1pg14.5))') b0q(0,0),b0q(0,1),b0q(0,2),  &
      bvar(0,0),bvar(1,0),bvar(2,0), bvar(0,1),bvar(1,1),bvar(2,1),  &
      bvar(0,2),bvar(1,2),bvar(2,2), n0q(0),n0q(1),n0q(2),  &
      nvar(0),nvar(1),nvar(2), pvar(0,1),pvar(1,1)
END IF

!     write the normalisation
!     -----------------------------------------------------------------
IF(tprecpri) THEN
!  WRITE(32,'(9(1pg22.12))') b0q(0,0),  &
!      nrm0(0,1),nrm0(1,1),          ! GCM and CHF  [Beta,P]  &
!      nrm0(0,0),nrm0(1,0),nrm0(2,0),! Multipole    [Beta,P]  &
!      nrm0(0,5),nrm0(1,5),nrm0(2,5) ! ATDHF        [Beta,P]
  WRITE(32,'(9(1pg22.12))') b0q(0,0),  &
      nrm0(0,1),nrm0(1,1),             &
      nrm0(0,0),nrm0(1,0),nrm0(2,0),   &
      nrm0(0,5),nrm0(1,5),nrm0(2,5)  
ELSE
  WRITE(32,'(9(1pg14.5))') b0q(0,0),  &
      nrm0(0,1),nrm0(1,1),            &
      nrm0(0,0),nrm0(1,0),nrm0(2,0),  &
      nrm0(0,5),nrm0(1,5),nrm0(2,5) 
END IF
CALL FLUSH(32)

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE linresp


!=======================================================================

!     cranking()

!     calculate dynamic Q-operator and P-operator by different cranking
!     methods.

!     Different Methods are supported:

!       0)  CHF:           [H,[H,Q]]             = -2*dlambda/dq*B Beta
!       1)  GCM:                 Q               =   -i P/(2*<P^2>)
!       2)  inglis:           [H,Q]              =   i  2B     P
!       3)  self-cons.:   LQ= [H,Q] + H^{resp}_Q =   i  2B     P
!       4)  ATDHF:            L^2Q               =   Omega^2   Q

!     H is the static Hamiltonian, Beta is the messurement operator
!     of the collective coordiante. L is the Liouvillian

!                      LA = [H,A] + H^{resp}_A

!     While CHF and ATDHF locally define the collective path by
!     themselves GCM, inglis and self-cons. require a given P-Operator.

!=======================================================================
!      subroutine cranking(type,imu,H2q,
!     &                    P2q,resP,Q2q,resQ,B2q,resB,
!     &                    B,BN,C,PVar,QVar,OmgSq,Nrm,Nrm0)
!c
!c     ------------------------------------------------------------------
!c
!      include 'paramr.inc'
!      include 'dynpar.inc'
!      include 'linresp.inc'
!      include 'levels.inc'
!c
!c     ------------------------------------------------------------------
!c
!      double precision H2q(1:MTD)
!      double precision P2q(1:MTD)
!      double precision Q2q(1:MTD)
!      double precision B2q(1:MTD)
!      double precision HP2q(1:MTD)
!      double precision HQ2q(1:MTD)
!#IF(FULLSPACE)
!      double precision ResP(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision ResQ(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision ResB(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision ResHP(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision ResHQ(0:MZD,0:MRD,MSD,MPSID,2)

!      double precision drRes(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision dzRes(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision dpRes(0:MZD,0:MRD,MSD,MPSID,2)
!      double precision ddRes(0:MZD,0:MRD,MSD,MPSID,2)
!#ELSE
!      double precision ResP
!      double precision ResQ
!      double precision ResB
!      double precision ResHP
!      double precision ResHQ

!      double precision drRes
!      double precision dzRes
!      double precision dpRes
!      double precision ddRes
!#ENDIF

!      double precision B,BN,C
!      double precision PVar,BVar,QVar,OmgSq
!      double precision Nrm,Nrm0
!      integer          imu,type

!      integer          max_citer
!      double precision stp_res_Q  , stp_qp_Q  , max_err_Q
!      double precision stp_res_P  , stp_qp_P  , max_err_P
!      double precision stp_res_Omg, stp_qp_Omg, max_err_Omg

!      integer          ir,iz,itr,ist,iter,buff2,lograte
!      integer          ibpar
!      double precision norm,error_P,error_Q,buffer
!      double precision var1,var2
!      double precision Pnorm(MPSID),Qnorm(MPSID),Bnorm(MPSID)
!      double precision HPnorm(MPSID),HQnorm(MPSID)
!      double precision rho,rhomax,r,z,buff,buff1

!c
!c     ------------------------------------------------------------------
!c
!      error_P  =0.0D0
!      error_Q  =0.0D0
!      error_Omg=0.0D0
!c
!c     iteration parameters
!c     ------------------------------------------------------------------
!      max_citer    = 1000   ! maximum number of cranking iterations
!      stp_res_P    = LRStep ! gradient steps in residual components
!      stp_res_Q    = LRStep
!      stp_res_Omg  = LRStep**2
!      stp_qp_P     = 0.8    ! gradient steps in qp-components
!      stp_qp_Q     = 0.8
!      stp_qp_Omg   = 0.8
!      max_err_P    = 1.0E-12 ! maximum relative error || P - LQ   ||/<P^2>
!      max_err_Q    = 1.0E-12 ! maximum relative error || Q - LP   ||/<Q^2>
!      max_err_Omg  = 1.0E-12 ! maximum relative error || Q - L^2Q ||/<Omg^2>
!      lograte      = 1       ! logging rate
!c
!c     check for unsupported cranking method
!c     ------------------------------------------------------------------
!      if (type.lt.0.or.type.gt.5) stop 'unsupported cranking method'
!c
!c     ==================================================================
!c     now the five different methods:
!c     ==================================================================
!c
!c     ------------------------------------------------------------------
!c     CHF + SELFCONSISTANT CRANKING                     (requires Beta)
!c     ------------------------------------------------------------------
!      if (type.eq.0) then
!c
!c        initialization
!c        ---------------------------------------------------------------
!         BVar   = var_oper(B2q,ResB,imu)
!c
!c        find Q-operator by {\cal L}^2 Q = 4 BC Q
!c        ---------------------------------------------------------------
!         do itr=1,TRlen
!            Q2q(itr)= B2q(itr)/H2q(itr)**2 ! Q_ij = B_ij/(e_i + e_j)**2
!         enddo
!c
!#IF(FULLSPACE)
!c        inversion of {\cal L}_0 in residual space
!c        ---------------------------------------------------------------
!c        initialize using the B-Operator
!         if (isresact.eq.1) then
!c
!c           initialize using the B-Operator
!            do ist=1,nst(1)+nst(2)
!               do ir=0,nr
!                  do iz=0,nz
!                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
!                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
!                   ResQ(iz,ir,1,ist,2)=ResB(iz,ir,1,ist,2)
!                   ResQ(iz,ir,2,ist,2)=ResB(iz,ir,2,ist,2)
!                  enddo
!               enddo
!            enddo
!c
!c           calculate norm of residual componentes B-Operator
!            do ist=1,nst(1)+nst(2)
!               Bnorm(ist)=
!     &              wfovlp(ResB(0,0,1,ist,1),ResB(0,0,1,ist,1))
!     &             +wfovlp(ResB(0,0,1,ist,2),ResB(0,0,1,ist,2))
!            enddo
!c
!         else
!c
!c           set residual components to zero
!            call zero_oper(Q2q,ResQ,2)
!            call zero_oper(HP2q,ResHP,2)
!c
!         endif
!#ENDIF
!c
!c        ===============================================================
!c        the cranking iteration (find Q-Operator)
!c        ===============================================================
!         do iter=1,max_citer
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           calculate derivatives of resiual part of Q
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!               call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                              ddRes(0,0,1,1,1),
!     &                              (-1)**imu,+imu,nst(1)+nst(2))
!               call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                              ddRes(0,0,1,1,2),
!     &                              (-1)**imu,-imu,nst(1)+nst(2))
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           calculate action of Liouville-Operator
!c           ------------------------------------------------------------
!            call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
!     &                    HQ2q,ResHQ,imu)
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           calculate derivatives of resiual part of LQ \propto P
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!               call upd_derive_spinor(ResHQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                              ddRes(0,0,1,1,1),
!     &                              (-1)**imu,+imu,nst(1)+nst(2))
!               call upd_derive_spinor(ResHQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                              ddRes(0,0,1,1,2),
!     &                              (-1)**imu,-imu,nst(1)+nst(2))
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           calculate action of Liouville-Operator
!c           ------------------------------------------------------------
!            call comm_HP(H2q,HQ2q,ResHQ,drRes,dzRes,dpRes,ddRes,
!     &                    HP2q,ResHP,imu)
!c
!c           ------------------------------------------------------------
!c           perform the gradient step in qp-space
!c           ------------------------------------------------------------
!            do itr=1,TRLen
!               HQ2q(itr) = HP2q(itr)-B2q(itr)
!               Q2q(itr)  =
!     &         Q2q(itr)  - stp_qp_Omg*(HP2q(itr)-B2q(itr))/H2q(itr)**2
!            enddo
!c
!c           ------------------------------------------------------------
!c           evaluate errors (only if residual parts are neglected)
!c           ------------------------------------------------------------
!            if (isresact.eq.0) then
!               buffer   = error_Q
!               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
!               derror_Q = error_Q-buffer
!            endif
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           perform the gradient step in residual space
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!c
!c              ---------------------------------------------------------
!c              calculate norm of residual components of [H,Q]-Operator
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  HPnorm(ist)=
!     &                 wfovlp(resHP(0,0,1,ist,1),resHP(0,0,1,ist,1))+
!     &                 wfovlp(resHP(0,0,1,ist,2),resHP(0,0,1,ist,2))
!               enddo
!c
!c              ---------------------------------------------------------
!c              calculate errors in residual space
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  do ir=0,nr
!                     do iz=0,nz
!                        ResHQ(iz,ir,1,ist,1) = ResHP(iz,ir,1,ist,1)-
!     &                                         ResB(iz,ir,1,ist,1)
!                        ResHQ(iz,ir,2,ist,1) = ResHP(iz,ir,2,ist,1)-
!     &                                         ResB(iz,ir,2,ist,1)
!                        ResHQ(iz,ir,1,ist,2) = ResHP(iz,ir,1,ist,2)-
!     &                                         ResB(iz,ir,1,ist,2)
!                        ResHQ(iz,ir,2,ist,2) = ResHP(iz,ir,2,ist,2)-
!     &                                         ResB(iz,ir,2,ist,2)
!                     enddo
!                  enddo
!               enddo
!c
!               buffer   = error_Q
!               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
!               derror_Q = error_Q-buffer
!c
!c              ---------------------------------------------------------
!c              apply damping
!c              ---------------------------------------------------------
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHQ(0,0,1,1,1),ResHP(0,0,1,1,1),
!     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHP(0,0,1,1,1),ResHQ(0,0,1,1,1),
!     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHQ(0,0,1,1,2),ResHP(0,0,1,1,2),
!     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHP(0,0,1,1,2),ResHQ(0,0,1,1,2),
!     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))

!c
!c              ---------------------------------------------------------
!c              the gradient step
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  do ir=0,nr
!                     do iz=0,nz
!                      ResQ(iz,ir,1,ist,1) =
!     &                ResQ(iz,ir,1,ist,1) - ResHQ(iz,ir,1,ist,1)
!                      ResQ(iz,ir,2,ist,1) =
!     &                ResQ(iz,ir,2,ist,1) - ResHQ(iz,ir,2,ist,1)
!                      ResQ(iz,ir,1,ist,2) =
!     &                ResQ(iz,ir,1,ist,2) - ResHQ(iz,ir,1,ist,2)
!                      ResQ(iz,ir,2,ist,2) =
!     &                ResQ(iz,ir,2,ist,2) - ResHQ(iz,ir,2,ist,2)
!                     enddo
!                  enddo
!               enddo
!c
!c              ---------------------------------------------------------
!c              norm stabilization
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  norm=sqrt(Bnorm(ist)/HPnorm(ist))
!                  do ir=0,nr
!                     do iz=0,nz
!                        ResQ(iz,ir,1,ist,1) = ResQ(iz,ir,1,ist,1)*norm
!                        ResQ(iz,ir,2,ist,1) = ResQ(iz,ir,2,ist,1)*norm
!                        ResQ(iz,ir,1,ist,2) = ResQ(iz,ir,1,ist,2)*norm
!                        ResQ(iz,ir,2,ist,2) = ResQ(iz,ir,2,ist,2)*norm
!                     enddo
!                  enddo
!               enddo
!c
!c              ---------------------------------------------------------
!c              project out sp-contributions
!c              ---------------------------------------------------------
!               call project_Omega(ResQ,(-1)**imu,imu,1,'Q')
!c
!c           ------------------------------------------------------------
!c           END of gradient step in residual space
!c           ------------------------------------------------------------
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           Logging and convergency-check
!c           ------------------------------------------------------------
!c
!            if (mod(iter,lograte).eq.0               .or.
!     &          error_Q.le.max_err_Q*sqrt(abs(BVar)) .or.
!     &          iter.eq.max_citer ) then
!c
!              write(21,'(a,i5,1(1pg16.4),2(e16.3))') 'CHF+SELFC  : ',
!     &              iter,error_Q,error_Q/sqrt(abs(BVar)),
!     &              derror_Q
!c
!               if (error_Q.le.max_err_Q*sqrt(abs(BVar))) goto 2100
!               if (abs(derror_Q).lt.1.0E-10) goto 2100
!c
!            endif
!c
!         enddo
!c        ===============================================================
!c        END of cranking iteration (find Q-Operator)
!c        ===============================================================
! 2100    continue
!c
!#IF(FULLSPACE)
!c        ------------------------------------------------------------
!c        calculate derivatives of residual part of Q
!c         ------------------------------------------------------------
!         if (isresact.eq.1) then
!            call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                           dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                           ddRes(0,0,1,1,1),
!     &                           (-1)**imu,+imu,nst(1)+nst(2))
!            call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                            dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                            ddRes(0,0,1,1,2),
!     &                            (-1)**imu,-imu,nst(1)+nst(2))
!         endif
!#ENDIF
!c
!c        ------------------------------------------------------------
!c        calculate action of Liouville-Operator to obtain P-Operator
!c        ------------------------------------------------------------
!         call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
!     &                P2q,ResP,imu)
!c
!c     ------------------------------------------------------------------
!c     GCM-formula:                                    (requires P,Beta)
!c     ------------------------------------------------------------------
!      else if (type.eq.1) then
!c
!c        GCM-Mass requires P-Variance
!#IF(FULLSPACE)
!         PVar   = var_oper(P2q,ResP,imu)
!#ELSE
!         call comm_and_var_QP(buff,PVar,buff1,CPbeta(0),
!     &                        (-1)**imu,imu,0)
!#ENDIF
!c
!c        the quasiparticle part of the GCM Q-Operator
!         do itr=1,TRlen
!            Q2q(itr)= P2q(itr)/(2.0D0*PVar)
!         enddo
!c
!#IF(FULLSPACE)
!c        the residual part of the GCM Q-Operator
!         if (isresact.eq.1) then
!            do ist=1,nst(1)+nst(2)
!               do ir=0,nr
!                  do iz=0,nz
!                  ResQ(iz,ir,1,ist,1)=ResP(iz,ir,1,ist,1)/(2.0D0*PVar)
!                  ResQ(iz,ir,2,ist,1)=ResP(iz,ir,2,ist,1)/(2.0D0*PVar)
!                  ResQ(iz,ir,1,ist,2)=ResP(iz,ir,1,ist,2)/(2.0D0*PVar)
!                  ResQ(iz,ir,2,ist,2)=ResP(iz,ir,2,ist,2)/(2.0D0*PVar)
!                  enddo
!               enddo
!            enddo
!         else
!            call zero_oper(Q2q,ResQ,2)  ! set residual part to zero
!         endif
!#ENDIF
!c
!c     ------------------------------------------------------------------
!c     inglis-formula:                                 (requires P,Beta)
!c     ------------------------------------------------------------------
!      else if (type.eq.2) then
!c
!c        temporary disable MF-response
!c        ---------------------------------------------------------------
!         buff2=isMFresp
!         isMFresp=0
!c
!c        initialization
!c        ---------------------------------------------------------------
!#IF(FULLSPACE)
!         PVar   = var_oper(P2q,ResP,imu)
!#ELSE
!         call comm_and_var_QP(buff,PVar,buff1,CPbeta(0),
!     &                        (-1)**imu,imu,0)
!#ENDIF
!c
!c        find Q-operator by {\cal L}_0 Q = [H,Q] = -2i B P
!c        ---------------------------------------------------------------
!c        initialize with inglis values for the 2qp-part
!         do itr=1,TRlen
!            Q2q(itr)= P2q(itr)/H2q(itr)   ! Q_ij = P_ij/(e_i + e_j)
!         enddo
!c
!#IF(FULLSPACE)
!c        inversion of {\cal L}_0 in residual space
!c        ---------------------------------------------------------------
!c        initialize using the B-Operator
!         if (isresact.eq.1) then
!            do ist=1,nst(1)+nst(2)
!               do ir=0,nr
!                  do iz=0,nz
!                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
!                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
!                   ResQ(iz,ir,1,ist,2)=0.0D0*ResB(iz,ir,1,ist,2)
!                   ResQ(iz,ir,2,ist,2)=0.0D0*ResB(iz,ir,2,ist,2)
!                  enddo
!               enddo
!            enddo
!         else
!            call zero_oper(Q2q,ResQ,2)  ! set residual part to zero
!         endif
!c
!c        calculate norm of residual componentes P-Operator
!         do ist=1,nst(1)+nst(2)
!            Pnorm(ist)=
!     &           wfovlp(ResP(0,0,1,ist,1),ResP(0,0,1,ist,1))
!c     &          +wfovlp(ResP(0,0,1,ist,2),ResP(0,0,1,ist,2))
!         enddo
!c
!c        find resQ by iteration
!         if (isresact.eq.1) then
!            do iter=1,max_citer
!c
!c            prepare derivatives of residual part of Q
!             call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                              ddRes(0,0,1,1,1),
!     &                              (-1)**imu,+imu,nst(1)+nst(2))
!c             call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!c     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!c     &                              ddRes(0,0,1,1,2),
!c     &                              (-1)**imu,-imu,nst(1)+nst(2))
!c
!c            calculate action of Liouville-Operator
!             call comm_HQ(H2q,Q2q,resQ,drRes,dzRes,dpRes,ddRes,
!     &                    HQ2q,resHQ,imu)
!c
!c            calculate variance <[H,Q]^2> and <P[H,Q]>
!             var1=var_oper(HQ2q,ResHQ,imu)
!             var2=0.5D0*comm_QP(HQ2q,ResHQ,P2q,ResP,imu)
!c
!c            calculate norm of residual components of [H,Q]-Operator
!             do ist=1,nst(1)+nst(2)
!                HQnorm(ist)=
!     &               wfovlp(resHQ(0,0,1,ist,1),resHQ(0,0,1,ist,1))
!c     &               +wfovlp(resHQ(0,0,1,ist,2),resHQ(0,0,1,ist,2))
!             enddo
!c
!c            calculate error in residual space
!             do ist=1,nst(1)+nst(2)
!                do ir=0,nr
!                   do iz=0,nz
!                      resHP(iz,ir,1,ist,1) = resHQ(iz,ir,1,ist,1)-
!     &                                       resP(iz,ir,1,ist,1)
!                      resHP(iz,ir,2,ist,1) = resHQ(iz,ir,2,ist,1)-
!     &                                       resP(iz,ir,2,ist,1)
!                      resHP(iz,ir,1,ist,2) = 0.0D0!resHQ(iz,ir,1,ist,2)-
!c     &                                       resP(iz,ir,1,ist,2)
!                      resHP(iz,ir,2,ist,2) = 0.0D0!resHQ(iz,ir,2,ist,2)-
!c     &                                       resP(iz,ir,2,ist,2)
!                   enddo
!                enddo
!             enddo
!c
!c            store error of [H,Q]-Operator
!             buffer   = error_P
!             error_P  = sqrt(var_oper(HP2q,resHP,imu))
!             derror_P = error_P-buffer
!c
!c            log components which have big errors
!             do ist=1,nst(1)+nst(2)
!                norm=wfovlp(resHP(0,0,1,ist,1),resHP(0,0,1,ist,1))
!c     &              +wfovlp(resHP(0,0,1,ist,2),resHP(0,0,1,ist,2))
!                if (norm.gt.1E-2) then
!                write (88,*) ist,norm,CPvsqf(ist,0)
!                endif
!             enddo
!c
!c            apply damping
!             call damp_spinor(LRDampr,LRDampz,
!     &                        resHP(0,0,1,1,1),resHQ(0,0,1,1,1),
!     &                        (-1)**imu,imu,stp_res_P,nst(1)+nst(2))
!c             call damp_spinor(LRDampr,LRDampz,
!c     &                        resHP(0,0,1,1,2),resHQ(0,0,1,1,2),
!c     &                        (-1)**imu,-imu,stp_res_P,nst(1)+nst(2))
!             call project_Omega(ResHQ,(-1)**imu,imu)
!c
!c            perform the gradient step in residual space
!             do ist=1,nst(1)+nst(2)
!                do ir=0,nr
!                   do iz=0,nz
!                      ResQ(iz,ir,1,ist,1) =
!     &                ResQ(iz,ir,1,ist,1) - resHQ(iz,ir,1,ist,1)
!                      ResQ(iz,ir,2,ist,1) =
!     &                ResQ(iz,ir,2,ist,1) - resHQ(iz,ir,2,ist,1)
!                      ResQ(iz,ir,1,ist,2) = 0.0D0!
!c     &                ResQ(iz,ir,1,ist,2) - resHQ(iz,ir,1,ist,2)
!                      ResQ(iz,ir,2,ist,2) = 0.0D0!
!c     &                ResQ(iz,ir,2,ist,2) - resHQ(iz,ir,2,ist,2)
!                      enddo
!                   enddo
!             enddo
!c
!c            norm stabilization
!             do ist=1,nst(1)+nst(2)
!                norm=sqrt(Pnorm(ist)/HQnorm(ist))
!                do ir=0,nr
!                   do iz=0,nz
!                      ResQ(iz,ir,1,ist,1) = ResQ(iz,ir,1,ist,1)*norm
!                      ResQ(iz,ir,2,ist,1) = ResQ(iz,ir,2,ist,1)*norm
!                      ResQ(iz,ir,1,ist,2) = 0.0*ResQ(iz,ir,1,ist,2)*norm
!                      ResQ(iz,ir,2,ist,2) = 0.0*ResQ(iz,ir,2,ist,2)*norm
!                   enddo
!                enddo
!             enddo
!c
!c            project out sp-contributions
!             call project_Omega(ResQ,(-1)**imu,imu)
!c
!c            Logging and convergency-check
!c            ---------------------------------------------------------
!             norm     = 1.0D0-var2**2/(Pvar*var1)
!c
!             if (mod(iter,lograte).eq.0 .or.
!     &           error_P.le.max_err_P*sqrt(abs(PVar)) .or.
!     &           iter.eq.max_citer ) then
!c
!               write(21,'(a,i5,2(1pg16.4),3(e16.3))') 'INGLIS   : ',
!     &              iter,error_P,error_P/sqrt(abs(PVar)),
!     &              derror_P,norm
!c
!               if ((error_P.le.max_err_P*sqrt(abs(PVar))).and.
!     &             (norm.lt.1E-4)) goto 2102
!               if (abs(derror_P).lt.1.0E-10) goto 2102
!c
!              endif
!c
!            enddo
! 2102       continue
!c
!         endif
!#ENDIF
!c
!c        reactivate MF-response
!c        ---------------------------------------------------------------
!         isMFresp=buff2
!c
!c     ------------------------------------------------------------------
!c     self-consistent cranking:                       (requires P,Beta)
!c     ------------------------------------------------------------------
!      else if (type.eq.3) then
!c
!#IF(FULLSPACE)
!         PVar   = var_oper(P2q,ResP,imu)
!#ELSE
!         call comm_and_var_QP(buff,PVar,buff1,CPbeta(0),
!     &                        (-1)**imu,imu,0)
!         write(666,*) PVar
!#ENDIF
!c
!c        find Q-operator by {\cal L}_0 Q = [H,Q] = -2i B P
!c        ---------------------------------------------------------------
!c        initialize with inglis values for the 2qp-part
!         do itr=1,TRlen
!            Q2q(itr)= P2q(itr)/H2q(itr)   ! Q_ij = P_ij/(e_i + e_j)
!         enddo
!c
!#IF(FULLSPACE)
!c        inversion of {\cal L}_0 in residual space
!c        ---------------------------------------------------------------
!c        initialize using the B-Operator
!         if (isresact.eq.1) then
!c
!c           initialize using the B-Operator
!            do ist=1,nst(1)+nst(2)
!               do ir=0,nr
!                  do iz=0,nz
!                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
!                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
!                   ResQ(iz,ir,1,ist,2)=ResB(iz,ir,1,ist,2)
!                   ResQ(iz,ir,2,ist,2)=ResB(iz,ir,2,ist,2)
!                  enddo
!               enddo
!            enddo
!c
!c           calculate norm of residual componentes P-Operator
!            do ist=1,nst(1)+nst(2)
!               Pnorm(ist)=
!     &              wfovlp(ResP(0,0,1,ist,1),ResP(0,0,1,ist,1))
!     &             +wfovlp(ResP(0,0,1,ist,2),ResP(0,0,1,ist,2))
!            enddo
!c
!         else
!c
!c           set residual components to zero
!            call zero_oper(Q2q,ResQ,2)
!            call zero_oper(HP2q,ResHP,2)
!c
!         endif
!#ENDIF
!c
!c        ===============================================================
!c        the cranking iteration (find Q-Operator)
!c        ===============================================================
!         do iter=1,max_citer
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           calculate derivatives of resiual part of Q
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!               call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                              ddRes(0,0,1,1,1),
!     &                              (-1)**imu,+imu,nst(1)+nst(2))
!               call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                              ddRes(0,0,1,1,2),
!     &                              (-1)**imu,-imu,nst(1)+nst(2))
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           calculate action of Liouville-Operator
!c           ------------------------------------------------------------
!            call comm_HQ(H2q,Q2q,resQ,drRes,dzRes,dpRes,ddRes,
!     &                    HQ2q,resHQ,imu)
!c
!c           ------------------------------------------------------------
!c           perform the gradient step in qp-space
!c           ------------------------------------------------------------
!            do itr=1,TRLen
!               HP2q(itr) = HQ2q(itr)-P2q(itr)
!               Q2q(itr)  =
!     &         Q2q(itr)  - stp_qp_P*(HQ2q(itr)-P2q(itr))/H2q(itr)
!            enddo
!c
!c           ------------------------------------------------------------
!c           evaluate errors (only if residual parts are neglected)
!c           ------------------------------------------------------------
!            if (isresact.eq.0) then
!               buffer   = error_P
!               var1     = var_oper(HQ2q,resHQ,imu)
!               var2     = var_oper( P2q, resP,imu)
!               var3     = 0.5D0*comm_QP(HQ2q,resHQ,P2q,resP,imu)
!               error_P  = sqrt(var_oper(HP2q,resHP,imu))
!               derror_P = error_P-buffer
!            endif
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           perform the gradient step in residual space
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!c
!c              ---------------------------------------------------------
!c              calculate norm of residual components of [H,Q]-Operator
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  HQnorm(ist)=
!     &                 wfovlp(resHQ(0,0,1,ist,1),resHQ(0,0,1,ist,1))+
!     &                 wfovlp(resHQ(0,0,1,ist,2),resHQ(0,0,1,ist,2))
!               enddo
!c
!c              ---------------------------------------------------------
!c              calculate errors in residual space
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  do ir=0,nr
!                     do iz=0,nz
!                        ResHP(iz,ir,1,ist,1) = ResHQ(iz,ir,1,ist,1)+
!     &                                         ResP(iz,ir,1,ist,1)
!                        ResHP(iz,ir,2,ist,1) = ResHQ(iz,ir,2,ist,1)+
!     &                                         ResP(iz,ir,2,ist,1)
!                        ResHP(iz,ir,1,ist,2) = ResHQ(iz,ir,1,ist,2)+
!     &                                         ResP(iz,ir,1,ist,2)
!                        ResHP(iz,ir,2,ist,2) = ResHQ(iz,ir,2,ist,2)+
!     &                                         ResP(iz,ir,2,ist,2)
!                     enddo
!                  enddo
!               enddo
!c
!               buffer   = error_P
!               error_P  = sqrt(var_oper(HP2q,resHP,imu))
!               derror_P = error_P-buffer
!c
!c              ---------------------------------------------------------
!c              apply damping
!c              ---------------------------------------------------------
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHP(0,0,1,1,1),ResHQ(0,0,1,1,1),
!     &                          (-1)**imu,+imu,stp_res_P,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHP(0,0,1,1,2),ResHQ(0,0,1,1,2),
!     &                          (-1)**imu,-imu,stp_res_P,nst(1)+nst(2))
!c
!c              ---------------------------------------------------------
!c              the gradient step
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  do ir=0,nr
!                     do iz=0,nz
!                      ResQ(iz,ir,1,ist,1) =
!     &                ResQ(iz,ir,1,ist,1) + resHQ(iz,ir,1,ist,1)
!                      ResQ(iz,ir,2,ist,1) =
!     &                ResQ(iz,ir,2,ist,1) + resHQ(iz,ir,2,ist,1)
!                      ResQ(iz,ir,1,ist,2) =
!     &                ResQ(iz,ir,1,ist,2) + resHQ(iz,ir,1,ist,2)
!                      ResQ(iz,ir,2,ist,2) =
!     &                ResQ(iz,ir,2,ist,2) + resHQ(iz,ir,2,ist,2)
!                     enddo
!                  enddo
!               enddo
!c
!c              ---------------------------------------------------------
!c              norm stabilization
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  norm=sqrt(Pnorm(ist)/HQnorm(ist))
!                  do ir=0,nr
!                     do iz=0,nz
!                        ResQ(iz,ir,1,ist,1) = ResQ(iz,ir,1,ist,1)*norm
!                        ResQ(iz,ir,2,ist,1) = ResQ(iz,ir,2,ist,1)*norm
!                        ResQ(iz,ir,1,ist,2) = ResQ(iz,ir,1,ist,2)*norm
!                        ResQ(iz,ir,2,ist,2) = ResQ(iz,ir,2,ist,2)*norm
!                     enddo
!                  enddo
!               enddo
!c
!c              ---------------------------------------------------------
!c              project out sp-contributions
!c              ---------------------------------------------------------
!               call project_qp(Q2q,ResQ,+1,imu,2)
!c               call project_Omega(ResQ,(-1)**imu,imu)
!c
!c           ------------------------------------------------------------
!c           END of gradient step in residual space
!c           ------------------------------------------------------------
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           Logging and convergency-check
!c           ------------------------------------------------------------
!            if (mod(iter,lograte).eq.0               .or.
!c     &          error_P.le.max_err_P*sqrt(abs(PVar)) .or.
!     &          iter.eq.max_citer ) then
!c
!               write(21,'(a,i5,4(f6.1))') 'SELFC    : ',
!     &              iter,
!     &              log10(MAX(abs(error_P                  ),1E-16)),
!     &              log10(MAX(abs(error_P/sqrt(abs(PVar))  ),1E-16)),
!     &              log10(MAX(abs(derror_P                 ),1E-16)),
!     &              log10(MAX(abs(1.0D0-var3**2/(var1*var2)),1E-16))
!c
!c               if (error_P.le.max_err_P*sqrt(abs(PVar))) goto 2103
!c
!            endif
!c
!            if (error_P.le.max_err_P) goto 2103
!c
!         enddo
!c        ===============================================================
!c        END of cranking iteration (find Q-Operator)
!c        ===============================================================
! 2103    continue
!c
!c
!c     ------------------------------------------------------------------
!c     ATDHF                                              (requires Beta)
!c     ------------------------------------------------------------------
!      else if (type.eq.4) then
!c
!c        initialization
!c        ---------------------------------------------------------------
!         BVar   = var_oper(B2q,ResB,imu)
!c
!c        find Q-operator by {\cal L}^2 Q = 4 BC Q
!c        ---------------------------------------------------------------
!         do itr=1,TRlen
!            Q2q(itr)= B2q(itr)/H2q(itr)**2 ! Q_ij = B_ij/(e_i + e_j)**2
!         enddo
!c
!#IF(FULLSPACE)
!c        inversion of {\cal L}_0 in residual space
!c        ---------------------------------------------------------------
!c        initialize using the B-Operator
!         if (isresact.eq.1) then
!c
!c           initialize using the B-Operator
!            do ist=1,nst(1)+nst(2)
!               do ir=0,nr
!                  do iz=0,nz
!                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
!                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
!                   ResQ(iz,ir,1,ist,2)=ResB(iz,ir,1,ist,2)
!                   ResQ(iz,ir,2,ist,2)=ResB(iz,ir,2,ist,2)
!                  enddo
!               enddo
!            enddo
!c
!c           calculate norm of residual componentes B-Operator
!            do ist=1,nst(1)+nst(2)
!               Bnorm(ist)=
!     &              wfovlp(ResB(0,0,1,ist,1),ResB(0,0,1,ist,1))
!     &             +wfovlp(ResB(0,0,1,ist,2),ResB(0,0,1,ist,2))
!            enddo
!c
!         else
!c
!c           set residual components to zero
!            call zero_oper(Q2q,ResQ,2)
!            call zero_oper(HP2q,ResHP,2)
!c
!         endif
!#ENDIF
!c
!c        ===============================================================
!c        the cranking iteration (find Q-Operator)
!c        ===============================================================
!         do iter=1,max_citer
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           calculate derivatives of resiual part of Q
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!               call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                              ddRes(0,0,1,1,1),
!     &                              (-1)**imu,+imu,nst(1)+nst(2))
!               call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                              ddRes(0,0,1,1,2),
!     &                              (-1)**imu,-imu,nst(1)+nst(2))
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           calculate action of Liouville-Operator
!c           ------------------------------------------------------------
!            call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
!     &                    HQ2q,ResHQ,imu)
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           calculate derivatives of resiual part of LQ \propto P
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!               call upd_derive_spinor(ResHQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                              ddRes(0,0,1,1,1),
!     &                              (-1)**imu,+imu,nst(1)+nst(2))
!               call upd_derive_spinor(ResHQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                              ddRes(0,0,1,1,2),
!     &                              (-1)**imu,-imu,nst(1)+nst(2))
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           calculate action of Liouville-Operator
!c           ------------------------------------------------------------
!            call comm_HP(H2q,HQ2q,ResHQ,drRes,dzRes,dpRes,ddRes,
!     &                    HP2q,ResHP,imu)
!c
!c           ------------------------------------------------------------
!c           calculate Omega^2
!c           ------------------------------------------------------------
!            OmgSq=comm_QP(Q2q,ResQ,HP2q,ResHP,imu)/
!     &            comm_QP(Q2q,ResQ, Q2q,ResQ, imu)
!c
!c           ------------------------------------------------------------
!c           perform the gradient step in qp-space
!c           ------------------------------------------------------------
!            do itr=1,TRLen
!               HQ2q(itr) = HP2q(itr)-OmgSq*Q2q(itr)
!               Q2q(itr)  = Q2q(itr) -stp_qp_Omg *
!     &                    (HP2q(itr)-OmgSq*Q2q(itr))/H2q(itr)**2
!            enddo
!c
!c           ------------------------------------------------------------
!c           evaluate errors (only if residual parts are neglected)
!c           ------------------------------------------------------------
!            if (isresact.eq.0) then
!               buffer   = error_Q
!               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
!               derror_Q = error_Q-buffer
!            endif
!c
!#IF(FULLSPACE)
!c           ------------------------------------------------------------
!c           perform the gradient step in residual space
!c           ------------------------------------------------------------
!            if (isresact.eq.1) then
!c
!c              ---------------------------------------------------------
!c              calculate norm of residual components of [H,Q]-Operator
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  HPnorm(ist)=
!     &                 wfovlp(resHP(0,0,1,ist,1),resHP(0,0,1,ist,1))+
!     &                 wfovlp(resHP(0,0,1,ist,2),resHP(0,0,1,ist,2))
!               enddo
!c
!c              ---------------------------------------------------------
!c              calculate errors in residual space
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  do ir=0,nr
!                     do iz=0,nz
!                        ResHQ(iz,ir,1,ist,1) = ResHP(iz,ir,1,ist,1)-
!     &                                         ResB(iz,ir,1,ist,1)
!                        ResHQ(iz,ir,2,ist,1) = ResHP(iz,ir,2,ist,1)-
!     &                                         ResB(iz,ir,2,ist,1)
!                        ResHQ(iz,ir,1,ist,2) = ResHP(iz,ir,1,ist,2)-
!     &                                         ResB(iz,ir,1,ist,2)
!                        ResHQ(iz,ir,2,ist,2) = ResHP(iz,ir,2,ist,2)-
!     &                                         ResB(iz,ir,2,ist,2)
!                     enddo
!                  enddo
!               enddo
!c
!               buffer   = error_Q
!               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
!               derror_Q = error_Q-buffer
!c
!c              ---------------------------------------------------------
!c              apply damping
!c              ---------------------------------------------------------
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHQ(0,0,1,1,1),ResHP(0,0,1,1,1),
!     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHP(0,0,1,1,1),ResHQ(0,0,1,1,1),
!     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHQ(0,0,1,1,2),ResHP(0,0,1,1,2),
!     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
!               call damp_spinor(LRDampr,LRDampz,
!     &                          ResHP(0,0,1,1,2),ResHQ(0,0,1,1,2),
!     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))

!c
!c              ---------------------------------------------------------
!c              the gradient step
!c              ---------------------------------------------------------
!               do ist=1,nst(1)+nst(2)
!                  do ir=0,nr
!                     do iz=0,nz
!                      ResQ(iz,ir,1,ist,1) =
!     &                ResQ(iz,ir,1,ist,1) - ResHQ(iz,ir,1,ist,1)
!                      ResQ(iz,ir,2,ist,1) =
!     &                ResQ(iz,ir,2,ist,1) - ResHQ(iz,ir,2,ist,1)
!                      ResQ(iz,ir,1,ist,2) =
!     &                ResQ(iz,ir,1,ist,2) - ResHQ(iz,ir,1,ist,2)
!                      ResQ(iz,ir,2,ist,2) =
!     &                ResQ(iz,ir,2,ist,2) - ResHQ(iz,ir,2,ist,2)
!                     enddo
!                  enddo
!               enddo
!c
!c              ---------------------------------------------------------
!c              project out sp-contributions
!c              ---------------------------------------------------------
!               call project_qp(Q2q,ResQ,+1,imu,2)
!c
!c           ------------------------------------------------------------
!c           END of gradient step in residual space
!c           ------------------------------------------------------------
!            endif
!#ENDIF
!c
!c           ------------------------------------------------------------
!c           normalize Q to 2*<QQ> = 1
!c           ------------------------------------------------------------
!            norm=1.0/sqrt(comm_QP(Q2q,ResQ, Q2q,ResQ, imu))
!            do itr=1,TRlen
!               Q2q(itr)=Q2q(itr)*norm
!            enddo
!c
!c
!c           ------------------------------------------------------------
!c           Logging and convergency-check
!c           ------------------------------------------------------------
!c
!            if (mod(iter,lograte).eq.0               .or.
!     &          error_Q.le.max_err_Q*sqrt(abs(QVar)) .or.
!     &          iter.eq.max_citer ) then
!c
!               write(21,'(a,i5,1(1pg16.4),2(e16.3))') 'ATDHF     : ',
!     &              iter,error_Q,error_Q/sqrt(abs(BVar)),
!     &              derror_Q
!c
!               if (error_Q.le.max_err_Q*sqrt(abs(QVar))) goto 2104
!               if (abs(derror_Q).lt.1.0E-10) goto 2104
!c
!            endif
!c
!         enddo
!c        ===============================================================
!c        END of cranking iteration (find Q-Operator)
!c        ===============================================================
! 2104    continue
!c
!#IF(FULLSPACE)
!c        ------------------------------------------------------------
!c        calculate derivatives of residual part of Q
!c         ------------------------------------------------------------
!         if (isresact.eq.1) then
!            call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                           dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                           ddRes(0,0,1,1,1),
!     &                           (-1)**imu,+imu,nst(1)+nst(2))
!            call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                            dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                            ddRes(0,0,1,1,2),
!     &                            (-1)**imu,-imu,nst(1)+nst(2))
!         endif
!#ENDIF
!c
!c        ------------------------------------------------------------
!c        calculate action of Liouville-Operator to obtain P-Operator
!c        ------------------------------------------------------------
!         call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
!     &                P2q,ResP,imu)
!c
!      endif
!c
!      write(88,'(a)') '--- cranking finished  ---'
!c
!c     ==================================================================
!c     end of the cranking methods
!c     ==================================================================
!c
!c     Normalization:
!c
!c     normalize P to Beta  and Q to P
!c
!c       1) GCM-formulas do not need any explicit normalization.
!c       2) for inglis and sc-cranking normalization is defined by
!c
!c             [H,Q] = -i 2B P  --> [H, Q/(2B)] = -iP
!c
!c          and restored by Q -> 2BQ
!c
!c       3) CHF-cranking and ATDHF-cranking need explicit normalisation
!c
!c     in case of fullspace calculations (isresact=1) Q,P could be
!c     normalized explicitly in any case.
!c
!c
!c     explicit normalization of P and Q to B
!c     ==================================================================
!      if (FULLSPACE.or.
!     & ((type.eq.0.or.type.eq.4.or.type.eq.5).or.(isresact.eq.1))) then
!c
!         B=-1.0D0   ! just set a flag -- this is not the real mass
!c
!#IF(FULLSPACE)
!c
!c       normalization of P to B       <[B,P]> = i
!c       ----------------------------------------------------------------
!         norm=1.0D0/comm_QP(B2q,resB,P2q,resP,imu)
!         do itr=1,TRlen
!            P2q(itr)= P2q(itr)*norm
!         enddo
!c
!         do ist=1,nst(1)+nst(2)
!            do ir=0,nr
!               do iz=0,nz
!                  resP(iz,ir,1,ist,1) = resP(iz,ir,1,ist,1)*norm
!                  resP(iz,ir,2,ist,1) = resP(iz,ir,2,ist,1)*norm
!                  resP(iz,ir,1,ist,2) = resP(iz,ir,1,ist,2)*norm
!                  resP(iz,ir,2,ist,2) = resP(iz,ir,2,ist,2)*norm
!               enddo
!            enddo
!         enddo
!c
!c        normalization of Q to P      <[Q,P]> = i
!c        ---------------------------------------------------------------
!         norm = 1.0D0/comm_QP(Q2q,resQ,P2q,resP,imu)
!         do itr=1,TRlen
!            Q2q(itr)= Q2q(itr)*norm
!         enddo
!c
!         do ist=1,nst(1)+nst(2)
!            do ir=0,nr
!               do iz=0,nz
!                  resQ(iz,ir,1,ist,1) = resQ(iz,ir,1,ist,1)*norm
!                  resQ(iz,ir,2,ist,1) = resQ(iz,ir,2,ist,1)*norm
!                  resQ(iz,ir,1,ist,2) = resQ(iz,ir,1,ist,2)*norm
!                  resQ(iz,ir,2,ist,2) = resQ(iz,ir,2,ist,2)*norm
!               enddo
!            enddo
!         enddo
!#ENDIF
!c
!c
!c     implicit normalization of Q to P:  i[H,Q] = 2BP
!c     ==================================================================
!      else if (type.eq.2.or.type.eq.3) then
!c
!         B    = 0.5/comm_QP(Q2q,resQ,P2q,resP,imu)
!         do itr=1,TRlen
!            Q2q(itr)= 2.0D0*B*Q2q(itr)
!         enddo
!c
!c     the GCM-mass with fullspace precission
!c     ==================================================================
!      else if (type.eq.1) then
!c


!c
!      endif
!c
!      write(88,'(a)') '--- normalization finished  ---'
!c
!c     Calculation of widths, masses and stiffness:
!c     ==================================================================
!c     calculate <Q^2>,<P^2>,<[Beta,P]>,<[Q,P]>
!c     ------------------------------------------------------------------
!      QVar   = var_oper(Q2q,resQ,imu)
!#IF(FULLSPACE)
!      Nrm0   = comm_QP(B2q,resB,P2q,resP,imu)
!      Nrm    = comm_QP(Q2q,resQ,P2q,resP,imu)
!      PVar   = var_oper(P2q,resP,imu)
!#ELSE
!      call comm_and_var_QP(BVar,PVar,Nrm0,CPbeta(0),
!     &                     (-1)**imu,imu,0)
!#ENDIF
!c
!c     calculate LQ = [H,Q] + H_Q
!c     and        B = 1/2*<[Q,LQ]]>
!c     ------------------------------------------------------------------
!c
!#IF(FULLSPACE)
!c     prepare derivatives of residual part of Q
!      call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                       dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                       ddRes(0,0,1,1,1),
!     &                       (-1)**imu,+imu,nst(1)+nst(2))
!      call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                       dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                       ddRes(0,0,1,1,2),
!     &                       (-1)**imu,-imu,nst(1)+nst(2))
!#ENDIF
!c
!c     if mass has not been calculated before
!      if (B.eq.-1.0D0) then
!c
!c        calculate action of Liouville-Operator
!         call comm_HQ(H2q,Q2q,resQ,drRes,dzRes,dpRes,ddRes,
!     &                HQ2q,resHQ,imu)
!         B = 0.5D0*comm_QP(Q2q,resQ,HQ2q,resHQ,imu)
!c
!      endif
!c
!c     calculate  LP = [H,P] + H_P
!c     and        C  = 1/2*<[P,LP]> = -1/2*<[LP,P]>
!c     ------------------------------------------------------------------
!#IF(FULLSPACE)
!      call upd_derive_spinor(ResP(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                       dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                       ddRes(0,0,1,1,1),
!     &                       (-1)**imu,+imu,nst(1)+nst(2))
!      call upd_derive_spinor(ResP(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                       dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                       ddRes(0,0,1,1,2),
!     &                       (-1)**imu,-imu,nst(1)+nst(2))
!#ENDIF
!c
!      call comm_HP(H2q,P2q,ResP,drRes,dzRes,dpRes,ddRes,
!     &             HP2q,resHP,imu)
!c
!      C = -0.5D0*comm_QP(HP2q,ResHP,P2q,ResP,imu)
!c
!c     calculate L^2Q    = [H,LQ] + H_{LQ}
!c     and       Omega^2 = <[P,L^2Q]>
!c     ------------------------------------------------------------------
!#IF(FULLSPACE)
!      call upd_derive_spinor(ResHQ(0,0,1,1,1),drRes(0,0,1,1,1),
!     &                       dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
!     &                       ddRes(0,0,1,1,1),
!     &                       (-1)**imu,+imu,nst(1)+nst(2))
!      call upd_derive_spinor(ResHQ(0,0,1,1,2),drRes(0,0,1,1,2),
!     &                       dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
!     &                       ddRes(0,0,1,1,2),
!     &                       (-1)**imu,-imu,nst(1)+nst(2))
!#ENDIF
!c
!      call comm_HP(H2q,HQ2q,ResHQ,drRes,dzRes,dpRes,ddRes,
!     &             HP2q,resHP,imu)
!c
!      OmgSq = comm_QP(HP2q,ResHP,P2q,ResP,imu)
!c
!c     calculate particle number mass B_N
!c     ------------------------------------------------------------------
!      call comm_NQ(Q2q,ResQ,HQ2q,ResHQ)
!      BN = 0.5D0*comm_QP(Q2q,ResQ,HQ2q,ResHQ,imu)
!c
!      write(88,'(a)') '--- evaluation finished  ---'
!c
!c     Logging
!c     ------------------------------------------------------------------
!      if (type.eq.0) then
!         write(21,'(a,8(1pe16.4))') 'CHF_SC_VAL: ',
!     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
!      elseif (type.eq.1) then
!         write(21,'(a,8(1pe16.4))') 'GCM_VAL   : ',
!     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
!      elseif (type.eq.2) then
!         write(21,'(a,8(1pe16.4))') 'INGLIS_VAL: ',
!     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
!      elseif (type.eq.3) then
!         write(21,'(a,8(1pe16.4))') 'SELFC_VAL : ',
!     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
!      elseif (type.eq.4) then
!         write(21,'(a,8(1pe16.4))') 'ATDHF_VAL : ',
!     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
!      endif
!c
!c     ------------------------------------------------------------------
!c
!      return
!      end


!=======================================================================

!     upd_path()

!       updates the collective path which is stored in trpsi.
!       The subroutne manages the storage of the last three sets
!       of wavefunctions and occupation amplitudes from the static
!       iteration and also some relevant observables.

!=======================================================================

SUBROUTINE upd_path()

!     ------------------------------------------------------------------

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!INCLUDE 'fields.inc'
!INCLUDE 'observ.inc' ! beta2,beta3,atnum,rms
!INCLUDE 'constr.inc' ! lqcon

!     ------------------------------------------------------------------
INTEGER :: ist,ir,iz,i
REAL(8) :: equasi,deltai,edifi


!     shift wavefunctions in path storage CPpsi (-1,0,+1) -> (0,+1,+1)
!     and copy current static wavefunction to last path storage index
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  DO ir=0,nr
    DO iz=0,nz
      cppsi(iz,ir,1,ist,-1) = cppsi(iz,ir,1,ist, 0)
      cppsi(iz,ir,2,ist,-1) = cppsi(iz,ir,2,ist, 0)
      cppsi(iz,ir,1,ist, 0) = cppsi(iz,ir,1,ist,+1)
      cppsi(iz,ir,2,ist, 0) = cppsi(iz,ir,2,ist,+1)
      cppsi(iz,ir,1,ist,+1) = psi(iz,ir,1,ist)
      cppsi(iz,ir,2,ist,+1) = psi(iz,ir,2,ist)
    END DO
  END DO
END DO

!     shift mean-fields in path storage   (-1,0,+1) -> (0,+1,+1)
!     and copy current mean-fields to last path storage index
!     ------------------------------------------------------------------
DO ir=0,nr
  DO iz=0,nz
    cprho(iz,ir,1,-1)     = cprho(iz,ir,1, 0)
    cprho(iz,ir,2,-1)     = cprho(iz,ir,2, 0)
    cprho(iz,ir,1, 0)     = cprho(iz,ir,1,+1)
    cprho(iz,ir,2, 0)     = cprho(iz,ir,2,+1)
    cprho(iz,ir,1,+1)     = rho(iz,ir,1)
    cprho(iz,ir,2,+1)     = rho(iz,ir,2)
    
    cpkin(iz,ir,1,-1)     = cpkin(iz,ir,1, 0)
    cpkin(iz,ir,2,-1)     = cpkin(iz,ir,2, 0)
    cpkin(iz,ir,1, 0)     = cpkin(iz,ir,1,+1)
    cpkin(iz,ir,2, 0)     = cpkin(iz,ir,2,+1)
    cpkin(iz,ir,1,+1)     = tau(iz,ir,1)
    cpkin(iz,ir,2,+1)     = tau(iz,ir,2)
    
    cpsocr(iz,ir,1,-1)    = cpsocr(iz,ir,1, 0)
    cpsocr(iz,ir,2,-1)    = cpsocr(iz,ir,2, 0)
    cpsocr(iz,ir,1, 0)    = cpsocr(iz,ir,1,+1)
    cpsocr(iz,ir,2, 0)    = cpsocr(iz,ir,2,+1)
    cpsocr(iz,ir,1,+1)    = rholsr(iz,ir,1)
    cpsocr(iz,ir,2,+1)    = rholsr(iz,ir,2)
    
    cpsocz(iz,ir,1,-1)    = cpsocz(iz,ir,1, 0)
    cpsocz(iz,ir,2,-1)    = cpsocz(iz,ir,2, 0)
    cpsocz(iz,ir,1, 0)    = cpsocz(iz,ir,1,+1)
    cpsocz(iz,ir,2, 0)    = cpsocz(iz,ir,2,+1)
    cpsocz(iz,ir,1,+1)    = rholsz(iz,ir,1)
    cpsocz(iz,ir,2,+1)    = rholsz(iz,ir,2)
    
    cplocpot(iz,ir,1,-1) = cplocpot(iz,ir,1, 0)
    cplocpot(iz,ir,2,-1) = cplocpot(iz,ir,2, 0)
    cplocpot(iz,ir,1, 0) = cplocpot(iz,ir,1,+1)
    cplocpot(iz,ir,2, 0) = cplocpot(iz,ir,2,+1)
    cplocpot(iz,ir,1,+1) = pot(iz,ir,1)
    cplocpot(iz,ir,2,+1) = pot(iz,ir,2)
    
    cpeffmas(iz,ir,1,-1) = cpeffmas(iz,ir,1, 0)
    cpeffmas(iz,ir,2,-1) = cpeffmas(iz,ir,2, 0)
    cpeffmas(iz,ir,1, 0) = cpeffmas(iz,ir,1,+1)
    cpeffmas(iz,ir,2, 0) = cpeffmas(iz,ir,2,+1)
    cpeffmas(iz,ir,1,+1) = bmass(iz,ir,1)
    cpeffmas(iz,ir,2,+1) = bmass(iz,ir,2)
    
    cpsocpotr(iz,ir,1,-1) = cpsocpotr(iz,ir,1, 0)
    cpsocpotr(iz,ir,2,-1) = cpsocpotr(iz,ir,2, 0)
    cpsocpotr(iz,ir,1, 0) = cpsocpotr(iz,ir,1,+1)
    cpsocpotr(iz,ir,2, 0) = cpsocpotr(iz,ir,2,+1)
    cpsocpotr(iz,ir,1,+1) = potls(iz,ir,1,1)
    cpsocpotr(iz,ir,2,+1) = potls(iz,ir,2,1)
    
    cpsocpotz(iz,ir,1,-1) = cpsocpotz(iz,ir,1, 0)
    cpsocpotz(iz,ir,2,-1) = cpsocpotz(iz,ir,2, 0)
    cpsocpotz(iz,ir,1, 0) = cpsocpotz(iz,ir,1,+1)
    cpsocpotz(iz,ir,2, 0) = cpsocpotz(iz,ir,2,+1)
    cpsocpotz(iz,ir,1,+1) = potls(iz,ir,1,2)
    cpsocpotz(iz,ir,2,+1) = potls(iz,ir,2,2)
    
    cpcoupot(iz,ir,-1)    = cpcoupot(iz,ir, 0)
    cpcoupot(iz,ir, 0)    = cpcoupot(iz,ir,+1)
    cpcoupot(iz,ir,+1)    = coulpt(iz,ir)
    
    cpconstr(iz,ir,-1)    = cpconstr(iz,ir, 0)
    cpconstr(iz,ir, 0)    = cpconstr(iz,ir,+1)
    cpconstr(iz,ir,+1)    = qfac2(iz,ir)           ! 'linresp' only for L=2
  END DO
END DO

!     get sp-observables from stationary iteration and shift path
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  cpvsqf(ist,-1)   = cpvsqf(ist, 0)
  cpvsqf(ist, 0)   = cpvsqf(ist,+1)
  cpvsqf(ist,+1)   = vsqf(ist)
  
  cpdeltaf(ist,-1) = cpdeltaf(ist, 0)
  cpdeltaf(ist, 0) = cpdeltaf(ist,+1)
  cpdeltaf(ist,+1) = deltaf(ist)
  
  cpwtpair(ist,-1) = cpwtpair(ist, 0)
  cpwtpair(ist, 0) = cpwtpair(ist,+1)
  cpwtpair(ist,+1) = wtpair(ist)
  
  cpspenrg(ist,-1) = cpspenrg(ist, 0)
  cpspenrg(ist, 0) = cpspenrg(ist,+1)
  cpspenrg(ist,+1) = spenrg(ist)
  
  cpqpenrg(ist,-1) = cpqpenrg(ist, 0)
  cpqpenrg(ist, 0) = cpqpenrg(ist,+1)
  
  deltai = deltaf(ist)*wtpair(ist)
  edifi  = spenrg(ist)-fer(icharg(ist))
  equasi = SQRT(edifi*edifi+deltai*deltai)
  
  IF (ipair /= -1) THEN
    cpqpenrg(ist,+1) = equasi-fer(icharg(ist))
  ELSE
    cpqpenrg(ist,+1) = -spenrg(ist)
  END IF
  
END DO

!     get global observables from stationary iteration and shift path
!     ------------------------------------------------------------------
cpenergy(-1)   = cpenergy( 0)
cpenergy( 0)   = cpenergy(+1)
cpenergy(+1)   = energy

cpenergyc(-1)   = cpenergyc( 0)
cpenergyc( 0)   = cpenergyc(+1)
cpenergyc(+1)   = qlamb*beta2

cplambda(-1)   = cplambda( 0)
cplambda( 0)   = cplambda(+1)
cplambda(+1)   = qlamb

cpzcm(-1)      = cpzcm( 0)
cpzcm( 0)      = cpzcm(+1)
cpzcm(+1)      = zcm

cpzpemod(-1)   = cpzpemod( 0)
cpzpemod( 0)   = cpzpemod(+1)
cpzpemod(+1)   = zpemod

cpradrms(-1)   = cpradrms( 0)      ! r      (RMS radius)
cpradrms( 0)   = cpradrms(+1)
cpradrms(+1)   = rchrg2

cpraddif(-1)   = cpraddif( 0)      ! R      (diff. radius)
cpraddif( 0)   = cpraddif(+1)
cpraddif(+1)   = rdmsc

cpradsrf(-1)   = cpradsrf( 0)      ! sigma  (surface. thickness)
cpradsrf( 0)   = cpradsrf(+1)
cpradsrf(+1)   = surfc

cpfer(1,-1)    = cpfer(1, 0)
cpfer(1, 0)    = cpfer(1,+1)
cpfer(1,+1)    = fer(1)

cpfer(2,-1)    = cpfer(2, 0)
cpfer(2, 0)    = cpfer(2,+1)
cpfer(2,+1)    = fer(2)

cppairfc(1,-1) = cppairfc(1, 0)
cppairfc(1, 0) = cppairfc(1, 1)
cppairfc(1, 1) = pairfc(1)

cppairfc(2,-1) = cppairfc(2, 0)
cppairfc(2, 0) = cppairfc(2, 1)
cppairfc(2, 1) = pairfc(2)

cpatnum(-1) = cpatnum( 0)
cpatnum( 0) = cpatnum(+1)
cpatnum(+1) = atnum

IF (lqcon == 2) THEN
  cpbeta(-1) = cpbeta( 0)
  cpbeta( 0) = cpbeta(+1)
  cpbeta(+1) = beta2
ELSE IF (lqcon == 3) THEN
  cpbeta(-1) = cpbeta( 0)
  cpbeta( 0) = cpbeta(+1)
  cpbeta(+1) = beta3
END IF

!     set current length of collective path
!     ------------------------------------------------------------------
cplen = cplen + 1

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,3f8.4)') 'UPD_PATH()    : collective path updated --'//  &
      ' beta= ',(cpbeta(i),i=-1,1)
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE upd_path


!=======================================================================

!     upd_derivatives()

!       updates the derivatives of the sp-wavefunctions and Mean-Fields
!       of the central element of the collective path.

!=======================================================================

SUBROUTINE upd_derivatives()

!     ------------------------------------------------------------------

USE paramr
USE constr
USE dynpar
USE observ, ONLY: beta2,beta3,atnum,rms
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'  !h2meff
!INCLUDE 'linresp.inc'

!INCLUDE 'observ.inc' ! beta2,beta3,atnum,rms
!INCLUDE 'constr.inc' ! lqcon

!      double precision HPsi(0:MZD,0:MRD,MSD,MPSID)
REAL(8) :: energy,var
INTEGER :: ist

!     ------------------------------------------------------------------

ipos= 0 ! position on path    ==> check with ini_matrix_p
!     calculate derivatives of single particle wavefunctions
!     ------------------------------------------------------------------
CALL upd_derive_spinor(cppsi(0,0,1,1,ipos), cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1), +1,0,nst(1)+nst(2))

!     derivatives requiere effmass to vanish on bound of box
!     ------------------------------------------------------------------
DO ir=0,nr
  DO iz=0,nz
    cpeffmas(iz,ir,1,0) = cpeffmas(iz,ir,1,0) - h2meff(1)
    cpeffmas(iz,ir,2,0) = cpeffmas(iz,ir,2,0) - h2meff(2)
  END DO
END DO

CALL nablaz_scalar(cpeffmas(0:nz,0:nr,1,0),cpdzeffmas(0:nz,0:nr,1),1,1,one)
CALL nablar_scalar(cpeffmas(0:nz,0:nr,1,0),cpdreffmas(0:nz,0:nr,1),0,1,one)
CALL nablaz_scalar(cpeffmas(0:nz,0:nr,2,0),cpdzeffmas(0:nz,0:nr,2),1,1,one)
CALL nablar_scalar(cpeffmas(0:nz,0:nr,2,0),cpdreffmas(0:nz,0:nr,2),0,1,one)

DO ir=0,nr
  DO iz=0,nz
    cpeffmas(iz,ir,1,0) = cpeffmas(iz,ir,1,0) + h2meff(1)
    cpeffmas(iz,ir,2,0) = cpeffmas(iz,ir,2,0) + h2meff(2)
  END DO
END DO

!     test action of Mean-Field-Hamiltonian on states
!     ------------------------------------------------------------------
!      if (1.eq.0) then
!      call act_mfh_even(CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),
!     &                  CPdzpsi(0,0,1,1),CPdppsi(0,0,1,1),
!     &                  CPddpsi(0,0,1,1),HPsi(0,0,1,1),
!     &                  CPlocpot(0,0,1,0),                    1,
!     &                  CPeffmas(0,0,1,0),CPdreffmas(0,0,1),
!     &                  CPdzeffmas(0,0,1),CPdzeffmas(0,0,1),  1,
!     &                  CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
!     &                  CPsocpotr(0,0,1,0),                   1,
!     &                  CPconstr(0,0,0),CPlambda(0),          1,
!     &                  CPspenrg(1,0),                        0,
!     &                  0,1.0D0,nst(1)+nst(2))
!c
!      do ist=1,nst(1)+nst(2)
!         energy = wfovlp(CPpsi(0,0,1,ist,0),HPsi(0,0,1,ist))
!         var    = wfovlp(HPsi(0,0,1,ist),HPsi(0,0,1,ist))-energy**2
!         write(88,'(a,3(1pg16.4),3i3)') 'PKtest:',energy,var,
!     &        CPspenrg(ist,0),mpsi(ist),nparz(1,ist),icharg(ist)
!      enddo
!      endif

!     write derivatives of single particle wavefunctions
!     ------------------------------------------------------------------
!      do ist=1,nst(1)+nst(2)
!         call wrt_2D_spin(CPpsi(0,0,1,ist,0) ,'psi_______',1,'tst')
!         call wrt_2D_spin(CPdrpsi(0,0,1,ist) ,'drpsi_____',0,'tst')
!         call wrt_2D_spin(CPdzpsi(0,0,1,ist) ,'dzpsi_____',0,'tst')
!         call wrt_2D_spin(CPdppsi(0,0,1,ist) ,'dppsi_____',0,'tst')
!         call wrt_2D_spin(CPddpsi(0,0,1,ist) ,'ddpsi_____',0,'tst')
!      enddo


!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a)') 'UPD_DERIVES() : derivatives of '// 'wavefunctions updated'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE upd_derivatives


!=======================================================================

!     ini_matrix_P(P2q,PVar,il,imu,iiso)

!       initialize the generator of axial deformation (mu=0) or
!       a rotation arround the y-axis (mu=1). In case of rotations
!       the operator is scaled by 1/(sqrt(3)*beta) to be equivalent
!       with the cartesean P_{1} (only at the axial configurations.

!       furthermore this function calculates

!       2qp-matrixelements    P^{2q}_ji = eta^(-)_ji <j|P|i>
!       variance              <P^2>

!=======================================================================

SUBROUTINE ini_matrix_p(p2q,pvar,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
USE wavefr
USE observ, ONLY: beta2
USE coulomb, ONLY: ylm
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: p2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! 2qp-matrixelements





INTEGER :: trzpar
REAL(8) :: pfac,eta,ovlp
REAL(8) :: z,r
INTEGER :: iz,ir
REAL(8) :: lastv,lastu,nextv,nextu,rcntv,rcntu
REAL(8) :: ui,uj,vi,vj
INTEGER :: itr,ist,jst,jstmin,jstmax

LOGICAL :: oldcrank=d_oldcrank
LOGICAL,PARAMETER :: toutside=.false.,ttestam=.FALSE.
!DATA             oldcrank/d_oldcrank/

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     some factors
!     ------------------------------------------------------------------
IF (imu == 0) THEN
  IF (oldcrank) THEN
    pfac     = 1.0D0/(cpbeta( 0)-cpbeta(-1))
  ELSE
    pfac     = 1.0D0/(cpbeta(+1)-cpbeta(-1))
  END IF
ELSE IF (imu == 1) THEN
  IF (ABS(cpbeta(0)) > 0.0D0) THEN
    pfac=1.0D0 !/(SQRT(3.0D0)*cpbeta(0))
  ELSE
    pfac=0.0D0
  END IF
ELSE
  WRITE(21,'(a,i2)') 'ERROR: only mu=0 and mu=1 possible -- not ',imu
  STOP 'in linresp'
END IF

!     set variance and 2q-matrixelements to zero
!     ------------------------------------------------------------------
pvar = 0D0
pvaroff=0D0
DO itr2=1,trntrtot
  p2q(itr2)=0D0
END DO

!     loop over all sp-states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
!        calculate action of P operator
!        ===============================================================
  
!        P_0 (generator of translation along the axial path, asym. def.)
!        --------------------------------------------------------------
  IF (imu == 0.AND.oldcrank) THEN
    
    DO ir=0,nr
      DO iz=0,nz
        lrbuff1(iz,ir,1,1)=pfac* (cppsi(iz,ir,1,ist, 0)-cppsi(iz,ir,1,ist,-1))
        lrbuff1(iz,ir,2,1)=pfac* (cppsi(iz,ir,2,ist, 0)-cppsi(iz,ir,2,ist,-1))
        lrbuff1(iz,ir,1,2)=0.0D0
        lrbuff1(iz,ir,2,2)=0.0D0
      END DO
    END DO
    
!        P_0 (generator of translation along the axial path, sym. def.)
!        --------------------------------------------------------------
  ELSE IF (imu == 0) THEN
    
    DO ir=0,nr
      DO iz=0,nz
        lrbuff1(iz,ir,1,1)=pfac* (cppsi(iz,ir,1,ist,+1)-cppsi(iz,ir,1,ist,-1))
        lrbuff1(iz,ir,2,1)=pfac* (cppsi(iz,ir,2,ist,+1)-cppsi(iz,ir,2,ist,-1))
        lrbuff1(iz,ir,1,2)=0.0D0
        lrbuff1(iz,ir,2,2)=0.0D0
      END DO
    END DO
    
!        P_1 (generator of rotation,  J_y = L_y + S_y)
!        --------------------------------------------------------------
  ELSE IF (imu == 1) THEN

    ipos= 0    ! position on path
  IF(tangmomsub) THEN
    CALL angmomop(cppsi(0,0,1,ist,ipos),lrbuff1(:,:,:,1),&
                  lrbuff1(:,:,:,2),nparz(1,ist),mpsi(ist))  
    lrbuff1=lrbuff1*0.5D0*pfac     ! to comply with normalization here
  ELSE
    CALL nablaz_spin(cppsi(0:nz,0:nr,1:2,ist,ipos),psi1(0:nz,0:nr,1:2),nparz(1,ist),1,1D0)
    CALL nablar_spin(cppsi(0:nz,0:nr,1:2,ist,ipos),psi2(0:nz,0:nr,1:2),   mpsi(ist),1,1D0)
!
! test derivatives:
!
    DO is=1,nsd
    mm=mpsi(ist)+is-1
    DO ir=0,nr
      r = (ir+sr)*dr
      IF(ir>0) THEN
        rinvact = mm/r
      ELSE
        rinvact = 0D0
      END IF
      DO iz=0,nz
        psi3(iz,ir,is) = rinvact *cppsi(iz,ir,is,ist,ipos)
      END DO
    END DO
    END DO
!    ovlpn = wfovlp(cpdppsi(0,0,1,ist),cpdppsi(0,0,1,ist))
!    ovlpo = wfovlp(psi3,psi3)
!    WRITE(6,'(a,2i4,6(1pg12.4))') 'DERIVATIVES:',ist,mpsi(ist), &
!      wfovlp(cpdzpsi(0,0,1,ist),psi1)/ &
!      1D0-&
!      SQRT(ABS(wfovlp(cpdzpsi(0,0,1,ist),cpdzpsi(0,0,1,ist))*wfovlp(psi1,psi1))), &
!      wfovlp(cpdrpsi(0,0,1,ist),psi2)/ &
!      1D0-&
!      SQRT(ABS(wfovlp(cpdrpsi(0,0,1,ist),cpdrpsi(0,0,1,ist))*wfovlp(psi2,psi2))),  &
!      1D0-wfovlp(cpdppsi(0,0,1,ist),psi3)/SQRT(ovlpn*ovlpo),ovlpn,ovlpo

! test path
!    ovlp = wfovlp(psi(0,0,1,ist),psi(0,0,1,ist))
!    WRITE(6,'(a,2i4,6(1pg15.7))') 'PATH:',ist,mpsi(ist), &
!      1D0-wfovlp(cppsi(0,0,1,ist,-1),psi(0,0,1,ist))/ &
!      SQRT(wfovlp(cppsi(0,0,1,ist,-1),cppsi(0,0,1,ist,-1))*ovlp),&
!      1D0-wfovlp(cppsi(0,0,1,ist,0),psi(0,0,1,ist))/ &
!      SQRT(wfovlp(cppsi(0,0,1,ist,0),cppsi(0,0,1,ist,0))*ovlp),&
!      1D0-wfovlp(cppsi(0,0,1,ist,1),psi(0,0,1,ist))/ &
!      SQRT(wfovlp(cppsi(0,0,1,ist,1),cppsi(0,0,1,ist,1))*ovlp)
    
    DO is=1,nsd
      mm=mpsi(ist)+is-1
      ir=0
      r = (ir+sr)*dr
      DO iz=0,nz
        z = (iz+sz-icm)*dz - cpzcm(0)
        IF(mm==1) THEN        
!                 ! upward laddering components of J_y
          lrbuff1(iz,ir,is,1) = 0D0 ! -  0.5D0*pfac*cppsi(iz,ir,2,ist,ipos)
        
!                 ! downward laddering components of J_y
          lrbuff1(iz,ir,is,2) =  -pfac*z*cpdrpsi(iz,ir,is,ist)
!                              +  0.5D0*pfac*cppsi(iz,ir,1,ist,ipos)
        ELSE
!                 ! upward laddering components of J_y
          lrbuff1(iz,ir,is,1) =  0.5D0*pfac*z*cpdrpsi(iz,ir,is,ist) ! &
!                              - 0.5D0*pfac*cppsi(iz,ir,2,ist,ipos)
        
!                 ! downward laddering components of J_y
          lrbuff1(iz,ir,is,2) =  -0.5D0*pfac*z*cpdrpsi(iz,ir,is,ist)  
!                               + 0.5D0*pfac*cppsi(iz,ir,1,ist,ipos)
        END IF
      END DO
    END DO
!                        spin flip ?

    DO ir=1,nr
      r = (ir+sr)*dr
      DO iz=0,nz
        z = (iz+sz-icm)*dz - cpzcm(ipos)
        
!                 ! upward laddering components of J_y
        lrbuff1(iz,ir,1,1) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,1,ist)  &
            +z*cpdppsi(iz,ir,1,ist) +r*cpdzpsi(iz,ir,1,ist)  &
            -  cppsi(iz,ir,2,ist,ipos))   ! - ??
        lrbuff1(iz,ir,2,1) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,2,ist)  &
            +z*cpdppsi(iz,ir,2,ist) +r*cpdzpsi(iz,ir,2,ist))
        
!                 ! downward laddering components of J_y
        lrbuff1(iz,ir,1,2) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,1,ist)  &
            -z*cpdppsi(iz,ir,1,ist) +r*cpdzpsi(iz,ir,1,ist))
        lrbuff1(iz,ir,2,2) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,2,ist)  &
            -z*cpdppsi(iz,ir,2,ist) +r*cpdzpsi(iz,ir,2,ist)  &
            +  cppsi(iz,ir,1,ist,ipos))
        
      END DO
    END DO
  END IF
! test J_i

  IF(tangmomsub) THEN
    CALL angmomop(cppsi(0,0,1,ist,ipos),psip,psin,nparz(1,ist),mpsi(ist))  
  ELSE
  DO is=1,nsd
    mm        = mpsi(ist)+is-1
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
        psip(iz,ir,is)=z*(psi2(iz,ir,is)-psi3(iz,ir,is)) -r*psi1(iz,ir,is)
        IF(is==1) psip(iz,ir,1) = psip(iz,ir,1)+psi(iz,ir,2,ist)
        psin(iz,ir,is)=-z*(psi2(iz,ir,is)+psi3(iz,ir,is))+r*psi1(iz,ir,is)
        IF(is==2) psin(iz,ir,2) = psin(iz,ir,2)+psi(iz,ir,1,ist)
      END DO
    END DO
  END DO
  END IF
  
  IF(ttestam) THEN
!  lrbuff1(:,:,:,1)=-lrbuff1(:,:,:,1)    !  ???
    ovlpn1 = wfovlp(lrbuff1(0:nz,0:nr,1:2,1),lrbuff1(0:nz,0:nr,1:2,1))
    ovlpn2 = wfovlp(lrbuff1(0:nz,0:nr,1:2,2),lrbuff1(0:nz,0:nr,1:2,2))
    ovlpo1 = wfovlp(psip(0:nz,0:nr,1:2),psip(0:nz,0:nr,1:2))
    ovlpo2 = wfovlp(psin(0:nz,0:nr,1:2),psin(0:nz,0:nr,1:2))
    ovlpp1 = wfovlp(ampsi(0:nz,0:nr,1:2,ist,1),ampsi(0:nr,0:nz,1:2,ist,1))
    ovlpp2 = wfovlp(ampsi(0:nz,0:nr,1:2,ist,2),ampsi(0:nz,0:nr,1:2,ist,2))
    WRITE(6,'(a,2i4,6(1pg12.4))') 'ANG.MOM:    ',ist,mpsi(ist), &
      1D0-wfovlp(lrbuff1(0:nz,0:nr,1:2,1),psip(0:nz,0:nr,1:2))/SQRT(ovlpn1*ovlpo1),&
      1D0-wfovlp(lrbuff1(0:nz,0:nr,1:2,2),psin(0:nz,0:nr,1:2))/SQRT(ovlpn2*ovlpo2),&
      ovlpp1,ovlpo1,ovlpp2,ovlpo2
    WRITE(6,'(a,2i4,6(1pg12.4))') '            ',ist,mpsi(ist), &
      1D0-wfovlp(ampsi(0:nz,0:nr,1:2,ist,1),psip(0:nz,0:nr,1:2))/SQRT(ovlpp1*ovlpo1),&
      1D0-wfovlp(ampsi(0:nz,0:nr,1:2,ist,2),psin(0:nz,0:nr,1:2))/SQRT(ovlpp2*ovlpo2),&
      ovlpp1,ovlpo1,ovlpp2,ovlpo2
    psip=ampsi(:,:,:,ist,1)
    psin=ampsi(:,:,:,ist,2)
    WRITE(6,'(a,2i4,6(1pg12.4))') '            ',ist,mpsi(ist), &
      1D0-wfovlp(lrbuff1(0:nz,0:nr,1:2,1),psip(0:nz,0:nr,1:2))/SQRT(ovlpn1*ovlpp1),&
      1D0-wfovlp(lrbuff1(0:nz,0:nr,1:2,2),psin(0:nz,0:nr,1:2))/SQRT(ovlpn2*ovlpp2),&
      ovlpn1,ovlpp1,ovlpn2,ovlpp2
    END IF
    
  END IF
  
!        adjust lower boundary in case of no z-symmetry
!        ---------------------------------------------------------------
  IF (isym == 0) THEN
    DO ir=0,nr
      lrbuff1(0,ir,1,1) = 0.0D0
      lrbuff1(0,ir,2,1) = 0.0D0
      lrbuff1(0,ir,1,2) = 0.0D0
      lrbuff1(0,ir,2,2) = 0.0D0
    END DO
  END IF
  
!        evaluate 2qp-matrixelements
!        ------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    itp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    eta  = vj*ui-uj*vi
    
!    ovlp1=wfovlp(cppsi(0,0,1,jst,0),lrbuff1(0,0,1,1))
!    ovlp2=wfovlp(cppsi(0,0,1,jst,0),lrbuff1(0,0,1,2))
    ovlp3=0D0
!    IF(ABS(mpsi(ist)-mpsi(jst))==1) THEN
    IF (itp == +1) THEN
!      IF(mpsi(ist)+1==mpsi(jst)) THEN
      IF (msgn >= 0) THEN
        ovlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,ipos),lrbuff1(0:nz,0:nr,1:2,1))
        ovlp3=wfovlp(cppsi(0:nz,0:nr,1:2,jst,ipos),psip(0:nz,0:nr,1:2))/2D0
        IF(toutside) CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,ipos),lrbuff1(0:nz,0:nr,1:2,1), ovlp)
        isw=3
      ELSE 
        ovlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,ipos),lrbuff1(0:nz,0:nr,1:2,2))
        ovlp3=wfovlp(cppsi(0:nz,0:nr,1:2,jst,ipos),psin(0:nz,0:nr,1:2))/2D0
        IF(toutside) CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,ipos),lrbuff1(0:nz,0:nr,1:2,2), ovlp)
        isw=2
      END IF
!    ELSE IF(mpsi(ist)==0 .AND. mpsi(jst)==0) THEN
    ELSE
      ovlp=wfovlptr(cppsi(0:nz,0:nr,1:2,jst,ipos),lrbuff1(0:nz,0:nr,1:2,2))
      ovlp3=wfovlptr(cppsi(0:nz,0:nr,1:2,jst,ipos),psin(0:nz,0:nr,1:2))/2D0
      IF(toutside) CALL wfprojtr(cppsi(0:nz,0:nr,1:2,jst,ipos),lrbuff1(0:nz,0:nr,1:2,2), ovlp)
      isw=1
    END IF
    
!           diagonal elements describe change of occupation amplitudes
    IF (ist == jst.AND.itp == +1) THEN
      
      IF (oldcrank) THEN
        lastv   = sqrtb(cpvsqf(ist,-1))
        lastu   = sqrtb(one-cpvsqf(ist,-1))
        nextv   = sqrtb(cpvsqf(ist,0))
        nextu   = sqrtb(one-cpvsqf(ist,0))
        rcntv   = 0.5D0*(lastv+nextv)
        rcntu   = 0.5D0*(lastu+nextu)
        IF (rcntv >= 0.5D0) THEN
          ovlp =  -(nextu-lastu)/rcntv*pfac
        ELSE
          ovlp =  +(nextv-lastv)/rcntu*pfac
        END IF
      ELSE
        lastv   = sqrtb(cpvsqf(ist,-1))
        lastu   = sqrtb(one-cpvsqf(ist,-1))
        rcntv   = sqrtb(cpvsqf(ist,0))
        rcntu   = sqrtb(one-cpvsqf(ist,0))
        nextv   = sqrtb(cpvsqf(ist,+1))
        nextu   = sqrtb(one-cpvsqf(ist,+1))
        IF (rcntv >= 0.5D0) THEN
          ovlp =  -(nextu-lastu)/rcntv*pfac
        ELSE
          ovlp =  +(nextv-lastv)/rcntu*pfac
        END IF
      END IF
      IF(imu==0) THEN            ! ??
        p2q(itr2) =               ovlp
        pvar      = pvar + (      ovlp)**2
      END IF
      
!           straight forward for off-diagonal elements
    ELSE
      
      p2q(itr2) =         eta * ovlp
      pvar      = pvar + (eta * ovlp)**2
      IF(ttestam .AND. eta**2>1D-5 .AND. imu==1) &
        WRITE(992,'(a,2f9.6,5i4,4(1pg13.5))') &
        'LINRES:',betampsi,cpbeta(0),ist,jst,mpsi(ist),mpsi(jst),isw,&
          ovlp,(eta*ovlp*2D0)**2,ovlp3
!      IF(eta**2>1D-5 .AND. imu==1) WRITE(992,'(a,3f9.6,6i4,4(1pg13.5))') &
!        'LINRES:',beta2,cpbeta(0),b0q(0,0),ist,jst,&
!          itp,msgn,mpsi(ist),mpsi(jst),ovlp,(eta*ovlp*2D0)**2
!    (eta*ovlp1*2D0)**2,(eta*ovlp2*2D0)**2,(eta*ovlp3*2D0)**2
    END IF
    
  END DO
  
!        accumulate contributions from outside qp-space
!        ---------------------------------------------------------------
  IF(toutside) THEN
  pvaradd = 2.0D0*cpvsqf(ist,0)*(  &
      wfovlp(lrbuff1(0:nz,0:nr,1:2,1),lrbuff1(0:nz,0:nr,1:2,1))  &
      + wfovlp(lrbuff1(0:nz,0:nr,1:2,2),lrbuff1(0:nz,0:nr,1:2,2)) )
  pvar = pvar + pvaradd
  pvaroff = pvaroff + pvaradd
  END IF  
END DO
IF(ttestam .AND. imu==1) THEN
  WRITE(992,'(a,4(1pg13.5))') 'total:',pvar
  CALL flush(992)
END IF


IF(imu==1 .AND. toutside) THEN
  WRITE(276,'(6(1pg13.5))') b0q(0,0),cpbeta(0),pvar*3.0D0*cpbeta(0)**2,&
    pvaroff*3.0D0*cpbeta(0)**2,&
    pexp1*3.0D0*cpbeta(0)**2,pexp2*3.0D0*cpbeta(0)**2
  CALL flush(276)
END IF
!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'INI_MATRIX_P(): P-operator initialized'//  &
      '  (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE ini_matrix_p


!=======================================================================

!     crank_test_resp(B2q,Q2q,P2q,H2q,il,imu)

!       checks the responses

!=======================================================================

SUBROUTINE crank_test_resp(b2q,q2q,p2q,h2q,mass,stiff,qvar,pvar, il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: p2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                     :: mass
REAL(8), INTENT(IN OUT)                     :: stiff
REAL(8), INTENT(IN OUT)                     :: qvar
REAL(8), INTENT(IN OUT)                     :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements




!REAL(8) :: mass                   ! collective mass B=1/(2M)
!REAL(8) :: qvar                   ! <Q^2>-<Q>^2
!REAL(8) :: pvar                   ! <P^2>
!REAL(8) :: stiff                  ! stiffness C

INTEGER :: itr
INTEGER :: iter
REAL(8) :: norm,norm0,omgsq
REAL(8) :: hresp2q(1:mtdtot)
REAL(8) :: hq2q(1:mtdtot)
REAL(8) :: d2q(1:mtdtot)
REAL(8) :: maxresp

INTEGER :: maxiter
REAL(8) :: stepq
REAL(8) :: m1,m2,m3,m4,m5,m6
REAL(8) :: s1,s2,s3,s4,s5,s6

maxiter = 200
stepq   = 0.8D0

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp A'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp B'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',lambdamax,')'
  STOP 'in linresp C'
END IF

!     initialize with Q-operator in inglis approximation
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  p2q(itr) = b2q(itr)
  q2q(itr) = b2q(itr)
END DO

!     check the Q-response
!     ------------------------------------------------------------------
CALL upd_response_qtest(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu, 1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
m1  =0.0D0
DO itr=1,trntrtot
  m1     = m1  + 2.0D0*q2q(itr)*hresp2q(itr)
END DO

CALL upd_response_qtest(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu, 0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
m2  =0.0D0
DO itr=1,trntrtot
  m2     = m2  + 2.0D0*q2q(itr)*hresp2q(itr)
END DO

CALL upd_response_qtest(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu, 0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
m3  =0.0D0
DO itr=1,trntrtot
  m3     = m3  + 2.0D0*q2q(itr)*hresp2q(itr)
END DO

CALL upd_response_qtest(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu, 0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
m4  =0.0D0
DO itr=1,trntrtot
  m4     = m4  + 2.0D0*q2q(itr)*hresp2q(itr)
END DO

CALL upd_response_qtest(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu, 0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
m5  =0.0D0
DO itr=1,trntrtot
  m5     = m5  + 2.0D0*q2q(itr)*hresp2q(itr)
END DO

CALL upd_response_qtest(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu, 0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
m6  =0.0D0
DO itr=1,trntrtot
  m6     = m6  + 2.0D0*q2q(itr)*hresp2q(itr)
END DO

!      write(89,'(i2,7g16.4)') imu,CPbeta(0),m1,m2,m3,m4,m5,m6
!      return

!     check the P-response
!     ------------------------------------------------------------------
CALL upd_response_ptest(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu,  &
    1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
s1  =0.0D0
DO itr=1,trntrtot
  s1     = s1  + 2.0D0*p2q(itr)*hresp2q(itr)
END DO

CALL upd_response_ptest(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu,  &
    0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
s2  =0.0D0
DO itr=1,trntrtot
  s2     = s2  + 2.0D0*p2q(itr)*hresp2q(itr)
END DO

CALL upd_response_ptest(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu,  &
    0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
s3  =0.0D0
DO itr=1,trntrtot
  s3     = s3  + 2.0D0*p2q(itr)*hresp2q(itr)
END DO

CALL upd_response_ptest(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu,  &
    0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
s4  =0.0D0
DO itr=1,trntrtot
  s4     = s4  + 2.0D0*p2q(itr)*hresp2q(itr)
END DO

CALL upd_response_ptest(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu,  &
    0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
s5  =0.0D0
DO itr=1,trntrtot
  s5     = s5  + 2.0D0*p2q(itr)*hresp2q(itr)
END DO

CALL upd_response_ptest(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu,  &
    0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
s6  =0.0D0
DO itr=1,trntrtot
  s6     = s6  + 2.0D0*p2q(itr)*hresp2q(itr)
END DO

WRITE(89,'(i2,7g16.4)') imu,cpbeta(0),s1,s2,s3,s4,s5,s6

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_test_resp


!=======================================================================

!     crank_inglis(B2q,Q2q,P2q,H2q,il,imu)

!       cranking in inglis approximation

!               i[H,(Q/2B)] = P            1/(2B) = <[(Q/2B),P]>

!=======================================================================

SUBROUTINE crank_inglis(b2q,q2q,p2q,h2q,mass,stiff,qvar,pvar, il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(IN)             :: p2q(1:mtdtot)
REAL(8), INTENT(IN)             :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                     :: mass
REAL(8), INTENT(OUT)                        :: stiff
REAL(8), INTENT(OUT)                        :: qvar
REAL(8), INTENT(OUT)                        :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements




!REAL(8) :: mass                   ! collective mass B=1/(2M)
!REAL(8) :: qvar                   ! <Q^2>-<Q>^2
!REAL(8) :: pvar                   ! <P^2>
!REAL(8) :: stiff                  ! stiffness C

INTEGER :: itr
REAL(8) :: norm,norm0,omgsq
REAL(8) :: hp2q(1:mtdtot)
REAL(8) :: hresp2q(1:mtdtot)

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     evaluate the P-response
!     ------------------------------------------------------------------
CALL upd_response_even(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
DO itr=1,trntrtot
  hp2q(itr) = h2q(itr)*p2q(itr) + hresp2q(itr)
END DO

!     set variables to zero
!     ------------------------------------------------------------------
mass  = 0.0D0
norm  = 0.0D0
qvar  = 0.0D0
pvar  = 0.0D0
stiff = 0.0D0

!     loop over all transitions
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  
  q2q(itr) =               p2q(itr)/h2q(itr)       ! 1/2B    Q
  mass     = mass  + 2.0D0*p2q(itr)**2/h2q(itr)    ! 1/2B
  qvar     = qvar  +       p2q(itr)**2/h2q(itr)**2 ! 1/(4BB) <Q^2>
  pvar     = pvar  +       p2q(itr)**2             !         <P^2>
  stiff    = stiff + 2.0D0*p2q(itr)*hp2q(itr  )    ! 2*C
  
END DO

!     normalize
!     ------------------------------------------------------------------
norm  = 1.0D0/mass
mass  = 0.5D0/mass
stiff = 0.5D0*stiff

qvar  =  norm**2*qvar
pvar  =          pvar

omgsq = mass*stiff

DO itr=1,trntrtot
  q2q(itr) = q2q(itr)*norm
END DO

!     check the normalisation
norm =0.0D0
norm0=0.0D0
DO itr=1,trntrtot
  norm0    = norm0 + 2.0D0*b2q(itr)*p2q(itr)
  norm     = norm  + 2.0D0*q2q(itr)*p2q(itr)
END DO

!     logging
WRITE(21,'(a,7(1pe16.4))') 'INGLIS_VAL: ',  &
    mass,stiff,pvar,qvar,omgsq,norm,norm0

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_inglis

!=======================================================================

!     crank_GCM(B2q,Q2q,P2q,H2q,il,imu)

!       cranking in static GCM approximation

!               Q = P/lambda

!=======================================================================

SUBROUTINE crank_gcm(b2q,q2q,p2q,h2q,mass,stiff,qvar,pvar, il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(IN)             :: p2q(1:mtdtot)
REAL(8), INTENT(IN)             :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                     :: mass
REAL(8), INTENT(OUT)                        :: stiff
REAL(8), INTENT(OUT)                        :: qvar
REAL(8), INTENT(OUT)                        :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements




!REAL(8) :: mass                   ! collective mass B=1/(2M)
!REAL(8) :: qvar                   ! <Q^2>-<Q>^2
!REAL(8) :: pvar                   ! <P^2>
!REAL(8) :: stiff                  ! stiffness C

INTEGER :: itr
REAL(8) :: norm,norm0,omgsq
REAL(8) :: hresp2q(1:mtdtot)
REAL(8) :: hq2q(1:mtdtot)
REAL(8) :: hp2q(1:mtdtot)


!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     calculate the Q-Operator in static GCM approximation
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  q2q(itr) = p2q(itr)/(2.0D0*pvar)
END DO

!     evaluate the Q-response
!     ------------------------------------------------------------------
CALL upd_response_odd(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
DO itr=1,trntrtot
  hq2q(itr) = h2q(itr)*q2q(itr) + hresp2q(itr)
END DO

!     loop over all transitions
!     ------------------------------------------------------------------
mass=0.0D0
qvar=0.0D0
DO itr=1,trntrtot
  mass     = mass  + 2.0D0*hq2q(itr)*q2q(itr) ! 2B
  qvar     = qvar  +        q2q(itr)*q2q(itr)
END DO
mass = 0.5D0*mass

!     evaluate the P-response
!     ------------------------------------------------------------------
CALL upd_response_even(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
DO itr=1,trntrtot
  hp2q(itr) = h2q(itr)*p2q(itr) + hresp2q(itr)
END DO

!     loop over all transitions
!     ------------------------------------------------------------------
stiff=0.0D0
pvar =0.0D0
norm =0.0D0
norm0=0.0D0
DO itr=1,trntrtot
  pvar      = pvar   +        p2q(itr)*p2q(itr)
  norm      = norm   + 2.0D0* q2q(itr)*p2q(itr)
  norm0     = norm0  + 2.0D0* b2q(itr)*p2q(itr)
  stiff     = stiff  + 2.0D0*hp2q(itr)*p2q(itr) ! 2C
END DO
stiff = 0.5D0*stiff

omgsq = mass*stiff

DO itr=1,trntrtot
  q2q(itr) = q2q(itr)*norm
END DO

!     logging
WRITE(21,'(a,7(1pe16.4))') 'GCM_VAL   : ',  &
    mass,stiff,pvar,qvar,omgsq,norm,norm0

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_gcm


!=======================================================================

!     crank_scaling(B2q,Q2q,P2q,H2q,il,imu)

!       cranking in scaling approximation

!          [H,[H,Q]] = 4 B C Beta      B=1/2*<[Q,[H,Q]]>
!                                      C=1/2*<[P,[H,P]]>

!=======================================================================

SUBROUTINE crank_scaling(b2q,q2q,p2q,h2q,mass,stiff, qvar,pvar,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: p2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: mass
REAL(8), INTENT(OUT)            :: stiff
REAL(8), INTENT(OUT)            :: qvar
REAL(8), INTENT(OUT)            :: pvar
INTEGER, INTENT(IN OUT)                  :: il
INTEGER, INTENT(IN OUT)                  :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements






INTEGER :: ist,jst
INTEGER :: trzpar
INTEGER :: itr
INTEGER :: iso,iiso
REAL(8) :: norm,norm0,normp,normq,omgsq

!     ------------------------------------------------------------------

!     set variables to zero
!     ------------------------------------------------------------------
mass  = 0.0D0
norm  = 0.0D0
norm0 = 0.0D0
qvar  = 0.0D0
pvar  = 0.0D0

!     loop over all transitions
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  
  p2q(itr) =               b2q(itr)/h2q(itr)       ! = 1/(2C)   P
  q2q(itr) =               b2q(itr)/h2q(itr)**2    ! = 1/(4BC) Q
  norm0    = norm0 + 2.0D0*b2q(itr)**2/h2q(itr)    ! = 1/(2C)
  pvar     = pvar  +       b2q(itr)**2/h2q(itr)**2 ! = 1/(4CC) <P^2>
  mass     = mass  + 2.0D0*b2q(itr)**2/h2q(itr)**3 ! = 1/(8BCC)
  qvar     = qvar  +       b2q(itr)**2/h2q(itr)**4 ! = 1/(16BBCC) Q
  
END DO

!     normalize
!     ------------------------------------------------------------------
normq = norm0/mass   ! factor to restore norm of Q operator
normp = 1.0D0/norm0  ! factor to restore norm of P operator

qvar  = normq**2       * qvar
pvar  = normp**2       * pvar
mass  = 0.5D0*norm0**2 / mass
stiff = 0.5D0*normp
omgsq = mass*stiff

DO itr=1,trntrtot
  p2q(itr) = p2q(itr)*normp
  q2q(itr) = q2q(itr)*normq
END DO

!     check the normalisation
norm =0.0D0
norm0=0.0D0
DO itr=1,trntrtot
  norm0    = norm0 + 2.0D0*b2q(itr)*p2q(itr)
  norm     = norm  + 2.0D0*q2q(itr)*p2q(itr)
END DO

!     logging
WRITE(21,'(a,7(1pe16.4))') 'SCALE_VAL : ',  &
    mass,stiff,pvar,qvar,omgsq,norm,norm0

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_scaling

!=======================================================================

!     crank_selfc(B2q,Q2q,P2q,H2q,il,imu)

!       ATDHF cranking with a given P-operator

!               i[H,(Q/2B)] = P            1/(2B) = <[(Q/2B),P]>

!=======================================================================

SUBROUTINE crank_selfc(b2q,q2q,p2q,h2q,mass,stiff,qvar,pvar, il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(IN)             :: p2q(1:mtdtot)
REAL(8), INTENT(IN)             :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                     :: mass
REAL(8), INTENT(OUT)                        :: stiff
REAL(8), INTENT(OUT)                        :: qvar
REAL(8), INTENT(OUT)                        :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements




!REAL(8) :: mass                   ! collective mass B=1/(2M)
!REAL(8) :: qvar                   ! <Q^2>-<Q>^2
!REAL(8) :: pvar                   ! <P^2>
!REAL(8) :: stiff                  ! stiffness C

INTEGER :: itr
INTEGER :: iter
REAL(8) :: norm,norm0,omgsq
REAL(8) :: hresp2q(1:mtdtot)
REAL(8) :: hq2q(1:mtdtot)
REAL(8) :: hp2q(1:mtdtot)
REAL(8) :: d2q(1:mtdtot)
REAL(8) :: maxresp

INTEGER :: maxiter
REAL(8) :: stepq

maxiter = 200
stepq   = 0.8D0

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     evaluate the P-response
!     ------------------------------------------------------------------
CALL upd_response_even(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
DO itr=1,trntrtot
  hp2q(itr) = h2q(itr)*p2q(itr) + hresp2q(itr)
END DO

!     initialize with Q-operator in inglis approximation
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  q2q(itr) = p2q(itr)/h2q(itr)
END DO

!     solve i[H,(Q/2B)] + h_{Q,resp} = P by iteration
!     ------------------------------------------------------------------
DO iter=1,maxiter
  
  CALL upd_response_odd(q2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
      redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
      respnmasr,respnmasp,respnmasz,il,imu)
  
  CALL upd_hresp_odd(hresp2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
      redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
      respnmasr,respnmasp,respnmasz,il,imu)
  
!        evaluate the Liouville operator i[H,Q] + h_{Q,resp}
  DO itr=1,trntrtot
    hq2q(itr) = h2q(itr)*q2q(itr) + hresp2q(itr)
  END DO
  
!        perform the gradient step
  DO itr=1,trntrtot
    d2q(itr) = hq2q(itr) - p2q(itr)
    q2q(itr) = q2q(itr)  - stepq*(hq2q(itr)-p2q(itr))/(h2q(itr))
  END DO
  
!        evaluate mass and errors
  mass  =0.0D0
  qvar  =0.0D0
  buff1 =0.0D0
  buff2 =0.0D0
  buff3 =0.0D0
  DO itr=1,trntrtot
    mass     = mass  + 2.0D0*q2q(itr)* p2q(itr) ! 1/2B
    qvar     = qvar  +       q2q(itr)**2        ! 1/(4BB) <Q^2>
    buff1    = buff1 +      hq2q(itr)* p2q(itr)
    buff2    = buff2 +      hq2q(itr)*hq2q(itr)
    buff3    = buff3 +       p2q(itr)* p2q(itr)
  END DO
  errp=1.0D0-buff1**2/(MAX(buff2*buff3,1D-10))
  IF(buff1**2/(MAX(buff2*buff3,1D-10)) < 1D-6) errp=1D-10
  
!        logging
  WRITE(22,'(a,2i4,3(1pg12.4))') 'SELFC: iter,imu,mass,qvar,errP:',  &
      iter,imu,mass,qvar,errp
  
!        check for convergency
  IF(ABS(errp) < 1E-6) GO TO 4000
  
END DO
4000 CONTINUE

!     set variables to zero
!     ------------------------------------------------------------------
mass  = 0.0D0
norm  = 0.0D0
qvar  = 0.0D0
pvar  = 0.0D0
stiff = 0.0D0

!     loop over all transitions
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  
  mass     = mass  + 2.0D0*q2q(itr)*p2q(itr)       ! 1/2B
  qvar     = qvar  +       q2q(itr)**2             ! 1/(4BB) <Q^2>
  pvar     = pvar  +       p2q(itr)**2             !         <P^2>
  stiff    = stiff + 2.0D0*p2q(itr)*hp2q(itr)      !  2C
  
END DO

!     normalize
!     ------------------------------------------------------------------
norm  = 1.0D0/mass
mass  = 0.5D0/mass
stiff = 0.5D0*stiff

qvar  =  norm**2*qvar
pvar  =          pvar

omgsq = mass*stiff

DO itr=1,trntrtot
  q2q(itr) = q2q(itr)*norm
END DO

!     check the normalisation
norm =0.0D0
norm0=0.0D0
DO itr=1,trntrtot
  norm0    = norm0 + 2.0D0*b2q(itr)*p2q(itr)
  norm     = norm  + 2.0D0*q2q(itr)*p2q(itr)
END DO

!     logging
WRITE(21,'(a,7(1pe16.4))') 'SELFC_VAL : ',  &
    mass,stiff,pvar,qvar,omgsq,norm,norm0

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_selfc


!=======================================================================

!     crank_CHF(B2q,Q2q,P2q,H2q,il,imu)

!       CHF + ATDHF-cranking for a given Beta-operator

!               [H,[H,(Q/4Omega^2)]] = Beta

!=======================================================================

SUBROUTINE crank_chf(b2q,q2q,p2q,h2q,mass,stiff,qvar,pvar, il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: p2q(1:mtdtot)
REAL(8), INTENT(IN)             :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                     :: mass
REAL(8), INTENT(OUT)                        :: stiff
REAL(8), INTENT(OUT)                        :: qvar
REAL(8), INTENT(OUT)                        :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements




!REAL(8) :: mass                   ! collective mass B=1/(2M)
!REAL(8) :: qvar                   ! <Q^2>-<Q>^2
!REAL(8) :: pvar                   ! <P^2>
!REAL(8) :: stiff                  ! stiffness C

INTEGER :: itr
INTEGER :: iter
REAL(8) :: norm,norm0,omgsq,normp,normq
REAL(8) :: hresp2q(1:mtdtot)
REAL(8) :: hq2q(1:mtdtot)
REAL(8) :: hp2q(1:mtdtot)
REAL(8) :: d2q(1:mtdtot)

INTEGER :: maxiter
REAL(8) :: stepq

maxiter = 400
stepq   = 0.8D0

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     initialize with Q-operator with saved matrixelements
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  q2q(itr) = q2qchf(itr,imu)
END DO

!     solve i[H,[H,Q] + h_{Q,resp}] + h_{P,resp} = Beta by iteration
!     ------------------------------------------------------------------
DO iter=1,maxiter
  
  CALL upd_response_odd(q2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
      redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
      respnmasr,respnmasp,respnmasz,il,imu)
  
  CALL upd_hresp_odd(hresp2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
      redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
      respnmasr,respnmasp,respnmasz,il,imu)
  
!        evaluate the Liouville operator i[H,Q] + h_{Q,resp}
  DO itr=1,trntrtot
    hq2q(itr) = h2q(itr)*q2q(itr) + hresp2q(itr)
  END DO
  
  CALL upd_response_even(hq2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
      redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)
  
  CALL upd_hresp_even(hresp2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
      resocpotr,resocpotp,resocpotz,il,imu)
  
!        evaluate the Liouville operator -i[H,P] + h_{P,resp}
  DO itr=1,trntrtot
    hp2q(itr) = h2q(itr)*hq2q(itr) + hresp2q(itr)
  END DO
  
!        perform the gradient step
  DO itr=1,trntrtot
    d2q(itr) = hp2q(itr) - b2q(itr)
    q2q(itr) = q2q(itr)  - stepq*(hp2q(itr)-b2q(itr))/ (h2q(itr)**2)
  END DO
  
!        evaluate mass and errors
  norm  =0.0D0
  norm0 =0.0D0
  qvar  =0.0D0
  pvar  =0.0D0
  buff1 =0.0D0
  buff2 =0.0D0
  buff3 =0.0D0
  DO itr=1,trntrtot
    norm     = norm  + 2.0D0*b2q(itr)*hq2q(itr) ! 1/(2C)
    norm0    = norm0 + 2.0D0*q2q(itr)*hq2q(itr) ! 1/(8BCC)
    qvar     = qvar  +       q2q(itr)**2        ! 1/(16BBCC)<Q^2>
    pvar     = pvar  +      hq2q(itr)**2        ! 1/(4CC)   <P^2>
    buff1    = buff1 +      hp2q(itr)* b2q(itr)
    buff2    = buff2 +      hp2q(itr)*hp2q(itr)
    buff3    = buff3 +       b2q(itr)* b2q(itr)
  END DO
  mass  = 0.5D0 *norm**2/norm0
  stiff = 0.5D0         /norm
  omgsq = 0.25D0*norm   /norm0
  pvar  = pvar          /norm**2
  qvar  = qvar*  (norm  /norm0)**2
  errp  = 1.0D0-buff1**2/(buff2*buff3)
  
!        logging
  WRITE(22,'(a,i4,6(1pg12.3))') 'CHF  : ',  &
      iter,mass,stiff,omgsq,qvar,pvar,errp
  
!        check for convergency
  IF(ABS(errp) < 1E-6) GO TO 4001
  
END DO
4001 CONTINUE

!     normalize
!     ------------------------------------------------------------------
normp  =    1.0D0/norm
normq  =     norm/norm0

DO itr=1,trntrtot
  q2q(itr) =  q2q(itr)*normq
  p2q(itr) = hq2q(itr)*normp
END DO

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
!     ------------------------------------------------------------------
CALL upd_response_odd(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

DO itr=1,trntrtot
  hq2q(itr) = h2q(itr)*q2q(itr) + hresp2q(itr)
END DO

!     evaluate the Liouville operator -i[H,P] + h_{P,resp}
!     ------------------------------------------------------------------
CALL upd_response_even(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

DO itr=1,trntrtot
  hp2q(itr) = h2q(itr)*p2q(itr) + hresp2q(itr)
END DO

!     check the norm and evaluate the variances, mass and stiffness
!     ------------------------------------------------------------------
mass =0.0D0
stiff=0.0D0
norm =0.0D0
norm0=0.0D0
qvar =0.0D0
pvar =0.0D0
DO itr=1,trntrtot
  norm0    = norm0 + 2.0D0*b2q(itr)*p2q(itr)
  norm     = norm  + 2.0D0*q2q(itr)*p2q(itr)
  mass     = mass  +       q2q(itr)*hq2q(itr)
  stiff    = stiff +       p2q(itr)*hp2q(itr)
  qvar     = qvar  +       q2q(itr)*q2q(itr)
  pvar     = pvar  +       p2q(itr)*p2q(itr)
END DO
omgsq = mass*stiff

!     logging
WRITE(21,'(a,7(1pe16.4))') 'CHF_VAL   : ',  &
    mass,stiff,pvar,qvar,omgsq,norm,norm0

!     save matrixelements of the Q- and P-Operator
!     ------------------------------------------------------------------
DO itr=1,trntrtot
  q2qchf(itr,imu) = q2q(itr)
  p2qchf(itr,imu) = p2q(itr)
END DO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_chf


!=======================================================================

!     crank_ATDHF(B2q,Q2q,P2q,H2q,il,imu)

!       full ATDHF

!               [H,[H,Q]] = 4 Omega^2 Q

!=======================================================================

SUBROUTINE crank_atdhf(b2q,q2q,p2q,h2q,mass,stiff,qvar,pvar, il,imu,pathmode)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: q2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: p2q(1:mtdtot)
REAL(8), INTENT(IN)             :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                     :: mass
REAL(8), INTENT(OUT)                        :: stiff
REAL(8), INTENT(OUT)                        :: qvar
REAL(8), INTENT(OUT)                        :: pvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
INTEGER, INTENT(IN)                         :: pathmode
!REAL(8) :: ! e_i + e_j
!REAL(8) :: ! 2qp-matrixelements




!REAL(8) :: mass                   ! collective mass B=1/(2M)
!REAL(8) :: qvar                   ! <Q^2>-<Q>^2
!REAL(8) :: pvar                   ! <P^2>
!REAL(8) :: stiff                  ! stiffness C
!INTEGER :: pathmode

INTEGER :: itr
INTEGER :: iter
REAL(8) :: norm,norm0,omgsq,normp,normq
REAL(8) :: hresp2q(1:mtdtot)
REAL(8) :: hq2q(1:mtdtot)
REAL(8) :: hp2q(1:mtdtot)
REAL(8) :: d2q(1:mtdtot)
REAL(8) :: dmp2q(1:mtdtot)
REAL(8) :: c2q(1:mtdtot)

INTEGER :: maxiter
REAL(8) :: stepq

maxiter = 1000
stepq   = 0.1D0

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
ELSE IF (pathmode /= 1.AND.pathmode /= 0) THEN
  WRITE(21,'(a,i2,a,i2,a)')  &
      'ERROR: pathmode must be 0 (ATHDF) or 1 (CHF) -- not', pathmode
  STOP 'in linresp'
END IF

!     initialize with Constraint-operator and Q-operator
!     ------------------------------------------------------------------
IF (pathmode == 0) THEN
  DO itr=1,trntrtot      ! ATDHF
    q2q(itr) = q2qatdhf(itr,imu)
    c2q(itr) = q2qatdhf(itr,imu)
  END DO
ELSE                      ! CHF
  DO itr=1,trntrtot
    q2q(itr) = q2qchf(itr,imu)
    c2q(itr) = b2q(itr)
  END DO
END IF

!     solve i[H,[H,Q] + h_{Q,resp}] + h_{P,resp} = Beta by iteration
!     ------------------------------------------------------------------
DO iter=1,maxiter
  
  CALL upd_response_odd(q2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
      redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
      respnmasr,respnmasp,respnmasz,il,imu)
  
  CALL upd_hresp_odd(hresp2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
      redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
      respnmasr,respnmasp,respnmasz,il,imu)
  
!        evaluate the Liouville operator i[H,Q] + h_{Q,resp}
  mass=0.0D0   ! B
  qvar=0.0D0
  DO itr=1,trntrtot
    hq2q(itr) = h2q(itr)*q2q(itr)    + hresp2q(itr)
    mass      = mass + h2q(itr)*q2q(itr)*b2q(itr) +  &
        hresp2q(itr)     *b2q(itr)
    qvar      = qvar +   q2q(itr)*q2q(itr)
  END DO
  
  CALL upd_response_even(hq2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
      redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)
  
  CALL upd_hresp_even(hresp2q,  &
      cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
      cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
      resocpotr,resocpotp,resocpotz,il,imu)
  
!        evaluate the Liouville operator -i[H,P] + h_{P,resp}
  omgsq=0.0D0   ! 4*Omega^2 = 4BC
  pvar =0.0D0
  DO itr=1,trntrtot
    hp2q(itr) = h2q(itr)*hq2q(itr) + hresp2q(itr)
    omgsq     = omgsq + h2q(itr)*hq2q(itr)*hq2q(itr) +  &
        hresp2q(itr)      *hq2q(itr)
    pvar      = pvar  +  hq2q(itr)*hq2q(itr)
  END DO
  omgsq = omgsq/mass   ! 4*Omega^2
  pvar  = pvar/(4.0D0*mass**2)
  stiff = omgsq/(4.0D0*mass)
  
!        perform the gradient step
  DO itr=1,trntrtot
    d2q(itr) = hp2q(itr) - omgsq*c2q(itr)
    q2q(itr) = q2q(itr)  - stepq*(hp2q(itr)-omgsq*c2q(itr))/  &
        (h2q(itr)**2+1.0D0)
!     &                                    Dmp2q(itr)
  END DO
  
!        evaluate variances, norms and errors
  norm  =0.0D0
  buff1 =0.0D0
  buff2 =0.0D0
  buff3 =0.0D0
  DO itr=1,trntrtot
    norm     = norm  + 2.0D0*q2q(itr)*hq2q(itr) ! should be 2B
    buff1    = buff1 +      hp2q(itr)* c2q(itr)
    buff2    = buff2 +      hp2q(itr)*hp2q(itr)
    buff3    = buff3 +       c2q(itr)* c2q(itr)
  END DO
  norm  = 2.0D0*mass/norm
  errp  = 1.0D0-buff1**2/(buff2*buff3)
  
  DO itr=1,trntrtot
    q2q(itr) = q2q(itr)*norm
  END DO
  
  IF (pathmode == 0) THEN
    DO itr=1,trntrtot
      c2q(itr) = q2q(itr)
    END DO
  END IF
  
!        logging
  IF (pathmode == 0) THEN
    WRITE(22,'(a,i4,6(1pg12.3))') 'ATDHF: ',  &
        iter,mass,stiff,omgsq,qvar,pvar,errp
  ELSE
    WRITE(22,'(a,i4,6(1pg12.3))') 'CHF  : ',  &
        iter,mass,stiff,omgsq,qvar,pvar,errp
  END IF
  
!        check for convergency
  IF(ABS(errp) < 1E-6) GO TO 4002
  
END DO
4002 CONTINUE

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
!     ------------------------------------------------------------------
CALL upd_response_odd(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

!     evaluate the normalization commutator <[beta,P]> and normalize P
!     ------------------------------------------------------------------
mass = 0.0D0  ! B
DO itr=1,trntrtot
  hq2q(itr) = h2q(itr)*q2q(itr) + hresp2q(itr)
  mass      = mass + h2q(itr)*q2q(itr)*b2q(itr) +  &
      hresp2q(itr)     *b2q(itr)
END DO

norm=1.0D0/(2.0D0*mass)
DO itr=1,trntrtot
  p2q(itr)  = hq2q(itr)*norm
END DO

!     evaluate the normalization commutator <[Q,P]> and normalize Q
!     ------------------------------------------------------------------
norm = 0.0D0
DO itr=1,trntrtot
  norm      = norm + 2.0D0*q2q(itr)*p2q(itr)
END DO

norm=1.0/norm
DO itr=1,trntrtot
  q2q(itr)  = q2q(itr)*norm
END DO

!     evaluate the Liouville operator i[H,Q] + h_{Q,resp} and mass
!     ------------------------------------------------------------------
CALL upd_response_odd(q2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

CALL upd_hresp_odd(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0),  &
    redynpotr,redynpotp,redynpotz, respnpotr,respnpotp,respnpotz,  &
    respnmasr,respnmasp,respnmasz,il,imu)

DO itr=1,trntrtot
  hq2q(itr) =  h2q(itr)*q2q(itr) + hresp2q(itr)
END DO

!     evaluate the Liouville operator -i[H,P] + h_{P,resp} and stiffness
!     ------------------------------------------------------------------
CALL upd_response_even(p2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    redreffmas,redpeffmas,redzeffmas, resocpotr,resocpotp,resocpotz,il,imu)

CALL upd_hresp_even(hresp2q,  &
    cppsi(0,0,1,1,0),cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1),  &
    cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),cpvsqf(1,0), relocpot,reeffmas,  &
    resocpotr,resocpotp,resocpotz,il,imu)

DO itr=1,trntrtot
  hp2q(itr) =  h2q(itr)*p2q(itr) + hresp2q(itr)
END DO

!     check the norm and evaluate the variances, mass and stiffness
!     ------------------------------------------------------------------
mass =0.0D0
stiff=0.0D0
norm =0.0D0
norm0=0.0D0
qvar =0.0D0
pvar =0.0D0
DO itr=1,trntrtot
  norm0    = norm0 + 2.0D0*b2q(itr)*p2q(itr)
  norm     = norm  + 2.0D0*q2q(itr)*p2q(itr)
  mass     = mass  +       q2q(itr)*hq2q(itr)
  stiff    = stiff +       p2q(itr)*hp2q(itr)
  qvar     = qvar  +       q2q(itr)*q2q(itr)
  pvar     = pvar  +       p2q(itr)*p2q(itr)
END DO
omgsq = mass*stiff

!     logging
IF (pathmode == 0) THEN
  WRITE(21,'(a,7(1pe16.4))') 'ATDHF_VAL : ',  &
      mass,stiff,pvar,qvar,omgsq,norm,norm0
ELSE
  WRITE(21,'(a,7(1pe16.4))') 'CHF_VAL   : ',  &
      mass,stiff,pvar,qvar,omgsq,norm,norm0
END IF

!     save matrixelements of the Q- and P-Operator
!     ------------------------------------------------------------------
IF (pathmode == 0) THEN
  DO itr=1,trntrtot
    q2qatdhf(itr,imu) = q2q(itr)
    p2qatdhf(itr,imu) = p2q(itr)
  END DO
ELSE
  DO itr=1,trntrtot
    q2qchf(itr,imu) = q2q(itr)
    p2qchf(itr,imu) = p2q(itr)
  END DO
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE crank_atdhf


!=======================================================================

!     eval_com_QHQ(Comm,il,imu,iiso)

!       evaluates the particle number mass B_N=1/2*<[Q,[N,Q]]>

!=======================================================================

SUBROUTINE eval_com_qnq(n2q,q2q,mass,il,imu,iiso)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: n2q(1:mtdtot)
REAL(8), INTENT(IN)             :: q2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: mass
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN OUT)                  :: imu
INTEGER, INTENT(IN)                      :: iiso
!REAL(8) :: ! qh-matrixelements of N
!REAL(8) :: ! 2q-matrixelements of Q





INTEGER :: ist,jst,jstmin,jstmax
INTEGER :: trzpar
INTEGER :: itr
INTEGER :: iso

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     set mass to zero
!     ------------------------------------------------------------------
mass  = 0.0D0

!     loop over all sp-states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
  IF (iso /= 3-iiso) THEN
    DO itr=1,trntr(ist)
      
      itr2 = tridx(ist,itr,4)
      
      mass = mass + q2q(itr2)**2*n2q(itr2)
      
    END DO
  END IF
  
END DO

!      mass = mass   ! factors from 0.5 from 1/M to B and 2 from the
!                    ! commutator compensate

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a,i2,a)')  &
      'EVAL_COM_QNQ(): QNQ-commutator evaluated'//  &
      ' (lambda=',il,' mu=',imu,' iso=', iiso,')'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE eval_com_qnq

!=======================================================================

!     eval_comm_BP(Comm,il,imu,iiso)

!       evaluates the commutator -i<[B,P]> including contributions
!       from outside the qp-space.

!=======================================================================

SUBROUTINE eval_comm_bp(comm,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
USE wavefr
USE coulomb, ONLY: ylm
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: comm
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN OUT)                      :: imu




REAL(8) :: qlmoper(0:mzd,0:mrd)   ! multipole operator
REAL(8) :: pfac,normfac,qlmfac
REAL(8) :: peta,qeta,povlp,qovlp
REAL(8) :: z,r
INTEGER :: iz,ir
REAL(8) :: lastv,lastu,nextv,nextu,rcntv,rcntu
REAL(8) :: ui,uj,vi,vj
INTEGER :: itr,ist,jst,jstmin,jstmax

LOGICAL :: oldcrank
DATA             oldcrank/d_oldcrank/

!     ------------------------------------------------------------------

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     some factors
!     ------------------------------------------------------------------
normfac = 0.6D0*1.2D0**2*  &
    (npart(1)+npart(2))**(2.0D0/3.0D0)*(npart(1)+npart(2))
IF (imu == 0) THEN
  qlmfac   = 0.8D0*pi/normfac
  IF (oldcrank) THEN
    pfac     = 1.0D0/(cpbeta( 0)-cpbeta(-1))
  ELSE
    pfac     = 1.0D0/(cpbeta(+1)-cpbeta(-1))
  END IF
ELSE IF (imu == 1) THEN
  qlmfac   = 0.8D0*pi/(SQRT(2.0D0)*normfac)
  IF (ABS(cpbeta(0)) > 0.0D0) THEN
    pfac=1.0D0/(SQRT(3D0)*cpbeta(0))
  ELSE
    pfac=0.0D0
  END IF
ELSE
  WRITE(21,'(a,i2)') 'ERROR: only mu=0 and mu=1 possible -- not ',imu
  STOP 'in linresp'
END IF

!     set commutator value to zero
!     ------------------------------------------------------------------
comm = 0.0D0
communocc = 0.0D0

!     initialize multipole operator
!     ------------------------------------------------------------------
DO ir=0,nr
  r    = dr*(ir+sr)
  DO iz=0,nz
    z    = ((iz)+sz-icm)*dz - cpzcm(0)
    qlmoper(iz,ir) = qlmfac*ylm(z,r,2,imu)
  END DO
END DO

!     loop over all sp-states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
!        calculate action of P operator
!        ===============================================================
  
!        P_0 (generator of translation along the axial path, asym. def.)
!        --------------------------------------------------------------
  IF (imu == 0.AND.oldcrank) THEN
    
    DO ir=0,nr
      DO iz=0,nz
        lrbuff1(iz,ir,1,1)=pfac* (cppsi(iz,ir,1,ist, 0)-cppsi(iz,ir,1,ist,-1))
        lrbuff1(iz,ir,2,1)=pfac* (cppsi(iz,ir,2,ist, 0)-cppsi(iz,ir,2,ist,-1))
        lrbuff1(iz,ir,1,2)=0.0D0
        lrbuff1(iz,ir,2,2)=0.0D0
      END DO
    END DO
    
!        P_0 (generator of translation along the axial path, sym. def.)
!        --------------------------------------------------------------
  ELSE IF (imu == 0) THEN
    
    DO ir=0,nr
      DO iz=0,nz
        lrbuff1(iz,ir,1,1)=pfac* (cppsi(iz,ir,1,ist,+1)-cppsi(iz,ir,1,ist,-1))
        lrbuff1(iz,ir,2,1)=pfac* (cppsi(iz,ir,2,ist,+1)-cppsi(iz,ir,2,ist,-1))
        lrbuff1(iz,ir,1,2)=0.0D0
        lrbuff1(iz,ir,2,2)=0.0D0
      END DO
    END DO
    
!        P_1 (generator of rotation,  J_y = L_y + S_y)
!        --------------------------------------------------------------
  ELSE IF (imu == 1) THEN
    
  IF(tangmomsub) THEN
    CALL angmomop(cppsi(0,0,1,ist,ipos),lrbuff1(:,:,:,1),&
                  lrbuff1(:,:,:,2),nparz(1,ist),mpsi(ist))  
    lrbuff1=lrbuff1*0.5D0*pfac     ! to comply with normalization here
  ELSE
    DO ir=0,nr
      r = (ir+sr)*dr
      DO iz=0,nz
        z = (iz+sz-icm)*dz - cpzcm(0)
        
!                 ! upward laddering components of J_y
        lrbuff1(iz,ir,1,1) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,1,ist)  &
            +z*cpdppsi(iz,ir,1,ist) +r*cpdzpsi(iz,ir,1,ist)  &
            -  cppsi(iz,ir,2,ist,0))
        lrbuff1(iz,ir,2,1) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,2,ist)  &
            +z*cpdppsi(iz,ir,2,ist) +r*cpdzpsi(iz,ir,2,ist))
        
!                 ! downward laddering components of J_y
        lrbuff1(iz,ir,1,2) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,1,ist)  &
            -z*cpdppsi(iz,ir,1,ist) +r*cpdzpsi(iz,ir,1,ist))
        lrbuff1(iz,ir,2,2) =  0.5D0*pfac*( -z*cpdrpsi(iz,ir,2,ist)  &
            -z*cpdppsi(iz,ir,2,ist) +r*cpdzpsi(iz,ir,2,ist)  &
            +  cppsi(iz,ir,1,ist,0))
        
      END DO
    END DO

  END IF    
  END IF
  
!        calculate action of Q operator
!        ===============================================================
  IF (imu == 0) THEN
    DO ir=0,nr
      DO iz=0,nz
        lrbuff2(iz,ir,1,1)= qlmoper(iz,ir)*cppsi(iz,ir,1,ist,0)
        lrbuff2(iz,ir,2,1)= qlmoper(iz,ir)*cppsi(iz,ir,2,ist,0)
        lrbuff2(iz,ir,1,2)= 0.0D0
        lrbuff2(iz,ir,2,2)= 0.0D0
      END DO
    END DO
  ELSE
    DO ir=0,nr
      DO iz=0,nz
        lrbuff2(iz,ir,1,1)= qlmoper(iz,ir)*cppsi(iz,ir,1,ist,0)
        lrbuff2(iz,ir,2,1)= qlmoper(iz,ir)*cppsi(iz,ir,2,ist,0)
        lrbuff2(iz,ir,1,2)= qlmoper(iz,ir)*cppsi(iz,ir,1,ist,0)
        lrbuff2(iz,ir,2,2)= qlmoper(iz,ir)*cppsi(iz,ir,2,ist,0)
      END DO
    END DO
  END IF
  
!        adjust lower boundary in case of no z-symmetry
!        ---------------------------------------------------------------
  IF (isym == 0) THEN
    DO ir=0,nr
      lrbuff1(0,ir,1,1) = 0.0D0
      lrbuff1(0,ir,2,1) = 0.0D0
      lrbuff1(0,ir,1,2) = 0.0D0
      lrbuff1(0,ir,2,2) = 0.0D0
      lrbuff2(0,ir,1,1) = 0.0D0
      lrbuff2(0,ir,2,1) = 0.0D0
      lrbuff2(0,ir,1,2) = 0.0D0
      lrbuff2(0,ir,2,2) = 0.0D0
    END DO
  END IF
  
!        evaluate 2q-matrixlelements and 0q-part
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    itp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    qeta = uj*vi + vj*ui
    peta = uj*vi - vj*ui
    
    IF (itp == +1) THEN
      IF (msgn >= 0) THEN
        povlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff1(0:nz,0:nr,1:2,1))
        CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff1(0:nz,0:nr,1:2,1), povlp)
        qovlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff2(0:nz,0:nr,1:2,1))
        CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff2(0:nz,0:nr,1:2,1), qovlp)
      ELSE
        povlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff1(0:nz,0:nr,1:2,2))
        CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff1(0:nz,0:nr,1:2,2), povlp)
        qovlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff2(0:nr,0:nz,1:2,2))
        CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff2(0:nz,0:nr,1:2,2), qovlp)
      END IF
    ELSE
      povlp=wfovlptr(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff1(0:nz,0:nr,1:2,2))
      CALL wfprojtr(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff1(0:nz,0:nr,1:2,2), povlp)
      qovlp=wfovlptr(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff2(0:nz,0:nr,1:2,2))
      CALL wfprojtr(cppsi(0:nz,0:nr,1:2,jst,0),lrbuff2(0:nz,0:nr,1:2,2), qovlp)
    END IF
    
    IF (ist == jst.AND.itp == +1) THEN
      
      IF (oldcrank) THEN
        lastv   = sqrtb(cpvsqf(ist,-1))
        lastu   = sqrtb(one-cpvsqf(ist,-1))
        nextv   = sqrtb(cpvsqf(ist,0))
        nextu   = sqrtb(one-cpvsqf(ist,0))
        rcntv   = 0.5D0*(lastv+nextv)
        rcntu   = 0.5D0*(lastu+nextu)
        IF (rcntv >= 0.5D0) THEN
          povlp =  -(nextu-lastu)/rcntv*pfac
        ELSE
          povlp =  +(nextv-lastv)/rcntu*pfac
        END IF
      ELSE
        lastv   = sqrtb(cpvsqf(ist,-1))
        lastu   = sqrtb(one-cpvsqf(ist,-1))
        rcntv   = sqrtb(cpvsqf(ist,0))
        rcntu   = sqrtb(one-cpvsqf(ist,0))
        nextv   = sqrtb(cpvsqf(ist,+1))
        nextu   = sqrtb(one-cpvsqf(ist,+1))
        IF (rcntv >= 0.5D0) THEN
          povlp =  -(nextu-lastu)/rcntv*pfac
        ELSE
          povlp =  +(nextv-lastv)/rcntu*pfac
        END IF
      END IF
      
      comm = comm +  qeta*     povlp*qovlp
      
    ELSE
      
      comm = comm +  qeta*peta*povlp*qovlp
      
    END IF
    
  END DO
  
!        accumulate contributions from outside qp-space
!        ---------------------------------------------------------------
  communocc = communocc + 2.0D0*cpvsqf(ist,0)*(  &
      wfovlp(lrbuff2(0:nz,0:nr,1:2,1),lrbuff1(0:nz,0:nr,1:2,1))  &
      +wfovlp(lrbuff2(0:nz,0:nr,1:2,2),lrbuff1(0:nz,0:nr,1:2,2)) )
  
END DO

comm=comm+communocc
comm = 2.0D0*comm
communocc=2D0*communocc
IF(imu==1) WRITE(994,'(a,f8.5,2(1pg13.5))') '[Q,J]:',cpbeta(0),comm,communocc
CALL flush(994)
!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'EVAL_COMM_QP(): QP-commutator evaluated'//  &
      ' (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE eval_comm_bp


!=======================================================================

!     ini_matrix_N(N2q,N0q,NVar,il,imu,iiso)

!       initialize the particle number operator

!       and calculate 2-quasiparticle transition matrixelements

!       2qp-matrixelements    N^{2q}_ji = eta^(+)_ji <j|N|i>
!       expectation           <N>
!       variance              <N^2> - <N>^2

!=======================================================================

SUBROUTINE ini_matrix_n(n2q,n0q,nvar,il,imu,iiso)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: n2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: n0q
REAL(8), INTENT(OUT)            :: nvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN OUT)                  :: imu
INTEGER, INTENT(IN)                      :: iiso
!REAL(8) :: ! N^{qh}_ii + N^{qh}_jj






INTEGER :: ist,jst
INTEGER :: itr,itr2
INTEGER :: iso,msgn
REAL(8) :: occ,eta,xi
REAL(8) :: vi,ui,vj,uj


!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
ELSE IF (iiso < 0.OR.iiso > 2) THEN
  WRITE(21,'(a,i2)') 'ERROR: isospin index has to be 0,1 or 2 -- not ',iiso
  STOP 'in linresp'
END IF

!     set 0qp-part and variance to zero
!     ------------------------------------------------------------------
nvar = 0.0D0
n0q  = 0.0D0
DO itr2=1,trntrtot
  n2q(itr2)=0.0D0
END DO

!     loop over all sp-states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        evaluate variance, 0q-part and N2q=(N^{qh}_{ii} + N^{qh}_{jj})
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    itp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    eta  = uj*vi + vj*ui
    occ  = 2.0D0*vi**2
    
!           only diagonal componets contribute
    IF (ist == jst.AND.itp == +1.AND.iso /= 3-iiso) THEN
      nvar      = nvar + (eta)**2
      n0q       = n0q  +  occ
    END IF
    
    n2q(itr2) = (ui*ui - vi*vi + uj*uj - vj*vj)
    
  END DO
  
END DO

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a,i2,a)') 'INI_MATRIX_N(): N-operator initialized'//  &
      '  (lambda=',il,' mu=',imu,' iso=',iiso,')'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE ini_matrix_n


!=======================================================================

!     ini_matrix_B(B2q,B0q,BVar,il,imu,iiso,ipath)

!       initialize the deformation measurement-operator -- beta

!                   4 pi    Q_{l,+m} + (-1)^m Q_{l,-m}
!             Q = --------  --------------------------
!                  5 A r^2              2

!                   4 pi    Q_{l,+m} + Q^*_{l,m}
!               = --------  --------------------
!                  5 A r^2           2

!       r^2  = RMS radius       =  3/5 * R^2
!       A    = nucleon number
!       Q_lm = R^l Y_lm

!       and calculate

!       2qp-matrixelements    Q^{2q}_ji = eta^(+)_ji <j|Q|i>
!       expectation           <Q>
!       variance              <Q^2> - <Q>^2

!=======================================================================

SUBROUTINE ini_matrix_b(b2q,b0q,bvar,il,imu,iiso,ipath)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
USE wavefr
USE coulomb, ONLY: ylm
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: b2q(1:mtdtot)
REAL(8), INTENT(OUT)            :: b0q
REAL(8), INTENT(OUT)            :: bvar
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN OUT)                  :: imu
INTEGER, INTENT(IN)                      :: iiso
INTEGER, INTENT(IN)                      :: ipath
!REAL(8) :: ! 2qp-matrixelements







REAL(8) :: qlmoper(0:mzd,0:mrd)   ! multipole operator
INTEGER :: ist,jst,jstmin,jstmax
INTEGER :: itr,itr2
REAL(8) :: qlmfac,normfac
REAL(8) :: r,z
INTEGER :: ir,iz,iso,msgn,itp
REAL(8) :: ovlp,occ,eta
REAL(8) :: vi,ui,vj,uj

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF
IF (ipath > 1.OR.ipath < -1) THEN
  WRITE(21,'(a,i1,a)') 'ERROR: invalid path index (',ipath,  &
      ') in ini_matrix_beta'
  STOP 'in linresp'
END IF

!     calculate normalization of Q_mu for given isospin (0: total 1:p 2:n)
!     ------------------------------------------------------------------
IF (iiso == 0) THEN
  normfac = 0.6D0*1.2D0**2*  &
      (npart(1)+npart(2))**(2.0D0/3.0D0)*(npart(1)+npart(2))
ELSE IF (iiso == 1) THEN
  normfac = 0.6D0*1.2D0**2* (npart(1)+npart(2))**(2.0D0/3.0D0)*npart(1)
ELSE IF (iiso == 2) THEN
  normfac = 0.6D0*1.2D0**2* (npart(1)+npart(2))**(2.0D0/3.0D0)*npart(2)
END IF

!     initialize multipole operator
!     ------------------------------------------------------------------
IF (imu == 0) THEN
  qlmfac   = 0.8D0*pi/normfac
ELSE
  qlmfac   = 0.8D0*pi/(SQRT(2.0D0)*normfac)
END IF

DO ir=0,nr
  r    = dr*(ir+sr)
  DO iz=0,nz
    z    = ((iz)+sz-icm)*dz - cpzcm(ipath)
    qlmoper(iz,ir) = qlmfac*ylm(z,r,il,imu)
  END DO
END DO

!     adjust lower boundary in case of no z-symmetry
!     ------------------------------------------------------------------
IF (isym == 0) THEN
  DO ir=0,nr
    qlmoper(0,ir) = 0.0D0
  END DO
END IF

!     set variance, 0q-part and 2q-matrixelements to zero
!     ------------------------------------------------------------------
bvar = 0.0D0
b0q  = 0.0D0
DO itr2=1,trntrtot
  b2q(itr2)=0.0D0
END DO

!     loop over all sp-states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        calculate action of Q operator
!        ---------------------------------------------------------------
  IF (iso /= 3-iiso) THEN
    IF (imu == 0) THEN
      DO ir=0,nr
        DO iz=0,nz
          lrbuff1(iz,ir,1,1)= qlmoper(iz,ir)*cppsi(iz,ir,1,ist,ipath)
          lrbuff1(iz,ir,2,1)= qlmoper(iz,ir)*cppsi(iz,ir,2,ist,ipath)
          lrbuff1(iz,ir,1,2)= 0.0D0
          lrbuff1(iz,ir,2,2)= 0.0D0
        END DO
      END DO
    ELSE
      DO ir=0,nr
        DO iz=0,nz
          lrbuff1(iz,ir,1,1)= qlmoper(iz,ir)*cppsi(iz,ir,1,ist,ipath)
          lrbuff1(iz,ir,2,1)= qlmoper(iz,ir)*cppsi(iz,ir,2,ist,ipath)
          lrbuff1(iz,ir,1,2)= qlmoper(iz,ir)*cppsi(iz,ir,1,ist,ipath)
          lrbuff1(iz,ir,2,2)= qlmoper(iz,ir)*cppsi(iz,ir,2,ist,ipath)
        END DO
      END DO
    END IF
  ELSE
    DO ir=0,nr
      DO iz=0,nz
        lrbuff1(iz,ir,1,1)= 0.0D0
        lrbuff1(iz,ir,2,1)= 0.0D0
        lrbuff1(iz,ir,1,2)= 0.0D0
        lrbuff1(iz,ir,2,2)= 0.0D0
      END DO
    END DO
  END IF
  
!        evaluate 2q-matrixlelements and 0q-part
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    itp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,ipath))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,ipath))
    uj   = sqrtb(one-vj*vj)
    eta  = uj*vi + vj*ui
    occ  = 2.0D0*vi**2
    
    IF (itp == +1) THEN
      IF (msgn >= 0) THEN
        ovlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,ipath),lrbuff1(0:nz,0:nr,1:2,1))
        CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,ipath),lrbuff1(0:nz,0:nr,1:2,1), ovlp)
      ELSE
        ovlp=wfovlp(cppsi(0:nz,0:nr,1:2,jst,ipath),lrbuff1(0:nz,0:nr,1:2,2))
        CALL wfproj(cppsi(0:nz,0:nr,1:2,jst,ipath),lrbuff1(0:nz,0:nz,1:2,2), ovlp)
      END IF
    ELSE
      ovlp=wfovlptr(cppsi(0:nz,0:nr,1:2,jst,ipath),lrbuff1(0:nz,0:nr,1:2,2))
      CALL wfprojtr(cppsi(0:nz,0:nr,1:2,jst,ipath),lrbuff1(0:nz,0:nr,1:2,2), ovlp)
    END IF
    
    IF (ist == jst.AND.itp == +1) THEN
      b2q(itr2) =         eta * ovlp
      bvar      = bvar + (eta * ovlp)**2
      b0q       = b0q  +  occ * ovlp
    ELSE
      b2q(itr2) =         eta * ovlp
      bvar      = bvar + (eta * ovlp)**2
    END IF
    
  END DO
  
!        accumulate contributions from outside qp-space
!        ---------------------------------------------------------------
  bvar = bvar + 2.0D0*cpvsqf(ist,ipath)*(  &
      wfovlp(lrbuff1(0:nz,0:nr,1:2,1),lrbuff1(0:nz,0:nr,1:2,1))  &
      + wfovlp(lrbuff1(0:nz,0:nr,1:2,2),lrbuff1(0:nz,0:nr,1:2,2)) )
  
END DO

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a,i2,a)') 'INI_MATRIX_B(): B-operator initialized'//  &
      '  (lambda=',il,' mu=',imu,' iso=',iiso,')'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE ini_matrix_b


!=======================================================================

!     ini_matrix_H(H2q,H0q,il,imu,iiso)

!     initialize the quasiparticle energies

!=======================================================================

SUBROUTINE ini_matrix_h(h2q,h0q,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: h2q(1:mtdtot)
REAL(8), INTENT(OUT)                        :: h0q
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
!REAL(8) :: ! e_i + e_j
!REAL(8) :: h0q                    ! 0qp-part of H (energy)



INTEGER :: ist,jst
INTEGER :: itr,itr2
INTEGER :: iso
REAL(8) :: deltai,edifi,equasi,deltaj,edifj,equasj

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     calculate 0qp and 2qp-matrixelements of beta-operator
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     set 0q-part and qh-matrixelements to zero
!     ------------------------------------------------------------------
h0q  = 0.0D0
DO itr=1,trntrtot
  h2q(itr) = 0.0D0
END DO

!     loop over all sp-states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso=1
  ELSE
    iso=2
  END IF
  
  DO itr=1,trntr(ist)
    jst =tridx(ist,itr,1)
    itr2=tridx(ist,itr,4)
    
    deltai = cpdeltaf(ist,0)*cpwtpair(ist,0)
    edifi  = cpspenrg(ist,0)-cpfer(iso,0)
    equasi = sqrtb(edifi*edifi+deltai*deltai)
    
    deltaj = cpdeltaf(jst,0)*cpwtpair(jst,0)
    edifj  = cpspenrg(jst,0)-cpfer(iso,0)
    equasj = sqrtb(edifj*edifj+deltaj*deltaj)
    
    IF (ist /= jst) THEN
      h2q(itr2) = equasi + equasj
    ELSE
      h2q(itr2) = 2.0D0*equasi
      h0q       = h0q + 2.0D0*equasi
    END IF
    
  END DO
  
END DO

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'INI_MATRIX_H(): H0-operator initialized'//  &
      ' (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE ini_matrix_h


!=======================================================================

!     ini_transitions(il,imu)

!     count the number of transitions (matrixelements) for a given
!     cranking mode

!=======================================================================

SUBROUTINE ini_transtab(il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN OUT)                  :: imu



INTEGER :: ist,jst,jstmin,jstmax
INTEGER :: trzpar,maxtr
INTEGER :: itr,itr2
INTEGER :: iso

!     calculate z-parity of transition
!     ------------------------------------------------------------------
trzpar=(-1)**(imu+il)

!     calculate 0qp and 2qp-matrixelements of beta-operator
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
  
!     the case mu=0
!     ------------------------------------------------------------------
ELSE IF (imu == 0) THEN
  
!        loop over all sp-states
!        ---------------------------------------------------------------
  itr   = 1
  itr2  = 1
  maxtr = 0
  DO ist=1,nst(1)+nst(2)
    
    IF (ist <= nst(1)) THEN
      jstmin = 1
      jstmax = nst(1)
      iso    = 1
    ELSE
      jstmin = nst(1)+1
      jstmax = nst(1)+nst(2)
      iso    = 2
    END IF
    
    DO jst=jstmin,jstmax
      IF (mpsi(jst)-mpsi(ist) == 0)  THEN
        IF (npar(jst)*npar(ist)*trzpar*isym /= -1) THEN
          IF (itr <= mtd) THEN
            tridx(ist,itr,1) = jst  ! final state
            tridx(ist,itr,2) = +1   ! from original space
            tridx(ist,itr,3) = 0    ! angular momentum transfer
            tridx(ist,itr,4) = itr2 ! global transition index
          END IF
          itr              = itr  + 1
          itr2             = itr2 + 1
        END IF
      END IF
    END DO
    
    trntr(ist) = itr - 1
    maxtr      = MAX(itr-1,maxtr)
    itr        = 1
    
  END DO
  
!     the case mu>0
!     ------------------------------------------------------------------
ELSE
  
!        loop over all sp-states
!        ---------------------------------------------------------------
  itr   = 1
  itr2  = 1
  maxtr = 0
  DO ist=1,nst(1)+nst(2)
    
    IF (ist <= nst(1)) THEN
      jstmin = 1
      jstmax = nst(1)
      iso    = 1
    ELSE
      jstmin = nst(1)+1
      jstmax = nst(1)+nst(2)
      iso    = 2
    END IF
    
    DO jst=jstmin,jstmax
!              the upward laddering transitions
      IF (mpsi(jst)-mpsi(ist) == imu)  THEN
        IF (npar(jst)*npar(ist)*trzpar*isym /= -1) THEN
          IF (itr <= mtd) THEN
            tridx(ist,itr,1) = jst  ! final state
            tridx(ist,itr,2) = +1   ! from original space
            tridx(ist,itr,3) = imu  ! angular momentum transfer
            tridx(ist,itr,4) = itr2 ! global transition index
          END IF
          itr              = itr  + 1
          itr2             = itr2 + 1
        END IF
!              the downward laddering transitions inside the K>0 space
      ELSE IF (mpsi(jst)-mpsi(ist) == -imu)  THEN
        IF (npar(jst)*npar(ist)*trzpar*isym /= -1) THEN
          IF (itr <= mtd) THEN
            tridx(ist,itr,1) = jst  ! final state
            tridx(ist,itr,2) = +1   ! from original space
            tridx(ist,itr,3) = -imu ! angular momentum transfer
            tridx(ist,itr,4) = itr2 ! global transition index
          END IF
          itr              = itr  + 1
          itr2             = itr2 + 1
        END IF
!              the downward laddering components to the K<0 space
      ELSE IF (mpsi(ist)+mpsi(jst)+1 == imu)  THEN
        IF (npar(jst)*npar(ist)*trzpar*isym /= +1) THEN
          IF (itr <= mtd) THEN
            tridx(ist,itr,1) = jst  ! final state
            tridx(ist,itr,2) = -1   ! from time reversed space
            tridx(ist,itr,3) = -imu ! angular momentum transfer
            tridx(ist,itr,4) = itr2 ! global transition index
          END IF
          itr              = itr  + 1
          itr2             = itr2 + 1
        END IF
      END IF
    END DO
    
    trntr(ist) = itr - 1
    maxtr      = MAX(itr-1,maxtr)
    itr        = 1
    
  END DO
END IF

trntrtot = itr2 - 1

!     ------------------------------------------------------------------

IF (maxtr > mtd) THEN
  WRITE(21,'(a,6i3,a,6i4)')  &
      'ERROR: number of transition too big (enhance MTD) ',  &
      ' maxtr=',maxtr,' MTD=',mtd
  STOP 'in linresp'
END IF

IF (trntrtot > mtdtot) THEN
  WRITE(21,'(a,6i3,a,6i6)')  &
      'ERROR: number of transition too big (enhance MTDTOT) ',  &
      ' TRntr=',trntrtot,' MTD=',mtdtot
  STOP 'in linresp'
END IF

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a,i2,a,i6,a)')  &
      'INI_TRANSTAB(): transitions initialized'//  &
      ' (lambda=',il,' mu=',imu,' pi= ',trzpar,' ntr=',trntrtot,')'
END IF

!     ------------------------------------------------------------------

RETURN

END SUBROUTINE ini_transtab


!=======================================================================

!     upd_response_odd(Q2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
!                      RESP_Ar,RESP_Ap,RESP_Az,
!                      RESP_Sr,RESP_Sp,RESP_Sz,
!                      RESP_Cr,RESP_Cp,RESP_Cz,il,imu)

!     updates the time-odd mean-field potentials (and densities)

!=======================================================================

SUBROUTINE upd_response_odd(q2q, phi,drphi,dzphi,dpphi,ddphi,occphi,  &
    resp_ar,resp_ap,resp_az, resp_sr,resp_sp,resp_sz,  &
    resp_cr,resp_cp,resp_cz,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: q2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: phi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: occphi(mpsid)
REAL(8), INTENT(OUT)            :: resp_ar(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_ap(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_az(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_sr(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_sp(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_sz(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_cr(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_cp(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_cz(0:mzd,0:mrd,2)
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu




















REAL(8) :: resp_delta(0:mzd,0:mrd,2)

REAL(8) :: resp_spnr(0:mzd,0:mrd,2)
REAL(8) :: resp_spnp(0:mzd,0:mrd,2)
REAL(8) :: resp_spnz(0:mzd,0:mrd,2)
REAL(8) :: resp_skdr(0:mzd,0:mrd,2)
REAL(8) :: resp_skdp(0:mzd,0:mrd,2)
REAL(8) :: resp_skdz(0:mzd,0:mrd,2)
REAL(8) :: resp_curr(0:mzd,0:mrd,2)
REAL(8) :: resp_curp(0:mzd,0:mrd,2)
REAL(8) :: resp_curz(0:mzd,0:mrd,2)
REAL(8) :: resp_paio(0:mzd,0:mrd,2)
REAL(8) :: resp_rotjr(0:mzd,0:mrd,2)
REAL(8) :: resp_rotjp(0:mzd,0:mrd,2)
REAL(8) :: resp_rotjz(0:mzd,0:mrd,2)
REAL(8) :: resp_rotsr(0:mzd,0:mrd,2)
REAL(8) :: resp_rotsp(0:mzd,0:mrd,2)
REAL(8) :: resp_rotsz(0:mzd,0:mrd,2)

INTEGER :: ist,jst,jtp,msgn
INTEGER :: itr,iz,ir
INTEGER :: iso
REAL(8) :: eta,xi,wgt1,wgt2,wgt3

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     set response densities to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      resp_spnr(iz,ir,iso)  = 0.0D0
      resp_spnp(iz,ir,iso)  = 0.0D0
      resp_spnz(iz,ir,iso)  = 0.0D0
      resp_skdr(iz,ir,iso)  = 0.0D0
      resp_skdp(iz,ir,iso)  = 0.0D0
      resp_skdz(iz,ir,iso)  = 0.0D0
      resp_curr(iz,ir,iso)  = 0.0D0
      resp_curp(iz,ir,iso)  = 0.0D0
      resp_curz(iz,ir,iso)  = 0.0D0
      resp_paio(iz,ir,iso)  = 0.0D0
      resp_rotjr(iz,ir,iso) = 0.0D0
      resp_rotjp(iz,ir,iso) = 0.0D0
      resp_rotjz(iz,ir,iso) = 0.0D0
      resp_rotsr(iz,ir,iso) = 0.0D0
      resp_rotsp(iz,ir,iso) = 0.0D0
      resp_rotsz(iz,ir,iso) = 0.0D0
    END DO
  END DO
END DO

!     loop all single particle states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        loop all transitions
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    jtp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    
!           if the final state is from the original space
    IF (jtp == +1) THEN
      
      IF (msgn == 0) THEN
        wgt1=   0.0D0                           ! for   e^{i*mu*phi}
        wgt2=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for -ie^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      ELSE IF (msgn > 0) THEN
        wgt1=  +2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for   e^{i*mu*phi}
        wgt2=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for -ie^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      ELSE
        wgt1=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for   e^{i*mu*phi}
        wgt2=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for -ie^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      END IF
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
          
          resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
          
          resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          
          resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) + 0.5D0*wgt1*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          
          resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          
          resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(  &
              - dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(  &
              + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist) )
          
          resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(  &
              + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
        END DO
      END DO
      
!           if the final state is from the time-reversed space
    ELSE
      
      wgt1=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for   e^{i*mu*phi}
      wgt2=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for -ie^{i*mu*phi}
      wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
          
          resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)  &
              - phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
          
          resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              - phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          
          resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) + 0.5D0*wgt1*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          
          resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          
          resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(  &
              - dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist) )
          
          resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(  &
              + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
        END DO
      END DO
      
    END IF
    
  END DO    ! end of loop over all transitions from state ist
END DO       ! end of loop over all states ist

!     recouple the Skyrme-Parameters
!     ------------------------------------------------------------------
b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
b4  = 0.500D0*       t4
b5  = 0.125D0*(   x1*t1 +       x2*t2)  ! 2*b5
b5p = 0.125D0*(      t1 -          t2)  ! 2*b5'

wgt_b1  = -2.0D0*(b1-b1p)
wgt_b1c = -2.0D0* b1
IF (b4p == 0.0D0) THEN    ! no isotopic spin-orbit force
  wgt_b4  = - 2.0D0*b4
  wgt_b4c = -       b4
ELSE
  wgt_b4  = -(b4+b4p)
  wgt_b4c = - b4
END IF
IF(ifls == 1) THEN
  wgt_b5  =  (b5-b5p)    ! b5 already includes the factor 2
  wgt_b5c =   b5
ELSE
  wgt_b5  = 0.0D0
  wgt_b5c = 0.0D0
END IF

!     set the mean-fields to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      resp_ar(iz,ir,iso) = 0.0D0
      resp_ap(iz,ir,iso) = 0.0D0
      resp_az(iz,ir,iso) = 0.0D0
      resp_sr(iz,ir,iso) = 0.0D0
      resp_sp(iz,ir,iso) = 0.0D0
      resp_sz(iz,ir,iso) = 0.0D0
      resp_cr(iz,ir,iso) = 0.0D0
      resp_cp(iz,ir,iso) = 0.0D0
      resp_cz(iz,ir,iso) = 0.0D0
    END DO
  END DO
END DO

!     construct the mean-fields from the densities
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      
      resp_ar(iz,ir,iso) = resp_ar(iz,ir,iso)  &
          + wgt_b1 *resp_curr(iz,ir,iso) + wgt_b1c*resp_curr(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotsr(iz,ir,iso) + wgt_b4c*resp_rotsr(iz,ir,3-iso)
      resp_ap(iz,ir,iso) = resp_ap(iz,ir,iso)  &
          + wgt_b1 *resp_curp(iz,ir,iso) + wgt_b1c*resp_curp(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotsp(iz,ir,iso) + wgt_b4c*resp_rotsp(iz,ir,3-iso)
      resp_az(iz,ir,iso) = resp_az(iz,ir,iso)  &
          + wgt_b1 *resp_curz(iz,ir,iso) + wgt_b1c*resp_curz(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotsz(iz,ir,iso) + wgt_b4c*resp_rotsz(iz,ir,3-iso)
      
      resp_sr(iz,ir,iso) = resp_sr(iz,ir,iso)  &
          + wgt_b5 *resp_skdr(iz,ir,iso) + wgt_b5c*resp_skdr(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotjr(iz,ir,iso) + wgt_b4c*resp_rotjr(iz,ir,3-iso)
      resp_sp(iz,ir,iso) = resp_sp(iz,ir,iso)  &
          + wgt_b5 *resp_skdp(iz,ir,iso) + wgt_b5c*resp_skdp(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotjp(iz,ir,iso) + wgt_b4c*resp_rotjp(iz,ir,3-iso)
      resp_sz(iz,ir,iso) = resp_sz(iz,ir,iso)  &
          + wgt_b5 *resp_skdz(iz,ir,iso) + wgt_b5c*resp_skdz(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotjz(iz,ir,iso) + wgt_b4c*resp_rotjz(iz,ir,3-iso)
      
      resp_cr(iz,ir,iso) = resp_cr(iz,ir,iso)  &
          + wgt_b5 *resp_spnr(iz,ir,iso) + wgt_b5c*resp_spnr(iz,ir,3-iso)
      resp_cp(iz,ir,iso) = resp_cp(iz,ir,iso)  &
          + wgt_b5 *resp_spnp(iz,ir,iso) + wgt_b5c*resp_spnp(iz,ir,3-iso)
      resp_cz(iz,ir,iso) = resp_cz(iz,ir,iso)  &
          + wgt_b5 *resp_spnz(iz,ir,iso) + wgt_b5c*resp_spnz(iz,ir,3-iso)
      
    END DO
  END DO
END DO

!     ALL FIELDS NOW CONTAIN THE (r,z)-PART OF THE RESPONSE FIELDS
!     THE ANGULAR DEPENDENCY READS

!     Cr,Cz,Sr,Sz,Ap  ~  sin(mu*phi)*X(r,z)
!     Cp,Sp,Ar,Az     ~  cos(mu*phi)*X(r,z)

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'UPD_RESP_ODD(): time-odd response updated'//  &
      '  (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN

END SUBROUTINE upd_response_odd


!=======================================================================

!     upd_response_even(P2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
!                       RESP_U,RESP_B,RESP_drB,RESP_dzB,RESP_Wz,RESP_Wr,


!     count the number of transitions (matrixelements) for a given
!     cranking mode

!=======================================================================

SUBROUTINE upd_response_ptest(p2q, phi,drphi,dzphi,dpphi,ddphi,occphi,  &
    resp_u, resp_b,resp_drb,resp_dpb,resp_dzb,  &
    resp_wr,resp_wp,resp_wz,il,imu, a1,a2,a3,a4,a5,a6)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: p2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: phi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: occphi(mpsid)
REAL(8), INTENT(OUT)            :: resp_u(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_b(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_drb(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_dpb(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_dzb(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_wr(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_wp(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_wz(0:mzd,0:mrd,2)
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
REAL(8), INTENT(IN)             :: a1
REAL(8), INTENT(IN)             :: a2
REAL(8), INTENT(IN)             :: a3
REAL(8), INTENT(IN)             :: a4
REAL(8), INTENT(IN)             :: a5
REAL(8), INTENT(IN)             :: a6



















REAL(8) :: resp_delta(0:mzd,0:mrd,2)

REAL(8) :: resp_rho(0:mzd,0:mrd,2)
REAL(8) :: resp_tau(0:mzd,0:mrd,2)
REAL(8) :: resp_socr(0:mzd,0:mrd,2)
REAL(8) :: resp_socp(0:mzd,0:mrd,2)
REAL(8) :: resp_socz(0:mzd,0:mrd,2)
REAL(8) :: resp_grhor(0:mzd,0:mrd,2)
REAL(8) :: resp_grhop(0:mzd,0:mrd,2)
REAL(8) :: resp_grhoz(0:mzd,0:mrd,2)
REAL(8) :: resp_drho(0:mzd,0:mrd,2)
REAL(8) :: resp_divsoc(0:mzd,0:mrd,2)
REAL(8) :: resp_pair(0:mzd,0:mrd,2)

INTEGER :: ist,jst,jtp
INTEGER :: itr,iz,ir
INTEGER :: iso
REAL(8) :: eta,xi,wgt1,wgt2,wgt3

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     set response densities to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO iz=0,nz
    DO ir=0,nr
      resp_rho(iz,ir,iso)    = 0.0D0
      resp_tau(iz,ir,iso)    = 0.0D0
      resp_socr(iz,ir,iso)   = 0.0D0
      resp_socp(iz,ir,iso)   = 0.0D0
      resp_socz(iz,ir,iso)   = 0.0D0
      resp_grhor(iz,ir,iso)  = 0.0D0
      resp_grhop(iz,ir,iso)  = 0.0D0
      resp_grhoz(iz,ir,iso)  = 0.0D0
      resp_drho(iz,ir,iso)   = 0.0D0
      resp_divsoc(iz,ir,iso) = 0.0D0
      resp_pair(iz,ir,iso)   = 0.0D0
    END DO
  END DO
END DO

!     loop all single particle states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        loop all transitions
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    jtp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    
!           if the final state is from the original space
    IF (jtp == +1) THEN
      
      IF (msgn == 0) THEN
        wgt1=   0.0D0                           ! for -ie^{i*mu*phi}
        wgt2=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      ELSE IF (msgn > 0) THEN
        wgt1=   2.0D0*a1*(uj*vi + vj*ui)*p2q(itr2) ! for -ie^{i*mu*phi}
        wgt2=   2.0D0*a2*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      ELSE
        wgt1=   2.0D0*a3*(uj*vi + vj*ui)*p2q(itr2) ! for -ie^{i*mu*phi}
        wgt2=   2.0D0*a4*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      END IF
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 0.5D0*wgt1*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          resp_grhor(iz,ir,iso) = resp_grhor(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhop(iz,ir,iso) = resp_grhop(iz,ir,iso) + wgt1*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhoz(iz,ir,iso) = resp_grhoz(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
!                 2*resp_tau has to be added later
          resp_drho(iz,ir,iso) = resp_drho(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*ddphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*ddphi(iz,ir,2,ist)  &
              + ddphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + ddphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(  &
              + drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
          
        END DO
      END DO
      
!           if the final state is from the time-reversed space
    ELSE
      
      wgt1=   2.0D0*a5*(uj*vi + vj*ui)*p2q(itr2) ! for -ie^{i*mu*phi}
      wgt2=   2.0D0*a6*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
      wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 0.5D0*wgt1*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          resp_grhor(iz,ir,iso) = resp_grhor(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhop(iz,ir,iso) = resp_grhop(iz,ir,iso) + wgt1*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhoz(iz,ir,iso) = resp_grhoz(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
!                 2*resp_tau has to be added later
          resp_drho(iz,ir,iso) = resp_drho(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*ddphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*ddphi(iz,ir,2,ist)  &
              - ddphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + ddphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(  &
              - drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
          
        END DO
      END DO
      
    END IF
    
  END DO    ! end of loop over all transitions from state ist
END DO       ! end of loop over all states ist

!     add 2*resp_tau to get the full Drho response
!     ------------------------------------------------------------------
DO iz=0,nz
  DO ir=0,nr
    resp_drho(iz,ir,1) = resp_drho(iz,ir,1) + 2.0D0*resp_tau(iz,ir,1)
    resp_drho(iz,ir,2) = resp_drho(iz,ir,2) + 2.0D0*resp_tau(iz,ir,2)
  END DO
END DO

!     recouple the Skyrme-Parameters
!     ------------------------------------------------------------------
b0  =                t0 + 0.5D0*t0*x0
b0p =          0.5D0*t0 +       t0*x0
b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
b2  = 0.125D0*(3.0D0*t1 + 1.5D0*t1*x1 -       t2 - 0.5D0*t2*x2)
b2p = 0.125D0*(1.5D0*t1 + 3.0D0*t1*x1 +   0.5*t2 +       t2*x2)
b3  = 0.250D0*(      t3 + 0.5D0*t3*x3)
b3p = 0.250D0*(0.5D0*t3 +       t3*x3)
b4  = 0.500D0*       t4
b5  = 0.125D0*(   x1*t1 +       x2*t2)  ! 2*b5
b5p = 0.125D0*(      t1 -          t2)  ! 2*b5'

wgt_b0  =  (b0-b0p)
wgt_b0c =   b0
wgt_b1  =  (b1-b1p)
wgt_b1c =   b1
wgt_b2  = -(b2-b2p)
wgt_b2c = - b2
IF (b4p == 0.0D0) THEN    ! no isotopic spin-orbit force
  wgt_b4  = - 2.0D0*b4
  wgt_b4c = -       b4
ELSE
  wgt_b4  = -(b4+b4p)
  wgt_b4c = - b4
END IF
IF(ifls == 1) THEN
  wgt_b5  =  -(b5-b5p)    ! b5 already includes the factor 2
  wgt_b5c =  - b5
ELSE
  wgt_b5  = 0.0D0
  wgt_b5c = 0.0D0
END IF

!     set the mean-fields to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      resp_u(iz,ir,iso) = 0.0D0
      resp_b(iz,ir,iso) = 0.0D0
      resp_wr(iz,ir,iso) = 0.0D0
      resp_wp(iz,ir,iso) = 0.0D0
      resp_wz(iz,ir,iso) = 0.0D0
!            RESP_drB(iz,ir,iso) = 0.0D0
!            RESP_dpB(iz,ir,iso) = 0.0D0
!            RESP_dzB(iz,ir,iso) = 0.0D0
    END DO
  END DO
END DO

!     construct the mean-fields from the densities
!     ------------------------------------------------------------------

DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      
!           the mean-field response without density-dependency (b3)
!           and coulomb
      resp_u(iz,ir,iso)  = resp_u(iz,ir,iso) + wgt_b0 *resp_rho(iz,ir,iso)  &
          + wgt_b0c*resp_rho(iz,ir,3-iso) + wgt_b1 *resp_tau(iz,ir,iso)  &
          + wgt_b1c*resp_tau(iz,ir,3-iso) + wgt_b2 *resp_drho(iz,ir,iso)  &
          + wgt_b2c*resp_drho(iz,ir,3-iso) + wgt_b4 *resp_divsoc(iz,ir,iso)  &
          + wgt_b4c*resp_divsoc(iz,ir,3-iso)
      
      resp_b(iz,ir,iso)  = resp_b(iz,ir,iso) + wgt_b1 *resp_rho(iz,ir,iso)  &
          + wgt_b1c*resp_rho(iz,ir,3-iso)
      
      resp_wr(iz,ir,iso) = resp_wr(iz,ir,iso)  &
          - wgt_b4 *resp_grhor(iz,ir,iso) - wgt_b4c*resp_grhor(iz,ir,3-iso)  &
          + wgt_b5 *resp_socr(iz,ir,iso) + wgt_b5c*resp_socr(iz,ir,3-iso)
      
      resp_wp(iz,ir,iso) = resp_wp(iz,ir,iso)  &
          - wgt_b4 *resp_grhop(iz,ir,iso) - wgt_b4c*resp_grhop(iz,ir,3-iso)  &
          + wgt_b5 *resp_socp(iz,ir,iso) + wgt_b5c*resp_socp(iz,ir,3-iso)
      
      resp_wz(iz,ir,iso) = resp_wz(iz,ir,iso)  &
          - wgt_b4 *resp_grhoz(iz,ir,iso) - wgt_b4c*resp_grhoz(iz,ir,3-iso)  &
          + wgt_b5 *resp_socz(iz,ir,iso) + wgt_b5c*resp_socz(iz,ir,3-iso)
      
!            RESP_drB(iz,ir,iso) = RESP_drB(iz,ir,iso)
!     &                         + wgt_b1 *resp_Grhor(iz,ir,iso)
!     &                         + wgt_b1c*resp_Grhor(iz,ir,3-iso)
!c
!            RESP_dpB(iz,ir,iso) = RESP_dpB(iz,ir,iso)
!     &                         + wgt_b1 *resp_Grhop(iz,ir,iso)
!     &                         + wgt_b1c*resp_Grhop(iz,ir,3-iso)
!c
!            RESP_dzB(iz,ir,iso) = RESP_dzB(iz,ir,iso)
!     &                         + wgt_b1 *resp_Grhoz(iz,ir,iso)
!     &                         + wgt_b1c*resp_Grhoz(iz,ir,3-iso)
!c
    END DO
  END DO
END DO

!     now the contribution from the density-dependency (b3)
!     ------------------------------------------------------------------
DO ir=0,nr
  DO iz=0,nz
    resp_u(iz,ir,1) = resp_u(iz,ir,1)  &
        + locpotwgt1(iz,ir,1)*resp_rho(iz,ir,1)  &
        + locpotwgt1(iz,ir,2)*resp_rho(iz,ir,2)
    resp_u(iz,ir,2) = resp_u(iz,ir,2)  &
        + locpotwgt2(iz,ir,1)*resp_rho(iz,ir,1)  &
        + locpotwgt2(iz,ir,2)*resp_rho(iz,ir,2)
  END DO
END DO

!     and from the coulomb
!     ------------------------------------------------------------------
IF (chg /= 0.0D0) THEN
  
END IF

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'UPD_RESP_EVEN(): time-even response updated'//  &
      '  (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN

END SUBROUTINE upd_response_ptest


!=======================================================================

!     upd_response_odd(Q2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
!                      RESP_Ar,RESP_Ap,RESP_Az,
!                      RESP_Sr,RESP_Sp,RESP_Sz,
!                      RESP_Cr,RESP_Cp,RESP_Cz,il,imu)

!     updates the time-odd mean-field potentials (and densities)

!=======================================================================

SUBROUTINE upd_response_qtest(q2q, phi,drphi,dzphi,dpphi,ddphi,occphi,  &
    resp_ar,resp_ap,resp_az, resp_sr,resp_sp,resp_sz,  &
    resp_cr,resp_cp,resp_cz,il,imu, a1,a2,a3,a4,a5,a6)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: q2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: phi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: occphi(mpsid)
REAL(8), INTENT(OUT)            :: resp_ar(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_ap(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_az(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_sr(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_sp(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_sz(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_cr(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_cp(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_cz(0:mzd,0:mrd,2)
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu
REAL(8), INTENT(IN)             :: a1
REAL(8), INTENT(IN)             :: a2
REAL(8), INTENT(IN)             :: a3
REAL(8), INTENT(IN)             :: a4
REAL(8), INTENT(IN)             :: a5
REAL(8), INTENT(IN)             :: a6




















REAL(8) :: resp_delta(0:mzd,0:mrd,2)

REAL(8) :: resp_spnr(0:mzd,0:mrd,2)
REAL(8) :: resp_spnp(0:mzd,0:mrd,2)
REAL(8) :: resp_spnz(0:mzd,0:mrd,2)
REAL(8) :: resp_skdr(0:mzd,0:mrd,2)
REAL(8) :: resp_skdp(0:mzd,0:mrd,2)
REAL(8) :: resp_skdz(0:mzd,0:mrd,2)
REAL(8) :: resp_curr(0:mzd,0:mrd,2)
REAL(8) :: resp_curp(0:mzd,0:mrd,2)
REAL(8) :: resp_curz(0:mzd,0:mrd,2)
REAL(8) :: resp_paio(0:mzd,0:mrd,2)
REAL(8) :: resp_rotjr(0:mzd,0:mrd,2)
REAL(8) :: resp_rotjp(0:mzd,0:mrd,2)
REAL(8) :: resp_rotjz(0:mzd,0:mrd,2)
REAL(8) :: resp_rotsr(0:mzd,0:mrd,2)
REAL(8) :: resp_rotsp(0:mzd,0:mrd,2)
REAL(8) :: resp_rotsz(0:mzd,0:mrd,2)

INTEGER :: ist,jst,jtp,msgn
INTEGER :: itr,iz,ir
INTEGER :: iso
REAL(8) :: eta,xi,wgt1,wgt2,wgt3

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     set response densities to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      resp_spnr(iz,ir,iso)  = 0.0D0
      resp_spnp(iz,ir,iso)  = 0.0D0
      resp_spnz(iz,ir,iso)  = 0.0D0
      resp_skdr(iz,ir,iso)  = 0.0D0
      resp_skdp(iz,ir,iso)  = 0.0D0
      resp_skdz(iz,ir,iso)  = 0.0D0
      resp_curr(iz,ir,iso)  = 0.0D0
      resp_curp(iz,ir,iso)  = 0.0D0
      resp_curz(iz,ir,iso)  = 0.0D0
      resp_paio(iz,ir,iso)  = 0.0D0
      resp_rotjr(iz,ir,iso) = 0.0D0
      resp_rotjp(iz,ir,iso) = 0.0D0
      resp_rotjz(iz,ir,iso) = 0.0D0
      resp_rotsr(iz,ir,iso) = 0.0D0
      resp_rotsp(iz,ir,iso) = 0.0D0
      resp_rotsz(iz,ir,iso) = 0.0D0
    END DO
  END DO
END DO

!     loop all single particle states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        loop all transitions
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    jtp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    
!           if the final state is from the original space
    IF (jtp == +1) THEN
      
      IF (msgn == 0) THEN
        wgt1=   0.0D0                           ! for   e^{i*mu*phi}
        wgt2=  -2.0D0*(uj*vi - vj*ui)*q2q(itr2) ! for -ie^{i*mu*phi}
        wgt3=  +2.0D0*(uj*ui + vj*vi)*q2q(itr2)
      ELSE IF (msgn > 0) THEN
        wgt1=   2.0D0*a1*(uj*vi - vj*ui)*q2q(itr2) ! for   e^{i*mu*phi}
        wgt2=   2.0D0*a2*(uj*vi - vj*ui)*q2q(itr2) ! for -ie^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      ELSE
        wgt1=   2.0D0*a3*(uj*vi - vj*ui)*q2q(itr2) !    e^{i*mu*phi}
        wgt2=   2.0D0*a4*(uj*vi - vj*ui)*q2q(itr2) !  -ie^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      END IF
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
          
          resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
          
          resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          
          resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) + 0.5D0*wgt1*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          
          resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          
          resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(  &
              - dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(  &
              + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist) )
          
          resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(  &
              + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
        END DO
      END DO
      
!           if the final state is from the time-reversed space
    ELSE
      
      wgt1=   2.0D0*a5*(uj*vi - vj*ui)*q2q(itr2) !    e^{i*mu*phi}
      wgt2=   2.0D0*a6*(uj*vi - vj*ui)*q2q(itr2) !  -ie^{i*mu*phi}
      wgt3=   0.0D0*(uj*ui + vj*vi)*q2q(itr2)
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
          
          resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)  &
              - phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
          
          resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              - phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          
          resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
          
          resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) + 0.5D0*wgt1*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          
          resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          
          resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(  &
              - dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist) )
          
          resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(  &
              + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          
          resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
        END DO
      END DO
      
    END IF
    
  END DO    ! end of loop over all transitions from state ist
END DO       ! end of loop over all states ist

!     recouple the Skyrme-Parameters
!     ------------------------------------------------------------------
b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
b4  = 0.500D0*       t4
b5  = 0.125D0*(   x1*t1 +       x2*t2)  ! 2*b5
b5p = 0.125D0*(      t1 -          t2)  ! 2*b5'

wgt_b1  = -2.0D0*(b1-b1p)
wgt_b1c = -2.0D0* b1
IF (b4p == 0.0D0) THEN    ! no isotopic spin-orbit force
  wgt_b4  = - 2.0D0*b4
  wgt_b4c = -       b4
ELSE
  wgt_b4  = -(b4+b4p)
  wgt_b4c = - b4
END IF
IF(ifls == 1) THEN
  wgt_b5  =  (b5-b5p)    ! b5 already includes the factor 2
  wgt_b5c =   b5
ELSE
  wgt_b5  = 0.0D0
  wgt_b5c = 0.0D0
END IF

!     set the mean-fields to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      resp_ar(iz,ir,iso) = 0.0D0
      resp_ap(iz,ir,iso) = 0.0D0
      resp_az(iz,ir,iso) = 0.0D0
      resp_sr(iz,ir,iso) = 0.0D0
      resp_sp(iz,ir,iso) = 0.0D0
      resp_sz(iz,ir,iso) = 0.0D0
      resp_cr(iz,ir,iso) = 0.0D0
      resp_cp(iz,ir,iso) = 0.0D0
      resp_cz(iz,ir,iso) = 0.0D0
    END DO
  END DO
END DO

!     construct the mean-fields from the densities
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      
      resp_ar(iz,ir,iso) = resp_ar(iz,ir,iso)  &
          + wgt_b1 *resp_curr(iz,ir,iso) + wgt_b1c*resp_curr(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotsr(iz,ir,iso) + wgt_b4c*resp_rotsr(iz,ir,3-iso)
      resp_ap(iz,ir,iso) = resp_ap(iz,ir,iso)  &
          + wgt_b1 *resp_curp(iz,ir,iso) + wgt_b1c*resp_curp(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotsp(iz,ir,iso) + wgt_b4c*resp_rotsp(iz,ir,3-iso)
      resp_az(iz,ir,iso) = resp_az(iz,ir,iso)  &
          + wgt_b1 *resp_curz(iz,ir,iso) + wgt_b1c*resp_curz(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotsz(iz,ir,iso) + wgt_b4c*resp_rotsz(iz,ir,3-iso)
      
      resp_sr(iz,ir,iso) = resp_sr(iz,ir,iso)  &
          + wgt_b5 *resp_skdr(iz,ir,iso) + wgt_b5c*resp_skdr(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotjr(iz,ir,iso) + wgt_b4c*resp_rotjr(iz,ir,3-iso)
      resp_sp(iz,ir,iso) = resp_sp(iz,ir,iso)  &
          + wgt_b5 *resp_skdp(iz,ir,iso) + wgt_b5c*resp_skdp(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotjp(iz,ir,iso) + wgt_b4c*resp_rotjp(iz,ir,3-iso)
      resp_sz(iz,ir,iso) = resp_sz(iz,ir,iso)  &
          + wgt_b5 *resp_skdz(iz,ir,iso) + wgt_b5c*resp_skdz(iz,ir,3-iso)  &
          + wgt_b4 *resp_rotjz(iz,ir,iso) + wgt_b4c*resp_rotjz(iz,ir,3-iso)
      
      resp_cr(iz,ir,iso) = resp_cr(iz,ir,iso)  &
          + wgt_b5 *resp_spnr(iz,ir,iso) + wgt_b5c*resp_spnr(iz,ir,3-iso)
      resp_cp(iz,ir,iso) = resp_cp(iz,ir,iso)  &
          + wgt_b5 *resp_spnp(iz,ir,iso) + wgt_b5c*resp_spnp(iz,ir,3-iso)
      resp_cz(iz,ir,iso) = resp_cz(iz,ir,iso)  &
          + wgt_b5 *resp_spnz(iz,ir,iso) + wgt_b5c*resp_spnz(iz,ir,3-iso)
      
    END DO
  END DO
END DO

!     ALL FIELDS NOW CONTAIN THE (r,z)-PART OF THE RESPONSE FIELDS
!     THE ANGULAR DEPENDENCY READS

!     Cr,Cz,Sr,Sz,Ap  ~  sin(mu*phi)*X(r,z)
!     Cp,Sp,Ar,Az     ~  cos(mu*phi)*X(r,z)

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'UPD_RESP_ODD(): time-odd response updated'//  &
      '  (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN

END SUBROUTINE upd_response_qtest


!=======================================================================

!     upd_response_even(P2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
!                       RESP_U,RESP_B,RESP_drB,RESP_dzB,RESP_Wz,RESP_Wr,


!     count the number of transitions (matrixelements) for a given
!     cranking mode

!=======================================================================

SUBROUTINE upd_response_even(p2q, phi,drphi,dzphi,dpphi,ddphi,occphi,  &
    resp_u, resp_b,resp_drb,resp_dpb,resp_dzb,  &
    resp_wr,resp_wp,resp_wz,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN)             :: p2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: phi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: occphi(mpsid)
REAL(8), INTENT(OUT)            :: resp_u(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_b(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_drb(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_dpb(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_dzb(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_wr(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_wp(0:mzd,0:mrd,2)
REAL(8), INTENT(OUT)            :: resp_wz(0:mzd,0:mrd,2)
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu



















REAL(8) :: resp_delta(0:mzd,0:mrd,2)

REAL(8) :: resp_rho(0:mzd,0:mrd,2)
REAL(8) :: resp_tau(0:mzd,0:mrd,2)
REAL(8) :: resp_socr(0:mzd,0:mrd,2)
REAL(8) :: resp_socp(0:mzd,0:mrd,2)
REAL(8) :: resp_socz(0:mzd,0:mrd,2)
REAL(8) :: resp_grhor(0:mzd,0:mrd,2)
REAL(8) :: resp_grhop(0:mzd,0:mrd,2)
REAL(8) :: resp_grhoz(0:mzd,0:mrd,2)
REAL(8) :: resp_drho(0:mzd,0:mrd,2)
REAL(8) :: resp_divsoc(0:mzd,0:mrd,2)
REAL(8) :: resp_pair(0:mzd,0:mrd,2)

INTEGER :: ist,jst,jtp
INTEGER :: itr,iz,ir
INTEGER :: iso
REAL(8) :: eta,xi,wgt1,wgt2,wgt3

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     set response densities to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO iz=0,nz
    DO ir=0,nr
      resp_rho(iz,ir,iso)    = 0.0D0
      resp_tau(iz,ir,iso)    = 0.0D0
      resp_socr(iz,ir,iso)   = 0.0D0
      resp_socp(iz,ir,iso)   = 0.0D0
      resp_socz(iz,ir,iso)   = 0.0D0
      resp_grhor(iz,ir,iso)  = 0.0D0
      resp_grhop(iz,ir,iso)  = 0.0D0
      resp_grhoz(iz,ir,iso)  = 0.0D0
      resp_drho(iz,ir,iso)   = 0.0D0
      resp_divsoc(iz,ir,iso) = 0.0D0
      resp_pair(iz,ir,iso)   = 0.0D0
    END DO
  END DO
END DO

!     loop all single particle states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        loop all transitions
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    jtp  = tridx(ist,itr,2)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    
!           if the final state is from the original space
    IF (jtp == +1) THEN
      
      IF (msgn == 0) THEN
        wgt1=   0.0D0                           ! for -ie^{i*mu*phi}
        wgt2=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      ELSE IF (msgn > 0) THEN
        wgt1=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for -ie^{i*mu*phi}
        wgt2=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      ELSE
        wgt1=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for -ie^{i*mu*phi}
        wgt2=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
        wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      END IF
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
          
          resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 0.5D0*wgt1*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
          
          resp_grhor(iz,ir,iso) = resp_grhor(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhop(iz,ir,iso) = resp_grhop(iz,ir,iso) + wgt1*(  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhoz(iz,ir,iso) = resp_grhoz(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
!                 2*resp_tau has to be added later
          resp_drho(iz,ir,iso) = resp_drho(iz,ir,iso) + wgt2*(  &
              +   phi(iz,ir,1,jst)*ddphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*ddphi(iz,ir,2,ist)  &
              + ddphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + ddphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(  &
              + drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
          
        END DO
      END DO
      
!           if the final state is from the time-reversed space
    ELSE
      
      wgt1=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for -ie^{i*mu*phi}
      wgt2=   2.0D0*(uj*vi + vj*ui)*p2q(itr2) ! for   e^{i*mu*phi}
      wgt3=   0.0D0*(uj*ui - vj*vi)*p2q(itr2)
      
      DO ir=0,nr
        DO iz=0,nz
          
          resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
          
          resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 0.5D0*wgt2*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 0.5D0*wgt1*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 0.5D0*wgt2*(  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
          
          resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
          
          resp_grhor(iz,ir,iso) = resp_grhor(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhop(iz,ir,iso) = resp_grhop(iz,ir,iso) + wgt1*(  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_grhoz(iz,ir,iso) = resp_grhoz(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
!                 2*resp_tau has to be added later
          resp_drho(iz,ir,iso) = resp_drho(iz,ir,iso) + wgt2*(  &
              -   phi(iz,ir,2,jst)*ddphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*ddphi(iz,ir,2,ist)  &
              - ddphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + ddphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
          
          resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(  &
              - drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
          
        END DO
      END DO
      
    END IF
    
  END DO    ! end of loop over all transitions from state ist
END DO       ! end of loop over all states ist

!     add 2*resp_tau to get the full Drho response
!     ------------------------------------------------------------------
DO iz=0,nz
  DO ir=0,nr
    resp_drho(iz,ir,1) = resp_drho(iz,ir,1) + 2.0D0*resp_tau(iz,ir,1)
    resp_drho(iz,ir,2) = resp_drho(iz,ir,2) + 2.0D0*resp_tau(iz,ir,2)
  END DO
END DO

!     recouple the Skyrme-Parameters
!     ------------------------------------------------------------------
b0  =                t0 + 0.5D0*t0*x0
b0p =          0.5D0*t0 +       t0*x0
b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
b2  = 0.125D0*(3.0D0*t1 + 1.5D0*t1*x1 -       t2 - 0.5D0*t2*x2)
b2p = 0.125D0*(1.5D0*t1 + 3.0D0*t1*x1 +   0.5*t2 +       t2*x2)
b3  = 0.250D0*(      t3 + 0.5D0*t3*x3)
b3p = 0.250D0*(0.5D0*t3 +       t3*x3)
b4  = 0.500D0*       t4
b5  = 0.125D0*(   x1*t1 +       x2*t2)  ! 2*b5
b5p = 0.125D0*(      t1 -          t2)  ! 2*b5'

wgt_b0  =  (b0-b0p)
wgt_b0c =   b0
wgt_b1  =  (b1-b1p)
wgt_b1c =   b1
wgt_b2  = -(b2-b2p)
wgt_b2c = - b2
IF (b4p == 0.0D0) THEN    ! no isotopic spin-orbit force
  wgt_b4  = - 2.0D0*b4
  wgt_b4c = -       b4
ELSE
  wgt_b4  = -(b4+b4p)
  wgt_b4c = - b4
END IF
IF(ifls == 1) THEN
  wgt_b5  =  -(b5-b5p)    ! b5 already includes the factor 2
  wgt_b5c =  - b5
ELSE
  wgt_b5  = 0.0D0
  wgt_b5c = 0.0D0
END IF

!     set the mean-fields to zero
!     ------------------------------------------------------------------
DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      resp_u(iz,ir,iso) = 0.0D0
      resp_b(iz,ir,iso) = 0.0D0
      resp_wr(iz,ir,iso) = 0.0D0
      resp_wp(iz,ir,iso) = 0.0D0
      resp_wz(iz,ir,iso) = 0.0D0
!            RESP_drB(iz,ir,iso) = 0.0D0
!            RESP_dpB(iz,ir,iso) = 0.0D0
!            RESP_dzB(iz,ir,iso) = 0.0D0
    END DO
  END DO
END DO

!     construct the mean-fields from the densities
!     ------------------------------------------------------------------

DO iso=1,2
  DO ir=0,nr
    DO iz=0,nz
      
!           the mean-field response without density-dependency (b3)
!           and coulomb
      resp_u(iz,ir,iso)  = resp_u(iz,ir,iso) + wgt_b0 *resp_rho(iz,ir,iso)  &
          + wgt_b0c*resp_rho(iz,ir,3-iso) + wgt_b1 *resp_tau(iz,ir,iso)  &
          + wgt_b1c*resp_tau(iz,ir,3-iso) + wgt_b2 *resp_drho(iz,ir,iso)  &
          + wgt_b2c*resp_drho(iz,ir,3-iso) + wgt_b4 *resp_divsoc(iz,ir,iso)  &
          + wgt_b4c*resp_divsoc(iz,ir,3-iso)
      
      resp_b(iz,ir,iso)  = resp_b(iz,ir,iso) + wgt_b1 *resp_rho(iz,ir,iso)  &
          + wgt_b1c*resp_rho(iz,ir,3-iso)
      
      resp_wr(iz,ir,iso) = resp_wr(iz,ir,iso)  &
          - wgt_b4 *resp_grhor(iz,ir,iso) - wgt_b4c*resp_grhor(iz,ir,3-iso)  &
          + wgt_b5 *resp_socr(iz,ir,iso) + wgt_b5c*resp_socr(iz,ir,3-iso)
      
      resp_wp(iz,ir,iso) = resp_wp(iz,ir,iso)  &
          - wgt_b4 *resp_grhop(iz,ir,iso) - wgt_b4c*resp_grhop(iz,ir,3-iso)  &
          + wgt_b5 *resp_socp(iz,ir,iso) + wgt_b5c*resp_socp(iz,ir,3-iso)
      
      resp_wz(iz,ir,iso) = resp_wz(iz,ir,iso)  &
          - wgt_b4 *resp_grhoz(iz,ir,iso) - wgt_b4c*resp_grhoz(iz,ir,3-iso)  &
          + wgt_b5 *resp_socz(iz,ir,iso) + wgt_b5c*resp_socz(iz,ir,3-iso)
      
!            RESP_drB(iz,ir,iso) = RESP_drB(iz,ir,iso)
!     &                         + wgt_b1 *resp_Grhor(iz,ir,iso)
!     &                         + wgt_b1c*resp_Grhor(iz,ir,3-iso)
!c
!            RESP_dpB(iz,ir,iso) = RESP_dpB(iz,ir,iso)
!     &                         + wgt_b1 *resp_Grhop(iz,ir,iso)
!     &                         + wgt_b1c*resp_Grhop(iz,ir,3-iso)
!c
!            RESP_dzB(iz,ir,iso) = RESP_dzB(iz,ir,iso)
!     &                         + wgt_b1 *resp_Grhoz(iz,ir,iso)
!     &                         + wgt_b1c*resp_Grhoz(iz,ir,3-iso)
!c
    END DO
  END DO
END DO

!     now the contribution from the density-dependency (b3)
!     ------------------------------------------------------------------
DO ir=0,nr
  DO iz=0,nz
    resp_u(iz,ir,1) = resp_u(iz,ir,1)  &
        + locpotwgt1(iz,ir,1)*resp_rho(iz,ir,1)  &
        + locpotwgt1(iz,ir,2)*resp_rho(iz,ir,2)
    resp_u(iz,ir,2) = resp_u(iz,ir,2)  &
        + locpotwgt2(iz,ir,1)*resp_rho(iz,ir,1)  &
        + locpotwgt2(iz,ir,2)*resp_rho(iz,ir,2)
  END DO
END DO

!     and from the coulomb
!     ------------------------------------------------------------------
IF (chg /= 0.0D0) THEN
  
END IF

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'UPD_RESP_EVEN(): time-even response updated'//  &
      '  (lambda=',il,' mu=',imu,')'
END IF

!     ------------------------------------------------------------------

RETURN

END SUBROUTINE upd_response_even


!=======================================================================

!     upd_Hresp_even(phi,drphi,dzphi,dpphi,ddphi,occphi,
!                   RESP_U,
!                   RESP_B,RESP_drB,RESP_dpB,RESP_dzB,
!                   RESP_Wr,RESP_Wp,RESP_Wz,il,imu)

!     evaluate the time-even response hamiltonian

!=======================================================================

SUBROUTINE upd_hresp_even(hresp2q, phi,drphi,dzphi,dpphi,ddphi,occphi,  &
    resp_u,resp_b, resp_wr,resp_wp,resp_wz,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: hresp2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: phi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: occphi(mpsid)
REAL(8), INTENT(IN OUT)         :: resp_u(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_b(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_wr(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_wp(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_wz(0:mzd,0:mrd,2)
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu

REAL(8) :: resp_delta(0:mzd,0:mrd,2)
REAL(8) :: resp_rho(0:mzd,0:mrd)
REAL(8) :: resp_tau(0:mzd,0:mrd)
REAL(8) :: resp_socr(0:mzd,0:mrd)
REAL(8) :: resp_socp(0:mzd,0:mrd)
REAL(8) :: resp_socz(0:mzd,0:mrd)
REAL(8) :: resp_pair(0:mzd,0:mrd)

INTEGER :: ist,jst,jtp
INTEGER :: itr,iz,ir
INTEGER :: iso
REAL(8) :: eta,xi,wgt1,wgt2

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     loop all single particle states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        loop all transitions
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    jtp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    xi   = uj*ui - vj*vi
    eta  = uj*vi + vj*ui
    
!           if the final state is from the original space
    IF (jtp == +1) THEN
      
      DO ir=0,nr
        DO iz=0,nz
             ! ~    exp(i*mu*phi)
          resp_rho(iz,ir) =               &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
           ! ~    exp(i*mu*phi)
          resp_tau(iz,ir) =                 &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
           ! ~    exp(i*mu*phi)
          resp_socr(iz,ir) = 0.5D0*(         &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist))
            ! ~  -i*exp(i*mu*phi)
          resp_socp(iz,ir) = 0.5D0*(       &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist))
            ! ~    exp(i*mu*phi)
          resp_socz(iz,ir) = 0.5D0*(       &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist))
           ! ~    exp(i*mu*phi) 
          resp_pair(iz,ir) =                &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
          
        END DO
      END DO
      
!           if the final state is from the time-reversed space
    ELSE
      
      DO ir=0,nr
        DO iz=0,nz
            ! ~    exp(i*mu*phi) 
          resp_rho(iz,ir) =               &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
            ! ~    exp(i*mu*phi)
          resp_tau(iz,ir) =                &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
           ! ~    exp(i*mu*phi)
          resp_socr(iz,ir) = 0.5D0*(       &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist))
           ! ~  i*exp(i*mu*phi)
          resp_socp(iz,ir) = 0.5D0*(       &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist))
            ! ~    exp(i*mu*phi)
          resp_socz(iz,ir) = 0.5D0*(      &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist))
           ! ~    exp(i*mu*phi) 
          resp_pair(iz,ir) =                &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
          
        END DO
      END DO
      
    END IF
    
!          evaluate the matrixelement of the transition
!          -------------------------------------------------------------
    IF (imu == 0) THEN
      hresp2q(itr2) = eta*( + f2int(resp_u(0:nz,0:nr,iso) ,resp_rho(0:nz,0:nr))  &
          + f2int(resp_b(0:nz,0:nr,iso) ,resp_tau(0:nz,0:nr))  &
          + f2int(resp_wr(0:nz,0:nr,iso),resp_socr(0:nz,0:nr))  &
          + f2int(resp_wz(0:nz,0:nr,iso),resp_socz(0:nz,0:nr)) )
    ELSE
      IF (msgn > 0) THEN
        hresp2q(itr2) = 0.5D0*eta*( + f2int(resp_u(0:nz,0:nr,iso) ,resp_rho(0:nz,0:nr))  &
            + f2int(resp_b(0:nz,0:nr,iso) ,resp_tau(0:nz,0:nr))  &
            + f2int(resp_wr(0:nz,0:nr,iso),resp_socr(0:nz,0:nr))  &
            + f2int(resp_wz(0:nz,0:nr,iso),resp_socz(0:nz,0:nr))  &
            - f2int(resp_wp(0:nz,0:nr,iso),resp_socp(0:nz,0:nr)) )
      ELSE
        hresp2q(itr2) =  0.5D0*eta*( + f2int(resp_u(0:nz,0:nr,iso) ,resp_rho(0:nz,0:nr))  &
            + f2int(resp_b(0:nz,0:nr,iso) ,resp_tau(0:nz,0:nr))  &
            + f2int(resp_wr(0:nz,0:nr,iso),resp_socr(0:nz,0:nr))  &
            + f2int(resp_wz(0:nz,0:nr,iso),resp_socz(0:nz,0:nr))  &
            + f2int(resp_wp(0:nz,0:nr,iso),resp_socp(0:nz,0:nr)) )
      END IF
    END IF
    
  END DO    ! end of loop over all transitions from state ist
END DO       ! end of loop over all states ist

RETURN
END SUBROUTINE upd_hresp_even


!=======================================================================

!     upd_Hresp_odd(phi,drphi,dzphi,dpphi,ddphi,occphi,
!                   RESP_Ar,RESP_Ap,RESP_Az,
!                   RESP_Sr,RESP_Sp,RESP_Sz,
!                   RESP_Cr,RESP_Cp,RESP_Cz,il,imu)

!     evaluate the time-odd response hamiltonian

!=======================================================================

SUBROUTINE upd_hresp_odd(hresp2q, phi,drphi,dzphi,dpphi,ddphi,occphi,  &
    resp_ar,resp_ap,resp_az, resp_sr,resp_sp,resp_sz,  &
    resp_cr,resp_cp,resp_cz,il,imu)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: hresp2q(1:mtdtot)
REAL(8), INTENT(IN OUT)         :: phi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: drphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dzphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: dpphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: ddphi(0:mzd,0:mrd,msd,mpsid)
REAL(8), INTENT(IN OUT)         :: occphi(mpsid)
REAL(8), INTENT(IN OUT)         :: resp_ar(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_ap(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_az(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_sr(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_sp(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_sz(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_cr(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_cp(0:mzd,0:mrd,2)
REAL(8), INTENT(IN OUT)         :: resp_cz(0:mzd,0:mrd,2)
INTEGER, INTENT(OUT)                     :: il
INTEGER, INTENT(IN)                      :: imu




















REAL(8) :: resp_delta(0:mzd,0:mrd,2)

REAL(8) :: resp_spnr(0:mzd,0:mrd)
REAL(8) :: resp_spnp(0:mzd,0:mrd)
REAL(8) :: resp_spnz(0:mzd,0:mrd)
REAL(8) :: resp_skdr(0:mzd,0:mrd)
REAL(8) :: resp_skdp(0:mzd,0:mrd)
REAL(8) :: resp_skdz(0:mzd,0:mrd)
REAL(8) :: resp_curr(0:mzd,0:mrd)
REAL(8) :: resp_curp(0:mzd,0:mrd)
REAL(8) :: resp_curz(0:mzd,0:mrd)
REAL(8) :: resp_paio(0:mzd,0:mrd)
REAL(8) :: resp_rotjr(0:mzd,0:mrd)
REAL(8) :: resp_rotjp(0:mzd,0:mrd)
REAL(8) :: resp_rotjz(0:mzd,0:mrd)
REAL(8) :: resp_rotsr(0:mzd,0:mrd)
REAL(8) :: resp_rotsp(0:mzd,0:mrd)
REAL(8) :: resp_rotsz(0:mzd,0:mrd)

INTEGER :: ist,jst,jtp
INTEGER :: itr,iz,ir
INTEGER :: iso,msgn
REAL(8) :: eta,xi,wgt1,wgt2

!     buffered SQRT as internal statement function
!     ------------------------------------------------------------------
!sqrtb(x) = DSQRT(MAX(x,smal))

!     ------------------------------------------------------------------

!     check parameters
!     ------------------------------------------------------------------
IF (imu < 0) THEN
  WRITE(21,'(a,i2)') 'ERROR: mu-component must be >=0 -- not ',imu
  STOP 'in linresp'
ELSE IF (imu > il) THEN
  WRITE(21,'(a,i2,a,i2)') 'ERROR: mu-component must be <=il (=',il,  &
      ')  -- not ',imu
  STOP 'in linresp'
ELSE IF (il > lambdamax) THEN
  WRITE(21,'(a,i2,a,i2,a)') 'ERROR: LAMBDAMAX too small (il=',il,  &
      ' LAMBDAMAX=',il,')'
  STOP 'in linresp'
END IF

!     loop all single particle states
!     ------------------------------------------------------------------
DO ist=1,nst(1)+nst(2)
  
  IF (ist <= nst(1)) THEN
    iso    = 1
  ELSE
    iso    = 2
  END IF
  
!        loop all transitions
!        ---------------------------------------------------------------
  DO itr=1,trntr(ist)
    
    jst  = tridx(ist,itr,1)
    jtp  = tridx(ist,itr,2)
    msgn = tridx(ist,itr,3)
    itr2 = tridx(ist,itr,4)
    
    vi   = sqrtb(cpvsqf(ist,0))
    ui   = sqrtb(one-vi*vi)
    vj   = sqrtb(cpvsqf(jst,0))
    uj   = sqrtb(one-vj*vj)
    
    xi   = uj*ui + vj*vi
    eta  = uj*vi - vj*ui
    
!           if the final state is from the original space
    IF (jtp == +1) THEN
      
      DO ir=0,nr
        DO iz=0,nz
            ! ~    exp(i*mu*phi)
          resp_spnr(iz,ir) =               &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
            ! ~ -i*exp(i*mu*phi)
          resp_spnp(iz,ir) =               &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)  &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
            ! ~    exp(i*mu*phi)
          resp_spnz(iz,ir) =               &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
          
           ! ~    exp(i*mu*phi) 
          resp_skdr(iz,ir) =               &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
           ! ~ -i*exp(i*mu*phi)
          resp_skdp(iz,ir) =                &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
           ! ~    exp(i*mu*phi) 
          resp_skdz(iz,ir) =               &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              - dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
          
           ! ~ -i*exp(i*mu*phi) 
          resp_curr(iz,ir) = 0.5D0*(              &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist))
           ! ~    exp(i*mu*phi)
          resp_curp(iz,ir) = 0.5D0*(               &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist))
           ! ~ -i*exp(i*mu*phi) 
          resp_curz(iz,ir) = 0.5D0*(              &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist))
          
           ! ~  i*exp(i*mu*phi)
          resp_paio(iz,ir) =                &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
          
        END DO
      END DO
      
!           if the final state is from the time-reversed space
    ELSE
      
      DO ir=0,nr
        DO iz=0,nz
           ! ~    exp(i*mu*phi)
          resp_spnr(iz,ir) =                &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
           ! ~ -i*exp(i*mu*phi)
          resp_spnp(iz,ir) =                &
              - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)  &
              - phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
           ! ~    exp(i*mu*phi)
          resp_spnz(iz,ir) =                &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              - phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
          
           ! ~    exp(i*mu*phi)
          resp_skdr(iz,ir) =                &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              - dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
           ! ~ -i*exp(i*mu*phi)
          resp_skdp(iz,ir) =                &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)  &
              - drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
           ! ~    exp(i*mu*phi)
          resp_skdz(iz,ir) =                &
              - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
          
           ! ~ -i*exp(i*mu*phi)
          resp_curr(iz,ir) = 0.5D0*(               &
              -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)  &
              + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist))
           ! ~    exp(i*mu*phi) 
          resp_curp(iz,ir) = 0.5D0*(              &
              -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)  &
              + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist))
            ! ~ -i*exp(i*mu*phi)
          resp_curz(iz,ir) = 0.5D0*(              &
              -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)  &
              +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)  &
              + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)  &
              - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist))
          
            ! ~  i*exp(i*mu*phi)
          resp_paio(iz,ir) =               &
              - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)  &
              + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
          
        END DO
      END DO
      
    END IF
    
!          -------------------------------------------------------------
    IF (imu == 0) THEN
      hresp2q(itr2) = eta*( - f2int(resp_ar(0:nz,0:nr,iso),resp_curr(0:nz,0:nr))  &
          - f2int(resp_az(0:nz,0:nr,iso),resp_curz(0:nz,0:nr))  &
          - f2int(resp_sp(0:nz,0:nr,iso),resp_spnp(0:nz,0:nr))  &
          - f2int(resp_cp(0:nz,0:nr,iso),resp_skdp(0:nz,0:nr)) )
    ELSE
      IF (msgn > 0) THEN
        hresp2q(itr2) = 0.5D0*eta*( - f2int(resp_ar(0:nz,0:nr,iso),resp_curr(0:nz,0:nr))  &
            - f2int(resp_az(0:nz,0:nr,iso),resp_curz(0:nz,0:nr))  &
            - f2int(resp_sp(0:nz,0:nr,iso),resp_spnp(0:nz,0:nr))  &
            - f2int(resp_cp(0:nz,0:nr,iso),resp_skdp(0:nz,0:nr))  &
            + f2int(resp_ap(0:nz,0:nr,iso),resp_curp(0:nz,0:nr))  &
            + f2int(resp_sr(0:nz,0:nr,iso),resp_spnr(0:nz,0:nr))  &
            + f2int(resp_sz(0:nz,0:nr,iso),resp_spnz(0:nz,0:nr))  &
            + f2int(resp_cr(0:nz,0:nr,iso),resp_skdr(0:nz,0:nr))  &
            + f2int(resp_cz(0:nz,0:nr,iso),resp_skdz(0:nz,0:nr)) )
      ELSE
        hresp2q(itr2) = 0.5D0*eta*( - f2int(resp_ar(0:nz,0:nr,iso),resp_curr(0:nz,0:nr))  &
            - f2int(resp_az(0:nz,0:nr,iso),resp_curz(0:nz,0:nr))  &
            - f2int(resp_sp(0:nz,0:nr,iso),resp_spnp(0:nz,0:nr))  &
            - f2int(resp_cp(0:nz,0:nr,iso),resp_skdp(0:nz,0:nr))  &
            - f2int(resp_ap(0:nz,0:nr,iso),resp_curp(0:nz,0:nr))  &
            - f2int(resp_sr(0:nz,0:nr,iso),resp_spnr(0:nz,0:nr))  &
            - f2int(resp_sz(0:nz,0:nr,iso),resp_spnz(0:nz,0:nr))  &
            - f2int(resp_cr(0:nz,0:nr,iso),resp_skdr(0:nz,0:nr))  &
            - f2int(resp_cz(0:nz,0:nr,iso),resp_skdz(0:nz,0:nr)) )
      END IF
    END IF
    
  END DO    ! end of loop over all transitions from state ist
END DO       ! end of loop over all states ist

RETURN
END SUBROUTINE upd_hresp_odd




!=======================================================================

!     upd_locpotwgt()

!       the response of the non linear contributions to the local
!       mean-field potential U is linearized

!              U^{NL}_q = U^{NL}_{q,0} +
!                         U^{NL}_{q,resp,p} rerho_p +
!                         U^{NL}_{q,resp,n} rerho_n

!       with

!              U^{NL}_{q,resp,q'} = dU_q/drho_q'

!       this function calculates the response wheight dU_q/drho_q' and
!       stores its values in locpotwgt:

!          locpotwgt1(0,0,1) = dU_p/drho_p
!          locpotwgt1(0,0,2) = dU_p/drho_n
!          locpotwgt2(0,0,1) = dU_n/drho_p
!          locpotwgt2(0,0,2) = dU_n/drho_n

!=======================================================================

SUBROUTINE upd_locpotwgt()

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'linresp.inc'

!     ------------------------------------------------------------------

REAL(8) :: fac1,fac2,fac3,fac4,fac5
REAL(8) :: buff1,buff2,buff3,buff4,buff5,buff6,buff7,buff8
REAL(8) :: wgtfac
INTEGER :: ir,iz,ist

!     ------------------------------------------------------------------

!     calculate densities rho             --> LRwrk1_spinor
!     ------------------------------------------------------------------
DO ir=0,nr
  DO iz=0,nz
    lrwrk1_spinor(iz,ir,1) = 0.0D0
    lrwrk1_spinor(iz,ir,2) = 0.0D0
  END DO
END DO

DO ist=1,nst(1)
  wgtfac = 2.0D0*cpvsqf(ist,0)
  DO ir=0,nr
    DO iz=0,nz
      lrwrk1_spinor(iz,ir,1) = lrwrk1_spinor(iz,ir,1)+  &
          wgtfac*(cppsi(iz,ir,1,ist, 0)*cppsi(iz,ir,1,ist, 0) +  &
          cppsi(iz,ir,2,ist, 0)*cppsi(iz,ir,2,ist, 0))
    END DO
  END DO
END DO

DO ist=nst(1)+1,nst(2)
  wgtfac = 2.0D0*cpvsqf(ist,0)
  DO ir=0,nr
    DO iz=0,nz
      lrwrk1_spinor(iz,ir,2) = lrwrk1_spinor(iz,ir,2)+  &
          wgtfac*(cppsi(iz,ir,1,ist, 0)*cppsi(iz,ir,1,ist, 0) +  &
          cppsi(iz,ir,2,ist, 0)*cppsi(iz,ir,2,ist, 0))
    END DO
  END DO
END DO

!     calculate response weight dU/drho   --> locpotwgt
!     ------------------------------------------------------------------
fac1 =  b3 *(power+2.0D0)*(power+1.0D0)/3.0D0
fac2 = -b3p*(2.0D0/3.0D0)
fac3 = -b3p*(2.0D0*power)/3.0D0
fac4 = -b3p*power*(power-1.0D0)/3.0D0

IF(chgx /= zero) THEN
  fac5 = 1.418033D0*chgx/3.0D0
ELSE
  fac5 = 0.0D0
END IF

DO ir=0,nr
  DO iz=0,nz
    locpotwgt1(iz,ir,1)=0.0D0
    locpotwgt1(iz,ir,2)=0.0D0
    locpotwgt2(iz,ir,1)=0.0D0
    locpotwgt2(iz,ir,2)=0.0D0
  END DO
END DO

DO ir=0,nr
  DO iz=0,nz
    buff1 = lrwrk1_spinor(iz,ir,1) ! rho_p
    buff2 = lrwrk1_spinor(iz,ir,2) ! rho_n
    buff3 = buff1 + buff2          ! rho
    IF (buff3 <= smal) THEN
      buff4 = 0.0D0               ! rho**alpha
      buff5 = 1.0D0               ! (rho_p^2 + rho_n^2)/rho^2
      buff6 = 1.0D0               ! rho_p/rho
      buff7 = 1.0D0               ! rho_n/rho
      buff8 = 0.0D0
    ELSE
      buff4 =  buff3**(power)
      buff5 = (buff1**2 + buff2**2)/buff3**2
      buff6 = (buff1/buff3)
      buff7 = (buff2/buff3)
      buff8 = 0.0D0
      buff8 = buff1**(-0.6666666666666667D0)
    END IF
    
    
!           dU_p/drho_p
    locpotwgt1(iz,ir,1)= buff4*(fac1 + fac2 +  2.0D0*fac3        *buff6  &
        +        fac4        *buff5) +buff8*fac5
    
!           dU_n/drho_n
    locpotwgt2(iz,ir,2)= buff4*(fac1 + fac2 +  2.0D0*fac3        *buff7  &
        +        fac4        *buff5)
    
!           dU_p/drho_n
    locpotwgt1(iz,ir,2)= buff4*(fac1        +       fac3  &
        +       fac4         *buff5)
    
!           dU_n/drho_p
    locpotwgt2(iz,ir,1)= buff4*(fac1        +       fac3  &
        +       fac4         *buff5)
    
  END DO
END DO

!     ------------------------------------------------------------------

IF (d_lrdebug) THEN
  WRITE(21,'(a)') 'UPD_LOCPOTWGT() : response weight dU/drho updated'
END IF

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE upd_locpotwgt



!=======================================================================

!     ini_linrespfiles()

!=======================================================================

SUBROUTINE ini_linrespfiles()

!     ------------------------------------------------------------------

USE paramr
USE constr
USE dynpar
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'
!INCLUDE 'constr.inc'
!INCLUDE 'linresp.inc'


!CHARACTER (LEN=78), INTENT(IN OUT)       :: title
!CHARACTER (LEN=3), INTENT(IN OUT)        :: outnam
INTEGER :: i



!     ------------------------------------------------------------------

OPEN(21,STATUS='unknown',FORM='formatted', FILE='linresp.'//outnam)
OPEN(22,STATUS='unknown',FORM='formatted', FILE='LRiter.'//outnam)
OPEN(29,STATUS='unknown',FORM='formatted', FILE='radii.'//outnam)
OPEN(30,STATUS='unknown',FORM='formatted', FILE='width.'//outnam)
OPEN(32,STATUS='unknown',FORM='formatted', FILE='norm.'//outnam)

WRITE(29,'(a/a,2i4/a,i1,3a)') title(1:78),  &
    ' A,Z=',npart(1)+npart(2),npart(1), ' beta_',lqcon,'   <P^2_beta>',  &
    '     r             R             sigma  ',  &
    '     (d/db)^2 r    (d/db)^2 R    (d/db)^2 sigma'

!     Kopfzeile der beiben Dateien norm und width vertauscht 30<>32

WRITE(32,'(2a/4a/4a)')  &
    '                                 qp-contributions           ',  &
    '                                 residual contributions     ',  &
    '             norm <[P,Beta]>                   ',  &
    'width  <P^2>             ',  &
    '             norm <[P,Beta]>                   ', 'width  <P^2> ',  &
    ' beta    mu= 0         1            2         ',  &
    ' 0           1            2 ',  &
    '       mu= 0         1            2         ', ' 0           1            2 '


WRITE(30,'(a/a,2i4/2a/2a)') title(1:78),  &
    ' A,Z=',npart(1)+npart(2),npart(1),  &
    '                      <P^2>                ',  &
    '                      <Q^2>',  &
    ' beta     mu= 0         1              2   ',  &
    '          mu= 0         1             2 '


WRITE(21,'(a)') ' '
WRITE(22,'(a)') ' '

OPEN(23,STATUS='unknown', FORM='formatted',FILE='colpar_mp.'//outnam)
OPEN(24,STATUS='unknown', FORM='formatted',FILE='colpar_GCM.'//outnam)
OPEN(25,STATUS='unknown', FORM='formatted',FILE='colpar_in.'//outnam)
OPEN(26,STATUS='unknown', FORM='formatted',FILE='colpar_sc.'//outnam)
OPEN(27,STATUS='unknown', FORM='formatted',FILE='colpar_CHF.'//outnam)
OPEN(28,STATUS='unknown', FORM='formatted',FILE='colpar_at.'//outnam)


WRITE(23,'(a)') ' parameters of quad.-dynamics (scaling mass)'
WRITE(24,'(a)') ' parameters of quad.-dynamics (static GCM)'
WRITE(25,'(a)') ' parameters of quad.-dynamics (inglis-cranking)'
WRITE(26,'(a)') ' parameters of quad.-dynamics (SC-cranking)'
WRITE(27,'(a)') ' parameters of quad.-dynamics (CHF + SC-Cranking)'
WRITE(28,'(a)') ' parameters of quad.-dynamics (full ATDHF)'


DO i=23,28
  WRITE(i,'(a/a,2i4/7a)') title(1:78),  &
      ' A,Z=',npart(1)+npart(2),npart(1),  &
      ' beta_20    ',  &
      'B_20        C_20    <P_20^2>    <Q_20^2>', '  <Omg_20^2>        ',  &
      'B_21        C_21    <P_21^2>    <Q_21^2>', '  <Omg_21^2>        ',  &
      'B_22        C_22    <P_22^2>    <Q_22^2>', '  <Omg_22^2>'
END DO


!     ------------------------------------------------------------------

RETURN
END SUBROUTINE ini_linrespfiles




!=======================================================================

!     CPshift()

!=======================================================================

SUBROUTINE cpshift(cpval,newval)

!     ------------------------------------------------------------------


REAL(8), INTENT(OUT)            :: cpval(-1:1)
REAL(8), INTENT(IN)             :: newval



!     ------------------------------------------------------------------
!     shift values in path storage CPpsi (-1,0,+1) -> (0,+1,+1)
!     and copy current value to last path storage index
!     ------------------------------------------------------------------
cpval(-1) = cpval( 0)
cpval( 0) = cpval(+1)
cpval(+1) = newval

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE cpshift


!=======================================================================

!     static_afterburn()

!=======================================================================

SUBROUTINE static_afterburn(ipath,varerr,maxiter)

!     ------------------------------------------------------------------

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'linresp.inc'
!INCLUDE 'dynpar.inc'


!     ------------------------------------------------------------------


INTEGER, INTENT(IN)         :: ipath
REAL(8), INTENT(IN)         :: varerr
INTEGER, INTENT(IN)         :: maxiter




INTEGER :: ir,iz,im
INTEGER :: ist,jst
REAL(8) :: energy(mpsid),var(mpsid),norm(mpsid)
REAL(8) :: maxvar,lastmaxvar,dmaxvar
REAL(8) :: hpsi(0:mzd,0:mrd,msd,mpsid)
REAL(8) :: hhpsi(0:mzd,0:mrd,msd,mpsid)
REAL(8) :: dzhpsi(0:mzd,0:mrd,msd,mpsid)
REAL(8) :: drhpsi(0:mzd,0:mrd,msd,mpsid)
REAL(8) :: dphpsi(0:mzd,0:mrd,msd,mpsid)
REAL(8) :: ddhpsi(0:mzd,0:mrd,msd,mpsid)
REAL(8) :: buffer,normstep,normerror
INTEGER :: touched,nlo,nhi,nadd,kpsic,iq
INTEGER :: subiter

REAL(8) :: mfd_rho(0:mzd,0:mrd,2)
REAL(8) :: mfd_tau(0:mzd,0:mrd,2)
REAL(8) :: mfd_jz(0:mzd,0:mrd,2)
REAL(8) :: mfd_jr(0:mzd,0:mrd,2)
REAL(8) :: mfd_grhor(0:mzd,0:mrd,2)
REAL(8) :: mfd_grhoz(0:mzd,0:mrd,2)
REAL(8) :: mfd_drho(0:mzd,0:mrd,2)
REAL(8) :: mfd_divj(0:mzd,0:mrd,2)

REAL(8) :: mfp_u(0:mzd,0:mrd,2)
REAL(8) :: mfp_unl(0:mzd,0:mrd,2)
REAL(8) :: mfp_b(0:mzd,0:mrd,2)
REAL(8) :: mfp_drb(0:mzd,0:mrd,2)
REAL(8) :: mfp_dzb(0:mzd,0:mrd,2)
REAL(8) :: mfp_wz(0:mzd,0:mrd,2)
REAL(8) :: mfp_wr(0:mzd,0:mrd,2)

INTEGER :: istart(2),iend(2)
REAL(8) :: occ1,occ2

!     ------------------------------------------------------------------

WRITE(21,'(a,g12.4,i8)')  &
    'STATIC_AFTERBURN(): improving sp-base: max. var, max. iter:' ,varerr,maxiter

!     calculate derivatives of single particle wavefunctions
!     ------------------------------------------------------------------
CALL upd_derive_spinor(cppsi(0,0,1,1,ipath),  &
    cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1), cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),  &
    +1,0,nst(1)+nst(2))

!     accumulate the densities and mean-fields
!     ------------------------------------------------------------------
!      call upd_static_MF(CPpsi(0,0,1,1,ipath),
!     &                   CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
!     &                   CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),
!     &                   CPvsqf(1,ipath),
!     &                   MFP_U,MFP_B,MFP_drB,MFP_dzB,
!     &                   MFP_Wz,MFP_Wr)


DO ir=0,nr
  DO iz=0,nz
    cpeffmas(iz,ir,1,ipath) = cpeffmas(iz,ir,1,ipath)-h2meff(1)
    cpeffmas(iz,ir,2,ipath) = cpeffmas(iz,ir,2,ipath)-h2meff(2)
  END DO
END DO

CALL nablaz_scalar(cpeffmas(0:nz,0:nr,1,ipath), mfp_dzb(0:nz,0:nr,1),1,1,one)
CALL nablar_scalar(cpeffmas(0:nz,0:nr,1,ipath), mfp_drb(0:nz,0:nr,1),0,1,one)
CALL nablaz_scalar(cpeffmas(0:nz,0:nr,2,ipath), mfp_dzb(0:nz,0:nr,2),1,1,one)
CALL nablar_scalar(cpeffmas(0:nz,0:nr,2,ipath), mfp_drb(0:nz,0:nr,2),0,1,one)

DO ir=0,nr
  DO iz=0,nz
    cpeffmas(iz,ir,1,ipath) = cpeffmas(iz,ir,1,ipath)+h2meff(1)
    cpeffmas(iz,ir,2,ipath) = cpeffmas(iz,ir,2,ipath)+h2meff(2)
  END DO
END DO

DO ir=0,nr
  DO iz=0,nz
    mfp_u(iz,ir,1)          = cplocpot(iz,ir,1,ipath)
    mfp_u(iz,ir,2)          = cplocpot(iz,ir,2,ipath)
    mfp_b(iz,ir,1)          = cpeffmas(iz,ir,1,ipath)
    mfp_b(iz,ir,2)          = cpeffmas(iz,ir,2,ipath)
    mfp_wr(iz,ir,1)         = cpsocpotr(iz,ir,1,ipath)
    mfp_wr(iz,ir,2)         = cpsocpotr(iz,ir,2,ipath)
    mfp_wz(iz,ir,1)         = cpsocpotz(iz,ir,1,ipath)
    mfp_wz(iz,ir,2)         = cpsocpotz(iz,ir,2,ipath)
  END DO
END DO

!     calculate error in Schroedinger-Eq.
!     ------------------------------------------------------------------
CALL act_mfh_even(cppsi(0,0,1,1,ipath),cpdrpsi(0,0,1,1),  &
    cpdzpsi(0,0,1,1),cpdppsi(0,0,1,1), cpddpsi(0,0,1,1),hpsi(0,0,1,1),  &
    mfp_u,                                        1,  &
    mfp_b,mfp_drb,mfp_dzb,    mfp_dzb,            1,  &
    mfp_wr,mfp_wz,            mfp_wz,             1,  &
    cpconstr(0,0,ipath),cplambda(ipath),          1,  &
    cpspenrg(1,ipath),                            0, 0,1.0D0,  &
    1,nst(1)+nst(2))

CALL upd_derive_spinor(hpsi(0,0,1,1), drhpsi(0,0,1,1),dzhpsi(0,0,1,1),  &
    dphpsi(0,0,1,1),ddhpsi(0,0,1,1), +1,0,nst(1)+nst(2))

CALL act_mfh_even(hpsi(0,0,1,1),drhpsi(0,0,1,1),  &
    dzhpsi(0,0,1,1),dphpsi(0,0,1,1), ddhpsi(0,0,1,1),hhpsi(0,0,1,1),  &
    mfp_u,                                        1,  &
    mfp_b,mfp_drb,mfp_dzb,    mfp_dzb,            1,  &
    mfp_wr,mfp_wz,            mfp_wz,             1,  &
    cpconstr(0,0,ipath),cplambda(ipath),          1,  &
    cpspenrg(1,ipath),                            0, 0,1.0D0,  &
    1,nst(1)+nst(2))


!     calculate variance and energy expectation
!     ------------------------------------------------------------------
maxvar=0.0D0
DO ist=1,nst(1)+nst(2)
  energy(ist) = wfovlp(cppsi(0:nz,0:nr,1:2,ist,ipath), hpsi(0:nz,0:nr,1:2,ist))
  var(ist)    = wfovlp(cppsi(0:nz,0:nr,1:2,ist,ipath),hhpsi(0:nz,0:nr,1:2,ist))
!         var(ist)    = wfovlp(HPsi(0,0,1,ist),HPsi(0,0,1,ist))
  var(ist)    = var(ist)-energy(ist)**2
  maxvar=MAX(maxvar,ABS(var(ist)))
  
!              write(88,'(a,3(1pg16.4),3i3)') 'PKtest1:',
!     &              energy(ist),var(ist),
!     &              CPspenrg(ist,ipath),
!     &              mpsi(ist),nparz(1,ist),icharg(ist)
  
END DO

WRITE(21,'(a,g12.3)')  &
    'STATIC_AFTERBURN(): started   ---  max. var       : ', maxvar

!     the static iteration (with fixed mean-fields, occ. numbers, etc.)
!     ------------------------------------------------------------------
lastmaxvar=0.0D0
dmaxvar   =1E5
DO iter=1,maxiter
  
  IF (maxvar < varerr.OR.iter == maxiter) THEN
    DO ist=1,nst(1)+nst(2)
      energy(ist) = wfovlp(cppsi(0:nz,0:nr,1:2,ist,ipath), hpsi(0:nz,0:nr,1:2,ist))
      var(ist)    = wfovlp(cppsi(0:nz,0:nr,1:2,ist,ipath), hhpsi(0:nz,0:nr,1:2,ist))
      var(ist)    = var(ist)-energy(ist)**2
      
!              write(88,'(a,3(1pg16.4),3i3)') 'PKtest:',
!     &              energy(ist),var(ist),
!     &              CPspenrg(ist,ipath),
!     &              mpsi(ist),nparz(1,ist),icharg(ist)
    END DO
    
    WRITE(21,'(a,g12.3,i8)')  &
        'STATIC_AFTERBURN(): finished  ---  max. var, niter: ', maxvar,iter
    GO TO 4000
  END IF
  
  dmaxvar=DABS(maxvar-lastmaxvar)
  lastmaxvar=maxvar
  IF (MOD(iter,100) == 0) THEN
    WRITE(21,'(a,i6,2g12.2)') 'afterburn: ',iter,maxvar,dmaxvar
  END IF
  IF (MOD(iter,20) == 0)  &
      WRITE(6,'(a,i6,2g12.2)') 'afterburn: ',iter,maxvar,dmaxvar
  
!        loop all m-blocks
!        ---------------------------------------------------------------
  DO im=1,nom
    touched = 0
    nhi     = mhi(im)
    nlo     = mlo(im)
    nadd    = nhi-nlo+1
    kpsic   = mim(im)
    iq      = miq(im)
    
!           now the gradient step
!           ------------------------------------------------------------
    DO ist=nlo,nhi
      IF (ABS(var(ist)) > varerr) THEN
        touched = 1
      END IF
    END DO
    
    IF (touched == 1) THEN
      
      DO ist=nlo,nhi
        
!                 calculate error in sp-states
!                 ---------------------------------------------
        DO ir=0,nr
          DO iz=0,nz
            hhpsi(iz,ir,1,ist)=hpsi(iz,ir,1,ist)  &
                -energy(ist)*cppsi(iz,ir,1,ist,ipath)
            hhpsi(iz,ir,2,ist)=hpsi(iz,ir,2,ist)  &
                -energy(ist)*cppsi(iz,ir,2,ist,ipath)
          END DO
        END DO
        
!                 apply damping
!                 ---------------------------------------------
        CALL damnuc(hhpsi(0:nz,0:nr,1:2,ist),hpsi(0:nz,0:nr,1:2,ist),ist)
        
      END DO
      
      DO ist=nlo,nhi
        
!                 perform the gradient step
!                 ---------------------------------------------
        DO ir=0,nr
          DO iz=0,nz
            cppsi(iz,ir,1,ist,ipath)= cppsi(iz,ir,1,ist,ipath)  &
                -hpsi(iz,ir,1,ist)
            cppsi(iz,ir,2,ist,ipath)= cppsi(iz,ir,2,ist,ipath)  &
                -hpsi(iz,ir,2,ist)
          END DO
        END DO
        
      END DO
      
!              orthonormalize m-block (only if touched=1)
!              --------------------------------------------
      CALL schmit(cppsi(0:nz,0:nr,1:2,nlo:nhi,ipath),nlo-1,nadd)
      
!              update derivatives
!              --------------------------------------------
      CALL upd_derive_spinor_range(cppsi(0,0,1,1,ipath),  &
          cpdrpsi(0,0,1,1),cpdzpsi(0,0,1,1), cpdppsi(0,0,1,1),cpddpsi(0,0,1,1),  &
          +1,0,nlo,nhi)
      
!              calculate H|psi>
!              ---------------------------------------------
      CALL act_mfh_even(cppsi(0,0,1,1,ipath),cpdrpsi(0,0,1,1),  &
          cpdzpsi(0,0,1,1),cpdppsi(0,0,1,1), cpddpsi(0,0,1,1),hpsi(0,0,1,1),  &
          mfp_u,                                        1,  &
          mfp_b,mfp_drb,mfp_dzb,mfp_dzb,                1,  &
          mfp_wr,mfp_wz,            mfp_wz,             1,  &
          cpconstr(0,0,ipath),cplambda(ipath),          1,  &
          cpspenrg(1,ipath),                            0, 0,1.0D0,  &
          nlo,nhi)
      
!              update derivatives
!              --------------------------------------------
      CALL upd_derive_spinor_range(hpsi(0,0,1,1),  &
          drhpsi(0,0,1,1),dzhpsi(0,0,1,1), dphpsi(0,0,1,1),ddhpsi(0,0,1,1),  &
          +1,0,nlo,nhi)
      
!              calculate HH|psi>
!              ---------------------------------------------
      CALL act_mfh_even(hpsi(0,0,1,1),drhpsi(0,0,1,1),  &
          dzhpsi(0,0,1,1),dphpsi(0,0,1,1), ddhpsi(0,0,1,1),hhpsi(0,0,1,1),  &
          mfp_u,                                        1,  &
          mfp_b,mfp_drb,mfp_dzb,mfp_dzb,                1,  &
          mfp_wr,mfp_wz,            mfp_wz,             1,  &
          cpconstr(0,0,ipath),cplambda(ipath),          1,  &
          cpspenrg(1,ipath),                            0, 0,1.0D0,  &
          nlo,nhi)
      
    END IF
    
  END DO
  
  maxvar   =0.0D0
  DO ist=1,nst(1)+nst(2)
    energy(ist) = wfovlp(cppsi(0:nz,0:nr,1:2,ist,ipath), hpsi(0:nz,0:nr,1:2,ist))
    var(ist)    = wfovlp(cppsi(0:nz,0:nr,1:2,ist,ipath), hhpsi(0:nz,0:nr,1:2,ist))
    var(ist)    = var(ist)-energy(ist)**2
    maxvar=MAX(maxvar,ABS(var(ist)))
    
  END DO
  
END DO

4000 CONTINUE

DO ist=1,nst(1)+nst(2)
  cpspenrg(ist,ipath) = energy(ist)
END DO

!      do ist=1,nst(1)+nst(2)
!         write(88,'(a,3(1pg16.4),3i3)') 'PKtest2:',
!     &              energy(ist),var(ist),
!     &              CPspenrg(ist,ipath),
!     &              mpsi(ist),nparz(1,ist),icharg(ist)
!      enddo
!      write(88,*)



!     ------------------------------------------------------------------

RETURN
END SUBROUTINE static_afterburn


SUBROUTINE return_wavefunctions(ipath)

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'linresp.inc'
!INCLUDE 'dynpar.inc'



INTEGER, INTENT(IN)                  :: ipath


INTEGER :: ist,ir,iz

DO ist=1,nst(1)+nst(2)
  DO ir=0,nr
    DO iz=0,nz
      psi(iz,ir,1,ist) = cppsi(iz,ir,1,ist,ipath)
      psi(iz,ir,2,ist) = cppsi(iz,ir,2,ist,ipath)
    END DO
  END DO
END DO

RETURN
END SUBROUTINE return_wavefunctions

SUBROUTINE linprint()

USE dynpar, ONLY: outnam
!CHARACTER (LEN=3), INTENT(IN)        :: outnam


  OPEN(75,STATUS='unknown',FORM='formatted',&
          POSITION='append',FILE='crank2.'//outnam)
    IF(cpbeta(0).NE.0D0 .AND. b(0,3).NE.0D0) &
    WRITE(75,'(1x,f9.5,1pg20.12,5(1pg12.4),0pf13.5,0pf7.3, 4(1pg12.4))')  &
    cpbeta(0),0.5D0/b(0,3),0.5D0/b(0,2),0D0,0.5D0/b(1,3),2D0*pvar(1,3),pvar(0,3),  &
    cpenergy(0),cpzpemod(0),2D0*bn(0,3),2D0*bn(0,2),zpetri,bn(1,3)!/b(1,3)**2 !/cpbeta(0)**2/1.5D0
  CLOSE(75)

  OPEN(176,STATUS='unknown',FORM='formatted',&
           POSITION='append',FILE='zpeco2.'//outnam)
  IF(cpbeta(0).NE.0D0 .AND. b(0,3).NE.0D0) &
    WRITE(176,'(1x,f9.5,1pg20.12,2(1pg12.4),1pg20.12,1pg12.4,  &
       1pg20.12,2(1pg12.4),0pf10.4,5f9.5)') cpbeta(0),b(0,3)*pvar(0,3), &
       b(0,2)*pvar(0,2),0D0,2D0*b(1,3)*pvar(1,3),  &
       2D0*pvar(1,3),cpenergy(0)-cpzpemod(0),zpetri,omgtri, &
       rchrg2,betawp,rdmsc,surfc,betafull,betafullp
   CLOSE(176)

END SUBROUTINE linprint

!=======================================================================

!     upd_derive_spinor(psia,drpsia,dzpsia,dppsia,ddpsia,zpar,imu,bpar)

!     calculates the derivatives of an array of sp-wavefunctions.

!     imu is an integer global shift in K component which e.g. could be
!     the result of the action of a multipole operator:

!         K(psia)    =  K(psi) + imu

!     zpar is an shift in z-parity

!         Pi_z(psia) =  zpar*Pi_z(psi)

!     where psi are the states of the stationary iteration.

!=======================================================================

SUBROUTINE upd_derive_spinor(psia,drpsia,dzpsia,dppsia,ddpsia,  &
    zpar,imu,nstates)

!     ------------------------------------------------------------------

USE paramr
USE levels
USE grid, ONLY: nablaz_spin,nablar_spin,nablap_spin,laplace_spin
IMPLICIT REAL(8) (a-h,o-z)

!     ------------------------------------------------------------------


REAL(8), INTENT(IN OUT)         :: psia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: drpsia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dzpsia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dppsia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: ddpsia(0:nzd,0:nrd,nsd,npsid)
INTEGER, INTENT(IN)             :: zpar
INTEGER, INTENT(IN)             :: imu
INTEGER, INTENT(IN)             :: nstates
INTEGER :: zparist





INTEGER :: ist

!     calculate derivatives of single particle wavefunctions
!     ------------------------------------------------------------------
DO ist=1,nstates
  
  IF (zpar == +1) THEN
    zparist = nparz(1,ist)
  ELSE
    zparist = 3-nparz(1,ist)
  END IF
  
!         write(88,*) ist, mpsi(ist), icharg(ist), nparz(1,ist), zpar
  
  CALL nablaz_spin(psia(0:nz,0:nr,1:2,ist),dzpsia(0:nz,0:nr,1:2,ist), zparist,1,one)
  CALL nablar_spin(psia(0:nz,0:nr,1:2,ist),drpsia(0:nz,0:nr,1:2,ist), mpsi(ist)+imu,1,one)
  CALL nablap_spin(psia(0:nz,0:nr,1:2,ist),dppsia(0:nz,0:nr,1:2,ist), mpsi(ist)+imu,1,one)
  CALL laplace_spin(psia(0:nz,0:nr,1:2,ist),ddpsia(0:nz,0:nr,1:2,ist),  &
      zparist,mpsi(ist)+imu,1,one)
END DO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE upd_derive_spinor

!=======================================================================

!     upd_derive_spinor_range

!     same for only one spinor at position "ist"

!=======================================================================

SUBROUTINE upd_derive_spinor_range( psia,drpsia,dzpsia,dppsia,ddpsia,  &
    zpar,imu,istart,iend)

!     ------------------------------------------------------------------

USE paramr
USE levels
USE grid, ONLY: nablaz_spin,nablar_spin,nablap_spin,laplace_spin
IMPLICIT REAL(8) (a-h,o-z)

!     ------------------------------------------------------------------


REAL(8), INTENT(IN OUT)         :: psia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: drpsia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dzpsia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dppsia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: ddpsia(0:nzd,0:nrd,nsd,npsid)
INTEGER, INTENT(IN)             :: zpar
INTEGER, INTENT(IN)             :: imu
INTEGER, INTENT(IN)             :: istart
INTEGER, INTENT(IN)             :: iend

INTEGER :: zparist
INTEGER :: ist

!     calculate derivatives of single particle wavefunctions
!     ------------------------------------------------------------------


DO ist=istart,iend
  
  IF (zpar == +1) THEN
    zparist = nparz(1,ist)
  ELSE
    zparist = 3-nparz(1,ist)
  END IF
  
!         write(88,*) ist, mpsi(ist), icharg(ist), nparz(1,ist), zpar
  
  CALL nablaz_spin(psia(0:nz,0:nr,1:2,ist),dzpsia(0:nz,0:nr,1:2,ist), zparist,1,one)
  CALL nablar_spin(psia(0:nz,0:nr,1:2,ist),drpsia(0:nz,0:nr,1:2,ist), mpsi(ist)+imu,1,one)
  CALL nablap_spin(psia(0:nz,0:nr,1:2,ist),dppsia(0:nz,0:nr,1:2,ist), mpsi(ist)+imu,1,one)
  CALL laplace_spin(psia(0:nz,0:nr,1:2,ist),ddpsia(0:nz,0:nr,1:2,ist),  &
      zparist,mpsi(ist)+imu,1,one)
  
END DO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE upd_derive_spinor_range

!=======================================================================

!     damp_spinor(dampr,dampz,psia,psib,zpar,imu,bpar)

!     kinetic energy damping of an array of sp-wavefunction.

!     imu is an integer global shift in K component which e.g. could be
!     the result of the action of a multipole operator:

!         K(psia)    =  K(psi) + imu

!     zpar is an shift in z-parity

!         Pi_z(psia) =  zpar*Pi_z(psi)

!     where psi are the states of the stationary iteration.

!=======================================================================

SUBROUTINE damp_spinor(dampr,dampz,psia,psib,zpar,imu,fac, nstates)

!     ------------------------------------------------------------------

USE paramr
USE levels
USE grid
USE dynpar !??
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'

!     ------------------------------------------------------------------


REAL(8), INTENT(IN OUT)         :: dampr(0:nrd,0:nrd,0:kmpsi)
REAL(8), INTENT(IN OUT)         :: dampz(0:nzd,0:nzd,2)
REAL(8), INTENT(IN OUT)         :: psia(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: psib(0:nzd,0:nrd,nsd,npsid)
INTEGER, INTENT(IN)                      :: zpar
INTEGER, INTENT(IN)                      :: imu
REAL(8), INTENT(IN OUT)         :: fac
INTEGER, INTENT(IN)                      :: nstates

INTEGER :: zparist
INTEGER :: isigma,izpact,imact
INTEGER :: ist,is


!     damp the field of wavefunctions
!     ------------------------------------------------------------------
DO ist=1,nstates
  
  IF (zpar == +1) THEN
    zparist = nparz(1,ist)
  ELSE
    zparist = 3-nparz(1,ist)
  END IF
  
  DO is=1,nsd
    
!           lower spin component has m-component imact+1
    imact  = ABS(mpsi(ist)+imu+is-1)
    
!           lower spin component has opposite z-parity if isym=1
!           if isym=0 izpact should be always +1
    isigma = 3-2*is
    izpact =(3*(1-isigma)/2+zparist*isigma)*isym+(1-isym)
    
    CALL mulr(dampr(0:nr,0:nr,imact), psia(0:nz,0:nr,is,ist),fldbuf(0:nz,0:nr),one,0)
    CALL mulz(dampz(0,0,izpact), fldbuf(0,0),psib(0,0,is,ist),fac,0)
    
  END DO
END DO

!     ------------------------------------------------------------------

RETURN
END SUBROUTINE damp_spinor

!---------------------------------------------------------------------------  
! DESCRIPTION: angmomop
!> @brief
!! Calculates the action of agular momentum operators \f$J_+\f$ and \f$J_-\f$
!>
!> @param[in] wfin
!> REAL(8), array, takes the input wavefunction field.
!> @param[out] amoutp, amoutn
!> REAL(8), arrays, returns action of  \f$J_+\f$ and \f$J_-\f$ on w.f.
!> @param[in] npz
!> INTEGER, takes the z-parity of upper spinor (1 for positive, 2 for 
!!negative, 1 and isym=0 for no good z-parity).
!> @param[in] mpsin
!> INTEGER, takes \f$L_z\f$ angular momentum of upper component of w.f.
!!(zero-derivative).
!> @param[in] factor
!> REAL(8), takes the weight factor.
!--------------------------------------------------------------------------- 

SUBROUTINE angmomop(wfin,amoutp,amoutm,npz,mpsin)

USE paramr, ONLY: nsd
USE grid

IMPLICIT NONE

REAL(8), INTENT(IN)        :: wfin(0:nzd,0:nrd,nsd)
REAL(8), INTENT(OUT)       :: amoutp(0:nzd,0:nrd,nsd)
REAL(8), INTENT(OUT)       :: amoutm(0:nzd,0:nrd,nsd)
INTEGER, INTENT(IN)        :: npz
INTEGER, INTENT(IN)        :: mpsin


INTEGER :: is,ir,iz,mm
REAL(8) :: z,r,rinv

!                                     first, the L+ and L-
  CALL nablaz_spin(wfin(0:nz,0:nr,1:2),psi1(0:nz,0:nr,1:2),npz,1,1D0)
  CALL nablar_spin(wfin(0:nz,0:nr,1:2),psi2(0:nz,0:nr,1:2),mpsin,1,1D0)
  
  DO is=1,nsd
    mm     = mpsin+is-1
!   first the special case r=0
    ir = 0
    IF(mm == 1) THEN
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
        amoutp(iz,ir,is) = 0D0
        amoutm(iz,ir,is) = -(z+z)*psi2(iz,ir,is)
      END DO
    ELSE
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
        amoutp(iz,ir,is) = z*psi2(iz,ir,is)
        amoutm(iz,ir,is) = -amoutp(iz,ir,is)
      END DO
    END IF
    DO ir=1,nr
      r    = (ir+sr)*dr
      rinv = mm/r
      DO iz=0,nz
        z     = (iz+sz-icm)*dz
        amoutp(iz,ir,is) = z*(psi2(iz,ir,is)-rinv *wfin(iz,ir,is))  &
                          -r*psi1(iz,ir,is)
        amoutm(iz,ir,is) =-z*(psi2(iz,ir,is)+rinv *wfin(iz,ir,is))  &
                          +r*psi1(iz,ir,is)
      END DO
    END DO
  END DO
!                                          the S+ and S-
  amoutp(0:nz,0:nr,1) = amoutp(0:nz,0:nr,1)+wfin(0:nz,0:nr,2)
  amoutm(0:nz,0:nr,2) = amoutm(0:nz,0:nr,2)+wfin(0:nz,0:nr,1)

RETURN
END SUBROUTINE angmomop

!=======================================================================

!     act_mfh_even(psi,drpsi,dzpsi,dppsi,Hpsi
!                  fldU,actU,
!                  fldB,drfldB,dzfldB,actB,
!                  fldQr,fldWz,actW,
!                  fldQ,lambda,actQ,
!                  diagE,actE)

!     computes the action of a Mean-Field-Hamiltonian on an array of
!     states Psi(1..nst).

!=======================================================================

SUBROUTINE act_mfh_even( psi,drpsi,dzpsi,dppsi,ddpsi,h_psi,  &
    fldu,                       actu, fldb,drfldb,dzfldb,dpfldb,  actb,  &
    fldwr,fldwz,fldwp,          actw, fldq,lambda,                actq,  &
    diage,                      acte, imu,fac,istart,iend)

!------------------------------------------------------------------------

USE paramr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'  ! nz,nr


REAL(8), INTENT(IN OUT)         :: psi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: drpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dzpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dppsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: ddpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(OUT)            :: h_psi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: fldu(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: actu
REAL(8), INTENT(IN OUT)         :: fldb(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: drfldb(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: dzfldb(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: dpfldb(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: actb
REAL(8), INTENT(IN OUT)         :: fldwr(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldwz(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldwp(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: actw
REAL(8), INTENT(IN OUT)         :: fldq(0:nzd,0:nrd)
REAL(8), INTENT(IN OUT)         :: lambda
INTEGER, INTENT(IN)                      :: actq
REAL(8), INTENT(IN OUT)         :: diage(npsid)
INTEGER, INTENT(IN)                      :: acte
INTEGER, INTENT(IN)                      :: imu
REAL(8), INTENT(IN)             :: fac
INTEGER, INTENT(IN)                      :: istart
INTEGER, INTENT(IN)                      :: iend

INTEGER :: iq,ist,ir,iz

!------------------------------------------------------------------------

!     loop all states and define isospin
!     (first all components which also contribute in axial symmetry)
!     ===================================================================
DO ist=istart,iend
  
!     define isospin
  
  iq=icharg(ist)
!         write(88,*) ist,iq
  
!     if accumulate switch is off --> zero HPsi
  
  DO ir=0,nr
    DO iz=0,nz
      h_psi(iz,ir,1,ist) = zero
      h_psi(iz,ir,2,ist) = zero
    END DO
  END DO
  
!     calculate action of hamiltonian on state psi...
  
  
!     the local potential
  
  IF (actu == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
            + fldu(iz,ir,iq) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
            + fldu(iz,ir,iq) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     the effective mass
  
  IF (actb == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
            - fldb(iz,ir,iq)   * ddpsi(iz,ir,1,ist)  &
            - drfldb(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
            - dzfldb(iz,ir,iq) * dzpsi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
            - fldb(iz,ir,iq)   * ddpsi(iz,ir,2,ist)  &
            - drfldb(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
            - dzfldb(iz,ir,iq) * dzpsi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     the spin-orbit potential
  
  IF (actw == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
            -  fldwz(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
            +  fldwr(iz,ir,iq) * dzpsi(iz,ir,2,ist)  &
            +  fldwr(iz,ir,iq) * dppsi(iz,ir,1,ist)  &
            -  fldwz(iz,ir,iq) * dppsi(iz,ir,2,ist)
        h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
            +  fldwz(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
            -  fldwr(iz,ir,iq) * dzpsi(iz,ir,1,ist)  &
            -  fldwr(iz,ir,iq) * dppsi(iz,ir,2,ist)  &
            -  fldwz(iz,ir,iq) * dppsi(iz,ir,1,ist)
      END DO
    END DO
  END IF
  
!     the deformation constraint
  
  IF (actq == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) =   h_psi(iz,ir,1,ist)  &
            - lambda*fldq(iz,ir) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) =   h_psi(iz,ir,2,ist)  &
            - lambda*fldq(iz,ir) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     the norm constraint (energy, diagonal)
  
  IF (acte == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) =   h_psi(iz,ir,1,ist)  &
            - diage(ist) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) =   h_psi(iz,ir,2,ist)  &
            - diage(ist) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
  IF (acte == -1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) =   h_psi(iz,ir,1,ist)  &
            + diage(ist) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) =   h_psi(iz,ir,2,ist)  &
            + diage(ist) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     end of loop over all states
  
END DO

!     loop all states and define isospin
!     (now all components which contribute in case of NO AXIAL symmetry)
!     ===================================================================
IF (imu /= 0) THEN
  
  DO ist=istart,iend
    
!     define isospin
    
    iq=icharg(ist)
    
!     the effective mass
    
    IF (actb == 1) THEN
      DO ir=0,nr
        DO iz=0,nz
          h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
              - dpfldb(iz,ir,iq) * dppsi(iz,ir,1,ist)
          h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
              - dpfldb(iz,ir,iq) * dppsi(iz,ir,2,ist)
        END DO
      END DO
    END IF
    
!     the spin-orbit potential
    
    IF (actw == 1) THEN
      DO ir=0,nr
        DO iz=0,nz
          h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
              -  fldwp(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
              +  fldwp(iz,ir,iq) * dzpsi(iz,ir,2,ist)
          h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
              +  fldwp(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
              +  fldwp(iz,ir,iq) * dzpsi(iz,ir,1,ist)
          
        END DO
      END DO
    END IF
    
!     end of loop over all states
    
  END DO
  
END IF

!     weight

DO ist=istart,iend
  DO ir=0,nr
    DO iz=0,nz
      h_psi(iz,ir,1,ist) = fac*h_psi(iz,ir,1,ist)
      h_psi(iz,ir,2,ist) = fac*h_psi(iz,ir,2,ist)
    END DO
  END DO
END DO

!------------------------------------------------------------------------

RETURN
END SUBROUTINE act_mfh_even


!=======================================================================

!     act_mfh_even(psi,drpsi,dzpsi,dppsi,Hpsi
!                  fldU,actU,
!                  fldB,drfldB,dzfldB,actB,
!                  fldWr,fldWz,actW,
!                  fldQ,lambda,actQ,
!                  diagE,actE)

!     computes the action of a Mean-Field-Hamiltonian on an array of
!     states Psi(1..nst).

!=======================================================================

SUBROUTINE act_mfh_even2( psi,drpsi,dzpsi,dppsi,ddpsi,h_psi,  &
    fldu,                       actu, fldb,drfldb,dzfldb,dpfldb,  actb,  &
    fldwr,fldwz,fldwp,          actw, fldq,lambda,                actq,  &
    diage,                      acte, imu,fac,nstates)

!------------------------------------------------------------------------

USE paramr
USE levels
USE grid, ONLY: nz,nr
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'  ! nz,nr


REAL(8), INTENT(IN OUT)         :: psi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: drpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dzpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dppsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: ddpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(OUT)            :: h_psi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: fldu(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: actu
REAL(8), INTENT(IN OUT)         :: fldb(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: drfldb(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: dzfldb(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: dpfldb(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: actb
REAL(8), INTENT(IN OUT)         :: fldwr(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldwz(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldwp(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: actw
REAL(8), INTENT(IN OUT)         :: fldq(0:nzd,0:nrd)
REAL(8), INTENT(IN OUT)         :: lambda
INTEGER, INTENT(IN)                      :: actq
REAL(8), INTENT(IN OUT)         :: diage(npsid)
INTEGER, INTENT(IN)                      :: acte
INTEGER, INTENT(IN)                      :: imu
REAL(8), INTENT(IN)             :: fac
INTEGER, INTENT(IN)                      :: nstates

INTEGER :: iq,ist,ir,iz

!------------------------------------------------------------------------

!     loop all states and define isospin
!     (first all components which also contribute in axial symmetry)
!     ===================================================================
DO ist=1,nstates
  
!     define isospin
  
  iq=icharg(ist)
!         write(88,*) ist,iq
  
!     if accumulate switch is off --> zero HPsi
  
  DO ir=0,nr
    DO iz=0,nz
      h_psi(iz,ir,1,ist) = zero
      h_psi(iz,ir,2,ist) = zero
    END DO
  END DO
  
!     calculate action of hamiltonian on state psi...
  
  
!     the local potential
  
  IF (actu == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
            + fldu(iz,ir,iq) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
            + fldu(iz,ir,iq) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     the effective mass
  
  IF (actb == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
            - fldb(iz,ir,iq)   * ddpsi(iz,ir,1,ist)  &
            - drfldb(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
            - dzfldb(iz,ir,iq) * dzpsi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
            - fldb(iz,ir,iq)   * ddpsi(iz,ir,2,ist)  &
            - drfldb(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
            - dzfldb(iz,ir,iq) * dzpsi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     the spin-orbit potential
  
  IF (actw == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
            -  fldwz(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
            +  fldwr(iz,ir,iq) * dzpsi(iz,ir,2,ist)  &
            +  fldwr(iz,ir,iq) * dppsi(iz,ir,1,ist)  &
            -  fldwz(iz,ir,iq) * dppsi(iz,ir,2,ist)
        h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
            +  fldwz(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
            -  fldwr(iz,ir,iq) * dzpsi(iz,ir,1,ist)  &
            -  fldwr(iz,ir,iq) * dppsi(iz,ir,2,ist)  &
            -  fldwz(iz,ir,iq) * dppsi(iz,ir,1,ist)
      END DO
    END DO
  END IF
  
!     the deformation constraint
  
  IF (actq == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) =   h_psi(iz,ir,1,ist)  &
            - lambda*fldq(iz,ir) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) =   h_psi(iz,ir,2,ist)  &
            - lambda*fldq(iz,ir) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     the norm constraint (energy, diagonal)
  
  IF (acte == 1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) =   h_psi(iz,ir,1,ist)  &
            - diage(ist) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) =   h_psi(iz,ir,2,ist)  &
            - diage(ist) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
  IF (acte == -1) THEN
    DO ir=0,nr
      DO iz=0,nz
        h_psi(iz,ir,1,ist) =   h_psi(iz,ir,1,ist)  &
            + diage(ist) * psi(iz,ir,1,ist)
        h_psi(iz,ir,2,ist) =   h_psi(iz,ir,2,ist)  &
            + diage(ist) * psi(iz,ir,2,ist)
      END DO
    END DO
  END IF
  
!     end of loop over all states
  
END DO

!     loop all states and define isospin
!     (now all components which contribute in case of NO AXIAL symmetry)
!     ===================================================================
IF (imu /= 0) THEN
  
  DO ist=1,nstates
    
!     define isospin
    
    iq=icharg(ist)
    
!     the effective mass
    
    IF (actb == 1) THEN
      DO ir=0,nr
        DO iz=0,nz
          h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
              - dpfldb(iz,ir,iq) * dppsi(iz,ir,1,ist)
          h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
              - dpfldb(iz,ir,iq) * dppsi(iz,ir,2,ist)
        END DO
      END DO
    END IF
    
!     the spin-orbit potential
    
    IF (actw == 1) THEN
      DO ir=0,nr
        DO iz=0,nz
          h_psi(iz,ir,1,ist) = h_psi(iz,ir,1,ist)  &
              -  fldwp(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
              +  fldwp(iz,ir,iq) * dzpsi(iz,ir,2,ist)
          h_psi(iz,ir,2,ist) = h_psi(iz,ir,2,ist)  &
              +  fldwp(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
              +  fldwp(iz,ir,iq) * dzpsi(iz,ir,1,ist)
          
        END DO
      END DO
    END IF
    
!     end of loop over all states
    
  END DO
  
END IF

!     weight

DO ist=1,nstates
  DO ir=0,nr
    DO iz=0,nz
      h_psi(iz,ir,1,ist) = fac*h_psi(iz,ir,1,ist)
      h_psi(iz,ir,2,ist) = fac*h_psi(iz,ir,2,ist)
    END DO
  END DO
END DO

!------------------------------------------------------------------------

RETURN
END SUBROUTINE act_mfh_even2


!=======================================================================

!     act_mfh_odd(psi,drpsi,dzpsi,dppsi,Hpsi
!                  fldU,actU,
!                  fldB,drfldB,dzfldB,actB,
!                  fldQr,fldWz,actW,
!                  fldQ,lambda,actQ,
!                  diagE,actE)

!     computes the action of a Mean-Field-Hamiltonian on an array of
!     states Psi(1..nst).

!=======================================================================

SUBROUTINE act_mfh_odd( psi,drpsi,dzpsi,dppsi,ddpsi,hpsi,  &
    fldar,fldaz,fldap,flddiva,   acta, fldsr,fldsz,fldsp,           acts,  &
    fldcr,fldcz,fldcp,           actc, imu,fac)

!------------------------------------------------------------------------

USE paramr
USE levels
USE grid, ONLY: nz,nr
IMPLICIT REAL(8) (a-h,o-z)
!INCLUDE 'dynpar.inc'  ! nz,nr


REAL(8), INTENT(IN OUT)         :: psi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: drpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dzpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: dppsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: ddpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(OUT)            :: hpsi(0:nzd,0:nrd,nsd,npsid)
REAL(8), INTENT(IN OUT)         :: fldar(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldaz(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldap(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)                     :: flddiva
INTEGER, INTENT(IN)                      :: acta
REAL(8), INTENT(IN OUT)         :: fldsr(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldsz(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldsp(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN)                      :: acts
REAL(8), INTENT(IN OUT)         :: fldcr(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldcz(0:nzd,0:nrd,kqm)
REAL(8), INTENT(IN OUT)         :: fldcp(0:nzd,0:nrd,kqm)
INTEGER, INTENT(IN OUT)                  :: actc
INTEGER, INTENT(IN)                      :: imu
REAL(8), INTENT(IN)             :: fac


REAL(8) :: divflda(0:nzd,0:nrd,kqm)
INTEGER :: iq,ist,ir,iz



!------------------------------------------------------------------------

!     loop all states and define isospin
!     (first all components which also contribute in axial symmetry)
!     ===================================================================
DO ist=1,nst(1)+nst(2)
  
!     define isospin
  
  IF (ist <= nst(1)) THEN
    iq=1
  ELSE
    iq=2
  END IF
  
!     if accumulate switch is off --> zero HPsi
  
  DO ir=0,nr
    DO iz=0,nz
      hpsi(iz,ir,1,ist) = zero
      hpsi(iz,ir,2,ist) = zero
    END DO
  END DO
  
!     calculate action of time-odd hamiltonian on state psi...
  
  DO ir=0,nr
    
!     the dynamic potential
    
    IF (acta == 1) THEN
      DO iz=0,nz
        hpsi(iz,ir,1,ist) = hpsi(iz,ir,1,ist)  &
            + 0.5D0*divflda(iz,ir,iq) *   psi(iz,ir,1,ist)  &
            +         fldar(iz,ir,iq) * drpsi(iz,ir,1,ist)  &
            +         fldaz(iz,ir,iq) * dzpsi(iz,ir,1,ist)
      END DO
      DO iz=0,nz
        hpsi(iz,ir,2,ist) = hpsi(iz,ir,2,ist)  &
            + 0.5D0*divflda(iz,ir,iq) *   psi(iz,ir,2,ist)  &
            +         fldar(iz,ir,iq) * drpsi(iz,ir,2,ist)  &
            +         fldaz(iz,ir,iq) * dzpsi(iz,ir,2,ist)
      END DO
    END IF
    
!     the spin-potential
    
    IF (acts == 1) THEN
      DO iz=0,nz
        hpsi(iz,ir,1,ist) = hpsi(iz,ir,1,ist)  &
            - fldsp(iz,ir,iq)  * psi(iz,ir,2,ist)
        hpsi(iz,ir,2,ist) = hpsi(iz,ir,2,ist)  &
            + fldsp(iz,ir,iq)  * psi(iz,ir,1,ist)
      END DO
    END IF
  END DO
  
!     the effective spin-mass
  
!         if (actC.eq.1) then
!            do ir=0,nr
!               do iz=0,nz
!                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist)
!     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist)
!     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist)
!     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist)
!     &                          +  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist)
!                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist)
!     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
!     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
!     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
!     &                          +  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist)
!               enddo
!            enddo
!         endif
  
!     end of loop over all states
  
END DO

!     loop all states and define isospin
!     (now all components which contribute in case of NO AXIAL symmetry)
!     ===================================================================
IF (imu > 0) THEN
  
  DO ist=1,nst(1)+nst(2)
    
!     define isospin
    
    IF (ist <= nst(1)) THEN
      iq=1
    ELSE
      iq=2
    END IF
    
    DO ir=0,nr
      
!     the dynamic potential
      
      IF (acta == 1) THEN
        DO iz=0,nz
          hpsi(iz,ir,1,ist) = hpsi(iz,ir,1,ist)  &
              +         fldap(iz,ir,iq) * dppsi(iz,ir,1,ist)
          hpsi(iz,ir,2,ist) = hpsi(iz,ir,2,ist)  &
              +         fldap(iz,ir,iq) * dppsi(iz,ir,2,ist)
        END DO
      END IF
      
!     the spin-potential
      
      IF (acts == 1) THEN
        DO iz=0,nz
          hpsi(iz,ir,1,ist) = hpsi(iz,ir,1,ist)  &
              + fldsr(iz,ir,iq)  * psi(iz,ir,2,ist)  &
              + fldsz(iz,ir,iq)  * psi(iz,ir,1,ist)
          hpsi(iz,ir,2,ist) = hpsi(iz,ir,2,ist)  &
              + fldsr(iz,ir,iq)  * psi(iz,ir,1,ist)  &
              - fldsz(iz,ir,iq)  * psi(iz,ir,2,ist)
        END DO
      END IF
    END DO
    
    
!     the effective spin-mass
    
!         if (actC.eq.1) then
!            do ir=0,nr
!               do iz=0,nz
!                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist)
!     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist)
!     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist)
!     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist)
!     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist)
!                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist)
!     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
!     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
!     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
!     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist)
!               enddo
!            enddo
!         endif
    
!     end of loop over all states
    
  END DO
  
ELSE IF (imu < 0) THEN
  
  DO ist=1,nst(1)+nst(2)
    
!     define isospin
    
    IF (ist <= nst(1)) THEN
      iq=1
    ELSE
      iq=2
    END IF
    DO ir=0,nr
      
!     the dynamic potential
      
      IF (acta == 1) THEN
        DO iz=0,nz
          hpsi(iz,ir,1,ist) = hpsi(iz,ir,1,ist)  &
              -         fldap(iz,ir,iq) * dppsi(iz,ir,1,ist)
          hpsi(iz,ir,2,ist) = hpsi(iz,ir,2,ist)  &
              -         fldap(iz,ir,iq) * dppsi(iz,ir,2,ist)
        END DO
      END IF
      
!     the spin-potential
      
      IF (acts == 1) THEN
        DO iz=0,nz
          hpsi(iz,ir,1,ist) = hpsi(iz,ir,1,ist)  &
              - fldsr(iz,ir,iq)  * psi(iz,ir,2,ist)  &
              - fldsz(iz,ir,iq)  * psi(iz,ir,1,ist)
          hpsi(iz,ir,2,ist) = hpsi(iz,ir,2,ist)  &
              - fldsr(iz,ir,iq)  * psi(iz,ir,1,ist)  &
              + fldsz(iz,ir,iq)  * psi(iz,ir,2,ist)
        END DO
      END IF
    END DO
    
!     the effective spin-mass
    
!         if (actC.eq.1) then
!            do ir=0,nr
!               do iz=0,nz
!                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist)
!     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist)
!     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist)
!     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist)
!     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist)
!                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist)
!     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
!     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
!     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
!     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist)
!               enddo
!            enddo
!         endif
    
!     end of loop over all states
    
  END DO
  
END IF

!     weight

DO ist=1,nst(1)+nst(2)
  DO ir=0,nr
    DO iz=0,nz
      hpsi(iz,ir,1,ist) = fac*hpsi(iz,ir,1,ist)
      hpsi(iz,ir,2,ist) = fac*hpsi(iz,ir,2,ist)
    END DO
  END DO
END DO

!------------------------------------------------------------------------

RETURN
END SUBROUTINE act_mfh_odd

SUBROUTINE hhpsi(psiin,psiout,ist)

!**********************************************************************
!                                                                     *
!     same as subroutine hpsi -- but evaluates the variances by       *
!     <ist|HHist> - <ist|Hist>**2                                     *
!                                                                     *
!**********************************************************************

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT REAL(8) (a-h,o-z)

REAL(8), INTENT(IN)                         :: psiin(0:nzd,0:nrd,nsd)
REAL(8), INTENT(OUT)                        :: psiout(0:nzd,0:nrd,nsd)
INTEGER, INTENT(IN OUT)                     :: ist


INTEGER,PARAMETER :: ifclean=1

!**********************************************************************

!CALL restore_boundary(psiin,0,+1,ist) !BS: Routine does not do anything anyway!

iq      = icharg(ist)
ebin    = spenrg(ist)

CALL laplace_spin(psiin(0:nz,0:nr,1:2),psi1(0:nz,0:nr,1:2),nparz(1,ist),mpsi(ist),1,one)
CALL nablar_spin(psiin(0:nz,0:nr,1:2),psi2(0:nz,0:nr,1:2),   mpsi(ist),1,one)
CALL nablaz_spin(psiin(0:nz,0:nr,1:2),psi3(0:nz,0:nr,1:2),nparz(1,ist),1,one)

mmup      = mpsi(ist)
mmdown    = mpsi(ist)+1
nparrup   = 1+MOD(mmup,2)
nparrdown = 1+MOD(mmdown,2)
sigis     = MOD(one*(is),two)*two - one

IF(ifstab == 1) THEN
  eferm2 = fer(iq)*oh
  IF(ebin > zero) THEN
    spefac = one
  ELSE IF(ebin > eferm2) THEN
    spefac = one-ebin/eferm2
  ELSE
    spefac = zero
  END IF
ELSE
  spefac = zero
END IF

IF(spefac <= zero) THEN
  DO ir=0,nr
    DO iz=0,nz
      psiout(iz,ir,1) =  &
          (pot(iz,ir,iq)-confield(iz,ir)-ebin)*psiin(iz,ir,1)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,1)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
      psiout(iz,ir,2) =  &
          (pot(iz,ir,iq)-confield(iz,ir)-ebin)*psiin(iz,ir,2)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,2)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
    END DO
  END DO
ELSE
  DO ir=0,nr
    DO iz=0,nz
      psiout(iz,ir,1) = (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)  &
          - confield(iz,ir)-ebin)*psiin(iz,ir,1)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,1)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
      psiout(iz,ir,2) = (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)  &
          - confield(iz,ir)-ebin)*psiin(iz,ir,2)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,2)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
    END DO
  END DO
END IF

IF(t4 /= zero) THEN
  
!     the spin-orbit part. action of the s.o. force is computed
!     simply from
!       W_ls |psi> ~ (nabla rho)*sigma sigma*nabla |psi>
!                    - (nabla rho) * nabla |psi>
  
!                               r=0 requires l'Hospital for 1/r*psiin
  IF(mmup == 1) THEN
    DO iz=0,nz
      psiout(iz,0,1) = +              potls(iz,0,iq,1)*psi2(iz,0,1)  &
          +              potls(iz,0,iq,1)*psi3(iz,0,2)  &
          -              potls(iz,0,iq,2)*psi2(iz,0,2) + psiout(iz,0,1)
      
      psiout(iz,0,2) = -              potls(iz,0,iq,1)*psi3(iz,0,1)  &
          +              potls(iz,0,iq,2)*psi2(iz,0,1)  &
          -              potls(iz,0,iq,2)*psi2(iz,0,1) + psiout(iz,0,2)
    END DO
  ELSE IF(mmdown == 1) THEN
    psiout(iz,0,1) = +              potls(iz,0,iq,1)*psi3(iz,0,2)  &
        -              potls(iz,0,iq,2)*psi2(iz,0,2)  &
        -              potls(iz,0,iq,2)*psi2(iz,0,2) + psiout(iz,0,1)
    
    psiout(iz,0,2) = -              potls(iz,0,iq,1)*psi2(iz,0,2)  &
        -              potls(iz,0,iq,1)*psi3(iz,0,1)  &
        +              potls(iz,0,iq,2)*psi2(iz,0,1) + psiout(iz,0,2)
  END IF
  
!                               now all other r
  DO ir=1,nr
    
    r       = (sr+ir)*dr
    rinv    = one/r
    xmmup   = mmup  *rinv
    xmmdown = mmdown*rinv
    
    DO iz=0,nz
      psiout(iz,ir,1) = +        xmmup*potls(iz,ir,iq,1)*psiin(iz,ir,1)  &
          +              potls(iz,ir,iq,1)*psi3(iz,ir,2)  &
          -              potls(iz,ir,iq,2)*psi2(iz,ir,2)  &
          -      xmmdown*potls(iz,ir,iq,2)*psiin(iz,ir,2) + psiout(iz,ir,1)
      
      psiout(iz,ir,2) = -      xmmdown*potls(iz,ir,iq,1)*psiin(iz,ir,2)  &
          -              potls(iz,ir,iq,1)*psi3(iz,ir,1)  &
          +              potls(iz,ir,iq,2)*psi2(iz,ir,1)  &
          -        xmmup*potls(iz,ir,iq,2)*psiin(iz,ir,1) + psiout(iz,ir,2)
      
    END DO
  END DO
END IF

!     compute expectation value

esub        = wfovlp(psiin(0:nz,0:nr,1:2),psiout(0:nz,0:nr,1:2))
spenrg(ist) = esub+ebin
spenrg2     = wfovlp(psiout(0:nz,0:nr,1:2),psiout(0:nz,0:nr,1:2))
spevar(ist) = SQRT(ABS(spenrg2-esub*esub))

!     clean up the wavefunction from diagonal components
!     ?? note: this may turn out to be unnecessary. Check!

IF(ifclean /= 1) RETURN
DO is=1,2
  DO ir=0,nr
    DO iz=0,nz
      psiout(iz,ir,is) = psiout(iz,ir,is) -esub*psiin(iz,ir,is)
    END DO
  END DO
END DO
RETURN

!     the action of H^2 on psi

CALL laplace_spin(psiout(0:nz,0:nr,1:2),psi1(0:nz,0:nr,1:2),nparz(1,ist),mpsi(ist),1,one)
CALL nablar_spin(psiout(0:nz,0:nr,1:2),psi2(0:nz,0:nr,1:2),   mpsi(ist),1,one)
CALL nablaz_spin(psiout(0:nz,0:nr,1:2),psi3(0:nz,0:nr,1:2),nparz(1,ist),1,one)

mmup      = mpsi(ist)
mmdown    = mpsi(ist)+1
nparrup   = 1+MOD(mmup,2)
nparrdown = 1+MOD(mmdown,2)
sigis     = MOD(one*(is),two)*two - one

IF(ifstab == 1) THEN
  eferm2 = fer(iq)*oh
  IF(ebin > zero) THEN
    spefac = one
  ELSE IF(ebin > eferm2) THEN
    spefac = one-ebin/eferm2
  ELSE
    spefac = zero
  END IF
ELSE
  spefac = zero
END IF

IF(spefac <= zero) THEN
  DO ir=0,nr
    DO iz=0,nz
      psi4(iz,ir,1) =  &
          (pot(iz,ir,iq)-confield(iz,ir)-ebin)*psiout(iz,ir,1)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,1)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
      psi4(iz,ir,2) =  &
          (pot(iz,ir,iq)-confield(iz,ir)-ebin)*psiout(iz,ir,2)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,2)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
    END DO
  END DO
ELSE
  DO ir=0,nr
    DO iz=0,nz
      psi4(iz,ir,1) = (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)  &
          -confield(iz,ir)-ebin)*psiout(iz,ir,1)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,1)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
      psi4(iz,ir,2) = (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)  &
          -confield(iz,ir)-ebin)*psiout(iz,ir,2)  &
          - bmass(iz,ir,iq)                      *psi1(iz,ir,2)  &
          - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)  &
          - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
    END DO
  END DO
END IF

IF(t4 /= zero) THEN
  
!     the spin-orbit part. action of the s.o. force is computed
!     simply from
!       W_ls |psi> ~ (nabla rho)*sigma sigma*nabla |psi>
!                    - (nabla rho) * nabla |psi>
  
!                               r=0 requires l'Hospital for 1/r*psiout
  IF(mmup == 1) THEN
    DO iz=0,nz
      psi4(iz,0,1) = +              potls(iz,0,iq,1)*psi2(iz,0,1)  &
          +              potls(iz,0,iq,1)*psi3(iz,0,2)  &
          -              potls(iz,0,iq,2)*psi2(iz,0,2) + psi4(iz,0,1)
      
      psi4(iz,0,2) = -              potls(iz,0,iq,1)*psi3(iz,0,1)  &
          +              potls(iz,0,iq,2)*psi2(iz,0,1)  &
          -              potls(iz,0,iq,2)*psi2(iz,0,1) + psi4(iz,0,2)
    END DO
  ELSE IF(mmdown == 1) THEN
    psi4(iz,0,1) = +              potls(iz,0,iq,1)*psi3(iz,0,2)  &
        -              potls(iz,0,iq,2)*psi2(iz,0,2)  &
        -              potls(iz,0,iq,2)*psi2(iz,0,2) + psi4(iz,0,1)
    
    psi4(iz,0,2) = -              potls(iz,0,iq,1)*psi2(iz,0,2)  &
        -              potls(iz,0,iq,1)*psi3(iz,0,1)  &
        +              potls(iz,0,iq,2)*psi2(iz,0,1) + psi4(iz,0,2)
  END IF
  
!                               now all other r
  DO ir=1,nr
    
    r       = (sr+ir)*dr
    rinv    = one/r
    xmmup   = mmup  *rinv
    xmmdown = mmdown*rinv
    
    DO iz=0,nz
      psi4(iz,ir,1) = +        xmmup*potls(iz,ir,iq,1)*psiout(iz,ir,1)  &
          +              potls(iz,ir,iq,1)*psi3(iz,ir,2)  &
          -              potls(iz,ir,iq,2)*psi2(iz,ir,2)  &
          -      xmmdown*potls(iz,ir,iq,2)*psiout(iz,ir,2) + psi4(iz,ir,1)
      
      psi4(iz,ir,2) = -      xmmdown*potls(iz,ir,iq,1)*psiout(iz,ir,2)  &
          -              potls(iz,ir,iq,1)*psi3(iz,ir,1)  &
          +              potls(iz,ir,iq,2)*psi2(iz,ir,1)  &
          -        xmmup*potls(iz,ir,iq,2)*psiout(iz,ir,1) + psi4(iz,ir,2)
      
    END DO
  END DO
END IF

51   CONTINUE

!     compute variance

spenrg2     = wfovlp(psiin(0:nz,0:nr,1:2),psi4(0:nz,0:nr,1:2))
spevar(ist) = SQRT(ABS(spenrg2-esub*esub))

RETURN
END SUBROUTINE hhpsi


END MODULE linresponse
