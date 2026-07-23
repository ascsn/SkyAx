c***********************************************************************
c
c  file        : linresp.inc
c  description : global variables of the linear-response module
c
c
c  static parameters
c
c     MZD,MRD              grid dimensions in r,z
c     MID,MSD              dimension of (iso)spin hilbert
c     MPSID                maximum number of sp-wavefunctions
c     MTD                  maximum number of transitions
c
c
c  common block "LRdens" -- response densities
c
c     rerho                response density           (rho)
c     rekin                response kinetic density   (tau)
c     recur[r,p,z]         response current           (j)
c     respn[r,p,z]         response spin density      (sigma)
c     resoc[r,p,z]         response spin-orbit dens.  (J)
c     reskd[r,p,z]         response kin. spin dens.   (tau) (vector)
c     repair               response pairing-density.  (chi)
c  
c
c  common block "LRflds" -- response fields
c
c     relocpot             response local potential   (U)
c     recoupot             response coulomb potential (U_C)
c     reeffmas             response effective mass    (B)
c     redynpot[r,p,z]      response dynamic potential (A)
c     respnpot[r,p,z]      response spin potential    (S)
c     resocpot[r,p,z]      response sp.-orb. pot.     (W)
c     respnmas[r,p,z]      response spin mass         (C)
c
c
c  common block "LRtrans" -- transition arrays
c
c     TRidx                transition table 
c     TRlen                length of transition table
c     TRxi                 occ. factor xi for transition
c     TReta                occ. factor eta for transition
c
c  common block "LRDamp"  --  damping matrices
c
c     E0                   damping parameter
c     LRDampr              damping matrix in r direction
c     LRDampz              damping matrix in z direction
c
c  common block "LRpath" -- storage of the collective path
c
c     CPrho                density
c     CPkin                kinetic density
c     CPsoc[r,z]           spin-orbit-current
c
c     CPlocpot             stationary local potential      (U)
c     CPcoupot             stationary coulomb potential    (U_c)
c     CPeffmas             stationary effective mass       (B)
c     CPsocpot[r,z]        stationary spin-orbit potential (W)
c
c     CPconstr             constraint field                (Q)
c
c     CPpsi                sp-wavefunctions of path
c     CPvsq                occ-amplitudes v*v
c     CPspeng              sp-energies of sp-states
c     CPlen                number of sp-states
c     
c     CPdeltaf             ????
c     CPwtpair             ????
c     CPpairfc             ????
c
c     CPenergy             total energy
c     CPlambda             constraint energy (lagrange-multiplier)
c     CPfer                fermi energies of path elements
c     CPbeta               total deformation of path elements
c     CPrms                total rms-radius of path elements
c     CPzcm                center of mass
c     CPatnum              total number of nucleons of path elements
c
c     CPdrpsi              derivatives of sp-wavefunctions
c     CPdppsi
c     CPdzpsi        
c
c
c  common block "LRbuff" -- buffers for linear-response module
c
c     LRwrk_scalar         buffer for scalar function
c     LRwrk1_spinor        buffers for spinor functions
c     LRwrk2_spinor        
c
c  setting switching parameters:
c     FULLSPACE 
c     D_LRDEBUG 
c     D_OLDCRANK
c
c***********************************************************************

c=======================================================================
c
c     include static parameters
c
c=======================================================================

      integer MZD               ! maximum number of gridpoints in z
      integer MRD               ! maximum number of gridpoints in r
      integer MID               ! dimension of isospin hilbert
      integer MSD               ! dimension of spin hilbert
      integer MPSID             ! number of sp-wavefunctions
      integer MTD               ! maximum number of transitions per state
      integer MTDTOT            ! maximum number of transitions
      integer LAMBDAMAX         ! highest collective Mode in cranking
      parameter(MZD=nzd,   MRD   =nrd, MID=2, MSD=nsd, MPSID=npsid,
     &          MTD=1000,  MTDTOT=100000,  LAMBDAMAX=3)


c=======================================================================
c
c     common block "damping matrices" 
c
c=======================================================================
      common /LRDamp/ LRDampR,LRDampZ,LRe0inv,LRx0,LRstep

      double precision LRe0inv,LRx0,LRstep
      double precision LRDampR(0:MRD,0:MRD,0:KMPSI)
      double precision LRDampZ(0:MZD,0:MZD,2)
      


c=======================================================================
c
c     common block "response densities" 
c
c=======================================================================
      common /LRdens/
     &     rerho, rekin,               
     &     recurr, recurp, recurz,     
     &     respnr, respnp, respnz,     
     &     resocr, resocp, resocz,     
     &     reskdr, reskdp, reskdz, repair  

      double precision rerho(0:MZD,0:MRD,MID)  ! response rho
      double precision rekin(0:MZD,0:MRD,MID)  ! response tau
      double precision recurr(0:MZD,0:MRD,MID) ! response j_r
      double precision recurp(0:MZD,0:MRD,MID) ! response j_phi
      double precision recurz(0:MZD,0:MRD,MID) ! response j_z
      double precision respnr(0:MZD,0:MRD,MID) ! response sigma_r
      double precision respnp(0:MZD,0:MRD,MID) ! response sigma_phi
      double precision respnz(0:MZD,0:MRD,MID) ! response sigma_z
      double precision resocr(0:MZD,0:MRD,MID) ! response J_r
      double precision resocp(0:MZD,0:MRD,MID) ! response J_p
      double precision resocz(0:MZD,0:MRD,MID) ! response J_z
      double precision reskdr(0:MZD,0:MRD,MID) ! response tau_r
      double precision reskdp(0:MZD,0:MRD,MID) ! response tau_p
      double precision reskdz(0:MZD,0:MRD,MID) ! response tau_z
      double precision repair(0:MZD,0:MRD,MID) ! response pair-dens.


c=======================================================================
c
c     common block "response fields" 
c
c=======================================================================
      common /LRflds/
     &     relocpot, recoupot, reeffmas, 
     &     redreffmas,redpeffmas,redzeffmas,  
     &     redynpotr, redynpotp, redynpotz,    
     &     respnpotr, respnpotp, respnpotz,   
     &     resocpotr, resocpotp, resocpotz,   
     &     respnmasr, respnmasp, respnmasz,
     &     locpotwgt1,locpotwgt2,repairpot

      double precision relocpot(0:MZD,0:MRD,MID)  ! response U
      double precision recoupot(0:MZD,0:MRD)      ! response U_C
      double precision reeffmas(0:MZD,0:MRD,MID)  ! response B
      double precision redzeffmas(0:MZD,0:MRD,MID)  ! response nabla_r B
      double precision redpeffmas(0:MZD,0:MRD,MID)  ! response nabla_p B
      double precision redreffmas(0:MZD,0:MRD,MID)  ! response nabla_z B
      double precision redynpotr(0:MZD,0:MRD,MID) ! response A_r
      double precision redynpotp(0:MZD,0:MRD,MID) ! response A_phi
      double precision redynpotz(0:MZD,0:MRD,MID) ! response A_z
      double precision respnpotr(0:MZD,0:MRD,MID) ! response S_r
      double precision respnpotp(0:MZD,0:MRD,MID) ! response S_phi
      double precision respnpotz(0:MZD,0:MRD,MID) ! response S_z
      double precision resocpotr(0:MZD,0:MRD,MID) ! response W_r
      double precision resocpotp(0:MZD,0:MRD,MID) ! response W_phi
      double precision resocpotz(0:MZD,0:MRD,MID) ! response W_z
      double precision respnmasr(0:MZD,0:MRD,MID) ! response C_r
      double precision respnmasp(0:MZD,0:MRD,MID) ! response C_phi
      double precision respnmasz(0:MZD,0:MRD,MID) ! response C_z
      double precision repairpot(0:MZD,0:MRD,MID) ! response pair-pot.
      
      double precision locpotwgt1(0:MZD,0:MRD,MID)! response weight dU_p/drho_q
      double precision locpotwgt2(0:MZD,0:MRD,MID)! response weight dU_n/drho_q

c=======================================================================
c
c     common block "transition table" 
c
c=======================================================================
      common /LRtrans/
     &     TRidx,
     &     TRxi, TReta, TRntr,TRmmax,TRntrtot,
     &     TRlen,nsttot,origin,partner,nstadd,
     &     isresact,isMFresp

      integer           TRntrtot                ! number of transitions
      integer           TRntr(1:MPSID)          ! dto. per state
      integer           TRidx(1:MPSID,1:MTD,8)  ! transition table
      double precision TRxi(1:MPSID,1:MTD,2)   ! occ. factor xi
      double precision TReta(1:MPSID,1:MTD,2)  ! occ. factor eta
      integer          origin(MPSID)      ! original WF in K>0 space
      integer          partner(MPSID)     ! Time-Rev. partner
      integer          TRmmax             ! maximum m-component
      integer          nsttot             ! total size of ext. sp-base
      integer          nstadd(2)          ! number of added time-reversed
      integer          isresact,isMFresp  !                        states


c=======================================================================
c
c     common block "collective path" 
c
c=======================================================================
      common /LRpath/
     &     CPrho, CPkin, CPsocr, CPsocz,
     &     CPlocpot,CPeffmas,CPsocpotr,CPsocpotz,CPcoupot,
     &     CPdreffmas,CPdzeffmas,CPconstr,
     &     CPpsi,CPdrpsi,CPdppsi,CPdzpsi,CPddpsi,
     &     P2qATDHF,Q2qATDHF,P2qCHF,Q2qCHF,
     &     CPdeltaf,CPwtpair,CPspenrg,CPqpenrg,CPfer,CPvsqf,
     &     CPenergy,CPlambda,CPbeta,CPzcm,
     &     CPatnum,CPradRMS,CPradDIF,CPradSRF,
     &     CPenergyC,CPlen
      
      double precision CPrho(0:MZD,0:MRD,MID,-1:1)       ! density 
      double precision CPkin(0:MZD,0:MRD,MID,-1:1)       ! kinetic dens.
      double precision CPsocr(0:MZD,0:MRD,MID,-1:1)      ! sp.-orb. current
      double precision CPsocz(0:MZD,0:MRD,MID,-1:1)      ! 
      double precision CPlocpot(0:MZD,0:MRD,MID,-1:1)    ! local potential
      double precision CPeffmas(0:MZD,0:MRD,MID,-1:1)    ! effective mass
      double precision CPsocpotr(0:MZD,0:MRD,MID,-1:1)   ! sp.-orb. potential
      double precision CPsocpotz(0:MZD,0:MRD,MID,-1:1)   ! 
      double precision CPcoupot(0:MZD,0:MRD,-1:1)        ! coulomb potential
      double precision CPdreffmas(0:MZD,0:MRD,MID)       ! r-deriv. of mass
      double precision CPdzeffmas(0:MZD,0:MRD,MID)       ! z-deriv. of mass

      double precision CPconstr(0:MZD,0:MRD,-1:1)        ! constraint field

      double precision CPpsi(0:MZD,0:MRD,MSD,MPSID,-1:1) ! sp-waves
      double precision CPdrpsi(0:MZD,0:MRD,MSD,MPSID)    ! r-derivatives
      double precision CPdppsi(0:MZD,0:MRD,MSD,MPSID)    ! z-derivatives
      double precision CPdzpsi(0:MZD,0:MRD,MSD,MPSID)    ! phi-derivatives
      double precision CPddpsi(0:MZD,0:MRD,MSD,MPSID)    ! laplace

      double precision P2qATDHF(MTDTOT,0:LAMBDAMAX)
      double precision Q2qATDHF(MTDTOT,0:LAMBDAMAX)
      double precision P2qCHF(MTDTOT,0:LAMBDAMAX)
      double precision Q2qCHF(MTDTOT,0:LAMBDAMAX)

      double precision CPbeta(-1:1)           ! deformations
      double precision CPzcm(-1:1)            ! center of mass
      double precision CPatnum(-1:1)          ! numbr. of nucl.
      double precision CPradRMS(-1:1)         ! rms-radius 
      double precision CPradDIF(-1:1)         ! diff.-radius 
      double precision CPradSRF(-1:1)         ! sigma  
      double precision CPenergy(-1:1)         ! total energy
      double precision CPenergyC(-1:1)        ! total constraint energy
      double precision CPlambda(-1:1)         ! constraint energy
 
      double precision CPdeltaf(MPSID,-1:1) ! ????          
      double precision CPwtpair(MPSID,-1:1) ! ????
      double precision CPpairfc(MPSID,-1:1) ! ????
      double precision CPspenrg(MPSID,-1:1) ! sp-energies
      double precision CPqpenrg(MPSID,-1:1) ! qp-energies
      double precision CPvsqf(MPSID,-1:1)   ! occ-amplitudes v**2
      double precision CPfer(MID,-1:1)        ! fermi energy
      integer          CPlen                  ! number of path elements


c=======================================================================
c
c     common block "workspaces" 
c
c=======================================================================
      common /LRbuff/
     &     LRRes,LRbuffer,LRwrk_scalar,LRwrk2_scalar,
     &     LRwrk1_spinor,LRwrk2_spinor,
     &     LRwrk3_spinor,LRwrk4_spinor
      
      double precision LRRes(0:MZD,0:MRD,MSD,MPSID,2)
      double precision LRbuff1(0:MZD,0:MRD,MSD,2)
      double precision LRbuff2(0:MZD,0:MRD,MSD,2)
      double precision LRwrk_scalar(0:MZD,0:MRD)      ! scalar buffer
      double precision LRwrk2_scalar(0:MZD,0:MRD)     ! scalar buffer
      double precision LRwrk1_spinor(0:MZD,0:MRD,MSD) ! spinor buffers
      double precision LRwrk2_spinor(0:MZD,0:MRD,MSD)
      double precision LRwrk3_spinor(0:MZD,0:MRD,MSD)
      double precision LRwrk4_spinor(0:MZD,0:MRD,MSD)


c=======================================================================
c
c     setting some switches
c
c=======================================================================
      
      logical FULLSPACE, D_LRDEBUG, D_OLDCRANK 

c     parameter (FULLSPACE=0)
      parameter (D_LRDEBUG=.false.)
      parameter (D_OLDCRANK=.false.)

