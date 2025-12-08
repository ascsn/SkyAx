cc#include "define.h"
c=======================================================================
c
c     linresp()
c
c       main routine for the cranking module of SkyAx. The subroutine
c       calculates the inverse masses B of a collective Hamiltonian
c       and the curvature (stiffness) of the potential C
c     
c                                          1   d^2 V
c            H = B P^2 + V(Q)         C =  - * -----
c                                          2   dQ^2
c
c       with
c
c            P  = 1/(2B) * [H,Q]       C = 1/2 * [P,[H,P]]
c            B  = 1/2 * [Q,[H,Q]]
c            C  = 1/2 * [P,[H,P]]
c
c       for different models
c       
c        1) simple mass from multipole moments
c        2) inglis cranking 
c        3) self-consistent crankting (P from CHF)
c        3) ATDHF cranking (self-consistant P,Q-operators)
c
c
c
c       the phonon energy of the mode described by P,Q is
c       defined by
c
c            [H,[H,Q]] = omega^2 Q     E_0 = 1/2 * omega
c
c
c       Q_ij are choosen real ==> P_ij == purely complex
c       the variable "x2q",etc. stores
c
c           q_ij = Q_ij          (if Q time-even)
c           p_ij = -i/hbar P_ij  (if P time-odd)
c
c       If the response of the MF is included the
c       commutators become
c
c           [H,Q] = [H_0,Q] + H^{resp}_Q
c           [H,P] = [H_0,P] + H^{resp}_P
c
c       where H_0 is diagonal in quasiparticle-space
c
c           H_0 = e_i alpha^+_i alpha_i
c       
c       and therefore the 2q-part becomes
c      
c           [H,Q]^{2q} = (e_i + e_j) Q_ij + (H^{resp}_Q)_ij
c           [H,P]^{2q} = (e_i + e_j) P_ij + (H^{resp}_P)_ij
c
c       usefull commutators:
c       
c        <Q^2>     = 1/2 * Q_ij Q_ij       <P^2> = 1/2 * P_ij P_ij
c           
c
c
c
c=======================================================================
      subroutine linresp(outnam)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'linresp.inc'
      include 'wavefr.inc'
      include 'levels.inc'
      include 'observ.inc'
      include 'constr.inc'
c
c     ------------------------------------------------------------------
c
      character*3      outnam
c
      double precision B2q(1:MTDTOT)
      double precision P2q(1:MTDTOT)
      double precision Q2q(1:MTDTOT)

      double precision H2q(1:MTDTOT)
      double precision N2q(1:MTDTOT)
c
      double precision B0q(0:2,0:2)         ! deformation and variance
      double precision Bvar(0:2,0:2)        ! (second index: isospin)
      double precision N0q(0:2)             ! particle number and variance
      double precision H0q
      double precision Nvar(0:2)            ! (second index: isospin)
c   
      double precision Pvar(0:2,0:5)        ! <P^2>  momentum variance
      double precision Qvar(0:2,0:5)        ! <Q^2>  Qdyn     variance

      double precision C(0:2,0:5)           ! Stiffness C 
c
      double precision B(0:2,0:5)           ! Mass  B     
      double precision BN(0:2,0:5)          ! Part. Number Mass  B_N 
c
      double precision Nrm0(0:2,0:5)        ! <[Beta,P]>
      double precision Nrm(0:2,0:5)         ! <[Q,P]>
      double precision OmgSq(0:2,0:5)       ! <[[H,[H,Q],P]>
c
      double precision norm0(0:2),norm0qp(0:2),norm0res(0:2)
      
      double precision d2energy,d2radRMS,d2radDIF,d2radSRF
      double precision d2energyC
      double precision buff(10),buffer
c      
      integer          imu,i,ibuff
      double precision betapow
c
      integer          il
      integer          PathLen
c
c     ------------------------------------------------------------------
c     switches
c     ------------------------------------------------------------------
      isMFresp=1          ! switch to turn on/off MF-response
      il      =2          ! collective mode (2=quadrupole)

      if (swicon.gt.oh) then
      PathLen = 1 ! 4    
      else
      PathLen = 1   ! 4 
      endif

c
c     update the collective path
c     ------------------------------------------------------------------
      call upd_path()
c
c     static afterburn to get a good sp-basis also for unoccupied states
c     ------------------------------------------------------------------
      call static_afterburn(+1,0.01D0,1000)
      call return_wavefunctions(+1)
c
c     update deformations off last path element (to keep consistence)
c     ------------------------------------------------------------------
c      call ini_transtab(il,0)
c      call ini_matrix_B(B2q,CPbeta(+1),buffer,il,0,0,+1)
c
c     output the deformations kept in the path storage 
c     ------------------------------------------------------------------
      write(21,'(a,3f8.4)')         
     &  'LINRESP()    : elements of collective path -- beta= ',
     &                  (CPbeta(i),i=-1,1)
c
c     check if path-length is to small
c     ------------------------------------------------------------------
      if (CPlen.lt.PathLen) then
c      if (CPlen.ne.6) then
         write(21,'(a)')         
     &  'LINRESP()    : SKIPPED linear response  -- '//
     &  'collective path incomplete'
         return
      endif
c
c     prepare derivatives of central element of path 
c     ------------------------------------------------------------------
      call upd_derivatives()
c
c     update non-linear response weight dU_q'/drho_q
c     ------------------------------------------------------------------
      call upd_locpotwgt()
c
c     calculate quadrupole mass for m=0,2
c     ------------------------------------------------------------------
      do imu=0,il,1
c      do imu=0,0,1
c
c        write header to cranking log
c        ------------------------------------------------------------------
         write(21,'(a)')   ' ------------------------------------------'
         write(21,'(a,i2,a,i2)') 
     &                     ' cranking for lambda = ',il,'  mu = ',imu
         write(21,'(a)')   ' ------------------------------------------'
         write(21,'(a,7(a16))') 'METHOD    : ',
     &        'B','C','<P^2>','<Q^2>-Q^2','Omega^2',
     &        '<[Q,P]>','<[alpha,P]>'
         CALL FLUSH(21)
c
c        initialize transitions
c        ---------------------------------------------------------------
         call ini_transtab(il,imu) 
c
c        initialize quasiparticle energies, beta-operator
c        ---------------------------------------------------------------
         call ini_matrix_H(H2q,H0q,il,imu)
         call ini_matrix_B(B2q,B0q(imu,0),BVar(imu,0),il,imu,0,0)
         if (imu.eq.0) then
            call ini_matrix_N(N2q, N0q(0),NVar(0),il,imu,0)
         else
            call ini_matrix_N(N2q,buff(1),buff(2),il,imu,0)
         endif
c
         if (CPLen.eq.PathLen) then
            do itr=1,TRntrtot
               Q2qCHF(itr,imu)  =B2q(itr)/H2q(itr)**2
               P2qCHF(itr,imu)  =B2q(itr)/H2q(itr)
               Q2qATDHF(itr,imu)=B2q(itr)/H2q(itr)**2
               P2qATDHF(itr,imu)=B2q(itr)/H2q(itr)
            enddo
         endif
c
c
c        calculate coll. mass using inglis and self-consistant cranking
c        ---------------------------------------------------------------
         if (imu.eq.0.or.imu.eq.1) then
c
c           initialize P-operator from CHF (SkyAx, constraint loop)
            call ini_matrix_P(P2q,Pvar(imu,1),il,imu)
            call eval_comm_BP(Nrm0(imu,1),il,imu)
            Pvar(imu,2)=Pvar(imu,1)    ! SC-, INGLIS- and static GCM
            Pvar(imu,3)=Pvar(imu,1)    ! have the same P-Variance and 
            Nrm0(imu,2)=Nrm0(imu,1)    ! <[beta,P]> commutator
            Nrm0(imu,3)=Nrm0(imu,1)
c
c           the GCM-mass
            call crank_GCM(B2q,Q2q,P2q,H2q,
     &                        B(imu,1),C(imu,1),
     &                        QVar(imu,1),Pvar(imu,1),
     &                        il,imu)
            call eval_com_QNQ(N2q,Q2q,BN(imu,1),il,imu,0)
            Pvar(imu,1)=Pvar(imu,2)
c
c           inglis cranking mass
            call crank_inglis(B2q,Q2q,P2q,H2q,
     &                        B(imu,2),C(imu,2),
     &                        QVar(imu,2),buff(1),
     &                        il,imu)
            call eval_com_QNQ(N2q,Q2q,BN(imu,2),il,imu,0)
c
c           selfconsistent cranking mass
            call crank_selfc(B2q,Q2q,P2q,H2q,
     &                       B(imu,3),C(imu,3),
     &                       QVar(imu,3),buff(1),
     &                       il,imu)
             call eval_com_QNQ(N2q,Q2q,BN(imu,3),il,imu,0)
c
c           debug subroutine to test the responses
c            call crank_test_resp(B2q,Q2q,P2q,H2q,
c     &                         B(imu,0),C(imu,0),
c     &                         QVar(imu,0),PVar(imu,0),
c     &                         il,imu)
c
c           the scaling mass
c            call crank_scaling(B2q,Q2q,P2q,H2q,
c     &                         B(imu,0),C(imu,0),
c     &                         QVar(imu,0),PVar(imu,0),
c     &                         il,imu)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,0),il,imu,0)
cc
cc           CHF in linear response
c            call crank_CHF(B2q,Q2q,P2q,H2q,
c     &                     B(imu,4),C(imu,4),
c     &                     QVar(imu,4),PVar(imu,4),
c     &                     il,imu)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
c            call crank_ATDHF(B2q,Q2q,P2q,H2q,
c     &                     B(imu,4),C(imu,4),
c     &                     QVar(imu,4),PVar(imu,4),
c     &                     il,imu,1)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
cc
cc           full ATDHF
c            call crank_ATDHF(B2q,Q2q,P2q,H2q,
c     &                     B(imu,5),C(imu,5),
c     &                     QVar(imu,5),PVar(imu,5),
c     &                     il,imu,0)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,1),il,imu,0)
c
         else if (imu.gt.1) then
c
c           debug subroutine to test the responses
c            call crank_test_resp(B2q,Q2q,P2q,H2q,
c     &                         B(imu,0),C(imu,0),
c     &                         QVar(imu,0),PVar(imu,0),
c     &                         il,imu)
c
c           the scaling mass
c            call crank_scaling(B2q,Q2q,P2q,H2q,
c     &                         B(imu,0),C(imu,0),
c     &                         QVar(imu,0),PVar(imu,0),
c     &                         il,imu)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,0),il,imu,0)
cc
cc           CHF in linear response
c            call crank_CHF(B2q,Q2q,P2q,H2q,
c     &                     B(imu,4),C(imu,4),
c     &                     QVar(imu,4),PVar(imu,4),
c     &                     il,imu)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
c            call crank_ATDHF(B2q,Q2q,P2q,H2q,
c     &                     B(imu,4),C(imu,4),
c     &                     QVar(imu,4),PVar(imu,4),
c     &                     il,imu,1)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,4),il,imu,0)
cc
cc           the GCM-mass
cc            call crank_GCM(B2q,Q2q,P2q,H2q,
cc     &                        B(imu,1),C(imu,1),
cc     &                        QVar(imu,1),buff(1),
cc     &                        il,imu)
cc            call eval_com_QNQ(N2q,Q2q,BN(imu,1),il,imu,0)
cc
cc           full ATDHF
c            call crank_ATDHF(B2q,Q2q,P2q,H2q,
c     &                     B(imu,5),C(imu,5),
c     &                     QVar(imu,5),PVar(imu,5),
c     &                     il,imu,0)
c            call eval_com_QNQ(N2q,Q2q,BN(imu,5),il,imu,0)
c
         endif
c
c        calculate widths (and deformations) for particle sorts seperatly 
c        ---------------------------------------------------------------
         call ini_matrix_B(B2q,B0q(imu,1),BVar(imu,1),il,imu,1,0) ! protons
         call ini_matrix_B(B2q,B0q(imu,2),BVar(imu,2),il,imu,2,0) ! neutrons
c
c         PVar(imu,2)=PVar(imu,1)     ! <P^2> is equal for GCM,SC,INGLIS
c         PVar(imu,3)=PVar(imu,1)     !
c         Nrm0(imu,2)=Nrm0(imu,1)     ! <[Q,P]> is equal for GCM,SC,INGLIS
c         Nrm0(imu,3)=Nrm0(imu,1)
c
c        calculate particle numbers and variances 
c        ---------------------------------------------------------------
         if (imu.eq.0) then
            call ini_matrix_N(N2q,N0q(0),NVar(0),il,imu,0)
            call ini_matrix_N(N2q,N0q(1),NVar(1),il,imu,1)
            call ini_matrix_N(N2q,N0q(2),NVar(2),il,imu,2)
         endif
c
      enddo
c
c     write collective parameters of Hamiltonian
c     ------------------------------------------------------------------

      do i=0,5
!         write(23+i,'(f8.5,21(1pg12.4))') B0q(0,0),
         write(23+i,'(f8.5,21(1pg12.4))') CPbeta(0),
     &           B(0,i), C(0,i), PVar(0,i), QVar(0,i), OmgSq(0,i), 
     &           B(1,i), C(1,i), PVar(1,i), QVar(1,i), OmgSq(1,i),
     &           B(2,i), C(2,i), PVar(2,i), QVar(2,i), OmgSq(2,i),
     &           BN(0,i),BN(1,i),BN(2,i),
     &           Nrm0(0,i),Nrm0(1,i),Nrm0(2,i)
         CALL FLUSH(23+i)
      enddo
c
c     calculate curvature directly from path
c     ------------------------------------------------------------------
      d2energy        = (CPenergy(+1) - 2.0*CPenergy( 0) +CPenergy(-1))
     &                 /(0.5*(CPbeta( +1) - CPbeta(-1)))**2/2.0
      d2energyC       = (CPenergyC(+1)-2.0*CPenergyC( 0)+CPenergyC(-1))
     &                 /(0.5*(CPbeta( +1) - CPbeta(-1)))**2/2.0
      d2radRMS        = (CPradRMS(+1) - 2.0*CPradRMS( 0) +CPradRMS(-1))
     &                 /(0.5*(CPbeta( +1) - CPbeta(-1)))**2/2.0
      d2radDIF        = (CPradDIF(+1) - 2.0*CPradDIF( 0) +CPradDIF(-1))
     &                 /(0.5*(CPbeta( +1) - CPbeta(-1)))**2/2.0
      d2radSRF        = (CPradSRF(+1) - 2.0*CPradSRF( 0) +CPradSRF(-1))
     &                 /(0.5*(CPbeta( +1) - CPbeta(-1)))**2/2.0
c
c     write the radii
c     -----------------------------------------------------------------
!         write(29,'(6(1pg14.5))') 
!     &     B0q(0,0),CPradRMS(0),CPradDIF(0),CPradSRF(0),
         write(29,'(6(1pg14.5))') 
     &     CPbeta(0),CPradRMS(0),CPradDIF(0),CPradSRF(0),
     &     d2energy, d2energyC
c
c     write the deformations and deformation-widths (and part. number)
c     -----------------------------------------------------------------
!         write(30,'(20(1pg14.5))') 
!     &     B0q(0,0),B0q(0,1),B0q(0,2),
         write(30,'(20(1pg14.5))') 
     &     CPbeta(0),B0q(0,1),B0q(0,2),
     &     Bvar(0,0),Bvar(1,0),Bvar(2,0),
     &     Bvar(0,1),Bvar(1,1),Bvar(2,1),
     &     Bvar(0,2),Bvar(1,2),Bvar(2,2),
     &     N0q(0),N0q(1),N0q(2),
     &     Nvar(0),Nvar(1),Nvar(2),
     &     Pvar(0,1),Pvar(1,1)
c
c     write the normalisation 
c     -----------------------------------------------------------------
!         write(32,'(9(1pg14.5))') 
!     &     B0q(0,0),
         write(32,'(9(1pg14.5))') 
     &     CPbeta(0),
     &     Nrm0(0,1),Nrm0(1,1),          ! GCM and CHF  [Beta,P]
     &     Nrm0(0,0),Nrm0(1,0),Nrm0(2,0),! Multipole    [Beta,P]
     &     Nrm0(0,5),Nrm0(1,5),Nrm0(2,5) ! ADTHF        [Beta,P]
c
c     ------------------------------------------------------------------
c
        CALL FLUSH(29)
        CALL FLUSH(30)
        CALL FLUSH(32)
      return
      end


c=======================================================================
c
c     cranking()
c
c     calculate dynamic Q-operator and P-operator by different cranking
c     methods.
c     
c     Different Methods are supported:
c
c       0)  CHF:           [H,[H,Q]]             = -2*dlambda/dq*B Beta
c       1)  GCM:                 Q               =   -i P/(2*<P^2>)
c       2)  inglis:           [H,Q]              =   i  2B     P
c       3)  self-cons.:   LQ= [H,Q] + H^{resp}_Q =   i  2B     P
c       4)  ATDHF:            L^2Q               =   Omega^2   Q
c
c     H is the static Hamiltonian, Beta is the messurement operator
c     of the collective coordiante. L is the Liouvillian
c
c                      LA = [H,A] + H^{resp}_A
c
c     While CHF and ATDHF locally define the collective path by
c     themselves GCM, inglis and self-cons. require a given P-Operator.
c
c=======================================================================
c      subroutine cranking(type,imu,H2q,
c     &                    P2q,resP,Q2q,resQ,B2q,resB,
c     &                    B,BN,C,PVar,QVar,OmgSq,Nrm,Nrm0)
cc
cc     ------------------------------------------------------------------
cc
c      include 'paramr.inc'
c      include 'dynpar.inc'
c      include 'linresp.inc'
c      include 'levels.inc'
cc
cc     ------------------------------------------------------------------
cc
c      double precision H2q(1:MTD)
c      double precision P2q(1:MTD) 
c      double precision Q2q(1:MTD) 
c      double precision B2q(1:MTD) 
c      double precision HP2q(1:MTD)
c      double precision HQ2q(1:MTD)
c#IF(FULLSPACE)
c      double precision ResP(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision ResQ(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision ResB(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision ResHP(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision ResHQ(0:MZD,0:MRD,MSD,MPSID,2)
c
c      double precision drRes(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision dzRes(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision dpRes(0:MZD,0:MRD,MSD,MPSID,2)
c      double precision ddRes(0:MZD,0:MRD,MSD,MPSID,2)
c#ELSE
c      double precision ResP
c      double precision ResQ
c      double precision ResB
c      double precision ResHP
c      double precision ResHQ
c
c      double precision drRes
c      double precision dzRes
c      double precision dpRes
c      double precision ddRes
c#ENDIF
c
c      double precision B,BN,C
c      double precision PVar,BVar,QVar,OmgSq
c      double precision Nrm,Nrm0
c      integer          imu,type
c
c      integer          max_citer
c      double precision stp_res_Q  , stp_qp_Q  , max_err_Q 
c      double precision stp_res_P  , stp_qp_P  , max_err_P
c      double precision stp_res_Omg, stp_qp_Omg, max_err_Omg
c
c      integer          ir,iz,itr,ist,iter,buff2,lograte
c      integer          ibpar
c      double precision norm,error_P,error_Q,buffer
c      double precision var1,var2
c      double precision Pnorm(MPSID),Qnorm(MPSID),Bnorm(MPSID)
c      double precision HPnorm(MPSID),HQnorm(MPSID)
c      double precision rho,rhomax,r,z,buff,buff1
c
cc
cc     ------------------------------------------------------------------
cc
c      error_P  =0.0D0
c      error_Q  =0.0D0
c      error_Omg=0.0D0
cc      
cc     iteration parameters
cc     ------------------------------------------------------------------
c      max_citer    = 1000   ! maximum number of cranking iterations
c      stp_res_P    = LRStep ! gradient steps in residual components
c      stp_res_Q    = LRStep 
c      stp_res_Omg  = LRStep**2  
c      stp_qp_P     = 0.8    ! gradient steps in qp-components
c      stp_qp_Q     = 0.8    
c      stp_qp_Omg   = 0.8
c      max_err_P    = 1.0E-12 ! maximum relative error || P - LQ   ||/<P^2>
c      max_err_Q    = 1.0E-12 ! maximum relative error || Q - LP   ||/<Q^2>
c      max_err_Omg  = 1.0E-12 ! maximum relative error || Q - L^2Q ||/<Omg^2>
c      lograte      = 1       ! logging rate
cc
cc     check for unsupported cranking method
cc     ------------------------------------------------------------------
c      if (type.lt.0.or.type.gt.5) stop 'unsupported cranking method'
cc
cc     ==================================================================
cc     now the five different methods:
cc     ==================================================================
cc
cc     ------------------------------------------------------------------
cc     CHF + SELFCONSISTANT CRANKING                     (requires Beta) 
cc     ------------------------------------------------------------------
c      if (type.eq.0) then
cc
cc        initialization
cc        ---------------------------------------------------------------    
c         BVar   = var_oper(B2q,ResB,imu)
cc     
cc        find Q-operator by {\cal L}^2 Q = 4 BC Q
cc        ---------------------------------------------------------------
c         do itr=1,TRlen
c            Q2q(itr)= B2q(itr)/H2q(itr)**2 ! Q_ij = B_ij/(e_i + e_j)**2
c         enddo
cc
c#IF(FULLSPACE)
cc        inversion of {\cal L}_0 in residual space 
cc        ---------------------------------------------------------------
cc        initialize using the B-Operator
c         if (isresact.eq.1) then
cc           
cc           initialize using the B-Operator
c            do ist=1,nst(1)+nst(2)
c               do ir=0,nr
c                  do iz=0,nz
c                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
c                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
c                   ResQ(iz,ir,1,ist,2)=ResB(iz,ir,1,ist,2)
c                   ResQ(iz,ir,2,ist,2)=ResB(iz,ir,2,ist,2)
c                  enddo
c               enddo
c            enddo
cc
cc           calculate norm of residual componentes B-Operator
c            do ist=1,nst(1)+nst(2)
c               Bnorm(ist)=
c     &              wfovlp(ResB(0,0,1,ist,1),ResB(0,0,1,ist,1))
c     &             +wfovlp(ResB(0,0,1,ist,2),ResB(0,0,1,ist,2))
c            enddo            
cc
c         else
cc
cc           set residual components to zero
c            call zero_oper(Q2q,ResQ,2)  
c            call zero_oper(HP2q,ResHP,2)           
cc
c         endif
c#ENDIF
cc
cc        ===============================================================
cc        the cranking iteration (find Q-Operator)
cc        ===============================================================
c         do iter=1,max_citer
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           calculate derivatives of resiual part of Q
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
c               call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                              ddRes(0,0,1,1,1),
c     &                              (-1)**imu,+imu,nst(1)+nst(2))
c               call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                              ddRes(0,0,1,1,2),
c     &                              (-1)**imu,-imu,nst(1)+nst(2))
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           calculate action of Liouville-Operator
cc           ------------------------------------------------------------
c            call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
c     &                    HQ2q,ResHQ,imu)
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           calculate derivatives of resiual part of LQ \propto P
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
c               call upd_derive_spinor(ResHQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                              ddRes(0,0,1,1,1),
c     &                              (-1)**imu,+imu,nst(1)+nst(2))
c               call upd_derive_spinor(ResHQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                              ddRes(0,0,1,1,2),
c     &                              (-1)**imu,-imu,nst(1)+nst(2))
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           calculate action of Liouville-Operator
cc           ------------------------------------------------------------
c            call comm_HP(H2q,HQ2q,ResHQ,drRes,dzRes,dpRes,ddRes,
c     &                    HP2q,ResHP,imu)
cc
cc           ------------------------------------------------------------ 
cc           perform the gradient step in qp-space
cc           ------------------------------------------------------------
c            do itr=1,TRLen
c               HQ2q(itr) = HP2q(itr)-B2q(itr)
c               Q2q(itr)  = 
c     &         Q2q(itr)  - stp_qp_Omg*(HP2q(itr)-B2q(itr))/H2q(itr)**2
c            enddo
cc
cc           ------------------------------------------------------------
cc           evaluate errors (only if residual parts are neglected)
cc           ------------------------------------------------------------
c            if (isresact.eq.0) then
c               buffer   = error_Q
c               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
c               derror_Q = error_Q-buffer
c            endif            
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           perform the gradient step in residual space
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
cc
cc              --------------------------------------------------------- 
cc              calculate norm of residual components of [H,Q]-Operator
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  HPnorm(ist)=
c     &                 wfovlp(resHP(0,0,1,ist,1),resHP(0,0,1,ist,1))+
c     &                 wfovlp(resHP(0,0,1,ist,2),resHP(0,0,1,ist,2))
c               enddo
cc
cc              --------------------------------------------------------- 
cc              calculate errors in residual space
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  do ir=0,nr
c                     do iz=0,nz
c                        ResHQ(iz,ir,1,ist,1) = ResHP(iz,ir,1,ist,1)-
c     &                                         ResB(iz,ir,1,ist,1)
c                        ResHQ(iz,ir,2,ist,1) = ResHP(iz,ir,2,ist,1)-
c     &                                         ResB(iz,ir,2,ist,1)
c                        ResHQ(iz,ir,1,ist,2) = ResHP(iz,ir,1,ist,2)-
c     &                                         ResB(iz,ir,1,ist,2)
c                        ResHQ(iz,ir,2,ist,2) = ResHP(iz,ir,2,ist,2)-
c     &                                         ResB(iz,ir,2,ist,2)
c                     enddo
c                  enddo
c               enddo
cc
c               buffer   = error_Q
c               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
c               derror_Q = error_Q-buffer
cc
cc              --------------------------------------------------------- 
cc              apply damping
cc              ---------------------------------------------------------
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHQ(0,0,1,1,1),ResHP(0,0,1,1,1),
c     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHP(0,0,1,1,1),ResHQ(0,0,1,1,1),
c     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHQ(0,0,1,1,2),ResHP(0,0,1,1,2),
c     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHP(0,0,1,1,2),ResHQ(0,0,1,1,2),
c     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))
c
cc
cc              --------------------------------------------------------- 
cc              the gradient step
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  do ir=0,nr
c                     do iz=0,nz
c                      ResQ(iz,ir,1,ist,1) = 
c     &                ResQ(iz,ir,1,ist,1) - ResHQ(iz,ir,1,ist,1)
c                      ResQ(iz,ir,2,ist,1) = 
c     &                ResQ(iz,ir,2,ist,1) - ResHQ(iz,ir,2,ist,1)
c                      ResQ(iz,ir,1,ist,2) = 
c     &                ResQ(iz,ir,1,ist,2) - ResHQ(iz,ir,1,ist,2)
c                      ResQ(iz,ir,2,ist,2) = 
c     &                ResQ(iz,ir,2,ist,2) - ResHQ(iz,ir,2,ist,2)
c                     enddo
c                  enddo
c               enddo
cc
cc              --------------------------------------------------------- 
cc              norm stabilization
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  norm=sqrt(Bnorm(ist)/HPnorm(ist))
c                  do ir=0,nr
c                     do iz=0,nz
c                        ResQ(iz,ir,1,ist,1) = ResQ(iz,ir,1,ist,1)*norm
c                        ResQ(iz,ir,2,ist,1) = ResQ(iz,ir,2,ist,1)*norm
c                        ResQ(iz,ir,1,ist,2) = ResQ(iz,ir,1,ist,2)*norm
c                        ResQ(iz,ir,2,ist,2) = ResQ(iz,ir,2,ist,2)*norm
c                     enddo
c                  enddo
c               enddo
cc
cc              --------------------------------------------------------- 
cc              project out sp-contributions
cc              ---------------------------------------------------------
c               call project_Omega(ResQ,(-1)**imu,imu,1,'Q')
cc
cc           ------------------------------------------------------------ 
cc           END of gradient step in residual space
cc           ------------------------------------------------------------
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           Logging and convergency-check
cc           ------------------------------------------------------------
cc     
c            if (mod(iter,lograte).eq.0               .or.
c     &          error_Q.le.max_err_Q*sqrt(abs(BVar)) .or.
c     &          iter.eq.max_citer ) then
cc     
c              write(21,'(a,i5,1(1pg16.4),2(e16.3))') 'CHF+SELFC  : ',
c     &              iter,error_Q,error_Q/sqrt(abs(BVar)),
c     &              derror_Q
cc
c               if (error_Q.le.max_err_Q*sqrt(abs(BVar))) goto 2100
c               if (abs(derror_Q).lt.1.0E-10) goto 2100
cc     
c            endif
cc
c         enddo
cc        ===============================================================
cc        END of cranking iteration (find Q-Operator)
cc        ===============================================================
c 2100    continue    
cc
c#IF(FULLSPACE)
cc        ------------------------------------------------------------ 
cc        calculate derivatives of residual part of Q
cc         ------------------------------------------------------------
c         if (isresact.eq.1) then
c            call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                           dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                           ddRes(0,0,1,1,1),
c     &                           (-1)**imu,+imu,nst(1)+nst(2))
c            call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                            dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                            ddRes(0,0,1,1,2),
c     &                            (-1)**imu,-imu,nst(1)+nst(2))
c         endif
c#ENDIF
cc
cc        ------------------------------------------------------------ 
cc        calculate action of Liouville-Operator to obtain P-Operator
cc        ------------------------------------------------------------
c         call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
c     &                P2q,ResP,imu)
cc
cc     ------------------------------------------------------------------
cc     GCM-formula:                                    (requires P,Beta) 
cc     ------------------------------------------------------------------
c      else if (type.eq.1) then
cc
cc        GCM-Mass requires P-Variance
c#IF(FULLSPACE)
c         PVar   = var_oper(P2q,ResP,imu)
c#ELSE
c         call comm_and_var_QP(buff,PVar,buff1,CPbeta(0),
c     &                        (-1)**imu,imu,0)
c#ENDIF
cc
cc        the quasiparticle part of the GCM Q-Operator
c         do itr=1,TRlen
c            Q2q(itr)= P2q(itr)/(2.0D0*PVar)
c         enddo
cc
c#IF(FULLSPACE)
cc        the residual part of the GCM Q-Operator
c         if (isresact.eq.1) then
c            do ist=1,nst(1)+nst(2)
c               do ir=0,nr
c                  do iz=0,nz
c                  ResQ(iz,ir,1,ist,1)=ResP(iz,ir,1,ist,1)/(2.0D0*PVar)
c                  ResQ(iz,ir,2,ist,1)=ResP(iz,ir,2,ist,1)/(2.0D0*PVar)
c                  ResQ(iz,ir,1,ist,2)=ResP(iz,ir,1,ist,2)/(2.0D0*PVar)
c                  ResQ(iz,ir,2,ist,2)=ResP(iz,ir,2,ist,2)/(2.0D0*PVar)
c                  enddo
c               enddo
c            enddo
c         else
c            call zero_oper(Q2q,ResQ,2)  ! set residual part to zero
c         endif
c#ENDIF
cc
cc     ------------------------------------------------------------------
cc     inglis-formula:                                 (requires P,Beta) 
cc     ------------------------------------------------------------------
c      else if (type.eq.2) then
cc
cc        temporary disable MF-response
cc        ---------------------------------------------------------------
c         buff2=isMFresp
c         isMFresp=0
cc
cc        initialization
cc        ---------------------------------------------------------------    
c#IF(FULLSPACE)
c         PVar   = var_oper(P2q,ResP,imu)
c#ELSE
c         call comm_and_var_QP(buff,PVar,buff1,CPbeta(0),
c     &                        (-1)**imu,imu,0)
c#ENDIF
cc     
cc        find Q-operator by {\cal L}_0 Q = [H,Q] = -2i B P
cc        ---------------------------------------------------------------   
cc        initialize with inglis values for the 2qp-part
c         do itr=1,TRlen
c            Q2q(itr)= P2q(itr)/H2q(itr)   ! Q_ij = P_ij/(e_i + e_j)
c         enddo
cc
c#IF(FULLSPACE)
cc        inversion of {\cal L}_0 in residual space 
cc        ---------------------------------------------------------------
cc        initialize using the B-Operator
c         if (isresact.eq.1) then
c            do ist=1,nst(1)+nst(2)
c               do ir=0,nr
c                  do iz=0,nz
c                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
c                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
c                   ResQ(iz,ir,1,ist,2)=0.0D0*ResB(iz,ir,1,ist,2)
c                   ResQ(iz,ir,2,ist,2)=0.0D0*ResB(iz,ir,2,ist,2)
c                  enddo
c               enddo
c            enddo
c         else
c            call zero_oper(Q2q,ResQ,2)  ! set residual part to zero
c         endif
cc
cc        calculate norm of residual componentes P-Operator
c         do ist=1,nst(1)+nst(2)
c            Pnorm(ist)=
c     &           wfovlp(ResP(0,0,1,ist,1),ResP(0,0,1,ist,1))
cc     &          +wfovlp(ResP(0,0,1,ist,2),ResP(0,0,1,ist,2))
c         enddo
cc
cc        find resQ by iteration
c         if (isresact.eq.1) then
c            do iter=1,max_citer
cc
cc            prepare derivatives of residual part of Q
c             call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                              ddRes(0,0,1,1,1),
c     &                              (-1)**imu,+imu,nst(1)+nst(2))
cc             call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
cc     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
cc     &                              ddRes(0,0,1,1,2),
cc     &                              (-1)**imu,-imu,nst(1)+nst(2))
cc           
cc            calculate action of Liouville-Operator
c             call comm_HQ(H2q,Q2q,resQ,drRes,dzRes,dpRes,ddRes,
c     &                    HQ2q,resHQ,imu)
cc
cc            calculate variance <[H,Q]^2> and <P[H,Q]>
c             var1=var_oper(HQ2q,ResHQ,imu)
c             var2=0.5D0*comm_QP(HQ2q,ResHQ,P2q,ResP,imu)
cc
cc            calculate norm of residual components of [H,Q]-Operator
c             do ist=1,nst(1)+nst(2)
c                HQnorm(ist)=
c     &               wfovlp(resHQ(0,0,1,ist,1),resHQ(0,0,1,ist,1))
cc     &               +wfovlp(resHQ(0,0,1,ist,2),resHQ(0,0,1,ist,2))
c             enddo
cc
cc            calculate error in residual space
c             do ist=1,nst(1)+nst(2)
c                do ir=0,nr
c                   do iz=0,nz
c                      resHP(iz,ir,1,ist,1) = resHQ(iz,ir,1,ist,1)-
c     &                                       resP(iz,ir,1,ist,1)
c                      resHP(iz,ir,2,ist,1) = resHQ(iz,ir,2,ist,1)-
c     &                                       resP(iz,ir,2,ist,1)
c                      resHP(iz,ir,1,ist,2) = 0.0D0!resHQ(iz,ir,1,ist,2)-
cc     &                                       resP(iz,ir,1,ist,2)
c                      resHP(iz,ir,2,ist,2) = 0.0D0!resHQ(iz,ir,2,ist,2)-
cc     &                                       resP(iz,ir,2,ist,2)
c                   enddo
c                enddo
c             enddo
cc
cc            store error of [H,Q]-Operator
c             buffer   = error_P
c             error_P  = sqrt(var_oper(HP2q,resHP,imu))
c             derror_P = error_P-buffer
cc
cc            log components which have big errors
c             do ist=1,nst(1)+nst(2)
c                norm=wfovlp(resHP(0,0,1,ist,1),resHP(0,0,1,ist,1))
cc     &              +wfovlp(resHP(0,0,1,ist,2),resHP(0,0,1,ist,2))
c                if (norm.gt.1E-2) then
c                write (88,*) ist,norm,CPvsqf(ist,0)
c                endif
c             enddo     
cc
cc            apply damping          
c             call damp_spinor(LRDampr,LRDampz,
c     &                        resHP(0,0,1,1,1),resHQ(0,0,1,1,1),
c     &                        (-1)**imu,imu,stp_res_P,nst(1)+nst(2))
cc             call damp_spinor(LRDampr,LRDampz,
cc     &                        resHP(0,0,1,1,2),resHQ(0,0,1,1,2),
cc     &                        (-1)**imu,-imu,stp_res_P,nst(1)+nst(2))
c             call project_Omega(ResHQ,(-1)**imu,imu)
cc
cc            perform the gradient step in residual space
c             do ist=1,nst(1)+nst(2)
c                do ir=0,nr
c                   do iz=0,nz
c                      ResQ(iz,ir,1,ist,1) = 
c     &                ResQ(iz,ir,1,ist,1) - resHQ(iz,ir,1,ist,1)
c                      ResQ(iz,ir,2,ist,1) = 
c     &                ResQ(iz,ir,2,ist,1) - resHQ(iz,ir,2,ist,1)
c                      ResQ(iz,ir,1,ist,2) = 0.0D0!
cc     &                ResQ(iz,ir,1,ist,2) - resHQ(iz,ir,1,ist,2)
c                      ResQ(iz,ir,2,ist,2) = 0.0D0!
cc     &                ResQ(iz,ir,2,ist,2) - resHQ(iz,ir,2,ist,2)
c                      enddo
c                   enddo
c             enddo
cc
cc            norm stabilization
c             do ist=1,nst(1)+nst(2)
c                norm=sqrt(Pnorm(ist)/HQnorm(ist))
c                do ir=0,nr
c                   do iz=0,nz
c                      ResQ(iz,ir,1,ist,1) = ResQ(iz,ir,1,ist,1)*norm
c                      ResQ(iz,ir,2,ist,1) = ResQ(iz,ir,2,ist,1)*norm
c                      ResQ(iz,ir,1,ist,2) = 0.0*ResQ(iz,ir,1,ist,2)*norm
c                      ResQ(iz,ir,2,ist,2) = 0.0*ResQ(iz,ir,2,ist,2)*norm
c                   enddo
c                enddo
c             enddo
cc
cc            project out sp-contributions
c             call project_Omega(ResQ,(-1)**imu,imu)
cc
cc            Logging and convergency-check
cc            ---------------------------------------------------------
c             norm     = 1.0D0-var2**2/(Pvar*var1)
cc     
c             if (mod(iter,lograte).eq.0 .or.
c     &           error_P.le.max_err_P*sqrt(abs(PVar)) .or.
c     &           iter.eq.max_citer ) then
cc     
c               write(21,'(a,i5,2(1pg16.4),3(e16.3))') 'INGLIS   : ',
c     &              iter,error_P,error_P/sqrt(abs(PVar)),
c     &              derror_P,norm
cc
c               if ((error_P.le.max_err_P*sqrt(abs(PVar))).and.
c     &             (norm.lt.1E-4)) goto 2102
c               if (abs(derror_P).lt.1.0E-10) goto 2102
cc     
c              endif
cc
c            enddo
c 2102       continue     
cc
c         endif
c#ENDIF
cc
cc        reactivate MF-response
cc        ---------------------------------------------------------------
c         isMFresp=buff2
cc
cc     ------------------------------------------------------------------
cc     self-consistent cranking:                       (requires P,Beta) 
cc     ------------------------------------------------------------------
c      else if (type.eq.3) then
cc
c#IF(FULLSPACE)
c         PVar   = var_oper(P2q,ResP,imu)
c#ELSE
c         call comm_and_var_QP(buff,PVar,buff1,CPbeta(0),
c     &                        (-1)**imu,imu,0)
c         write(666,*) PVar
c#ENDIF
cc     
cc        find Q-operator by {\cal L}_0 Q = [H,Q] = -2i B P
cc        ---------------------------------------------------------------   
cc        initialize with inglis values for the 2qp-part
c         do itr=1,TRlen
c            Q2q(itr)= P2q(itr)/H2q(itr)   ! Q_ij = P_ij/(e_i + e_j)
c         enddo
cc
c#IF(FULLSPACE)
cc        inversion of {\cal L}_0 in residual space 
cc        ---------------------------------------------------------------
cc        initialize using the B-Operator
c         if (isresact.eq.1) then
cc           
cc           initialize using the B-Operator
c            do ist=1,nst(1)+nst(2)
c               do ir=0,nr
c                  do iz=0,nz
c                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
c                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
c                   ResQ(iz,ir,1,ist,2)=ResB(iz,ir,1,ist,2)
c                   ResQ(iz,ir,2,ist,2)=ResB(iz,ir,2,ist,2)
c                  enddo
c               enddo
c            enddo
cc
cc           calculate norm of residual componentes P-Operator
c            do ist=1,nst(1)+nst(2)
c               Pnorm(ist)=
c     &              wfovlp(ResP(0,0,1,ist,1),ResP(0,0,1,ist,1))
c     &             +wfovlp(ResP(0,0,1,ist,2),ResP(0,0,1,ist,2))
c            enddo            
cc
c         else
cc
cc           set residual components to zero
c            call zero_oper(Q2q,ResQ,2)  
c            call zero_oper(HP2q,ResHP,2)           
cc
c         endif
c#ENDIF
cc
cc        ===============================================================
cc        the cranking iteration (find Q-Operator)
cc        ===============================================================
c         do iter=1,max_citer
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           calculate derivatives of resiual part of Q
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
c               call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                              ddRes(0,0,1,1,1),
c     &                              (-1)**imu,+imu,nst(1)+nst(2))
c               call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                              ddRes(0,0,1,1,2),
c     &                              (-1)**imu,-imu,nst(1)+nst(2))
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           calculate action of Liouville-Operator
cc           ------------------------------------------------------------
c            call comm_HQ(H2q,Q2q,resQ,drRes,dzRes,dpRes,ddRes,
c     &                    HQ2q,resHQ,imu)
cc
cc           ------------------------------------------------------------ 
cc           perform the gradient step in qp-space
cc           ------------------------------------------------------------
c            do itr=1,TRLen
c               HP2q(itr) = HQ2q(itr)-P2q(itr)
c               Q2q(itr)  = 
c     &         Q2q(itr)  - stp_qp_P*(HQ2q(itr)-P2q(itr))/H2q(itr)
c            enddo
cc
cc           ------------------------------------------------------------
cc           evaluate errors (only if residual parts are neglected)
cc           ------------------------------------------------------------
c            if (isresact.eq.0) then
c               buffer   = error_P
c               var1     = var_oper(HQ2q,resHQ,imu)
c               var2     = var_oper( P2q, resP,imu)
c               var3     = 0.5D0*comm_QP(HQ2q,resHQ,P2q,resP,imu)
c               error_P  = sqrt(var_oper(HP2q,resHP,imu))
c               derror_P = error_P-buffer
c            endif            
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           perform the gradient step in residual space
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
cc
cc              --------------------------------------------------------- 
cc              calculate norm of residual components of [H,Q]-Operator
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  HQnorm(ist)=
c     &                 wfovlp(resHQ(0,0,1,ist,1),resHQ(0,0,1,ist,1))+
c     &                 wfovlp(resHQ(0,0,1,ist,2),resHQ(0,0,1,ist,2))
c               enddo
cc
cc              --------------------------------------------------------- 
cc              calculate errors in residual space
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  do ir=0,nr
c                     do iz=0,nz
c                        ResHP(iz,ir,1,ist,1) = ResHQ(iz,ir,1,ist,1)+
c     &                                         ResP(iz,ir,1,ist,1)
c                        ResHP(iz,ir,2,ist,1) = ResHQ(iz,ir,2,ist,1)+
c     &                                         ResP(iz,ir,2,ist,1)
c                        ResHP(iz,ir,1,ist,2) = ResHQ(iz,ir,1,ist,2)+
c     &                                         ResP(iz,ir,1,ist,2)
c                        ResHP(iz,ir,2,ist,2) = ResHQ(iz,ir,2,ist,2)+
c     &                                         ResP(iz,ir,2,ist,2)
c                     enddo
c                  enddo
c               enddo
cc
c               buffer   = error_P
c               error_P  = sqrt(var_oper(HP2q,resHP,imu))
c               derror_P = error_P-buffer
cc
cc              --------------------------------------------------------- 
cc              apply damping
cc              ---------------------------------------------------------
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHP(0,0,1,1,1),ResHQ(0,0,1,1,1),
c     &                          (-1)**imu,+imu,stp_res_P,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHP(0,0,1,1,2),ResHQ(0,0,1,1,2),
c     &                          (-1)**imu,-imu,stp_res_P,nst(1)+nst(2))
cc
cc              --------------------------------------------------------- 
cc              the gradient step
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  do ir=0,nr
c                     do iz=0,nz
c                      ResQ(iz,ir,1,ist,1) = 
c     &                ResQ(iz,ir,1,ist,1) + resHQ(iz,ir,1,ist,1)
c                      ResQ(iz,ir,2,ist,1) = 
c     &                ResQ(iz,ir,2,ist,1) + resHQ(iz,ir,2,ist,1)
c                      ResQ(iz,ir,1,ist,2) = 
c     &                ResQ(iz,ir,1,ist,2) + resHQ(iz,ir,1,ist,2)
c                      ResQ(iz,ir,2,ist,2) =
c     &                ResQ(iz,ir,2,ist,2) + resHQ(iz,ir,2,ist,2)
c                     enddo
c                  enddo
c               enddo
cc
cc              --------------------------------------------------------- 
cc              norm stabilization
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  norm=sqrt(Pnorm(ist)/HQnorm(ist))
c                  do ir=0,nr
c                     do iz=0,nz
c                        ResQ(iz,ir,1,ist,1) = ResQ(iz,ir,1,ist,1)*norm
c                        ResQ(iz,ir,2,ist,1) = ResQ(iz,ir,2,ist,1)*norm
c                        ResQ(iz,ir,1,ist,2) = ResQ(iz,ir,1,ist,2)*norm
c                        ResQ(iz,ir,2,ist,2) = ResQ(iz,ir,2,ist,2)*norm
c                     enddo
c                  enddo
c               enddo
cc
cc              --------------------------------------------------------- 
cc              project out sp-contributions
cc              ---------------------------------------------------------
c               call project_qp(Q2q,ResQ,+1,imu,2)
cc               call project_Omega(ResQ,(-1)**imu,imu)
cc
cc           ------------------------------------------------------------ 
cc           END of gradient step in residual space
cc           ------------------------------------------------------------
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           Logging and convergency-check
cc           ------------------------------------------------------------    
c            if (mod(iter,lograte).eq.0               .or.
cc     &          error_P.le.max_err_P*sqrt(abs(PVar)) .or.
c     &          iter.eq.max_citer ) then
cc     
c               write(21,'(a,i5,4(f6.1))') 'SELFC    : ',
c     &              iter,
c     &              log10(MAX(abs(error_P                  ),1E-16)),
c     &              log10(MAX(abs(error_P/sqrt(abs(PVar))  ),1E-16)),
c     &              log10(MAX(abs(derror_P                 ),1E-16)),
c     &              log10(MAX(abs(1.0D0-var3**2/(var1*var2)),1E-16))
cc
cc               if (error_P.le.max_err_P*sqrt(abs(PVar))) goto 2103
cc     
c            endif
cc
c            if (error_P.le.max_err_P) goto 2103
cc
c         enddo
cc        ===============================================================
cc        END of cranking iteration (find Q-Operator)
cc        ===============================================================
c 2103    continue     
cc
cc
cc     ------------------------------------------------------------------
cc     ATDHF                                              (requires Beta) 
cc     ------------------------------------------------------------------
c      else if (type.eq.4) then
cc
cc        initialization
cc        ---------------------------------------------------------------    
c         BVar   = var_oper(B2q,ResB,imu)
cc     
cc        find Q-operator by {\cal L}^2 Q = 4 BC Q
cc        ---------------------------------------------------------------
c         do itr=1,TRlen
c            Q2q(itr)= B2q(itr)/H2q(itr)**2 ! Q_ij = B_ij/(e_i + e_j)**2
c         enddo
cc
c#IF(FULLSPACE)
cc        inversion of {\cal L}_0 in residual space 
cc        ---------------------------------------------------------------
cc        initialize using the B-Operator
c         if (isresact.eq.1) then
cc           
cc           initialize using the B-Operator
c            do ist=1,nst(1)+nst(2)
c               do ir=0,nr
c                  do iz=0,nz
c                   ResQ(iz,ir,1,ist,1)=ResB(iz,ir,1,ist,1)
c                   ResQ(iz,ir,2,ist,1)=ResB(iz,ir,2,ist,1)
c                   ResQ(iz,ir,1,ist,2)=ResB(iz,ir,1,ist,2)
c                   ResQ(iz,ir,2,ist,2)=ResB(iz,ir,2,ist,2)
c                  enddo
c               enddo
c            enddo
cc
cc           calculate norm of residual componentes B-Operator
c            do ist=1,nst(1)+nst(2)
c               Bnorm(ist)=
c     &              wfovlp(ResB(0,0,1,ist,1),ResB(0,0,1,ist,1))
c     &             +wfovlp(ResB(0,0,1,ist,2),ResB(0,0,1,ist,2))
c            enddo            
cc
c         else
cc
cc           set residual components to zero
c            call zero_oper(Q2q,ResQ,2)  
c            call zero_oper(HP2q,ResHP,2)           
cc
c         endif
c#ENDIF
cc
cc        ===============================================================
cc        the cranking iteration (find Q-Operator)
cc        ===============================================================
c         do iter=1,max_citer
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           calculate derivatives of resiual part of Q
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
c               call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                              ddRes(0,0,1,1,1),
c     &                              (-1)**imu,+imu,nst(1)+nst(2))
c               call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                              ddRes(0,0,1,1,2),
c     &                              (-1)**imu,-imu,nst(1)+nst(2))
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           calculate action of Liouville-Operator
cc           ------------------------------------------------------------
c            call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
c     &                    HQ2q,ResHQ,imu)
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           calculate derivatives of resiual part of LQ \propto P
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
c               call upd_derive_spinor(ResHQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                              dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                              ddRes(0,0,1,1,1),
c     &                              (-1)**imu,+imu,nst(1)+nst(2))
c               call upd_derive_spinor(ResHQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                              dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                              ddRes(0,0,1,1,2),
c     &                              (-1)**imu,-imu,nst(1)+nst(2))
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           calculate action of Liouville-Operator
cc           ------------------------------------------------------------
c            call comm_HP(H2q,HQ2q,ResHQ,drRes,dzRes,dpRes,ddRes,
c     &                    HP2q,ResHP,imu)
cc
cc           ------------------------------------------------------------ 
cc           calculate Omega^2
cc           ------------------------------------------------------------
c            OmgSq=comm_QP(Q2q,ResQ,HP2q,ResHP,imu)/
c     &            comm_QP(Q2q,ResQ, Q2q,ResQ, imu)
cc
cc           ------------------------------------------------------------ 
cc           perform the gradient step in qp-space
cc           ------------------------------------------------------------
c            do itr=1,TRLen
c               HQ2q(itr) = HP2q(itr)-OmgSq*Q2q(itr)
c               Q2q(itr)  = Q2q(itr) -stp_qp_Omg *
c     &                    (HP2q(itr)-OmgSq*Q2q(itr))/H2q(itr)**2
c            enddo
cc
cc           ------------------------------------------------------------
cc           evaluate errors (only if residual parts are neglected)
cc           ------------------------------------------------------------
c            if (isresact.eq.0) then
c               buffer   = error_Q
c               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
c               derror_Q = error_Q-buffer
c            endif            
cc
c#IF(FULLSPACE)
cc           ------------------------------------------------------------ 
cc           perform the gradient step in residual space
cc           ------------------------------------------------------------
c            if (isresact.eq.1) then
cc
cc              --------------------------------------------------------- 
cc              calculate norm of residual components of [H,Q]-Operator
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  HPnorm(ist)=
c     &                 wfovlp(resHP(0,0,1,ist,1),resHP(0,0,1,ist,1))+
c     &                 wfovlp(resHP(0,0,1,ist,2),resHP(0,0,1,ist,2))
c               enddo
cc
cc              --------------------------------------------------------- 
cc              calculate errors in residual space
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  do ir=0,nr
c                     do iz=0,nz
c                        ResHQ(iz,ir,1,ist,1) = ResHP(iz,ir,1,ist,1)-
c     &                                         ResB(iz,ir,1,ist,1)
c                        ResHQ(iz,ir,2,ist,1) = ResHP(iz,ir,2,ist,1)-
c     &                                         ResB(iz,ir,2,ist,1)
c                        ResHQ(iz,ir,1,ist,2) = ResHP(iz,ir,1,ist,2)-
c     &                                         ResB(iz,ir,1,ist,2)
c                        ResHQ(iz,ir,2,ist,2) = ResHP(iz,ir,2,ist,2)-
c     &                                         ResB(iz,ir,2,ist,2)
c                     enddo
c                  enddo
c               enddo
cc
c               buffer   = error_Q
c               error_Q  = sqrt(var_oper(HQ2q,resHQ,imu))
c               derror_Q = error_Q-buffer
cc
cc              --------------------------------------------------------- 
cc              apply damping
cc              ---------------------------------------------------------
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHQ(0,0,1,1,1),ResHP(0,0,1,1,1),
c     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHP(0,0,1,1,1),ResHQ(0,0,1,1,1),
c     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHQ(0,0,1,1,2),ResHP(0,0,1,1,2),
c     &                          (-1)**imu,imu,1.0D0,nst(1)+nst(2))
c               call damp_spinor(LRDampr,LRDampz,
c     &                          ResHP(0,0,1,1,2),ResHQ(0,0,1,1,2),
c     &                          (-1)**imu,imu,stp_res_Omg,nst(1)+nst(2))
c
cc
cc              --------------------------------------------------------- 
cc              the gradient step
cc              ---------------------------------------------------------
c               do ist=1,nst(1)+nst(2)
c                  do ir=0,nr
c                     do iz=0,nz
c                      ResQ(iz,ir,1,ist,1) = 
c     &                ResQ(iz,ir,1,ist,1) - ResHQ(iz,ir,1,ist,1)
c                      ResQ(iz,ir,2,ist,1) = 
c     &                ResQ(iz,ir,2,ist,1) - ResHQ(iz,ir,2,ist,1)
c                      ResQ(iz,ir,1,ist,2) = 
c     &                ResQ(iz,ir,1,ist,2) - ResHQ(iz,ir,1,ist,2)
c                      ResQ(iz,ir,2,ist,2) = 
c     &                ResQ(iz,ir,2,ist,2) - ResHQ(iz,ir,2,ist,2)
c                     enddo
c                  enddo
c               enddo
cc
cc              --------------------------------------------------------- 
cc              project out sp-contributions
cc              ---------------------------------------------------------
c               call project_qp(Q2q,ResQ,+1,imu,2)
cc
cc           ------------------------------------------------------------ 
cc           END of gradient step in residual space
cc           ------------------------------------------------------------
c            endif
c#ENDIF
cc
cc           ------------------------------------------------------------ 
cc           normalize Q to 2*<QQ> = 1
cc           ------------------------------------------------------------
c            norm=1.0/sqrt(comm_QP(Q2q,ResQ, Q2q,ResQ, imu))
c            do itr=1,TRlen
c               Q2q(itr)=Q2q(itr)*norm
c            enddo
cc
cc
cc           ------------------------------------------------------------ 
cc           Logging and convergency-check
cc           ------------------------------------------------------------
cc     
c            if (mod(iter,lograte).eq.0               .or.
c     &          error_Q.le.max_err_Q*sqrt(abs(QVar)) .or.
c     &          iter.eq.max_citer ) then
cc     
c               write(21,'(a,i5,1(1pg16.4),2(e16.3))') 'ATDHF     : ',
c     &              iter,error_Q,error_Q/sqrt(abs(BVar)),
c     &              derror_Q
cc
c               if (error_Q.le.max_err_Q*sqrt(abs(QVar))) goto 2104
c               if (abs(derror_Q).lt.1.0E-10) goto 2104
cc     
c            endif
cc
c         enddo
cc        ===============================================================
cc        END of cranking iteration (find Q-Operator)
cc        ===============================================================
c 2104    continue    
cc
c#IF(FULLSPACE)
cc        ------------------------------------------------------------ 
cc        calculate derivatives of residual part of Q
cc         ------------------------------------------------------------
c         if (isresact.eq.1) then
c            call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                           dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                           ddRes(0,0,1,1,1),
c     &                           (-1)**imu,+imu,nst(1)+nst(2))
c            call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                            dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                            ddRes(0,0,1,1,2),
c     &                            (-1)**imu,-imu,nst(1)+nst(2))
c         endif
c#ENDIF
cc
cc        ------------------------------------------------------------ 
cc        calculate action of Liouville-Operator to obtain P-Operator
cc        ------------------------------------------------------------
c         call comm_HQ(H2q,Q2q,ResQ,drRes,dzRes,dpRes,ddRes,
c     &                P2q,ResP,imu)
cc
c      endif
cc
c      write(88,'(a)') '--- cranking finished  ---'
cc
cc     ==================================================================
cc     end of the cranking methods
cc     ==================================================================
cc
cc     Normalization:
cc
cc     normalize P to Beta  and Q to P
cc
cc       1) GCM-formulas do not need any explicit normalization.
cc       2) for inglis and sc-cranking normalization is defined by 
cc
cc             [H,Q] = -i 2B P  --> [H, Q/(2B)] = -iP
cc
cc          and restored by Q -> 2BQ
cc
cc       3) CHF-cranking and ATDHF-cranking need explicit normalisation
cc
cc     in case of fullspace calculations (isresact=1) Q,P could be
cc     normalized explicitly in any case.     
cc     
cc
cc     explicit normalization of P and Q to B 
cc     ==================================================================
c      if (FULLSPACE.or.
c     & ((type.eq.0.or.type.eq.4.or.type.eq.5).or.(isresact.eq.1))) then
cc
c         B=-1.0D0   ! just set a flag -- this is not the real mass
cc
c#IF(FULLSPACE)
cc
cc       normalization of P to B       <[B,P]> = i
cc       ----------------------------------------------------------------
c         norm=1.0D0/comm_QP(B2q,resB,P2q,resP,imu)
c         do itr=1,TRlen
c            P2q(itr)= P2q(itr)*norm
c         enddo
cc     
c         do ist=1,nst(1)+nst(2)
c            do ir=0,nr
c               do iz=0,nz
c                  resP(iz,ir,1,ist,1) = resP(iz,ir,1,ist,1)*norm
c                  resP(iz,ir,2,ist,1) = resP(iz,ir,2,ist,1)*norm
c                  resP(iz,ir,1,ist,2) = resP(iz,ir,1,ist,2)*norm
c                  resP(iz,ir,2,ist,2) = resP(iz,ir,2,ist,2)*norm
c               enddo
c            enddo
c         enddo
cc
cc        normalization of Q to P      <[Q,P]> = i
cc        ---------------------------------------------------------------
c         norm = 1.0D0/comm_QP(Q2q,resQ,P2q,resP,imu)
c         do itr=1,TRlen
c            Q2q(itr)= Q2q(itr)*norm
c         enddo
cc     
c         do ist=1,nst(1)+nst(2)
c            do ir=0,nr
c               do iz=0,nz
c                  resQ(iz,ir,1,ist,1) = resQ(iz,ir,1,ist,1)*norm
c                  resQ(iz,ir,2,ist,1) = resQ(iz,ir,2,ist,1)*norm
c                  resQ(iz,ir,1,ist,2) = resQ(iz,ir,1,ist,2)*norm
c                  resQ(iz,ir,2,ist,2) = resQ(iz,ir,2,ist,2)*norm
c               enddo
c            enddo
c         enddo
c#ENDIF
cc
cc
cc     implicit normalization of Q to P:  i[H,Q] = 2BP
cc     ==================================================================
c      else if (type.eq.2.or.type.eq.3) then
cc
c         B    = 0.5/comm_QP(Q2q,resQ,P2q,resP,imu) 
c         do itr=1,TRlen
c            Q2q(itr)= 2.0D0*B*Q2q(itr)
c         enddo
cc
cc     the GCM-mass with fullspace precission
cc     ==================================================================
c      else if (type.eq.1) then
cc
c         
c         
cc
c      endif
cc
c      write(88,'(a)') '--- normalization finished  ---'
cc
cc     Calculation of widths, masses and stiffness:
cc     ==================================================================
cc     calculate <Q^2>,<P^2>,<[Beta,P]>,<[Q,P]>
cc     ------------------------------------------------------------------
c      QVar   = var_oper(Q2q,resQ,imu)
c#IF(FULLSPACE)
c      Nrm0   = comm_QP(B2q,resB,P2q,resP,imu) 
c      Nrm    = comm_QP(Q2q,resQ,P2q,resP,imu)
c      PVar   = var_oper(P2q,resP,imu)
c#ELSE
c      call comm_and_var_QP(BVar,PVar,Nrm0,CPbeta(0),
c     &                     (-1)**imu,imu,0)
c#ENDIF
cc     
cc     calculate LQ = [H,Q] + H_Q  
cc     and        B = 1/2*<[Q,LQ]]>
cc     ------------------------------------------------------------------
cc
c#IF(FULLSPACE)
cc     prepare derivatives of residual part of Q
c      call upd_derive_spinor(ResQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                       dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                       ddRes(0,0,1,1,1),
c     &                       (-1)**imu,+imu,nst(1)+nst(2))
c      call upd_derive_spinor(ResQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                       dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                       ddRes(0,0,1,1,2),
c     &                       (-1)**imu,-imu,nst(1)+nst(2))
c#ENDIF
cc
cc     if mass has not been calculated before
c      if (B.eq.-1.0D0) then
cc           
cc        calculate action of Liouville-Operator
c         call comm_HQ(H2q,Q2q,resQ,drRes,dzRes,dpRes,ddRes,
c     &                HQ2q,resHQ,imu)
c         B = 0.5D0*comm_QP(Q2q,resQ,HQ2q,resHQ,imu)
cc
c      endif
cc     
cc     calculate  LP = [H,P] + H_P  
cc     and        C  = 1/2*<[P,LP]> = -1/2*<[LP,P]>
cc     ------------------------------------------------------------------
c#IF(FULLSPACE)
c      call upd_derive_spinor(ResP(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                       dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                       ddRes(0,0,1,1,1),
c     &                       (-1)**imu,+imu,nst(1)+nst(2))
c      call upd_derive_spinor(ResP(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                       dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                       ddRes(0,0,1,1,2),
c     &                       (-1)**imu,-imu,nst(1)+nst(2))
c#ENDIF
cc     
c      call comm_HP(H2q,P2q,ResP,drRes,dzRes,dpRes,ddRes,
c     &             HP2q,resHP,imu)
cc
c      C = -0.5D0*comm_QP(HP2q,ResHP,P2q,ResP,imu)
cc
cc     calculate L^2Q    = [H,LQ] + H_{LQ}
cc     and       Omega^2 = <[P,L^2Q]>
cc     ------------------------------------------------------------------
c#IF(FULLSPACE)
c      call upd_derive_spinor(ResHQ(0,0,1,1,1),drRes(0,0,1,1,1),
c     &                       dzRes(0,0,1,1,1),dpRes(0,0,1,1,1),
c     &                       ddRes(0,0,1,1,1),
c     &                       (-1)**imu,+imu,nst(1)+nst(2))
c      call upd_derive_spinor(ResHQ(0,0,1,1,2),drRes(0,0,1,1,2),
c     &                       dzRes(0,0,1,1,2),dpRes(0,0,1,1,2),
c     &                       ddRes(0,0,1,1,2),
c     &                       (-1)**imu,-imu,nst(1)+nst(2))
c#ENDIF
cc     
c      call comm_HP(H2q,HQ2q,ResHQ,drRes,dzRes,dpRes,ddRes,
c     &             HP2q,resHP,imu)
cc
c      OmgSq = comm_QP(HP2q,ResHP,P2q,ResP,imu)
cc
cc     calculate particle number mass B_N
cc     ------------------------------------------------------------------
c      call comm_NQ(Q2q,ResQ,HQ2q,ResHQ)
c      BN = 0.5D0*comm_QP(Q2q,ResQ,HQ2q,ResHQ,imu)
cc
c      write(88,'(a)') '--- evaluation finished  ---'
cc
cc     Logging
cc     ------------------------------------------------------------------
c      if (type.eq.0) then
c         write(21,'(a,8(1pe16.4))') 'CHF_SC_VAL: ',
c     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
c      elseif (type.eq.1) then
c         write(21,'(a,8(1pe16.4))') 'GCM_VAL   : ',
c     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
c      elseif (type.eq.2) then
c         write(21,'(a,8(1pe16.4))') 'INGLIS_VAL: ',
c     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
c      elseif (type.eq.3) then
c         write(21,'(a,8(1pe16.4))') 'SELFC_VAL : ',
c     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
c      elseif (type.eq.4) then
c         write(21,'(a,8(1pe16.4))') 'ATDHF_VAL : ',
c     &        B,C,PVar,QVar,OmgSq,BN,Nrm,Nrm0
c      endif
cc
cc     ------------------------------------------------------------------
cc
c      return
c      end


c=======================================================================
c
c     upd_path()
c     
c       updates the collective path which is stored in trpsi. 
c       The subroutine manages the storage of the last three sets 
c       of wavefunctions and occupation amplitudes from the static
c       iteration and also some relevant observables. 
c
c=======================================================================
      subroutine upd_path()
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'observ.inc' ! beta2,beta3,atnum,rms
      include 'constr.inc' ! lqcon
c
c     ------------------------------------------------------------------
      integer           ist,ir,iz,i
      double precision  equasi,deltai,edifi
c
c
c     shift wavefunctions in path storage CPpsi (-1,0,+1) -> (0,+1,+1)
c     and copy current static wavefunction to last path storage index 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               CPpsi(iz,ir,1,ist,-1) = CPpsi(iz,ir,1,ist, 0)
               CPpsi(iz,ir,2,ist,-1) = CPpsi(iz,ir,2,ist, 0)
               CPpsi(iz,ir,1,ist, 0) = CPpsi(iz,ir,1,ist,+1)
               CPpsi(iz,ir,2,ist, 0) = CPpsi(iz,ir,2,ist,+1)
               CPpsi(iz,ir,1,ist,+1) = psi(iz,ir,1,ist)
               CPpsi(iz,ir,2,ist,+1) = psi(iz,ir,2,ist)
            enddo
         enddo
      enddo
c
c     shift mean-fields in path storage   (-1,0,+1) -> (0,+1,+1)
c     and copy current mean-fields to last path storage index 
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            CPrho(iz,ir,1,-1)     = CPrho(iz,ir,1, 0)
            CPrho(iz,ir,2,-1)     = CPrho(iz,ir,2, 0)
            CPrho(iz,ir,1, 0)     = CPrho(iz,ir,1,+1)
            CPrho(iz,ir,2, 0)     = CPrho(iz,ir,2,+1)
            CPrho(iz,ir,1,+1)     = rho(iz,ir,1)
            CPrho(iz,ir,2,+1)     = rho(iz,ir,2)

            CPkin(iz,ir,1,-1)     = CPkin(iz,ir,1, 0)
            CPkin(iz,ir,2,-1)     = CPkin(iz,ir,2, 0)
            CPkin(iz,ir,1, 0)     = CPkin(iz,ir,1,+1)
            CPkin(iz,ir,2, 0)     = CPkin(iz,ir,2,+1)
            CPkin(iz,ir,1,+1)     = tau(iz,ir,1)
            CPkin(iz,ir,2,+1)     = tau(iz,ir,2)
            
            CPsocr(iz,ir,1,-1)    = CPsocr(iz,ir,1, 0)
            CPsocr(iz,ir,2,-1)    = CPsocr(iz,ir,2, 0)
            CPsocr(iz,ir,1, 0)    = CPsocr(iz,ir,1,+1)
            CPsocr(iz,ir,2, 0)    = CPsocr(iz,ir,2,+1)
            CPsocr(iz,ir,1,+1)    = rholsr(iz,ir,1)
            CPsocr(iz,ir,2,+1)    = rholsr(iz,ir,2)

            CPsocz(iz,ir,1,-1)    = CPsocz(iz,ir,1, 0)
            CPsocz(iz,ir,2,-1)    = CPsocz(iz,ir,2, 0)
            CPsocz(iz,ir,1, 0)    = CPsocz(iz,ir,1,+1)
            CPsocz(iz,ir,2, 0)    = CPsocz(iz,ir,2,+1)
            CPsocz(iz,ir,1,+1)    = rholsz(iz,ir,1)
            CPsocz(iz,ir,2,+1)    = rholsz(iz,ir,2)

            CPlocpot(iz,ir,1,-1) = CPlocpot(iz,ir,1, 0)
            CPlocpot(iz,ir,2,-1) = CPlocpot(iz,ir,2, 0)
            CPlocpot(iz,ir,1, 0) = CPlocpot(iz,ir,1,+1)
            CPlocpot(iz,ir,2, 0) = CPlocpot(iz,ir,2,+1)
            CPlocpot(iz,ir,1,+1) = pot(iz,ir,1)
            CPlocpot(iz,ir,2,+1) = pot(iz,ir,2)

            CPeffmas(iz,ir,1,-1) = CPeffmas(iz,ir,1, 0)
            CPeffmas(iz,ir,2,-1) = CPeffmas(iz,ir,2, 0)
            CPeffmas(iz,ir,1, 0) = CPeffmas(iz,ir,1,+1)
            CPeffmas(iz,ir,2, 0) = CPeffmas(iz,ir,2,+1)
            CPeffmas(iz,ir,1,+1) = bmass(iz,ir,1)
            CPeffmas(iz,ir,2,+1) = bmass(iz,ir,2)
            
            CPsocpotr(iz,ir,1,-1) = CPsocpotr(iz,ir,1, 0)
            CPsocpotr(iz,ir,2,-1) = CPsocpotr(iz,ir,2, 0)
            CPsocpotr(iz,ir,1, 0) = CPsocpotr(iz,ir,1,+1)
            CPsocpotr(iz,ir,2, 0) = CPsocpotr(iz,ir,2,+1)
            CPsocpotr(iz,ir,1,+1) = potls(iz,ir,1,1)
            CPsocpotr(iz,ir,2,+1) = potls(iz,ir,2,1)

            CPsocpotz(iz,ir,1,-1) = CPsocpotz(iz,ir,1, 0)
            CPsocpotz(iz,ir,2,-1) = CPsocpotz(iz,ir,2, 0)
            CPsocpotz(iz,ir,1, 0) = CPsocpotz(iz,ir,1,+1)
            CPsocpotz(iz,ir,2, 0) = CPsocpotz(iz,ir,2,+1)
            CPsocpotz(iz,ir,1,+1) = potls(iz,ir,1,2)
            CPsocpotz(iz,ir,2,+1) = potls(iz,ir,2,2)

            CPcoupot(iz,ir,-1)    = CPcoupot(iz,ir, 0)
            CPcoupot(iz,ir, 0)    = CPcoupot(iz,ir,+1)
            CPcoupot(iz,ir,+1)    = coulpt(iz,ir)

            CPconstr(iz,ir,-1)    = CPconstr(iz,ir, 0)
            CPconstr(iz,ir, 0)    = CPconstr(iz,ir,+1)
            CPconstr(iz,ir,+1)    = qfac(iz,ir)
         enddo
      enddo
c      
c     get sp-observables from stationary iteration and shift path
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
         CPvsqf(ist,-1)   = CPvsqf(ist, 0)
         CPvsqf(ist, 0)   = CPvsqf(ist,+1)
         CPvsqf(ist,+1)   = vsqf(ist)

         CPdeltaf(ist,-1) = CPdeltaf(ist, 0)
         CPdeltaf(ist, 0) = CPdeltaf(ist,+1)
         CPdeltaf(ist,+1) = deltaf(ist)
         
         CPwtpair(ist,-1) = CPwtpair(ist, 0)
         CPwtpair(ist, 0) = CPwtpair(ist,+1)
         CPwtpair(ist,+1) = wtpair(ist)

         CPspenrg(ist,-1) = CPspenrg(ist, 0)
         CPspenrg(ist, 0) = CPspenrg(ist,+1)
         CPspenrg(ist,+1) = spenrg(ist)

         CPqpenrg(ist,-1) = CPqpenrg(ist, 0)
         CPqpenrg(ist, 0) = CPqpenrg(ist,+1)

         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-fer(icharg(ist))
         equasi = sqrt(edifi*edifi+deltai*deltai)
         
         if (ipair.ne.-1) then
            CPqpenrg(ist,+1) = equasi-fer(icharg(ist))
         else
            CPqpenrg(ist,+1) = -spenrg(ist)
         endif

      enddo
c      
c     get global observables from stationary iteration and shift path
c     ------------------------------------------------------------------
      CPenergy(-1)   = CPenergy( 0)
      CPenergy( 0)   = CPenergy(+1)
      CPenergy(+1)   = energy 

      CPenergyC(-1)   = CPenergyC( 0)
      CPenergyC( 0)   = CPenergyC(+1)
      CPenergyC(+1)   = qlamb*beta2

      CPlambda(-1)   = CPlambda( 0)
      CPlambda( 0)   = CPlambda(+1)
      CPlambda(+1)   = qlamb

      CPzcm(-1)      = CPzcm( 0)
      CPzcm( 0)      = CPzcm(+1)
      CPzcm(+1)      = zcm

      CPradRMS(-1)   = CPradRMS( 0)      ! r      (RMS radius)
      CPradRMS( 0)   = CPradRMS(+1)
      CPradRMS(+1)   = rchrg2

      CPradDIF(-1)   = CPradDIF( 0)      ! R      (diff. radius)
      CPradDIF( 0)   = CPradDIF(+1)
      CPradDIF(+1)   = rdmsc

      CPradSRF(-1)   = CPradSRF( 0)      ! sigma  (surface. thickness)
      CPradSRF( 0)   = CPradSRF(+1)
      CPradSRF(+1)   = surfc

      CPfer(1,-1)    = CPfer(1, 0)
      CPfer(1, 0)    = CPfer(1,+1)
      CPfer(1,+1)    = fer(1)

      CPfer(2,-1)    = CPfer(2, 0)
      CPfer(2, 0)    = CPfer(2,+1)
      CPfer(2,+1)    = fer(2)

      CPpairfc(1,-1) = CPpairfc(1, 0) 
      CPpairfc(1, 0) = CPpairfc(1, 1) 
      CPpairfc(1, 1) = pairfc(1)

      CPpairfc(2,-1) = CPpairfc(2, 0) 
      CPpairfc(2, 0) = CPpairfc(2, 1) 
      CPpairfc(2, 1) = pairfc(2)

      CPatnum(-1) = CPatnum( 0)
      CPatnum( 0) = CPatnum(+1)
      CPatnum(+1) = atnum
      
      if (lqcon.eq.2) then 
         CPbeta(-1) = CPbeta( 0)
         CPbeta( 0) = CPbeta(+1)
         CPbeta(+1) = beta2
      elseif (lqcon.eq.3) then
         CPbeta(-1) = CPbeta( 0)
         CPbeta( 0) = CPbeta(+1)
         CPbeta(+1) = beta3
      endif
c      
c     set current length of collective path
c     ------------------------------------------------------------------
      CPlen = CPlen + 1
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,3f8.4)')         
     &  'UPD_PATH()    : collective path updated --'//
     &  ' beta= ',(CPbeta(i),i=-1,1)
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     upd_derivatives()
c     
c       updates the derivatives of the sp-wavefunctions and Mean-Fields
c       of the central element of the collective path.
c
c=======================================================================
      subroutine upd_derivatives()
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'  !h2meff
      include 'linresp.inc'
      include 'wavefr.inc'
      include 'levels.inc'
      include 'observ.inc' ! beta2,beta3,atnum,rms
      include 'constr.inc' ! lqcon

c      double precision HPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision energy,var
      integer          ist
c
c     ------------------------------------------------------------------

c
c     calculate derivatives of single particle wavefunctions
c     ------------------------------------------------------------------
      call upd_derive_spinor(CPpsi(0,0,1,1,0),
     &                       CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &                       CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),
     &                       +1,0,nst(1)+nst(2))
c
c     derivatives requiere effmass to vanish on bound of box
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            CPeffmas(iz,ir,1,0) = CPeffmas(iz,ir,1,0) - h2meff(1)
            CPeffmas(iz,ir,2,0) = CPeffmas(iz,ir,2,0) - h2meff(2)
         enddo
      enddo
      
      call nablaz_scalar(CPeffmas(0,0,1,0),CPdzeffmas(0,0,1),1,1,one)
      call nablar_scalar(CPeffmas(0,0,1,0),CPdreffmas(0,0,1),0,1,one) 
      call nablaz_scalar(CPeffmas(0,0,2,0),CPdzeffmas(0,0,2),1,1,one)
      call nablar_scalar(CPeffmas(0,0,2,0),CPdreffmas(0,0,2),0,1,one) 

      do ir=0,nr
         do iz=0,nz
            CPeffmas(iz,ir,1,0) = CPeffmas(iz,ir,1,0) + h2meff(1)
            CPeffmas(iz,ir,2,0) = CPeffmas(iz,ir,2,0) + h2meff(2)
         enddo
      enddo
c
c     test action of Mean-Field-Hamiltonian on states
c     ------------------------------------------------------------------
c      if (1.eq.0) then
c      call act_mfh_even(CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),
c     &                  CPdzpsi(0,0,1,1),CPdppsi(0,0,1,1),
c     &                  CPddpsi(0,0,1,1),HPsi(0,0,1,1),
c     &                  CPlocpot(0,0,1,0),                    1,
c     &                  CPeffmas(0,0,1,0),CPdreffmas(0,0,1),
c     &                  CPdzeffmas(0,0,1),CPdzeffmas(0,0,1),  1,
c     &                  CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
c     &                  CPsocpotr(0,0,1,0),                   1,
c     &                  CPconstr(0,0,0),CPlambda(0),          1,
c     &                  CPspenrg(1,0),                        0,
c     &                  0,1.0D0,nst(1)+nst(2))
cc
c      do ist=1,nst(1)+nst(2)
c         energy = wfovlp(CPpsi(0,0,1,ist,0),HPsi(0,0,1,ist))
c         var    = wfovlp(HPsi(0,0,1,ist),HPsi(0,0,1,ist))-energy**2
c         write(88,'(a,3(1pg16.4),3i3)') 'PKtest:',energy,var, 
c     &        CPspenrg(ist,0),mpsi(ist),nparz(1,ist),icharg(ist)
c      enddo
c      endif
c
c     write derivatives of single particle wavefunctions
c     ------------------------------------------------------------------
c      do ist=1,nst(1)+nst(2)
c         call wrt_2D_spin(CPpsi(0,0,1,ist,0) ,'psi_______',1,'tst')
c         call wrt_2D_spin(CPdrpsi(0,0,1,ist) ,'drpsi_____',0,'tst')
c         call wrt_2D_spin(CPdzpsi(0,0,1,ist) ,'dzpsi_____',0,'tst')
c         call wrt_2D_spin(CPdppsi(0,0,1,ist) ,'dppsi_____',0,'tst')
c         call wrt_2D_spin(CPddpsi(0,0,1,ist) ,'ddpsi_____',0,'tst')
c      enddo         
       
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a)') 'UPD_DERIVES() : derivatives of '//
     &                'wavefunctions updated'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     ini_matrix_P(P2q,PVar,il,imu,iiso)
c
c       initialize the generator of axial deformation (mu=0) or
c       a rotation arround the y-axis (mu=1). In case of rotations
c       the operator is scaled by 1/(sqrt(3)*beta) to be equivalent
c       with the cartesean P_{1} (only at the axial configurations.
c
c       furthermore this function calculates
c
c       2qp-matrixelements    P^{2q}_ji = eta^(-)_ji <j|P|i>
c       variance              <P^2>
c
c=======================================================================
      subroutine ini_matrix_P(P2q,PVar,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision P2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision PVar                   ! variance

      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      integer          trzpar
      double precision pfac,eta,ovlp
      double precision z,r
      integer          iz,ir
      double precision lastv,lastu,nextv,nextu,rcntv,rcntu
      double precision ui,uj,vi,vj
      integer          itr,ist,jst,jstmin,jstmax
c     
      logical          oldcrank
      data             oldcrank/D_OLDCRANK/
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     some factors
c     ------------------------------------------------------------------
      if (imu.eq.0) then
         if (oldcrank) then
            pfac     = 1.0D0/(CPbeta( 0)-CPbeta(-1))
         else
            pfac     = 1.0D0/(CPbeta(+1)-CPbeta(-1))
         endif
      else if (imu.eq.1) then
         if (abs(CPbeta(0)).gt.0.0D0) then
            pfac=1.0D0/(sqrt(3.0D0)*CPbeta(0))
         else
            pfac=0.0D0
         endif
      else
         write(21,'(a,i2)') 
     &            'ERROR: only mu=0 and mu=1 possible -- not ',imu
         stop 'in linresp'
      endif
c
c     set variance and 2q-matrixelements to zero
c     ------------------------------------------------------------------
      PVar = 0.0D0
      do itr2=1,TRntrtot
         P2q(itr2)=0.0D0
      enddo
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c     
c        calculate action of P operator 
c        ===============================================================
c
c        P_0 (generator of translation along the axial path, asym. def.)
c        --------------------------------------------------------------
         if (imu.eq.0.and.oldcrank) then   
c
            do ir=0,nr
               do iz=0,nz
                  LRbuff1(iz,ir,1,1)=pfac*
     &               (CPpsi(iz,ir,1,ist, 0)-CPpsi(iz,ir,1,ist,-1))
                  LRbuff1(iz,ir,2,1)=pfac*
     &               (CPpsi(iz,ir,2,ist, 0)-CPpsi(iz,ir,2,ist,-1))
                  LRbuff1(iz,ir,1,2)=0.0D0
                  LRbuff1(iz,ir,2,2)=0.0D0
               enddo
            enddo
c
c        P_0 (generator of translation along the axial path, sym. def.)
c        --------------------------------------------------------------
         else if (imu.eq.0) then 
c
            do ir=0,nr
               do iz=0,nz
                  LRbuff1(iz,ir,1,1)=pfac*
     &               (CPpsi(iz,ir,1,ist,+1)-CPpsi(iz,ir,1,ist,-1))
                  LRbuff1(iz,ir,2,1)=pfac*
     &               (CPpsi(iz,ir,2,ist,+1)-CPpsi(iz,ir,2,ist,-1))
                  LRbuff1(iz,ir,1,2)=0.0D0
                  LRbuff1(iz,ir,2,2)=0.0D0
               enddo
            enddo 
c
c        P_1 (generator of rotation,  J_y = L_y + S_y)
c        --------------------------------------------------------------
         else if (imu.eq.1) then
c
            do ir=0,nr
               r = (ir+sr)*dr
               do iz=0,nz 
                  z = (iz+sz-icm)*dz - CPzcm(0)
c     
c                 ! upward laddering components of J_y
                  LRbuff1(iz,ir,1,1) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,1,ist) 
     &                                  +z*CPdppsi(iz,ir,1,ist) 
     &                                  +r*CPdzpsi(iz,ir,1,ist)
     &                                  -  CPpsi(iz,ir,2,ist,0))
                  LRbuff1(iz,ir,2,1) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,2,ist) 
     &                                  +z*CPdppsi(iz,ir,2,ist) 
     &                                  +r*CPdzpsi(iz,ir,2,ist))
c
c                 ! downward laddering components of J_y
                  LRbuff1(iz,ir,1,2) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,1,ist) 
     &                                  -z*CPdppsi(iz,ir,1,ist) 
     &                                  +r*CPdzpsi(iz,ir,1,ist))
                  LRbuff1(iz,ir,2,2) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,2,ist) 
     &                                  -z*CPdppsi(iz,ir,2,ist) 
     &                                  +r*CPdzpsi(iz,ir,2,ist)
     &                                  +  CPpsi(iz,ir,1,ist,0))
c
               enddo
            enddo
c
         endif           
c
c        adjust lower boundary in case of no z-symmetry
c        ---------------------------------------------------------------
         if (isym.eq.0) then
            do ir=0,nr
               LRbuff1(0,ir,1,1) = 0.0D0
               LRbuff1(0,ir,2,1) = 0.0D0
               LRbuff1(0,ir,1,2) = 0.0D0
               LRbuff1(0,ir,2,2) = 0.0D0
            enddo
         endif
c
c        evaluate 2qp-matrixelements
c        ------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            itp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)
            eta  = uj*vi - vj*ui            
c
            if (itp.eq.+1) then
               if (msgn.ge.0) then
                  ovlp=wfovlp(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,1))  
                  call wfproj(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,1),
     &                        ovlp)
               else       
                  ovlp=wfovlp(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2))  
                  call wfproj(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2),
     &                        ovlp)
               endif
            else
               ovlp=wfovlptr(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2))  
               call wfprojtr(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2),
     &                       ovlp)
            endif
c
c           diagonal elements describe change of occupation amplitudes
            if (ist.eq.jst.and.itp.eq.+1) then
c  
               if (oldcrank) then
                  lastv   = sqrtb(CPvsqf(ist,-1))
                  lastu   = sqrtb(one-CPvsqf(ist,-1))
                  nextv   = sqrtb(CPvsqf(ist,0))
                  nextu   = sqrtb(one-CPvsqf(ist,0))
                  rcntv   = 0.5D0*(lastv+nextv)
                  rcntu   = 0.5D0*(lastu+nextu)
                  if (rcntv.ge.0.5D0) then
                     ovlp =  -(nextu-lastu)/rcntv*pfac
                  else
                     ovlp =  +(nextv-lastv)/rcntu*pfac
                  endif
               else
                  lastv   = sqrtb(CPvsqf(ist,-1))
                  lastu   = sqrtb(one-CPvsqf(ist,-1))
                  rcntv   = sqrtb(CPvsqf(ist,0)) 
                  rcntu   = sqrtb(one-CPvsqf(ist,0))
                  nextv   = sqrtb(CPvsqf(ist,+1))
                  nextu   = sqrtb(one-CPvsqf(ist,+1))
                  if (rcntv.ge.0.5D0) then
                     ovlp =  -(nextu-lastu)/rcntv*pfac
                  else
                     ovlp =  +(nextv-lastv)/rcntu*pfac
                  endif
               endif
               P2q(itr2) =               ovlp
               PVar      = PVar + (      ovlp)**2
c
c           straight forward for off-diagonal elements
            else
c
               P2q(itr2) =         eta * ovlp
               PVar      = PVar + (eta * ovlp)**2
c
            endif
c
         enddo
c
c        accumulate contributions from outside qp-space
c        ---------------------------------------------------------------
         PVar = PVar + 2.0D0*CPvsqf(ist,0)*(
     &                       wfovlp(LRbuff1(0,0,1,1),LRbuff1(0,0,1,1))
     &                     + wfovlp(LRbuff1(0,0,1,2),LRbuff1(0,0,1,2)) )
c
      enddo
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'INI_MATRIX_P(): P-operator initialized'//
     &  '  (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     crank_test_resp(B2q,Q2q,P2q,H2q,il,imu)
c
c       checks the responses
c
c=======================================================================
      subroutine crank_test_resp(B2q,Q2q,P2q,H2q,mass,stiff,qvar,pvar,
     &                       il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass                   ! collective mass B=1/(2M)
      double precision qvar                   ! <Q^2>-<Q>^2
      double precision pvar                   ! <P^2>
      double precision stiff                  ! stiffness C
c
      integer          itr
      integer          iter
      double precision norm,norm0,OmgSq
      double precision Hresp2q(1:MTDTOT)
      double precision HQ2q(1:MTDTOT)
      double precision D2q(1:MTDTOT)
      double precision maxresp
c
      integer          maxiter
      double precision stepQ
      double precision m1,m2,m3,m4,m5,m6
      double precision s1,s2,s3,s4,s5,s6
c
      maxiter = 200
      stepQ   = 0.8D0
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp A'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp B'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',LAMBDAMAX,')'
         stop 'in linresp C'
      endif
c
c     initialize with Q-operator in inglis approximation
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
         P2q(itr) = B2q(itr)
         Q2q(itr) = B2q(itr)
      enddo
c
c     check the Q-response
c     ------------------------------------------------------------------
      call upd_response_Qtest(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu,
     &     1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      m1  =0.0D0
      do itr=1,TRntrtot
         m1     = m1  + 2.0D0*Q2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Qtest(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu,
     &     0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      m2  =0.0D0
      do itr=1,TRntrtot
         m2     = m2  + 2.0D0*Q2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Qtest(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu,
     &     0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      m3  =0.0D0
      do itr=1,TRntrtot
         m3     = m3  + 2.0D0*Q2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Qtest(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu,
     &     0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      m4  =0.0D0
      do itr=1,TRntrtot
         m4     = m4  + 2.0D0*Q2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Qtest(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu,
     &     0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      m5  =0.0D0
      do itr=1,TRntrtot
         m5     = m5  + 2.0D0*Q2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Qtest(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      m6  =0.0D0
      do itr=1,TRntrtot
         m6     = m6  + 2.0D0*Q2q(itr)*Hresp2q(itr) 
      enddo
c
c      write(89,'(i2,7g16.4)') imu,CPbeta(0),m1,m2,m3,m4,m5,m6 
c      return
c
c     check the P-response
c     ------------------------------------------------------------------
      call upd_response_Ptest(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu,
     &     1.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      s1  =0.0D0
      do itr=1,TRntrtot
         s1     = s1  + 2.0D0*P2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Ptest(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu,
     &     0.0D0,1.0D0,0.0D0,0.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      s2  =0.0D0
      do itr=1,TRntrtot
         s2     = s2  + 2.0D0*P2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Ptest(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu,
     &     0.0D0,0.0D0,1.0D0,0.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      s3  =0.0D0
      do itr=1,TRntrtot
         s3     = s3  + 2.0D0*P2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Ptest(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu,
     &     0.0D0,0.0D0,0.0D0,1.0D0,0.0D0,0.0D0)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      s4  =0.0D0
      do itr=1,TRntrtot
         s4     = s4  + 2.0D0*P2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Ptest(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu,
     &     0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,0.0D0)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      s5  =0.0D0
      do itr=1,TRntrtot
         s5     = s5  + 2.0D0*P2q(itr)*Hresp2q(itr) 
      enddo
c
      call upd_response_Ptest(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu,
     &     0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      s6  =0.0D0
      do itr=1,TRntrtot
         s6     = s6  + 2.0D0*P2q(itr)*Hresp2q(itr) 
      enddo
c
      write(89,'(i2,7g16.4)') imu,CPbeta(0),s1,s2,s3,s4,s5,s6 
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     crank_inglis(B2q,Q2q,P2q,H2q,il,imu)
c
c       cranking in inglis approximation
c 
c               i[H,(Q/2B)] = P            1/(2B) = <[(Q/2B),P]>
c
c=======================================================================
      subroutine crank_inglis(B2q,Q2q,P2q,H2q,mass,stiff,qvar,pvar,
     &                        il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass                   ! collective mass B=1/(2M)
      double precision qvar                   ! <Q^2>-<Q>^2
      double precision pvar                   ! <P^2>
      double precision stiff                  ! stiffness C
c
      integer          itr
      double precision norm,norm0,OmgSq
      double precision HP2q(1:MTDTOT)
      double precision Hresp2q(1:MTDTOT)
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c        
c     evaluate the P-response
c     ------------------------------------------------------------------
      call upd_response_even(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      do itr=1,TRntrtot
         HP2q(itr) = H2q(itr)*P2q(itr) + Hresp2q(itr)
      enddo
c
c     set variables to zero
c     ------------------------------------------------------------------
      mass  = 0.0D0
      norm  = 0.0D0
      qvar  = 0.0D0
      pvar  = 0.0D0
      stiff = 0.0D0
c        
c     loop over all transitions
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
c
         Q2q(itr) =               P2q(itr)/H2q(itr)       ! 1/2B    Q
         mass     = mass  + 2.0D0*P2q(itr)**2/H2q(itr)    ! 1/2B
         qvar     = qvar  +       P2q(itr)**2/H2q(itr)**2 ! 1/(4BB) <Q^2>
         pvar     = pvar  +       P2q(itr)**2             !         <P^2>
         stiff    = stiff + 2.0D0*P2q(itr)*HP2q(itr  )    ! 2*C 
c
      enddo
c
c     normalize 
c     ------------------------------------------------------------------
      norm  = 1.0D0/mass
      mass  = 0.5D0/mass
      stiff = 0.5D0*stiff
c
      qvar  =  norm**2*qvar
      pvar  =          pvar
c
      OmgSq = mass*stiff
c
      do itr=1,TRntrtot
         Q2q(itr) = Q2q(itr)*norm
      enddo
c
c     check the normalisation
      norm =0.0D0
      norm0=0.0D0
      do itr=1,TRntrtot
         norm0    = norm0 + 2.0D0*B2q(itr)*P2q(itr)
         norm     = norm  + 2.0D0*Q2q(itr)*P2q(itr)
      enddo    
c
c     logging
      write(21,'(a,7(1pe16.4))') 'INGLIS_VAL: ',
     &        mass,stiff,PVar,QVar,OmgSq,norm,norm0
c
c     ------------------------------------------------------------------
c
      return
      end

c=======================================================================
c
c     crank_GCM(B2q,Q2q,P2q,H2q,il,imu)
c
c       cranking in static GCM approximation
c 
c               Q = P/lambda
c
c=======================================================================
      subroutine crank_GCM(B2q,Q2q,P2q,H2q,mass,stiff,qvar,pvar,
     &                        il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass                   ! collective mass B=1/(2M)
      double precision qvar                   ! <Q^2>-<Q>^2
      double precision pvar                   ! <P^2>
      double precision stiff                  ! stiffness C
c
      integer          itr
      double precision norm,norm0,OmgSq
      double precision Hresp2q(1:MTDTOT)
      double precision HQ2q(1:MTDTOT)
      double precision HP2q(1:MTDTOT)

c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c        
c     calculate the Q-Operator in static GCM approximation
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
         Q2q(itr) = P2q(itr)/(2.0D0*PVar)
      enddo
c        
c     evaluate the Q-response
c     ------------------------------------------------------------------
      call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      do itr=1,TRntrtot
         HQ2q(itr) = H2q(itr)*Q2q(itr) + Hresp2q(itr)
      enddo
c        
c     loop over all transitions
c     ------------------------------------------------------------------
      mass=0.0D0
      qvar=0.0D0
      do itr=1,TRntrtot
         mass     = mass  + 2.0D0*HQ2q(itr)*Q2q(itr) ! 2B
         qvar     = qvar  +        Q2q(itr)*Q2q(itr)
      enddo
      mass = 0.5D0*mass
c        
c     evaluate the P-response
c     ------------------------------------------------------------------
      call upd_response_even(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
      do itr=1,TRntrtot
         HP2q(itr) = H2q(itr)*P2q(itr) + Hresp2q(itr)
      enddo
c        
c     loop over all transitions
c     ------------------------------------------------------------------
      stiff=0.0D0
      pvar =0.0D0
      norm =0.0D0
      norm0=0.0D0
      do itr=1,TRntrtot
         pvar      = pvar   +        P2q(itr)*P2q(itr)
         norm      = norm   + 2.0D0* Q2q(itr)*P2q(itr)
         norm0     = norm0  + 2.0D0* B2q(itr)*P2q(itr)
         stiff     = stiff  + 2.0D0*HP2q(itr)*P2q(itr) ! 2C
      enddo
      stiff = 0.5D0*stiff
c
      OmgSq = mass*stiff
c
      do itr=1,TRntrtot
         Q2q(itr) = Q2q(itr)*norm
      enddo
c
c     logging
      write(21,'(a,7(1pe16.4))') 'GCM_VAL   : ',
     &        mass,stiff,PVar,QVar,OmgSq,norm,norm0
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     crank_scaling(B2q,Q2q,P2q,H2q,il,imu)
c
c       cranking in scaling approximation
c 
c          [H,[H,Q]] = 4 B C Beta      B=1/2*<[Q,[H,Q]]>
c                                      C=1/2*<[P,[H,P]]>
c
c=======================================================================
      subroutine crank_scaling(B2q,Q2q,P2q,H2q,mass,stiff,
     &                         qvar,pvar,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass,stiff
      double precision qvar,pvar
c
      integer          ist,jst
      integer          trzpar
      integer          itr
      integer          iso,iiso
      double precision norm,norm0,normP,normQ,OmgSq
c
c     ------------------------------------------------------------------
c
c     set variables to zero
c     ------------------------------------------------------------------
      mass  = 0.0D0
      norm  = 0.0D0
      norm0 = 0.0D0
      qvar  = 0.0D0
      pvar  = 0.0D0
c        
c     loop over all transitions
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
c
         P2q(itr) =               B2q(itr)/H2q(itr)       ! = 1/(2C)   P
         Q2q(itr) =               B2q(itr)/H2q(itr)**2    ! = 1/(4BC) Q
         norm0    = norm0 + 2.0D0*B2q(itr)**2/H2q(itr)    ! = 1/(2C) 
         pvar     = pvar  +       B2q(itr)**2/H2q(itr)**2 ! = 1/(4CC) <P^2>
         mass     = mass  + 2.0D0*B2q(itr)**2/H2q(itr)**3 ! = 1/(8BCC)
         qvar     = qvar  +       B2q(itr)**2/H2q(itr)**4 ! = 1/(16BBCC) Q
c
      enddo
c
c     normalize 
c     ------------------------------------------------------------------
      normQ = norm0/mass   ! factor to restore norm of Q operator
      normP = 1.0D0/norm0  ! factor to restore norm of P operator
c
      qvar  = normQ**2       * qvar
      pvar  = normP**2       * pvar
      mass  = 0.5D0*norm0**2 / mass
      stiff = 0.5D0*normP
      OmgSq = mass*stiff  
c
      do itr=1,TRntrtot
         P2q(itr) = P2q(itr)*normP
         Q2q(itr) = Q2q(itr)*normQ
      enddo
c
c     check the normalisation
      norm =0.0D0
      norm0=0.0D0
      do itr=1,TRntrtot
         norm0    = norm0 + 2.0D0*B2q(itr)*P2q(itr)
         norm     = norm  + 2.0D0*Q2q(itr)*P2q(itr)
      enddo      
c
c     logging
      write(21,'(a,7(1pe16.4))') 'SCALE_VAL : ',
     &        mass,stiff,PVar,QVar,OmgSq,norm,norm0
c
c     ------------------------------------------------------------------
c
      return
      end

c=======================================================================
c
c     crank_selfc(B2q,Q2q,P2q,H2q,il,imu)
c
c       ATDHF cranking with a given P-operator
c 
c               i[H,(Q/2B)] = P            1/(2B) = <[(Q/2B),P]>
c
c=======================================================================
      subroutine crank_selfc(B2q,Q2q,P2q,H2q,mass,stiff,qvar,pvar,
     &                       il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass                   ! collective mass B=1/(2M)
      double precision qvar                   ! <Q^2>-<Q>^2
      double precision pvar                   ! <P^2>
      double precision stiff                  ! stiffness C
c
      integer          itr
      integer          iter
      double precision norm,norm0,OmgSq
      double precision Hresp2q(1:MTDTOT)
      double precision HQ2q(1:MTDTOT)
      double precision HP2q(1:MTDTOT)
      double precision D2q(1:MTDTOT)
      double precision maxresp
c
      integer          maxiter
      double precision stepQ
c
      maxiter = 200
      stepQ   = 0.8D0
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      WRITE(21,'(a,i2)') 'CRANK_SELFC with IMU=',imu
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c        
c     evaluate the P-response
c     ------------------------------------------------------------------
      call upd_response_even(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
      do itr=1,TRntrtot
         HP2q(itr) = H2q(itr)*P2q(itr) + Hresp2q(itr)
      enddo
c
c     initialize with Q-operator in inglis approximation
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
         Q2q(itr) = P2q(itr)/H2q(itr)
      enddo
c
c     solve i[H,(Q/2B)] + h_{Q,resp} = P by iteration
c     ------------------------------------------------------------------
      do iter=1,maxiter
c
         call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
         call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c        evaluate the Liouville operator i[H,Q] + h_{Q,resp}
         do itr=1,TRntrtot
            HQ2q(itr) = H2q(itr)*Q2q(itr) + Hresp2q(itr)
         enddo
c
c        perform the gradient step
         do itr=1,TRntrtot
            D2q(itr) = HQ2q(itr) - P2q(itr)
            Q2q(itr) = Q2q(itr)  - stepQ*(HQ2q(itr)-P2q(itr))/(H2q(itr))
         enddo
c
c        evaluate mass and errors
         mass  =0.0D0
         qvar  =0.0D0
         buff1 =0.0D0
         buff2 =0.0D0
         buff3 =0.0D0
         do itr=1,TRntrtot
            mass     = mass  + 2.0D0*Q2q(itr)* P2q(itr) ! 1/2B
            qvar     = qvar  +       Q2q(itr)**2        ! 1/(4BB) <Q^2>
            buff1    = buff1 +      HQ2q(itr)* P2q(itr)
            buff2    = buff2 +      HQ2q(itr)*HQ2q(itr)
            buff3    = buff3 +       P2q(itr)* P2q(itr)
         enddo
         errP=1.0D0-buff1**2/(buff2*buff3)
c
c        logging
         write(22,'(a,2i4,3(1pg12.4))')
     &      'SELFC: imu,iter,mass,qvar,errP:',imu,iter,mass,qvar,errP
         CALL FLUSH(22)
c
c        check for convergency
         if(abs(errP).lt.1D-8) goto 4000
c         
      enddo
 4000 continue
c
c     set variables to zero
c     ------------------------------------------------------------------
      mass  = 0.0D0
      norm  = 0.0D0
      qvar  = 0.0D0
      pvar  = 0.0D0
      stiff = 0.0D0
c        
c     loop over all transitions
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
c
         mass     = mass  + 2.0D0*Q2q(itr)*P2q(itr)       ! 1/2B
         qvar     = qvar  +       Q2q(itr)**2             ! 1/(4BB) <Q^2>
         pvar     = pvar  +       P2q(itr)**2             !         <P^2>
         stiff    = stiff + 2.0D0*P2q(itr)*HP2q(itr)      !  2C 
c
      enddo
c
c     normalize 
c     ------------------------------------------------------------------
      norm  = 1.0D0/mass
      mass  = 0.5D0/mass
      stiff = 0.5D0*stiff
c
      qvar  =  norm**2*qvar
      pvar  =          pvar
c
      OmgSq = mass*stiff
c
      do itr=1,TRntrtot
         Q2q(itr) = Q2q(itr)*norm
      enddo
c
c     check the normalisation
      norm =0.0D0
      norm0=0.0D0
      do itr=1,TRntrtot
         norm0    = norm0 + 2.0D0*B2q(itr)*P2q(itr)
         norm     = norm  + 2.0D0*Q2q(itr)*P2q(itr)
      enddo    
c
c     logging
      write(21,'(a,i3,7(1pe16.4))') 'SELFC_VAL : ',
     &        imu,mass,stiff,PVar,QVar,OmgSq,norm,norm0
      CALL FLUSH(21)
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     crank_CHF(B2q,Q2q,P2q,H2q,il,imu)
c
c       CHF + ATDHF-cranking for a given Beta-operator
c 
c               [H,[H,(Q/4Omega^2)]] = Beta
c
c=======================================================================
      subroutine crank_CHF(B2q,Q2q,P2q,H2q,mass,stiff,qvar,pvar,
     &                       il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass                   ! collective mass B=1/(2M)
      double precision qvar                   ! <Q^2>-<Q>^2
      double precision pvar                   ! <P^2>
      double precision stiff                  ! stiffness C
c
      integer          itr
      integer          iter
      double precision norm,norm0,OmgSq,normP,normQ
      double precision Hresp2q(1:MTDTOT)
      double precision HQ2q(1:MTDTOT)
      double precision HP2q(1:MTDTOT)
      double precision D2q(1:MTDTOT)
c
      integer          maxiter
      double precision stepQ
c
      maxiter = 400
      stepQ   = 0.8D0
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     initialize with Q-operator with saved matrixelements
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
         Q2q(itr) = Q2qCHF(itr,imu)
      enddo
c
c     solve i[H,[H,Q] + h_{Q,resp}] + h_{P,resp} = Beta by iteration
c     ------------------------------------------------------------------
      do iter=1,maxiter
c
         call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
         call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c        evaluate the Liouville operator i[H,Q] + h_{Q,resp}
         do itr=1,TRntrtot
            HQ2q(itr) = H2q(itr)*Q2q(itr) + Hresp2q(itr)
         enddo
c
         call upd_response_even(HQ2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
         call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c        evaluate the Liouville operator -i[H,P] + h_{P,resp}
         do itr=1,TRntrtot
            HP2q(itr) = H2q(itr)*HQ2q(itr) + Hresp2q(itr)
         enddo
c
c        perform the gradient step
         do itr=1,TRntrtot
            D2q(itr) = HP2q(itr) - B2q(itr)
            Q2q(itr) = Q2q(itr)  - stepQ*(HP2q(itr)-B2q(itr))/
     &                                   (H2q(itr)**2)
         enddo
c
c        evaluate mass and errors
         norm  =0.0D0
         norm0 =0.0D0
         qvar  =0.0D0
         pvar  =0.0D0
         buff1 =0.0D0
         buff2 =0.0D0
         buff3 =0.0D0
         do itr=1,TRntrtot
            norm     = norm  + 2.0D0*B2q(itr)*HQ2q(itr) ! 1/(2C)
            norm0    = norm0 + 2.0D0*Q2q(itr)*HQ2q(itr) ! 1/(8BCC)
            qvar     = qvar  +       Q2q(itr)**2        ! 1/(16BBCC)<Q^2>
            pvar     = pvar  +      HQ2q(itr)**2        ! 1/(4CC)   <P^2>
            buff1    = buff1 +      HP2q(itr)* B2q(itr)
            buff2    = buff2 +      HP2q(itr)*HP2q(itr)
            buff3    = buff3 +       B2q(itr)* B2q(itr)
         enddo
         mass  = 0.5D0 *norm**2/norm0
         stiff = 0.5D0         /norm
         OmgSq = 0.25D0*norm   /norm0
         pvar  = pvar          /norm**2
         qvar  = qvar*  (norm  /norm0)**2
         errP  = 1.0D0-buff1**2/(buff2*buff3)
c
c        logging
         write(22,'(a,i4,6(1pg12.3))') 'CHF  : ',
     &            iter,mass,stiff,OmgSq,qvar,pvar,errP
c
c        check for convergency
         if(abs(errP).lt.1E-6) goto 4001
c         
      enddo
 4001 continue
c
c     normalize 
c     ------------------------------------------------------------------
      normP  =    1.0D0/norm
      normQ  =     norm/norm0
c
      do itr=1,TRntrtot
         Q2q(itr) =  Q2q(itr)*normQ
         P2q(itr) = HQ2q(itr)*normP
      enddo
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
c     ------------------------------------------------------------------
      call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
      do itr=1,TRntrtot
         HQ2q(itr) = H2q(itr)*Q2q(itr) + Hresp2q(itr)
      enddo
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp}
c     ------------------------------------------------------------------
      call upd_response_even(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      do itr=1,TRntrtot
         HP2q(itr) = H2q(itr)*P2q(itr) + Hresp2q(itr)
      enddo
c
c     check the norm and evaluate the variances, mass and stiffness
c     ------------------------------------------------------------------
      mass =0.0D0
      stiff=0.0D0
      norm =0.0D0
      norm0=0.0D0
      qvar =0.0D0
      pvar =0.0D0
      do itr=1,TRntrtot
         norm0    = norm0 + 2.0D0*B2q(itr)*P2q(itr)
         norm     = norm  + 2.0D0*Q2q(itr)*P2q(itr)
         mass     = mass  +       Q2q(itr)*HQ2q(itr)
         stiff    = stiff +       P2q(itr)*HP2q(itr)
         qvar     = qvar  +       Q2q(itr)*Q2q(itr)
         pvar     = pvar  +       P2q(itr)*P2q(itr)
      enddo    
      OmgSq = mass*stiff
c
c     logging
      write(21,'(a,7(1pe16.4))') 'CHF_VAL   : ',
     &        mass,stiff,PVar,QVar,OmgSq,norm,norm0
c
c     save matrixelements of the Q- and P-Operator
c     ------------------------------------------------------------------
      do itr=1,TRntrtot
         Q2qCHF(itr,imu) = Q2q(itr)
         P2qCHF(itr,imu) = P2q(itr)
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     crank_ATDHF(B2q,Q2q,P2q,H2q,il,imu)
c
c       full ATDHF
c 
c               [H,[H,Q]] = 4 Omega^2 Q
c
c=======================================================================
      subroutine crank_ATDHF(B2q,Q2q,P2q,H2q,mass,stiff,qvar,pvar,
     &                       il,imu,pathmode)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision P2q(1:MTDTOT),Q2q(1:MTDTOT)
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      double precision mass                   ! collective mass B=1/(2M)
      double precision qvar                   ! <Q^2>-<Q>^2
      double precision pvar                   ! <P^2>
      double precision stiff                  ! stiffness C
      integer          pathmode
c
      integer          itr
      integer          iter
      double precision norm,norm0,OmgSq,normP,normQ
      double precision Hresp2q(1:MTDTOT)
      double precision HQ2q(1:MTDTOT)
      double precision HP2q(1:MTDTOT)
      double precision D2q(1:MTDTOT)
      double precision Dmp2q(1:MTDTOT)
      double precision C2q(1:MTDTOT)
c
      integer          maxiter
      double precision stepQ
c
      maxiter = 1000
      stepQ   = 0.1D0
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      else if (pathmode.ne.1.and.pathmode.ne.0) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: pathmode must be 0 (ATHDF) or 1 (CHF) -- not',
     &             pathmode
         stop 'in linresp'
      endif
c
c     initialize with Constraint-operator and Q-operator
c     ------------------------------------------------------------------
      if (pathmode.eq.0) then
         do itr=1,TRntrtot      ! ATDHF
            Q2q(itr) = Q2qATDHF(itr,imu)
            C2q(itr) = Q2qATDHF(itr,imu)
         enddo
      else                      ! CHF
         do itr=1,TRntrtot
            Q2q(itr) = Q2qCHF(itr,imu)
            C2q(itr) = B2q(itr)
         enddo
      endif
c
c     solve i[H,[H,Q] + h_{Q,resp}] + h_{P,resp} = Beta by iteration
c     ------------------------------------------------------------------
      do iter=1,maxiter
c
         call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
         call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c        evaluate the Liouville operator i[H,Q] + h_{Q,resp}
         mass=0.0D0   ! B
         qvar=0.0D0
         do itr=1,TRntrtot
            HQ2q(itr) = H2q(itr)*Q2q(itr)    + Hresp2q(itr)
            mass      = mass +       
     &                  H2q(itr)*Q2q(itr)*B2q(itr) + 
     &                  Hresp2q(itr)     *B2q(itr)
            qvar      = qvar +   Q2q(itr)*Q2q(itr)
         enddo
c
         call upd_response_even(HQ2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
         call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
c        evaluate the Liouville operator -i[H,P] + h_{P,resp}
         OmgSq=0.0D0   ! 4*Omega^2 = 4BC
         pvar =0.0D0
         do itr=1,TRntrtot
            HP2q(itr) = H2q(itr)*HQ2q(itr) + Hresp2q(itr)
            OmgSq     = OmgSq + 
     &                  H2q(itr)*HQ2q(itr)*HQ2q(itr) + 
     &                  Hresp2q(itr)      *HQ2q(itr)
            pvar      = pvar  +  HQ2q(itr)*HQ2q(itr)
         enddo
         OmgSq = OmgSq/mass   ! 4*Omega^2
         pvar  = pvar/(4.0D0*mass**2)
         stiff = OmgSq/(4.0D0*mass)
c
c        perform the gradient step
         do itr=1,TRntrtot
            D2q(itr) = HP2q(itr) - OmgSq*C2q(itr)
            Q2q(itr) = Q2q(itr)  - stepQ*(HP2q(itr)-OmgSq*C2q(itr))/
     &                                    (H2q(itr)**2+1.0D0)
c     &                                    Dmp2q(itr)
         enddo
c
c        evaluate variances, norms and errors
         norm  =0.0D0
         buff1 =0.0D0
         buff2 =0.0D0
         buff3 =0.0D0
         do itr=1,TRntrtot
            norm     = norm  + 2.0D0*Q2q(itr)*HQ2q(itr) ! should be 2B
            buff1    = buff1 +      HP2q(itr)* C2q(itr)
            buff2    = buff2 +      HP2q(itr)*HP2q(itr)
            buff3    = buff3 +       C2q(itr)* C2q(itr)
         enddo
         norm  = 2.0D0*mass/norm
         errP  = 1.0D0-buff1**2/(buff2*buff3)
c
         do itr=1,TRntrtot
            Q2q(itr) = Q2q(itr)*norm
         enddo
c
         if (pathmode.eq.0) then
            do itr=1,TRntrtot
               C2q(itr) = Q2q(itr)
            enddo
         endif
c
c        logging
         if (pathmode.eq.0) then
            write(22,'(a,i4,6(1pg12.3))') 'ATDHF: ',
     &            iter,mass,stiff,OmgSq,qvar,pvar,errP
         else
            write(22,'(a,i4,6(1pg12.3))') 'CHF  : ',
     &            iter,mass,stiff,OmgSq,qvar,pvar,errP
         endif
c
c        check for convergency
         if(abs(errP).lt.1E-6) goto 4002
c         
      enddo
 4002 continue
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp}
c     ------------------------------------------------------------------
      call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
c     evaluate the normalization commutator <[beta,P]> and normalize P
c     ------------------------------------------------------------------
      mass = 0.0D0  ! B
      do itr=1,TRntrtot
         HQ2q(itr) = H2q(itr)*Q2q(itr) + Hresp2q(itr)
         mass      = mass +       
     &               H2q(itr)*Q2q(itr)*B2q(itr) + 
     &               Hresp2q(itr)     *B2q(itr)
      enddo 
c
      norm=1.0D0/(2.0D0*mass)
      do itr=1,TRntrtot
         P2q(itr)  = HQ2q(itr)*norm
      enddo
c
c     evaluate the normalization commutator <[Q,P]> and normalize Q
c     ------------------------------------------------------------------
      norm = 0.0D0  
      do itr=1,TRntrtot
         norm      = norm + 2.0D0*Q2q(itr)*P2q(itr)
      enddo 
c
      norm=1.0/norm
      do itr=1,TRntrtot
         Q2q(itr)  = Q2q(itr)*norm
      enddo
c
c     evaluate the Liouville operator i[H,Q] + h_{Q,resp} and mass
c     ------------------------------------------------------------------
      call upd_response_odd(Q2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
      call upd_Hresp_odd(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     redynpotr,redynpotp,redynpotz,
     &     respnpotr,respnpotp,respnpotz,
     &     respnmasr,respnmasp,respnmasz,il,imu)
c
      do itr=1,TRntrtot
         HQ2q(itr) =  H2q(itr)*Q2q(itr) + Hresp2q(itr)
      enddo
c
c     evaluate the Liouville operator -i[H,P] + h_{P,resp} and stiffness
c     ------------------------------------------------------------------
      call upd_response_even(P2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     redreffmas,redpeffmas,redzeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      call upd_Hresp_even(Hresp2q,
     &     CPpsi(0,0,1,1,0),CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &     CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),CPvsqf(1,0),
     &     relocpot,reeffmas,
     &     resocpotr,resocpotp,resocpotz,il,imu)
c
      do itr=1,TRntrtot
         HP2q(itr) =  H2q(itr)*P2q(itr) + Hresp2q(itr)
      enddo
c
c     check the norm and evaluate the variances, mass and stiffness
c     ------------------------------------------------------------------
      mass =0.0D0
      stiff=0.0D0
      norm =0.0D0
      norm0=0.0D0
      qvar =0.0D0
      pvar =0.0D0
      do itr=1,TRntrtot
         norm0    = norm0 + 2.0D0*B2q(itr)*P2q(itr)
         norm     = norm  + 2.0D0*Q2q(itr)*P2q(itr)
         mass     = mass  +       Q2q(itr)*HQ2q(itr)
         stiff    = stiff +       P2q(itr)*HP2q(itr)
         qvar     = qvar  +       Q2q(itr)*Q2q(itr)
         pvar     = pvar  +       P2q(itr)*P2q(itr)
      enddo    
      OmgSq = mass*stiff
c
c     logging
      if (pathmode.eq.0) then
         write(21,'(a,7(1pe16.4))') 'ATDHF_VAL : ',
     &        mass,stiff,PVar,QVar,OmgSq,norm,norm0
      else
         write(21,'(a,7(1pe16.4))') 'CHF_VAL   : ',
     &        mass,stiff,PVar,QVar,OmgSq,norm,norm0
      endif
c
c     save matrixelements of the Q- and P-Operator
c     ------------------------------------------------------------------
      if (pathmode.eq.0) then
         do itr=1,TRntrtot
            Q2qATDHF(itr,imu) = Q2q(itr)
            P2qATDHF(itr,imu) = P2q(itr)
         enddo
      else
         do itr=1,TRntrtot
            Q2qCHF(itr,imu) = Q2q(itr)
            P2qCHF(itr,imu) = P2q(itr)
         enddo
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     eval_com_QHQ(Comm,il,imu,iiso)
c
c       evaluates the particle number mass B_N=1/2*<[Q,[N,Q]]>
c
c=======================================================================
      subroutine eval_com_QNQ(N2q,Q2q,mass,il,imu,iiso)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision N2q(1:MTDTOT)          ! qh-matrixelements of N
      double precision Q2q(1:MTDTOT)          ! 2q-matrixelements of Q
      double precision mass                   ! particle number mass
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      integer          iiso                   ! iso-spin-projection
c
      integer          ist,jst,jstmin,jstmax
      integer          trzpar
      integer          itr
      integer          iso 
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     set mass to zero
c     ------------------------------------------------------------------
      mass  = 0.0D0
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif
c
         if (iso.ne.3-iiso) then
            do itr=1,TRntr(ist)
c     
               itr2 = TRidx(ist,itr,4)
c     
               mass = mass + Q2q(itr2)**2*N2q(itr2)
c
            enddo
         endif
c
      enddo
c
c      mass = mass   ! factors from 0.5 from 1/M to B and 2 from the
c                    ! commutator compensate
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a,i2,a)') 
     &  'EVAL_COM_QNQ(): QNQ-commutator evaluated'//
     &  ' (lambda=',il,' mu=',imu,' iso=', iiso,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end

c=======================================================================
c
c     eval_comm_BP(Comm,il,imu,iiso)
c
c       evaluates the commutator -i<[B,P]> including contributions
c       from outside the qp-space.
c
c=======================================================================
      subroutine eval_comm_BP(Comm,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision Comm                   ! variance
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision QlmOper(0:MZD,0:MRD)   ! multipole operator
      double precision pfac,normfac,qlmfac
      double precision peta,qeta,povlp,qovlp
      double precision z,r
      integer          iz,ir
      double precision lastv,lastu,nextv,nextu,rcntv,rcntu
      double precision ui,uj,vi,vj
      integer          itr,ist,jst,jstmin,jstmax
c
      logical          oldcrank
      data             oldcrank/D_OLDCRANK/
c
c     ------------------------------------------------------------------
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     some factors
c     ------------------------------------------------------------------
      normfac = 0.6D0*1.2D0**2*
     &           (npart(1)+npart(2))**(2.0D0/3.0D0)*(npart(1)+npart(2))
      if (imu.eq.0) then
         qlmfac   = 0.8D0*pi/normfac
         if (oldcrank) then
            pfac     = 1.0D0/(CPbeta( 0)-CPbeta(-1))
         else
            pfac     = 1.0D0/(CPbeta(+1)-CPbeta(-1))
         endif
      else if (imu.eq.1) then
         qlmfac   = 0.8D0*pi/(sqrt(2.0D0)*normfac)
         if (abs(CPbeta(0)).gt.0.0D0) then
            pfac=1.0D0/(sqrt(3.0D0)*CPbeta(0))
         else
            pfac=0.0D0
         endif
      else
         write(21,'(a,i2)') 
     &            'ERROR: only mu=0 and mu=1 possible -- not ',imu
         stop 'in linresp'
      endif
c
c     set commutator value to zero
c     ------------------------------------------------------------------
      Comm = 0.0D0
c        
c     initialize multipole operator
c     ------------------------------------------------------------------
      do ir=0,nr
         r    = dr*(ir+sr)
         do iz=0,nz
            z    = ((iz)+sz-icm)*dz - CPzcm(0)
            QlmOper(iz,ir) = qlmfac*ylm(z,r,2,imu) 
         enddo
      enddo
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c     
c        calculate action of P operator 
c        ===============================================================
c
c        P_0 (generator of translation along the axial path, asym. def.)
c        --------------------------------------------------------------
         if (imu.eq.0.and.oldcrank) then   
c
            do ir=0,nr
               do iz=0,nz
                  LRbuff1(iz,ir,1,1)=pfac*
     &               (CPpsi(iz,ir,1,ist, 0)-CPpsi(iz,ir,1,ist,-1))
                  LRbuff1(iz,ir,2,1)=pfac*
     &               (CPpsi(iz,ir,2,ist, 0)-CPpsi(iz,ir,2,ist,-1))
                  LRbuff1(iz,ir,1,2)=0.0D0
                  LRbuff1(iz,ir,2,2)=0.0D0
               enddo
            enddo
c
c        P_0 (generator of translation along the axial path, sym. def.)
c        --------------------------------------------------------------
         else if (imu.eq.0) then 
c
            do ir=0,nr
               do iz=0,nz
                  LRbuff1(iz,ir,1,1)=pfac*
     &               (CPpsi(iz,ir,1,ist,+1)-CPpsi(iz,ir,1,ist,-1))
                  LRbuff1(iz,ir,2,1)=pfac*
     &               (CPpsi(iz,ir,2,ist,+1)-CPpsi(iz,ir,2,ist,-1))
                  LRbuff1(iz,ir,1,2)=0.0D0
                  LRbuff1(iz,ir,2,2)=0.0D0
               enddo
            enddo 
c
c        P_1 (generator of rotation,  J_y = L_y + S_y)
c        --------------------------------------------------------------
         else if (imu.eq.1) then
c
            do ir=0,nr
               r = (ir+sr)*dr
               do iz=0,nz 
                  z = (iz+sz-icm)*dz - CPzcm(0)
c     
c                 ! upward laddering components of J_y
                  LRbuff1(iz,ir,1,1) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,1,ist) 
     &                                  +z*CPdppsi(iz,ir,1,ist) 
     &                                  +r*CPdzpsi(iz,ir,1,ist)
     &                                  -  CPpsi(iz,ir,2,ist,0))
                  LRbuff1(iz,ir,2,1) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,2,ist) 
     &                                  +z*CPdppsi(iz,ir,2,ist) 
     &                                  +r*CPdzpsi(iz,ir,2,ist))
c
c                 ! downward laddering components of J_y
                  LRbuff1(iz,ir,1,2) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,1,ist) 
     &                                  -z*CPdppsi(iz,ir,1,ist) 
     &                                  +r*CPdzpsi(iz,ir,1,ist))
                  LRbuff1(iz,ir,2,2) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,2,ist) 
     &                                  -z*CPdppsi(iz,ir,2,ist) 
     &                                  +r*CPdzpsi(iz,ir,2,ist)
     &                                  +  CPpsi(iz,ir,1,ist,0))
c
               enddo
            enddo
c
         endif           
c     
c        calculate action of Q operator 
c        ===============================================================
         if (imu.eq.0) then
            do ir=0,nr
               do iz=0,nz
                  LRbuff2(iz,ir,1,1)=
     &                 QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  LRbuff2(iz,ir,2,1)=
     &                 QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
                  LRbuff2(iz,ir,1,2)= 0.0D0
                  LRbuff2(iz,ir,2,2)= 0.0D0
               enddo
            enddo
         else
            do ir=0,nr
               do iz=0,nz
                  LRbuff2(iz,ir,1,1)=
     &                 QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  LRbuff2(iz,ir,2,1)=
     &                 QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
                  LRbuff2(iz,ir,1,2)=
     &                 QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  LRbuff2(iz,ir,2,2)=
     &                 QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
               enddo
            enddo
         endif
c
c        adjust lower boundary in case of no z-symmetry
c        ---------------------------------------------------------------
         if (isym.eq.0) then
            do ir=0,nr
               LRbuff1(0,ir,1,1) = 0.0D0
               LRbuff1(0,ir,2,1) = 0.0D0
               LRbuff1(0,ir,1,2) = 0.0D0
               LRbuff1(0,ir,2,2) = 0.0D0
               LRbuff2(0,ir,1,1) = 0.0D0
               LRbuff2(0,ir,2,1) = 0.0D0
               LRbuff2(0,ir,1,2) = 0.0D0
               LRbuff2(0,ir,2,2) = 0.0D0
            enddo
         endif
c     
c        evaluate 2q-matrixlelements and 0q-part
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            itp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)
            qeta = uj*vi + vj*ui
            peta = uj*vi - vj*ui
c
            if (itp.eq.+1) then
               if (msgn.ge.0) then
                  povlp=wfovlp(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,1))  
                  call wfproj(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,1),
     &                        povlp)
                  qovlp=wfovlp(CPpsi(0,0,1,jst,0),LRbuff2(0,0,1,1))  
                  call wfproj(CPpsi(0,0,1,jst,0),LRbuff2(0,0,1,1),
     &                        qovlp)
               else       
                  povlp=wfovlp(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2))  
                  call wfproj(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2),
     &                        povlp)
                  qovlp=wfovlp(CPpsi(0,0,1,jst,0),LRbuff2(0,0,1,2))  
                  call wfproj(CPpsi(0,0,1,jst,0),LRbuff2(0,0,1,2),
     &                        qovlp)
               endif
            else
               povlp=wfovlptr(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2))  
               call wfprojtr(CPpsi(0,0,1,jst,0),LRbuff1(0,0,1,2),
     &                       povlp)
               qovlp=wfovlptr(CPpsi(0,0,1,jst,0),LRbuff2(0,0,1,2))  
               call wfprojtr(CPpsi(0,0,1,jst,0),LRbuff2(0,0,1,2),
     &                       qovlp)
            endif
c
            if (ist.eq.jst.and.itp.eq.+1) then
c  
               if (oldcrank) then
                  lastv   = sqrtb(CPvsqf(ist,-1))
                  lastu   = sqrtb(one-CPvsqf(ist,-1))
                  nextv   = sqrtb(CPvsqf(ist,0))
                  nextu   = sqrtb(one-CPvsqf(ist,0))
                  rcntv   = 0.5D0*(lastv+nextv)
                  rcntu   = 0.5D0*(lastu+nextu)
                  if (rcntv.ge.0.5D0) then
                     povlp =  -(nextu-lastu)/rcntv*pfac
                  else
                     povlp =  +(nextv-lastv)/rcntu*pfac
                  endif
               else
                  lastv   = sqrtb(CPvsqf(ist,-1))
                  lastu   = sqrtb(one-CPvsqf(ist,-1))
                  rcntv   = sqrtb(CPvsqf(ist,0)) 
                  rcntu   = sqrtb(one-CPvsqf(ist,0))
                  nextv   = sqrtb(CPvsqf(ist,+1))
                  nextu   = sqrtb(one-CPvsqf(ist,+1))
                  if (rcntv.ge.0.5D0) then
                     povlp =  -(nextu-lastu)/rcntv*pfac
                  else
                     povlp =  +(nextv-lastv)/rcntu*pfac
                  endif
               endif
c
               Comm = Comm +  qeta*     povlp*qovlp
c
            else
c
               Comm = Comm +  qeta*peta*povlp*qovlp   
c
            endif
c
         enddo
c
c        accumulate contributions from outside qp-space
c        ---------------------------------------------------------------
         Comm = Comm + 2.0D0*CPvsqf(ist,0)*(
     &                     wfovlp(LRbuff2(0,0,1,1),LRbuff1(0,0,1,1))
     &                    +wfovlp(LRbuff2(0,0,1,2),LRbuff1(0,0,1,2)) )
c
      enddo
c
      Comm = 2.0D0*Comm
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'EVAL_COMM_QP(): QP-commutator evaluated'//
     &  ' (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     ini_matrix_N(N2q,N0q,NVar,il,imu,iiso)
c
c       initialize the particle number operator
c
c       and calculate 2-quasiparticle transition matrixelements
c
c       2qp-matrixelements    N^{2q}_ji = eta^(+)_ji <j|N|i>
c       expectation           <N>
c       variance              <N^2> - <N>^2
c
c=======================================================================
      subroutine ini_matrix_N(N2q,N0q,NVar,il,imu,iiso)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision N2q(1:MTDTOT)          ! N^{qh}_ii + N^{qh}_jj
      double precision N0q                    ! 0q-part of N
      double precision NVar                   ! variance
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      integer          iiso                   ! iso-spin-projection
c
      integer          ist,jst
      integer          itr,itr2
      integer          iso,msgn
      double precision occ,eta,xi
      double precision vi,ui,vj,uj
c
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      else if (iiso.lt.0.or.iiso.gt.2) then
         write(21,'(a,i2)') 
     &            'ERROR: isospin index has to be 0,1 or 2 -- not ',iiso
         stop 'in linresp' 
      endif
c
c     set 0qp-part and variance to zero
c     ------------------------------------------------------------------
      NVar = 0.0D0
      N0q  = 0.0D0
      do itr2=1,TRntrtot
         N2q(itr2)=0.0D0
      enddo
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif    
c     
c        evaluate variance, 0q-part and N2q=(N^{qh}_{ii} + N^{qh}_{jj})
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            itp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)
            eta  = uj*vi + vj*ui
            occ  = 2.0D0*vi**2
c
c           only diagonal componets contribute
            if (ist.eq.jst.and.itp.eq.+1.and.iso.ne.3-iiso) then 
               NVar      = NVar + (eta)**2
               N0q       = N0q  +  occ
            endif

            N2q(itr2) = (ui*ui - vi*vi + uj*uj - vj*vj) 
c
         enddo
c
      enddo
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a,i2,a)') 
     &  'INI_MATRIX_N(): N-operator initialized'//
     &  '  (lambda=',il,' mu=',imu,' iso=',iiso,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     ini_matrix_B(B2q,B0q,BVar,il,imu,iiso,ipath)
c
c       initialize the deformation measurement-operator -- beta
c
c                   4 pi    Q_{l,+m} + (-1)^m Q_{l,-m}
c             Q = --------  -------------------------- 
c                  5 A r^2              2
c
c                   4 pi    Q_{l,+m} + Q^*_{l,m}
c               = --------  -------------------- 
c                  5 A r^2           2
c
c       r^2  = RMS radius       =  3/5 * R^2
c       A    = nucleon number
c       Q_lm = R^l Y_lm
c
c       and calculate
c
c       2qp-matrixelements    Q^{2q}_ji = eta^(+)_ji <j|Q|i>
c       expectation           <Q>
c       variance              <Q^2> - <Q>^2
c
c=======================================================================
      subroutine ini_matrix_B(B2q,B0q,BVar,il,imu,iiso,ipath)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision B2q(1:MTDTOT)          ! 2qp-matrixelements
      double precision B0q                    ! 0qp-part
      double precision BVar                   ! variance
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
      integer          iiso                   ! iso-spin-projection
      integer          ipath                  ! path index
c 
      double precision QlmOper(0:MZD,0:MRD)   ! multipole operator
      integer          ist,jst,jstmin,jstmax
      integer          itr,itr2
      double precision qlmfac,normfac
      double precision r,z
      integer          ir,iz,iso,msgn,itp
      double precision ovlp,occ,eta
      double precision vi,ui,vj,uj
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
      if (ipath.gt.1.or.ipath.lt.-1) then
         write(21,'(a,i5,a)') 
     &            'ERROR: invalid path index (',ipath,
     &            ') in ini_matrix_beta'
         stop 'in linresp'
      endif
c
c     calculate normalization of Q_mu for given isospin (0: total 1:p 2:n)
c     ------------------------------------------------------------------
      if (iiso.eq.0) then
         normfac = 0.6D0*1.2D0**2*
     &           (npart(1)+npart(2))**(2.0D0/3.0D0)*(npart(1)+npart(2))
      else if (iiso.eq.1) then
         normfac = 0.6D0*1.2D0**2*
     &           (npart(1)+npart(2))**(2.0D0/3.0D0)*npart(1)
      else if (iiso.eq.2) then
         normfac = 0.6D0*1.2D0**2*
     &           (npart(1)+npart(2))**(2.0D0/3.0D0)*npart(2)
      endif
c        
c     initialize multipole operator
c     ------------------------------------------------------------------
      if (imu.eq.0) then
         qlmfac   = 0.8D0*pi/normfac
      else
         qlmfac   = 0.8D0*pi/(sqrt(2.0D0)*normfac)
      endif
c
      do ir=0,nr
         r    = dr*(ir+sr)
         do iz=0,nz
            z    = ((iz)+sz-icm)*dz - CPzcm(ipath)
            QlmOper(iz,ir) = qlmfac*ylm(z,r,il,imu) 
         enddo
      enddo
c
c     adjust lower boundary in case of no z-symmetry
c     ------------------------------------------------------------------
      if (isym.eq.0) then
         do ir=0,nr
            QlmOper(0,ir) = 0.0D0
         enddo
      endif
c
c     set variance, 0q-part and 2q-matrixelements to zero
c     ------------------------------------------------------------------
      BVar = 0.0D0
      B0q  = 0.0D0
      do itr2=1,TRntrtot
         B2q(itr2)=0.0D0
      enddo
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif         
c     
c        calculate action of Q operator
c        ---------------------------------------------------------------
         if (iso.ne.3-iiso) then
            if (imu.eq.0) then
               do ir=0,nr
                  do iz=0,nz
                     LRbuff1(iz,ir,1,1)=
     &                      QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,ipath)
                     LRbuff1(iz,ir,2,1)=
     &                      QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,ipath)
                     LRbuff1(iz,ir,1,2)= 0.0D0
                     LRbuff1(iz,ir,2,2)= 0.0D0
                  enddo
               enddo
            else
               do ir=0,nr
                  do iz=0,nz
                     LRbuff1(iz,ir,1,1)=
     &                      QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,ipath)
                     LRbuff1(iz,ir,2,1)=
     &                      QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,ipath)
                     LRbuff1(iz,ir,1,2)=
     &                      QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,ipath)
                     LRbuff1(iz,ir,2,2)=
     &                      QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,ipath)
                  enddo
               enddo
            endif
         else
            do ir=0,nr
               do iz=0,nz
                  LRbuff1(iz,ir,1,1)= 0.0D0
                  LRbuff1(iz,ir,2,1)= 0.0D0
                  LRbuff1(iz,ir,1,2)= 0.0D0
                  LRbuff1(iz,ir,2,2)= 0.0D0
               enddo
            enddo
         endif
c     
c        evaluate 2q-matrixlelements and 0q-part
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            itp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,ipath))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,ipath))
            uj   = sqrtb(one-vj*vj)
            eta  = uj*vi + vj*ui
            occ  = 2.0D0*vi**2
c
            if (itp.eq.+1) then
               if (msgn.ge.0) then
                  ovlp=wfovlp(CPpsi(0,0,1,jst,ipath),LRbuff1(0,0,1,1))  
                  call wfproj(CPpsi(0,0,1,jst,ipath),LRbuff1(0,0,1,1),
     &                        ovlp)
               else       
                  ovlp=wfovlp(CPpsi(0,0,1,jst,ipath),LRbuff1(0,0,1,2))  
                  call wfproj(CPpsi(0,0,1,jst,ipath),LRbuff1(0,0,1,2),
     &                        ovlp)
               endif
            else
               ovlp=wfovlptr(CPpsi(0,0,1,jst,ipath),LRbuff1(0,0,1,2))  
               call wfprojtr(CPpsi(0,0,1,jst,ipath),LRbuff1(0,0,1,2),
     &                       ovlp)
            endif
c
            if (ist.eq.jst.and.itp.eq.+1) then 
               B2q(itr2) =         eta * ovlp
               BVar      = BVar + (eta * ovlp)**2
               B0q       = B0q  +  occ * ovlp
            else
               B2q(itr2) =         eta * ovlp
               BVar      = BVar + (eta * ovlp)**2
            endif
c
         enddo
c
c        accumulate contributions from outside qp-space
c        ---------------------------------------------------------------
         BVar = BVar + 2.0D0*CPvsqf(ist,ipath)*(
     &                      wfovlp(LRbuff1(0,0,1,1),LRbuff1(0,0,1,1))
     &                    + wfovlp(LRbuff1(0,0,1,2),LRbuff1(0,0,1,2)) )
c
      enddo
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a,i2,a)') 
     &  'INI_MATRIX_B(): B-operator initialized'//
     &  '  (lambda=',il,' mu=',imu,' iso=',iiso,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     ini_matrix_H(H2q,H0q,il,imu,iiso)
c
c     initialize the quasiparticle energies
c
c=======================================================================
      subroutine ini_matrix_H(H2q,H0q,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTDTOT)          ! e_i + e_j
      double precision H0q                    ! 0qp-part of H (energy)
      integer          il                     ! angular momentum
      integer          imu 

      integer          ist,jst
      integer          itr,itr2
      integer          iso 
      double precision deltai,edifi,equasi,deltaj,edifj,equasj
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     calculate 0qp and 2qp-matrixelements of beta-operator
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     set 0q-part and qh-matrixelements to zero
c     ------------------------------------------------------------------
      H0q  = 0.0D0
      do itr=1,TRntrtot
         H2q(itr) = 0.0D0
      enddo
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso=1
         else
            iso=2
         endif
c
         do itr=1,TRntr(ist)
            jst =TRidx(ist,itr,1)
            itr2=TRidx(ist,itr,4)
c
            deltai = CPdeltaf(ist,0)*CPwtpair(ist,0)
            edifi  = CPspenrg(ist,0)-CPfer(iso,0)
            equasi = sqrtb(edifi*edifi+deltai*deltai)
c
            deltaj = CPdeltaf(jst,0)*CPwtpair(jst,0)
            edifj  = CPspenrg(jst,0)-CPfer(iso,0)
            equasj = sqrtb(edifj*edifj+deltaj*deltaj)
c                  
            if (ist.ne.jst) then
               H2q(itr2) = equasi + equasj
            else
               H2q(itr2) = 2.0D0*equasi
               H0q       = H0q + 2.0D0*equasi
            endif
c     
         enddo
c     
      enddo  
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'INI_MATRIX_H(): H0-operator initialized'//
     &  ' (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     ini_transitions(il,imu)
c
c     count the number of transitions (matrixelements) for a given
c     cranking mode
c
c=======================================================================
      subroutine ini_transtab(il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K

      integer          ist,jst,jstmin,jstmax
      integer          trzpar,maxtr
      integer          itr,itr2
      integer          iso
c
c     calculate z-parity of transition 
c     ------------------------------------------------------------------
      trzpar=(-1)**(imu+il)      
c
c     calculate 0qp and 2qp-matrixelements of beta-operator
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
c        
c     the case mu=0
c     ------------------------------------------------------------------
      else if (imu.eq.0) then
c        
c        loop over all sp-states 
c        ---------------------------------------------------------------
         itr   = 1
         itr2  = 1
         maxtr = 0
         do ist=1,nst(1)+nst(2)
c
            if (ist.le.nst(1)) then
               jstmin = 1
               jstmax = nst(1)
               iso    = 1
            else
               jstmin = nst(1)+1
               jstmax = nst(1)+nst(2)
               iso    = 2
            endif
c
            do jst=jstmin,jstmax
               if (mpsi(jst)-mpsi(ist).eq.0)  then
               if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
                  if (itr.le.MTD) then
                     TRidx(ist,itr,1) = jst  ! final state
                     TRidx(ist,itr,2) = +1   ! from original space
                     TRidx(ist,itr,3) = 0    ! angular momentum transfer
                     TRidx(ist,itr,4) = itr2 ! global transition index
                  endif
                  itr              = itr  + 1
                  itr2             = itr2 + 1
               endif
               endif
            enddo    
c
            TRntr(ist) = itr - 1
            maxtr      = MAX(itr-1,maxtr)
            itr        = 1
c
         enddo
c        
c     the case mu>0
c     ------------------------------------------------------------------
      else
c        
c        loop over all sp-states 
c        ---------------------------------------------------------------
         itr   = 1
         itr2  = 1
         maxtr = 0
         do ist=1,nst(1)+nst(2)
c
            if (ist.le.nst(1)) then
               jstmin = 1
               jstmax = nst(1)
               iso    = 1
            else
               jstmin = nst(1)+1
               jstmax = nst(1)+nst(2)
               iso    = 2
            endif
c
            do jst=jstmin,jstmax
c              the upward laddering transitions
               if (mpsi(jst)-mpsi(ist).eq.imu)  then
               if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
                  if (itr.le.MTD) then
                     TRidx(ist,itr,1) = jst  ! final state
                     TRidx(ist,itr,2) = +1   ! from original space
                     TRidx(ist,itr,3) = imu  ! angular momentum transfer
                     TRidx(ist,itr,4) = itr2 ! global transition index
                  endif
                  itr              = itr  + 1
                  itr2             = itr2 + 1
               endif
c              the downward laddering transitions inside the K>0 space
               else if (mpsi(jst)-mpsi(ist).eq.-imu)  then
               if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
                  if (itr.le.MTD) then
                     TRidx(ist,itr,1) = jst  ! final state
                     TRidx(ist,itr,2) = +1   ! from original space
                     TRidx(ist,itr,3) = -imu ! angular momentum transfer
                     TRidx(ist,itr,4) = itr2 ! global transition index
                  endif
                  itr              = itr  + 1
                  itr2             = itr2 + 1
               endif
c              the downward laddering components to the K<0 space
               else if (mpsi(ist)+mpsi(jst)+1.eq.imu)  then
               if (npar(jst)*npar(ist)*trzpar*isym.ne.+1) then
                  if (itr.le.MTD) then
                     TRidx(ist,itr,1) = jst  ! final state
                     TRidx(ist,itr,2) = -1   ! from time reversed space
                     TRidx(ist,itr,3) = -imu ! angular momentum transfer
                     TRidx(ist,itr,4) = itr2 ! global transition index
                  endif
                  itr              = itr  + 1
                  itr2             = itr2 + 1
               endif
               endif
            enddo    
c
            TRntr(ist) = itr - 1
            maxtr      = MAX(itr-1,maxtr)
            itr        = 1
c
         enddo
      endif
c
      TRntrtot = itr2 - 1
c
c     ------------------------------------------------------------------
c
      if (maxtr.gt.MTD) then
         write(21,'(a,6i5,a,6i5)') 
     &        'ERROR: number of transition too big (enhance MTD) ',
     &        ' maxtr=',maxtr,' MTD=',MTD
         stop 'in linresp'    
      endif
c
      if (TRntrtot.gt.MTDTOT) then
         write(21,'(a,6i5,a,6i5)') 
     &        'ERROR: number of transition too big (enhance MTDTOT) ',
     &        ' TRntr=',TRntrtot,' MTD=',MTDTOT
         stop 'in linresp'    
      endif
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a,i2,a,i6,a)') 
     &  'INI_TRANSTAB(): transitions initialized'//
     &  ' (lambda=',il,' mu=',imu,' pi= ',trzpar,' ntr=',TRntrtot,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
c
      end


c=======================================================================
c
c     upd_response_odd(Q2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
c                      RESP_Ar,RESP_Ap,RESP_Az,
c                      RESP_Sr,RESP_Sp,RESP_Sz,
c                      RESP_Cr,RESP_Cp,RESP_Cz,il,imu)
c                       
c     updates the time-odd mean-field potentials (and densities)
c
c=======================================================================
      subroutine upd_response_odd(Q2q,
     &                            phi,drphi,dzphi,dpphi,ddphi,occphi,
     &                            RESP_Ar,RESP_Ap,RESP_Az,
     &                            RESP_Sr,RESP_Sp,RESP_Sz,
     &                            RESP_Cr,RESP_Cp,RESP_Cz,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision  Q2q(1:MTDTOT)
      double precision  phi(0:MZD,0:MRD,MSD,MPSID)
      double precision  drphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dzphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dpphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  ddphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  occphi(MPSID)
c
      double precision  RESP_Ar(0:MZD,0:MRD,2)
      double precision  RESP_Ap(0:MZD,0:MRD,2)
      double precision  RESP_Az(0:MZD,0:MRD,2)
      double precision  RESP_Sr(0:MZD,0:MRD,2)
      double precision  RESP_Sp(0:MZD,0:MRD,2)
      double precision  RESP_Sz(0:MZD,0:MRD,2)
      double precision  RESP_Cr(0:MZD,0:MRD,2)
      double precision  RESP_Cp(0:MZD,0:MRD,2)
      double precision  RESP_Cz(0:MZD,0:MRD,2)
      double precision  RESP_Delta(0:MZD,0:MRD,2)
c
      double precision  resp_spnr(0:MZD,0:MRD,2)
      double precision  resp_spnp(0:MZD,0:MRD,2)
      double precision  resp_spnz(0:MZD,0:MRD,2)
      double precision  resp_skdr(0:MZD,0:MRD,2)
      double precision  resp_skdp(0:MZD,0:MRD,2)
      double precision  resp_skdz(0:MZD,0:MRD,2)
      double precision  resp_curr(0:MZD,0:MRD,2)
      double precision  resp_curp(0:MZD,0:MRD,2)
      double precision  resp_curz(0:MZD,0:MRD,2)
      double precision  resp_paio(0:MZD,0:MRD,2)
      double precision  resp_rotjr(0:MZD,0:MRD,2)
      double precision  resp_rotjp(0:MZD,0:MRD,2)
      double precision  resp_rotjz(0:MZD,0:MRD,2)
      double precision  resp_rotsr(0:MZD,0:MRD,2)
      double precision  resp_rotsp(0:MZD,0:MRD,2)
      double precision  resp_rotsz(0:MZD,0:MRD,2)
c
      integer          ist,jst,jtp,msgn
      integer          itr,iz,ir
      integer          iso
      double precision eta,xi,wgt1,wgt2,wgt3
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     set response densities to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do ir=0,nr
            do iz=0,nz
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
            enddo
         enddo
      enddo
c
c     loop all single particle states
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif                  
c       
c        loop all transitions
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            jtp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0)) 
            uj   = sqrtb(one-vj*vj)
c
c           if the final state is from the original space
            if (jtp.eq.+1) then
c
               if (msgn.eq.0) then
                  wgt1=   0.0D0                           ! for   e^{i*mu*phi}   
                  wgt2=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for -ie^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2) 
               else if (msgn.gt.0) then
                  wgt1=  +2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for   e^{i*mu*phi}
                  wgt2=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for -ie^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2)
               else 
                  wgt1=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for   e^{i*mu*phi}  
                  wgt2=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for -ie^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2) 
               endif
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
c
c
                  resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) 
     &                                                  + 0.5D0*wgt1*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
c
                  resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) ) 
c
c
                  resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(
     &                 - dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(
     &                 + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist) )
c     
                  resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(
     &                 + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist) 
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist) 
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist) 
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
               enddo
               enddo
c
c           if the final state is from the time-reversed space
            else
c
               wgt1=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for   e^{i*mu*phi}  
               wgt2=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for -ie^{i*mu*phi} 
               wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2) 
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 - phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 - phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
c
c
                  resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) 
     &                                                  + 0.5D0*wgt1*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
c
                  resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) ) 
c
c
                  resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(
     &                 - dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist) )
c     
                  resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(
     &                 + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist) 
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist) 
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist) 
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
               enddo
               enddo
c
            endif
c
         enddo    ! end of loop over all transitions from state ist 
      enddo       ! end of loop over all states ist
c
c     recouple the Skyrme-Parameters
c     ------------------------------------------------------------------
      b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
      b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
      b4  = 0.500D0*       t4
      b5  = 0.125D0*(   x1*t1 +       x2*t2)  ! 2*b5
      b5p = 0.125D0*(      t1 -          t2)  ! 2*b5'
c
      wgt_b1  = -2.0D0*(b1-b1p)  
      wgt_b1c = -2.0D0* b1  
      if (b4p.eq.0.0D0) then    ! no isotopic spin-orbit force
         wgt_b4  = - 2.0D0*b4
         wgt_b4c = -       b4 
      else
         wgt_b4  = -(b4+b4p)
         wgt_b4c = - b4
      endif
      if(ifls.eq.1) then
         wgt_b5  =  (b5-b5p)    ! b5 already includes the factor 2
         wgt_b5c =   b5
      else
         wgt_b5  = 0.0D0
         wgt_b5c = 0.0D0
      endif
c
c     set the mean-fields to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
            RESP_Ar(iz,ir,iso) = 0.0D0           
            RESP_Ap(iz,ir,iso) = 0.0D0
            RESP_Az(iz,ir,iso) = 0.0D0
            RESP_Sr(iz,ir,iso) = 0.0D0
            RESP_Sp(iz,ir,iso) = 0.0D0
            RESP_Sz(iz,ir,iso) = 0.0D0
            RESP_Cr(iz,ir,iso) = 0.0D0
            RESP_Cp(iz,ir,iso) = 0.0D0
            RESP_Cz(iz,ir,iso) = 0.0D0
         enddo   
         enddo
      enddo
c
c     construct the mean-fields from the densities
c     ------------------------------------------------------------------  
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
c
            RESP_Ar(iz,ir,iso) = RESP_Ar(iz,ir,iso) 
     &                         + wgt_b1 *resp_curr(iz,ir,iso)
     &                         + wgt_b1c*resp_curr(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotsr(iz,ir,iso)
     &                         + wgt_b4c*resp_rotsr(iz,ir,3-iso)
            RESP_Ap(iz,ir,iso) = RESP_Ap(iz,ir,iso) 
     &                         + wgt_b1 *resp_curp(iz,ir,iso)
     &                         + wgt_b1c*resp_curp(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotsp(iz,ir,iso)
     &                         + wgt_b4c*resp_rotsp(iz,ir,3-iso)
            RESP_Az(iz,ir,iso) = RESP_Az(iz,ir,iso) 
     &                         + wgt_b1 *resp_curz(iz,ir,iso)
     &                         + wgt_b1c*resp_curz(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotsz(iz,ir,iso)
     &                         + wgt_b4c*resp_rotsz(iz,ir,3-iso)
c
            RESP_Sr(iz,ir,iso) = RESP_Sr(iz,ir,iso) 
     &                         + wgt_b5 *resp_skdr(iz,ir,iso)
     &                         + wgt_b5c*resp_skdr(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotjr(iz,ir,iso)
     &                         + wgt_b4c*resp_rotjr(iz,ir,3-iso)
            RESP_Sp(iz,ir,iso) = RESP_Sp(iz,ir,iso) 
     &                         + wgt_b5 *resp_skdp(iz,ir,iso)
     &                         + wgt_b5c*resp_skdp(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotjp(iz,ir,iso)
     &                         + wgt_b4c*resp_rotjp(iz,ir,3-iso)
            RESP_Sz(iz,ir,iso) = RESP_Sz(iz,ir,iso) 
     &                         + wgt_b5 *resp_skdz(iz,ir,iso)
     &                         + wgt_b5c*resp_skdz(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotjz(iz,ir,iso)
     &                         + wgt_b4c*resp_rotjz(iz,ir,3-iso)
c
            RESP_Cr(iz,ir,iso) = RESP_Cr(iz,ir,iso) 
     &                         + wgt_b5 *resp_spnr(iz,ir,iso)
     &                         + wgt_b5c*resp_spnr(iz,ir,3-iso)
            RESP_Cp(iz,ir,iso) = RESP_Cp(iz,ir,iso) 
     &                         + wgt_b5 *resp_spnp(iz,ir,iso)
     &                         + wgt_b5c*resp_spnp(iz,ir,3-iso)
            RESP_Cz(iz,ir,iso) = RESP_Cz(iz,ir,iso) 
     &                         + wgt_b5 *resp_spnz(iz,ir,iso)
     &                         + wgt_b5c*resp_spnz(iz,ir,3-iso)
c
         enddo   
         enddo
      enddo  
c
c     ALL FIELDS NOW CONTAIN THE (r,z)-PART OF THE RESPONSE FIELDS
c     THE ANGULAR DEPENDENCY READS
c
c     Cr,Cz,Sr,Sz,Ap  ~  sin(mu*phi)*X(r,z)
c     Cp,Sp,Ar,Az     ~  cos(mu*phi)*X(r,z)
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'UPD_RESP_ODD(): time-odd response updated'//
     &  '  (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
c
      end


c=======================================================================
c
c     upd_response_even(P2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
c                       RESP_U,RESP_B,RESP_drB,RESP_dzB,RESP_Wz,RESP_Wr,
c                       
c
c     count the number of transitions (matrixelements) for a given
c     cranking mode
c
c=======================================================================
      subroutine upd_response_Ptest(P2q,
     &                             phi,drphi,dzphi,dpphi,ddphi,occphi,
     &                             RESP_U,
     &                             RESP_B,RESP_drB,RESP_dpB,RESP_dzB,
     &                             RESP_Wr,RESP_Wp,RESP_Wz,il,imu,
     &                             a1,a2,a3,a4,a5,a6)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision  P2q(1:MTDTOT)
      double precision  phi(0:MZD,0:MRD,MSD,MPSID)
      double precision  drphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dzphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dpphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  ddphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  occphi(MPSID)
c
      double precision  RESP_U(0:MZD,0:MRD,2)
      double precision  RESP_B(0:MZD,0:MRD,2)
      double precision  RESP_drB(0:MZD,0:MRD,2)
      double precision  RESP_dpB(0:MZD,0:MRD,2)
      double precision  RESP_dzB(0:MZD,0:MRD,2)
      double precision  RESP_Wr(0:MZD,0:MRD,2)
      double precision  RESP_Wp(0:MZD,0:MRD,2)
      double precision  RESP_Wz(0:MZD,0:MRD,2)
      double precision  RESP_Delta(0:MZD,0:MRD,2)
c
      double precision  resp_rho(0:MZD,0:MRD,2)
      double precision  resp_tau(0:MZD,0:MRD,2)
      double precision  resp_socr(0:MZD,0:MRD,2)
      double precision  resp_socp(0:MZD,0:MRD,2)
      double precision  resp_socz(0:MZD,0:MRD,2)
      double precision  resp_Grhor(0:MZD,0:MRD,2)
      double precision  resp_Grhop(0:MZD,0:MRD,2)
      double precision  resp_Grhoz(0:MZD,0:MRD,2)
      double precision  resp_Drho(0:MZD,0:MRD,2)
      double precision  resp_divsoc(0:MZD,0:MRD,2)
      double precision  resp_pair(0:MZD,0:MRD,2)
c
      integer          ist,jst,jtp
      integer          itr,iz,ir
      integer          iso
      double precision eta,xi,wgt1,wgt2,wgt3,a1,a2,a3,a4,a5,a6
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     set response densities to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do iz=0,nz
            do ir=0,nr
               resp_rho(iz,ir,iso)    = 0.0D0
               resp_tau(iz,ir,iso)    = 0.0D0
               resp_socr(iz,ir,iso)   = 0.0D0
               resp_socp(iz,ir,iso)   = 0.0D0
               resp_socz(iz,ir,iso)   = 0.0D0
               resp_Grhor(iz,ir,iso)  = 0.0D0
               resp_Grhop(iz,ir,iso)  = 0.0D0
               resp_Grhoz(iz,ir,iso)  = 0.0D0
               resp_Drho(iz,ir,iso)   = 0.0D0
               resp_divsoc(iz,ir,iso) = 0.0D0
               resp_pair(iz,ir,iso)   = 0.0D0
            enddo
         enddo
      enddo
c
c     loop all single particle states
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif                  
c       
c        loop all transitions
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            jtp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)
c
c           if the final state is from the original space
            if (jtp.eq.+1) then
c
               if (msgn.eq.0) then
                  wgt1=   0.0D0                           ! for -ie^{i*mu*phi}   
                  wgt2=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2) 
               else if (msgn.gt.0) then
                  wgt1=   2.0D0*a1*(uj*vi + vj*ui)*P2q(itr2) ! for -ie^{i*mu*phi}
                  wgt2=   2.0D0*a2*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2)
               else 
                  wgt1=   2.0D0*a3*(uj*vi + vj*ui)*P2q(itr2) ! for -ie^{i*mu*phi}  
                  wgt2=   2.0D0*a4*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2) 
               endif
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
c
                  resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 
     &                                                   0.5D0*wgt2*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 
     &                                                 0.5D0*wgt1*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 
     &                                                     0.5D0*wgt2*(
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )                  
c
                  resp_Grhor(iz,ir,iso) = resp_Grhor(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )                   
c
                  resp_Grhop(iz,ir,iso) = resp_Grhop(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_Grhoz(iz,ir,iso) = resp_Grhoz(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
c                 2*resp_tau has to be added later
                  resp_Drho(iz,ir,iso) = resp_Drho(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*ddphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*ddphi(iz,ir,2,ist)
     &                 + ddphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + ddphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(
     &                 + drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist) 
     &                 - dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
c
               enddo
               enddo
c
c           if the final state is from the time-reversed space
            else
c
               wgt1=   2.0D0*a5*(uj*vi + vj*ui)*P2q(itr2) ! for -ie^{i*mu*phi}  
               wgt2=   2.0D0*a6*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
               wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2) 
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
c
                  resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 
     &                                                     0.5D0*wgt2*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 
     &                                                     0.5D0*wgt1*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 
     &                                                     0.5D0*wgt2*(
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )                  
c
                  resp_Grhor(iz,ir,iso) = resp_Grhor(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )                   
c
                  resp_Grhop(iz,ir,iso) = resp_Grhop(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_Grhoz(iz,ir,iso) = resp_Grhoz(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
c                 2*resp_tau has to be added later
                  resp_Drho(iz,ir,iso) = resp_Drho(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*ddphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*ddphi(iz,ir,2,ist)
     &                 - ddphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + ddphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(
     &                 - drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist) 
     &                 + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
c
               enddo
               enddo
c
            endif
c
         enddo    ! end of loop over all transitions from state ist 
      enddo       ! end of loop over all states ist
c
c     add 2*resp_tau to get the full Drho response
c     ------------------------------------------------------------------
      do iz=0,nz
         do ir=0,nr
            resp_Drho(iz,ir,1) = resp_Drho(iz,ir,1) + 
     &                           2.0D0*resp_tau(iz,ir,1)
            resp_Drho(iz,ir,2) = resp_Drho(iz,ir,2) + 
     &                           2.0D0*resp_tau(iz,ir,2)
         enddo
      enddo
c
c     recouple the Skyrme-Parameters
c     ------------------------------------------------------------------
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
c
      wgt_b0  =  (b0-b0p)
      wgt_b0c =   b0
      wgt_b1  =  (b1-b1p)  
      wgt_b1c =   b1
      wgt_b2  = -(b2-b2p)
      wgt_b2c = - b2  
      if (b4p.eq.0.0D0) then    ! no isotopic spin-orbit force
         wgt_b4  = - 2.0D0*b4
         wgt_b4c = -       b4 
      else
         wgt_b4  = -(b4+b4p)
         wgt_b4c = - b4
      endif
      if(ifls.eq.1) then
         wgt_b5  =  -(b5-b5p)    ! b5 already includes the factor 2
         wgt_b5c =  - b5
      else
         wgt_b5  = 0.0D0
         wgt_b5c = 0.0D0
      endif
c
c     set the mean-fields to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
            RESP_U(iz,ir,iso) = 0.0D0           
            RESP_B(iz,ir,iso) = 0.0D0
            RESP_Wr(iz,ir,iso) = 0.0D0
            RESP_Wp(iz,ir,iso) = 0.0D0
            RESP_Wz(iz,ir,iso) = 0.0D0
c            RESP_drB(iz,ir,iso) = 0.0D0
c            RESP_dpB(iz,ir,iso) = 0.0D0
c            RESP_dzB(iz,ir,iso) = 0.0D0
         enddo   
         enddo
      enddo
c
c     construct the mean-fields from the densities
c     ------------------------------------------------------------------
  
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
c
c           the mean-field response without density-dependency (b3)
c           and coulomb
            RESP_U(iz,ir,iso)  = RESP_U(iz,ir,iso) 
     &                         + wgt_b0 *resp_rho(iz,ir,iso)
     &                         + wgt_b0c*resp_rho(iz,ir,3-iso)
     &                         + wgt_b1 *resp_tau(iz,ir,iso)
     &                         + wgt_b1c*resp_tau(iz,ir,3-iso)
     &                         + wgt_b2 *resp_Drho(iz,ir,iso)
     &                         + wgt_b2c*resp_Drho(iz,ir,3-iso)
     &                         + wgt_b4 *resp_divsoc(iz,ir,iso)
     &                         + wgt_b4c*resp_divsoc(iz,ir,3-iso)
c
            RESP_B(iz,ir,iso)  = RESP_B(iz,ir,iso) 
     &                         + wgt_b1 *resp_rho(iz,ir,iso)
     &                         + wgt_b1c*resp_rho(iz,ir,3-iso)
c
            RESP_Wr(iz,ir,iso) = RESP_Wr(iz,ir,iso) 
     &                         - wgt_b4 *resp_Grhor(iz,ir,iso)
     &                         - wgt_b4c*resp_Grhor(iz,ir,3-iso)
     &                         + wgt_b5 *resp_socr(iz,ir,iso)
     &                         + wgt_b5c*resp_socr(iz,ir,3-iso)        
c
            RESP_Wp(iz,ir,iso) = RESP_Wp(iz,ir,iso) 
     &                         - wgt_b4 *resp_Grhop(iz,ir,iso)
     &                         - wgt_b4c*resp_Grhop(iz,ir,3-iso)
     &                         + wgt_b5 *resp_socp(iz,ir,iso)
     &                         + wgt_b5c*resp_socp(iz,ir,3-iso)  
c
            RESP_Wz(iz,ir,iso) = RESP_Wz(iz,ir,iso) 
     &                         - wgt_b4 *resp_Grhoz(iz,ir,iso)
     &                         - wgt_b4c*resp_Grhoz(iz,ir,3-iso)
     &                         + wgt_b5 *resp_socz(iz,ir,iso)
     &                         + wgt_b5c*resp_socz(iz,ir,3-iso)  
c
c            RESP_drB(iz,ir,iso) = RESP_drB(iz,ir,iso) 
c     &                         + wgt_b1 *resp_Grhor(iz,ir,iso)
c     &                         + wgt_b1c*resp_Grhor(iz,ir,3-iso)
cc
c            RESP_dpB(iz,ir,iso) = RESP_dpB(iz,ir,iso) 
c     &                         + wgt_b1 *resp_Grhop(iz,ir,iso)
c     &                         + wgt_b1c*resp_Grhop(iz,ir,3-iso)
cc
c            RESP_dzB(iz,ir,iso) = RESP_dzB(iz,ir,iso) 
c     &                         + wgt_b1 *resp_Grhoz(iz,ir,iso)
c     &                         + wgt_b1c*resp_Grhoz(iz,ir,3-iso)
cc
         enddo   
         enddo
      enddo  
c
c     now the contribution from the density-dependency (b3)
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            RESP_U(iz,ir,1) = RESP_U(iz,ir,1) 
     &                      + locpotwgt1(iz,ir,1)*resp_rho(iz,ir,1)
     &                      + locpotwgt1(iz,ir,2)*resp_rho(iz,ir,2)
            RESP_U(iz,ir,2) = RESP_U(iz,ir,2) 
     &                      + locpotwgt2(iz,ir,1)*resp_rho(iz,ir,1)
     &                      + locpotwgt2(iz,ir,2)*resp_rho(iz,ir,2)
         enddo
      enddo
c
c     and from the coulomb
c     ------------------------------------------------------------------
      if (chg.ne.0.0D0) then
         
      endif
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'UPD_RESP_EVEN(): time-even response updated'//
     &  '  (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
c
      end


c=======================================================================
c
c     upd_response_odd(Q2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
c                      RESP_Ar,RESP_Ap,RESP_Az,
c                      RESP_Sr,RESP_Sp,RESP_Sz,
c                      RESP_Cr,RESP_Cp,RESP_Cz,il,imu)
c                       
c     updates the time-odd mean-field potentials (and densities)
c
c=======================================================================
      subroutine upd_response_Qtest(Q2q,
     &                            phi,drphi,dzphi,dpphi,ddphi,occphi,
     &                            RESP_Ar,RESP_Ap,RESP_Az,
     &                            RESP_Sr,RESP_Sp,RESP_Sz,
     &                            RESP_Cr,RESP_Cp,RESP_Cz,il,imu,
     &                            a1,a2,a3,a4,a5,a6)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision  Q2q(1:MTDTOT)
      double precision  phi(0:MZD,0:MRD,MSD,MPSID)
      double precision  drphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dzphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dpphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  ddphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  occphi(MPSID)
c
      double precision  RESP_Ar(0:MZD,0:MRD,2)
      double precision  RESP_Ap(0:MZD,0:MRD,2)
      double precision  RESP_Az(0:MZD,0:MRD,2)
      double precision  RESP_Sr(0:MZD,0:MRD,2)
      double precision  RESP_Sp(0:MZD,0:MRD,2)
      double precision  RESP_Sz(0:MZD,0:MRD,2)
      double precision  RESP_Cr(0:MZD,0:MRD,2)
      double precision  RESP_Cp(0:MZD,0:MRD,2)
      double precision  RESP_Cz(0:MZD,0:MRD,2)
      double precision  RESP_Delta(0:MZD,0:MRD,2)
c
      double precision  resp_spnr(0:MZD,0:MRD,2)
      double precision  resp_spnp(0:MZD,0:MRD,2)
      double precision  resp_spnz(0:MZD,0:MRD,2)
      double precision  resp_skdr(0:MZD,0:MRD,2)
      double precision  resp_skdp(0:MZD,0:MRD,2)
      double precision  resp_skdz(0:MZD,0:MRD,2)
      double precision  resp_curr(0:MZD,0:MRD,2)
      double precision  resp_curp(0:MZD,0:MRD,2)
      double precision  resp_curz(0:MZD,0:MRD,2)
      double precision  resp_paio(0:MZD,0:MRD,2)
      double precision  resp_rotjr(0:MZD,0:MRD,2)
      double precision  resp_rotjp(0:MZD,0:MRD,2)
      double precision  resp_rotjz(0:MZD,0:MRD,2)
      double precision  resp_rotsr(0:MZD,0:MRD,2)
      double precision  resp_rotsp(0:MZD,0:MRD,2)
      double precision  resp_rotsz(0:MZD,0:MRD,2)
c
      integer          ist,jst,jtp,msgn
      integer          itr,iz,ir
      integer          iso
      double precision eta,xi,wgt1,wgt2,wgt3,a1,a2,a3,a4,a5,a6
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     set response densities to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do ir=0,nr
            do iz=0,nz
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
            enddo
         enddo
      enddo
c
c     loop all single particle states
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif                  
c       
c        loop all transitions
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            jtp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0)) 
            uj   = sqrtb(one-vj*vj)
c
c           if the final state is from the original space
            if (jtp.eq.+1) then
c
               if (msgn.eq.0) then
                  wgt1=   0.0D0                           ! for   e^{i*mu*phi}   
                  wgt2=  -2.0D0*(uj*vi - vj*ui)*Q2q(itr2) ! for -ie^{i*mu*phi} 
                  wgt3=  +2.0D0*(uj*ui + vj*vi)*Q2q(itr2) 
               else if (msgn.gt.0) then
                  wgt1=   2.0D0*a1*(uj*vi - vj*ui)*Q2q(itr2) ! for   e^{i*mu*phi}
                  wgt2=   2.0D0*a2*(uj*vi - vj*ui)*Q2q(itr2) ! for -ie^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2)
               else 
                  wgt1=   2.0D0*a3*(uj*vi - vj*ui)*Q2q(itr2) !    e^{i*mu*phi}  
                  wgt2=   2.0D0*a4*(uj*vi - vj*ui)*Q2q(itr2) !  -ie^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2) 
               endif
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
c
c
                  resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) 
     &                                                  + 0.5D0*wgt1*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
c
                  resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) ) 
c
c
                  resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(
     &                 - dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(
     &                 + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist) )
c     
                  resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(
     &                 + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist) 
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist) 
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist) 
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
               enddo
               enddo
c
c           if the final state is from the time-reversed space
            else
c
               wgt1=   2.0D0*a5*(uj*vi - vj*ui)*Q2q(itr2) !    e^{i*mu*phi}  
               wgt2=   2.0D0*a6*(uj*vi - vj*ui)*Q2q(itr2) !  -ie^{i*mu*phi} 
               wgt3=   0.0D0*(uj*ui + vj*vi)*Q2q(itr2) 
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt2*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 - phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnz(iz,ir,iso) = resp_spnz(iz,ir,iso) + wgt1*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 - phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
c
c
                  resp_skdr(iz,ir,iso) = resp_skdr(iz,ir,iso) + wgt1*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt2*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
c
                  resp_skdz(iz,ir,iso) = resp_skdz(iz,ir,iso) + wgt1*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) 
     &                                                  + 0.5D0*wgt1*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) 
     &                                                  + 0.5D0*wgt2*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
c
                  resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt3*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) ) 
c
c
                  resp_rotjr(iz,ir,iso) = resp_rotjr(iz,ir,iso) + wgt1*(
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt2*(
     &                 - dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist) )
c     
                  resp_rotjz(iz,ir,iso) = resp_rotjz(iz,ir,iso) + wgt1*(
     &                 + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
c
                  resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist) 
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_rotsp(iz,ir,iso) = resp_rotsp(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist) 
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist) 
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
               enddo
               enddo
c
            endif
c
         enddo    ! end of loop over all transitions from state ist 
      enddo       ! end of loop over all states ist
c
c     recouple the Skyrme-Parameters
c     ------------------------------------------------------------------
      b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
      b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
      b4  = 0.500D0*       t4
      b5  = 0.125D0*(   x1*t1 +       x2*t2)  ! 2*b5
      b5p = 0.125D0*(      t1 -          t2)  ! 2*b5'
c
      wgt_b1  = -2.0D0*(b1-b1p)  
      wgt_b1c = -2.0D0* b1  
      if (b4p.eq.0.0D0) then    ! no isotopic spin-orbit force
         wgt_b4  = - 2.0D0*b4
         wgt_b4c = -       b4 
      else
         wgt_b4  = -(b4+b4p)
         wgt_b4c = - b4
      endif
      if(ifls.eq.1) then
         wgt_b5  =  (b5-b5p)    ! b5 already includes the factor 2
         wgt_b5c =   b5
      else
         wgt_b5  = 0.0D0
         wgt_b5c = 0.0D0
      endif
c
c     set the mean-fields to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
            RESP_Ar(iz,ir,iso) = 0.0D0           
            RESP_Ap(iz,ir,iso) = 0.0D0
            RESP_Az(iz,ir,iso) = 0.0D0
            RESP_Sr(iz,ir,iso) = 0.0D0
            RESP_Sp(iz,ir,iso) = 0.0D0
            RESP_Sz(iz,ir,iso) = 0.0D0
            RESP_Cr(iz,ir,iso) = 0.0D0
            RESP_Cp(iz,ir,iso) = 0.0D0
            RESP_Cz(iz,ir,iso) = 0.0D0
         enddo   
         enddo
      enddo
c
c     construct the mean-fields from the densities
c     ------------------------------------------------------------------  
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
c
            RESP_Ar(iz,ir,iso) = RESP_Ar(iz,ir,iso) 
     &                         + wgt_b1 *resp_curr(iz,ir,iso)
     &                         + wgt_b1c*resp_curr(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotsr(iz,ir,iso)
     &                         + wgt_b4c*resp_rotsr(iz,ir,3-iso)
            RESP_Ap(iz,ir,iso) = RESP_Ap(iz,ir,iso) 
     &                         + wgt_b1 *resp_curp(iz,ir,iso)
     &                         + wgt_b1c*resp_curp(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotsp(iz,ir,iso)
     &                         + wgt_b4c*resp_rotsp(iz,ir,3-iso)
            RESP_Az(iz,ir,iso) = RESP_Az(iz,ir,iso) 
     &                         + wgt_b1 *resp_curz(iz,ir,iso)
     &                         + wgt_b1c*resp_curz(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotsz(iz,ir,iso)
     &                         + wgt_b4c*resp_rotsz(iz,ir,3-iso)
c
            RESP_Sr(iz,ir,iso) = RESP_Sr(iz,ir,iso) 
     &                         + wgt_b5 *resp_skdr(iz,ir,iso)
     &                         + wgt_b5c*resp_skdr(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotjr(iz,ir,iso)
     &                         + wgt_b4c*resp_rotjr(iz,ir,3-iso)
            RESP_Sp(iz,ir,iso) = RESP_Sp(iz,ir,iso) 
     &                         + wgt_b5 *resp_skdp(iz,ir,iso)
     &                         + wgt_b5c*resp_skdp(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotjp(iz,ir,iso)
     &                         + wgt_b4c*resp_rotjp(iz,ir,3-iso)
            RESP_Sz(iz,ir,iso) = RESP_Sz(iz,ir,iso) 
     &                         + wgt_b5 *resp_skdz(iz,ir,iso)
     &                         + wgt_b5c*resp_skdz(iz,ir,3-iso)
     &                         + wgt_b4 *resp_rotjz(iz,ir,iso)
     &                         + wgt_b4c*resp_rotjz(iz,ir,3-iso)
c
            RESP_Cr(iz,ir,iso) = RESP_Cr(iz,ir,iso) 
     &                         + wgt_b5 *resp_spnr(iz,ir,iso)
     &                         + wgt_b5c*resp_spnr(iz,ir,3-iso)
            RESP_Cp(iz,ir,iso) = RESP_Cp(iz,ir,iso) 
     &                         + wgt_b5 *resp_spnp(iz,ir,iso)
     &                         + wgt_b5c*resp_spnp(iz,ir,3-iso)
            RESP_Cz(iz,ir,iso) = RESP_Cz(iz,ir,iso) 
     &                         + wgt_b5 *resp_spnz(iz,ir,iso)
     &                         + wgt_b5c*resp_spnz(iz,ir,3-iso)
c
         enddo   
         enddo
      enddo  
c
c     ALL FIELDS NOW CONTAIN THE (r,z)-PART OF THE RESPONSE FIELDS
c     THE ANGULAR DEPENDENCY READS
c
c     Cr,Cz,Sr,Sz,Ap  ~  sin(mu*phi)*X(r,z)
c     Cp,Sp,Ar,Az     ~  cos(mu*phi)*X(r,z)
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'UPD_RESP_ODD(): time-odd response updated'//
     &  '  (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
c
      end


c=======================================================================
c
c     upd_response_even(P2q,phi,drphi,dzphi,dpphi,ddphi,occphi,
c                       RESP_U,RESP_B,RESP_drB,RESP_dzB,RESP_Wz,RESP_Wr,
c                       
c
c     count the number of transitions (matrixelements) for a given
c     cranking mode
c
c=======================================================================
      subroutine upd_response_even(P2q,
     &                             phi,drphi,dzphi,dpphi,ddphi,occphi,
     &                             RESP_U,
     &                             RESP_B,RESP_drB,RESP_dpB,RESP_dzB,
     &                             RESP_Wr,RESP_Wp,RESP_Wz,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision  P2q(1:MTDTOT)
      double precision  phi(0:MZD,0:MRD,MSD,MPSID)
      double precision  drphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dzphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dpphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  ddphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  occphi(MPSID)
c
      double precision  RESP_U(0:MZD,0:MRD,2)
      double precision  RESP_B(0:MZD,0:MRD,2)
      double precision  RESP_drB(0:MZD,0:MRD,2)
      double precision  RESP_dpB(0:MZD,0:MRD,2)
      double precision  RESP_dzB(0:MZD,0:MRD,2)
      double precision  RESP_Wr(0:MZD,0:MRD,2)
      double precision  RESP_Wp(0:MZD,0:MRD,2)
      double precision  RESP_Wz(0:MZD,0:MRD,2)
      double precision  RESP_Delta(0:MZD,0:MRD,2)
c
      double precision  resp_rho(0:MZD,0:MRD,2)
      double precision  resp_tau(0:MZD,0:MRD,2)
      double precision  resp_socr(0:MZD,0:MRD,2)
      double precision  resp_socp(0:MZD,0:MRD,2)
      double precision  resp_socz(0:MZD,0:MRD,2)
      double precision  resp_Grhor(0:MZD,0:MRD,2)
      double precision  resp_Grhop(0:MZD,0:MRD,2)
      double precision  resp_Grhoz(0:MZD,0:MRD,2)
      double precision  resp_Drho(0:MZD,0:MRD,2)
      double precision  resp_divsoc(0:MZD,0:MRD,2)
      double precision  resp_pair(0:MZD,0:MRD,2)
c
      integer          ist,jst,jtp
      integer          itr,iz,ir
      integer          iso
      double precision eta,xi,wgt1,wgt2,wgt3
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     set response densities to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do iz=0,nz
            do ir=0,nr
               resp_rho(iz,ir,iso)    = 0.0D0
               resp_tau(iz,ir,iso)    = 0.0D0
               resp_socr(iz,ir,iso)   = 0.0D0
               resp_socp(iz,ir,iso)   = 0.0D0
               resp_socz(iz,ir,iso)   = 0.0D0
               resp_Grhor(iz,ir,iso)  = 0.0D0
               resp_Grhop(iz,ir,iso)  = 0.0D0
               resp_Grhoz(iz,ir,iso)  = 0.0D0
               resp_Drho(iz,ir,iso)   = 0.0D0
               resp_divsoc(iz,ir,iso) = 0.0D0
               resp_pair(iz,ir,iso)   = 0.0D0
            enddo
         enddo
      enddo
c
c     loop all single particle states
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif                  
c       
c        loop all transitions
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            jtp  = TRidx(ist,itr,2)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)
c
c           if the final state is from the original space
            if (jtp.eq.+1) then
c
               if (msgn.eq.0) then
                  wgt1=   0.0D0                           ! for -ie^{i*mu*phi}   
                  wgt2=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2) 
               else if (msgn.gt.0) then
                  wgt1=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for -ie^{i*mu*phi}
                  wgt2=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2)
               else 
                  wgt1=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for -ie^{i*mu*phi}  
                  wgt2=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
                  wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2) 
               endif
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )
c
                  resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 
     &                                                   0.5D0*wgt2*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 
     &                                                 0.5D0*wgt1*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 
     &                                                     0.5D0*wgt2*(
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) )                  
c
                  resp_Grhor(iz,ir,iso) = resp_Grhor(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )                   
c
                  resp_Grhop(iz,ir,iso) = resp_Grhop(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_Grhoz(iz,ir,iso) = resp_Grhoz(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
c                 2*resp_tau has to be added later
                  resp_Drho(iz,ir,iso) = resp_Drho(iz,ir,iso) + wgt2*(
     &                 +   phi(iz,ir,1,jst)*ddphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*ddphi(iz,ir,2,ist)
     &                 + ddphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + ddphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(
     &                 + drphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist) 
     &                 - dzphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) )
c
               enddo
               enddo
c
c           if the final state is from the time-reversed space
            else
c
               wgt1=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for -ie^{i*mu*phi}  
               wgt2=   2.0D0*(uj*vi + vj*ui)*P2q(itr2) ! for   e^{i*mu*phi} 
               wgt3=   0.0D0*(uj*ui - vj*vi)*P2q(itr2) 
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_rho(iz,ir,iso) = resp_rho(iz,ir,iso) + wgt2*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )
c
                  resp_tau(iz,ir,iso) = resp_tau(iz,ir,iso) + wgt2*(
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) )
c
                  resp_socr(iz,ir,iso) = resp_socr(iz,ir,iso) + 
     &                                                     0.5D0*wgt2*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socp(iz,ir,iso) = resp_socp(iz,ir,iso) + 
     &                                                     0.5D0*wgt1*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_socz(iz,ir,iso) = resp_socz(iz,ir,iso) + 
     &                                                     0.5D0*wgt2*(
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_pair(iz,ir,iso) = resp_pair(iz,ir,iso) + wgt3*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) )                  
c
                  resp_Grhor(iz,ir,iso) = resp_Grhor(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )                   
c
                  resp_Grhop(iz,ir,iso) = resp_Grhop(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_Grhoz(iz,ir,iso) = resp_Grhoz(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
c                 2*resp_tau has to be added later
                  resp_Drho(iz,ir,iso) = resp_Drho(iz,ir,iso) + wgt2*(
     &                 -   phi(iz,ir,2,jst)*ddphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*ddphi(iz,ir,2,ist)
     &                 - ddphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + ddphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_divsoc(iz,ir,iso)=resp_divsoc(iz,ir,iso)+wgt2*(
     &                 - drphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist) 
     &                 + dzphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) )
c
               enddo
               enddo
c
            endif
c
         enddo    ! end of loop over all transitions from state ist 
      enddo       ! end of loop over all states ist
c
c     add 2*resp_tau to get the full Drho response
c     ------------------------------------------------------------------
      do iz=0,nz
         do ir=0,nr
            resp_Drho(iz,ir,1) = resp_Drho(iz,ir,1) + 
     &                           2.0D0*resp_tau(iz,ir,1)
            resp_Drho(iz,ir,2) = resp_Drho(iz,ir,2) + 
     &                           2.0D0*resp_tau(iz,ir,2)
         enddo
      enddo
c
c     recouple the Skyrme-Parameters
c     ------------------------------------------------------------------
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
c
      wgt_b0  =  (b0-b0p)
      wgt_b0c =   b0
      wgt_b1  =  (b1-b1p)  
      wgt_b1c =   b1
      wgt_b2  = -(b2-b2p)
      wgt_b2c = - b2  
      if (b4p.eq.0.0D0) then    ! no isotopic spin-orbit force
         wgt_b4  = - 2.0D0*b4
         wgt_b4c = -       b4 
      else
         wgt_b4  = -(b4+b4p)
         wgt_b4c = - b4
      endif
      if(ifls.eq.1) then
         wgt_b5  =  -(b5-b5p)    ! b5 already includes the factor 2
         wgt_b5c =  - b5
      else
         wgt_b5  = 0.0D0
         wgt_b5c = 0.0D0
      endif
c
c     set the mean-fields to zero
c     ------------------------------------------------------------------
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
            RESP_U(iz,ir,iso) = 0.0D0           
            RESP_B(iz,ir,iso) = 0.0D0
            RESP_Wr(iz,ir,iso) = 0.0D0
            RESP_Wp(iz,ir,iso) = 0.0D0
            RESP_Wz(iz,ir,iso) = 0.0D0
c            RESP_drB(iz,ir,iso) = 0.0D0
c            RESP_dpB(iz,ir,iso) = 0.0D0
c            RESP_dzB(iz,ir,iso) = 0.0D0
         enddo   
         enddo
      enddo
c
c     construct the mean-fields from the densities
c     ------------------------------------------------------------------
  
      do iso=1,2
         do ir=0,nr
         do iz=0,nz
c
c           the mean-field response without density-dependency (b3)
c           and coulomb
            RESP_U(iz,ir,iso)  = RESP_U(iz,ir,iso) 
     &                         + wgt_b0 *resp_rho(iz,ir,iso)
     &                         + wgt_b0c*resp_rho(iz,ir,3-iso)
     &                         + wgt_b1 *resp_tau(iz,ir,iso)
     &                         + wgt_b1c*resp_tau(iz,ir,3-iso)
     &                         + wgt_b2 *resp_Drho(iz,ir,iso)
     &                         + wgt_b2c*resp_Drho(iz,ir,3-iso)
     &                         + wgt_b4 *resp_divsoc(iz,ir,iso)
     &                         + wgt_b4c*resp_divsoc(iz,ir,3-iso)
c
            RESP_B(iz,ir,iso)  = RESP_B(iz,ir,iso) 
     &                         + wgt_b1 *resp_rho(iz,ir,iso)
     &                         + wgt_b1c*resp_rho(iz,ir,3-iso)
c
            RESP_Wr(iz,ir,iso) = RESP_Wr(iz,ir,iso) 
     &                         - wgt_b4 *resp_Grhor(iz,ir,iso)
     &                         - wgt_b4c*resp_Grhor(iz,ir,3-iso)
     &                         + wgt_b5 *resp_socr(iz,ir,iso)
     &                         + wgt_b5c*resp_socr(iz,ir,3-iso)        
c
            RESP_Wp(iz,ir,iso) = RESP_Wp(iz,ir,iso) 
     &                         - wgt_b4 *resp_Grhop(iz,ir,iso)
     &                         - wgt_b4c*resp_Grhop(iz,ir,3-iso)
     &                         + wgt_b5 *resp_socp(iz,ir,iso)
     &                         + wgt_b5c*resp_socp(iz,ir,3-iso)  
c
            RESP_Wz(iz,ir,iso) = RESP_Wz(iz,ir,iso) 
     &                         - wgt_b4 *resp_Grhoz(iz,ir,iso)
     &                         - wgt_b4c*resp_Grhoz(iz,ir,3-iso)
     &                         + wgt_b5 *resp_socz(iz,ir,iso)
     &                         + wgt_b5c*resp_socz(iz,ir,3-iso)  
c
c            RESP_drB(iz,ir,iso) = RESP_drB(iz,ir,iso) 
c     &                         + wgt_b1 *resp_Grhor(iz,ir,iso)
c     &                         + wgt_b1c*resp_Grhor(iz,ir,3-iso)
cc
c            RESP_dpB(iz,ir,iso) = RESP_dpB(iz,ir,iso) 
c     &                         + wgt_b1 *resp_Grhop(iz,ir,iso)
c     &                         + wgt_b1c*resp_Grhop(iz,ir,3-iso)
cc
c            RESP_dzB(iz,ir,iso) = RESP_dzB(iz,ir,iso) 
c     &                         + wgt_b1 *resp_Grhoz(iz,ir,iso)
c     &                         + wgt_b1c*resp_Grhoz(iz,ir,3-iso)
cc
         enddo   
         enddo
      enddo  
c
c     now the contribution from the density-dependency (b3)
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            RESP_U(iz,ir,1) = RESP_U(iz,ir,1) 
     &                      + locpotwgt1(iz,ir,1)*resp_rho(iz,ir,1)
     &                      + locpotwgt1(iz,ir,2)*resp_rho(iz,ir,2)
            RESP_U(iz,ir,2) = RESP_U(iz,ir,2) 
     &                      + locpotwgt2(iz,ir,1)*resp_rho(iz,ir,1)
     &                      + locpotwgt2(iz,ir,2)*resp_rho(iz,ir,2)
         enddo
      enddo
c
c     and from the coulomb
c     ------------------------------------------------------------------
      if (chg.ne.0.0D0) then
         
      endif
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a,i2,a,i2,a)') 
     &  'UPD_RESP_EVEN(): time-even response updated'//
     &  '  (lambda=',il,' mu=',imu,')'
      endif
c
c     ------------------------------------------------------------------
c
      return
c
      end


c=======================================================================
c
c     upd_Hresp_even(phi,drphi,dzphi,dpphi,ddphi,occphi,
c                   RESP_U,
c                   RESP_B,RESP_drB,RESP_dpB,RESP_dzB,
c                   RESP_Wr,RESP_Wp,RESP_Wz,il,imu)
c                       
c     evaluate the time-even response hamiltonian
c
c=======================================================================
      subroutine upd_Hresp_even(Hresp2q,
     &                          phi,drphi,dzphi,dpphi,ddphi,occphi,
     &                          RESP_U,RESP_B,
     &                          RESP_Wr,RESP_Wp,RESP_Wz,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision  Hresp2q(1:MTDTOT)
      double precision  phi(0:MZD,0:MRD,MSD,MPSID)
      double precision  drphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dzphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dpphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  ddphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  occphi(MPSID)
c
      double precision  RESP_U(0:MZD,0:MRD,2)
      double precision  RESP_B(0:MZD,0:MRD,2)
      double precision  RESP_Wr(0:MZD,0:MRD,2)
      double precision  RESP_Wp(0:MZD,0:MRD,2)
      double precision  RESP_Wz(0:MZD,0:MRD,2)
      double precision  RESP_Delta(0:MZD,0:MRD,2)
c
      double precision  resp_rho(0:MZD,0:MRD)
      double precision  resp_tau(0:MZD,0:MRD)
      double precision  resp_socr(0:MZD,0:MRD)
      double precision  resp_socp(0:MZD,0:MRD)
      double precision  resp_socz(0:MZD,0:MRD)
      double precision  resp_pair(0:MZD,0:MRD)
c
      integer          ist,jst,jtp
      integer          itr,iz,ir
      integer          iso
      double precision eta,xi,wgt1,wgt2
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     loop all single particle states
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif                  
c       
c        loop all transitions
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            jtp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)     
            xi   = uj*ui - vj*vi
            eta  = uj*vi + vj*ui   
c
c           if the final state is from the original space
            if (jtp.eq.+1) then
c
               do ir=0,nr
               do iz=0,nz
c
                  resp_rho(iz,ir) =                ! ~    exp(i*mu*phi)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) 
c
                  resp_tau(iz,ir) =                ! ~    exp(i*mu*phi)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) 
c
                  resp_socr(iz,ir) = 0.5D0*(        ! ~    exp(i*mu*phi)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_socp(iz,ir) = 0.5D0*(       ! ~  -i*exp(i*mu*phi)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_socz(iz,ir) = 0.5D0*(       ! ~    exp(i*mu*phi)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_pair(iz,ir) =                ! ~    exp(i*mu*phi)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist)       
c
               enddo
               enddo
c
c           if the final state is from the time-reversed space
            else
c
               do ir=0,nr
               do iz=0,nz
c
                  resp_rho(iz,ir) =                ! ~    exp(i*mu*phi)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
c
                  resp_tau(iz,ir) =                ! ~    exp(i*mu*phi)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) 
c
                  resp_socr(iz,ir) = 0.5D0*(      ! ~    exp(i*mu*phi)
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_socp(iz,ir) = 0.5D0*(      ! ~  i*exp(i*mu*phi)
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_socz(iz,ir) = 0.5D0*(      ! ~    exp(i*mu*phi)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 -   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_pair(iz,ir) =                ! ~    exp(i*mu*phi)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)            
c
               enddo
               enddo
c
            endif
c       
c          evaluate the matrixelement of the transition
c          -------------------------------------------------------------
            if (imu.eq.0) then
               Hresp2q(itr2) = eta*(
     &              + f2int(RESP_U(0,0,iso) ,resp_rho)
     &              + f2int(RESP_B(0,0,iso) ,resp_tau)
     &              + f2int(RESP_Wr(0,0,iso),resp_socr)
     &              + f2int(RESP_Wz(0,0,iso),resp_socz)
     &                             )
            else 
               if (msgn.gt.0) then
                  Hresp2q(itr2) = 0.5D0*eta*(
     &              + f2int(RESP_U(0,0,iso) ,resp_rho)
     &              + f2int(RESP_B(0,0,iso) ,resp_tau)
     &              + f2int(RESP_Wr(0,0,iso),resp_socr)
     &              + f2int(RESP_Wz(0,0,iso),resp_socz)
     &              - f2int(RESP_Wp(0,0,iso),resp_socp)
     &                                      )
               else
                  Hresp2q(itr2) =  0.5D0*eta*(
     &              + f2int(RESP_U(0,0,iso) ,resp_rho)
     &              + f2int(RESP_B(0,0,iso) ,resp_tau)
     &              + f2int(RESP_Wr(0,0,iso),resp_socr)
     &              + f2int(RESP_Wz(0,0,iso),resp_socz)
     &              + f2int(RESP_Wp(0,0,iso),resp_socp)
     &                                       )
               endif
            endif
c
         enddo    ! end of loop over all transitions from state ist 
      enddo       ! end of loop over all states ist
c
      return
      end


c=======================================================================
c
c     upd_Hresp_odd(phi,drphi,dzphi,dpphi,ddphi,occphi,
c                   RESP_Ar,RESP_Ap,RESP_Az,
c                   RESP_Sr,RESP_Sp,RESP_Sz,
c                   RESP_Cr,RESP_Cp,RESP_Cz,il,imu)
c                       
c     evaluate the time-odd response hamiltonian
c
c=======================================================================
      subroutine upd_Hresp_odd(Hresp2q,
     &                   phi,drphi,dzphi,dpphi,ddphi,occphi,
     &                   RESP_Ar,RESP_Ap,RESP_Az,
     &                   RESP_Sr,RESP_Sp,RESP_Sz,
     &                   RESP_Cr,RESP_Cp,RESP_Cz,il,imu)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
c     ------------------------------------------------------------------
c
      integer          il                     ! angular momentum
      integer          imu                    ! delta K
c
      double precision  Hresp2q(1:MTDTOT)
      double precision  phi(0:MZD,0:MRD,MSD,MPSID)
      double precision  drphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dzphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  dpphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  ddphi(0:MZD,0:MRD,MSD,MPSID)
      double precision  occphi(MPSID)
c
      double precision  RESP_Ar(0:MZD,0:MRD,2)
      double precision  RESP_Ap(0:MZD,0:MRD,2)
      double precision  RESP_Az(0:MZD,0:MRD,2)
      double precision  RESP_Sr(0:MZD,0:MRD,2)
      double precision  RESP_Sp(0:MZD,0:MRD,2)
      double precision  RESP_Sz(0:MZD,0:MRD,2)
      double precision  RESP_Cr(0:MZD,0:MRD,2)
      double precision  RESP_Cp(0:MZD,0:MRD,2)
      double precision  RESP_Cz(0:MZD,0:MRD,2)
      double precision  RESP_Delta(0:MZD,0:MRD,2)
c
      double precision  resp_spnr(0:MZD,0:MRD)
      double precision  resp_spnp(0:MZD,0:MRD)
      double precision  resp_spnz(0:MZD,0:MRD)
      double precision  resp_skdr(0:MZD,0:MRD)
      double precision  resp_skdp(0:MZD,0:MRD)
      double precision  resp_skdz(0:MZD,0:MRD)
      double precision  resp_curr(0:MZD,0:MRD)
      double precision  resp_curp(0:MZD,0:MRD)
      double precision  resp_curz(0:MZD,0:MRD)
      double precision  resp_paio(0:MZD,0:MRD)
      double precision  resp_rotjr(0:MZD,0:MRD)
      double precision  resp_rotjp(0:MZD,0:MRD)
      double precision  resp_rotjz(0:MZD,0:MRD)
      double precision  resp_rotsr(0:MZD,0:MRD)
      double precision  resp_rotsp(0:MZD,0:MRD)
      double precision  resp_rotsz(0:MZD,0:MRD)
c
      integer          ist,jst,jtp
      integer          itr,iz,ir
      integer          iso,msgn
      double precision eta,xi,wgt1,wgt2
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     ------------------------------------------------------------------
c
c     check parameters
c     ------------------------------------------------------------------
      if (imu.lt.0) then
         write(21,'(a,i2)') 
     &            'ERROR: mu-component must be >=0 -- not ',imu
         stop 'in linresp'
      else if (imu.gt.il) then
         write(21,'(a,i2,a,i2)') 
     &            'ERROR: mu-component must be <=il (=',il,
     &            ')  -- not ',imu
         stop 'in linresp'
      else if (il.gt.LAMBDAMAX) then
         write(21,'(a,i2,a,i2,a)') 
     &            'ERROR: LAMBDAMAX too small (il=',il,
     &            ' LAMBDAMAX=',il,')'
         stop 'in linresp'
      endif
c
c     loop all single particle states
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (ist.le.nst(1)) then
            iso    = 1
         else
            iso    = 2
         endif                  
c       
c        loop all transitions
c        ---------------------------------------------------------------
         do itr=1,TRntr(ist)
c
            jst  = TRidx(ist,itr,1)
            jtp  = TRidx(ist,itr,2)
            msgn = TRidx(ist,itr,3)
            itr2 = TRidx(ist,itr,4)
c
            vi   = sqrtb(CPvsqf(ist,0))  
            ui   = sqrtb(one-vi*vi)
            vj   = sqrtb(CPvsqf(jst,0))
            uj   = sqrtb(one-vj*vj)     
c
            xi   = uj*ui + vj*vi
            eta  = uj*vi - vj*ui   
c
c           if the final state is from the original space
            if (jtp.eq.+1) then
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir) =               ! ~    exp(i*mu*phi)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,1,ist) 
c
                  resp_spnp(iz,ir) =               ! ~ -i*exp(i*mu*phi)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist) 
c
                  resp_spnz(iz,ir) =               ! ~    exp(i*mu*phi) 
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist) 
c
c
                  resp_skdr(iz,ir) =               ! ~    exp(i*mu*phi)  
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) 
c
                  resp_skdp(iz,ir) =               ! ~ -i*exp(i*mu*phi) 
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist) 
c
                  resp_skdz(iz,ir) =               ! ~    exp(i*mu*phi) 
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist) 
c
c
                  resp_curr(iz,ir) = 0.5D0*(              ! ~ -i*exp(i*mu*phi) 
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_curp(iz,ir) = 0.5D0*(              ! ~    exp(i*mu*phi)  
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)) 
c     
                  resp_curz(iz,ir) = 0.5D0*(              ! ~ -i*exp(i*mu*phi)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist)) 
c
c
                  resp_paio(iz,ir) =               ! ~  i*exp(i*mu*phi)  
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,2,ist) 
c
               enddo
               enddo
c
c           if the final state is from the time-reversed space
            else
c
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir) =               ! ~    exp(i*mu*phi) 
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist) 
c
                  resp_spnp(iz,ir) =               ! ~ -i*exp(i*mu*phi) 
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 - phi(iz,ir,1,jst)*phi(iz,ir,1,ist) 
c
                  resp_spnz(iz,ir) =               ! ~    exp(i*mu*phi)  
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 - phi(iz,ir,1,jst)*phi(iz,ir,2,ist) 
c
c
                  resp_skdr(iz,ir) =               ! ~    exp(i*mu*phi) 
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 - dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) 
c
                  resp_skdp(iz,ir) =               ! ~ -i*exp(i*mu*phi) 
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,1,ist) 
c
                  resp_skdz(iz,ir) =               ! ~    exp(i*mu*phi)  
     &                 - drphi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,1,jst)*dpphi(iz,ir,2,ist) 
c
c
                  resp_curr(iz,ir) = 0.5D0*(              ! ~ -i*exp(i*mu*phi)  
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)) 
c
                  resp_curp(iz,ir) = 0.5D0*(              ! ~    exp(i*mu*phi)  
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)) 
c     
                  resp_curz(iz,ir) = 0.5D0*(              ! ~ -i*exp(i*mu*phi) 
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist)) 
c
c
                  resp_paio(iz,ir) =               ! ~  i*exp(i*mu*phi)  
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,1,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist) 
c
               enddo
               enddo           
c
            endif
c
c          -------------------------------------------------------------
            if (imu.eq.0) then
               Hresp2q(itr2) = eta*( 
     &              - f2int(RESP_Ar(0,0,iso),resp_curr)
     &              - f2int(RESP_Az(0,0,iso),resp_curz)
     &              - f2int(RESP_Sp(0,0,iso),resp_spnp)
     &              - f2int(RESP_Cp(0,0,iso),resp_skdp)
     &                             )
            else 
               if (msgn.gt.0) then          
                  Hresp2q(itr2) = 0.5D0*eta*(
     &              - f2int(RESP_Ar(0,0,iso),resp_curr)
     &              - f2int(RESP_Az(0,0,iso),resp_curz)
     &              - f2int(RESP_Sp(0,0,iso),resp_spnp)
     &              - f2int(RESP_Cp(0,0,iso),resp_skdp)
     &              + f2int(RESP_Ap(0,0,iso),resp_curp)
     &              + f2int(RESP_Sr(0,0,iso),resp_spnr)
     &              + f2int(RESP_Sz(0,0,iso),resp_spnz)
     &              + f2int(RESP_Cr(0,0,iso),resp_skdr)
     &              + f2int(RESP_Cz(0,0,iso),resp_skdz)
     &                 )
               else   
                  Hresp2q(itr2) = 0.5D0*eta*(
     &              - f2int(RESP_Ar(0,0,iso),resp_curr)
     &              - f2int(RESP_Az(0,0,iso),resp_curz)
     &              - f2int(RESP_Sp(0,0,iso),resp_spnp)
     &              - f2int(RESP_Cp(0,0,iso),resp_skdp)
     &              - f2int(RESP_Ap(0,0,iso),resp_curp)
     &              - f2int(RESP_Sr(0,0,iso),resp_spnr)
     &              - f2int(RESP_Sz(0,0,iso),resp_spnz)
     &              - f2int(RESP_Cr(0,0,iso),resp_skdr)
     &              - f2int(RESP_Cz(0,0,iso),resp_skdz)
     &                                )
               endif
            endif
c
         enddo    ! end of loop over all transitions from state ist 
      enddo       ! end of loop over all states ist
c
      return
      end




c=======================================================================
c
c     upd_locpotwgt()
c     
c       the response of the non linear contributions to the local 
c       mean-field potential U is linearized
c
c              U^{NL}_q = U^{NL}_{q,0} + 
c                         U^{NL}_{q,resp,p} rerho_p +
c                         U^{NL}_{q,resp,n} rerho_n
c
c       with
c
c              U^{NL}_{q,resp,q'} = dU_q/drho_q'
c
c       this function calculates the response wheight dU_q/drho_q' and
c       stores its values in locpotwgt:
c       
c          locpotwgt1(0,0,1) = dU_p/drho_p
c          locpotwgt1(0,0,2) = dU_p/drho_n
c          locpotwgt2(0,0,1) = dU_n/drho_p
c          locpotwgt2(0,0,2) = dU_n/drho_n
c
c=======================================================================
      subroutine upd_locpotwgt()
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'
c
c     ------------------------------------------------------------------
c
      double precision b3,b3p,fac1,fac2,fac3,fac4,fac5
      double precision buff1,buff2,buff3,buff4,buff5,buff6,buff7,buff8
      double precision wgtfac
      integer          ir,iz,ist
c
c     ------------------------------------------------------------------
c
c     calculate densities rho             --> LRwrk1_spinor
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            LRwrk1_spinor(iz,ir,1) = 0.0D0
            LRwrk1_spinor(iz,ir,2) = 0.0D0
         enddo
      enddo
c
      do ist=1,nst(1)
         wgtfac = 2.0D0*CPvsqf(ist,0)  
         do ir=0,nr
            do iz=0,nz
               LRwrk1_spinor(iz,ir,1) = LRwrk1_spinor(iz,ir,1)+
     &            wgtfac*(CPpsi(iz,ir,1,ist, 0)*CPpsi(iz,ir,1,ist, 0) +
     &                    CPpsi(iz,ir,2,ist, 0)*CPpsi(iz,ir,2,ist, 0))
            enddo
         enddo
      enddo
c
      do ist=nst(1)+1,nst(2)
         wgtfac = 2.0D0*CPvsqf(ist,0) 
         do ir=0,nr
            do iz=0,nz
               LRwrk1_spinor(iz,ir,2) = LRwrk1_spinor(iz,ir,2)+
     &            wgtfac*(CPpsi(iz,ir,1,ist, 0)*CPpsi(iz,ir,1,ist, 0) +
     &                    CPpsi(iz,ir,2,ist, 0)*CPpsi(iz,ir,2,ist, 0))
            enddo
         enddo
      enddo
c      
c     calculate response weight dU/drho   --> locpotwgt
c     ------------------------------------------------------------------
      b3   = 0.250D0*(      t3 + 0.5D0*t3*x3                         )
      b3p  = 0.250D0*(0.5D0*t3 +       t3*x3                         )
      fac1 =  b3 *(power+2.0D0)*(power+1.0D0)/3.0D0
      fac2 = -b3p*(2.0D0/3.0D0)
      fac3 = -b3p*(2.0D0*power)/3.0D0
      fac4 = -b3p*power*(power-1.0D0)/3.0D0

      if(chgx.ne.zero) then
         fac5 = 1.418033D0*chgx/3.0D0
      else
         fac5 = 0.0D0
      endif

      do ir=0,nr
         do iz=0,nz
            locpotwgt1(iz,ir,1)=0.0D0
            locpotwgt1(iz,ir,2)=0.0D0
            locpotwgt2(iz,ir,1)=0.0D0
            locpotwgt2(iz,ir,2)=0.0D0
         enddo
      enddo
c
      do ir=0,nr
         do iz=0,nz
            buff1 = LRwrk1_spinor(iz,ir,1) ! rho_p
            buff2 = LRwrk1_spinor(iz,ir,2) ! rho_n
            buff3 = buff1 + buff2          ! rho
            if (buff3.le.smal) then
               buff4 = 0.0D0               ! rho**alpha
               buff5 = 1.0D0               ! (rho_p^2 + rho_n^2)/rho^2
               buff6 = 1.0D0               ! rho_p/rho
               buff7 = 1.0D0               ! rho_n/rho
               buff8 = 0.0D0
            else
               buff4 =  buff3**(power)   
               buff5 = (buff1**2 + buff2**2)/buff3**2
               buff6 = (buff1/buff3)               
               buff7 = (buff2/buff3)
               buff8 = 0.0D0
               buff8 = buff1**(-0.6666666666666667D0)
            endif
               
c
c           dU_p/drho_p
            locpotwgt1(iz,ir,1)=
     &           buff4*(fac1 + fac2 +  2.0D0*fac3        *buff6 
     &                              +        fac4        *buff5)
     &          +buff8*fac5
c     
c           dU_n/drho_n
            locpotwgt2(iz,ir,2)=
     &           buff4*(fac1 + fac2 +  2.0D0*fac3        *buff7   
     &                              +        fac4        *buff5)
c
c           dU_p/drho_n
            locpotwgt1(iz,ir,2)=
     &           buff4*(fac1        +       fac3          
     &                              +       fac4         *buff5)
c
c           dU_n/drho_p
            locpotwgt2(iz,ir,1)=
     &           buff4*(fac1        +       fac3          
     &                              +       fac4         *buff5)

         enddo
      enddo
c
c     ------------------------------------------------------------------
c
      if (D_LRDEBUG) then
         write(21,'(a)')         
     &  'UPD_LOCPOTWGT() : response weight dU/drho updated'
      endif
c
c     ------------------------------------------------------------------
c
      return
      end



c=======================================================================
c
c     ini_files()
c
c=======================================================================
      subroutine ini_files(outnam)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'constr.inc'
      include 'linresp.inc'
      include 'levels.inc'
c
      integer          i
      character*3      outnam
!      character*78     title
c
c     ------------------------------------------------------------------
c
      open(21,status='unknown',form='formatted',
     &        file='linresp.'//outnam)
      open(22,status='unknown',form='formatted',
     &        file='LRiter.'//outnam)
      open(29,status='unknown',form='formatted',
     &        file='radii.'//outnam)      
      open(30,status='unknown',form='formatted',      
     &        file='width.'//outnam)
      open(32,status='unknown',form='formatted',
     &        file='norm.'//outnam)

      write(29,'(a/a,2i4/a,i1,3a)')
     &   title(1:78),
     &   ' A,Z=',npart(1)+npart(2),npart(1),
     &   ' beta_',lqcon,'   <P^2_beta>',
     &   '     r             R             sigma  ',
     &   '     (d/db)^2 r    (d/db)^2 R    (d/db)^2 sigma'
c
c     Kopfzeile der beiben Dateien norm und width vertauscht 30<>32
c
      write(32,'(2a/4a/4a)')
     &  '                                 qp-contributions           ',
     &  '                                 residual contributions     ',    
     &  '             norm <[P,Beta]>                   ',
     &  'width  <P^2>             ',
     &  '             norm <[P,Beta]>                   ',
     &  'width  <P^2> ',
     &  ' beta    mu= 0         1            2         ',
     &  ' 0           1            2 ',
     &  '       mu= 0         1            2         ',
     &  ' 0           1            2 '
c

      write(30,'(a/a,2i4/2a/2a)')
     &     title(1:78),
     &     ' A,Z=',npart(1)+npart(2),npart(1),
     &     '                      <P^2>                ',
     &     '                      <Q^2>',
     &     ' beta     mu= 0         1              2   ',
     &     '          mu= 0         1             2 '
      
            
      write(21,'(a)') ' '
      write(22,'(a)') ' '
c
      open(23,status='unknown',
     &        form='formatted',file='colpar_mp.'//outnam)
      open(24,status='unknown',
     &        form='formatted',file='colpar_GCM.'//outnam)
      open(25,status='unknown',
     &        form='formatted',file='colpar_in.'//outnam)
      open(26,status='unknown',
     &        form='formatted',file='colpar_sc.'//outnam)
      open(27,status='unknown',
     &        form='formatted',file='colpar_CHF.'//outnam)
      open(28,status='unknown',
     &        form='formatted',file='colpar_at.'//outnam)

c
      write(23,'(a)')
     &     ' parameters of quad.-dynamics (scaling mass)'
      write(24,'(a)')
     &     ' parameters of quad.-dynamics (static GCM)'
      write(25,'(a)')
     &     ' parameters of quad.-dynamics (inglis-cranking)'
      write(26,'(a)')
     &     ' parameters of quad.-dynamics (SC-cranking)'
      write(27,'(a)')
     &     ' parameters of quad.-dynamics (CHF + SC-Cranking)'
      write(28,'(a)')
     &     ' parameters of quad.-dynamics (full ATDHF)'
      
c
      do i=23,28
         write(i,'(a/a,2i4/7a)')
     &        title(1:78),
     &        ' A,Z=',npart(1)+npart(2),npart(1),
     &        ' beta_20     B_20     ',
     &        '  C_20    <P_20^2>    <Q_20^2>',
     &        '  <Omg_20^2>        ',
     &        'B_21        C_21    <P_21^2>    <Q_21^2>',
     &        '  <Omg_21^2>        ',
     &        'B_22        C_22    <P_22^2>    <Q_22^2>',
     &        '  <Omg_22^2>'
      enddo

c
c     ------------------------------------------------------------------
c
      return
      end




c=======================================================================
c
c     CPshift()
c
c=======================================================================
      subroutine CPshift(CPval,newval)
c
c     ------------------------------------------------------------------
c      
      double precision CPval(-1:1)
      double precision newval
c
c     ------------------------------------------------------------------
c     shift values in path storage CPpsi (-1,0,+1) -> (0,+1,+1)
c     and copy current value to last path storage index 
c     ------------------------------------------------------------------
      CPval(-1) = CPval( 0)
      CPval( 0) = CPval(+1)
      CPval(+1) = newval
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     static_afterburn()
c
c=======================================================================
      subroutine static_afterburn(ipath,varerr,maxiter)
c
c     ------------------------------------------------------------------
c 
      include 'paramr.inc'
      include 'linresp.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'levels.inc'
c
c     ------------------------------------------------------------------
c 
      integer          ipath
      double precision varerr
      integer          maxiter
c
      integer          ir,iz,im
      integer          ist,jst
      double precision energy(MPSID),var(MPSID),norm(MPSID)
      double precision maxvar,lastmaxvar,dmaxvar
      double precision HPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision HHPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision dzHPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision drHPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision dpHPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision ddHPsi(0:MZD,0:MRD,MSD,MPSID)
      double precision buffer,normStep,normError
      integer          touched,nlo,nhi,nadd,kpsic,iq
      integer          subiter
c
      double precision MFD_rho(0:MZD,0:MRD,2)
      double precision MFD_tau(0:MZD,0:MRD,2)
      double precision MFD_Jz(0:MZD,0:MRD,2)
      double precision MFD_Jr(0:MZD,0:MRD,2)
      double precision MFD_Grhor(0:MZD,0:MRD,2)
      double precision MFD_Grhoz(0:MZD,0:MRD,2)
      double precision MFD_Drho(0:MZD,0:MRD,2)
      double precision MFD_divJ(0:MZD,0:MRD,2)
c
      double precision MFP_U(0:MZD,0:MRD,2)
      double precision MFP_UNL(0:MZD,0:MRD,2)
      double precision MFP_B(0:MZD,0:MRD,2)
      double precision MFP_drB(0:MZD,0:MRD,2)
      double precision MFP_dzB(0:MZD,0:MRD,2)
      double precision MFP_Wz(0:MZD,0:MRD,2)
      double precision MFP_Wr(0:MZD,0:MRD,2)
c
      integer          istart(2),iend(2)
      double precision occ1,occ2
c
c     ------------------------------------------------------------------
c     
      write(21,'(a,g12.4,i8)')
     &  'STATIC_AFTERBURN(): improving sp-base: max. var, max. iter:'
     &  ,varerr,maxiter
c
c     calculate derivatives of single particle wavefunctions
c     ------------------------------------------------------------------
      call upd_derive_spinor(CPpsi(0,0,1,1,ipath),
     &                       CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &                       CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),
     &                       +1,0,nst(1)+nst(2))
c
c     accumulate the densities and mean-fields
c     ------------------------------------------------------------------
c      call upd_static_MF(CPpsi(0,0,1,1,ipath),
c     &                   CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
c     &                   CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),
c     &                   CPvsqf(1,ipath),
c     &                   MFP_U,MFP_B,MFP_drB,MFP_dzB,
c     &                   MFP_Wz,MFP_Wr)

c
      do ir=0,nr
         do iz=0,nz
            CPeffmas(iz,ir,1,ipath) = CPeffmas(iz,ir,1,ipath)-h2meff(1)
            CPeffmas(iz,ir,2,ipath) = CPeffmas(iz,ir,2,ipath)-h2meff(2) 
         enddo
      enddo
c      
      call nablaz_scalar(CPeffmas(0,0,1,ipath),
     &                   MFP_dzB(0,0,1),1,1,one)
      call nablar_scalar(CPeffmas(0,0,1,ipath),
     &                   MFP_drB(0,0,1),0,1,one) 
      call nablaz_scalar(CPeffmas(0,0,2,ipath),
     &                   MFP_dzB(0,0,2),1,1,one)
      call nablar_scalar(CPeffmas(0,0,2,ipath),
     &                   MFP_drB(0,0,2),0,1,one) 
c
      do ir=0,nr
         do iz=0,nz
            CPeffmas(iz,ir,1,ipath) = CPeffmas(iz,ir,1,ipath)+h2meff(1)
            CPeffmas(iz,ir,2,ipath) = CPeffmas(iz,ir,2,ipath)+h2meff(2)
         enddo
      enddo      
c
      do ir=0,nr
         do iz=0,nz
            MFP_U(iz,ir,1)          = CPlocpot(iz,ir,1,ipath)
            MFP_U(iz,ir,2)          = CPlocpot(iz,ir,2,ipath)
            MFP_B(iz,ir,1)          = CPeffmas(iz,ir,1,ipath)
            MFP_B(iz,ir,2)          = CPeffmas(iz,ir,2,ipath)
            MFP_Wr(iz,ir,1)         = CPsocpotr(iz,ir,1,ipath)
            MFP_Wr(iz,ir,2)         = CPsocpotr(iz,ir,2,ipath)
            MFP_Wz(iz,ir,1)         = CPsocpotz(iz,ir,1,ipath)
            MFP_Wz(iz,ir,2)         = CPsocpotz(iz,ir,2,ipath)
         enddo
      enddo
c
c     calculate error in Schroedinger-Eq.
c     ------------------------------------------------------------------
      call act_mfh_even(CPpsi(0,0,1,1,ipath),CPdrpsi(0,0,1,1),
     &                  CPdzpsi(0,0,1,1),CPdppsi(0,0,1,1),
     &                  CPddpsi(0,0,1,1),HPsi(0,0,1,1),
     &                  MFP_U,                                        1,
     &                  MFP_B,MFP_drB,MFP_dzB,    MFP_dzB,            1,
     &                  MFP_Wr,MFP_Wz,            MFP_Wz,             1,
     &                  CPconstr(0,0,ipath),CPlambda(ipath),          1,
     &                  CPspenrg(1,ipath),                            0,
     &                  0,1.0D0,
     &                  1,nst(1)+nst(2))
c
      call upd_derive_spinor(Hpsi(0,0,1,1),
     &                       drHpsi(0,0,1,1),dzHpsi(0,0,1,1),
     &                       dpHpsi(0,0,1,1),ddHpsi(0,0,1,1),
     &                       +1,0,nst(1)+nst(2))
c
      call act_mfh_even(Hpsi(0,0,1,1),drHpsi(0,0,1,1),
     &                  dzHpsi(0,0,1,1),dpHpsi(0,0,1,1),
     &                  ddHpsi(0,0,1,1),HHPsi(0,0,1,1),
     &                  MFP_U,                                        1,
     &                  MFP_B,MFP_drB,MFP_dzB,    MFP_dzB,            1,
     &                  MFP_Wr,MFP_Wz,            MFP_Wz,             1,
     &                  CPconstr(0,0,ipath),CPlambda(ipath),          1,
     &                  CPspenrg(1,ipath),                            0,
     &                  0,1.0D0,
     &                  1,nst(1)+nst(2))       

c
c     calculate variance and energy expectation
c     ------------------------------------------------------------------
      maxvar=0.0D0
      do ist=1,nst(1)+nst(2)
         energy(ist) = wfovlp(CPPsi(0,0,1,ist,ipath), HPsi(0,0,1,ist))
         var(ist)    = wfovlp(CPPsi(0,0,1,ist,ipath),HHPsi(0,0,1,ist))
c         var(ist)    = wfovlp(HPsi(0,0,1,ist),HPsi(0,0,1,ist))
         var(ist)    = var(ist)-energy(ist)**2
         maxvar=MAX(maxvar,abs(var(ist)))
c
c              write(88,'(a,3(1pg16.4),3i3)') 'PKtest1:',
c     &              energy(ist),var(ist), 
c     &              CPspenrg(ist,ipath),
c     &              mpsi(ist),nparz(1,ist),icharg(ist)
c
      enddo
c
      write(21,'(a,g12.3)') 
     & 'STATIC_AFTERBURN(): started   ---  max. var       : ',
     &             maxvar
c
c     the static iteration (with fixed mean-fields, occ. numbers, etc.)
c     ------------------------------------------------------------------
      lastmaxvar=0.0D0
      dmaxvar   =1E5
      do iter=1,maxiter
c
         if (maxvar.lt.varerr.or.iter.eq.maxiter) then
            do ist=1,nst(1)+nst(2)
               energy(ist) = wfovlp(CPPsi(0,0,1,ist,ipath),
     &                              HPsi(0,0,1,ist))
               var(ist)    = wfovlp(CPPsi(0,0,1,ist,ipath),
     &                              HHPsi(0,0,1,ist))
               var(ist)    = var(ist)-energy(ist)**2
c
c              write(88,'(a,3(1pg16.4),3i3)') 'PKtest:',
c     &              energy(ist),var(ist), 
c     &              CPspenrg(ist,ipath),
c     &              mpsi(ist),nparz(1,ist),icharg(ist)
            enddo
c
            write(21,'(a,g12.3,i8)') 
     &      'STATIC_AFTERBURN(): finished  ---  max. var, niter: ',
     &       maxvar,iter
            goto 4000
         endif
c
         dmaxvar=dabs(maxvar-lastmaxvar)
         lastmaxvar=maxvar
         if (mod(iter,100).eq.0) then
            write(21,'(a,i6,2g12.2)') 'afterburn: ',iter,maxvar,dmaxvar
         endif
c
c        loop all m-blocks
c        ---------------------------------------------------------------
         do im=1,nom
            touched = 0
            nhi     = mhi(im)
            nlo     = mlo(im)
            nadd    = nhi-nlo+1
            kpsic   = mim(im)
            iq      = miq(im)
c
c           now the gradient step 
c           ------------------------------------------------------------
            do ist=nlo,nhi
               if (abs(var(ist)).gt.varerr) then
                  touched = 1
               endif
            enddo
c
            if (touched.eq.1) then
c
               do ist=nlo,nhi
c
c                 calculate error in sp-states
c                 ---------------------------------------------
                  do ir=0,nr
                     do iz=0,nz
                        HHpsi(iz,ir,1,ist)=Hpsi(iz,ir,1,ist)
     &                   -energy(ist)*CPpsi(iz,ir,1,ist,ipath)
                        HHpsi(iz,ir,2,ist)=Hpsi(iz,ir,2,ist)
     &                   -energy(ist)*CPpsi(iz,ir,2,ist,ipath)
                     enddo
                  enddo
c
c                 apply damping
c                 ---------------------------------------------
                  call damnuc(HHPsi(0,0,1,ist),HPsi(0,0,1,ist),ist)
c
               enddo
c
               do ist=nlo,nhi
c
c                 perform the gradient step
c                 ---------------------------------------------
                  do ir=0,nr
                     do iz=0,nz
                        CPpsi(iz,ir,1,ist,ipath)=
     &                       CPpsi(iz,ir,1,ist,ipath)
     &                       -HPsi(iz,ir,1,ist)
                        CPpsi(iz,ir,2,ist,ipath)=
     &                       CPpsi(iz,ir,2,ist,ipath)
     &                       -HPsi(iz,ir,2,ist)
                     enddo
                  enddo
c                  
               enddo
c
c              orthonormalize m-block (only if touched=1)
c              --------------------------------------------
               call schmit(CPpsi(0,0,1,nlo,ipath),nlo-1,nadd)
c
c              update derivatives
c              --------------------------------------------
               call upd_derive_spinor_range(CPpsi(0,0,1,1,ipath),
     &                       CPdrpsi(0,0,1,1),CPdzpsi(0,0,1,1),
     &                       CPdppsi(0,0,1,1),CPddpsi(0,0,1,1),
     &                       +1,0,nlo,nhi)
c
c              calculate H|psi>
c              ---------------------------------------------
               call act_mfh_even(CPpsi(0,0,1,1,ipath),CPdrpsi(0,0,1,1),
     &                  CPdzpsi(0,0,1,1),CPdppsi(0,0,1,1),
     &                  CPddpsi(0,0,1,1),HPsi(0,0,1,1),
     &                  MFP_U,                                        1,
     &                  MFP_B,MFP_drB,MFP_dzB,MFP_dzB,                1,
     &                  MFP_Wr,MFP_Wz,            MFP_Wz,             1,
     &                  CPconstr(0,0,ipath),CPlambda(ipath),          1,
     &                  CPspenrg(1,ipath),                            0,
     &                  0,1.0D0,
     &                  nlo,nhi)
c
c              update derivatives
c              --------------------------------------------
               call upd_derive_spinor_range(Hpsi(0,0,1,1),
     &                       drHpsi(0,0,1,1),dzHpsi(0,0,1,1),
     &                       dpHpsi(0,0,1,1),ddHpsi(0,0,1,1),
     &                       +1,0,nlo,nhi)
c
c              calculate HH|psi>
c              ---------------------------------------------
               call act_mfh_even(Hpsi(0,0,1,1),drHpsi(0,0,1,1),
     &                  dzHpsi(0,0,1,1),dpHpsi(0,0,1,1),
     &                  ddHpsi(0,0,1,1),HHPsi(0,0,1,1),
     &                  MFP_U,                                        1,
     &                  MFP_B,MFP_drB,MFP_dzB,MFP_dzB,                1,
     &                  MFP_Wr,MFP_Wz,            MFP_Wz,             1,
     &                  CPconstr(0,0,ipath),CPlambda(ipath),          1,
     &                  CPspenrg(1,ipath),                            0,
     &                  0,1.0D0,
     &                  nlo,nhi)
c
            endif
c
         enddo
c
         maxvar   =0.0D0
         do ist=1,nst(1)+nst(2)
            energy(ist) = wfovlp(CPPsi(0,0,1,ist,ipath),
     &                           HPsi(0,0,1,ist))
            var(ist)    = wfovlp(CPPsi(0,0,1,ist,ipath),
     &                           HHPsi(0,0,1,ist))
            var(ist)    = var(ist)-energy(ist)**2
            maxvar=MAX(maxvar,abs(var(ist)))
c
      enddo
c     
      enddo
c
 4000 continue
c
      do ist=1,nst(1)+nst(2)
         CPspenrg(ist,ipath) = energy(ist)
      enddo
c
c      do ist=1,nst(1)+nst(2)
c         write(88,'(a,3(1pg16.4),3i3)') 'PKtest2:',
c     &              energy(ist),var(ist), 
c     &              CPspenrg(ist,ipath),
c     &              mpsi(ist),nparz(1,ist),icharg(ist)
c      enddo
c      write(88,*)
c

c
c     ------------------------------------------------------------------
c
      return
      end

c
      subroutine return_wavefunctions(ipath)
c
      include 'paramr.inc'
      include 'linresp.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'levels.inc'
c     
      integer ipath
c
      integer ist,ir,iz
c
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               psi(iz,ir,1,ist) = CPpsi(iz,ir,1,ist,ipath)
               psi(iz,ir,2,ist) = CPpsi(iz,ir,2,ist,ipath)
            enddo
         enddo
      enddo
c
      return
      end
