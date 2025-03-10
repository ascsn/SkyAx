#INCLUDE "define.h"
c======================================================================
c
c     comm_HQ(H_2q,Q_2q,Q_res,drQ_res,dzQ_res,dpQ_res,ddQ_res,
c                  HQ_2q,HQ_res,imu)
c
c     calculates the liouvillian {\cal L}Q = [H,Q] + H^{resp,Q} 
c     
c        imu    : angular-momentum transfer of Q-Operator (\mu)
c
c=======================================================================
      subroutine comm_HQ(H2q,
     &                   Q2q,ResQ,drResQ,dzResQ,dpResQ,ddResQ,
     &                   HQ2q,ResHQ,
     &                   imu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc' ! nst(1/2)
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTD)
      double precision Q2q(1:MTD)
      double precision HQ2q(1:MTD)
      double precision HResp2q(1:MTD)

#IF(FULLSPACE)
      double precision ResQ(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ResHQ(0:MZD,0:MRD,MSD,MPSID,2) 
      double precision ResHResp(0:MZD,0:MRD,MSD,MPSID,2)

      double precision drResQ(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dzResQ(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dpResQ(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ddResQ(0:MZD,0:MRD,MSD,MPSID,2)
#ELSE
      double precision ResQ
      double precision ResHQ
      double precision ResHResp
      double precision drResQ
      double precision dzResQ
      double precision dpResQ
      double precision ddResQ
#ENDIF

      integer          ist,itr,imu,ir,iz
c
c     ------------------------------------------------------------------
c
c     2qp-matrixelements of  i[H,Q] + H_{Q,resp} = 2 B (-iP)
c     ------------------------------------------------------------------
      if (isMFresp.eq.1) then
c     
c        calculate the Q-response of the mean-field
c        ---------------------------------------------------------------
         call upd_respdens_odd(Q2q,ResQ,drResQ,dzResQ,dpResQ,imu)
         call upd_respflds_odd(imu)
c         call wrt_response()
         call upd_resp_H_odd(HResp2q,ResHResp,imu) 
c  
         do itr=1,TRlen
            HQ2q(itr) =  H2q(itr)*Q2q(itr) - HResp2q(itr)
         enddo
c
      else
c
         do itr=1,TRlen
            HQ2q(itr) =  H2q(itr)*Q2q(itr)
         enddo
c
      endif
c
#IF(FULLSPACE)
c     if no residual action --> all done
c     ------------------------------------------------------------------
      if (isresact.eq.0) then
         return
      endif
c
c     action of i[H,Q] in residual space
c     ( CPsocz and CPsocr only act as unaccessed buffers )
c     ------------------------------------------------------------------
      call act_mfh_even(
     &            ResQ(0,0,1,1,1),  drResQ(0,0,1,1,1),
     &            dzResQ(0,0,1,1,1),dpResQ(0,0,1,1,1),
     &            ddResQ(0,0,1,1,1), ResHQ(0,0,1,1,1),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas(0,0,1),CPdzeffmas(0,0,1),     
     &            CPsocz(0,0,1,0),                         1, 
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocr(0,0,1,0),                         1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                          -1,
     &            0,1.0D0,
     &            1,nst(1)+nst(2))
      if (imu.ne.0) then
         call act_mfh_even(
     &            ResQ(0,0,1,1,2),  drResQ(0,0,1,1,2),
     &            dzResQ(0,0,1,1,2),dpResQ(0,0,1,1,2),
     &            ddResQ(0,0,1,1,2), ResHQ(0,0,1,1,2),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas(0,0,1),CPdzeffmas(0,0,1),     
     &            CPsocz(0,0,1,0),                         1, 
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocr(0,0,1,0),                         1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                          -1,
     &            0,1.0D0,
     &            1,nst(1)+nst(2))
      endif
c
c     add contributions from H_{Q,resp}
c     ------------------------------------------------------------------
      if (isMFresp.eq.1) then
         do ist=1,nst(1)+nst(2)
            do ir=0,nr
               do iz=0,nz
                  ResHQ(iz,ir,1,ist,1) = ResHQ(iz,ir,1,ist,1) 
     &                                  -ResHResp(iz,ir,1,ist,1)
                  ResHQ(iz,ir,2,ist,1) = ResHQ(iz,ir,2,ist,1) 
     &                                  -ResHResp(iz,ir,2,ist,1)
                  ResHQ(iz,ir,1,ist,2) = ResHQ(iz,ir,1,ist,2) 
     &                                  -ResHResp(iz,ir,1,ist,2)
                  ResHQ(iz,ir,2,ist,2) = ResHQ(iz,ir,2,ist,2) 
     &                                  -ResHResp(iz,ir,2,ist,2)
               enddo
            enddo
         enddo
      endif
c
c     project out contributions which lead back to qp-space 
c     ------------------------------------------------------------------
      call project_Omega(ResHQ,(-1)**imu,imu)
#ENDIF
c
c     ------------------------------------------------------------------
c  
      return
      end




c======================================================================
c
c     comm_HP(H_2q,P_2q,P_res,drP_res,dzP_res,dpP_res,ddP_res,
c                  HP_2q,HP_res,isresp,imu)
c
c     calculate the liouvillian {\cal L}P = [H,P] + H^{resp,P} 
c     
c        isresp : enable/disable response of mean-field
c        imu    : angular-momentum transfer of P-Operator (\mu)
c
c=======================================================================
      subroutine comm_HP(H2q,
     &                   P2q,ResP,drResP,dzResP,dpResP,ddResP,
     &                   HP2q,ResHP,
     &                   imu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc' ! nst(1/2)
c
c     ------------------------------------------------------------------
c
      double precision H2q(1:MTD)
      double precision P2q(1:MTD)
      double precision HP2q(1:MTD)
      double precision HResp2q(1:MTD)

#IF(FULLSPACE)
      double precision ResP(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ResHP(0:MZD,0:MRD,MSD,MPSID,2) 
      double precision ResHResp(0:MZD,0:MRD,MSD,MPSID,2)
      double precision drResP(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dzResP(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dpResP(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ddResP(0:MZD,0:MRD,MSD,MPSID,2)
#ELSE
      double precision ResP    
      double precision ResHP    
      double precision ResHResp
      double precision drResP
      double precision dzResP
      double precision dpResP
      double precision ddResP
#ENDIF

      integer          ist,itr,imu,ir,iz
c
c     ------------------------------------------------------------------
c
c     2qp-matrixelements of  -i[H,P] + H_{P,resp} = 2 C Q
c     ------------------------------------------------------------------
      if (isMFresp.eq.1) then
c     
c        calculate the P-response of the mean-field
c        ---------------------------------------------------------------
c        call upd_respdens_even_CHF()
c        call upd_respflds_even_CHF()
c        call wrt_response()        
        call upd_respdens_even(P2q,ResP,drResP,dzResP,dpResP,imu)
        call upd_respflds_even(imu)
c        call wrt_response()
         call upd_resp_H_even(HResp2q,ResHResp,imu)     
c
         do itr=1,TRlen
            HP2q(itr) = H2q(itr)*P2q(itr) - HResp2q(itr)
         enddo
c     
      else
c
         do itr=1,TRlen
            HP2q(itr) = H2q(itr)*P2q(itr)
         enddo
c
      endif
c
#IF(FULLSPACE)
c     if no residual action --> all done
c     ------------------------------------------------------------------
      if (isresact.eq.0) then
         return
      endif
c
c     action of -i[H,P] in residual space
c     ------------------------------------------------------------------
      call act_mfh_even(
     &            ResP(0,0,1,1,1),  drResP(0,0,1,1,1),
     &            dzResP(0,0,1,1,1),dpResP(0,0,1,1,1),
     &            ddResP(0,0,1,1,1), ResHP(0,0,1,1,1),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0,
     &            1,nst(1)+nst(2))
      call act_mfh_even(
     &            ResP(0,0,1,1,2),  drResP(0,0,1,1,2),
     &            dzResP(0,0,1,1,2),dpResP(0,0,1,1,2),
     &            ddResP(0,0,1,1,2), ResHP(0,0,1,1,2),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0,
     &            1,nst(1)+nst(2))
c
c     add contributions from H_{P,resp}
c     ------------------------------------------------------------------
      if (isMFresp.eq.1) then
         do ist=1,nst(1)+nst(2)
            do ir=0,nr
               do iz=0,nz
                  ResHP(iz,ir,1,ist,1) = ResHP(iz,ir,1,ist,1) 
     &                                 - ResHResp(iz,ir,1,ist,1)
                  ResHP(iz,ir,2,ist,1) = ResHP(iz,ir,2,ist,1) 
     &                                 - ResHResp(iz,ir,2,ist,1)
                  ResHP(iz,ir,1,ist,2) = ResHP(iz,ir,1,ist,2) 
     &                                 - ResHResp(iz,ir,1,ist,2)
                  ResHP(iz,ir,2,ist,2) = ResHP(iz,ir,2,ist,2) 
     &                                 - ResHResp(iz,ir,2,ist,2)
               enddo
            enddo
         enddo
      endif
c
c     project out contributions which lead back to qp-space 
c     ------------------------------------------------------------------
      call project_Omega(ResHP,(-1)**imu,imu)
#ENDIF
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     comm_NQ(H_2q,Q_2q,Q_res,drQ_res,dzQ_res,dpQ_res,ddQ_res,
c                  HQ_2q,HQ_res,imu)
c
c     calculates the liouvillian {\cal L}Q = [H,Q] + H^{resp,Q} 
c     
c        imu    : angular-momentum transfer of Q-Operator (\mu)
c
c=======================================================================
      subroutine comm_NQ(Q_2q,ResQ,NQ_2q,ResNQ)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc' ! nst(1/2)
c
c     ------------------------------------------------------------------
c
      double precision Q_2q(1:MTD)
      double precision NQ_2q(1:MTD)

#IF(FULLSPACE)
      double precision ResQ(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ResNQ(0:MZD,0:MRD,MSD,MPSID,2)
#ELSE
      double precision ResQ
      double precision ResNQ
#ENDIF

      integer          ist,jst,itr,ir,iz
c
c     ------------------------------------------------------------------
c
c     2qp-matrixelements of  [N,Q] 
c     ------------------------------------------------------------------
      do itr=1,TRlen
c
         ist=TRidx(itr,1)   ! initial state index
         jst=TRidx(itr,2)   ! final state index

         NQ_2q(itr) =  2.0D0*
     &         (1.0D0 - CPvsqf(ist,0) - CPvsqf(jst,0))*Q_2q(itr)
c
      enddo
c
#IF(FULLSPACE)
c     residual action == zero
c     ------------------------------------------------------------------
      do ist=1,nsttot
         do ir=0,nr
            do iz=0,nz
               ResNQ(iz,ir,1,ist,1) = 0.0D0
               ResNQ(iz,ir,2,ist,1) = 0.0D0
               ResNQ(iz,ir,1,ist,2) = 0.0D0
               ResNQ(iz,ir,2,ist,2) = 0.0D0
            enddo
         enddo
      enddo
#ENDIF
c
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     comm_QP(Q_2q,Q_res,P_2q,P_res,imu)
c
c     calculate the commutator expectation <[Q,P]>
c     
c        imu    : angular-momentum transfer of Q and P-Operator (\mu)
c
c=======================================================================
      function comm_QP(Q2q,ResQ,P2q,ResP,imu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'linresp.inc'  ! CPvsqf,ntr,MTR
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision Q2q(1:MTD)
      double precision P2q(1:MTD)

#IF(FULLSPACE)
      double precision ResQ(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ResP(0:MZD,0:MRD,MSD,MPSID,2)
#ELSE
      double precision ResQ
      double precision ResP
#ENDIF

      integer          ist,itr
      double precision comm_QP,comm_QP_res,comm_QP_qp
c
c     ------------------------------------------------------------------
c     
      comm_QP_res=0.0D0
      comm_QP_qp =0.0D0
c
c     contributions from 2qp-matrixelements
c     ------------------------------------------------------------------
      do itr=1,TRlen
         comm_QP_qp = comm_QP_qp + TRidx(itr,6)*Q2q(itr)*P2q(itr)
      enddo
c
#IF(FULLSPACE)
c     contributions from residual action
c     ( factor 4 results from occupation (2) and commutator (2) )
c     ------------------------------------------------------------------
      if (isresact.ne.0) then
         if (imu.eq.0) then
            do ist=1,nst(1)+nst(2)
               comm_QP_res = comm_QP_res + 
     &              4.0D0*CPvsqf(ist,0)*                
     &              wfovlp(ResQ(0,0,1,ist,1) ,ResP(0,0,1,ist,1))
            enddo
         else
            do ist=1,nst(1)+nst(2)
               comm_QP_res = comm_QP_res + 
     &              4.0D0*CPvsqf(ist,0)*(                
     &              wfovlp(ResQ(0,0,1,ist,1) ,ResP(0,0,1,ist,1))+
     &              wfovlp(ResQ(0,0,1,ist,2) ,ResP(0,0,1,ist,2)) )
            enddo
         endif
      endif
c
      comm_QP = comm_QP_res + comm_QP_qp
#ELSE
      comm_QP = comm_QP_qp     
#ENDIF
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     var_oper(X_2q,X_res,imu)      
c
c     calculate the variance of an operator <X^2>  X=Q,P
c     
c        imu    : angular-momentum transfer of X-Operator (\mu)
c
c=======================================================================
      function var_oper(X2q,ResX,imu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'linresp.inc'  ! CPvsqf,TRlen,MTR,TRidx
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision X2q(1:MTD)
      double precision ResX(0:MZD,0:MRD,MSD,MPSID,2)

      integer          ist,itr
      double precision var_oper,var_X_res,var_X_qp
c
c     ------------------------------------------------------------------
c     
      var_X_res=0.0D0
      var_X_qp =0.0D0
c
c     contributions from 2qp-matrixelements
c     ------------------------------------------------------------------
      do itr=1,TRlen
         var_X_qp = var_X_qp + TRidx(itr,6)*X2q(itr)*X2q(itr)
      enddo
c
#IF(FULLSPACE)
c     contributions from residual action
c     ------------------------------------------------------------------
      if (isresact.ne.0) then
         if (imu.eq.0) then
            do ist=1,nst(1)+nst(2)
               var_X_res = var_X_res + 
     &              CPvsqf(ist,0)*                
     &              wfovlp(ResX(0,0,1,ist,1) ,ResX(0,0,1,ist,1))
            enddo
         else
            do ist=1,nst(1)+nst(2)
               var_X_res = var_X_res + 
     &              CPvsqf(ist,0)*(                
     &              wfovlp(ResX(0,0,1,ist,1) ,ResX(0,0,1,ist,1))+
     &              wfovlp(ResX(0,0,1,ist,2) ,ResX(0,0,1,ist,2)) )
            enddo
         endif
      endif
c
      var_oper = 2.0D0*var_X_res + 0.5D0*var_X_qp
#ELSE
      var_oper = 0.5D0*var_X_qp
#ENDIF
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     zero_oper(X_2q,ResX,ipart)      
c
c     set operator matrixelements (ipart=1) residual part (ipart=2) or
c     both (ipart=0) to zero
c     
c=======================================================================
      subroutine zero_oper(X2q,ResX,ipart)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'dynpar.inc'   ! nz,nr, etc.
      include 'linresp.inc'  ! CPvsqf,TRlen,MTR,TRidx
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision X2q(1:MTD)

#IF(FULLSPACE)
      double precision ResX(0:MZD,0:MRD,MSD,MPSID,2)
#ELSE
      double precision ResX
#ENDIF

      integer          ipart

      integer          ist,ir,iz,itr
c
c     ------------------------------------------------------------------
c
c     contributions from 2qp-matrixelements
c     ------------------------------------------------------------------
      if (ipart.ne.2) then
         do itr=1,TRlen
            X2q(itr)=0.0D0
         end do
      endif
c
#IF(FULLSPACE)
      if (ipart.ne.1) then
         do ist=1,nst(1)+nst(2)
            do ir=0,nr
               do iz=0,nz 
                  ResX(iz,ir,1,ist,1)=0.0D0
                  ResX(iz,ir,2,ist,1)=0.0D0
                  ResX(iz,ir,1,ist,2)=0.0D0
                  ResX(iz,ir,2,ist,2)=0.0D0
               enddo
            enddo
         enddo
      endif
#ENDIF
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     project_qp(X_2q,ResX,Xtpar,imu,iaccu)
c
c     project out components in sp-space from an operator-wavefunction
c     (=Pi_Omega)
c     
c        X_2q   : 2qp-matrixelements of ResX
c        ResX  : projection of ResX on residual space (output)
c                 or unprojected operator-wavefunctions (input)
c        imu    : angular-momentum transfer of X-Operator (\mu)
c        iaccu  : switch if matrixlements should be accumulated in X_2q
c                 or     if X_2q should be zeroed (iaccu=0)
c
c=======================================================================
      subroutine project_qp(X2q,ResX,tpar,imu,iaccu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'dynpar.inc'   ! nr,nz
      include 'linresp.inc'  ! CPvsqf,ntr,MTR
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision X2q(1:MTD)
      double precision ResX(0:MZD,0:MRD,MSD,MPSID,2)
      integer          imu,tpar,iaccu

      integer          ist,jst,jtp,musgn,tpidx,iz,ir,itr
      double precision ovlp,ovlp2
      double precision fac
c
c     ------------------------------------------------------------------
c
      tpidx = (3-tpar)/2          ! +1 |--> 1,  -1 |--> 2
c
      if (iaccu.eq.0) then
         do itr=1,TRlen
            X2q(itr)=0.0D0
         enddo
      endif
c
c     ==================================================================
c     if residual action is included
c     ==================================================================
      if (isresact.eq.1) then
c
c     the case imu == 0      ==>  only ResX(0,0,1,1,1) is used !!!
c     ==================================================================
      if (imu.eq.0) then
c
c     loop all non-vanishing matrixelements 
c     ------------------------------------------------------------------
      do itr=1,TRlen
         ist       = TRidx(itr,1)
         jst       = TRidx(itr,2)
c
c        calculate <jst| X |ist> and
c        project out |jst> from |Xist> and |ist> from |Xjst>
c        ---------------------------------------------------------------
         if (ist.ne.jst) then
           ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1))  
           call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1),ovlp)
           call wfproj(CPpsi(0,0,1,ist,0),ResX(0,0,1,jst,1),tpar*ovlp)
         else
           ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1))  
           call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1),ovlp)
         endif
c         
c        transform matrixelement into qp-space
c        ---------------------------------------------------------------
         X2q(itr) = X2q(itr) + TReta(itr,tpidx)*ovlp
c
      enddo
c
      else
c
c     the case imu != 0  ==>  both components of ResX are used !!!
c     ==================================================================
c
c     loop all non-vanishing matrixelements 
c     ------------------------------------------------------------------
      do itr=1,TRlen
         ist       = TRidx(itr,1)
         jst       = TRidx(itr,2)
         jtp       = TRidx(itr,3)
         musgn     = TRidx(itr,5)
c
c        calculate <jst| X |ist> and
c        project out |jst> from |Xist> and |ist> from |Xjst>
c        ---------------------------------------------------------------
         if (jtp.eq.+1) then
         if (musgn.eq.+1) then
            ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1))  
            call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1),ovlp)
            call wfproj(CPpsi(0,0,1,ist,0),ResX(0,0,1,jst,2),tpar*ovlp)
         else
            ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2))  
            call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),ovlp)
            call wfproj(CPpsi(0,0,1,ist,0),ResX(0,0,1,jst,1),tpar*ovlp)
         endif
         else
            if (ist.ne.jst) then
               ovlp=wfovlptr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2)) 
               call wfprojtr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),ovlp)
               call wfprojtr(CPpsi(0,0,1,ist,0),ResX(0,0,1,jst,2),
     &                                                       -tpar*ovlp)
            else
               ovlp=wfovlptr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2)) 
               call wfprojtr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),ovlp)
            endif
         endif
c         
c        transform matrixelement into qp-space
c        ---------------------------------------------------------------
         X2q(itr) = X2q(itr) + TReta(itr,tpidx)*ovlp
c
      enddo
c     
      endif
c
c     ==================================================================
c     if residual action is NOT included
c     ==================================================================
      else
c
c     the case imu == 0      ==>  only ResX(0,0,1,1,1) is used !!!
c     ==================================================================
      if (imu.eq.0) then
c
c     loop all non-vanishing matrixelements 
c     ------------------------------------------------------------------
      do itr=1,TRlen
         ist       = TRidx(itr,1)
         jst       = TRidx(itr,2)
c
c        calculate <jst| X |ist> and
c        project out |jst> from |Xist> and |ist> from |Xjst>
c        ---------------------------------------------------------------
         ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1))  
c         
c        transform matrixelement into qp-space
c        ---------------------------------------------------------------
         X2q(itr) = X2q(itr) + TReta(itr,tpidx)*ovlp
c
      enddo
c
      else
c
c     the case imu != 0  ==>  both components of ResX are used !!!
c     ==================================================================
c
c     loop all non-vanishing matrixelements 
c     ------------------------------------------------------------------
      do itr=1,TRlen
         ist       = TRidx(itr,1)
         jst       = TRidx(itr,2)
         jtp       = TRidx(itr,3)
         musgn     = TRidx(itr,5)
c
c        calculate <jst| X |ist> and
c        project out |jst> from |Xist> and |ist> from |Xjst>
c        ---------------------------------------------------------------
         if (jtp.eq.+1) then
         if (musgn.eq.+1) then
            ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1))  
         else
            ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2))  
         endif
         else
            ovlp=wfovlptr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2)) 
         endif
c         
c        transform matrixelement into qp-space
c        ---------------------------------------------------------------
         X2q(itr) = X2q(itr) + TReta(itr,tpidx)*ovlp
c
      enddo
c     
      endif
      endif
c
#IF(FULLSPACE)
c     if residual action is neglected ==> zero residual action
c     ------------------------------------------------------------------
      if (isresact.eq.0) then
c
         do ist=1,nst(1)+nst(2)
            do ir=0,nr
               do iz=0,nz
                  ResX(iz,ir,1,ist,1) = 0.0D0
                  ResX(iz,ir,2,ist,1) = 0.0D0
                  ResX(iz,ir,1,ist,2) = 0.0D0
                  ResX(iz,ir,2,ist,2) = 0.0D0
               enddo
            enddo  
         enddo
c
      else
c
         call project_Omega(ResX,(-1)**imu,imu) 
c
      endif
#ENDIF
c
c
c     ------------------------------------------------------------------
c  
      return
      end


#IF(FULLSPACE)
c======================================================================
c
c     project_Omega(X_2q,X_res,Xtpar,imu,iaccu)
c
c     checks the sp-contributions of the residual action of an operator
c     
c        X_res  : projection of X_res on residual space (output)
c                 or unprojected operator-wavefunctions (input)
c        imu    : angular-momentum transfer of X-Operator (\mu)
c
c=======================================================================
      subroutine project_Omega(ResX,trzpar,trmcomp)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'dynpar.inc'   ! nr,nz
      include 'linresp.inc'  ! CPvsqf,ntr,MTR
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision ResX(0:MZD,0:MRD,MSD,MPSID,2)
      integer          trmcomp,trzpar

      integer          ist,jst,jtp,itp,jst2
      double precision ovlp
c
c     project out sp-components of ResX
c     ------------------------------------------------------------------
      do ist=1,nst(1)
         do jst=1,nst(1)
c
            if (mpsi(ist)+trmcomp.eq.mpsi(jst)) then
               if (npar(ist)*npar(jst)*trzpar.ne.-1) then
                  ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1))  
                  call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1),
     &                        ovlp)
               endif
            elseif (mpsi(ist)-trmcomp.eq.mpsi(jst)) then
               if (npar(ist)*npar(jst)*trzpar.ne.-1) then
                  ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2))  
                  call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),
     &                        ovlp)
               endif
            elseif (mpsi(ist)-trmcomp.eq.(-mpsi(jst)-1)) then
               if (npar(ist)*npar(jst)*trzpar.ne.+1) then
                  ovlp=wfovlptr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2))  
                  call wfprojtr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),
     &                          ovlp)
               endif
            endif
c
         enddo
      enddo
c
      do ist=nst(1)+1,nst(2)
         do jst=nst(1)+1,nst(2)
c
            if (mpsi(ist)+trmcomp.eq.mpsi(jst)) then
               if (npar(ist)*npar(jst)*trzpar.ne.-1) then
                  ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1)) 
                  call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,1),
     &                        ovlp) 
               endif
            elseif (mpsi(ist)-trmcomp.eq.mpsi(jst)) then
               if (npar(ist)*npar(jst)*trzpar.ne.-1) then
                  ovlp=wfovlp(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2))  
                  call wfproj(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),
     &                        ovlp)
               endif
            elseif (mpsi(ist)-trmcomp.eq.(-mpsi(jst)-1)) then
               if (npar(ist)*npar(jst)*trzpar.ne.+1) then
                  ovlp=wfovlptr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2)) 
                  call wfprojtr(CPpsi(0,0,1,jst,0),ResX(0,0,1,ist,2),
     &                          ovlp) 
               endif
            endif
c     
         enddo
      enddo
c
c     ------------------------------------------------------------------
c  
      return
      end
#ENDIF


c======================================================================
c
c     upd_respflds_even_CHF()
c
c       updates the time-even response fields from the collective 
c       path by explicit evaluation of the change of the MF-potentials
c
c=======================================================================
      subroutine upd_respflds_even_CHF()
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
      integer          ir,iz
      double precision invdb
c
c     ------------------------------------------------------------------
c
      invdb = 1.0D0/(CPbeta(+1)-CPbeta(-1)) 
c
      do ir=0,nr
         do iz=0,nz
c
            relocpot(iz,ir,1)=
     &           (CPlocpot(iz,ir,1,+1)-CPlocpot(iz,ir,1,-1))*invdb
            relocpot(iz,ir,2)=
     &           (CPlocpot(iz,ir,2,+1)-CPlocpot(iz,ir,2,-1))*invdb
            reeffmas(iz,ir,1)=
     &           (CPeffmas(iz,ir,1,+1)-CPeffmas(iz,ir,1,-1))*invdb
            reeffmas(iz,ir,2)=
     &           (CPeffmas(iz,ir,2,+1)-CPeffmas(iz,ir,2,-1))*invdb
            resocpotr(iz,ir,1)=
     &           (CPsocpotr(iz,ir,1,+1)-CPsocpotr(iz,ir,1,-1))*invdb
            resocpotr(iz,ir,2)=
     &           (CPsocpotr(iz,ir,2,+1)-CPsocpotr(iz,ir,2,-1))*invdb
            resocpotz(iz,ir,1)=
     &           (CPsocpotz(iz,ir,1,+1)-CPsocpotz(iz,ir,1,-1))*invdb
            resocpotz(iz,ir,2)=
     &           (CPsocpotz(iz,ir,2,+1)-CPsocpotz(iz,ir,2,-1))*invdb
c
         enddo
      enddo
c
      do ir=0,nr
         do iz=0,nz
c
            relocpot(iz,ir,1)=relocpot(iz,ir,1) -
     &           ( CPlambda(+1)*CPconstr(iz,ir,+1)
     &            -CPlambda(-1)*CPconstr(iz,ir,-1))*invdb
            relocpot(iz,ir,2)=relocpot(iz,ir,2) -
     &           ( CPlambda(+1)*CPconstr(iz,ir,+1)
     &            -CPlambda(-1)*CPconstr(iz,ir,-1))*invdb
c 
         enddo
      enddo
c
c
c     ------------------------------------------------------------------
c
c      call wrt_2D_spin(relocpot(0,0,1) ,'relocpot_C',1,'tst')
c      call wrt_2D_spin(resocpotr(0,0,1),'resocpotrC',0,'tst')
c      call wrt_2D_spin(resocpotz(0,0,1),'resocpotzC',0,'tst')
c      call wrt_2D_spin(resocpotp(0,0,1),'resocpotpC',0,'tst')
c      call wrt_2D_spin(reeffmas(0,0,1), 'reeffmas_C',0,'tst')
c
c     ------------------------------------------------------------------
c
      return
      end

c======================================================================
c
c     upd_respdens_even_CHF()
c
c       updates the time-even response densities from the collective 
c       path by explicit evaluation of the change of the MF-potentials
c
c=======================================================================
      subroutine upd_respdens_even_CHF()
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
      integer          ir,iz
      double precision invdb
c
c     ------------------------------------------------------------------
c
      invdb = 1.0D0/(CPbeta(+1)-CPbeta(-1)) 
c
      do ir=0,nr
         do iz=0,nz
c
            rerho(iz,ir,1)=
     &           (CPrho(iz,ir,1,+1)-CPrho(iz,ir,1,-1))*invdb
            rerho(iz,ir,2)=
     &           (CPrho(iz,ir,2,+1)-CPrho(iz,ir,2,-1))*invdb
            rekin(iz,ir,1)=
     &           (CPkin(iz,ir,1,+1)-CPkin(iz,ir,1,-1))*invdb
            rekin(iz,ir,2)=
     &           (CPkin(iz,ir,2,+1)-CPkin(iz,ir,2,-1))*invdb
            resocr(iz,ir,1)=
     &           (CPsocr(iz,ir,1,+1)-CPsocr(iz,ir,1,-1))*invdb
            resocr(iz,ir,2)=
     &           (CPsocr(iz,ir,2,+1)-CPsocr(iz,ir,2,-1))*invdb
            resocz(iz,ir,1)=
     &           (CPsocz(iz,ir,1,+1)-CPsocz(iz,ir,1,-1))*invdb
            resocz(iz,ir,2)=
     &           (CPsocz(iz,ir,2,+1)-CPsocz(iz,ir,2,-1))*invdb
c
         enddo
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end


c======================================================================
c
c     upd_respdens_even(P_2q,P_res,drP_res,dzP_res,dpP_res)
c
c       updates the time-even response densities resulting from 
c       a small change described by the time-odd generator P with 
c       2q-matrixelements 'P_2q' and residual action 'P_res'. 
c
c     imu: angular-momentum transfer of P-operator (mu)
c
c=======================================================================
#IF(FULLSPACE)
      subroutine upd_respdens_even(P_2q,P_res,
     &                             drP_res,dzP_res,dpP_res,
     &                             imu)
#ELSE
      subroutine upd_respdens_even(P_2q,dummy1,
     &                             dummy2,dummy3,dummy4,
     &                             imu)
#ENDIF
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
      double precision P_2q(1:MTD)
#IF(FULLSPACE)
      double precision P_res(0:MZD,0:MRD,MSD,MPSID)
      double precision drP_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dzP_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dpP_res(0:MZD,0:MRD,MSD,MPSID)
#ELSE
      double precision dummy1,dummy2,dummy3,dummy4
#ENDIF
      double precision wgtfac,wgtfac0,wgtfac2,multi,symfac
      double precision buff(0:9)
      integer          itr,ist,jst,ir,iz,jtp,iq,imu
c
c     ------------------------------------------------------------------
c
      symfac=1.0D0
      if (imu.eq.0) symfac=0.0D0
c
c     set time-even transition densities to zero
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            rerho(iz,ir,1)  = zero
            rerho(iz,ir,2)  = zero
            rekin(iz,ir,1)  = zero
            rekin(iz,ir,2)  = zero
            resocr(iz,ir,1) = zero
            resocr(iz,ir,2) = zero
            resocp(iz,ir,1) = zero
            resocp(iz,ir,2) = zero
            resocz(iz,ir,1) = zero
            resocz(iz,ir,2) = zero
            repair(iz,ir,1) = zero
            repair(iz,ir,2) = zero
         enddo
      enddo
c                  
c     accumulate time-even response densities from 2q-matrixelements
c     ------------------------------------------------------------------
      do itr=1,TRlen
c
         ist     = TRidx(itr,1)
         jst     = TRidx(itr,2)
         jtp     = TRidx(itr,3)
         iq      = TRidx(itr,4)
         multi   = TRidx(itr,6)

c        
         wgtfac  = multi*TReta(itr,1)*P_2q(itr) ! HF-densities
         wgtfac0 = multi*TReta(itr,1)*P_2q(itr)*symfac
         wgtfac2 = multi* TRxi(itr,1)*P_2q(itr) ! Pairing-densities
c                  
c        if state jst has positive "T-parity": T|jst>=|Tjst>
c        ---------------------------------------------------------------
         if (jtp.eq.+1) then
c
            do ir=0,nr
               do iz=0,nz
c                  
c                 propability density                                 
c                 ------------------------------------------------------
                  rerho(iz,ir,iq)  = wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0) 
     &                +CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + rerho(iz,ir,iq)

c                  
c                 kinetic density                                     
c                 ------------------------------------------------------
                  rekin(iz,ir,iq)  = wgtfac*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist) 
     &                +CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist) 
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist) 
     &                +CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist) 
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist) 
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)
     &                                      ) + rekin(iz,ir,iq)
c                  
c                 spin-orbit current                                  
c                 ------------------------------------------------------
                  resocr(iz,ir,iq) = wgtfac*0.5D0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,2,ist) 
     &                +CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,1,ist) 
     &                +CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                +CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                                      ) + resocr(iz,ir,iq)
                  
                  resocp(iz,ir,iq) = wgtfac0*0.5D0*(
     &                -CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,2,ist) 
     &                +CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                                      ) + resocp(iz,ir,iq)
c                  
                  resocz(iz,ir,iq) = wgtfac*0.5D0*(
     &                -CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,2,ist) 
     &                +CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,1,ist) 
     &                -CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                -CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                                      ) + resocz(iz,ir,iq)
c                  
c                 pairing density                                 
c                 ------------------------------------------------------
                  repair(iz,ir,iq)  = wgtfac2*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0) 
     &                +CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + repair(iz,ir,iq)
c
               enddo
            enddo
c                  
c           if state jst has negative "T-parity": T|jst>=-|Tjst>
c           ------------------------------------------------------------
            else 
c              
            do ir=0,nr
               do iz=0,nz
c                  
c                 propability density                                 
c                 ------------------------------------------------------
                  rerho(iz,ir,iq)  = wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0) 
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + rerho(iz,ir,iq)

c                  
c                 kinetic density                                     
c                 ------------------------------------------------------
                  rekin(iz,ir,iq)  = wgtfac*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist) 
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist) 
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist) 
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist) 
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist) 
     &                -CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)
     &                                      ) + rekin(iz,ir,iq)
c                  
c                 spin-orbit current                                  
c                 ------------------------------------------------------
                  resocr(iz,ir,iq) = wgtfac*0.5D0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,1,ist) 
     &                +CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                +CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                                      ) + resocr(iz,ir,iq)
                  
                  resocp(iz,ir,iq) = wgtfac0*0.5D0*(
     &                +CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,2,ist) 
     &                +CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                +CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                                      ) + resocp(iz,ir,iq)
c                  
                  resocz(iz,ir,iq) = wgtfac*0.5D0*(
     &                +CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,2,ist) 
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,1,ist) 
     &                +CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                                      ) + resocz(iz,ir,iq)
c                  
c                 pairing density                                 
c                 ------------------------------------------------------
                  repair(iz,ir,iq)  = wgtfac2*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0) 
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + repair(iz,ir,iq)
c
               enddo
            enddo 
c
         endif
c
      enddo
c          
#IF(FULLSPACE)        
c     accumulate time-even response densities from residual-components
c     ------------------------------------------------------------------
      if (isresact.eq.1) then
c
c      do ist=1,nsttot
cc     
c         wgtfac = 4.0D0*CPvsqf(ist,0)  ! comm.(x2), time-reversal (x2)
c         wgtfac0= 4.0D0*CPvsqf(ist,0)*symfac
cc
c         if (ist.le.nst(1)) then
c            iq = 1
c         else
c            iq = 2
c         endif
cc
c         do ir=0,nr
c            do iz=0,nz
ccc                  
cc           propability density                                 
cc           ------------------------------------------------------------
c               rerho(iz,ir,iq)  = wgtfac*(
c     &             +CPpsi(iz,ir,1,ist,0) * P_res(iz,ir,1,ist,1) 
c     &             +CPpsi(iz,ir,2,ist,0) * P_res(iz,ir,2,ist,1)
c     &             +CPpsi(iz,ir,1,ist,0) * P_res(iz,ir,1,ist,2) 
c     &             +CPpsi(iz,ir,2,ist,0) * P_res(iz,ir,2,ist,2)
c     &                                   ) + rerho(iz,ir,iq)
c
cc                  
cc                 kinetic density                                     
cc                 ------------------------------------------------------
c                  rekin(iz,ir,iq)  = wgtfac*(
c     &                +CPdrpsi(iz,ir,1,ist) * drP_res(iz,ir,1,ist,1)
c     &                +CPdrpsi(iz,ir,2,ist) * drP_res(iz,ir,2,ist,1)
c     &                +CPdzpsi(iz,ir,1,ist) * dzP_res(iz,ir,1,ist,1)
c     &                +CPdzpsi(iz,ir,2,ist) * dzP_res(iz,ir,2,ist,1)
c     &                +CPdppsi(iz,ir,1,ist) * dpP_res(iz,ir,1,ist,1)
c     &                +CPdppsi(iz,ir,2,ist) * dpP_res(iz,ir,2,ist,1)
c     &                +CPdrpsi(iz,ir,1,ist) * drP_res(iz,ir,1,ist,2)
c     &                +CPdrpsi(iz,ir,2,ist) * drP_res(iz,ir,2,ist,2)
c     &                +CPdzpsi(iz,ir,1,ist) * dzP_res(iz,ir,1,ist,2)
c     &                +CPdzpsi(iz,ir,2,ist) * dzP_res(iz,ir,2,ist,2)
c     &                +CPdppsi(iz,ir,1,ist) * dpP_res(iz,ir,1,ist,2)
c     &                +CPdppsi(iz,ir,2,ist) * dpP_res(iz,ir,2,ist,2)
c     &                                      ) + rekin(iz,ir,iq)
cc                  
cc                 spin-orbit current                                  
cc                 ------------------------------------------------------
c                  resocr(iz,ir,iq) = wgtfac*0.5D0*(
c     &                +CPpsi(iz,ir,1,ist,0) * dpP_res(iz,ir,1,ist,1) 
c     &                -CPpsi(iz,ir,2,ist,0) * dpP_res(iz,ir,2,ist,1) 
c     &                +CPpsi(iz,ir,1,ist,0) * dzP_res(iz,ir,2,ist,1)   
c     &                -CPpsi(iz,ir,2,ist,0) * dzP_res(iz,ir,1,ist,1)
c     &                +CPdppsi(iz,ir,1,ist) *   P_res(iz,ir,1,ist,1) 
c     &                -CPdppsi(iz,ir,2,ist) *   P_res(iz,ir,2,ist,1) 
c     &                +CPdzpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,1)   
c     &                -CPdzpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,1)
c     &                +CPpsi(iz,ir,1,ist,0) * dpP_res(iz,ir,1,ist,2) 
c     &                -CPpsi(iz,ir,2,ist,0) * dpP_res(iz,ir,2,ist,2) 
c     &                +CPpsi(iz,ir,1,ist,0) * dzP_res(iz,ir,2,ist,2)   
c     &                -CPpsi(iz,ir,2,ist,0) * dzP_res(iz,ir,1,ist,2)
c     &                +CPdppsi(iz,ir,1,ist) *   P_res(iz,ir,1,ist,2) 
c     &                -CPdppsi(iz,ir,2,ist) *   P_res(iz,ir,2,ist,2) 
c     &                +CPdzpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,2)   
c     &                -CPdzpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,2)
c     &                                      ) + resocr(iz,ir,iq)
cc                  
c                  resocp(iz,ir,iq) = wgtfac0*0.5D0*(
c     &                -CPpsi(iz,ir,1,ist,0) * dzP_res(iz,ir,2,ist,1) 
c     &                -CPpsi(iz,ir,2,ist,0) * dzP_res(iz,ir,1,ist,1) 
c     &                +CPpsi(iz,ir,1,ist,0) * drP_res(iz,ir,1,ist,1)   
c     &                -CPpsi(iz,ir,2,ist,0) * drP_res(iz,ir,2,ist,1)
c     &                +CPdzpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,1) 
c     &                +CPdzpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,1) 
c     &                -CPdrpsi(iz,ir,1,ist) *   P_res(iz,ir,1,ist,1)   
c     &                +CPdrpsi(iz,ir,2,ist) *   P_res(iz,ir,2,ist,1)
cc     &                -CPpsi(iz,ir,1,ist,0) * dzP_res(iz,ir,2,ist,2) 
cc     &                -CPpsi(iz,ir,2,ist,0) * dzP_res(iz,ir,1,ist,2) 
cc     &                +CPpsi(iz,ir,1,ist,0) * drP_res(iz,ir,1,ist,2)   
cc     &                -CPpsi(iz,ir,2,ist,0) * drP_res(iz,ir,2,ist,2)
cc     &                +CPdzpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,2) 
cc     &                +CPdzpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,2) 
cc     &                -CPdrpsi(iz,ir,1,ist) *   P_res(iz,ir,1,ist,2)   
cc     &                +CPdrpsi(iz,ir,2,ist) *   P_res(iz,ir,2,ist,2)
c     &                                      ) + resocp(iz,ir,iq)
cc                  
c                  resocz(iz,ir,iq) = wgtfac*0.5D0*(
c     &                -CPpsi(iz,ir,1,ist,0) * dpP_res(iz,ir,2,ist,1) 
c     &                -CPpsi(iz,ir,2,ist,0) * dpP_res(iz,ir,1,ist,1) 
c     &                -CPpsi(iz,ir,1,ist,0) * drP_res(iz,ir,2,ist,1)   
c     &                +CPpsi(iz,ir,2,ist,0) * drP_res(iz,ir,1,ist,1)
c     &                -CPdppsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,1) 
c     &                -CPdppsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,1) 
c     &                -CPdrpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,1)   
c     &                +CPdrpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,1)
c     &                -CPpsi(iz,ir,1,ist,0) * dpP_res(iz,ir,2,ist,2) 
c     &                -CPpsi(iz,ir,2,ist,0) * dpP_res(iz,ir,1,ist,2) 
c     &                -CPpsi(iz,ir,1,ist,0) * drP_res(iz,ir,2,ist,2)   
c     &                +CPpsi(iz,ir,2,ist,0) * drP_res(iz,ir,1,ist,2)
c     &                -CPdppsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,2) 
c     &                -CPdppsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,2) 
c     &                -CPdrpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist,2)   
c     &                +CPdrpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist,2)
c     &                                      ) + resocz(iz,ir,iq)
cc
c               enddo
c            enddo
c      enddo
c
      endif
#ENDIF
c
c     calculate norm of currents and densities
c     ------------------------------------------------------------------
c      buff(0)=wfovoc(rerho(0,0,1),rerho(0,0,1))
c      buff(1)=wfovoc(rekin(0,0,1),rekin(0,0,1))
c      buff(2)=wfovoc(resocr(0,0,1),resocr(0,0,1))
c      buff(3)=wfovoc(resocp(0,0,1),resocp(0,0,1))
c      buff(4)=wfovoc(resocz(0,0,1),resocz(0,0,1))
c     
c      write(99,'(i,5(g12.3))') 
c     &   imu,buff(0),buff(1),buff(2),buff(3),buff(4)
c
c
c     ------------------------------------------------------------------
c
      return
      end


c======================================================================
c
c     upd_respdens_odd(Q_2q,Q_res,drQ_res,dzQ_res,dpQ_res)
c
c       updates the time-odd response densities resulting from 
c       a small change described by the time-even generator Q with 
c       2q-matrixelements 'Q_2q' and residual action 'Q_res'. 
c
c       imu: mu-component of response fields
c
c=======================================================================
#IF(FULLSPACE)
      subroutine upd_respdens_odd(Q_2q,Q_res,drQ_res,dzQ_res,dpQ_res,
     &                            imu)
#ELSE
      subroutine upd_respdens_odd(Q_2q,dummy1,
     &                             dummy2,dummy3,dummy4,
     &                             imu)
#ENDIF
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
      double precision Q_2q(1:MTD)
#IF(FULLSPACE)
      double precision Q_res(0:MZD,0:MRD,MSD,MPSID)
      double precision drQ_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dzQ_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dpQ_res(0:MZD,0:MRD,MSD,MPSID)
#ELSE
      double precision dummy1,dummy2,dummy3,dummy4
#ENDIF
      double precision wgtfac,wgtfac0,wgtfac2,symfac,multi
      double precision buff(0:10)
      integer          itr,ist,jst,ir,iz,jtp,iq,imu
c
c     ------------------------------------------------------------------
c
      symfac=1.0D0
      if (imu.eq.0) symfac=0.0D0
c
c     set time-odd response densities to zero
c     ------------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
            recurr(iz,ir,1) = zero
            recurr(iz,ir,2) = zero
            recurp(iz,ir,1) = zero
            recurp(iz,ir,2) = zero
            recurz(iz,ir,1) = zero
            recurz(iz,ir,2) = zero
            respnr(iz,ir,1) = zero
            respnr(iz,ir,2) = zero
            respnp(iz,ir,1) = zero
            respnp(iz,ir,2) = zero
            respnz(iz,ir,1) = zero
            respnz(iz,ir,2) = zero
            reskdr(iz,ir,1) = zero
            reskdr(iz,ir,2) = zero
            reskdp(iz,ir,1) = zero
            reskdp(iz,ir,2) = zero
            reskdz(iz,ir,1) = zero
            reskdz(iz,ir,2) = zero
         enddo
      enddo
c                  
c     accumulate time-odd response densities from 2q-matrixelements
c     ------------------------------------------------------------------
      do itr=1,TRlen
c
         ist     = TRidx(itr,1)
         jst     = TRidx(itr,2)
         jtp     = TRidx(itr,3)
         iq      = TRidx(itr,4)
         multi   = TRidx(itr,6)

         wgtfac  = multi*TReta(itr,2)*Q_2q(itr) ! HF-densities
         wgtfac0 = multi*TReta(itr,2)*Q_2q(itr)*symfac
         wgtfac2 = multi* TRxi(itr,2)*Q_2q(itr) ! Pairing-densities
c                  
c        if state jst has positive "T-parity": T|jst>=|Tjst>
c        ---------------------------------------------------------------
         if (jtp.eq.+1) then
c              
            do ir=0,nr
               do iz=0,nz
c 
c                 current
c                 ------------------------------------------------------
                  recurr(iz,ir,iq) = -0.5D0*wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + recurr(iz,ir,iq)

                  recurp(iz,ir,iq) = 0.5D0*wgtfac0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0)  
     &                +CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + recurp(iz,ir,iq)
                  recurz(iz,ir,iq) = -0.5D0*wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + recurz(iz,ir,iq)
c                  
c                 spin density
c                 ------------------------------------------------------
                  respnr(iz,ir,iq) = wgtfac0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                +CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + respnr(iz,ir,iq)
                  
                  respnp(iz,ir,iq) = -wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + respnp(iz,ir,iq)
                  
                  respnz(iz,ir,iq) = wgtfac0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0)  
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + respnz(iz,ir,iq) 
c                  
c                 kinetic spin density
c                 ------------------------------------------------------
                  reskdr(iz,ir,iq) = wgtfac0*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + reskdr(iz,ir,iq)
                  
                  reskdp(iz,ir,iq) = -wgtfac*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)  
     &                -CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + reskdp(iz,ir,iq)
                  
                  reskdz(iz,ir,iq) = wgtfac0*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist)  
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist)
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist)  
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist)
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist)  
     &                -CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)
     &                                      ) + reskdz(iz,ir,iq) 
c
               enddo
            enddo
c                  
c           if state jst has negative "T-parity": T|jst>=-|Tjst>
c           ------------------------------------------------------------
         else 
c              
            do ir=0,nr
               do iz=0,nz
c 
c                 current
c                 ------------------------------------------------------
                  recurr(iz,ir,iq) = -0.5D0*wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + recurr(iz,ir,iq)

                  recurp(iz,ir,iq) = 0.5D0*wgtfac0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + recurp(iz,ir,iq)
                  recurz(iz,ir,iq) = -0.5D0*wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + recurz(iz,ir,iq)
c                  
c                 spin density
c                 ------------------------------------------------------
                  respnr(iz,ir,iq) = wgtfac0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + respnr(iz,ir,iq)
                  
                  respnp(iz,ir,iq) = -wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                -CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + respnp(iz,ir,iq)
                  
                  respnz(iz,ir,iq) = wgtfac0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0)  
     &                -CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + respnz(iz,ir,iq) 
c                  
c                 kinetic spin density
c                 ------------------------------------------------------
                  reskdr(iz,ir,iq) = wgtfac0*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist)
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)  
     &                -CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + reskdr(iz,ir,iq)
                  
                  reskdp(iz,ir,iq) = -wgtfac*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist)
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + reskdp(iz,ir,iq)
                  
                  reskdz(iz,ir,iq) = wgtfac0*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist)  
     &                -CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist)
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist)  
     &                -CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist)
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist)  
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)
     &                                      ) + reskdz(iz,ir,iq) 
c                  
               enddo
            enddo 
c
         endif
c     
      enddo
c 
#IF(FULLSPACE)                 
c     accumulate time-odd response densities from residual-components
c     ------------------------------------------------------------------
      if (isresact.eq.1) then
c      do ist=1,nsttot
cc     
c         wgtfac  = 4.0D0*CPvsqf(ist,0)
c         wgtfac0 = 4.0D0*CPvsqf(ist,0)*symfac
cc
c         if (ist.le.nst(1)) then
c            iq = 1
c         else
c            iq = 2
c         endif
cc
c         do ir=0,nr
c            do iz=0,nz
cc 
cc                 current
cc                 ------------------------------------------------------
c                  recurr(iz,ir,iq) = -0.5D0*wgtfac*(
c     &                +CPpsi(iz,ir,1,ist,0) * drQ_res(iz,ir,1,ist,1) 
c     &                +CPpsi(iz,ir,2,ist,0) * drQ_res(iz,ir,2,ist,1) 
c     &                -CPdrpsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist,1) 
c     &                -CPdrpsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist,1)
c     &                -CPpsi(iz,ir,1,ist,0) * drQ_res(iz,ir,1,ist,2) 
c     &                -CPpsi(iz,ir,2,ist,0) * drQ_res(iz,ir,2,ist,2) 
c     &                +CPdrpsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist,2) 
c     &                +CPdrpsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist,2)
c     &                                      ) + recurr(iz,ir,iq)
c                  
c                  recurp(iz,ir,iq) = 0.5D0*wgtfac0*(
c     &                +CPpsi(iz,ir,1,ist,0) * dpQ_res(iz,ir,1,ist,1) 
c     &                +CPpsi(iz,ir,2,ist,0) * dpQ_res(iz,ir,2,ist,1) 
c     &                +CPdppsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist,1) 
c     &                +CPdppsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist,1)
c     &                +CPpsi(iz,ir,1,ist,0) * dpQ_res(iz,ir,1,ist,2) 
c     &                +CPpsi(iz,ir,2,ist,0) * dpQ_res(iz,ir,2,ist,2) 
c     &                +CPdppsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist,2) 
c     &                +CPdppsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist,2)
c     &                                      ) + recurp(iz,ir,iq)
c                  
c                  recurz(iz,ir,iq) = -0.5D0*wgtfac*(
c     &                +CPpsi(iz,ir,1,ist,0) * dzQ_res(iz,ir,1,ist,1) 
c     &                +CPpsi(iz,ir,2,ist,0) * dzQ_res(iz,ir,2,ist,1) 
c     &                -CPdzpsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist,1) 
c     &                -CPdzpsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist,1)
c     &                -CPpsi(iz,ir,1,ist,0) * dzQ_res(iz,ir,1,ist,2) 
c     &                -CPpsi(iz,ir,2,ist,0) * dzQ_res(iz,ir,2,ist,2) 
c     &                +CPdzpsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist,2) 
c     &                +CPdzpsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist,2)
c     &                                      ) + recurz(iz,ir,iq)
cc                  
cc                 spin density
cc                 ------------------------------------------------------
c                  respnr(iz,ir,iq) = wgtfac0*(
c     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,2,ist,1) 
c     &                +CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,1,ist,1)
c     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,2,ist,2) 
c     &                +CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,1,ist,2)
c     &                                      ) + respnr(iz,ir,iq)
c                  
c                  respnp(iz,ir,iq) = -wgtfac*(
c     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,2,ist,1) 
c     &                -CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,1,ist,1)
c     &                -CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,2,ist,2)
c     &                +CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,1,ist,2)
c     &                                      ) + respnp(iz,ir,iq)
c                  
c                  respnz(iz,ir,iq) = wgtfac0*(
c     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,1,ist,1) 
c     &                -CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,2,ist,1)
c     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,1,ist,2) 
c     &                -CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,2,ist,2)
c     &                                      ) + respnz(iz,ir,iq) 
cc                  
cc                 kinetic spin density
cc                 ------------------------------------------------------
c                  reskdr(iz,ir,iq) = wgtfac0*(
c     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,2,ist,1)  
c     &                +CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,1,ist,1)
c     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,2,ist,1)  
c     &                +CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,1,ist,1)
c     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,2,ist,1)  
c     &                +CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,1,ist,1)
c     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,2,ist,2)  
c     &                +CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,1,ist,2)
c     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,2,ist,2)  
c     &                +CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,1,ist,2)
c     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,2,ist,2)  
c     &                +CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,1,ist,2)
c     &                                      ) + reskdr(iz,ir,iq)
c                  
c                  reskdp(iz,ir,iq) = -wgtfac*(
c     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,2,ist,1)  
c     &                -CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,1,ist,1)
c     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,2,ist,1)  
c     &                -CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,1,ist,1)
c     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,2,ist,1)  
c     &                -CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,1,ist,1)
c     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,2,ist,2)  
c     &                -CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,1,ist,2)
c     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,2,ist,2)  
c     &                -CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,1,ist,2)
c     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,2,ist,2)  
c     &                -CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,1,ist,2)
c     &                                      ) + reskdp(iz,ir,iq)
c                  
c                  reskdz(iz,ir,iq) = wgtfac0*(
c     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,1,ist,1)  
c     &                -CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,2,ist,1)
c     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,1,ist,1)  
c     &                -CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,2,ist,1)
c     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,1,ist,1)  
c     &                -CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,2,ist,1)
c     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,1,ist,2)  
c     &                -CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,2,ist,2)
c     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,1,ist,2)  
c     &                -CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,2,ist,2)
c     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,1,ist,2)  
c     &                -CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,2,ist,2)
c     &                                      ) + reskdz(iz,ir,iq) 
c            enddo
c         enddo
cc
c      enddo
      endif
#ENDIF
c
c     calculate norm of currents and densities
c     ------------------------------------------------------------------
c      buff(0)=wfovoc(recurr(0,0,1),recurr(0,0,1))+
c     &        wfovoc(recurr(0,0,2),recurr(0,0,2))
c      buff(1)=wfovoc(recurp(0,0,1),recurp(0,0,1))+
c     &        wfovoc(recurp(0,0,2),recurp(0,0,2))
c      buff(2)=wfovoc(recurz(0,0,1),recurz(0,0,1))+
c     &        wfovoc(recurz(0,0,2),recurz(0,0,2))
c      buff(3)=wfovoc(respnr(0,0,1),respnr(0,0,1))+
c     &        wfovoc(respnr(0,0,2),respnr(0,0,2))
c      buff(4)=wfovoc(respnp(0,0,1),respnp(0,0,1))+
c     &        wfovoc(respnp(0,0,2),respnp(0,0,2))
c      buff(5)=wfovoc(respnz(0,0,1),respnz(0,0,1))+
c     &        wfovoc(respnz(0,0,2),respnz(0,0,2))
c      buff(6)=wfovoc(reskdr(0,0,1),reskdr(0,0,1))+     
c     &        wfovoc(reskdr(0,0,2),reskdr(0,0,2))
c      buff(7)=wfovoc(reskdp(0,0,1),reskdp(0,0,1))+
c     &        wfovoc(reskdp(0,0,2),reskdp(0,0,2))
c      buff(8)=wfovoc(reskdz(0,0,1),reskdz(0,0,1))+
c     &        wfovoc(reskdz(0,0,2),reskdz(0,0,2))
cc     
c      write(100+imu,'(i,10(g12.3))') 
c     &   imu,CPbeta(0),buff(0),buff(1),buff(2),buff(3),buff(4),
c     &       buff(5),buff(6),buff(7),buff(8)
c
         call wrt_2D_spin(recurr(0,0,1),'recurr____',1,'tst')
         call wrt_2D_spin(recurz(0,0,1),'recurz____',0,'tst')
         call wrt_2D_spin(recurp(0,0,1),'recurp____',0,'tst')
         call wrt_2D_spin(respnr(0,0,1),'respnr____',0,'tst')
         call wrt_2D_spin(respnz(0,0,1),'respnz____',0,'tst')
         call wrt_2D_spin(respnp(0,0,1),'respnp____',0,'tst')
         call wrt_2D_spin(reskdr(0,0,1),'reskdr____',0,'tst')
         call wrt_2D_spin(reskdz(0,0,1),'reskdz____',0,'tst')
         call wrt_2D_spin(reskdp(0,0,1),'reskdp____',0,'tst')
c         call wrt_2D_spin(repairo(0,0,1),'repairo___',0,'tst')

c
c     ------------------------------------------------------------------
c
      return
      end


c======================================================================
c
c     upd_respflds_even()
c
c       calculate the time-even response fields from the time-even 
c       response densities.
c
c       imfld: mu-component of response fields
c
c=======================================================================
      subroutine upd_respflds_even(imfld)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc' ! skyrme parameters
      include 'linresp.inc'
      include 'levels.inc'
c
c     ------------------------------------------------------------------
c
      integer          Gtpar,imfld,zpidx
      double precision buff1,buff2,buff3,buff4,buff5
      double precision buff6,buff7,buff8,buff9,buff10
      double precision buffr1,buffr2,buffr3,buffr4
      double precision buffp1,buffp2,buffp3,buffp4
      double precision buffz1,buffz2,buffz3,buffz4
      double precision b0,b0p,b1,b1p,b2,b2p,b3,b3p,b4,b5,b5p
      double precision wgt_b0,wgt_b0p,wgt_b0c
      double precision wgt_b1,wgt_b1p,wgt_b1c
      double precision wgt_b2,wgt_b2p,wgt_b2c
      double precision wgt_b3,wgt_b3p,wgt_b3c
      double precision wgt_b4,wgt_b4p,wgt_b4c
      double precision wgt_b5,wgt_b5p,wgt_b5c
      double precision wgt_pair1,wgt_pair2
c
c     ------------------------------------------------------------------
c
c     prepare parameters
c     ------------------------------------------------------------------
      b0  =                t0 + 0.5D0*t0*x0 
      b0p =          0.5D0*t0 +       t0*x0 
      b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
      b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
      b2  = 0.125D0*(3.0D0*t1 + 1.5D0*t1*x1 -       t2 - 0.5D0*t2*x2)
      b2p = 0.125D0*(1.5D0*t1 + 3.0D0*t1*x1 + 0.5D0*t2 +       t2*x2)
      b3  = 0.250D0*(      t3 + 0.5D0*t3*x3                         )
      b3p = 0.250D0*(0.5D0*t3 +       t3*x3                         )
      b4  = 0.500D0*       t4
c      b4p = 0.500D0*       t4    ! b4p is a parameter in dynpar.inc
      b5  = 0.125D0*(                 t1*x1 +                  t2*x2)
      b5p = 0.125D0*(      t1               -       t2              )
c
      wgt_b0  = b0-b0p
      wgt_b0c = b0
      wgt_b1  = b1-b1p
      wgt_b1c = b1
      wgt_b2  = -(b2-b2p)
      wgt_b2c = -b2
c
      if (b4p.eq.0.0D0) then
         wgt_b4  = (b4+b4)    ! assume b4p = b4 (W_std)
         wgt_b4c = b4
      else
         wgt_b4  = (b4+b4p)   ! assume b4p != b4 (W_ext)
         wgt_b4c = b4
      endif
c
      if (mod(imfld,2).eq.0) then
         zpidx=1
      else
         zpidx=2
      endif
c                   
c     local potential U
c     -----------------------------------------------------------------
      if (imfld.eq.0) then
         call div_vector(resocr(0,0,1),resocp(0,0,1),resocz(0,0,1),
     &                   LRwrk1_spinor(0,0,1),
     &                   1,imfld,-1, 
     &                   0,1,0, 1.0D0)
         call div_vector(resocr(0,0,2),resocp(0,0,2),resocz(0,0,2),
     &                   LRwrk1_spinor(0,0,2),
     &                   1,imfld,-1, 
     &                   0,1,0, 1.0D0)
      else
         call div_vector(resocr(0,0,1),resocp(0,0,1),resocz(0,0,1),
     &                   LRwrk1_spinor(0,0,1),
     &                   zpidx,imfld,-1, 
     &                   0,0,0, 1.0D0)
         call div_vector(resocr(0,0,2),resocp(0,0,2),resocz(0,0,2),
     &                   LRwrk1_spinor(0,0,2),
     &                   zpidx,imfld,-1, 
     &                   0,0,0, 1.0D0)
      endif
c
      call laplace_scalar(rerho(0,0,1),LRwrk2_spinor(0,0,1),
     &                    zpidx,imfld,1,
     &                    1.0D0)
      call laplace_scalar(rerho(0,0,2),LRwrk2_spinor(0,0,2),
     &                    zpidx,imfld,1,
     &                    1.0D0)

c      call colrun(rerho(0,0,1),recoupot(0,0),0,imfld)
c
c      call wrt_2D_spin(LRwrk1_spinor,'divsoc____',1,'tst')
c      call wrt_2D_spin(LRwrk2_spinor,'laprho____',0,'tst')
c
      do ir=0,nr
         do iz=0,nz
c
            relocpot(iz,ir,1) = wgt_b0              * rerho(iz,ir,1)
     &                        + wgt_b0c             * rerho(iz,ir,2)
     &                        + locpotwgt1(iz,ir,1) * rerho(iz,ir,1)
     &                        + locpotwgt1(iz,ir,2) * rerho(iz,ir,2)
     &                        + wgt_b1              * rekin(iz,ir,1)
     &                        + wgt_b1c             * rekin(iz,ir,2)
     &                        + wgt_b2     *  LRwrk2_spinor(iz,ir,1) 
     &                        + wgt_b2c    *  LRwrk2_spinor(iz,ir,2)
     &                        - wgt_b4     *  LRwrk1_spinor(iz,ir,1)
     &                        - wgt_b4c    *  LRwrk1_spinor(iz,ir,2)
c     &                        - recoupot(iz,ir)
             
            relocpot(iz,ir,2) = wgt_b0              * rerho(iz,ir,2)
     &                        + wgt_b0c             * rerho(iz,ir,1)
     &                        + locpotwgt2(iz,ir,1) * rerho(iz,ir,1)
     &                        + locpotwgt2(iz,ir,2) * rerho(iz,ir,2)
     &                        + wgt_b1              * rekin(iz,ir,2)
     &                        + wgt_b1c             * rekin(iz,ir,1)
     &                        + wgt_b2     *  LRwrk2_spinor(iz,ir,2) 
     &                        + wgt_b2c    *  LRwrk2_spinor(iz,ir,1)
     &                        - wgt_b4     *  LRwrk1_spinor(iz,ir,2)
     &                        - wgt_b4c    *  LRwrk1_spinor(iz,ir,1)
c 
         enddo
      enddo
c                  
c     effective mass B
c     ---------------------------------------------------------------
      do ir=0,nr
         do iz=0,nz
c                
            reeffmas(iz,ir,1) = +wgt_b1  * rerho(iz,ir,1)
     &                          +wgt_b1c * rerho(iz,ir,2)
            reeffmas(iz,ir,2) = +wgt_b1  * rerho(iz,ir,2) 
     &                          +wgt_b1c * rerho(iz,ir,1)
c         
         enddo
      enddo
c                  
c     spin-orbit potential W 
c     ---------------------------------------------------------------
      call nablar_scalar(rerho(0,0,1),LRwrk1_spinor(0,0,1),
     &                   imfld,1,1.0D0) 
      call nablar_scalar(rerho(0,0,2),LRwrk1_spinor(0,0,2),
     &                   imfld,1,1.0D0) 
      call nablaz_scalar(rerho(0,0,1),LRwrk2_spinor(0,0,1),
     &                   zpidx,1,1.0D0)
      call nablaz_scalar(rerho(0,0,2),LRwrk2_spinor(0,0,2),
     &                   zpidx,1,1.0D0)
c
      do ir=0,nr
         do iz=0,nz
c            
            resocpotr(iz,ir,1) = wgt_b4  * LRwrk1_spinor(iz,ir,1)
     &                         + wgt_b4c * LRwrk1_spinor(iz,ir,2)
            resocpotr(iz,ir,2) = wgt_b4  * LRwrk1_spinor(iz,ir,2)
     &                         + wgt_b4c * LRwrk1_spinor(iz,ir,1)
            resocpotz(iz,ir,1) = wgt_b4  * LRwrk2_spinor(iz,ir,1)
     &                         + wgt_b4c * LRwrk2_spinor(iz,ir,2)
            resocpotz(iz,ir,2) = wgt_b4  * LRwrk2_spinor(iz,ir,2)
     &                         + wgt_b4c * LRwrk2_spinor(iz,ir,1)
c
         enddo
      enddo
c
      if (imfld.ne.0) then
c
      call nablap_scalar(rerho(0,0,1),LRwrk3_spinor(0,0,1),
     &                   imfld,+1,1,1.0D0) 
      call nablap_scalar(rerho(0,0,2),LRwrk3_spinor(0,0,2),
     &                   imfld,+1,1,1.0D0)
c
      do ir=0,nr
         do iz=0,nz
c   
            resocpotp(iz,ir,1) = wgt_b4  * LRwrk3_spinor(iz,ir,1)
     &                         + wgt_b4c * LRwrk3_spinor(iz,ir,2)
            resocpotp(iz,ir,2) = wgt_b4  * LRwrk3_spinor(iz,ir,2)
     &                         + wgt_b4c * LRwrk3_spinor(iz,ir,1)
        
         enddo
      enddo
c
      endif
c
c      call wrt_2D_spin(LRwrk1_spinor,'drrho_____',0,'tst')
c      call wrt_2D_spin(LRwrk2_spinor,'dzrho_____',0,'tst')
c      call wrt_2D_spin(LRwrk3_spinor,'dprho_____',0,'tst')
c                  
c     pairing potential
c     ---------------------------------------------------------------
c         wgt_pair1 = -CPpairfc(1,0)
c         wgt_pair2 = -CPpairfc(2,0)
c         do ir=0,nr
c            do iz=0,nz
c
c               buff1 = repair(iz,ir,1)
c               buff2 = repair(iz,ir,2)
c                
c               repairpot(iz,ir,1) = wgt_pair1 * buff1
c               repairpot(iz,ir,2) = wgt_pair2 * buff2
c         
c            enddo
c         enddo
c
c     ------------------------------------------------------------------
c
c      call wrt_2D_spin(relocpot(0,0,1) ,'relocpot__',0,'tst')
c      call wrt_2D_spin(resocpotr(0,0,1),'resocpotr_',0,'tst')
c      call wrt_2D_spin(resocpotz(0,0,1),'resocpotz_',0,'tst')
c      call wrt_2D_spin(resocpotp(0,0,1),'resocpotp_',0,'tst')
c      call wrt_2D_spin(reeffmas(0,0,1), 'reeffmas__',0,'tst')
c
c     ------------------------------------------------------------------
c
c      write(21,'(a)') 'UPD_RESPFLDS(): response fields updated'
c
c     ------------------------------------------------------------------
c            
      return
      end


c======================================================================
c
c     upd_respflds_odd()
c
c       calculate the response fields from the response densities for
c       a small change described by generator G with time-parity Gtpar
c       and m-component imfld. The response densities must be updated
c       before.
c
c=======================================================================
      subroutine upd_respflds_odd(imfld)
c
c     ------------------------------------------------------------------
c      
      include 'paramr.inc'
      include 'dynpar.inc' ! skyrme parameters
      include 'linresp.inc'
      include 'levels.inc'
c
c     ------------------------------------------------------------------
c
      integer          imfld,zpidx
      double precision buff1,buff2,buff3,buff4,buff5
      double precision buff6,buff7,buff8,buff9,buff10
      double precision buffr1,buffr2,buffr3,buffr4
      double precision buffp1,buffp2,buffp3,buffp4
      double precision buffz1,buffz2,buffz3,buffz4
      double precision wgt_b1,wgt_b1p,wgt_b1c
      double precision wgt_b2,wgt_b2p,wgt_b2c
      double precision wgt_b3,wgt_b3p,wgt_b3c
      double precision wgt_b4,wgt_b4p,wgt_b4c
      double precision wgt_b5,wgt_b5p,wgt_b5c
      double precision wgt_pair1,wgt_pair2
c
c     ------------------------------------------------------------------
c
c     prepare parameters
c     ------------------------------------------------------------------
      b0  =                t0 + 0.5D0*t0*x0 
      b0p =          0.5D0*t0 +       t0*x0 
      b1  = 0.250D0*(      t1 + 0.5D0*t1*x1 +       t2 + 0.5D0*t2*x2)
      b1p = 0.250D0*(0.5D0*t1 +       t1*x1 - 0.5D0*t2 -       t2*x2)
      b2  = 0.125D0*(3.0D0*t1 + 1.5D0*t1*x1 -       t2 - 0.5D0*t2*x2)
      b2p = 0.125D0*(1.5D0*t1 + 3.0D0*t1*x1 + 0.5D0*t2 +       t2*x2)
      b3  = 0.250D0*(      t3 + 0.5D0*t3*x3                         )
      b3p = 0.250D0*(0.5D0*t3 +       t3*x3                         )
      b4  = 0.500D0*       t4
c      b4p = 0.500D0*       t4    ! b4p is a parameter in dynpar.inc
      b5  = 0.125D0*(                 t1*x1 +                  t2*x2)
      b5p = 0.125D0*(      t1               -       t2              )
c
      wgt_b1  = -2.0D0 *(b1-b1p)
      wgt_b1c = -2.0D0 * b1
      wgt_b5  = 0.125D0*( t1*x1 + t2*x2 - t1 + t2 )
      wgt_b5c = 0.125D0*( t1*x1 + t2*x2 )
c
      if (b4p.eq.0.0D0) then
         wgt_b4  = -(b4+b4)    ! assume b4p = b4 (W_std)
         wgt_b4c = -b4
      else
         wgt_b4  = -(b4+b4p)   ! assume b4p != b4 (W_ext)
         wgt_b4c = -b4
      endif

c
      if (mod(imfld,2).eq.0) then
         zpidx=1
      else
         zpidx=2
      endif
c                  
c     dynamic potential A
c     ------------------------------------------------------------------
      if (imfld.eq.0) then
         call rotr_vector(respnp(0,0,1),respnz(0,0,1),redynpotr(0,0,1),
     &                    3-zpidx,imfld,-1, 0,1, 1.0D0)
         call rotr_vector(respnp(0,0,2),respnz(0,0,2),redynpotr(0,0,2),
     &                    3-zpidx,imfld,-1, 0,1, 1.0D0)
         call rotp_vector(respnz(0,0,1),respnr(0,0,1),redynpotp(0,0,1),
     &                    3-zpidx,imfld,-1, 1,1, 0.0D0)
         call rotp_vector(respnz(0,0,2),respnr(0,0,2),redynpotp(0,0,2),
     &                    3-zpidx,imfld,-1, 1,1, 0.0D0)
         call rotz_vector(respnr(0,0,1),respnp(0,0,1),redynpotz(0,0,1),
     &                    3-zpidx,imfld,-1, 1,0, 1.0D0)
         call rotz_vector(respnr(0,0,2),respnp(0,0,2),redynpotz(0,0,2),
     &                    3-zpidx,imfld,-1, 1,0, 1.0D0)
      else
         call rotr_vector(respnp(0,0,1),respnz(0,0,1),redynpotr(0,0,1),
     &                    3-zpidx,imfld,-1, 0,0, 1.0D0)
         call rotr_vector(respnp(0,0,2),respnz(0,0,2),redynpotr(0,0,2),
     &                    3-zpidx,imfld,-1, 0,0, 1.0D0)
         call rotp_vector(respnz(0,0,1),respnr(0,0,1),redynpotp(0,0,1),
     &                    3-zpidx,imfld,-1, 0,0, 1.0D0)
         call rotp_vector(respnz(0,0,2),respnr(0,0,2),redynpotp(0,0,2),
     &                    3-zpidx,imfld,-1, 0,0, 1.0D0)
         call rotz_vector(respnr(0,0,1),respnp(0,0,1),redynpotz(0,0,1),
     &                    3-zpidx,imfld,-1, 0,0, 1.0D0)
         call rotz_vector(respnr(0,0,2),respnp(0,0,2),redynpotz(0,0,2),
     &                    3-zpidx,imfld,-1, 0,0, 1.0D0)
      endif
c     
         do ir=0,nr
            do iz=0,nz
c
c              buffr1 = redynpotr(iz,ir,1)
c              buffr2 = redynpotr(iz,ir,2)
c              buffp1 = redynpotp(iz,ir,1)
c              buffp2 = redynpotp(iz,ir,2)
c              buffz1 = redynpotz(iz,ir,1)
c              buffz2 = redynpotz(iz,ir,2)
              buffr3 = recurr(iz,ir,1)
              buffr4 = recurr(iz,ir,2)
              buffp3 = recurp(iz,ir,1)
              buffp4 = recurp(iz,ir,2)
              buffz3 = recurz(iz,ir,1)
              buffz4 = recurz(iz,ir,2)
c                
              redynpotr(iz,ir,1) = wgt_b4  * buffr1 + wgt_b4c  * buffr2
     &                            +wgt_b1  * buffr3 + wgt_b1c  * buffr4
              redynpotr(iz,ir,2) = wgt_b4  * buffr2 + wgt_b4c  * buffr1
     &                            +wgt_b1  * buffr4 + wgt_b1c  * buffr3
              redynpotp(iz,ir,1) = wgt_b4  * buffp1 + wgt_b4c  * buffp2
     &                            +wgt_b1  * buffp3 + wgt_b1c  * buffp4
              redynpotp(iz,ir,2) = wgt_b4  * buffp2 + wgt_b4c  * buffp1
     &                            +wgt_b1  * buffp4 + wgt_b1c  * buffp3
              redynpotz(iz,ir,1) = wgt_b4  * buffz1 + wgt_b4c  * buffz2
     &                            +wgt_b1  * buffz3 + wgt_b1c  * buffz4
              redynpotz(iz,ir,2) = wgt_b4  * buffz2 + wgt_b4c  * buffz1
     &                            +wgt_b1  * buffz4 + wgt_b1c  * buffz3
c         
            enddo
         enddo
c                  
c        spin potential S
c        ---------------------------------------------------------------
         if (imfld.eq.0) then
         call rotr_vector(recurp(0,0,1),recurz(0,0,1),respnpotr(0,0,1),
     &                    zpidx,imfld,-1, 1,0,1.0D0)
         call rotr_vector(recurp(0,0,2),recurz(0,0,2),respnpotr(0,0,2),
     &                    zpdix,imfld,-1, 1,0,1.0D0)
         call rotp_vector(recurz(0,0,1),recurr(0,0,1),respnpotp(0,0,1),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotp_vector(recurz(0,0,2),recurr(0,0,2),respnpotp(0,0,2),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotz_vector(recurr(0,0,1),recurp(0,0,1),respnpotz(0,0,1),
     &                    zpidx,imfld,-1, 0,1,1.0D0)
         call rotz_vector(recurr(0,0,2),recurp(0,0,2),respnpotz(0,0,2),
     &                    zpidx,imfld,-1, 0,1,1.0D0)
         else
         call rotr_vector(recurp(0,0,1),recurz(0,0,1),respnpotr(0,0,1),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotr_vector(recurp(0,0,2),recurz(0,0,2),respnpotr(0,0,2),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotp_vector(recurz(0,0,1),recurr(0,0,1),respnpotp(0,0,1),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotp_vector(recurz(0,0,2),recurr(0,0,2),respnpotp(0,0,2),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotz_vector(recurr(0,0,1),recurp(0,0,1),respnpotz(0,0,1),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         call rotz_vector(recurr(0,0,2),recurp(0,0,2),respnpotz(0,0,2),
     &                    zpidx,imfld,-1, 0,0,1.0D0)
         endif
c
         do ir=0,nr
            do iz=0,nz
c
               buffr1 = respnpotr(iz,ir,1)
               buffr2 = respnpotr(iz,ir,2)
               buffp1 = respnpotp(iz,ir,1)
               buffp2 = respnpotp(iz,ir,2)
               buffz1 = respnpotz(iz,ir,1)
               buffz2 = respnpotz(iz,ir,2)
                
               respnpotr(iz,ir,1) = wgt_b4  * buffr1 + wgt_b4c  * buffr2
               respnpotr(iz,ir,2) = wgt_b4  * buffr2 + wgt_b4c  * buffr1
               respnpotp(iz,ir,1) = wgt_b4  * buffp1 + wgt_b4c  * buffp2
               respnpotp(iz,ir,2) = wgt_b4  * buffp2 + wgt_b4c  * buffp1
               respnpotz(iz,ir,1) = wgt_b4  * buffz1 + wgt_b4c  * buffz2
               respnpotz(iz,ir,2) = wgt_b4  * buffz2 + wgt_b4c  * buffz1
c                   
            enddo
         enddo
c                  
c        effective spin mass C    
c        ---------------------------------------------------------------
         do ir=0,nr
            do iz=0,nz
c
c               buffr1 = respnr(iz,ir,1)
c               buffr2 = respnr(iz,ir,2)
c               buffp1 = respnp(iz,ir,1)
c               buffp2 = respnp(iz,ir,2)
c               buffz1 = respnz(iz,ir,1)
c               buffz2 = respnz(iz,ir,2)
c                
c               respnmasr(iz,ir,1) = wgt_b5  * buffr1 + wgt_b5c  * buffr2
c               respnmasr(iz,ir,2) = wgt_b5  * buffr2 + wgt_b5c  * buffr1
c               respnmasp(iz,ir,1) = wgt_b5  * buffp1 + wgt_b5c  * buffp2
c               respnmasp(iz,ir,2) = wgt_b5  * buffp2 + wgt_b5c  * buffp1
c               respnmasz(iz,ir,1) = wgt_b5  * buffz1 + wgt_b5c  * buffz2
c               respnmasz(iz,ir,2) = wgt_b5  * buffz2 + wgt_b5c  * buffz1
c                   
            enddo
         enddo
c
c     in case of imfld=1 "div A" and "S_z" are singular at r=0. The
c     singularity is of 1. order. We use the r=0 components of the
c     field to store the regularized values "r * div A" and "r * S_z"
c     -----------------------------------------------------------------
      if (imfld.eq.1) then
c       ... to be coded ...         
      endif
c
c     ------------------------------------------------------------------
c
c      write(21,'(a)') 'UPD_RESPFLDS(): response fields updated'
c
c     ------------------------------------------------------------------
c            
      return
      end





c=======================================================================
c
c     upd_resp_H_even(Hresp_2q,Hresp_res)
c
c       updates the matrixelements of the time-even response hamiltonian
c
c     imu: mu-component of response fields
c
c=======================================================================
#IF(FULLSPACE)
      subroutine upd_resp_H_even(Hresp_2q,Hresp_res,imu)
#ELSE
      subroutine upd_resp_H_even(Hresp_2q,dummy,imu)
#ENDIF
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
      double precision Hresp_2q(1:MTD)
      double precision Hresp_res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dummy
c
      integer          itr,ist,jst,ir,iz,jtp,iq,imu,zpidx
      double precision ovlp
c
c     ------------------------------------------------------------------
c
c     the case imu = 0
c     ==================================================================
      if (imu.eq.0) then
c
c        prepare derivatives of mean-fields
c
         call nablar_scalar(reeffmas(0,0,1),LRwrk1_spinor(0,0,1),
     &                      0,1,one) 
         call nablar_scalar(reeffmas(0,0,2),LRwrk1_spinor(0,0,2),
     &                      0,1,one) 
         call nablaz_scalar(reeffmas(0,0,1),LRwrk2_spinor(0,0,1),
     &                      1,1,one)
         call nablaz_scalar(reeffmas(0,0,2),LRwrk2_spinor(0,0,2),
     &                      1,1,one)
c
c         call wrt_2D_spin(LRwrk2_spinor,'dzeffmas__',1,'tst')
c         call wrt_2D_spin(LRwrk1_spinor,'dreffmas__',0,'tst')
c
c        calculate the action of response Hamiltonian on sp-states
c        (the second appearence of reeffmas and resocpotz are just 
c         dummy fields because phi-derivative of reeffmas and 
c         resocpotp vanish in case of axial symmetry)
c
      call act_mfh_even(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,1),
     &            relocpot(0,0,1),                         1,
     &            reeffmas,LRwrk1_spinor,LRwrk2_spinor,
     &            reeffmas,                                1,
     &            resocpotr,resocpotz,resocpotz,           1,
     &            CPconstr(0,0,0),CPLambda(0),             0,
     &            CPspenrg(1,0),                           0,
     &            0,1.0D0,
     &            1,nst(1)+nst(2)) 
c
c     the case imu != 0
c     ==================================================================
      else
c
c        prepare derivatives of mean-fields
c
c         zpidx=1+mod(imu,2)
c
c         call nablar_scalar(reeffmas(0,0,1),LRwrk1_spinor(0,0,1),
c     &                      imu,1,one) 
c         call nablar_scalar(reeffmas(0,0,2),LRwrk1_spinor(0,0,2),
c     &                      imu,1,one) 
c         call nablaz_scalar(reeffmas(0,0,1),LRwrk2_spinor(0,0,1),
c     &                      zpidx,1,one)
c         call nablaz_scalar(reeffmas(0,0,2),LRwrk2_spinor(0,0,2),
c     &                      zpidx,1,one)
c         call nablap_scalar(reeffmas(0,0,1),LRwrk3_spinor(0,0,1),
c     &                      imu,+1,1,one) 
c         call nablap_scalar(reeffmas(0,0,2),LRwrk3_spinor(0,0,2),
c     &                      imu,+1,1,one) 
c
c         call wrt_2D_spin(LRwrk3_spinor,'dpeffmas__',1,'tst')
c         call wrt_2D_spin(LRwrk2_spinor,'dzeffmas__',0,'tst')
c         call wrt_2D_spin(LRwrk1_spinor,'dreffmas__',0,'tst')
c
c        calculate the action of response Hamiltonian on sp-states
c
c         call act_mfh_even(
c     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
c     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
c     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,1),
c     &            relocpot(0,0,1),                          1,
c     &            reeffmas,LRwrk1_spinor,LRwrk2_spinor,
c     &            LRwrk3_spinor,                            1,
c     &            resocpotr,resocpotz,resocpotp,            1,
c     &            CPconstr(0,0,0),CPlambda(0),              0,
c     &            CPspenrg(1,0),                            0,
c     &            imu,0.5D0)  
cc
      endif
c
c     project qp-components
c     ==================================================================
      call project_qp(Hresp_2q,Hresp_res,+1,imu,0)
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     upd_resp_H_odd(Hresp_2q,Hresp_res,imu)
c
c       updates the matrixelements of the time-odd response hamiltonian
c
c     imu: mu-component of response fields
c
c=======================================================================
#IF(FULLSPACE)
      subroutine upd_resp_H_odd(Hresp_2q,Hresp_res,imu)
#ELSE
      subroutine upd_resp_H_odd(Hresp_2q,dummy,imu)
#ENDIF
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
      double precision Hresp_2q(1:MTD)
      double precision Hresp_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dummy
c
      integer          itr,ist,jst,ir,iz,jtp,iq,imu
      integer          izpidx
c
c     ------------------------------------------------------------------
c
c     the case imu = 0
c     ==================================================================
      if (imu.eq.0) then
c
c        prepare derivatives of mean-fields
c
         call div_vector(redynpotr(0,0,1),redynpotp(0,0,1),
     &                   redynpotz(0,0,1),LRwrk1_spinor(0,0,1),
     &                   1,0,-1,
     &                   0,1,0,1.0D0)
         call div_vector(redynpotr(0,0,2),redynpotp(0,0,2),
     &                   redynpotz(0,0,2),LRwrk1_spinor(0,0,2),
     &                   1,0,-1,
     &                   0,1,0,1.0D0)
c
         call wrt_2D_spin(LRwrk1_spinor(0,0,1),'divdyn____',1,'tst')
         call wrt_2D_spin(redynpotr(0,0,1),'dynpotr___',0,'tst')
         call wrt_2D_spin(redynpotz(0,0,1),'dynpotz___',0,'tst')
         call wrt_2D_spin(redynpotp(0,0,1),'dynpotp___',0,'tst')
         call wrt_2D_spin(respnpotr(0,0,1),'spnpotr___',0,'tst')
         call wrt_2D_spin(respnpotz(0,0,1),'spnpotz___',0,'tst')
         call wrt_2D_spin(respnpotp(0,0,1),'spnpotp___',0,'tst')
c
c        calculate the action of response Hamiltonian on sp-states
c
         call act_mfh_odd(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1),
     &            redynpotr,redynpotz,redynpotp,LRwrk1_spinor, 1,
     &            respnpotr,respnpotz,respnpotp,               1,
     &            respnmasr,respnmasz,respnmasp,               0,
     &            0,1.0D0)  
c
c     the case imu != 0
c     ==================================================================
      else
c
c        prepare derivatives of mean-fields
c
         izpidx =  1 + mod(abs(imu),2)
c
         call div_vector(redynpotr(0,0,1),redynpotp(0,0,1),
     &                   redynpotz(0,0,1),LRwrk1_spinor(0,0,1),
     &                   izpidx,imu,-1,
     &                   0,0,0,1.0D0)
         call div_vector(redynpotr(0,0,2),redynpotp(0,0,2),
     &                   redynpotz(0,0,2),LRwrk1_spinor(0,0,2),
     &                   izpidx,imu,-1,
     &                   0,0,0,1.0D0)
c
         call wrt_2D_spin(LRwrk1_spinor(0,0,1),'divdyn____',1,'ts1')
         call wrt_2D_spin(redynpotr(0,0,1),'dynpotr___',0,'ts1')
         call wrt_2D_spin(redynpotz(0,0,1),'dynpotz___',0,'ts1')
         call wrt_2D_spin(redynpotp(0,0,1),'dynpotp___',0,'ts1')
         call wrt_2D_spin(respnpotr(0,0,1),'spnpotr___',0,'ts1')
         call wrt_2D_spin(respnpotz(0,0,1),'spnpotz___',0,'ts1')
         call wrt_2D_spin(respnpotp(0,0,1),'spnpotp___',0,'ts1')
c
c      NOTE: the div A and S_z are singular at r=0 for mu=1. 
c            This is OK. In case of mu=1 at least one state 
c            has |m|>1 and therefore a term propto r appears 
c            in all matrixelements. So the total integration 
c            kernel keeps regular. Nevertheless
c            the value at r=0 is numerically lost and has to be
c            recovered by l.Hospital or by using the alternative
c            representation of
c
c            1/2 {nabla, A} = 1/2 div A + A*nabla = 
c                                   <--              -->
c                             -1/2 nabla A + 1/2  A nabla
c
c            which is only applicable in the calculation of
c            matrixelements.
c
c        calculate the action of response Hamiltonian on sp-states
c
c
         call act_mfh_odd(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1),
     &            redynpotr,redynpotz,redynpotp,LRwrk1_spinor, 1,
     &            respnpotr,respnpotz,respnpotp,               1,
     &            respnmasr,respnmasz,respnmasp,               0,
     &            imu,0.5D0)  
c
      endif
c
c     project qp-components
c     ==================================================================
      call project_qp(Hresp_2q,Hresp_res,-1,imu,0)
c
c     ------------------------------------------------------------------
c
      return
      end


c======================================================================
c
c     ini_damping(E0,X0)      
c
c     initialize the damping matrices for the residual space iteration
c     
c=======================================================================
      subroutine ini_damping(einv,Dampr,Dampz)
c
      include 'paramr.inc' ! one, two
      include 'dynpar.inc' ! h2meff
      include 'derivs.inc'
c
      double precision einv
      double precision Dampr(0:NRD,0:NRD,0:kmpsi)
      double precision Dampz(0:NZD,0:NZD,2)
c
      integer ipi,i,j,impsi
      double precision befac
      double precision zaux(0:nzd,0:nzd),raux(0:nrd,0:nrd)
c
c      befac  = h2meff/(two*einv)
      befac  = h2mav/(two*einv)
c
c     kinetic energy damping in z-direction
c     ---------------------------------------------------------------
      do ipi=1,2
c
        do i=0,nz
          do j=0,nz
            Dampz(j,i,ipi) = -befac * derz2f(j,i,ipi)
          enddo
          Dampz(i,i,ipi) = one + Dampz(i,i,ipi)
        enddo
c
c       inversion yields final result for zdnuc
        call minv(Dampz(0,0,ipi),zaux,nzd+1,nz+1)
c
c       ensure correct boundary values
        if (ipi.eq.2) then
          do i=0,nz
            Dampz(0,i,ipi) = 0.0D0
            Dampz(i,0,ipi) = 0.0D0  
          enddo
        endif
c
      enddo
c
c     kinetic energy damping in r-direction
c     ---------------------------------------------------------------
      do impsi=0,kmpsi
c
         do j=0,nr
            do i=0,nr
               Dampr(i,j,impsi) = -befac * derr2af(i,j,impsi,1)
            enddo
         enddo
         do i=0,nr
           Dampr(i,i,impsi) =  one + Dampr(i,i,impsi)
         enddo
c
         call minv(Dampr(0,0,impsi),raux,nrd+1,nr+1)
c
         if (impsi.gt.0) then
            do i=0,nr
             Dampr(0,i,impsi) = zero
             Dampr(i,0,impsi) = zero 
           enddo
         endif
c
      enddo
c
c     ------------------------------------------------------------------
c
      write(21,'(a,g12.4)') 
     &  'INI_DAMPING() : damping matrices initialized with E0=',einv
c
c     ------------------------------------------------------------------
c
      return
      end


c======================================================================
c
c     wrt_response()
c
c     print responses
c
c=======================================================================
      subroutine wrt_response()
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'
c
      call wrt_2D_spin(rerho(0,0,1) ,'rerho_____',1,'ts3') ! OK for m=0
      call wrt_2D_spin(rekin(0,0,1) ,'rekin_____',0,'ts3') ! OK for m=0
      call wrt_2D_spin(resocz(0,0,1),'resocz____',0,'ts3') ! OK for m=0
      call wrt_2D_spin(resocr(0,0,1),'resocr____',0,'ts3') ! OK for m=0
      call wrt_2D_spin(resocp(0,0,1),'resocp____',0,'ts3') 
      call wrt_2D_spin(recurr(0,0,1),'recurr____',0,'ts3') ! OK for m=0
      call wrt_2D_spin(recurp(0,0,1),'recurp____',0,'ts3') ! OK for m=0
      call wrt_2D_spin(recurz(0,0,1),'recurz____',0,'ts3')
      call wrt_2D_spin(respnr(0,0,1),'respnr____',0,'ts3')
      call wrt_2D_spin(respnp(0,0,1),'respnp____',0,'ts3')
      call wrt_2D_spin(respnz(0,0,1),'respnz____',0,'ts3')
      call wrt_2D_spin(reskdr(0,0,1),'reskdr____',0,'ts3')
      call wrt_2D_spin(reskdp(0,0,1),'reskdp____',0,'ts3')
      call wrt_2D_spin(reskdz(0,0,1),'reskdz____',0,'ts3')
      call wrt_2D_spin(relocpot(0,0,1) ,'relocpot__',0,'ts3') ! ? 
      call wrt_2D_spin(reeffmas(0,0,1) ,'reeffmas__',0,'ts3') ! OK for m=0
      call wrt_2D_spin(redynpotr(0,0,1),'redynpotr_',0,'ts3') ! OK for m=0
      call wrt_2D_spin(redynpotz(0,0,1),'redynpotz_',0,'ts3') ! OK for m=0
      call wrt_2D_spin(redynpotp(0,0,1),'redynpotp_',0,'ts3')
      call wrt_2D_spin(respnpotr(0,0,1),'respnpotr_',0,'ts3')
      call wrt_2D_spin(respnpotz(0,0,1),'respnpotz_',0,'ts3')
      call wrt_2D_spin(respnpotp(0,0,1),'respnpotp_',0,'ts3')
      call wrt_2D_spin(resocpotr(0,0,1),'resocpotr_',0,'ts3') ! OK for m=0
      call wrt_2D_spin(resocpotz(0,0,1),'resocpotz_',0,'ts3') ! OK for m=0
      call wrt_2D_spin(resocpotp(0,0,1),'resocpotp_',0,'ts3')
      call wrt_2D_spin(respnmasr(0,0,1),'respnmasr_',0,'ts3')
      call wrt_2D_spin(respnmasz(0,0,1),'respnmasz_',0,'ts3')
      call wrt_2D_spin(respnmasp(0,0,1),'respnmasp_',0,'ts3')
c
      return
      end


c======================================================================
c
c     comm_and_var_QP(imu,BVar,PVar,CommBP)
c
c     memory saving evaluation of fullspace <B^2>-B^2, <P^2> and <[B,P]>
c
c=======================================================================
      subroutine comm_and_var_QP(BVar,PVar,CommBP,beta0,trzpar,imu,iiso)
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'linresp.inc'
      include 'levels.inc'   !nst
c
      integer iiso,imu,trzpar
      real*8  BVar,PVar,CommBP,beta0
c
      real*8  normfac,qlmfac,pfac
      real*8  lastv,lastu,nextv,nextu,rcntv,rcntu
      real*8  ui,vi,uj,vj,peta,qeta,povlp,qovlp
      real*8  QlmOper(0:MZD,0:MRD)
      real*8  Qbuff(0:MZD,0:MRD,MSD,2)
      real*8  Pbuff(0:MZD,0:MRD,MSD,2)
      real*8  r,z
      integer ir,iz,ist,jst,jstmin,jstmax
c
c     buffered SQRT as internal statement function
c     ------------------------------------------------------------------
      sqrtb(x) = dsqrt(max(x,smal))
c
c     zero return values
c     ------------------------------------------------------------------
      BVar  =0.0D0
      PVar  =0.0D0
      CommBP=0.0D0
c
c     calc. normalization of Q_mu for given isosp. (0: total 1:p 2:n)
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
      end if
c        
c     initialize multipole operator and scale of momentum operators
c     ------------------------------------------------------------------
      if (imu.eq.0) then
         qlmfac   = 0.8D0*pi/normfac
         pfac     = 1.0D0/(CPbeta(+1)-CPbeta(-1))
      else
         qlmfac   = 0.8D0*pi/(sqrt(2.0D0)*normfac)
         if (beta0.ne.0.0D0) then
            pfac     = 1.0D0/(sqrt(3.0D0)*beta0)
         else
            pfac     = 0.0D0
         endif
      endif
c
      do ir=0,nr
         r    = dr*(ir+sr)
         do iz=0,nz
            z    = ((iz)+sz-icm)*dz - CPzcm(0)
            QlmOper(iz,ir) = qlmfac*ylm(z,r,2,imu) 
         enddo
      enddo
c
c     operators with imu = 0 
c     ------------------------------------------------------------------
      if (imu.eq.0) then
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
         do ist=1,nst(1)+nst(2)
c     
c           calculate action of Q and P operator for mu=0
c           ------------------------------------------------------------
            do ir=0,nr
               do iz=0,nz
                  Qbuff(iz,ir,1,1)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  Qbuff(iz,ir,2,1)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
                  Pbuff(iz,ir,1,1)=pfac*
     &               (CPpsi(iz,ir,1,ist,+1)-CPpsi(iz,ir,1,ist,-1))
                  Pbuff(iz,ir,2,1)=pfac*
     &               (CPpsi(iz,ir,2,ist,+1)-CPpsi(iz,ir,2,ist,-1))
               enddo
            enddo
c
c           evaluate 2qp-matrixelements
c           ------------------------------------------------------------
            if (ist.le.nst(1)) then
               jstmin=1
               jstmax=nst(1)
            else
               jstmin=nst(1)+1
               jstmax=nst(1)+nst(2)
            endif
c
            do jst=jstmin,jstmax
               if (mpsi(jst)-mpsi(ist).eq.0)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     if (ist.ne.jst) then
c
                        povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1)) 
                        call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1),
     &                                                          povlp)
c
                        qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
                        call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
     &                                                          qovlp)
c
                        BVar   = BVar   + (qeta*qovlp)**2
                        PVar   = PVar   + (peta*povlp)**2
                        CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                     else
c
                        lastv   = sqrtb(CPvsqf(ist,-1))
                        lastu   = sqrtb(one-CPvsqf(ist,-1))
                        rcntv   = sqrtb(CPvsqf(ist,0)) 
                        rcntu   = sqrtb(one-CPvsqf(ist,0))
                        nextv   = sqrtb(CPvsqf(ist,+1))
                        nextu   = sqrtb(one-CPvsqf(ist,+1))
c
                        if (rcntv.ge.0.5D0) then
                           povlp =  -(nextu-lastu)/rcntv*pfac
                        else
                           povlp =  +(nextv-lastv)/rcntu*pfac
                        endif
c
                        qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
                        call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
     &                                                          qovlp)
c
                        BVar   = BVar   +       (qeta*qovlp)**2
                        PVar   = PVar   +       (povlp)**2
                        CommBP = CommBP + 2.0D0*qeta*povlp*qovlp
c
                     endif  
c
                  endif
               endif
c     
            enddo
c
c           project out diagonal components of momentum operator
c           ------------------------------------------------------------
            ovlp=wfovlp(CPpsi(0,0,1,ist,0),Pbuff)
            call wfproj(CPpsi(0,0,1,ist,0),Pbuff,ovlp)
c     
c           accumulate residual components
c           ------------------------------------------------------------
            BVar   = BVar   + 2.0D0*CPvsqf(ist,0)*
     &                        wfovlp(Qbuff(0,0,1,1) ,Qbuff(0,0,1,1))
            PVar   = PVar   + 2.0D0*CPvsqf(ist,0)*
     &                        wfovlp(Pbuff(0,0,1,1) ,Pbuff(0,0,1,1))
            CommBP = CommBP + 4.0D0*CPvsqf(ist,0)*
     &                        wfovlp(Qbuff(0,0,1,1) ,Pbuff(0,0,1,1))            
c
         enddo
c
c     operators with imu = 1
c     ------------------------------------------------------------------
      else if (imu.eq.1) then
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
         do ist=1,nst(1)+nst(2)
c     
c           calculate action of Q and P operator for mu=1
c           ------------------------------------------------------------
            do ir=0,nr
               r = (ir+sr)*dr
               do iz=0,nz
                  z = (iz+sz-icm)*dz - CPzcm(0)
c
                  Qbuff(iz,ir,1,1)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  Qbuff(iz,ir,2,1)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
                  Qbuff(iz,ir,1,2)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  Qbuff(iz,ir,2,2)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
c     
c                 ! upward laddering components of J_y
                  Pbuff(iz,ir,1,1) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,1,ist) 
     &                                  +z*CPdppsi(iz,ir,1,ist) 
     &                                  +r*CPdzpsi(iz,ir,1,ist)
     &                                  -  CPpsi(iz,ir,2,ist,0))
                  Pbuff(iz,ir,2,1) =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,2,ist) 
     &                                  +z*CPdppsi(iz,ir,2,ist) 
     &                                  +r*CPdzpsi(iz,ir,2,ist))
c
c                ! downward laddering components of J_y
                 Pbuff(iz,ir,1,2)  =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,1,ist) 
     &                                  -z*CPdppsi(iz,ir,1,ist) 
     &                                  +r*CPdzpsi(iz,ir,1,ist))
                 Pbuff(iz,ir,2,2)  =  0.5D0*pfac*(
     &                                  -z*CPdrpsi(iz,ir,2,ist) 
     &                                  -z*CPdppsi(iz,ir,2,ist) 
     &                                  +r*CPdzpsi(iz,ir,2,ist)
     &                                  +  CPpsi(iz,ir,1,ist,0))
               enddo
            enddo
c
c           evaluate 2qp-matrixelements
c           ------------------------------------------------------------
            if (ist.le.nst(1)) then
               jstmin=1
               jstmax=nst(1)
            else
               jstmin=nst(1)+1
               jstmax=nst(1)+nst(2)
            endif
c
            do jst=jstmin,jstmax
c
c              first the upward transitions
               if (mpsi(jst)-mpsi(ist).eq.+1)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1)) 
                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1),
     &                                                          povlp)
c
                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
     &                                                          qovlp)
                     BVar   = BVar   + (qeta*qovlp)**2
                     PVar   = PVar   + (peta*povlp)**2
                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                  endif
               endif
c
c              now the downward transitions
               if (mpsi(jst)-mpsi(ist).eq.-1)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
     &                                                          povlp)
c
                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
     &                                                          qovlp)
                     BVar   = BVar   + (qeta*qovlp)**2
                     PVar   = PVar   + (peta*povlp)**2
                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                  endif
               endif
c     
            enddo
c
c           and finally the downward transitions into the T-mirrored
c           states
            do jst=jstmin,jstmax
c
               if (mpsi(ist)+mpsi(jst)+1.eq.+1)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.+1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     povlp=wfovlptr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
                     call wfprojtr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
     &                                                          povlp)
c
                     qovlp=wfovlptr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
                     call wfprojtr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
     &                                                          qovlp)
                     BVar   = BVar   + (qeta*qovlp)**2
                     PVar   = PVar   + (peta*povlp)**2
                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                  endif
               endif
c     
            enddo
c     
c           accumulate residual components
c           ------------------------------------------------------------
            BVar   = BVar   + 2.0D0*CPvsqf(ist,0)*(
     &                        wfovlp(Qbuff(0,0,1,1) ,Qbuff(0,0,1,1))
     &                       +wfovlp(Qbuff(0,0,1,2) ,Qbuff(0,0,1,2)))
            PVar   = PVar   + 2.0D0*CPvsqf(ist,0)*(
     &                        wfovlp(Pbuff(0,0,1,1) ,Pbuff(0,0,1,1))
     &                       +wfovlp(Pbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
            CommBP = CommBP + 4.0D0*CPvsqf(ist,0)*(
     &                        wfovlp(Qbuff(0,0,1,1) ,Pbuff(0,0,1,1))
     &                       +wfovlp(Qbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
c
         enddo
c
c     operators with imu = 2
c     ------------------------------------------------------------------
      else if (imu.eq.2) then
c        
c     loop over all sp-states 
c     ------------------------------------------------------------------
         do ist=1,nst(1)+nst(2)
c     
c           calculate action of Q and P operator for mu=1
c           ------------------------------------------------------------
            do ir=0,nr
               r = (ir+sr)*dr
               do iz=0,nz
                  z = (iz+sz-icm)*dz - CPzcm(0)
c
                  Qbuff(iz,ir,1,1)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  Qbuff(iz,ir,2,1)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
                  Qbuff(iz,ir,1,2)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
                  Qbuff(iz,ir,2,2)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
c     
                  Pbuff(iz,ir,1,1)  =  0.0D0
                  Pbuff(iz,ir,2,1)  =  0.0D0
                  Pbuff(iz,ir,1,2)  =  0.0D0
                  Pbuff(iz,ir,2,2)  =  0.0D0
               enddo
            enddo
c
c           evaluate 2qp-matrixelements
c           ------------------------------------------------------------
            if (ist.le.nst(1)) then
               jstmin=1
               jstmax=nst(1)
            else
               jstmin=nst(1)+1
               jstmax=nst(1)+nst(2)
            endif
c
            do jst=jstmin,jstmax
c
c              first the upward transitions
               if (mpsi(jst)-mpsi(ist).eq.+2)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1)) 
                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1),
     &                                                          povlp)
c
                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
     &                                                          qovlp)
                     BVar   = BVar   + (qeta*qovlp)**2
                     PVar   = PVar   + (peta*povlp)**2
                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                  endif
               endif
c
c              now the downward transitions
               if (mpsi(jst)-mpsi(ist).eq.-2)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
     &                                                          povlp)
c
                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
     &                                                          qovlp)
                     BVar   = BVar   + (qeta*qovlp)**2
                     PVar   = PVar   + (peta*povlp)**2
                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                  endif
               endif
c     
            enddo
c
c           and finally the downward transitions into the T-mirrored
c           states
            do jst=jstmin,jstmax
c
               if (mpsi(ist)+mpsi(jst)+1.eq.+2)  then
                  if (npar(jst)*npar(ist)*trzpar*isym.ne.+1) then
c
                     vi   = sqrtb(CPvsqf(ist,0))  
                     ui   = sqrtb(one-vi*vi)
                     vj   = sqrtb(CPvsqf(jst,0))
                     uj   = sqrtb(one-vj*vj)
                     qeta = uj*vi + vj*ui
                     peta = uj*vi - vj*ui
c
                     povlp=wfovlptr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
                     call wfprojtr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
     &                                                          povlp)
c
                     qovlp=wfovlptr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
                     call wfprojtr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
     &                                                          qovlp)
                     BVar   = BVar   + (qeta*qovlp)**2
                     PVar   = PVar   + (peta*povlp)**2
                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
c
                  endif
               endif
c     
            enddo
c     
c           accumulate residual components
c           ------------------------------------------------------------
            BVar   = BVar   + 2.0D0*CPvsqf(ist,0)*(
     &                        wfovlp(Qbuff(0,0,1,1) ,Qbuff(0,0,1,1))
     &                       +wfovlp(Qbuff(0,0,1,2) ,Qbuff(0,0,1,2)))
            PVar   = PVar   + 2.0D0*CPvsqf(ist,0)*(
     &                        wfovlp(Pbuff(0,0,1,1) ,Pbuff(0,0,1,1))
     &                       +wfovlp(Pbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
            CommBP = CommBP + 4.0D0*CPvsqf(ist,0)*(
     &                        wfovlp(Qbuff(0,0,1,1) ,Pbuff(0,0,1,1))
     &                       +wfovlp(Qbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
c
         enddo
c
      endif
c
      return
      end


cc======================================================================
cc
cc     GCM_mass(imu,BVar,PVar,CommBP)
cc
cc     memory saving evaluation of fullspace GCM-mass
cc
cc=======================================================================
c      subroutine GCM_mass(mass,beta0,trzpar,imu,iiso)
cc
c      include 'paramr.inc'
c      include 'dynpar.inc'
c      include 'linresp.inc'
c      include 'levels.inc'   !nst
cc
c      integer iiso,imu,trzpar
c      real*8  mass,beta0
cc
c      real*8  normfac,qlmfac,pfac
c      real*8  lastv,lastu,nextv,nextu,rcntv,rcntu
c      real*8  ui,vi,uj,vj,peta,qeta,povlp,qovlp
c      real*8  HPbuff(0:MZD,0:MRD,MSD,2)
c      real*8  Pbuff(0:MZD,0:MRD,MSD,2)
c      real*8  drPbuff(0:MZD,0:MRD,MSD,2)
c      real*8  dzPbuff(0:MZD,0:MRD,MSD,2)
c      real*8  dpPbuff(0:MZD,0:MRD,MSD,2)
c      real*8  ddPbuff(0:MZD,0:MRD,MSD,2)
c      real*8  r,z
c      integer ir,iz,ist,jst,jstmin,jstmax
cc
cc     buffered SQRT as internal statement function
cc     ------------------------------------------------------------------
c      sqrtb(x) = dsqrt(max(x,smal))
cc
cc     zero return value
cc     ------------------------------------------------------------------
c      mass  =0.0D0
cc
c      if (imu.eq.0) then
c         qlmfac   = 0.8D0*pi/normfac
c         pfac     = 1.0D0/(CPbeta(+1)-CPbeta(-1))
c      else
c         qlmfac   = 0.8D0*pi/(sqrt(2.0D0)*normfac)
c         if (beta0.ne.0.0D0) then
c            pfac     = 1.0D0/(sqrt(3.0D0)*beta0)
c         else
c            pfac     = 0.0D0
c         endif
c      endif
cc
cc     operators with imu = 0 
cc     ------------------------------------------------------------------
c      if (imu.eq.0) then
cc        
cc     loop over all sp-states 
cc     ------------------------------------------------------------------
c         do ist=1,nst(1)+nst(2)
cc     
cc           calculate action of P operator for mu=0
cc           ------------------------------------------------------------
c            do ir=0,nr
c               do iz=0,nz
c                  Pbuff(iz,ir,1,1)=pfac*
c     &               (CPpsi(iz,ir,1,ist,+1)-CPpsi(iz,ir,1,ist,-1))
c                  Pbuff(iz,ir,2,1)=pfac*
c     &               (CPpsi(iz,ir,2,ist,+1)-CPpsi(iz,ir,2,ist,-1))
c                  HPbuff(iz,ir,1,1)=pfac*
c     &               (CPpsi(iz,ir,1,ist,+1)-CPpsi(iz,ir,1,ist,-1))
c                  HPbuff(iz,ir,2,1)=pfac*
c     &               (CPpsi(iz,ir,2,ist,+1)-CPpsi(iz,ir,2,ist,-1))
c               enddo
c            enddo
cc
cc           evaluate 2qp-matrixelements
cc           ------------------------------------------------------------
c            if (ist.le.nst(1)) then
c               jstmin=1
c               jstmax=nst(1)
c            else
c               jstmin=nst(1)+1
c               jstmax=nst(1)+nst(2)
c            endif
cc
c            do jst=jstmin,jstmax
c               if (mpsi(jst)-mpsi(ist).eq.0)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     if (ist.ne.jst) then
cc
c                        povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1)) 
c                        call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1),
c     &                                                          povlp)
cc
c                        qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
c                        call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
c     &                                                          qovlp)
cc
c                        BVar   = BVar   + (qeta*qovlp)**2
c                        PVar   = PVar   + (peta*povlp)**2
c                        CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                     else
cc
c                        lastv   = sqrtb(CPvsqf(ist,-1))
c                        lastu   = sqrtb(one-CPvsqf(ist,-1))
c                        rcntv   = sqrtb(CPvsqf(ist,0)) 
c                        rcntu   = sqrtb(one-CPvsqf(ist,0))
c                        nextv   = sqrtb(CPvsqf(ist,+1))
c                        nextu   = sqrtb(one-CPvsqf(ist,+1))
cc
c                        if (rcntv.ge.0.5D0) then
c                           povlp =  -(nextu-lastu)/rcntv*pfac
c                        else
c                           povlp =  +(nextv-lastv)/rcntu*pfac
c                        endif
cc
c                        qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
c                        call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
c     &                                                          qovlp)
cc
c                        BVar   = BVar   +       (qeta*qovlp)**2
c                        PVar   = PVar   +       (povlp)**2
c                        CommBP = CommBP + 2.0D0*qeta*povlp*qovlp
cc
c                     endif  
cc
c                  endif
c               endif
cc     
c            enddo
cc
cc           project out diagonal components of momentum operator
cc           ------------------------------------------------------------
c            ovlp=wfovlp(CPpsi(0,0,1,ist,0),Pbuff)
c            call wfproj(CPpsi(0,0,1,ist,0),Pbuff,ovlp)
cc     
cc           update derivatives of P|ist>
cc           ------------------------------------------------------------
c            call upd_derive_spinor_range(Pbuff(0,0,1,1),
c     &                       drPbuff(0,0,1,1),dzPbuff(0,0,1,1),
c     &                       dpPbuff(0,0,1,1),ddPbuff(0,0,1,1),
c     &                       (-1)**imu,imu,1,1)
cc     
cc           calculate action of H operator on P|ist> 
cc           ------------------------------------------------------------
c            call act_mfh_even(Pbuff(0,0,1,1),drPbuff(0,0,1,1),
c     &                  dzPbuff(0,0,1,1),dpPbuff(0,0,1,1),
c     &                  ddPbuff(0,0,1,1),HPbuff(0,0,1,1),
c     &                  CPlocpot(0,0,1,0),                        1,
c     &                  CPeffmas(0,0,1,0),CPdreffmas(0,0,1),
c     &                  CPdzeffmas(0,0,1),CPdzeffmas(0,0,1),      1,
c     &                  CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
c     &                  CPsocpotr(0,0,1,0),                       1,
c     &                  CPconstr(0,0,0),CPlambda(0),              1,
c     &                  CPspenrg(1,0),                            0,
c     &                  imu,1.0D0,
c     &                  1,1)
cc     
cc           accumulate residual components
cc           ------------------------------------------------------------
c            BVar   = BVar   + 2.0D0*CPvsqf(ist,0)*
c     &                        wfovlp(Qbuff(0,0,1,1) ,Qbuff(0,0,1,1))
c            PVar   = PVar   + 2.0D0*CPvsqf(ist,0)*
c     &                        wfovlp(Pbuff(0,0,1,1) ,Pbuff(0,0,1,1))
c            CommBP = CommBP + 4.0D0*CPvsqf(ist,0)*
c     &                        wfovlp(Qbuff(0,0,1,1) ,Pbuff(0,0,1,1))            
cc
c         enddo
cc
cc     operators with imu = 1
cc     ------------------------------------------------------------------
c      else if (imu.eq.1) then
cc        
cc     loop over all sp-states 
cc     ------------------------------------------------------------------
c         do ist=1,nst(1)+nst(2)
cc     
cc           calculate action of Q and P operator for mu=1
cc           ------------------------------------------------------------
c            do ir=0,nr
c               r = (ir+sr)*dr
c               do iz=0,nz
c                  z = (iz+sz-icm)*dz - CPzcm(0)
cc
c                  Qbuff(iz,ir,1,1)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
c                  Qbuff(iz,ir,2,1)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
c                  Qbuff(iz,ir,1,2)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
c                  Qbuff(iz,ir,2,2)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
cc     
cc                 ! upward laddering components of J_y
c                  Pbuff(iz,ir,1,1) =  0.5D0*pfac*(
c     &                                  -z*CPdrpsi(iz,ir,1,ist) 
c     &                                  +z*CPdppsi(iz,ir,1,ist) 
c     &                                  +r*CPdzpsi(iz,ir,1,ist)
c     &                                  -  CPpsi(iz,ir,2,ist,0))
c                  Pbuff(iz,ir,2,1) =  0.5D0*pfac*(
c     &                                  -z*CPdrpsi(iz,ir,2,ist) 
c     &                                  +z*CPdppsi(iz,ir,2,ist) 
c     &                                  +r*CPdzpsi(iz,ir,2,ist))
cc
cc                ! downward laddering components of J_y
c                 Pbuff(iz,ir,1,2)  =  0.5D0*pfac*(
c     &                                  -z*CPdrpsi(iz,ir,1,ist) 
c     &                                  -z*CPdppsi(iz,ir,1,ist) 
c     &                                  +r*CPdzpsi(iz,ir,1,ist))
c                 Pbuff(iz,ir,2,2)  =  0.5D0*pfac*(
c     &                                  -z*CPdrpsi(iz,ir,2,ist) 
c     &                                  -z*CPdppsi(iz,ir,2,ist) 
c     &                                  +r*CPdzpsi(iz,ir,2,ist)
c     &                                  +  CPpsi(iz,ir,1,ist,0))
c               enddo
c            enddo
cc
cc           evaluate 2qp-matrixelements
cc           ------------------------------------------------------------
c            if (ist.le.nst(1)) then
c               jstmin=1
c               jstmax=nst(1)
c            else
c               jstmin=nst(1)+1
c               jstmax=nst(1)+nst(2)
c            endif
cc
c            do jst=jstmin,jstmax
cc
cc              first the upward transitions
c               if (mpsi(jst)-mpsi(ist).eq.+1)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1)) 
c                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1),
c     &                                                          povlp)
cc
c                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
c                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
c     &                                                          qovlp)
c                     BVar   = BVar   + (qeta*qovlp)**2
c                     PVar   = PVar   + (peta*povlp)**2
c                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                  endif
c               endif
cc
cc              now the downward transitions
c               if (mpsi(jst)-mpsi(ist).eq.-1)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
c                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
c     &                                                          povlp)
cc
c                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
c                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
c     &                                                          qovlp)
c                     BVar   = BVar   + (qeta*qovlp)**2
c                     PVar   = PVar   + (peta*povlp)**2
c                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                  endif
c               endif
cc     
c            enddo
cc
cc           and finally the downward transitions into the T-mirrored
cc           states
c            do jst=jstmin,jstmax
cc
c               if (mpsi(ist)+mpsi(jst)+1.eq.+1)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.+1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     povlp=wfovlptr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
c                     call wfprojtr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
c     &                                                          povlp)
cc
c                     qovlp=wfovlptr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
c                     call wfprojtr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
c     &                                                          qovlp)
c                     BVar   = BVar   + (qeta*qovlp)**2
c                     PVar   = PVar   + (peta*povlp)**2
c                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                  endif
c               endif
cc     
c            enddo
cc     
cc           accumulate residual components
cc           ------------------------------------------------------------
c            BVar   = BVar   + 2.0D0*CPvsqf(ist,0)*(
c     &                        wfovlp(Qbuff(0,0,1,1) ,Qbuff(0,0,1,1))
c     &                       +wfovlp(Qbuff(0,0,1,2) ,Qbuff(0,0,1,2)))
c            PVar   = PVar   + 2.0D0*CPvsqf(ist,0)*(
c     &                        wfovlp(Pbuff(0,0,1,1) ,Pbuff(0,0,1,1))
c     &                       +wfovlp(Pbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
c            CommBP = CommBP + 4.0D0*CPvsqf(ist,0)*(
c     &                        wfovlp(Qbuff(0,0,1,1) ,Pbuff(0,0,1,1))
c     &                       +wfovlp(Qbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
cc
c         enddo
cc
cc     operators with imu = 2
cc     ------------------------------------------------------------------
c      else if (imu.eq.2) then
cc        
cc     loop over all sp-states 
cc     ------------------------------------------------------------------
c         do ist=1,nst(1)+nst(2)
cc     
cc           calculate action of Q and P operator for mu=1
cc           ------------------------------------------------------------
c            do ir=0,nr
c               r = (ir+sr)*dr
c               do iz=0,nz
c                  z = (iz+sz-icm)*dz - CPzcm(0)
cc
c                  Qbuff(iz,ir,1,1)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
c                  Qbuff(iz,ir,2,1)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
c                  Qbuff(iz,ir,1,2)=QlmOper(iz,ir)*CPpsi(iz,ir,1,ist,0)
c                  Qbuff(iz,ir,2,2)=QlmOper(iz,ir)*CPpsi(iz,ir,2,ist,0)
cc     
c                  Pbuff(iz,ir,1,1)  =  0.0D0
c                  Pbuff(iz,ir,2,1)  =  0.0D0
c                  Pbuff(iz,ir,1,2)  =  0.0D0
c                  Pbuff(iz,ir,2,2)  =  0.0D0
c               enddo
c            enddo
cc
cc           evaluate 2qp-matrixelements
cc           ------------------------------------------------------------
c            if (ist.le.nst(1)) then
c               jstmin=1
c               jstmax=nst(1)
c            else
c               jstmin=nst(1)+1
c               jstmax=nst(1)+nst(2)
c            endif
cc
c            do jst=jstmin,jstmax
cc
cc              first the upward transitions
c               if (mpsi(jst)-mpsi(ist).eq.+2)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1)) 
c                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,1),
c     &                                                          povlp)
cc
c                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1))  
c                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,1),
c     &                                                          qovlp)
c                     BVar   = BVar   + (qeta*qovlp)**2
c                     PVar   = PVar   + (peta*povlp)**2
c                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                  endif
c               endif
cc
cc              now the downward transitions
c               if (mpsi(jst)-mpsi(ist).eq.-2)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.-1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     povlp=wfovlp(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
c                     call wfproj(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
c     &                                                          povlp)
cc
c                     qovlp=wfovlp(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
c                     call wfproj(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
c     &                                                          qovlp)
c                     BVar   = BVar   + (qeta*qovlp)**2
c                     PVar   = PVar   + (peta*povlp)**2
c                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                  endif
c               endif
cc     
c            enddo
cc
cc           and finally the downward transitions into the T-mirrored
cc           states
c            do jst=jstmin,jstmax
cc
c               if (mpsi(ist)+mpsi(jst)+1.eq.+2)  then
c                  if (npar(jst)*npar(ist)*trzpar*isym.ne.+1) then
cc
c                     vi   = sqrtb(CPvsqf(ist,0))  
c                     ui   = sqrtb(one-vi*vi)
c                     vj   = sqrtb(CPvsqf(jst,0))
c                     uj   = sqrtb(one-vj*vj)
c                     qeta = uj*vi + vj*ui
c                     peta = uj*vi - vj*ui
cc
c                     povlp=wfovlptr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2)) 
c                     call wfprojtr(CPpsi(0,0,1,jst,0),Pbuff(0,0,1,2),
c     &                                                          povlp)
cc
c                     qovlp=wfovlptr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2))  
c                     call wfprojtr(CPpsi(0,0,1,jst,0),Qbuff(0,0,1,2),
c     &                                                          qovlp)
c                     BVar   = BVar   + (qeta*qovlp)**2
c                     PVar   = PVar   + (peta*povlp)**2
c                     CommBP = CommBP + 2.0D0*qeta*peta*povlp*qovlp
cc
c                  endif
c               endif
cc     
c            enddo
cc     
cc           accumulate residual components
cc           ------------------------------------------------------------
c            BVar   = BVar   + 2.0D0*CPvsqf(ist,0)*(
c     &                        wfovlp(Qbuff(0,0,1,1) ,Qbuff(0,0,1,1))
c     &                       +wfovlp(Qbuff(0,0,1,2) ,Qbuff(0,0,1,2)))
c            PVar   = PVar   + 2.0D0*CPvsqf(ist,0)*(
c     &                        wfovlp(Pbuff(0,0,1,1) ,Pbuff(0,0,1,1))
c     &                       +wfovlp(Pbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
c            CommBP = CommBP + 4.0D0*CPvsqf(ist,0)*(
c     &                        wfovlp(Qbuff(0,0,1,1) ,Pbuff(0,0,1,1))
c     &                       +wfovlp(Qbuff(0,0,1,2) ,Pbuff(0,0,1,2)))
cc
c         enddo
cc
c      endif
cc
c      return
c      end
c


