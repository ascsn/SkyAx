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
      subroutine comm_HQ(H_2q,
     &                   Q_2q,Q_Res,drQ_Res,dzQ_Res,dpQ_Res,ddQ_Res,
     &                   HQ_2q,HQ_Res,
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
      double precision H_2q(0:MTD)
      double precision Q_2q(0:MTD)
      double precision Hresp_2q(0:MTD)
      double precision HQ_2q(0:MTD)

      double precision Q_res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision HQ_Res(0:MZD,0:MRD,MSD,MPSID,2) 
      double precision Hresp_Res(0:MZD,0:MRD,MSD,MPSID,2)

      double precision drQ_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dzQ_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dpQ_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ddQ_Res(0:MZD,0:MRD,MSD,MPSID,2)

      integer          ist,itr,imu,ir,iz
c
c     ------------------------------------------------------------------
c
c     2qp-matrixelements of  ([H,Q] + H_{Q,resp}) = 2 B (-iP)
c     ------------------------------------------------------------------
      if (isMFresp.eq.1) then
c     
c        calculate the Q-response of the mean-field
c       ----------------------------------------------------------------
         call upd_respdens_odd(Q_2q,Q_Res,drQ_Res,dzQ_Res,dpQ_Res,imu)
         call upd_respflds_odd(imu)
c         call wrt_response()
         call upd_resp_H_odd(Hresp_2q,Hresp_res,imu) 
c  
         do itr=1,TRlen
            HQ_2q(itr) =  H_2q(itr)*Q_2q(itr) - Hresp_2q(itr)
         enddo
c
      else
c
         do itr=1,TRlen
            HQ_2q(itr) =  H_2q(itr)*Q_2q(itr)
         enddo
c
      endif
c
c     if no residual action --> all done
c     ------------------------------------------------------------------
      if (isresact.ne.1) then
         return
      endif
c
c     action of [H,Q] + H_{Q,resp} in residual space
c     ------------------------------------------------------------------
c     the case mu=0
      if (imu.eq.0) then
c
c     contributions from H_0 
      call act_mfh_even(
     &            Q_res(0,0,1,1,1),  drQ_res(0,0,1,1,1),
     &            dzQ_res(0,0,1,1,1),dpQ_res(0,0,1,1,1),
     &            ddQ_res(0,0,1,1,1), HQ_res(0,0,1,1,1),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0)
c
c     project out contributions which lead back to qp-space and 
c     accumulate on HQ_2q
c     ---> should be zero <---- !!!????
c      call project_qp(HQ_2q,HQ_res,-1,imu,1) !!!!!PK!!!!!
      call project_qp(HQ_2q,HQ_res,-1,imu,2)
c
c     add contributions from H_{Q,resp}
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               HQ_res(iz,ir,1,ist,1) = HQ_res(iz,ir,1,ist,1) 
     &                                -Hresp_res(iz,ir,1,ist,1)
               HQ_res(iz,ir,2,ist,1) = HQ_res(iz,ir,2,ist,1) 
     &                                -Hresp_res(iz,ir,2,ist,1)
               HQ_res(iz,ir,1,ist,2) = 0.0D0
               HQ_res(iz,ir,2,ist,2) = 0.0D0
            enddo
         enddo
      enddo
c
c     the case mu!=0
      else
c
c     contributions from H_0 
      call act_mfh_even(
     &            Q_res(0,0,1,1,1),  drQ_res(0,0,1,1,1),
     &            dzQ_res(0,0,1,1,1),dpQ_res(0,0,1,1,1),
     &            ddQ_res(0,0,1,1,1), HQ_res(0,0,1,1,1),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0)
      call act_mfh_even(
     &            Q_res(0,0,1,1,2),  drQ_res(0,0,1,1,2),
     &            dzQ_res(0,0,1,1,2),dpQ_res(0,0,1,1,2),
     &            ddQ_res(0,0,1,1,2), HQ_res(0,0,1,1,2),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0)
c
c     project out contributions which lead back to qp-space and 
c     accumulate on HQ_2q
c      call project_qp(HQ_2q,HQ_res,-1,imu,1) !!!!!PK!!!!!
      call project_qp(HQ_2q,HQ_res,-1,imu,2)
c
c     add contributions from H_{Q,resp}
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               HQ_res(iz,ir,1,ist,1) = HQ_res(iz,ir,1,ist,1) 
     &                                -Hresp_res(iz,ir,1,ist,1)
               HQ_res(iz,ir,2,ist,1) = HQ_res(iz,ir,2,ist,1) 
     &                                -Hresp_res(iz,ir,2,ist,1)
               HQ_res(iz,ir,1,ist,2) = HQ_res(iz,ir,1,ist,2) 
     &                                -Hresp_res(iz,ir,1,ist,2)
               HQ_res(iz,ir,2,ist,2) = HQ_res(iz,ir,2,ist,2) 
     &                                -Hresp_res(iz,ir,2,ist,2)
            enddo
         enddo
      enddo
c
      endif
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
      subroutine comm_HP(H_2q,
     &                   P_2q,P_Res,drP_Res,dzP_Res,dpP_Res,ddP_Res,
     &                   HP_2q,HP_Res,
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
      double precision H_2q(0:MTD)
      double precision P_2q(0:MTD)
      double precision Hresp_2q(0:MTD)
      double precision HP_2q(0:MTD)

      double precision P_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision HP_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision Hresp_Res(0:MZD,0:MRD,MSD,MPSID,2)

      double precision drP_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dzP_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision dpP_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision ddP_Res(0:MZD,0:MRD,MSD,MPSID,2)

      integer          ist,itr,imu
c
c     ------------------------------------------------------------------
c
c     2qp-matrixelements of ([H,(-iP)] + (-i)H_{P,resp}) = 2 C Q
c     ------------------------------------------------------------------
      if (isMFresp.eq.1) then
c     
c        calculate the P-response of the mean-field
c        ---------------------------------------------------------------
c        call upd_respdens_even_CHF()
c        call upd_respflds_even_CHF()
c        call wrt_response()        
         call upd_respdens_even(P_2q,P_Res,drP_Res,dzP_Res,dpP_Res,imu)
         call upd_respflds_even(imu)
c        call wrt_response()
c
         call upd_resp_H_even(Hresp_2q,Hresp_Res,imu)     
c
         do itr=1,TRlen
            HP_2q(itr) = H_2q(itr)*P_2q(itr) - Hresp_2q(itr)
         enddo
c     
      else
c
         do itr=1,TRlen
            HP_2q(itr) = H_2q(itr)*P_2q(itr)
         enddo
c
      endif
c
c     if no residual action --> all done
c     ------------------------------------------------------------------
      if (isresact.ne.1) then
         return
      endif        
c
c     action of -i([H,P] + H_{P,resp}) in residual space
c     ------------------------------------------------------------------
c     the case mu=0
      if (imu.eq.0) then
c
c     contributions from H_0
c     ( the second appearence of CPdzeffmas and CPsocpotz are just 
c       dummy fields because CPsocpotp and CPdpeffmas vanish in case
c       of axial symmetry)
      call act_mfh_even(
     &            P_Res(0,0,1,1,1),  drP_res(0,0,1,1,1),
     &            dzP_res(0,0,1,1,1),dpP_res(0,0,1,1,1),
     &            ddP_res(0,0,1,1,1), HP_res(0,0,1,1,1),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0)
c
c     project out contributions which lead back to qp-space 
c     ( should be zero ) 
c      call project_qp(HP_2q,HP_res,+1,imu,1) !!!!!PK!!!!!
      call project_qp(HP_2q,HP_res,+1,imu,2)
c
c     add contributions from H_{P,resp}
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               HP_res(iz,ir,1,ist,1) =-HP_res(iz,ir,1,ist,1) 
     &                                -Hresp_res(iz,ir,1,ist,1)
               HP_res(iz,ir,2,ist,1) =-HP_res(iz,ir,2,ist,1) 
     &                                -Hresp_res(iz,ir,2,ist,1)
               HP_res(iz,ir,1,ist,2) = 0.0D0
               HP_res(iz,ir,2,ist,2) = 0.0D0
            enddo
         enddo
      enddo
c
c     the case mu!=0
      else
c
c     contributions from H_0
c     ( the second appearence of CPdzeffmas and CPsocpotz are just 
c       dummy fields because CPsocpotp and CPdpeffmas vanish in case
c       of axial symmetry)
      call act_mfh_even(
     &            P_res(0,0,1,1,1),  drP_res(0,0,1,1,1),
     &            dzP_res(0,0,1,1,1),dpP_res(0,0,1,1,1),
     &            ddP_res(0,0,1,1,1), HP_res(0,0,1,1,1),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0)

      call act_mfh_even(
     &            P_res(0,0,1,1,2),  drP_res(0,0,1,1,2),
     &            dzP_res(0,0,1,1,2),dpP_res(0,0,1,1,2),
     &            ddP_res(0,0,1,1,2), HP_res(0,0,1,1,2),
     &            CPlocpot(0,0,1,0),                       1,
     &            CPeffmas(0,0,1,0),
     &            CPdreffmas,CPdzeffmas,CPdzeffmas,        1,
     &            CPsocpotr(0,0,1,0),CPsocpotz(0,0,1,0),
     &            CPsocpotz(0,0,1,0),                      1,
     &            CPconstr(0,0,0),CPlambda(0),             1,
     &            CPqpenrg(1,0),                           -1,
     &            0,1.0D0)
c
c     project out contributions which lead back to qp-space 
c     ( should be zero ) 
c      call project_qp(HP_2q,HP_res,+1,imu,1) !!!!!PK!!!!!
      call project_qp(HP_2q,HP_res,+1,imu,2)
c
c     add contributions from H_{P,resp}
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               HP_res(iz,ir,1,ist,1) =-HP_res(iz,ir,1,ist,1) 
     &                                -Hresp_res(iz,ir,1,ist,1)
               HP_res(iz,ir,2,ist,1) =-HP_res(iz,ir,2,ist,1) 
     &                                -Hresp_res(iz,ir,2,ist,1)
               HP_res(iz,ir,1,ist,2) =-HP_res(iz,ir,1,ist,2) 
     &                                -Hresp_res(iz,ir,1,ist,2)
               HP_res(iz,ir,2,ist,2) =-HP_res(iz,ir,2,ist,2) 
     &                                -Hresp_res(iz,ir,2,ist,2)
            enddo
         enddo
      enddo
c
      endif
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
      subroutine comm_NQ(Q_2q,Q_Res,NQ_2q,NQ_Res)
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
      double precision Q_2q(0:MTD)
      double precision NQ_2q(0:MTD)
      double precision Q_Res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision NQ_Res(0:MZD,0:MRD,MSD,MPSID,2)

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
     &            (1.0D0 - CPvsqf(ist,0) - CPvsqf(jst,0))*Q_2q(itr)
c
      enddo
c
c     residual action == zero
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               NQ_Res(iz,ir,1,ist,1) = 0.0D0
               NQ_Res(iz,ir,2,ist,1) = 0.0D0
               NQ_Res(iz,ir,1,ist,2) = 0.0D0
               NQ_Res(iz,ir,2,ist,2) = 0.0D0
            enddo
         enddo
      enddo
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
      function comm_QP(Q_2q,Q_res,P_2q,P_res,imu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'linresp.inc'  ! CPvsqf,ntr,MTR
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision Q_2q(0:MTD)
      double precision P_2q(0:MTD)
      double precision Q_res(0:MZD,0:MRD,MSD,MPSID,2)
      double precision P_res(0:MZD,0:MRD,MSD,MPSID,2)
      integer          imu

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
c
c         multi   = TRidx(itr,6)
c
         comm_QP_qp = comm_QP_qp + TRidx(itr,6)*Q_2q(itr)*P_2q(itr)
c
      enddo
c
c     contributions from residual action
c     ( factor 4 results from occupation (2) and commutator (2) )
c     ------------------------------------------------------------------
      if (isresact.ne.0) then
c
      if (imu.eq.0) then
c
c        only upper component of X_res = X_res(0,0,1,ist,1) is used
         do ist=1,nst(1)+nst(2)
            comm_QP_res = comm_QP_res +
     &           4.0D0*CPvsqf(ist,0)*
     &           wfovlp(Q_res(0,0,1,ist,1),P_res(0,0,1,ist,1))
         enddo
c
      else
c
c        upper & lower component of X_res contribute
         do ist=1,nst(1)+nst(2)
            comm_QP_res = comm_QP_res +
     &           4.0D0*CPvsqf(ist,0)*(
     &           wfovlp(Q_res(0,0,1,ist,1),P_res(0,0,1,ist,1))+
     &           wfovlp(Q_res(0,0,1,ist,2),P_res(0,0,1,ist,2)) )
         enddo
c
      endif
c     
      endif
c
      comm_QP = comm_QP_res + comm_QP_qp
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
      function var_oper(X_2q,X_res,imu)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'linresp.inc'  ! CPvsqf,TRlen,MTR,TRidx
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision X_2q(0:MTD)
      double precision X_res(0:MZD,0:MRD,MSD,MPSID,2)
      integer          imu

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
c
c         multi   = TRidx(itr,6)
c
         var_X_qp = var_X_qp + TRidx(itr,6)*X_2q(itr)*X_2q(itr)
c
      enddo
c
      var_X_qp=0.5D0*var_X_qp
c
c     contributions from residual action
c     ( factor 2 results from occupation (2) )
c     ------------------------------------------------------------------
      if (isresact.ne.0) then
c
      if (imu.eq.0) then
c
c        only upper component of X_res = X_res(0,0,1,ist,1) is used
         do ist=1,nst(1)+nst(2)
            var_X_res = var_X_res + 
     &           2.0D0*CPvsqf(ist,0)*
     &           wfovlp(X_res(0,0,1,ist,1),X_res(0,0,1,ist,1))
         enddo
c
      else
c
c        upper & lower component of X_res contribute
         do ist=1,nst(1)+nst(2)
            var_X_res = var_X_res +
     &           2.0D0*CPvsqf(ist,0)*(
     &           wfovlp(X_res(0,0,1,ist,1),X_res(0,0,1,ist,1))+
     &           wfovlp(X_res(0,0,1,ist,2),X_res(0,0,1,ist,2)) )
         enddo
c
      endif
c
      endif
c
      var_oper = var_X_res + var_X_qp
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     zero_oper(X_2q,X_res,ipart)      
c
c     set operator matrixelements (ipart=1) residual part (ipart=2) or
c     both (ipart=0) to zero
c     
c=======================================================================
      subroutine zero_oper(X_2q,X_res,ipart)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'   ! MZD,MRD,MSD,MPSID
      include 'linresp.inc'  ! CPvsqf,TRlen,MTR,TRidx
      include 'levels.inc'   ! nst
c
c     ------------------------------------------------------------------
c
      double precision X_2q(0:MTD)
      double precision X_res(0:MZD,0:MRD,MSD,MPSID,2)
      integer          ipart

      integer          ist,ir,iz,itr
c
c     ------------------------------------------------------------------
c
c     contributions from 2qp-matrixelements
c     ------------------------------------------------------------------
      if (ipart.ne.2) then
         do itr=1,TRlen
            X_2q(itr)=0.0D0
         end do
      endif

      if (ipart.ne.1) then
         do ist=1,nst(1)+nst(2)
            do ir=0,nr
               do iz=0,nz 
                  X_res(iz,ir,1,ist,1)=0.0D0
                  X_res(iz,ir,2,ist,1)=0.0D0
                  X_res(iz,ir,1,ist,2)=0.0D0
                  X_res(iz,ir,2,ist,2)=0.0D0
               enddo
            enddo
         enddo
      endif
c
c     ------------------------------------------------------------------
c  
      return
      end


c======================================================================
c
c     project_qp(X_2q,X_res,Xtpar,imu,iaccu)
c
c     project out components in sp-space from an operator-wavefunction
c     (=Pi_Omega)
c     
c        X_2q   : 2qp-matrixelements of X_res
c        X_res  : projection of X_res on residual space (output)
c                 or unprojected operator-wavefunctions (input)
c        imu    : angular-momentum transfer of X-Operator (\mu)
c        iaccu  : switch if matrixlements should be accumulated in X_2q
c                 or     if X_2q should be zeroed (iaccu=0)
c
c=======================================================================
      subroutine project_qp(X_2q,X_res,tpar,imu,iaccu)
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
      double precision X_2q(0:MTD)
      double precision X_res(0:MZD,0:MRD,MSD,MPSID,2)
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
            X_2q(itr)=0.0D0
         enddo
      endif
c
      if (iaccu.eq.2) then
         fac=0.0D0
      else
         fac=1.0D0
      endif
c
c     the case imu == 0      ==>  only X_res(0,0,1,1,1) is used !!!
c     ==================================================================
      if (imu.eq.0) then
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
         if (ist.ne.jst) then
           ovlp=wfovlp(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,1))  
           call wfproj(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,1),ovlp)
           call wfproj(CPpsi(0,0,1,ist,0),X_res(0,0,1,jst,1),tpar*ovlp)
         else
           ovlp=wfovlp(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,1))  
           call wfproj(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,1),ovlp)
         endif
c         
c        transform matrixelement into qp-space
c        ---------------------------------------------------------------
         X_2q(itr) = X_2q(itr) + fac*TReta(itr,tpidx)*ovlp
c
      enddo
c
      else
c
c     the case imu != 0  ==>  both components of X_res are used !!!
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
            ovlp=wfovlp(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,1))  
            call wfproj(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,1),ovlp)
            call wfproj(CPpsi(0,0,1,ist,0),X_res(0,0,1,jst,2),tpar*ovlp)
         else
            ovlp=wfovlp(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,2))  
            call wfproj(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,2),ovlp)
            call wfproj(CPpsi(0,0,1,ist,0),X_res(0,0,1,jst,1),tpar*ovlp)
         endif
         else
            if (ist.ne.jst) then
               ovlp=wfovlptr(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,2)) 
               call wfprojtr(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,2),ovlp)
               call wfprojtr(CPpsi(0,0,1,ist,0),X_res(0,0,1,jst,2),
     &                                                       -tpar*ovlp)
            else
               ovlp=wfovlptr(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,2)) 
               call wfprojtr(CPpsi(0,0,1,jst,0),X_res(0,0,1,ist,2),ovlp)
            endif
         endif
c         
c        transform matrixelement into qp-space
c        ---------------------------------------------------------------
         X_2q(itr) = X_2q(itr) + fac*TReta(itr,tpidx)*ovlp
c
      enddo
c     
      endif
c
c     if residual action is neglected ==> zero residual action
c     ------------------------------------------------------------------
      if (isresact.eq.0) then
c
         do ist=1,nst(1)+nst(2)
            do ir=0,nr
               do iz=0,nz
                  X_res(iz,ir,1,ist,1) = 0.0D0
                  X_res(iz,ir,2,ist,1) = 0.0D0
                  X_res(iz,ir,1,ist,2) = 0.0D0
                  X_res(iz,ir,2,ist,2) = 0.0D0
               enddo
            enddo  
         enddo
c
      endif
c
c     ------------------------------------------------------------------
c  
      return
      end


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
      subroutine upd_respdens_even(P_2q,P_res,
     &                             drP_res,dzP_res,dpP_res,
     &                             imu)
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
      double precision P_2q(0:MTD)
      double precision P_res(0:MZD,0:MRD,MSD,MPSID)
      double precision drP_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dzP_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dpP_res(0:MZD,0:MRD,MSD,MPSID)
      double precision wgtfac,wgtfac0,wgtfac2,multi,symfac
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
c            repair(iz,ir,1) = zero
c            repair(iz,ir,2) = zero
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
     &                +CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                +CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                +CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,2,ist) 
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
c                  repair(iz,ir,iq)  = wgtfac2*(
c     &                 CPpsi(iz,ir,1,ist,0) * CPpsi(iz,ir,1,jst,0) +
c     &                 CPpsi(iz,ir,2,ist,0) * CPpsi(iz,ir,2,jst,0)
c     &                                      ) + repair(iz,ir,iq)
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
     &                +CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0) 
     &                -CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0) 
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0) 
     &                +CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,2,ist) 
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
c                  repair(iz,ir,iq)  = wgtfac2*(
c     &                 CPpsi(iz,ir,1,ist,0) * CPpsi(iz,ir,1,jst,0) +
c     &                 CPpsi(iz,ir,2,ist,0) * CPpsi(iz,ir,2,jst,0)
c     &                                      ) + repair(iz,ir,iq)
c
               enddo
            enddo 
c
         endif
c
      enddo
c                  
c     accumulate time-even response densities from residual-components
c     ------------------------------------------------------------------
      if (isresact.eq.1) then
c
      do ist=1,nst(1)+nst(2)
c     
         wgtfac = 4.0D0*CPvsqf(ist,0)  ! comm.(x2), time-reversal (x2)
         wgtfac0= 4.0D0*CPvsqf(ist,0)*symfac
c
         if (ist.le.nst(1)) then
            iq = 1
         else
            iq = 2
         endif
c
         do ir=0,nr
            do iz=0,nz
c                  
c           propability density                                 
c           ------------------------------------------------------------
               rerho(iz,ir,iq)  = wgtfac*(
     &             +CPpsi(iz,ir,1,ist,0) * P_res(iz,ir,1,ist) 
     &             +CPpsi(iz,ir,2,ist,0) * P_res(iz,ir,2,ist)
     &                                   ) + rerho(iz,ir,iq)

c                  
c                 kinetic density                                     
c                 ------------------------------------------------------
                  rekin(iz,ir,iq)  = wgtfac*(
     &                +CPdrpsi(iz,ir,1,ist) * drP_res(iz,ir,1,ist)
     &                +CPdrpsi(iz,ir,2,ist) * drP_res(iz,ir,2,ist)
     &                +CPdzpsi(iz,ir,1,ist) * dzP_res(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,2,ist) * dzP_res(iz,ir,2,ist)
     &                +CPdppsi(iz,ir,1,ist) * dpP_res(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,2,ist) * dpP_res(iz,ir,2,ist)
     &                                      ) + rekin(iz,ir,iq)
c                  
c                 spin-orbit current                                  
c                 ------------------------------------------------------
                  resocr(iz,ir,iq) = wgtfac*0.5D0*(
     &                +CPpsi(iz,ir,1,ist,0) * dpP_res(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,2,ist,0) * dpP_res(iz,ir,2,ist) 
     &                +CPpsi(iz,ir,1,ist,0) * dzP_res(iz,ir,2,ist)   
     &                -CPpsi(iz,ir,2,ist,0) * dzP_res(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,ist) *   P_res(iz,ir,1,ist) 
     &                -CPdppsi(iz,ir,2,ist) *   P_res(iz,ir,2,ist) 
     &                +CPdzpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist)   
     &                -CPdzpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist)
     &                                      ) + resocr(iz,ir,iq)
c                  
                  resocp(iz,ir,iq) = wgtfac0*0.5D0*(
     &                +CPdzpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist) 
     &                +CPdzpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist) 
     &                +CPdrpsi(iz,ir,1,ist) *   P_res(iz,ir,1,ist)   
     &                -CPdrpsi(iz,ir,2,ist) *   P_res(iz,ir,2,ist)
     &                -CPpsi(iz,ir,1,ist,0) * dzP_res(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,ist,0) * dzP_res(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,1,ist,0) * drP_res(iz,ir,1,ist)   
     &                +CPpsi(iz,ir,2,ist,0) * drP_res(iz,ir,2,ist)
     &                                      ) + resocp(iz,ir,iq)
c                  
                  resocz(iz,ir,iq) = wgtfac*0.5D0*(
     &                -CPpsi(iz,ir,1,ist,0) * dpP_res(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,ist,0) * dpP_res(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,1,ist,0) * drP_res(iz,ir,2,ist)   
     &                +CPpsi(iz,ir,2,ist,0) * drP_res(iz,ir,1,ist)
     &                -CPdppsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist) 
     &                -CPdppsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist) 
     &                -CPdrpsi(iz,ir,2,ist) *   P_res(iz,ir,1,ist)   
     &                +CPdrpsi(iz,ir,1,ist) *   P_res(iz,ir,2,ist)
     &                                      ) + resocz(iz,ir,iq)
c
               enddo
            enddo
      enddo
c
      endif
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
      subroutine upd_respdens_odd(Q_2q,Q_res,drQ_res,dzQ_res,dpQ_res,
     &                            imu)
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
      double precision Q_2q(0:MTD)
      double precision Q_res(0:MZD,0:MRD,MSD,MPSID)
      double precision drQ_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dzQ_res(0:MZD,0:MRD,MSD,MPSID)
      double precision dpQ_res(0:MZD,0:MRD,MSD,MPSID)
      double precision trecurr(0:MZD,0:MRD,MID)
      double precision trecurp(0:MZD,0:MRD,MID)
      double precision trecurz(0:MZD,0:MRD,MID)
      double precision trespnr(0:MZD,0:MRD,MID)
      double precision trespnp(0:MZD,0:MRD,MID)
      double precision trespnz(0:MZD,0:MRD,MID)
      double precision treskdr(0:MZD,0:MRD,MID)
      double precision treskdp(0:MZD,0:MRD,MID)
      double precision treskdz(0:MZD,0:MRD,MID)
      double precision wgtfac,wgtfac0,wgtfac2,symfac,multi
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
C            recurr(iz,ir,1) = zero
c            recurr(iz,ir,2) = zero
c            recurp(iz,ir,1) = zero
c            recurp(iz,ir,2) = zero
c            recurz(iz,ir,1) = zero
c            recurz(iz,ir,2) = zero
c            respnr(iz,ir,1) = zero
c            respnr(iz,ir,2) = zero
c            respnp(iz,ir,1) = zero
c            respnp(iz,ir,2) = zero
c            respnz(iz,ir,1) = zero
c            respnz(iz,ir,2) = zero
c            reskdr(iz,ir,1) = zero
c            reskdr(iz,ir,2) = zero
c            reskdp(iz,ir,1) = zero
c            reskdp(iz,ir,2) = zero
c            reskdz(iz,ir,1) = zero
c            reskdz(iz,ir,2) = zero
            trecurr(iz,ir,1) = zero
            trecurr(iz,ir,2) = zero
            trecurp(iz,ir,1) = zero
            trecurp(iz,ir,2) = zero
            trecurz(iz,ir,1) = zero
            trecurz(iz,ir,2) = zero
            trespnr(iz,ir,1) = zero
            trespnr(iz,ir,2) = zero
            trespnp(iz,ir,1) = zero
            trespnp(iz,ir,2) = zero
            trespnz(iz,ir,1) = zero
            trespnz(iz,ir,2) = zero
            treskdr(iz,ir,1) = zero
            treskdr(iz,ir,2) = zero
            treskdp(iz,ir,1) = zero
            treskdp(iz,ir,2) = zero
            treskdz(iz,ir,1) = zero
            treskdz(iz,ir,2) = zero
         enddo
      enddo
c                  
c     accumulate time-odd response densities from 2q-matrixelements
c     ------------------------------------------------------------------
cGH
c      write(*,*) 'TRlen=',TRlen
!$OMP PARALLEL DO private(ist,jst,jtp,iq,multi,wgtfac,wgtfac0,
!$OMP+ ir,iz,wgtfac2)
!$OMP+ reduction(+:trecurr,trecurp,trecurz,trespnr,trespnp,trespnz,
!$OMP+ treskdr,treskdp,treskdz)
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
                  trecurr(iz,ir,iq) = -0.5D0*wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trecurr(iz,ir,iq)
                  
                  trecurp(iz,ir,iq) = 0.5D0*wgtfac0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0)  
     &                +CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trecurp(iz,ir,iq)
                  
                  trecurz(iz,ir,iq) = -0.5D0*wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trecurz(iz,ir,iq)
c                  
c                 spin density
c                 ------------------------------------------------------
                  trespnr(iz,ir,iq) = wgtfac0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                +CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + trespnr(iz,ir,iq)
                  
                  trespnp(iz,ir,iq) = -wgtfac*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + trespnp(iz,ir,iq)
                  
                  trespnz(iz,ir,iq) = wgtfac0*(
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0)  
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trespnz(iz,ir,iq) 
c                  
c                 kinetic spin density
c                 ------------------------------------------------------
                  treskdr(iz,ir,iq) = wgtfac0*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + treskdr(iz,ir,iq)
                  
                  treskdp(iz,ir,iq) = -wgtfac*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)  
     &                -CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + treskdp(iz,ir,iq)
                  
                  treskdz(iz,ir,iq) = wgtfac0*(
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist)  
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist)
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist)  
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist)
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist)  
     &                -CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)
     &                                      ) + treskdz(iz,ir,iq) 
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
                  trecurr(iz,ir,iq) = -0.5D0*wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdrpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,1,jst,0) * CPdrpsi(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdrpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trecurr(iz,ir,iq)
                  
                  trecurp(iz,ir,iq) = 0.5D0*wgtfac0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdppsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,1,jst,0) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdppsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trecurp(iz,ir,iq)
                  
                  trecurz(iz,ir,iq) = -0.5D0*wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPdzpsi(iz,ir,1,ist)  
     &                +CPpsi(iz,ir,1,jst,0) * CPdzpsi(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,2,jst) * CPpsi(iz,ir,1,ist,0)  
     &                -CPdzpsi(iz,ir,1,jst) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trecurz(iz,ir,iq)
c                  
c                 spin density
c                 ------------------------------------------------------
                  trespnr(iz,ir,iq) = wgtfac0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                +CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + trespnr(iz,ir,iq)
                  
                  trespnp(iz,ir,iq) = -wgtfac*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,2,ist,0)  
     &                -CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,1,ist,0)
     &                                      ) + trespnp(iz,ir,iq)
                  
                  trespnz(iz,ir,iq) = wgtfac0*(
     &                -CPpsi(iz,ir,2,jst,0) * CPpsi(iz,ir,1,ist,0)  
     &                -CPpsi(iz,ir,1,jst,0) * CPpsi(iz,ir,2,ist,0)
     &                                      ) + trespnz(iz,ir,iq) 
c                  
c                 kinetic spin density
c                 ------------------------------------------------------
                  treskdr(iz,ir,iq) = wgtfac0*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist)
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)  
     &                -CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + treskdr(iz,ir,iq)
                  
                  treskdp(iz,ir,iq) = -wgtfac*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,1,ist)
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,1,ist)
     &                                      ) + treskdp(iz,ir,iq)
                  
                  treskdz(iz,ir,iq) = wgtfac0*(
     &                -CPdrpsi(iz,ir,2,jst) * CPdrpsi(iz,ir,1,ist)  
     &                -CPdrpsi(iz,ir,1,jst) * CPdrpsi(iz,ir,2,ist)
     &                -CPdzpsi(iz,ir,2,jst) * CPdzpsi(iz,ir,1,ist)  
     &                -CPdzpsi(iz,ir,1,jst) * CPdzpsi(iz,ir,2,ist)
     &                +CPdppsi(iz,ir,2,jst) * CPdppsi(iz,ir,1,ist)  
     &                +CPdppsi(iz,ir,1,jst) * CPdppsi(iz,ir,2,ist)
     &                                      ) + treskdz(iz,ir,iq) 
c                  
               enddo
            enddo 
c
         endif
c     
      enddo
!$OMP END PARALLEL DO


      do ir=0,nr
         do iz=0,nz
            recurr(iz,ir,1) = trecurr(iz,ir,1)
            recurr(iz,ir,2) = trecurr(iz,ir,2)
            recurp(iz,ir,1) = trecurp(iz,ir,1)
            recurp(iz,ir,2) = trecurp(iz,ir,2)
            recurz(iz,ir,1) = trecurz(iz,ir,1)
            recurz(iz,ir,2) = trecurz(iz,ir,2)
            respnr(iz,ir,1) = trespnr(iz,ir,1)
            respnr(iz,ir,2) = trespnr(iz,ir,2)
            respnp(iz,ir,1) = trespnp(iz,ir,1)
            respnp(iz,ir,2) = trespnp(iz,ir,2)
            respnz(iz,ir,1) = trespnz(iz,ir,1)
            respnz(iz,ir,2) = trespnz(iz,ir,2)
            reskdr(iz,ir,1) = treskdr(iz,ir,1)
            reskdr(iz,ir,2) = treskdr(iz,ir,2)
            reskdp(iz,ir,1) = treskdp(iz,ir,1)
            reskdp(iz,ir,2) = treskdp(iz,ir,2)
            reskdz(iz,ir,1) = treskdz(iz,ir,1)
            reskdz(iz,ir,2) = treskdz(iz,ir,2)
         enddo
      enddo






c                  
c     accumulate time-odd response densities from residual-components
c     ------------------------------------------------------------------
      if (isresact.eq.1) then
C GH
C      write(*,*) 'nst(1)=',nst(1),' nst(2)=',nst(2)
      do ist=1,nst(1)+nst(2)
c     
         wgtfac  = 4.0D0*CPvsqf(ist,0)
         wgtfac0 = 4.0D0*CPvsqf(ist,0)*symfac
c
         if (ist.le.nst(1)) then
            iq = 1
         else
            iq = 2
         endif
c
         do ir=0,nr
            do iz=0,nz
c 
c                 current
c                 ------------------------------------------------------
                  recurr(iz,ir,iq) = -0.5D0*wgtfac*(
     &                +CPpsi(iz,ir,1,ist,0) * drQ_res(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,2,ist,0) * drQ_res(iz,ir,2,ist) 
     &                -CPdrpsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist) 
     &                -CPdrpsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist)
     &                                      ) + recurr(iz,ir,iq)
                  
                  recurp(iz,ir,iq) = 0.5D0*wgtfac0*(
     &                +CPpsi(iz,ir,1,ist,0) * dpQ_res(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,2,ist,0) * dpQ_res(iz,ir,2,ist) 
     &                +CPdppsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist) 
     &                +CPdppsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist)
     &                                      ) + recurp(iz,ir,iq)
                  
                  recurz(iz,ir,iq) = -0.5D0*wgtfac*(
     &                +CPpsi(iz,ir,1,ist,0) * dzQ_res(iz,ir,1,ist) 
     &                +CPpsi(iz,ir,2,ist,0) * dzQ_res(iz,ir,2,ist) 
     &                -CPdzpsi(iz,ir,1,ist) *   Q_res(iz,ir,1,ist) 
     &                -CPdzpsi(iz,ir,2,ist) *   Q_res(iz,ir,2,ist)
     &                                      ) + recurz(iz,ir,iq)
c                  
c                 spin density
c                 ------------------------------------------------------
                  respnr(iz,ir,iq) = wgtfac0*(
     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,2,ist) 
     &                +CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,1,ist)
     &                                      ) + respnr(iz,ir,iq)
                  
                  respnp(iz,ir,iq) = -wgtfac*(
     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,2,ist) 
     &                -CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,1,ist)
     &                                      ) + respnp(iz,ir,iq)
                  
                  respnz(iz,ir,iq) = wgtfac0*(
     &                +CPpsi(iz,ir,1,ist,0) *   Q_res(iz,ir,1,ist) 
     &                -CPpsi(iz,ir,2,ist,0) *   Q_res(iz,ir,2,ist)
     &                                      ) + respnz(iz,ir,iq) 
c                  
c                 kinetic spin density
c                 ------------------------------------------------------
                  reskdr(iz,ir,iq) = wgtfac0*(
     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,2,ist)  
     &                +CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,2,ist)  
     &                +CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,2,ist)  
     &                +CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,1,ist)
     &                                      ) + reskdr(iz,ir,iq)
                  
                  reskdp(iz,ir,iq) = -wgtfac*(
     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,2,ist)  
     &                -CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,1,ist)
     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,2,ist)  
     &                -CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,1,ist)
     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,2,ist)  
     &                -CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,1,ist)
     &                                      ) + reskdp(iz,ir,iq)
                  
                  reskdz(iz,ir,iq) = wgtfac0*(
     &                +CPdrpsi(iz,ir,1,ist) * drQ_res(iz,ir,1,ist)  
     &                -CPdrpsi(iz,ir,2,ist) * drQ_res(iz,ir,2,ist)
     &                +CPdzpsi(iz,ir,1,ist) * dzQ_res(iz,ir,1,ist)  
     &                -CPdzpsi(iz,ir,2,ist) * dzQ_res(iz,ir,2,ist)
     &                +CPdppsi(iz,ir,1,ist) * dpQ_res(iz,ir,1,ist)  
     &                -CPdppsi(iz,ir,2,ist) * dpQ_res(iz,ir,2,ist)
     &                                      ) + reskdz(iz,ir,iq) 
            enddo
         enddo
c
      enddo
      endif
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
c        dynamic potential A
c        ---------------------------------------------------------------
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
              buffr1 = redynpotr(iz,ir,1)
              buffr2 = redynpotr(iz,ir,2)
              buffp1 = redynpotp(iz,ir,1)
              buffp2 = redynpotp(iz,ir,2)
              buffz1 = redynpotz(iz,ir,1)
              buffz2 = redynpotz(iz,ir,2)
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
c      write(21,'(a)')         
c     &  'UPD_LOCPOTWGT() : response weight dU/drho updated'
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
      subroutine upd_resp_H_even(Hresp_2q,Hresp_res,imu)
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
      double precision Hresp_2q(0:MTD)
      double precision Hresp_res(0:MZD,0:MRD,MSD,MPSID,2)
      integer          itr,ist,jst,ir,iz,jtp,iq,imu,zpidx
      double precision ovlp
c
c     ------------------------------------------------------------------
c     the case imu==0
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
c      call wrt_2D_spin(LRwrk2_spinor,'dzeffmas__',1,'tst')
c      call wrt_2D_spin(LRwrk1_spinor,'dreffmas__',0,'tst')
c
c        calculate the action of response Hamiltonian on sp-states
c       ( the second appearence of reeffmas and resocpotz are just 
c         dummy fields because reeffmas and resocpotz vanish in case
c         of axial symmetry)
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
     &            0,1.0D0) 
c
         call project_qp(Hresp_2q,Hresp_res,+1,imu,0)
c
c     the case imu!=0
      else
c
c
c        prepare derivatives of mean-fields
c
         zpidx=1+mod(imu,2)
c
      call nablar_scalar(reeffmas(0,0,1),LRwrk1_spinor(0,0,1),
     &                   imu,1,one) 
      call nablar_scalar(reeffmas(0,0,2),LRwrk1_spinor(0,0,2),
     &                   imu,1,one) 
      call nablaz_scalar(reeffmas(0,0,1),LRwrk2_spinor(0,0,1),
     &                   zpidx,1,one)
      call nablaz_scalar(reeffmas(0,0,2),LRwrk2_spinor(0,0,2),
     &                   zpidx,1,one)
      call nablap_scalar(reeffmas(0,0,1),LRwrk3_spinor(0,0,1),
     &                   imu,+1,1,one) 
      call nablap_scalar(reeffmas(0,0,2),LRwrk3_spinor(0,0,2),
     &                   imu,+1,1,one) 
c
c      call wrt_2D_spin(LRwrk3_spinor,'dpeffmas__',1,'tst')
c      call wrt_2D_spin(LRwrk2_spinor,'dzeffmas__',0,'tst')
c      call wrt_2D_spin(LRwrk1_spinor,'dreffmas__',0,'tst')
c
c        calculate the action of response Hamiltonian on sp-states
c
         call act_mfh_even(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,1),
     &            relocpot(0,0,1),                          1,
     &            reeffmas,LRwrk1_spinor,LRwrk2_spinor,
     &            LRwrk3_spinor,                            1,
     &            resocpotr,resocpotz,resocpotp,            1,
     &            CPconstr(0,0,0),CPlambda(0),              0,
     &            CPspenrg(1,0),                            0,
     &            imu,0.5D0)  
c
         call act_mfh_even(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,2),
     &            relocpot(0,0,1),                          1,
     &            reeffmas,LRwrk1_spinor,LRwrk2_spinor,
     &            LRwrk3_spinor,                            1,
     &            resocpotr,resocpotz,resocpotp,            1,
     &            CPconstr(0,0,0),CPlambda(0),              0,
     &            CPspenrg(1,0),                            0,
     &            -imu,0.5D0)  
c     &            -imu,(+1.0D0)**imu)  
c
         call project_qp(Hresp_2q,Hresp_res,+1,imu,0)
c
      endif
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
      subroutine upd_resp_H_odd(Hresp_2q,Hresp_res,imu)
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
      double precision Hresp_2q(0:MTD)
      double precision Hresp_res(0:MZD,0:MRD,MSD,MPSID,2)
      integer          itr,ist,jst,ir,iz,jtp,iq,imu
      integer          izpidx
c
c     ------------------------------------------------------------------
c     the case imu==0
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
c      call wrt_2D_spin(LRwrk1_spinor(0,0,1),'divdyn____',1,'tst')
c      call wrt_2D_spin(redynpotr(0,0,1),'dynpotr___',0,'tst')
c      call wrt_2D_spin(redynpotz(0,0,1),'dynpotz___',0,'tst')
c      call wrt_2D_spin(redynpotp(0,0,1),'dynpotp___',0,'tst')
c       call wrt_2D_spin(respnpotr(0,0,1),'spnpotr___',0,'tst')
c       call wrt_2D_spin(respnpotz(0,0,1),'spnpotz___',0,'tst')
c       call wrt_2D_spin(respnpotp(0,0,1),'spnpotp___',0,'tst')
c
c        calculate the action of response Hamiltonian on sp-states
c
         call act_mfh_odd(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,1),
     &            redynpotr,redynpotz,redynpotp,LRwrk1_spinor, 1,
     &            respnpotr,respnpotz,respnpotp,               1,
     &            respnmasr,respnmasz,respnmasp,               0,
     &            0,1.0D0)  
c
         call project_qp(Hresp_2q,Hresp_res,-1,imu,0)
c
c     the case imu!=0
      else
c
c        set parities of dynamic potential A
         izpidx =  1 + mod(abs(imu),2)
c
c        prepare derivatives of mean-fields
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
c      call wrt_2D_spin(LRwrk1_spinor(0,0,1),'divdyn____',1,'tst')
c      call wrt_2D_spin(redynpotr(0,0,1),'dynpotr___',0,'tst')
c      call wrt_2D_spin(redynpotz(0,0,1),'dynpotz___',0,'tst')
c      call wrt_2D_spin(redynpotp(0,0,1),'dynpotp___',0,'tst')
c       call wrt_2D_spin(respnpotr(0,0,1),'spnpotr___',0,'tst')
c       call wrt_2D_spin(respnpotz(0,0,1),'spnpotz___',0,'tst')
c       call wrt_2D_spin(respnpotp(0,0,1),'spnpotp___',0,'tst')
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
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,1),
     &            redynpotr,redynpotz,redynpotp,LRwrk1_spinor, 1,
     &            respnpotr,respnpotz,respnpotp,               1,
     &            respnmasr,respnmasz,respnmasp,               0,
     &            imu,0.5D0)  
c
         call act_mfh_odd(
     &            CPpsi(0,0,1,1,0),CPdrPsi(0,0,1,1),
     &            CPdzpsi(0,0,1,1),CPdpPsi(0,0,1,1),
     &            CPddpsi(0,0,1,1), Hresp_res(0,0,1,1,2),
     &            redynpotr,redynpotz,redynpotp,LRwrk1_spinor, 1,
     &            respnpotr,respnpotz,respnpotp,               1,
     &            respnmasr,respnmasz,respnmasp,               0,
     &            -imu,0.5D0)  
c
         call project_qp(Hresp_2q,Hresp_res,-1,imu,0)
c
      endif
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
      call wrt_2D_spin(rerho(0,0,1) ,'rerho_____',1,'ts3')
      call wrt_2D_spin(rekin(0,0,1) ,'rekin_____',0,'ts3')
      call wrt_2D_spin(resocz(0,0,1),'resocz____',0,'ts3')
      call wrt_2D_spin(resocr(0,0,1),'resocr____',0,'ts3')
      call wrt_2D_spin(resocp(0,0,1),'resocp____',0,'ts3')
      call wrt_2D_spin(recurr(0,0,1),'recurr____',0,'ts3')
      call wrt_2D_spin(recurp(0,0,1),'recurp____',0,'ts3')
      call wrt_2D_spin(recurz(0,0,1),'recurz____',0,'ts3')
      call wrt_2D_spin(respnr(0,0,1),'respnr____',0,'ts3')
      call wrt_2D_spin(respnp(0,0,1),'respnp____',0,'ts3')
      call wrt_2D_spin(respnz(0,0,1),'respnz____',0,'ts3')
      call wrt_2D_spin(reskdr(0,0,1),'reskdr____',0,'ts3')
      call wrt_2D_spin(reskdp(0,0,1),'reskdp____',0,'ts3')
      call wrt_2D_spin(reskdz(0,0,1),'reskdz____',0,'ts3')
      call wrt_2D_spin(relocpot(0,0,1) ,'relocpot__',0,'ts3')
      call wrt_2D_spin(reeffmas(0,0,1) ,'reeffmas__',0,'ts3')
      call wrt_2D_spin(redynpotr(0,0,1),'redynpotr_',0,'ts3')
      call wrt_2D_spin(redynpotz(0,0,1),'redynpotz_',0,'ts3')
      call wrt_2D_spin(redynpotp(0,0,1),'redynpotp_',0,'ts3')
      call wrt_2D_spin(respnpotr(0,0,1),'respnpotr_',0,'ts3')
      call wrt_2D_spin(respnpotz(0,0,1),'respnpotz_',0,'ts3')
      call wrt_2D_spin(respnpotp(0,0,1),'respnpotp_',0,'ts3')
      call wrt_2D_spin(resocpotr(0,0,1),'resocpotr_',0,'ts3')
      call wrt_2D_spin(resocpotz(0,0,1),'resocpotz_',0,'ts3')
      call wrt_2D_spin(resocpotp(0,0,1),'resocpotp_',0,'ts3')
      call wrt_2D_spin(respnmasr(0,0,1),'respnmasr_',0,'ts3')
      call wrt_2D_spin(respnmasz(0,0,1),'respnmasz_',0,'ts3')
      call wrt_2D_spin(respnmasp(0,0,1),'respnmasp_',0,'ts3')
c
      return
      end

