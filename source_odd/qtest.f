c-----hcPT--------------------------------------------------------------
c
      subroutine hcPT
c
c     Test-routine for Commutator [H,P] 
c     using P-operator 'presp' from routine 'crankin'
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
      write(7,'(a,1pg12.4)') 
     &  ' now test HCOMMP: initial <P**2>=',opovlp(presp,presp)
c
      call hcommP(presp,qresp)
c
c     test
c
      write(7,'(a,3(1pg12.4))')
     &     ' hcommP: <P**2>,<[P,[H,P]]>,<Q**2>=',
     &     opovlp(presp,presp),opovlp(presp,qresp),opovlp(qresp,qresp)
c
      call hcommQ(qresp,presp,ptest,0)
      write(7,'(a,3(1pg12.4))')
     &     ' hcommQ: <[Q,[H,Q]]>,<P_new**2>=',
     &     opovlp(ptest,qresp),opovlp(ptest,ptest)
c      
      return
      end 
c
c-----hcQT--------------------------------------------------------------
c
      subroutine hcQT(icase)
c
c     Test-routine for Commutator [H,Q] of Hamiltonian 
c     with cases {b1=0,b1 not 0}x{Q1=2z^2-r^2, Q2=2z^2, Q3=-r^2}
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
c     0. step:
c     
c
      if(itest.eq.0) then       ! only for testing
c
c     test particle number
c
         do ir=0,nr
            do iz=0,nz
c               wrk1(iz,ir)  = 0.0
c               wrk2(iz,ir)  = 0.0
               wrk4(iz,ir) = rho(iz,ir,1) + rho(iz,ir,2)
            enddo
         enddo
         write(17,'(/a,1(1pg12.4))')
     &        '> hcQT: test output : Tr{1}  = #particles =',
     &        f1int(wrk4)
c     
         do i=1,2               ! i=1 -> b1=0, i=2 -> b1 not 0
c
c     setting force parametres
c
            if (i.eq.1)  then
               b1  = 0.0
               b1p = 0.0
            else
               b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
               b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
            endif
            b1c = b1-b1p
c
            write(17,'(a,2(1pg12.4))')
     &           '> hcQT: test case: b1, b1p  = ',b1,b1p
c
c
c     calculating test masses icases: 3 , 1 , 2
c 
c     mass1 \       / Q1 = 2 z^2 - r^2 => (Nabla Q1)^2 =  16 z^2 +  4 r^2
c     mass2  | <-> |  Q2 = 2 z^2       => (Nabla Q2)^2 =  16 z^2
c     mass3 /       \ Q3 = - r^2       => (Nabla Q3)^2 =   4 r^2
c
c     making integrant in workspaces
c
            do ir=0,nr
               do iz=0,nz
c     
c     coordinates
c
                  r  = dr*(ir+sr)
                  r2 = r*r
                  z  = dz*(iz+sz-icm)
                  z2 = z*z
c     
c    h2meffq ??? factor 2.0 : +2 b1 ???
c
                  h2meffp 
     &                 =(h2m(1)+b1c*rho(iz,ir,1)+b1*rho(iz,ir,2))
                  h2meffn 
     &                 =(h2m(2)+b1c*rho(iz,ir,2)+b1*rho(iz,ir,1)) 
c
                  wh2meff = h2meffp+h2meffn
c
                  wrk1(iz,ir) = wh2meff * ( z2 ) * rho(iz,ir,1)
                  wrk2(iz,ir) = wh2meff * ( r2 ) * rho(iz,ir,2)
               enddo
            enddo
            
c
c     final integrating, with factor
c
            dummy1 = f1int(wrk1)
            dummy2 = f1int(wrk2)
            xmass1 =  8.0  * (4.0*dummy1+dummy2)
            xmass2 = 32.0  * (dummy1)
            xmass3 =  8.0  * (dummy2)
c
CR250599
c
c            write(7,'(a/,a,3(1pg12.4))')
c     &           '> hcQT: 2z^2-r^2,  2z^2,   -r^2',
c     &           '> hcQT: xmass1,  xmass2, xmass3 (inv) = ',
c     &           xmass1,xmass2,xmass3     
c
            if(icase.eq.1) then
               write(17,'(a/,a,1(1pg12.4))')
     &              '> hcQT: 2z^2 ',
     &              '> hcQT: xmass2 (inv) = ', xmass2 
            else if(icase.eq.2) then
               write(17,'(a/,a,1(1pg12.4))')
     &              '> hcQT: -r^2 ',
     &              '> hcQT: xmass3 (inv) = ', xmass3 
            else               
               write(17,'(a/,a,1(1pg12.4))')
     &              '> hcQT: 2z^2-r^2 ',
     &              '> hcQT: xmass1 (inv) = ', xmass1
            endif
c
c     all done
c
         enddo
      endif
      return
      end 
c
c-----Qabtst-----------------------------------------------------------
c
      subroutine Qabtst(icase)
c
c     calculates the matrix-elements Q_\alpha\beta of Q=2z^2-r^2 for 
c     testing and stores results on 'qresp' where it can be processed
c     further.
c     Three cases are provided:
c      icase = 1  --> Q=2z^2
c      icase = 2  --> Q=-r^2
c      icase = 3  --> Q=2z^2-r^2
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
c
      if(icase.eq.1) then
        write(7,'(/a)') ' test HCOMMQ with Q=2*z**2:'
        write(17,'(/a)') '# Qabtst: test HCOMMQ with Q=2*z^2:'
      else if(icase.eq.2) then
        write(7,'(/a)') ' test HCOMMQ with Q=-r**2:'
        write(17,'(/a)') '# Qabtst: test HCOMMQ with Q=-r^2:'
      else
        write(7,'(/a)') ' test HCOMMQ with Q=2*z**2-r**2:'
        write(17,'(/a)') '# Qabtst: test HCOMMQ with Q=2*z^2-r^2:'
      endif
c
      sumq   = 0.0
      sumq2  = 0.0
      ipostr = 0
c
      do ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
c
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         vi    = sqrt(vsqf(ist))
         ui    = sqrt(one-vi*vi)
c
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
            jstmax = nst(1)+nst(2)
         endif
c
         do jst=ist,jstmax      ! <-> beta
           if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
             ipostr = 1 + ipostr
             vj    = sqrt(vsqf(jst))
             uj    = sqrt(one-vj*vj)
             if(icase.eq.1) then
               do ir=0,nr
                 do iz=0,nz
                   z  = dz*(iz+sz-icm)
                   z2 = z*z
                   wrk1(iz,ir) =                                 ! Q=2z^2
     &                     2.0*z2*(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                           + psi(iz,ir,2,ist)*psi(iz,ir,2,jst))
                 enddo
               enddo
             elseif(icase.eq.2) then
               do ir=0,nr
                 r  = dr*(ir+sr)
                 r2 = r*r
                 do iz=0,nz
                   wrk1(iz,ir) =                                 ! Q=2z^2
     &                     -r2*(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                         + psi(iz,ir,2,ist)*psi(iz,ir,2,jst))
                 enddo
               enddo
             else
               do ir=0,nr
                 r  = dr*(ir+sr)
                 r2 = r*r
                 do iz=0,nz
                   z  = dz*(iz+sz-icm)
                   wrk1(iz,ir) =  (2.0*(z*z)-r2)               ! quadrupole
     &                     *(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                       + psi(iz,ir,2,ist)*psi(iz,ir,2,jst))
                 enddo
               enddo
            endif
c
c     integrating and storing
c
            qmatre = f1int(wrk1)
            qresp(ipostr) = qmatre*(ui*vj+uj*vi)
CR---------------------------------------------------------------------
cr fehlt faktor 2?
            if(ist.eq.jst) then
              sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
              sumq  = qmatre*vi*vi + sumq
            else
              sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
            endif
CR---------------------------------------------------------------------
c
          endif
        enddo
      enddo
      iposmx = ipostr
      write(7,'(a,2f8.3,a,2(1pg12.4))')  
     &  ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
      write(17,'(a,2f8.3,a,2(1pg12.4))')  
     &  '# Qabtst: for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
c
c     testing by calling commutator
c
      call hcommQ(qresp,presp,ptest,icase)
c     
      return
      end 
c
c-----Sabtst-----------------------------------------------------------
CR040399
c
      subroutine Sabtst(icase)
c
c     calculates the matrix-elements Q_\alpha\beta of Q=sigma_r for 
c     testing and stores results on 'qresp' where it can be processed
c     further. (more or less the same as Qabtst)
c      icase = 3  --> Q=2z^2-r^2
c
c     icase is only used for the call of hcommQ
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
c
      write(7,'(/a)') ' test HCOMMQ with Q=sigma_r:'
      write(17,'(/a)') ' test HCOMMQ with Q=sigma_r:'
c
      sumq   = 0.0
      sumq2  = 0.0
      ipostr = 0
c
      do ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
c
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         vi     = sqrt(vsqf(ist))
         ui     = sqrt(one-vi*vi)
c
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
            jstmax = nst(1)+nst(2)
         endif
c
         do jst=ist,jstmax      ! <-> beta
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               vj    = sqrt(vsqf(jst))
               uj    = sqrt(one-vj*vj)
c
               do ir=0,nr
                  do iz=0,nz
c                     z  = dz*(iz+sz-icm)
c                     z2 = z*z
                     wrk1(iz,ir) = ! Q=sigma_r
     &                    ( psi(iz,ir,1,ist)*psi(iz,ir,2,jst)
     &                    + psi(iz,ir,2,ist)*psi(iz,ir,1,jst))
                  enddo
               enddo
c
c     integrating and storing
c
               qmatre = f1int(wrk1)
               qresp(ipostr) = qmatre*(ui*vj+uj*vi)
Cr--------------------------------------------------------
CR ??? fehlt da nicht ein Faktor 2 ?
               if(ist.eq.jst) then
                  sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
                  sumq  = qmatre*vi*vi + sumq
               else
                  sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
               endif
CR
Cr--------------------------------------------------------
c
            endif
         enddo
      enddo
      iposmx = ipostr
      write(7,'(a,2f8.3,a,2(1pg12.4))')  
     &  ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
      write(17,'(a,2f8.3,a,2(1pg12.4))')  
     &  '# Sabtst: for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
c
c     testing by calling commutator
c
      call hcommQ(qresp,presp,ptest,icase)
c     
      return
      end 
c
c-----------------------------------------------------------------------
CR210499
c
      subroutine hcDT
c
c     Test-routine for Commutator [H,Q] of Hamiltonian 
c     with the case Q=z
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
c     0. step:
c     
c
      if(itest.eq.0) then       ! only for testing
c
c     test particle number
c
         do ir=0,nr
            do iz=0,nz
c               wrk1(iz,ir)  = 0.0
c               wrk2(iz,ir)  = 0.0
               wrk4(iz,ir) = rho(iz,ir,1) + rho(iz,ir,2)
            enddo
         enddo
         write(17,'(/a,1(1pg12.4))')
     &        '> hcDT: test output : Tr{1}  = #particles =',
     &        f1int(wrk4)
c     
         do i=1,2               ! i=1 -> b1=0, i=2 -> b1 not 0
c
c     setting force parametres
c
            if (i.eq.1)  then
               b1  = 0.0
               b1p = 0.0
               b1c = 0.0
            else
               b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
               b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
               b1c = b1-b1p
            endif
c
            write(17,'(a,2(1pg12.4))')
     &           '> hcDT: test case: b1, b1p  = ',b1,b1p
c
c
c     calculating test masse
c 
c     mass2   <->   Q =  z        => (Nabla Q)^2 =  1
c
c     making integrant in workspaces
c
            do ir=0,nr
               do iz=0,nz
c     
c     coordinates
c
                  r  = dr*(ir+sr)
                  z  = dz*(iz+sz-icm)
c     
c    h2meffq +2.0 in 2.zeile ????
c
c
                  h2meffp =
     &                 h2m(1)+b1c*rho(iz,ir,1)+b1*rho(iz,ir,2) 
                  h2meffn =
     &                 h2m(2)+b1c*rho(iz,ir,2)+b1*rho(iz,ir,1)   
c
                  wrk1(iz,ir) = 
     &                 h2meffp*rho(iz,ir,1) + h2meffn*rho(iz,ir,2)
c
               enddo
            enddo
            
c
c     final integrating, with factor
c
            xmass = f1int(wrk1)
c
            write(17,'(a,1(1pg12.4))')
     &           '> hcDT: xmass (inv) = ',
     &           xmass     
c
c     all done
c
         enddo
      endif
      return
      end 
c
c-----Dabtst-----------------------------------------------------------
CR200499
c
      subroutine Dabtst(icase)
c
c     calculates the matrix-elements Q_\alpha\beta of Q=z for 
c     testing and stores results on 'qresp' where it can be processed
c     further. 
c
c     icase is only used for the call of hcommQ
c
c     Dipole-Test
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
c
      write(7,'(/a)') ' test HCOMMQ with Q=z:'
      write(17,'(/a)') ' test HCOMMQ with Q=z:'
c
      sumq   = 0.0
      sumq2  = 0.0
      ipostr = 0
c
      do ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
c
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         vi     = sqrt(vsqf(ist))
         ui     = sqrt(one-vi*vi)
c
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
            jstmax = nst(1)+nst(2)
         endif
c
         do jst=ist,jstmax      ! <-> beta
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               vj    = sqrt(vsqf(jst))
               uj    = sqrt(one-vj*vj)
c
               do ir=0,nr
                  do iz=0,nz
                     z  = dz*(iz+sz-icm)
c  
                     wrk1(iz,ir) = ! Q=z
     &                    ( psi(iz,ir,1,ist)*z*psi(iz,ir,1,jst)
     &                    + psi(iz,ir,2,ist)*z*psi(iz,ir,2,jst))
                  enddo
               enddo
c
c     integrating and storing
c
               qmatre = f1int(wrk1)
               qresp(ipostr) = qmatre*(ui*vj+uj*vi)
Cr----------------------------------------------------------
CR fehlt da nicht ein Faktor 2?
               if(ist.eq.jst) then
                  sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
                  sumq  = qmatre*vi*vi + sumq
               else
                  sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
               endif
CR----------------------------------------------------------
c
            endif
         enddo
      enddo
      iposmx = ipostr
      write(7,'(a,2f8.3,a,2(1pg12.4))')  
     &  ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
      write(17,'(a,2f8.3,a,2(1pg12.4))')  
     &  '# Dabtst for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
c
c     testing by calling commutator
c
      call hcommQ(qresp,presp,ptest,icase)
c     
      return
      end 
C#######################################################################
CR250599
c-----PSabtst-----------------------------------------------------------
c
      subroutine PSabtst(icase)
c
c     calculates the matrix-elements Q_\alpha\beta of Q= p_r sigma_r for 
c     testing and stores results on 'qresp' where it can be processed
c     further. 
c
c     icase is only used for the call of hcommQ
c
c     p_r*sigma_r Test
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c-----------------------------------------------------------------------
c
c
      write(7,'(/a)') ' test HCOMMQ with Q= p_r sigma_r:'
      write(17,'(/a)') ' test HCOMMQ with Q= p_r sigma_r:'
c
      sumq   = 0.0
      sumq2  = 0.0
      ipostr = 0
c
      do ist=1,nst(1)+nst(2)    ! over all particles <-> alpha
c
         mist   = mpsi(ist)
         npist  = nparz(1,ist)
         vi     = sqrt(vsqf(ist))
         ui     = sqrt(one-vi*vi)
c
         if(ist.le.nst(1)) then ! protons
            iq     = 1
            jstmax = nst(1)
         else                   ! neutrons
            iq     = 2
            jstmax = nst(1)+nst(2)
         endif
c
         do jst=ist,jstmax      ! <-> beta
            if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
               ipostr = 1 + ipostr
               vj    = sqrt(vsqf(jst))
               uj    = sqrt(one-vj*vj)
c
               do ir=0,nr
                  do iz=0,nz
                     z  = dz*(iz+sz-icm)
c  
                     wrk1(iz,ir) = ! Q= p_r sigma_r
c p_r = -i d_r
     &                    -0.5*( 
     &                      dpsir(iz,ir,1,ist) * psi(iz,ir,2,jst)
     &                    + psi(iz,ir,1,ist)   * dpsir(iz,ir,2,jst)
     &                    + dpsir(iz,ir,2,ist) * psi(iz,ir,1,jst)
     &                    + psi(iz,ir,2,ist)   * dpsir(iz,ir,1,jst)
     &                    )
                  enddo
               enddo
c
c     integrating and storing
c
               qmatre = f1int(wrk1)
               qresp(ipostr) = qmatre*(ui*vj+uj*vi)
Cr----------------------------------------------------------
CR fehlt da nicht ein Faktor 2?
               if(ist.eq.jst) then
                  sumq2 = qresp(ipostr)*qresp(ipostr) + sumq2
                  sumq  = qmatre*vi*vi + sumq
               else
                  sumq2 = 2.0*qresp(ipostr)*qresp(ipostr) + sumq2
               endif
CR----------------------------------------------------------
c
            endif
         enddo
      enddo
      iposmx = ipostr
      write(7,'(a,2f8.3,a,2(1pg12.4))')  
     &  ' for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
      write(17,'(a,2f8.3,a,2(1pg12.4))')  
     &  '# PSabtst for t1,t2=',t1,t2,': <Q>,<Q**2>=',sumq,sumq2
c
c     testing by calling commutator
c
      call hcommQ(qresp,presp,ptest,icase)
c     
      return
      end 
C#######################################################################
CR250599
c-----Oabtst-----------------------------------------------------------
c
      subroutine Oabtst(icase)
c
c     Major-Routine to calculate Test-Operator-Matrixelements
c 
c     with cases 
c
c     Qabtst:   icase = 1,2,3
c
c                     1: Q = 2z^2
c                     2: Q = -r^2
c                     3: Q = 2z^2-r^2
c
c     Dabtst:   icase = 4: Q = z
c
c     Sabtst:   icase = 5: Q = sigma_r
c     PSabtest: icase = 6: Q = p_r sigma_r
c
      include 'paramr.inc'
c
      if((icase.ge.1).and.(icase.le.3)) then
         call Qabtst(icase)
      else if(icase.eq.4) then
         call Dabtst(icase)
      else if(icase.eq.5) then
         call Sabtst(icase)
      else if(icase.eq.6) then
         call PSabtst(icase)
      else if(icase.gt.6) then
         STOP 'Oabtst: unkown icase' 
      endif
c
      return
      end 
C#######################################################################
CR250599
c-----hcOT-----------------------------------------------------------
c
      subroutine hcOT(icase)
c
c     Major-Routine to call hcommQ-Test-Routines
c 
c     with cases 
c
c     hcQT: icase = 1,2,3
c
c                     1: Q = 2z^2
c                     2: Q = -r^2
c                     3: Q = 2z^2-r^2
c
c     hcDT: icase = 4  : Q = z
c     icase 5: sigma_r
c     icase 6: pr_sigma_r
c
      include 'paramr.inc'
c
      if((icase.ge.1).and.(icase.le.3)) then
         call hcQT(icase)
      else if(icase.eq.4) then
         call hcDT
      else if(icase.gt.6) then
         STOP 'hcOT: unkown icase' 
      endif
c
      return
      end 
C#######################################################################
CR250599
c-----ptstcur-----------------------------------------------------------
c
      subroutine ptstcur(icase)
c
c     subroutine of hcommQ: output of test currents on unit 17
c 
c     with cases 
c
c     ptstcur: icase = 1,2,3
c
c                     1: Q = 2z^2
c                     2: Q = -r^2
c                     3: Q = 2z^2-r^2
c
c     ptstcur: icase = 4  Q = z
c     ptstcur: icase = 5  Q = sigma_r
c     ptstcur: icase = 6  Q = p_r sigma_r ! todo
c
c
      include 'paramr.inc'      
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
      include 'crankin.inc'
c
c----------------------------------------------------------------------
c
      if(icase.le.3) then       ! Qadrupol
c
         if(icase.eq.1) then    ! plots of currents, sigma_Q,\phi==0
            write(17,'(//a)') 
     &           '# ptstcur: plot of currents for Q=2z^2'
         else if(icase.eq.2) then 
            write(17,'(//a)') 
     &           '# ptstcur: plot of currents for Q=-r^2'
         else
            write(17,'(//a)') 
     &           '# ptstcur: plot of currents for Q=2z^2-r^2' 
         endif
c     
         iz0 = (nz+1)/2+1       ! this is z=0+dz
c
         write(17,'(/a)')
     &        '# test curqz for ir=1 :',
     &        '# iz  z  curqz_p  curqz/rho_p  curqz_n  curqz/rho_n '
         do iz=0,nz
            write(17,'(1x,i3,f8.2,4(1pg12.4))')
     &           iz,(iz+sz-icm)*dz,curqz(iz,1,1),
     &           curqz(iz,1,1)/rho(iz,1,1),
     &           curqz(iz,1,2),curqz(iz,1,2)/rho(iz,1,2)
         enddo
crold         write(17,'(/a)')
crold     &        '# test curqr for iz=0 :',
crold     &        '# ir  r  curqr  curqr/rho '
crold         do ir=0,nr
crold            write(17,'(1x,i3,f8.2,4(1pg12.4))')
crold     &           ir,(ir+sr)*dr,curqr(0,ir,1),curqr(0,ir,1)/rho(0,ir,1),
crold     &           curqr(0,ir,2),curqr(0,ir,2)/rho(0,ir,2)
crold         enddo
         write(17,'(/a)')
     &        '# test curqr for z=0+dz :',
     &        '# ir  r  curqr_p  curqr/rho_p  curqr_n  curqr/rho_n '
         do ir=0,nr
            write(17,'(1x,i3,f8.2,4(1pg12.4))')
     &           ir,(ir+sr)*dr,curqr(iz0,ir,1),
     &           curqr(iz0,ir,1)/rho(iz0,ir,1),
     &           curqr(0,ir,2),curqr(iz0,ir,2)/rho(iz0,ir,2)
         enddo
C     
CR030399 test currents j_Q,q Q=2z^2-r^2
CR210499
crold         write(17,'(/a/)') '# Q=2z^2-r^2 '
crold         write(17,'(/a,2g12.6)') '# dr,dz= ',
crold     &        dr,dz
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,p;r /(rho_p*dr))|iz=0,ir=1 = -2r * ',
crold     &        curqr(0,iact,1) / (-2.0 * rho(0,iact,1) * dr)
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,p;z /(rho_p*dz))|ir=0,iz=1 =  4z * ',
crold     &        curqz(iact,0,1) / ( 4.0 * rho(iact,0,1) * dz)
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,n;r /(rho_n*dr))|iz=0,ir=1 = -2r * ',
crold     &        curqr(0,iact,2) / (-2.0 * rho(0,iact,2) * dr)
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,n;z /(rho_n*dz))|ir=0,iz=1 =  4z * ',
crold     &        curqz(iact,0,2) / ( 4.0 * rho(iact,0,2) * dz)
crold
         if(icase.eq.1) then    ! plots of currents, sigma_Q,\phi==0
            write(17,'(//a)') 
     &           '# ptstcur: Q=2z^2'
         else if(icase.eq.2) then 
            write(17,'(//a)') 
     &           '# ptstcur: Q=-r^2'
         else
            write(17,'(//a)') 
     &           '# ptstcur: Q=2z^2-r^2' 
         endif
c
         write(17,'(/a,2g12.6)') '# dr,dz= ',
     &        dr,dz
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 = -2r * ',
     &        curqr(iz0,1,1) / (-2.0 * rho(iz0,1,1) * dr)
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  4z * ',
     &        curqz(iz0,1,1) / ( 4.0 * rho(iz0,1,1) * dz)
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 = -2r * ',
     &        curqr(iz0,1,2) / (-2.0 * rho(iz0,1,2) * dr)
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  4z * ',
     &        curqz(iz0,1,2) / ( 4.0 * rho(iz0,1,2) * dz)
c         
      else if(icase.eq.4) then ! Q=z
         iz0 = (nz+1)/2+1 ! this is z=0+dz
crold         write(17,'(/a)') 
crold     &        '# ptstcur: plot of currents for Q=z'
croldc     
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,p;r /(rho_p*dr))|ir=1,iz=0 =  0 ',
crold     &        curqr(0,1,1) 
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,p;z /(rho_p*dz))|ir=0,iz=1 =  rho /2 ',
crold     &        curqz(iact,0,1) / ( rho(iact,0,1) * dz) * 1.0
crold         write(17,'(a,1g12.6)') 
crold     &        '#(j_Q,n;r /(rho_n*dr))|ir=1,iz=0 =  0 ',
crold     &        curqr(0,1,2) 
crold         write(17,'(a,1g12.6)') 
crold     &        '# (j_Q,n;z /(rho_n*dz))|ir=0,iz=1 =  rho /2 ',
crold     &        curqz(iact,0,2) / ( rho(iact,0,2) * dz) * 1.0
         write(17,'(/a)') 
     &        '# ptstcur: plot of currents for Q=z'
c     
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 =  0 ',
     &        curqr(iz0,1,1) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  rho /2 ',
     &        curqz(iz0,1,1) / ( rho(iz0,1,1) * dz) * 1.0
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 =  0 ',
     &        curqr(iz0,1,2) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  rho /2 ',
     &        curqz(iz0,1,2) / ( rho(iz0,1,2) * dz) * 1.0
c     
      else if(icase.eq.5) then
         iz0 = (nz+1)/2+1 ! this is z=0+dz
         write(17,'(/a)') 
     &        '# ptstcur: plot of currents for Q=sigma_r'
CR040399 test currents sigma_Q,q Q=sigma_r 
c     for protons ->wrk1, neutons ->wrk2
c
c     cleaning up
         do ir=0,nr
            do iz=0,nz
               wrk1(iz,ir) = 0.0
               wrk2(iz,ir) = 0.0
               wrk3(iz,ir) = 0.0
               wrk4(iz,ir) = 0.0
            enddo
         enddo
c     
         do ist=1,nst(1)+nst(2)
            vi     = sqrt(vsqf(ist))
            ui     = sqrt(one-vi*vi)
            mist   = mpsi(ist)
            npist  = nparz(1,ist)
            if(ist.le.nst(1)) then ! protons
               iq     = 1
               eferm  = fer(1)
               jstmax = nst(1)
            else                ! neutrons
               iq     = 2
               eferm  = fer(2)
               jstmax = nst(1)+nst(2)
            endif
            do jst=ist,jstmax
               if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
                  vj     = sqrt(vsqf(jst))
                  uj     = sqrt(one-vj*vj)
                  if(ist.ne.jst) then 
                     factor = (ui*ui*vj*vj-vi*vi*uj*uj) 
                     do ir=0,nr
                        do iz=0,nz
                           if(iq.eq.1) then ! protons
                              wrk1(iz,ir) = -2.0  
     &                             *(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                             - psi(iz,ir,2,ist)*psi(iz,ir,2,jst))
     &                             *factor  +  wrk1(iz,ir)
                           else ! neutrons
                              wrk1(iz,ir) = -2.0  
     &                             *(psi(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                             - psi(iz,ir,2,ist)*psi(iz,ir,2,jst))
     &                             *factor  +  wrk2(iz,ir)
                           endif
                        enddo
                     enddo
                  endif
               endif
            enddo
         enddo
C     
C     
         write(17,'(/a)')
     &        '# test spinph for ir=1 :',
     &        '# iz  z  spinph_p  exact_p spinph_n  exact_n '
         do iz=0,nz
            write(17,'(1x,i3,f8.2,4(1pg12.4))')
     &           iz,(iz+sz-icm)*dz,spinph(iz,1,1),wrk1(iz,1),
     &           spinph(iz,1,2),wrk2(iz,1)
         enddo
         write(17,'(a)')
     &        '# test curqr for z=0+dz :',
     &        '# ir  r  curqr_p  curqr/rho_p  curqr_n  curqr/rho_n '
         do ir=0,nr
            write(17,'(1x,i3,f8.2,4(1pg12.4))')
     &           ir,(ir+sr)*dr,spinph(iz0,ir,1),wrk1(iz0,ir),
     &              spinph(iz0,ir,2),wrk2(iz0,ir)
         enddo
c     
         do ir=0,nr
            do iz=0,nz
               wrk3(iz,ir)=(wrk1(iz,ir)-spinph(iz,ir,1))**2
               wrk4(iz,ir)=(wrk2(iz,ir)-spinph(iz,ir,2))**2
            enddo
         enddo
C
         write(17,'(/a)') 
     &        '# ptstcur: plot of currents for Q=sigma_r'
c     
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 =  0 ',
     &        curqr(iz0,1,1) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  0 ',
     &        curqz(iz0,1,1) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 =  0 ',
     &        curqr(iz0,1,2) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  0 ',
     &        curqz(iz0,1,2)     
C
         xst1 = f1int(wrk3)
         xst2 = f1int(wrk4)
C     
         write(17,'(/a,1g12.6)') '# sigma_Q=sigma_r,p: 0?=xst1 ', xst1
         write(17,'(a,1g12.6)') '# sigma_Q=sigma_r,n: 0?=xst2 ', xst2
c     
c     cleaning up   (not important, but more clean)
         do ir=0,nr
            do iz=0,nz
               wrk1(iz,ir) = 0.0
               wrk2(iz,ir) = 0.0
               wrk3(iz,ir) = 0.0
               wrk4(iz,ir) = 0.0
            enddo
         enddo
C     
      else if(icase.eq.6) then
         iz0 = (nz+1)/2+1 ! this is z=0+dz
         write(17,'(/a)') 
     &        '# ptstcur: plot of currents for Q=p_r sigma_r'
CR250599 test currents sigma_Q,q Q=p_r sigma_r 
c     for protons ->wrk1, neutons ->wrk2
c
c     cleaning up
         do ir=0,nr
            do iz=0,nz
               wrk1(iz,ir) = 0.0
               wrk2(iz,ir) = 0.0
               wrk3(iz,ir) = 0.0
               wrk4(iz,ir) = 0.0
            enddo
         enddo
c     
         do ist=1,nst(1)+nst(2)
            vi     = sqrt(vsqf(ist))
            ui     = sqrt(one-vi*vi)
            mist   = mpsi(ist)
            npist  = nparz(1,ist)
            if(ist.le.nst(1)) then ! protons
               iq     = 1
               eferm  = fer(1)
               jstmax = nst(1)
            else                ! neutrons
               iq     = 2
               eferm  = fer(2)
               jstmax = nst(1)+nst(2)
            endif
            do jst=ist,jstmax
               if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
                  vj     = sqrt(vsqf(jst))
                  uj     = sqrt(one-vj*vj)
                  if(ist.ne.jst) then 
                     factor = (ui*ui*vj*vj-vi*vi*uj*uj) 
                     do ir=0,nr
                        do iz=0,nz
c
c sigma_Q_q = - 2 i sum phi^+ p_r sigma_z phi 
c
                           if(iq.eq.1) then ! protons
                              wrk1(iz,ir) = -2.0*0.5  
     &                             *(dpsir(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                             - dpsir(iz,ir,2,ist)*psi(iz,ir,2,jst)
     &                             + psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                             - psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                             )
     &                             *factor  +  wrk1(iz,ir)
                           else ! neutrons
                              wrk1(iz,ir) = -2.0*0.5
     &                             *(dpsir(iz,ir,1,ist)*psi(iz,ir,1,jst)
     &                             - dpsir(iz,ir,2,ist)*psi(iz,ir,2,jst)
     &                             + psi(iz,ir,1,ist)*dpsir(iz,ir,1,jst)
     &                             - psi(iz,ir,2,ist)*dpsir(iz,ir,2,jst)
     &                             )
     &                             *factor  +  wrk2(iz,ir)
                           endif
                        enddo
                     enddo
                  endif
               endif
            enddo
         enddo
C     
C     
         write(17,'(/a)')
     &        '# test spinph for ir=1 :',
     &        '# iz  z  spinph_p  exact_p spinph_n  exact_n '
         do iz=0,nz
            write(17,'(1x,i3,f8.2,4(1pg12.4))')
     &           iz,(iz+sz-icm)*dz,spinph(iz,1,1),wrk1(iz,1),
     &           spinph(iz,1,2),wrk2(iz,1)
         enddo
         write(17,'(a)')
     &        '# test curqr for z=0+dz :',
     &        '# ir  r  curqr_p  curqr/rho_p  curqr_n  curqr/rho_n '
         do ir=0,nr
            write(17,'(1x,i3,f8.2,4(1pg12.4))')
     &           ir,(ir+sr)*dr,spinph(iz0,ir,1),wrk1(iz0,ir),
     &              spinph(iz0,ir,2),wrk2(iz0,ir)
         enddo
c     
         do ir=0,nr
            do iz=0,nz
               wrk3(iz,ir)=(wrk1(iz,ir)-spinph(iz,ir,1))**2
               wrk4(iz,ir)=(wrk2(iz,ir)-spinph(iz,ir,2))**2
            enddo
         enddo
C
         write(17,'(/a)') 
     &        '# ptstcur: plot of currents for Q=p_r sigma_r'
c     
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;r /(rho_p*dr))|z=0+dz,ir=1 =  0 ',
     &        curqr(iz0,1,1) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,p;z /(rho_p*dz))|z=0+dz,ir=1 =  0 ',
     &        curqz(iz0,1,1) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;r /(rho_n*dr))|z=0+dz,ir=1 =  0 ',
     &        curqr(iz0,1,2) 
         write(17,'(a,1g12.6)') 
     &        '# (j_Q,n;z /(rho_n*dz))|z=0+dz,ir=1 =  0 ',
     &        curqz(iz0,1,2)     
C
         xst1 = f1int(wrk3)
         xst2 = f1int(wrk4)
C     
         write(17,'(/a,1g12.6)') '# sigma_Q=pr sigma_r,p: 0?=xst1 ',xst1
         write(17,'(a,1g12.6)') '# sigma_Q=pr sigma_r,n: 0?=xst2 ', xst2
c     
c     cleaning up   (not important, but more clean)
         do ir=0,nr
            do iz=0,nz
               wrk1(iz,ir) = 0.0
               wrk2(iz,ir) = 0.0
               wrk3(iz,ir) = 0.0
               wrk4(iz,ir) = 0.0
            enddo
         enddo
C     
      endif
      return
      end
CR250599 end of currents test output
c############################################################################
