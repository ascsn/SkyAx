c
c        find [H,Q] = i B P^{res} by iteration
         do iter=1,2000

         call upd_derive_spinor(resQ(0,0,1,1,1),drresQ(0,0,1,1,1),
     &                          dzresQ(0,0,1,1,1),dpresQ(0,0,1,1,1),
c     &                          ddresQ(0,0,1,1,1),+1,0)
c         call upd_derive_spinor(resQ(0,0,1,1,2),drresQ(0,0,1,1,2),
c     &                          dzresQ(0,0,1,1,2),dpresQ(0,0,1,1,2),
c     &                          ddresQ(0,0,1,1,2),+1,0)

c         call comm_HP_res(resQ(0,0,1,1,1),drresQ(0,0,1,1,1),
c     &                    dzresQ(0,0,1,1,1),dpresQ(0,0,1,1,1),
c     &                    ddresQ(0,0,1,1,1),resHQ(0,0,1,1,1))
c         call comm_HP_res(resQ(0,0,1,1,2),drresQ(0,0,1,1,2),
c     &                    dzresQ(0,0,1,1,2),dpresQ(0,0,1,1,2),
c     &                    ddresQ(0,0,1,1,2),resHQ(0,0,1,1,2))

c         buff=0.0D0
c        do ist=1,nst(1)+nst(2)
c
c            do ir=0,nr
c               do iz=0,nz
c                  LRwrk1_spinor(iz,ir,1) = damp(iz,ir)*(
c     &               resHQ(iz,ir,1,ist,1) - resP(iz,ir,1,ist,1))
c                  LRwrk1_spinor(iz,ir,2) = damp(iz,ir)*(
c     &               resHQ(iz,ir,2,ist,1) - resP(iz,ir,2,ist,1))
c                  LRwrk2_spinor(iz,ir,1) = damp(iz,ir)*(
c     &               resHQ(iz,ir,1,ist,2) - resP(iz,ir,1,ist,2))
c                  LRwrk2_spinor(iz,ir,2) = damp(iz,ir)*(
c     &               resHQ(iz,ir,2,ist,2) - resP(iz,ir,2,ist,2))
c               enddo
c            enddo
ccc
c            buff = buff+
c     &             wfovlp(LRwrk1_spinor(0,0,1),LRwrk1_spinor(0,0,1))+
c     &             wfovlp(LRwrk2_spinor(0,0,1),LRwrk2_spinor(0,0,1))
cc     
c           call damp_nucleon(LRwrk1_spinor(0,0,1),resHQ(0,0,1,ist,1),
c     &                        nparz(1,ist),mpsi(ist),0.005D0)
c           call damp_nucleon(LRwrk2_spinor(0,0,1),resHQ(0,0,1,ist,2),
c     &                        nparz(1,ist),mpsi(ist),0.005D0)
c
c             do ir=0,nr
c               do iz=0,nz
c                  resQ(iz,ir,1,ist,1) =
c     &                 resQ(iz,ir,1,ist,1) - resHQ(iz,ir,1,ist,1)
c                  resQ(iz,ir,2,ist,1) =
c     &                 resQ(iz,ir,2,ist,1) - resHQ(iz,ir,2,ist,1)
c                  resQ(iz,ir,1,ist,2) =
c     &                 resQ(iz,ir,1,ist,2) - resHQ(iz,ir,1,ist,2)
c                  resQ(iz,ir,2,ist,2) =
c     &                 resQ(iz,ir,2,ist,2) - resHQ(iz,ir,2,ist,2)
c               enddo
c            enddo
c            
c         enddo
c        
c
c         write(7,'(a,i4,1(1pg15.5))') 'CRNK_ITER: '
c     &                                ,iter,buff/(nst(1)+nst(2))
c         if (buff/(nst(1)+nst(2)).le.1E-2) goto 2001
c            
c      enddo
c 2001 continue
c      write(7,'(a,i4,1(1pg15.5))') '##########  PKtest:',
c     &                             iter,buff/(nst(1)+nst(2))
cc
