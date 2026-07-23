                  resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,2,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt1*(
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
                  resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt1*(
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
                  resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*drphi(iz,ir,2,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
                  resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) + wgt1*(
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,2,jst)*dzphi(iz,ir,2,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,2,jst)*  phi(iz,ir,2,ist) )
c
c
                  resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt2*(
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
                  resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt1*(
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
                  resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt1*(
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
                  resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt1*(
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
               do ir=0,nr
               do iz=0,nz
c     
                  resp_spnr(iz,ir,iso) = resp_spnr(iz,ir,iso) + wgt1*(
     &                 - phi(iz,ir,2,jst)*phi(iz,ir,2,ist)
     &                 + phi(iz,ir,1,jst)*phi(iz,ir,1,ist) )
c
                  resp_spnp(iz,ir,iso) = resp_spnp(iz,ir,iso) + wgt1*(
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
                  resp_skdp(iz,ir,iso) = resp_skdp(iz,ir,iso) + wgt1*(
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
                  resp_curr(iz,ir,iso) = resp_curr(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*drphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*drphi(iz,ir,2,ist)
     &                 + drphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - drphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
                  resp_curp(iz,ir,iso) = resp_curp(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*dpphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dpphi(iz,ir,2,ist)
     &                 + dpphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dpphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c     
                  resp_curz(iz,ir,iso) = resp_curz(iz,ir,iso) + wgt1*(
     &                 -   phi(iz,ir,2,jst)*dzphi(iz,ir,1,ist)
     &                 +   phi(iz,ir,1,jst)*dzphi(iz,ir,2,ist)
     &                 + dzphi(iz,ir,2,jst)*  phi(iz,ir,1,ist)
     &                 - dzphi(iz,ir,1,jst)*  phi(iz,ir,2,ist) )
c
c
                  resp_paio(iz,ir,iso) = resp_paio(iz,ir,iso) + wgt2*(
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
                  resp_rotjp(iz,ir,iso) = resp_rotjp(iz,ir,iso) + wgt1*(
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
                  resp_rotsr(iz,ir,iso) = resp_rotsr(iz,ir,iso) + wgt1*(
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
                  resp_rotsz(iz,ir,iso) = resp_rotsz(iz,ir,iso) + wgt1*(
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
