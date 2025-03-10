

c-----hcQT-----------------------------------------------------------
c
      subroutine hcQT
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
c     cleaning up workspaces (nicht noetig, ist ok)
c
         do ir=0,nr
            do iz=0,nz
               wrk1(iz,ir)  = 0.0
               wrk2(iz,ir)  = 0.0
c test ob er dichte kriegt, der code (ja)
c               write(7,'(a,2I4,2(1pg12.4))')
c     &           '> hcQT: test output : rho_p/n (iz,ir),  = ',
c     &              iz,ir, rho(iz,ir,1), rho(iz,ir,2)               
               wrk4(iz,ir) = rho(iz,ir,1) + rho(iz,ir,2)
            enddo
         enddo
         write(7,'(a,1(1pg12.4))')
     &        '> hcQT: test output : Tr{1}  = ',
     &        f1int(wrk4)
c     
         do i=1,2               ! i=1 -> b1=0, i=2 -> b1 not 0
c
c     setting force parametres
c
            if (i.eq.1)  then
               b1  = 0.0;
               b1p = 0.0;
            else
               b1  = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
               b1p = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
            endif
            b1c = b1-b1p
c
            write(7,'(a,2(1pg12.4))')
     &           '> hcQT: test case: b1, b1p  = ',b1,b1p
c
c
c     calculating test masses
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
c    h2meffq
                  h2meffp = (h2m(1) + b1c*rho(iz,ir,1) + b1*rho(iz,ir,2)) 
                  h2meffn = (h2m(2) + b1c*rho(iz,ir,2) + b1*rho(iz,ir,1)) 
c
                  wh2meff = h2meffp*rho(iz,ir,1)+h2meffn*rho(iz,ir,2)
c
                  wrk1(iz,ir) = wh2meff * ( z2 )
                  wrk2(iz,ir) = wh2meff * ( r2 )
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
            write(7,'(a,3(1pg12.4))')
     &           '> hcQT: mass1, mass2, mass3 (inv) = ',
     &           xmass1,xmass2,xmass3     
c
c     all done
c
         enddo
      endif
      return
      end 

