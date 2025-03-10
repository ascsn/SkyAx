c-----wspot -----------------------------------------------------------
c
      subroutine wspot()
*
***********************************************************************
*                                                                     *
*        compute the mean field hamiltonian .                         *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'levels.inc'
c
      dimension h2mnist(2) 
      data h2mnist/20.74982,20.72126/
c
      logical tfirst
      data tfirst/.true./
c
c---------------------------------------------------------------------
c
      if(tfirst) then
        write(6,'(a)')
     &  ' WS potential with the parameters:',
     &  ' rad0   strenght  kapp  width  b0  b2  b4 '
        do iq=1,2
          write(6,'(1x,7f9.4)') wsrad0(iq),ws0(iq),wskap(iq),
     &       wswid(iq),wsb0(iq),wsb2(iq),wsb4(iq)
        enddo
      endif
c
c     the local potentials, the (constant) inverse mass, 
c     and the charge distribution
c
      yf2  = 0.25 * sqrt(5.0/pi)
      yf4  = 3.0 / 16.0 / sqrt(pi)
      anum13 = anumb**0.3333333D0
      do ir=0,nr
        r     = (sr+ir)*dr
        r2    = r*r
        do iz=0,nz
          z     = (sz+iz-icm)*dz
          z2    = z*z
          rad2  = r2+z2
          rad   = sqrt(rad2)
          rad2  = max(rad2,smal)
          y20p  = yf2*(two*z2-r2)/rad2
          y40p  = yf4*(8.0d0*z2*z2 + 3.0d0 *r2*r2 - 24.0d0 *z2*r2)
     &           /(rad2*rad2)
          do iq=1,2
            wsrad = wsrad0(iq)*anum13
            wsradef = wsrad*(one+wsb0(iq)+wsb2(iq)*y20p+wsb4(iq)*y40p)
            pot(iz,ir,iq) = ws0(iq)/(one+exp((rad-wsradef)/wswid(iq)))
            bmass(iz,ir,iq) = h2mnist(iq)
          enddo
          wsrad = wsrad0(1)*anum13
          wsradef = wsrad*(one+wsb0(1)+wsb2(1)*y20p+wsb4(1)*y40p)
          crho(iz,ir) = one/(one+exp((rad-wsradef)/wswid(1)))
        enddo
      enddo
      if(tfirst) write(6,'(a)') ' WS potentials done '
      call testpr(pot(0,0,1),1,10,1,10,' pot_p')
      call testpr(pot(0,0,2),1,10,1,10,' pot_n')
c
c     the spin-orbit potentials
c
      call nablaz_scalar(pot(0,0,1),potls(0,0,1,1), 1,1,wskap(1))
      call nablaz_scalar(pot(0,0,2),potls(0,0,2,1), 1,1,wskap(2))
      call nablar_scalar(pot(0,0,1),potls(0,0,1,2), 0,1,wskap(1))
      call nablar_scalar(pot(0,0,2),potls(0,0,2,2), 0,1,wskap(2))
      
cPK   do iq=1,2
cPK     call mulr(derr1(0,0,1),pot(0,0,iq),potls(0,0,iq,1),one,0)
cPK     call mulz(derz1(0,0,1),pot(0,0,iq),potls(0,0,iq,2),one,0)
cPK     do ir=0,nr
cPK     do iz=0,nz
cPK       potls(iz,ir,iq,1) = wskap(iq)*potls(iz,ir,iq,1)
cPK       potls(iz,ir,iq,2) = wskap(iq)*potls(iz,ir,iq,2)
cPK     enddo
cPK     enddo
cPK   enddo

      if(tfirst) write(6,'(a)') ' WS-ls potentials done '
c
c     the Coulomb potential
c
ctest      return
      xnorm = (npart(1)-one)/f1int(crho)
      write(7,'(a,1pg12.4)') ' norm of CRHO=',xnorm
      call testpr(crho,1,10,1,10,' crho')
      call coulmb(0,0)
      do ir=0,nr
      do iz=0,nz
        coulpt(iz,ir) = xnorm*coulpt(iz,ir)
        pot(iz,ir,1) = pot(iz,ir,1) + coulpt(iz,ir)
      enddo
      enddo
      write(6,'(a)') ' WS Coulomb done '
c
      call testpr(coulpt,1,10,1,10,' Coulomb')
      call testpr(potls(0,0,1,1),1,10,1,10,' potls_p1')
      call testpr(potls(0,0,1,2),1,10,1,10,' potls_p2')
      call testpr(potls(0,0,2,1),1,10,1,10,' potls_n1')
      call testpr(potls(0,0,2,2),1,10,1,10,' potls_n1')
      tfirst = .false.
c
      return
      end
