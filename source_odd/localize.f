c-----localize------------------------------------------------------------
c
      subroutine localize(iwunit)
c
c     computes localization criterion
c
      include 'paramr.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'observ.inc'
c
      logical ttest
      data ttest/.false./
c
c
c----------------------------------------------------------------------
c
c     interpolate to standard grid points and sum up
c 
      call nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,one)
      call nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,one)
      call nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,one)
      call nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,one)

      if(ttest) write(6,'(a)')   'localization:',
     &      '  iz  ir    z     r      c      c_tilde '
      write(iwunit,'(a,2f8.3/a)')  
     &      '# beta2,beta3=',beta2,beta3,
     &      '#  iz  ir    z     r      c      c_tilde '
      ttffac = 0.6*(3.0*pi**2)**0.6666666667
      do iq=1,2
      do ir=0,nrd
      do iz=0,nzd
        rhoact = rho(iz,ir,iq)
        if(rhoact.gt.SMAL) then
          clocal = tau(iz,ir,iq) 
     &            -(psip(iz,ir,iq)**2+psin(iz,ir,iq)**2)/rhoact
          chom   = ttffac*rhoact**1.666666667
          ctilde = 1.0D0/(1.0D0+(clocal/chom)**2)
        else
          clocal = 1.0D-20
          chom   = 1.0D-10
          curvt  = 0.0D0
          ctilde = 0.0D0
        endif
        zpoint  = ((iz)+sz-icm)*dz
        if(ttest) write(6,'(1x,2i5,2f8.3,3(1pg13.5))')
     &      iz,ir,zpoint,(ir-0.5D0)*dr,clocal,ctilde
        write(iwunit,'(1x,2i5,2f8.3,3(1pg13.5))')
     &      iz,ir,zpoint,(ir-0.5D0)*dr,clocal,ctilde,rhoact
      enddo
      write(iwunit,'(5x,a)') ' '
      enddo
      write(iwunit,'(5x,a)') ' '
      write(iwunit,'(5x,a)') ' '
      enddo
c
c
      return
      end

