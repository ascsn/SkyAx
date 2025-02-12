c-----localize------------------------------------------------------------
c
      subroutine localize(outnam)
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
      character*3 outnam,tnum
      logical ttest
      data ttest/.false./
      data numcall/0/
c
c
c----------------------------------------------------------------------
c
      numcall = 1 + numcall
      open(78,status='scratch')
      write(78,'(i4)') 1000+numcall
      rewind(78)
      read(78,'(1x,a)') tnum
      close(78)
      iwunit = 78
c
c     outer loop over isospin
c
      do iq=1,2
      if(iq.eq.1) then
        istlow = 1
        istup  = nst(1)
      else
        istlow = nst(1)+1
        istup  = nst(1)+nst(2)
      endif
c
c     output for 'plotmtv' of density
c
      if(iq.eq.1) then
        open(iwunit,file='densip.'//tnum//'.'//outnam)
        write(iwunit,'(a)') '# plot density for protons:'
      else
        open(iwunit,file='densin.'//tnum//'.'//outnam)
        write(iwunit,'(a)') '# plot density for neutrons:'
      endif
      write(iwunit,'(a,2i4,2f8.3)')  
     &      '# nprot,nneut,beta2,beta3=',npart,beta2,beta3
      write(iwunit,'(a/a,i3,a,f7.2/a,i3,a,i3/a,i3,a,f7.2/a,i3)')
     &  '$ data=contour',
     &  '% nx=',nz+1,' xmin= 0.0 xmax=',nz*dz,
     &  '% vxmin=',int(nz*dz*(0.25-0.25*isym)),
     &      ' vxmax=',int(nz*dz*(0.75-0.25*isym)),
     &  '% ny=',nr+1,' ymin= 0.0 ymax=',nr*dr,
     &  '% vymax=',int(nr*dr/2.0)
      write(iwunit,'(a)')
     &  '% contstyle=2',
     &  '% nsteps=50',
     &  '% interp=2',
     &  '% xlabel="z [fm]"',
     &  '% ylabel="r [fm]"'
      if(iq.eq.1) then
        write(iwunit,'(a,2i4,a/a,2f6.3,a)')
     &  '% toplabel ="proton density for Z,N=',npart,'" '
      else
        write(iwunit,'(a,2i4,a/a,2f6.3,a)')
     &  '% toplabel ="neutron density for Z,N=',npart,'" '
      endif
      write(iwunit,'(a,2f6.3,a)')
     &  '% subtitle ="deformation beta_2,beta_3=',beta2,beta3,'"'
      write(iwunit,'(a)')
     &  '% grid=true',
     &  '% equalscale=true',
     &  '% fitpage=false',
     &  '% clip=true'
      write(iwunit,'(a,f6.3)')
     &  '% xyratio=',1.0*nr/nz
      write(iwunit,'(a)')
     &  '% cmin=0',
     &  '% cmax=0.1',
     &  '% cstep=0.01   '
      do ir=0,nr
        write(iwunit,'(10f8.5)') (rho(iz,ir,iq),iz=0,nz)
        write(iwunit,'(a)') '     '
      enddo
      close(iwunit)
c
c     prepare header for 'plotmtv' file
c
      if(iq.eq.1) then
        open(iwunit,file='localp.'//tnum//'.'//outnam)
        write(iwunit,'(a)') '# plot localization for protons:'
      else
        open(iwunit,file='localn.'//tnum//'.'//outnam)
        write(iwunit,'(a)') '# plot localization for neutrons:'
      endif
      write(iwunit,'(a,2i4,2f8.3)')  
     &      '# nprot,nneut,beta2,beta3=',npart,beta2,beta3
      write(iwunit,'(a/a,i3,a,f7.2/a,i3,a,i3/a,i3,a,f7.2/a,i3)')
     &  '$ data=contour',
     &  '% nx=',nz+1,' xmin= 0.0 xmax=',nz*dz,
     &  '% vxmin=',int(nz*dz*(0.25-0.25*isym)),
     &      ' vxmax=',int(nz*dz*(0.75-0.25*isym)),
     &  '% ny=',nr+1,' ymin= 0.0 ymax=',nr*dr,
     &  '% vymax=',int(nr*dr/2.0)
      write(iwunit,'(a)')
     &  '% contstyle=2',
     &  '% nsteps=50',
     &  '% interp=2',
     &  '% xlabel="z [fm]"',
     &  '% ylabel="r [fm]"'
      if(iq.eq.1) then
        write(iwunit,'(a,2i4,a/a,2f6.3,a)')
     &  '% toplabel ="proton localization for Z,N=',npart,'" '
      else
        write(iwunit,'(a,2i4,a/a,2f6.3,a)')
     &  '% toplabel ="neutron localization for Z,N=',npart,'" '
      endif
      write(iwunit,'(a,2f6.3,a)')
     &  '% subtitle ="deformation beta_2,beta_3=',beta2,beta3,'"'
      write(iwunit,'(a)')
     &  '% grid=true',
     &  '% equalscale=true',
     &  '% fitpage=false',
     &  '% clip=true'
      write(iwunit,'(a,f6.3)')
     &  '% xyratio=',1.0*nr/nz
      write(iwunit,'(a)')
     &  '% cmin=0',
     &  '% cmax=1',
     &  '% cstep=0.1   '
c      
c     compute derivative of density consistently
c
      do ir=0,nr
      do iz=0,nz
        wrk1(iz,ir) = 0.0D0                 !  d/dr rho
        wrk2(iz,ir) = 0.0D0                 !  d/dz rho
        rho(iz,ir,iq) = 0.0D0
        tau(iz,ir,iq) = 0.0D0
      enddo
      enddo
      do ist=istlow,istup
        fpair  = vsqf(ist)*wate(ist)
        mm     = mpsi(ist)
        xm2    = mm*mm
        xm2l   = (mm+1)*(mm+1)
        nparr  = 1+mod(mm,2)
        nparrc = 3-nparr
        if(ttest) write(iwunit,'(a,3i5,1pg13.5)')
     &   ' localize: iq,ist,mm,fpair=',iq,ist,mm,fpair
        call nablaz_spin(psi(0,0,1,ist),psin,nparz(1,ist),1,one)
        call nablar_spin(psi(0,0,1,ist),psip,   mpsi(ist),1,one)
cPK        do is=1,2
cPK          call mulz(derz1(0,0,nparz(1,ist)),psi(0,0,is,ist),
cPK     &              psin(0,0,is),one,0)
cPK          call mulr(derr1(0,0,nparr),psi(0,0,is,ist),
cPK     &              psip(0,0,is),one,0)
cPK        enddo        
        call densit(psi(0,0,1,ist),vsqf(ist),ist,rho(0,0,iq),
     &              tau(0,0,iq),wrk3(0,0),wrk3(0,0),wrk3(0,0))
        do ir=0,nr
        do iz=0,nz
        do is=1,2
          psiact = psi(iz,ir,is,ist)
ctest          rho(iz,ir,iq) = fpair*(psiact*psiact) + rho(iz,ir,iq)
ctest          tau(iz,ir,iq) = fpair*(psip(iz,ir,is)**2+psin(iz,ir,is)**2)
ctest     &                   + tau(iz,ir,iq)
c          write(78,'(3i5,3(1pg13.5))')
c     &     iz,ir,is,psiact,rho(iz,ir,iq),tau(iz,ir,iq)
          psiact = psiact*fpair
          wrk1(iz,ir) = psip(iz,ir,is)*psiact + wrk1(iz,ir)
          wrk2(iz,ir) = psin(iz,ir,is)*psiact + wrk2(iz,ir)
        enddo
        enddo
        enddo
      enddo
c     
c     interpolate to standard grid points and sum up
c 
      if(ttest) write(6,'(a)')   'localization:',
     &      '  iz  ir    z     r      c      c_tilde   rho'
      write(iwunit,'(a,i3/a)')  
     &      '# for nucleon type iq=',iq,
     &      '#  iz  ir    z     r      c      c_tilde   rho'
      ttffac = 0.6*(3.0*pi**2)**0.6666666667
      do ir=0,nr
      do iz=0,nz
        rhoact = rho(iz,ir,iq)
c          write(78,'(2i5,3(1pg13.5))')
c     &     iz,ir,rho(iz,ir,iq),tau(iz,ir,iq)
        if(rhoact.gt.SMAL) then
          clocal = tau(iz,ir,iq) 
     &            -(wrk1(iz,ir)**2+wrk2(iz,ir)**2)/rhoact
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
c        write(iwunit,'(1x,2i5,2f8.3,3(1pg13.5))')
c     &      iz,ir,zpoint,(ir-0.5D0)*dr,clocal,ctilde,rhoact
         write(iwunit,'(f8.5)') ctilde
      enddo
      write(iwunit,'(5x,a)') ' '
      enddo
      close(78)
c
c     end of outer loop
c
      enddo
c
c
      return
      end

