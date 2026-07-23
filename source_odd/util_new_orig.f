

c=======================================================================
c
c
c     block #5: utilities (printing, math subroutines, etc)
c
c
c=======================================================================
c
c
c======================================================================
c
c     wrt_spprop()
c
c     print sinlge-particle-properties
c
c=======================================================================
      subroutine wrt_spprop(outnam,isinit)
c
c----------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'levels.inc'
c
c----------------------------------------------------------------------
c
      integer           ist
      integer           isinit
      character         outnum*5
      character         outnam*3
c
c----------------------------------------------------------------------
c
c     if initialization should be done (isinit=1)
c----------------------------------------------------------------------
      if (isinit.eq.1) then
         do ist=1,nst(1)+nst(2)
c
c     generate a string containing the index of the file
c----------------------------------------------------------------------
            if (ist.le.nst(1)) then
               WRITE(UNIT=outnum, FMT='(f5.3)') 
     &               1.0D0+ist/1000.0D0
            else
               WRITE(UNIT=outnum, FMT='(f5.3)') 
     &               2.0D0+(ist-nst(1))/1000.D0
            endif
c     
c     open file and write header
c----------------------------------------------------------------------
            open(77,file='levels.'//outnum//'.'//outnam)
            write(77,'(//a)') 
     &           '# beta2, epsilon_i, Var_i, w_i, Delta_i, wtpair_i'
            close(77)
c     
         enddo
c
      else
c
c     standard call (isinit=0)
c----------------------------------------------------------------------
         do ist=1,nst(1)+nst(2)
c
c     generate a string containing the index of the file
c----------------------------------------------------------------------
            if (ist.le.nst(1)) then
               WRITE(UNIT=outnum, FMT='(f5.3)') 
     &               1.0D0+ist/1000.0D0
            else
               WRITE(UNIT=outnum, FMT='(f5.3)') 
     &               2.0D0+(ist-nst(1))/1000.D0
            endif
c     
c     open file and write header
c----------------------------------------------------------------------
            open(77,file='levels.'//outnum//'.'//outnam,access='append')
            write(77,'(f7.4,1(1pg14.7),1(1pg14.5),3(1pg14.5))')
     &           beta2,
     &           spenrg(ist),spevar(ist),
     &           vsqf(ist),deltaf(ist),wtpair(ist)
            close(77)
c     
         enddo
      endif
c
c----------------------------------------------------------------------
c
      return
      end


c-----testpr------------------------------------------------------------
c
      subroutine testpr(field,nzmin,nzmax,nrmin,nrmax,header)
*
***********************************************************************
*                                                                     *
*        testprint of real 'field' in the range nzmin<iz<nzmax and    *
*        nrmin<ir<nrmax. maximally 10 ir are printed.                 *
*        print is named by 'header'.                                  *
*        maximum and minimum value of 'field' is evaluated and        *
*        printed (in whole range).                                    *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
*
      character*(*) header
      dimension field(0:nzd,0:nrd)
cxx      dimension field(0:nrd,0:nrd)
*
      write(7,'(a,a)') ' testprint of field:',header
      do 10 irmin=nrmin,nrmax,10
      nrpmax = min0(nrmax,irmin+9)
      write(7,'(a,3x,a,i2,9(7x,a,i2))')
     x        ' nz',('ir=',ir,ir=irmin,nrpmax)
      do 10 iz=nzmin,nzmax
      write(7,'(i3,10g12.4)') iz,(field(iz,ir),ir=irmin,nrpmax)
   10 continue
*
      fldmin = +1.e30
      fldmax = -1.e30
      do 20 ir=0,nr
      do 20 iz=0,nz
      fldmin = min(abs(field(iz,ir)),fldmin)
      fldmax = max(abs(field(iz,ir)),fldmax)
   20 continue
      write(7,'(/a,g12.4,a,g12.4)') ' abs maximum value = ',fldmax,
     x                              ' abs minimum value = ',fldmin
*
      return
      end
c
c-----f2int------------------------------------------------------------
c
      function f2int(fua,fub)
c
c     computes the integral over fua*fub
c
      include 'paramr.inc'
      include 'dynpar.inc'
c
      dimension fua(0:nzd,0:nrd),fub(0:nzd,0:nrd)
c
c     summing up the wavefunction
c
      drezw  = one/12.0
c
      sum1   = oh*fua(0,0)*fub(0,0)
      do 10 iz=1,nz
       sum1  = fua(iz,0)*fub(iz,0)+sum1
10    continue
      sum1   =  drezw*sum1               ! Euler-Mc'Laurin correction at r=0
c
      sum2   = zero                      ! accumulator for boundary z=0
      sum3   = zero                      ! accumulator for the interior
      dir    = sr
      do 30 ir=1,nr
      dir   = one + dir                  ! the actual value for 'r'
      sum2  = fua(0,ir)*fub(0,ir)*dir+sum2
      do 30 iz=1,nz
       sum3  = fua(iz,ir)*fub(iz,ir)*dir + sum3
30    continue
      sum2   = oh*sum2
c
      sum = sum1+sum2+sum3
c
c     the final result
c
      f2int = sum*(sym+one)*dvol
c
      return
      end
c
c-----f1int------------------------------------------------------------
c
      function f1int(fua)
c
c     computes the integral over fua
c
      include 'paramr.inc'
      include 'dynpar.inc'
c
      dimension fua(0:nzd,0:nrd)
c
      drezw  = one/12.0
c
      sum1   = oh*fua(0,0)
      do 10 iz=1,nz
       sum1  = fua(iz,0)+sum1
10    continue
      sum1   =  drezw*sum1            ! Euler-Mc'Laurin correction at r=0
c
      sum2   = zero                   ! accumulator for boundary z=0
      sum3   = zero                   ! accumulator for all interior
      dir    = sr                     ! set value of r to left boundary 'sr'
      do 30 ir=1,nr
      dir   = one + dir
      sum2  = fua(0,ir)*dir+sum2
      do 30 iz=1,nz
       sum3  = fua(iz,ir)*dir + sum3
30    continue
      sum2   = oh*sum2
c
      sum = sum1+sum2+sum3
c
c     the final result
c
      f1int = sum*(sym+one)*dvol
c
      return
      end

c=======================================================================
c
c     upd_derive_spinor(psia,drpsia,dzpsia,dppsia,ddpsia,zpar,imu,bpar)
c
c     calculates the derivatives of an array of sp-wavefunction.
c
c     imu is an integer global shift in K component which e.g. could be
c     the result of the action of a multipole operator:
c
c         K(psia)    =  K(psi) + imu
c
c     zpar is an shift in z-parity
c
c         Pi_z(psia) =  zpar*Pi_z(psi)
c
c     where psi are the states of the stationary iteration.
c     
c=======================================================================
      subroutine upd_derive_spinor(psia,drpsia,dzpsia,dppsia,ddpsia,
     &                             zpar,imu,bpar)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'  ! one
      include 'levels.inc'  ! nst(1/2)
c
c     ------------------------------------------------------------------
c
      integer          imu,zparist,zpar,bpar
      double precision psia(0:nzd,0:nrd,nsd,npsid)
      double precision drpsia(0:nzd,0:nrd,nsd,npsid)
      double precision dzpsia(0:nzd,0:nrd,nsd,npsid)
      double precision dppsia(0:nzd,0:nrd,nsd,npsid)
      double precision ddpsia(0:nzd,0:nrd,nsd,npsid)
      integer          ist
c
c     calculate derivatives of single particle wavefunctions
c     ------------------------------------------------------------------
c
c      imu=abs(imu)  !! this should be removed
c
      do ist=1,nst(1)+nst(2)
c
         if (zpar.eq.+1) then
            zparist = nparz(1,ist)
         else
            zparist = 3-nparz(1,ist)
         endif

         call nablaz_spin(psia(0,0,1,ist),dzpsia(0,0,1,ist),
     &                    zparist,bpar,one)
         call nablar_spin(psia(0,0,1,ist),drpsia(0,0,1,ist),
     &                    mpsi(ist)+imu,bpar,one) 
         call nablap_spin(psia(0,0,1,ist),dppsia(0,0,1,ist),
     &                    mpsi(ist)+imu,bpar,one) 
         call laplace_spin(psia(0,0,1,ist),ddpsia(0,0,1,ist),
     &                     zparist,mpsi(ist)+imu,bpar,
     &                     one) 
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     damp_spinor(psia,psib,zpar,imu,bpar)
c
c     kinetic energy damping of an array of sp-wavefunction.
c
c     imu is an integer global shift in K component which e.g. could be
c     the result of the action of a multipole operator:
c
c         K(psia)    =  K(psi) + imu
c
c     zpar is an shift in z-parity
c
c         Pi_z(psia) =  zpar*Pi_z(psi)
c
c     where psi are the states of the stationary iteration.
c     
c=======================================================================
      subroutine damp_spinor(psia,psib,zpar,imu,bpar,fac)
c
c     ------------------------------------------------------------------
c
      include 'paramr.inc'  ! one
      include 'levels.inc'  ! nst(1/2)
c
c     ------------------------------------------------------------------
c
      integer          imu,zparist,zpar,bpar
      double precision psia(0:nzd,0:nrd,nsd,npsid)
      double precision psib(0:nzd,0:nrd,nsd,npsid)
      integer          ist
      double precision fac
c
c     calculate derivatives of single particle wavefunctions
c     ------------------------------------------------------------------
      do ist=1,nst(1)+nst(2)
c
         if (zpar.eq.+1) then
            zparist = nparz(1,ist)
         else
            zparist = 3-nparz(1,ist)
         endif
 
         call damp_nucleon(psia(0,0,1,ist),psib(0,0,1,ist),
     &                     zparist,mpsi(ist)+imu,fac)
c
      enddo
c
c     ------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     act_mfh_even(psi,drpsi,dzpsi,dppsi,Hpsi
c                  fldU,actU,
c                  fldB,drfldB,dzfldB,actB,
c                  fldQr,fldWz,actW,
c                  fldQ,lambda,actQ,
c                  diagE,actE)
c     
c     computes the action of a Mean-Field-Hamiltonian on an array of
c     states Psi(1..nst). 
c
c=======================================================================
      subroutine act_mfh_even(
     &        Psi,drPsi,dzPsi,dpPsi,ddPsi,H_Psi,
     &        fldU,                       actU,
     &        fldB,drfldB,dzfldB,dpfldB,  actB,
     &        fldWr,fldWz,fldWp,          actW,
     &        fldQ,lambda,                actQ,
     &        diagE,                      actE,
     &        imu,fac)
c
c------------------------------------------------------------------------
c
      include 'paramr.inc'  ! nzd,nrd,nsd,npsid
      include 'dynpar.inc'  ! nz,nr
      include 'levels.inc'  ! nst(1/2)
c
      double precision Psi(0:nzd,0:nrd,nsd,npsid)
      double precision drPsi(0:nzd,0:nrd,nsd,npsid)
      double precision dzPsi(0:nzd,0:nrd,nsd,npsid)
      double precision dpPsi(0:nzd,0:nrd,nsd,npsid)
      double precision ddPsi(0:nzd,0:nrd,nsd,npsid)
      double precision H_Psi(0:nzd,0:nrd,nsd,npsid)
c
      double precision fldU(0:nzd,0:nrd,kqm)
      double precision fldB(0:nzd,0:nrd,kqm)
      double precision fldWr(0:nzd,0:nrd,kqm)
      double precision fldWz(0:nzd,0:nrd,kqm)
      double precision fldWp(0:nzd,0:nrd,kqm)
      double precision fldQ(0:nzd,0:nrd)

      double precision diagE(npsid)
      double precision drfldB(0:nzd,0:nrd,kqm)
      double precision dzfldB(0:nzd,0:nrd,kqm)
      double precision dpfldB(0:nzd,0:nrd,kqm)
c
      double precision lambda,fac
c
      integer          actU,actB,actW,actQ,actE
      integer          iq,ist,ir,iz,imu
c
c------------------------------------------------------------------------
c
c     loop all states and define isospin 
c     (first all components which also contribute in axial symmetry)
c     ===================================================================
      do ist=1,nst(1)+nst(2)
c         
c     define isospin
c
         if (ist.le.nst(1)) then
            iq=1
         else              
            iq=2
         endif
c
c     if accumulate switch is off --> zero HPsi
c
         do ir=0,nr
            do iz=0,nz
               H_Psi(iz,ir,1,ist) = zero
               H_Psi(iz,ir,2,ist) = zero
            enddo
         enddo
c      
c     calculate action of hamiltonian on state psi...
c
c
c     the local potential
c
         if (actU.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) = H_Psi(iz,ir,1,ist) 
     &                              + fldU(iz,ir,iq) * Psi(iz,ir,1,ist)   
                  H_Psi(iz,ir,2,ist) = H_Psi(iz,ir,2,ist) 
     &                              + fldU(iz,ir,iq) * Psi(iz,ir,2,ist)
               enddo
            enddo
         endif
c
c     the effective mass
c
         if (actB.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) = H_Psi(iz,ir,1,ist) 
     &                          - fldB(iz,ir,iq)   * ddPsi(iz,ir,1,ist) 
     &                          - drfldB(iz,ir,iq) * drPsi(iz,ir,1,ist) 
     &                          - dzfldB(iz,ir,iq) * dzPsi(iz,ir,1,ist) 
                  H_Psi(iz,ir,2,ist) = H_Psi(iz,ir,2,ist) 
     &                          - fldB(iz,ir,iq)   * ddPsi(iz,ir,2,ist) 
     &                          - drfldB(iz,ir,iq) * drPsi(iz,ir,2,ist) 
     &                          - dzfldB(iz,ir,iq) * dzPsi(iz,ir,2,ist) 
               enddo
            enddo
         endif
c
c     the spin-orbit potential
c
         if (actW.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) = H_Psi(iz,ir,1,ist) 
     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist) 
     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist) 
     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist) 
                  H_Psi(iz,ir,2,ist) = H_Psi(iz,ir,2,ist) 
     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
               enddo
            enddo
         endif
c
c     the deformation constraint 
c     
         if (actQ.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) =   H_Psi(iz,ir,1,ist) 
     &                          - lambda*fldQ(iz,ir) * Psi(iz,ir,1,ist)
                  H_Psi(iz,ir,2,ist) =   H_Psi(iz,ir,2,ist) 
     &                          - lambda*fldQ(iz,ir) * Psi(iz,ir,2,ist)
               enddo
            enddo
         endif
c
c     the norm constraint (energy, diagonal)
c     
         if (actE.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) =   H_Psi(iz,ir,1,ist) 
     &                              - diagE(ist) * Psi(iz,ir,1,ist)
                  H_Psi(iz,ir,2,ist) =   H_Psi(iz,ir,2,ist) 
     &                              - diagE(ist) * Psi(iz,ir,2,ist)
               enddo
            enddo
         endif
c
         if (actE.eq.-1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) =   H_Psi(iz,ir,1,ist) 
     &                              + diagE(ist) * Psi(iz,ir,1,ist)
                  H_Psi(iz,ir,2,ist) =   H_Psi(iz,ir,2,ist) 
     &                              + diagE(ist) * Psi(iz,ir,2,ist)
               enddo
            enddo
         endif
c
c     end of loop over all states
c     
      enddo
c
c     loop all states and define isospin 
c     (now all components which contribute in case of NO AXIAL symmetry)
c     ===================================================================
      if (imu.ne.0) then
c
      do ist=1,nst(1)+nst(2)
c         
c     define isospin
c
         if (ist.le.nst(1)) then
            iq=1
         else              
            iq=2
         endif
c
c     the effective mass
c
         if (actB.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) = H_Psi(iz,ir,1,ist) 
     &                          - dpfldB(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
                  H_Psi(iz,ir,2,ist) = H_Psi(iz,ir,2,ist) 
     &                          - dpfldB(iz,ir,iq) * dpPsi(iz,ir,2,ist) 
               enddo
            enddo
         endif
c
c     the spin-orbit potential
c 
         if (actW.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  H_Psi(iz,ir,1,ist) = H_Psi(iz,ir,1,ist) 
     &                          -  fldWp(iz,ir,iq) * drPsi(iz,ir,1,ist) 
     &                          +  fldWp(iz,ir,iq) * dzPsi(iz,ir,2,ist) 
                  H_Psi(iz,ir,2,ist) = H_Psi(iz,ir,2,ist) 
     &                          +  fldWp(iz,ir,iq) * drPsi(iz,ir,2,ist)
     &                          +  fldWp(iz,ir,iq) * dzPsi(iz,ir,1,ist) 
                  
               enddo
            enddo
         endif
c
c     end of loop over all states
c     
      enddo
c
      endif
c
c     weight
c
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               H_Psi(iz,ir,1,ist) = fac*H_Psi(iz,ir,1,ist)
               H_Psi(iz,ir,2,ist) = fac*H_Psi(iz,ir,2,ist)
            enddo
         enddo
      enddo
c
c------------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     act_mfh_odd(psi,drpsi,dzpsi,dppsi,Hpsi
c                  fldU,actU,
c                  fldB,drfldB,dzfldB,actB,
c                  fldQr,fldWz,actW,
c                  fldQ,lambda,actQ,
c                  diagE,actE)
c     
c     computes the action of a Mean-Field-Hamiltonian on an array of
c     states Psi(1..nst). 
c
c=======================================================================
      subroutine act_mfh_odd(
     &        Psi,drPsi,dzPsi,dpPsi,ddPsi,HPsi,
     &        fldAr,fldAz,fldAp,flddivA,   actA,
     &        fldSr,fldSz,fldSp,           actS,
     &        fldCr,fldCz,fldCp,           actC,
     &        imu,fac)
c
c------------------------------------------------------------------------
c
      include 'paramr.inc'  ! nzd,nrd,nsd,npsid
      include 'dynpar.inc'  ! nz,nr
      include 'levels.inc'  ! nst(1/2)
c
      double precision Psi(0:nzd,0:nrd,nsd,npsid)
      double precision drPsi(0:nzd,0:nrd,nsd,npsid)
      double precision dzPsi(0:nzd,0:nrd,nsd,npsid)
      double precision dpPsi(0:nzd,0:nrd,nsd,npsid)
      double precision ddPsi(0:nzd,0:nrd,nsd,npsid)
      double precision HPsi(0:nzd,0:nrd,nsd,npsid)
c
      double precision fldAr(0:nzd,0:nrd,kqm)
      double precision fldAz(0:nzd,0:nrd,kqm)
      double precision fldAp(0:nzd,0:nrd,kqm)
      double precision fldSr(0:nzd,0:nrd,kqm)
      double precision fldSz(0:nzd,0:nrd,kqm)
      double precision fldSp(0:nzd,0:nrd,kqm)
      double precision fldCr(0:nzd,0:nrd,kqm)
      double precision fldCz(0:nzd,0:nrd,kqm)
      double precision fldCp(0:nzd,0:nrd,kqm)

      double precision divfldA(0:nzd,0:nrd,kqm)
c
      integer          actA,actS,actC
      integer          iq,ist,ir,iz
c
      double precision fac
c
c------------------------------------------------------------------------
c
c     loop all states and define isospin 
c     (first all components which also contribute in axial symmetry)
c     ===================================================================
      do ist=1,nst(1)+nst(2)
c
c     define isospin
c       
         if (ist.le.nst(1)) then
            iq=1
         else              
            iq=2
         endif
c
c     if accumulate switch is off --> zero HPsi
c
         do ir=0,nr
            do iz=0,nz
               HPsi(iz,ir,1,ist) = zero
               HPsi(iz,ir,2,ist) = zero
            enddo
         enddo
c      
c     calculate action of time-odd hamiltonian on state psi...
c
c
c     the dynamic potential
c     
         if (actA.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
     &                 + 0.5D0*divfldA(iz,ir,iq) *   Psi(iz,ir,1,ist)
     &                 +         fldAr(iz,ir,iq) * drPsi(iz,ir,1,ist)
     &                 +         fldAz(iz,ir,iq) * dzPsi(iz,ir,1,ist)
!GH
               enddo
               do iz=0,nz
                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
     &                 + 0.5D0*divfldA(iz,ir,iq) *   Psi(iz,ir,2,ist)
     &                 +         fldAr(iz,ir,iq) * drPsi(iz,ir,2,ist)
     &                 +         fldAz(iz,ir,iq) * dzPsi(iz,ir,2,ist)
               enddo
            enddo
         endif
c     
c     the spin-potential
c
         if (actS.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
     &                          - fldSp(iz,ir,iq)  * Psi(iz,ir,2,ist) 
                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
     &                          + fldSp(iz,ir,iq)  * Psi(iz,ir,1,ist) 
               enddo
            enddo
         endif
c
c     the effective spin-mass
c
c         if (actC.eq.1) then
c            do ir=0,nr
c               do iz=0,nz
c                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
c     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist) 
c     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist) 
c     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
c     &                          +  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist) 
c                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
c     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
c     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
c     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
c     &                          +  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
c               enddo
c            enddo
c         endif
c
c     end of loop over all states
c     
      enddo
c
c     loop all states and define isospin 
c     (now all components which contribute in case of NO AXIAL symmetry)
c     ===================================================================
      if (imu.gt.0) then
c
      do ist=1,nst(1)+nst(2)
c
c     define isospin
c       
         if (ist.le.nst(1)) then
            iq=1
         else              
            iq=2
         endif
c
c     the dynamic potential
c
         if (actA.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
     &                    +         fldAp(iz,ir,iq) * dpPsi(iz,ir,1,ist)
                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
     &                    +         fldAp(iz,ir,iq) * dpPsi(iz,ir,2,ist)
               enddo
            enddo
         endif
c
c     the spin-potential
c
         if (actS.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
     &                          + fldSr(iz,ir,iq)  * Psi(iz,ir,2,ist) 
     &                          + fldSz(iz,ir,iq)  * Psi(iz,ir,1,ist) 
                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
     &                          + fldSr(iz,ir,iq)  * Psi(iz,ir,1,ist) 
     &                          - fldSz(iz,ir,iq)  * Psi(iz,ir,2,ist) 
               enddo
            enddo
         endif
c
c     the effective spin-mass
c
c         if (actC.eq.1) then
c            do ir=0,nr
c               do iz=0,nz
c                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
c     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist) 
c     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist) 
c     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
c     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist) 
c                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
c     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
c     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
c     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
c     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
c               enddo
c            enddo
c         endif
c
c     end of loop over all states
c     
      enddo
c
      else if (imu.lt.0) then
c
      do ist=1,nst(1)+nst(2)
c
c     define isospin
c       
         if (ist.le.nst(1)) then
            iq=1
         else              
            iq=2
         endif
c
c     the dynamic potential
c
         if (actA.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
     &                    -         fldAp(iz,ir,iq) * dpPsi(iz,ir,1,ist)
                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
     &                    -         fldAp(iz,ir,iq) * dpPsi(iz,ir,2,ist)
               enddo
            enddo
         endif
c
c     the spin-potential
c
         if (actS.eq.1) then
            do ir=0,nr
               do iz=0,nz
                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
     &                          - fldSr(iz,ir,iq)  * Psi(iz,ir,2,ist) 
     &                          - fldSz(iz,ir,iq)  * Psi(iz,ir,1,ist) 
                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
     &                          - fldSr(iz,ir,iq)  * Psi(iz,ir,1,ist) 
     &                          + fldSz(iz,ir,iq)  * Psi(iz,ir,2,ist) 
               enddo
            enddo
         endif
c
c     the effective spin-mass
c
c         if (actC.eq.1) then
c            do ir=0,nr
c               do iz=0,nz
c                  HPsi(iz,ir,1,ist) = HPsi(iz,ir,1,ist) 
c     &                          -  fldWz(iz,ir,iq) * drPsi(iz,ir,2,ist) 
c     &                          +  fldWr(iz,ir,iq) * dzPsi(iz,ir,2,ist) 
c     &                          +  fldWr(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
c     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,2,ist) 
c                  HPsi(iz,ir,2,ist) = HPsi(iz,ir,2,ist) 
c     &                          +  fldWz(iz,ir,iq) * drPsi(iz,ir,1,ist)
c     &                          -  fldWr(iz,ir,iq) * dzPsi(iz,ir,1,ist)
c     &                          -  fldWr(iz,ir,iq) * dpPsi(iz,ir,2,ist)
c     &                          -  fldWz(iz,ir,iq) * dpPsi(iz,ir,1,ist) 
c               enddo
c            enddo
c         endif
c
c     end of loop over all states
c     
      enddo
c
      endif
c
c     weight
c
      do ist=1,nst(1)+nst(2)
         do ir=0,nr
            do iz=0,nz
               HPsi(iz,ir,1,ist) = fac*HPsi(iz,ir,1,ist)
               HPsi(iz,ir,2,ist) = fac*HPsi(iz,ir,2,ist)
            enddo
         enddo
      enddo
c
c------------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     wfovlp(psia,psib)
c
c     computes the overlap between psia and psib , including spin
c
c=======================================================================
      function wfovlp(psia,psib)
c
c------------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'levels.inc'
c
      dimension psia(0:nzd,0:nrd,nsd),psib(0:nzd,0:nrd,nsd)
c
      parameter         (drezw=1.0D0/12.0D0)
      double precision  wfovlp
c
c------------------------------------------------------------------------
c
c     summing up the wavefunction
c
      sum1   = oh*(psia(0,0,1)*psib(0,0,1)+psia(0,0,2)*psib(0,0,2))
      do iz=1,nz
       sum1 = sum1 + psia(iz,0,1)*psib(iz,0,1)
     &             + psia(iz,0,2)*psib(iz,0,2)
      enddo
      sum1 = sum1*drezw              ! Euler-Mc'Laurin correction at r=0
c
      sum2 = zero                      
      r    = sr
      do ir=1,nr
        sum3 = oh*(psia(0,ir,1)*psib(0,ir,1)
     &             +psia(0,ir,2)*psib(0,ir,2))
        do iz=1,nz
          sum3 = sum3 + psia(iz,ir,1)*psib(iz,ir,1)
     &                + psia(iz,ir,2)*psib(iz,ir,2)
        enddo
        r    = one + r
        sum2 = sum2 + sum3*r
      enddo
c
c     the final result
c
      wfovlp = (sum1+sum2)*(sym+one)*dvol
c
c------------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     wfovlptr(psia,psib)
c
c     computes the overlap between time-reversed state of psia and 
c     psib , including spin
c
c=======================================================================
      function wfovlptr(psia,psib)
c
c------------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'levels.inc'
c
      dimension psia(0:nzd,0:nrd,nsd),psib(0:nzd,0:nrd,nsd)
c
      parameter (drezw=1.0D0/12.0D0)
      double precision  wfovlptr
c
c------------------------------------------------------------------------
c
c     summing up the wavefunction
c
      sum1   = oh*(-psia(0,0,2)*psib(0,0,1)+psia(0,0,1)*psib(0,0,2))
      do iz=1,nz
       sum1 = sum1 - psia(iz,0,2)*psib(iz,0,1)
     &             + psia(iz,0,1)*psib(iz,0,2)
      enddo
      sum1 = sum1*drezw             ! Euler-Mc'Laurin correction at r=0
c
      sum2 = zero                      
      r    = sr
      do ir=1,nr
        sum3 = oh*(-psia(0,ir,2)*psib(0,ir,1)
     &             +psia(0,ir,1)*psib(0,ir,2))
        do iz=1,nz
          sum3 = sum3 - psia(iz,ir,2)*psib(iz,ir,1)
     &                + psia(iz,ir,1)*psib(iz,ir,2)
        enddo
        r    = one + r
        sum2 = sum2 + sum3*r
      enddo
c
c     the final result
c
      wfovlptr = (sum1+sum2)*(sym+one)*dvol
c
c------------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     wfproj(psia,psib)
c
c     projects components of spinor psia out of psib: 
c     |b> = |b> - fac*|a>
c
c=======================================================================
      subroutine wfproj(psia,psib,fac)
c
c------------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
c
      dimension psia(0:nzd,0:nrd,nsd),psib(0:nzd,0:nrd,nsd)
c
      double precision   fac
      integer            ir,iz
c
c------------------------------------------------------------------------
c
c     project out fac*|a> from |b>
c     
      do ir=0,nr
         do iz=0,nz
            psib(iz,ir,1) = psib(iz,ir,1)-fac*psia(iz,ir,1)
            psib(iz,ir,2) = psib(iz,ir,2)-fac*psia(iz,ir,2)
         enddo
      enddo
c
c------------------------------------------------------------------------
c
      return
      end


c=======================================================================
c
c     wfprojtr(psia,psib)
c
c     projects components of time reversed state of psia out of psib: 
c     |b> = |b> - fac*|Ta>
c
c=======================================================================
      subroutine wfprojtr(psia,psib,fac)
c
c------------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'dynpar.inc'
c
      dimension psia(0:nzd,0:nrd,nsd),psib(0:nzd,0:nrd,nsd)
c
      double precision   fac
      integer            ir,iz
c
c------------------------------------------------------------------------
c
c     project out fac*|Ta> from |b>
c     
      do ir=0,nr
         do iz=0,nz
            psib(iz,ir,1) = psib(iz,ir,1)+fac*psia(iz,ir,2)
            psib(iz,ir,2) = psib(iz,ir,2)-fac*psia(iz,ir,1)
         enddo
      enddo
c
c------------------------------------------------------------------------
c
      return
      end

c
c-----wfovpc-----------------------------------------------------------
c
      function wfovpc(psia)
c
c     computes the integral over psia*psia , no spin summation !
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'levels.inc'
c
      dimension psia(0:nzd,0:nrd)
c
      drezw  = one/12.0
      sum1   = oh*psia(0,0)*psia(0,0)
      do 10 iz=1,nz
       sum1 = sum1 + psia(iz,0)*psia(iz,0)
 10   continue
      sum1 = sum1*drezw                  ! Euler-Mc'Laurin correction at r=0
c
      sum2 = zero                        ! accumulator for boundary z=0
      sum3 = zero                        ! accumulator for the interior
      r    = sr
      do 30 ir=1,nr
      r    = one+ r
      sum2 = sum2 + psia(0,ir)*psia(0,ir)*r
      do 30 iz=1,nz
       sum3 = sum3 + psia(iz,ir)*psia(iz,ir)*r
 30   continue
      sum2 = sum2*oh
c
      sum = sum1+sum2+sum3
c
c     the final result
c
      wfovpc = sum*(sym+one)*dvol
c
      return
      end
c
c-----wfovoc-----------------------------------------------------------
c
      function wfovoc(psia,psib)
c
c     computes the integral over psia*psia , no spin summation !
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'levels.inc'
c
      dimension psia(0:nzd,0:nrd),psib(0:nzd,0:nrd)
c
      drezw  = one/12.0
      sum1   = oh*psia(0,0)*psib(0,0)
      do 10 iz=1,nz
       sum1 = sum1 + psia(iz,0)*psib(iz,0)
 10   continue
      sum1 = sum1*drezw                  ! Euler-Mc'Laurin correction at r=0
c
      sum2 = zero                        ! accumulator for boundary z=0
      sum3 = zero                        ! accumulator for the interior
      r    = sr
      do 30 ir=1,nr
      r    = one+ r
      sum2 = sum2 + psia(0,ir)*psib(0,ir)*r
      do 30 iz=1,nz
       sum3 = sum3 + psia(iz,ir)*psib(iz,ir)*r
 30   continue
      sum2 = sum2*oh
c
      sum = sum1+sum2+sum3
c
c     the final result
c
      wfovoc = sum*(sym+one)*dvol
c
      return
      end
c
c

c-----div   -----------------------------------------------------------
c
      subroutine div(divout,fldrin,fldzin)
c
c     divergence of vector field:
c        fldrin  =  r-component of vector field, input
c        fldzin  =  z-component of vector field, input
c        divout  =  divergence, output
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'wavefr.inc'
*
      dimension divout(0:nzd,0:nrd)
      dimension fldrin(0:nzd,0:nrd)
      dimension fldzin(0:nzd,0:nrd)
*
*-----------------------------------------------------------
*
*
*---- calculate derivatives of fields
*
      call nablaz_scalar(fldzin(0,0),wrk1(0,0), 1,1,one)
      call nablar_scalar(fldrin(0,0),wrk2(0,0), 0,1,one)
*
*---- value at 'r=0'
*
      do iz=0,nz
          divout(iz,0) = 2.0*wrk2(iz,0)+wrk1(iz,0) 
      enddo
*
*---- all other 'r'
*
      do ir=1,nr
       r      = (sr+ir)*dr
       rinv = one/r
       do iz=0,nz
        divout(iz,ir) = wrk2(iz,ir)+rinv*fldrin(iz,ir)+wrk1(iz,ir) 
       enddo
      enddo
c
      return
      end

c-----testls-----------------------------------------------------------
c
      subroutine testls
c
c     test of l*s current by checking vector field against its divergence
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'wavefr.inc'
*
*
*-----------------------------------------------------------
*
      call div(wrk3,rholsr(0,0,1),rholsz(0,0,1))
      write(6,'(/a/a)')
     & ' test spin-orbit fields for protons:',
     & '  ir   iz   div-direct    div-deduced '
      do ir=1,nr
       do iz=0,nz
         write(6,'(2i5,3(1pg13.4))')
     &   ir,iz,drhols(iz,ir,1),wrk3(iz,ir),drhols(iz,ir,1)-wrk3(iz,ir)
       enddo
      enddo

      call div(wrk3,rholsr(0,0,2),rholsz(0,0,2))
      write(6,'(/a/a)')
     & ' test spin-orbit fields for neutrons:',
     & '  ir   iz   div-direct    div-deduced '
      do ir=1,nr
       do iz=0,nz
         write(6,'(2i5,3(1pg13.4))')
     &   ir,iz,drhols(iz,ir,2),wrk3(iz,ir),drhols(iz,ir,2)-wrk3(iz,ir)
       enddo
      enddo
c
      return
      end

c-----pripot-----------------------------------------------------------
c
      subroutine pripot(title,outnam)
c
c     print various potentials on file 'potentials.<name>'
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'fields.inc'
      include 'wavefr.inc'
      include 'levels.inc'
c
      character*3  outnam
      character*80 title
c
c---------------------------------------------------------------------------
c
      open(31,status='unknown',form='formatted',
     &     file='potentials.'//outnam)
      do ispool=1,69999
       read(31,'(1x)',end=9999)
      enddo
      stop " file count exhausted for print of potentials "
 9999 continue
c
      write(31,'(a/2a/a,2i5/a,2f6.1,2i5/a,3f7.3)') 
     &     ' potentials:',' title=',title,
     &     ' A,Z=',npart(1)+npart(2),npart(1),
     &     ' dz,dr,nz,nr=',dz,dr,nz,nr,
     &     ' beta2,beta2p,beta2n=',beta2,beta2p,beta2n
      write(31,'(2a)')
     &  '  iz ir    z     r     V_proton      V_neutron   ',
     &  '    V_Coulomb   V_prot-V_coul '
      do ir=1,nr
       do iz=0,nz
         r  = (ir+sr)*dr
         z  = (iz+sz-icm)*dz
         write(31,'(1x,2i3,2f6.1,4g15.7)')
     &   iz,ir,z,r,pot(iz,ir,1),pot(iz,ir,2),coulpt(iz,ir),
     &   pot(iz,ir,1)-coulpt(iz,ir)
       enddo
      enddo
c
      close(31)
      write(6,'(a)') ' potentials written on file.'
c
      return
      end
ctestc-----yuktest-----------------------------------------------------------
ctestc
ctest      subroutine yuktest
ctestc
ctestc**********************************************************************
ctestc                                                                     *
ctestc     tests the yukawa solver                                         *
ctestc                                                                     *
ctestc**********************************************************************
ctestc
ctest      include 'paramr.inc'
ctest      include 'fields.inc'
ctest      include 'dynpar.inc'
ctestc
ctest      dimension rhs(0:nzd,0:nrd)
ctestc
ctestc
ctest      a      = 4.0
ctest      b      = sqrt(h2meff*a**(one/3)/41.0)
ctestc
ctest      do 5 ir=0,nr
ctest      rsqr          = ((ir)+sr)**2 * drsqr
ctest      do 5 iz=0,nz
ctest      zsqr          = ((iz)+sz)**2 * dzsqr
ctest      rhosqr        = rsqr + zsqr
ctestcc    exarg1        = (rho-rad)/sur
ctestcc    exarg2        =-(rho+rad)/sur
ctestc
ctestcc    if (earg1.ge.-eighty) then
ctestcc      exfac1      = exp(exarg1)
ctestcc    else
ctestcc      exfac1      = 0.0
ctestcc    endif
ctestcc    if (exarg2.ge.-eighty) then
ctestcc      exfac2      = exp(exarg2)
ctestcc    else
ctestcc      exfac1      = 0.0
ctestcc    endif
ctestc
ctestcc    opex1         = (one+exfac1)
ctestcc    opex2         = (one+exfac2)
ctestcc    fac1          =-exfac1/(sur*opex1**2*opex2)
ctestcc    fac2          = exfac2/(sur*opex2**2*opex1)
ctestcc    fac3          = 2.0*exfac1**2/(sur*opex1**3*opex2)
ctestcc    fac4          =-2.0*erd/(sur*opex1**2*opex2**2)
ctestcc    fac5          = 2.0*exfac2**2/(sur*opex1*opex2**3)
ctestc
ctestcc    td2           = (fac1-fac2+fac3+fac4+fac5)/sur
ctestcc    d1            =  fac1+fac2
ctestcc    td1           = 2.0*d1/rho
ctestcc    f             = one/(opex1*opex2)
ctestcc    tmass         = xmass**2*f
ctestc
ctestcc    rhs(iz,ir)    = (-td2-td1+tmass)*v0
ctestcc    smeson(iz,ir) = f*v0
ctest      pfac    =  4.0*vcoupl/(pi**1.5*b**3)
ctest      earg    =  (-rhosqr/b**2)
ctest      if (abs(earg).le.50.0*one) then
ctest         rhs(iz,ir) = pfac*exp(earg)
ctest      else
ctest         rhs(iz,ir) = 0.0
ctest      endif
ctest 5    continue
ctest      call mesini(vmeson,rhs,vmshc)
ctest      call mesrun(vmeson,rhs,vmshc,2,1)
ctestc     call mesini(vmeson,rhs,vaccel,vmass)
ctestc     call mesrun(vmeson,rhs,vaccel,vmass)
ctestc
ctestcc    errmx = -one
ctestcc    do 6 ir=0,nr
ctestcc    do 6 iz=0,nz
ctestcc     err= abs(vmeson(iz,ir)-smeson(iz,ir))
ctestcc     if (err.gt.errmx) then
ctestcc       errmx = err
ctestcc       irmx  = ir
ctestcc       izmx  = iz
ctestcc     endif
ctestc6    continue
ctestc     tsum = 0.0
ctestc     do 6 ir=0,nr
ctestc     r  = (sr+ir)
ctestc     do 6 iz=0,nz
ctestc       tsum = tsum+vmeson(iz,ir)*r
ctestc6    continue
ctestc     tsum = tsum*dz*dr
ctestc     write(7,'(a,g15.5)') ' testsum=',tsum
ctest      call testpr(rhs,0,9,0,9,' density')
ctestcc    write(7,'(/a,g15.5,2(a,i3)/)') ' error(ana-num)=',errmx,' ir=',
ctestcc   &                             irmx,' iz=',izmx
ctestc
ctestc      call testpr(vmeson,0,9,0,9,' vmeson')
ctestc     call testpr(smeson,0,9,0,9,' ana result')
ctestc
ctestc     call testpr(vmeson,17,25,26,34,' num result')
ctestc     call testpr(smeson,17,25,26,34,' ana result')
ctestc
ctestc     call testpr(vmeson,0,9,11,20,' num result')
ctestc     call testpr(smeson,0,9,11,20,' ana result')
ctestc
ctestc     call testpr(vmeson,11,20,11,20,' num result')
ctestc     call testpr(smeson,11,20,11,20,' ana result')
ctestc
ctest      stop
ctest      end
