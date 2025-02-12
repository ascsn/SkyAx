c
c=======================================================================
c
c     block #4: the main routines for the propagation
c
c=======================================================================
c-----hpsi_old  -----------------------------------------------------------
c
      subroutine hpsi(psiin,psiout,ist)
*
***********************************************************************
*                                                                     *
*        action of the mean field hamiltonian on the wave function    *
*        'psiin'. result is returned on 'psiout'.                     *
*        at the same time hpsi evaluates the single particle energy   *
*        and the variances of the state 'ist' and returns it on the   *
*        common /levels/.                                             *
*                                                                     *
*        'psiin' is supposed to be related to the state 'ist' in      *
*        the global ordering.                                         *
*                                                                     *
*        notes:                                                       *
*                                                                     *
*          ist    =  absolute adress of the state, to be used in      *
*                    all fields running in the range (1:npsid)        *
*                                                                     *
*          calculation of the variance is done by <Hist|Hist> here.   *
*          this turned out to be bad for states with M=0 and M=1      *
*          due to the non hermicity of the hamiltonian                *
*          if you need a improved calculation of the variance you     *
*          should use the subroutine hhpsi insted.                    *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'observ.inc'
*
      dimension psiin(0:nzd,0:nrd,nsd),psiout(0:nzd,0:nrd,nsd)
      dimension diffpsi(0:nzd,0:nrd,nsd)
      dimension psic(0:nzd,0:nrd,nsd)
*
      double precision xifac

      data ifclean/1/
*
***********************************************************************
*     
      IF(b_strength.NE.0D0) b_str_act(:)=SQRT(e2)*h2m(:)/hc*b_strength
!      IF(b_strength.NE.0D0) WRITE(*,*) 'b_str_act=',b_str_act(:)

c      if(((econ.lt.(0.1D0*errenr)).and.(devh2old.lt.2.0D0*errvarr))
c     & .or.(econ.le.1E-9 .and. devh2old.lt.10.0D0*errvar)
c     & .or.(itlam.gt.(itmax-5)))then
c         xifac = qlamb*0.0D0
c         write(7,*)itlam 
c      else
c         xifac = 1/(qw*qw)*qlamb*0.0D0
c      endif    
	
      iq      = icharg(ist)
      ebin    = spenrg(ist)
c
c     first, zero the output wavfunction because 'psiout' is
c     accumulated during the next spin loop.
c     ?? note: it may be more efficient to introduce workspaces
c              which distinguish spin 'is'. This would allow
c              to avoid costly zeroing and to avoid the spin loop.
c
      do is=1,2
      do ir=0,nr
      do iz=0,nz
        psiout(iz,ir,is) = zero
      enddo
      enddo
      enddo
c
      call laplace_spin(psiin,psi1,nparz(1,ist),mpsi(ist),1,one)
      call nablar_spin(psiin,psi2,mpsi(ist),1,one)
      call nablaz_spin(psiin,psi3,nparz(1,ist),1,one)
c
      do 50 is=1,2
       isi       = 3-is
       mm        = mpsi(ist)+is-1
       nparr     = 1+mod(mm,2)
       sigis     = mod(one*(is),two)*two - one
c
c       if (mm.ne.0) then
c          irmin=1
c       else      
          irmin=0
c       endif
cc
c       if (nparz(is,ist).eq.2.or.isym.eq.0) then
c          izmin=1
c       else
          izmin=0
c       endif
c
       if(ifstab.eq.1) then
         eferm2 = fer(iq)*oh
         if(ebin.gt.zero) then
           spefac = one
         elseif(ebin.gt.eferm2) then
           spefac = one-ebin/eferm2
         else
           spefac = zero
         endif
       else
         spefac = zero
       endif
       if(spefac.le.zero) then
         do ir=irmin,nr
         do iz=izmin,nz
            
c       diffpsi(iz,ir,is) = psiin(iz,ir,is) - psiold(iz,ir,is,ist)   


           psiout(iz,ir,is) =
     &      (pot(iz,ir,iq)-qlamb*qfac(iz,ir)-ebin)*psiin(iz,ir,is)
c     &     -xifac*diffpsi(iz,ir,is)
     &     - bmass(iz,ir,iq)*psi1(iz,ir,is)
     &     - drbmas(iz,ir,iq)*psi2(iz,ir,is)
     &     - dzbmas(iz,ir,iq)*psi3(iz,ir,is)
     &     + psiout(iz,ir,is)
         enddo
         enddo
       else
c         write(6,'(a,i5,3g12.5)')
c     &   ' ist,ebin,eferm2,spefac=',ist,ebin,eferm2,spefac

         do ir=irmin,nr
         do iz=izmin,nz
            
c       diffpsi(iz,ir,is) = psiin(iz,ir,is) - psiold(iz,ir,is,ist)   

           psiout(iz,ir,is) =
     &      (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)
     &       -qlamb*qfac(iz,ir)-ebin)*psiin(iz,ir,is)
c     &     -xifac*diffpsi(iz,ir,is)
     &     - bmass(iz,ir,iq)*psi1(iz,ir,is)
     &     - drbmas(iz,ir,iq)*psi2(iz,ir,is)
     &     - dzbmas(iz,ir,iq)*psi3(iz,ir,is)
     &     + psiout(iz,ir,is)
         enddo
         enddo
       endif

!       optionally for additional quadr. constr.
       IF(qlambquad.NE.0D0)   
     &    psiout(:,:,is)= psiout(:,:,is)+qfacquad(:,:)*psiin(:,:,is)

c       optionally magnetic field
       IF(b_strength.NE.0D0) THEN
!        WRITE(*,*) 'add B field:',b_str_act(iq)
        IF(g_magnetic(iq)>0D0) THEN
          facm=1D0
        ELSE
          facm=0D0
        END IF
        do ir=irmin,nr
        do iz=izmin,nz
         psiout(iz,ir,is) =
     &     -(facm*mm+sigis*0.5D0*g_magnetic(iq))*b_str_act(iq)*
     &     psiin(iz,ir,is)   + psiout(iz,ir,is)
!old         psiout(iz,ir,is) =
!old     &     -((2-iq)*mm+sigis*0.5D0*g_magnetic(iq))*b_str_act(iq)*
!old     &     psiin(iz,ir,is)   + psiout(iz,ir,is)
        enddo
        enddo
       END IF
c
       if(t4.ne.zero .OR. ioptsk==13) then        
c
c     the spin-orbit part. action of the s.o. force is computed
c     simply from
c       W_ls |psi> ~ (nabla rho)*sigma sigma*nabla |psi>
c                    - (nabla rho) * nabla |psi>
c
c                               r=0 requires l'Hospital for 1/r*psiin
       if(mm.eq.1) then
        do iz=izmin,nz
          psiout(iz,0,is) =
     &      +sigis*potls(iz,0,iq,1)*psi2(iz,0,is)     ! 1 = r   2 = z
     &      + psiout(iz,0,is)
          psiout(iz,0,isi)=
     &      -sigis*potls(iz,0,iq,1)*psi3(iz,0,is)
     &      -potls(iz,0,iq,2)*(-sigis*psi2(iz,0,is)
     &                 +psi2(iz,0,is) )
     &      + psiout(iz,0,isi)
        enddo
       endif
codd
       if(mm.eq.(-1)) then
        do iz=izmin,nz
          psiout(iz,0,is) =
     &      -sigis*potls(iz,0,iq,1)*psi2(iz,0,is)     ! 1 = r   2 = z
     &      + psiout(iz,0,is)
          psiout(iz,0,isi)=
     &      -sigis*potls(iz,0,iq,1)*psi3(iz,0,is)
     &      +potls(iz,0,iq,2)*(sigis*psi2(iz,0,is)
     &                 +psi2(iz,0,is) )
     &      + psiout(iz,0,isi)
        enddo
       endif
codd
c                               now all other r

       do ir=1,nr
         r   = (sr+ir)*dr
         rinv = one/r
         xmm  = mm*rinv
         do iz=izmin,nz
           psiout(iz,ir,is) =
     &      +sigis*xmm*potls(iz,ir,iq,1)*psiin(iz,ir,is)
     &      + psiout(iz,ir,is)
           psiout(iz,ir,isi)=
     &      -sigis*potls(iz,ir,iq,1)*psi3(iz,ir,is)
     &      -potls(iz,ir,iq,2)*(-sigis*psi2(iz,ir,is)
     &                 +xmm*psiin(iz,ir,is) )
     &      + psiout(iz,ir,isi)
         enddo
       enddo
      endif


   50 continue
      if((isTSym.ge.2)) then   ! time-odd part of Hamiltonian only if required
         call hpsiodd(psiin,psic,ist)
         do is=1,nsd
          do ir=0,nr
          do iz=0,nz
             psiout(iz,ir,is) = psiout(iz,ir,is)+psic(iz,ir,is)
          enddo
          enddo
         enddo
       endif 

! c
! c     compute expectation value and variance
! c
      esub        = wfovlp(psiin,psiout)
      spenrg2     = wfovlp(psiout,psiout)
      spevar(ist) = sqrt(abs(spenrg2-esub*esub))
      spenrg(ist) = esub+ebin
!      IF(ist==9) WRITE(7,'(a,i4,2(1pg13.5))') 'HPSI:',ist,
!     &    ebin,esub
c
c     clean up the wavefunction from diagonal components
c     ?? note: this may turn out to be unnecessary. Check!
c
      if(ifclean.ne.1) return

      do is=1,2
      do ir=0,nr
      do iz=0,nz
        psiout(iz,ir,is) = psiout(iz,ir,is) -esub*psiin(iz,ir,is)
      enddo
      enddo
      enddo
c

codd
      return
      end
c
c

      subroutine hpsiodd(psiin,psiout,ist)
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'observ.inc'
*
      dimension psiin(0:nzd,0:nrd,nsd),psiout(0:nzd,0:nrd,nsd)
      dimension drpotCz(0:nzd,0:nrd)
      dimension drpotCr(0:nzd,0:nrd)
      dimension dzpotCz(0:nzd,0:nrd)
      dimension dzpotCr(0:nzd,0:nrd)
*
      data ifclean/1/
*
      iq      = icharg(ist)
      ebin    = spenrg(ist)
c
c     first, zero the output wavfunction because 'psiout' is
c     accumulated during the next spin loop.
c     ?? note: it may be more efficient to introduce workspaces
c              which distinguish spin 'is'. This would allow
c              to avoid costly zeroing and to avoid the spin loop.
c

c
      call laplace_spin(psiin,psi1,nparz(1,ist),mpsi(ist),1,one)
      call nablar_spin(psiin,psi2,   mpsi(ist),1,one)
      call nablaz_spin(psiin,psi3,nparz(1,ist),1,one)
      call nablaz_scalar(potC(0,0,iq,1), dzpotCz,1,1,one)
      call nablaz_scalar(potC(0,0,iq,2), dzpotCr,1,1,one)
      call nablar_scalar(potC(0,0,iq,1), drpotCz,0,1,one)
      call nablar_scalar(potC(0,0,iq,2), drpotCr,0,1,one)


****************************************************
*  time-odd part of the hamiltonian
*
****************************************************      

	do is=1,2
	do ir=0,nr
	do iz=0,nz
	   psiout(iz,ir,is) = 0.0D0
	enddo
	enddo
	enddo




      do is=1,2      
      mm        = mpsi(ist)+is-1
      mmm      =  mm + 3 - 2*is
      if(is.eq.1) then
       sigis = one
      else
       sigis = -1.0
      endif
      do ir=0,nr
        r   = (sr+ir)*dr
        rinv = one/r
	do iz=0,nz
         isi = 3-is


*************************************
*        S*sigma 
*************************************

 	psiout(iz,ir,is)=psiout(iz,ir,is)
     &        +(potS(iz,ir,iq,2)*psiin(iz,ir,isi)
     & +sigis*potS(iz,ir,iq,1)*psiin(iz,ir,is))
ctest
c 	psiout(iz,ir,is)=psiout(iz,ir,is)
c     &        +(sigis*potS(iz,ir,iq,2)*psiin(iz,ir,isi)
c     & -potS(iz,ir,iq,1)*psiin(iz,ir,is))
ctest


*************************************
*        A_q 
*************************************
         if(ir.eq.0) then

* at r = 0 nothing has to be done
* if mm = 0, then the product is zero, else 
* A and psi are at least linear in r

         else

           fac = (one)*potA(iz,ir,iq)*mm*rinv	
           psiout(iz,ir,is)=psiout(iz,ir,is)
     &             + fac*psiin(iz,ir,is)

         endif
*************************************
* C_q 
*************************************
          psiout(iz,ir,is)=psiout(iz,ir,is) 
     & + (3-2*is)*(potC(iz,ir,iq,1)*psi1(iz,ir,is)
     & + drpotCz(ir,ir)*psi2(iz,ir,is)
     & + dzpotCz(iz,ir)*psi3(iz,ir,is))
     & + potC(iz,ir,iq,2)*psi1(iz,ir,isi)
     & + drpotCr(iz,ir)*psi2(iz,ir,isi)
     & + dzpotCr(iz,ir)*psi3(iz,ir,isi)
*     & + last term, two cases

        if(ir.gt.zero) then
	   psiout(iz,ir,is)=psiout(iz,ir,is) 
     & + sigis*(rinv*potC(iz,ir,iq,1)*(psi2(iz,ir,is)
     &  - rinv*mm*mm*psiin(iz,ir,is)))
     &  + sigis*mmm*potC(iz,ir,iq,2)*psiin(iz,ir,isi)*rinv*rinv
     &  + potC(iz,ir,iq,2)*rinv*psi2(iz,ir,isi)
     &  - potC(iz,ir,iq,2)*rinv*rinv*mmm*mmm*
     &  psiin(iz,ir,isi)
        else
          if(abs(mm).eq.1)then
          psiout(iz,ir,isi)=psiout(iz,ir,isi)
     & + sigis*drpotCr(iz,ir)*psi2(iz,is,is)
     & * sign(1,mpsi(ist))
          endif
        endif

        enddo
      enddo
      enddo

      return
      end

c
      subroutine hhpsi(psiin,psiout,ist)
*
***********************************************************************
*                                                                     *
*     same as subroutine hpsi -- but evaluates the variances by       *
*     <ist|HHist> - <ist|Hist>**2                                     *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'derivs.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'observ.inc'
*
      dimension psiin(0:nzd,0:nrd,nsd),psiout(0:nzd,0:nrd,nsd)
*
      data ifclean/1/
*
***********************************************************************
*
      call restore_boundary(psiin,0,+1,ist)
c
      iq      = icharg(ist)
      ebin    = spenrg(ist)
c
      call laplace_spin(psiin,psi1,nparz(1,ist),mpsi(ist),1,one)
      call nablar_spin(psiin,psi2,   mpsi(ist),1,one)
      call nablaz_spin(psiin,psi3,nparz(1,ist),1,one)
c
      mmup      = mpsi(ist)
      mmdown    = mpsi(ist)+1
      nparrup   = 1+mod(mmup,2)
      nparrdown = 1+mod(mmdown,2)
      sigis     = mod(one*(is),two)*two - one
c     
      if(ifstab.eq.1) then
         eferm2 = fer(iq)*oh
         if(ebin.gt.zero) then
            spefac = one
         elseif(ebin.gt.eferm2) then
            spefac = one-ebin/eferm2
         else
            spefac = zero
         endif
      else
         spefac = zero
      endif
c     
      if(spefac.le.zero) then
         do ir=0,nr
         do iz=0,nz
           psiout(iz,ir,1) =
     &      (pot(iz,ir,iq)-qlamb*qfac(iz,ir)-ebin)*psiin(iz,ir,1)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,1)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
           psiout(iz,ir,2) =
     &      (pot(iz,ir,iq)-qlamb*qfac(iz,ir)-ebin)*psiin(iz,ir,2)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,2)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
        enddo
        enddo
      else
         do ir=0,nr
         do iz=0,nz
           psiout(iz,ir,1) =
     &      (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)
     &                    -qlamb*qfac(iz,ir)-ebin)*psiin(iz,ir,1)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,1)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
           psiout(iz,ir,2) =
     &      (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)
     &                    -qlamb*qfac(iz,ir)-ebin)*psiin(iz,ir,2)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,2)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
        enddo
        enddo
      endif
c
      if(t4.ne.zero .OR. ioptsk==13) then        
c
c     the spin-orbit part. action of the s.o. force is computed
c     simply from
c       W_ls |psi> ~ (nabla rho)*sigma sigma*nabla |psi>
c                    - (nabla rho) * nabla |psi>
c
c                               r=0 requires l'Hospital for 1/r*psiin
       if(mmup.eq.1) then
          do iz=0,nz
             psiout(iz,0,1) =
     &      +              potls(iz,0,iq,1)*psi2(iz,0,1)
     &      +              potls(iz,0,iq,1)*psi3(iz,0,2)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,2)
     &      + psiout(iz,0,1)

           psiout(iz,0,2) =
     &      -              potls(iz,0,iq,1)*psi3(iz,0,1)
     &      +              potls(iz,0,iq,2)*psi2(iz,0,1)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,1)
     &      + psiout(iz,0,2)
          enddo
       elseif(mmdown.eq.1) then
           psiout(iz,0,1) =
     &      +              potls(iz,0,iq,1)*psi3(iz,0,2)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,2)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,2)
     &      + psiout(iz,0,1)

           psiout(iz,0,2) =
     &      -              potls(iz,0,iq,1)*psi2(iz,0,2)
     &      -              potls(iz,0,iq,1)*psi3(iz,0,1)
     &      +              potls(iz,0,iq,2)*psi2(iz,0,1)
     &      + psiout(iz,0,2)
       endif
          
c                               now all other r
       do ir=1,nr
c
         r       = (sr+ir)*dr
         rinv    = one/r
         xmmup   = mmup  *rinv
         xmmdown = mmdown*rinv
c
         do iz=0,nz
           psiout(iz,ir,1) =
     &      +        xmmup*potls(iz,ir,iq,1)*psiin(iz,ir,1)
     &      +              potls(iz,ir,iq,1)*psi3(iz,ir,2)
     &      -              potls(iz,ir,iq,2)*psi2(iz,ir,2)
     &      -      xmmdown*potls(iz,ir,iq,2)*psiin(iz,ir,2)
     &      + psiout(iz,ir,1)

           psiout(iz,ir,2) =
     &      -      xmmdown*potls(iz,ir,iq,1)*psiin(iz,ir,2)
     &      -              potls(iz,ir,iq,1)*psi3(iz,ir,1)
     &      +              potls(iz,ir,iq,2)*psi2(iz,ir,1)
     &      -        xmmup*potls(iz,ir,iq,2)*psiin(iz,ir,1)
     &      + psiout(iz,ir,2)
c
         enddo
       enddo
      endif
c
c     compute expectation value 
c
      esub        = wfovlp(psiin,psiout)
      spenrg(ist) = esub+ebin
      spenrg2     = wfovlp(psiout,psiout)
      spevar(ist) = sqrt(abs(spenrg2-esub*esub))
c
c     clean up the wavefunction from diagonal components
c     ?? note: this may turn out to be unnecessary. Check!
c
      if(ifclean.ne.1) return
      do is=1,2
      do ir=0,nr
      do iz=0,nz
        psiout(iz,ir,is) = psiout(iz,ir,is) -esub*psiin(iz,ir,is)
      enddo
      enddo
      enddo
      return
c
c     the action of H^2 on psi
c
      call laplace_spin(psiout,psi1,nparz(1,ist),mpsi(ist),1,one)
      call nablar_spin(psiout,psi2,   mpsi(ist),1,one)
      call nablaz_spin(psiout,psi3,nparz(1,ist),1,one)
c
      mmup      = mpsi(ist)
      mmdown    = mpsi(ist)+1
      nparrup   = 1+mod(mmup,2)
      nparrdown = 1+mod(mmdown,2)
      sigis     = mod(one*(is),two)*two - one
c     
      if(ifstab.eq.1) then
         eferm2 = fer(iq)*oh
         if(ebin.gt.zero) then
            spefac = one
         elseif(ebin.gt.eferm2) then
            spefac = one-ebin/eferm2
         else
            spefac = zero
         endif
      else
         spefac = zero
      endif
c     
      if(spefac.le.zero) then
         do ir=0,nr
         do iz=0,nz
           psi4(iz,ir,1) =
     &      (pot(iz,ir,iq)-qlamb*qfac(iz,ir)-ebin)*psiout(iz,ir,1)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,1)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
           psi4(iz,ir,2) =
     &      (pot(iz,ir,iq)-qlamb*qfac(iz,ir)-ebin)*psiout(iz,ir,2)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,2)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
        enddo
        enddo
      else
         do ir=0,nr
         do iz=0,nz
           psi4(iz,ir,1) =
     &      (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)
     &                    -qlamb*qfac(iz,ir)-ebin)*psiout(iz,ir,1)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,1)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,1)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,1)
           psi4(iz,ir,2) =
     &      (pot(iz,ir,iq)+spefac*potadd(iz,ir,iq)
     &                    -qlamb*qfac(iz,ir)-ebin)*psiout(iz,ir,2)
     &     - bmass(iz,ir,iq)                      *psi1(iz,ir,2)
     &     - drbmas(iz,ir,iq)                     *psi2(iz,ir,2)
     &     - dzbmas(iz,ir,iq)                     *psi3(iz,ir,2)
        enddo
        enddo
      endif
c
      if(t4.ne.zero .OR. ioptsk==13) then        
c
c     the spin-orbit part. action of the s.o. force is computed
c     simply from
c       W_ls |psi> ~ (nabla rho)*sigma sigma*nabla |psi>
c                    - (nabla rho) * nabla |psi>
c
c                               r=0 requires l'Hospital for 1/r*psiout
       if(mmup.eq.1) then
          do iz=0,nz
             psi4(iz,0,1) =
     &      +              potls(iz,0,iq,1)*psi2(iz,0,1)
     &      +              potls(iz,0,iq,1)*psi3(iz,0,2)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,2)
     &      + psi4(iz,0,1)

           psi4(iz,0,2) =
     &      -              potls(iz,0,iq,1)*psi3(iz,0,1)
     &      +              potls(iz,0,iq,2)*psi2(iz,0,1)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,1)
     &      + psi4(iz,0,2)
          enddo
       elseif(mmdown.eq.1) then
           psi4(iz,0,1) =
     &      +              potls(iz,0,iq,1)*psi3(iz,0,2)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,2)
     &      -              potls(iz,0,iq,2)*psi2(iz,0,2)
     &      + psi4(iz,0,1)

           psi4(iz,0,2) =
     &      -              potls(iz,0,iq,1)*psi2(iz,0,2)
     &      -              potls(iz,0,iq,1)*psi3(iz,0,1)
     &      +              potls(iz,0,iq,2)*psi2(iz,0,1)
     &      + psi4(iz,0,2)
       endif
          
c                               now all other r
       do ir=1,nr
c
         r       = (sr+ir)*dr
         rinv    = one/r
         xmmup   = mmup  *rinv
         xmmdown = mmdown*rinv
c
         do iz=0,nz
           psi4(iz,ir,1) =
     &      +        xmmup*potls(iz,ir,iq,1)*psiout(iz,ir,1)
     &      +              potls(iz,ir,iq,1)*psi3(iz,ir,2)
     &      -              potls(iz,ir,iq,2)*psi2(iz,ir,2)
     &      -      xmmdown*potls(iz,ir,iq,2)*psiout(iz,ir,2)
     &      + psi4(iz,ir,1)

           psi4(iz,ir,2) =
     &      -      xmmdown*potls(iz,ir,iq,1)*psiout(iz,ir,2)
     &      -              potls(iz,ir,iq,1)*psi3(iz,ir,1)
     &      +              potls(iz,ir,iq,2)*psi2(iz,ir,1)
     &      -        xmmup*potls(iz,ir,iq,2)*psiout(iz,ir,1)
     &      + psi4(iz,ir,2)
c
         enddo
       enddo
      endif
c
 51   continue
c
c     compute variance
c
      spenrg2     = wfovlp(psiin,psi4)
      spevar(ist) = sqrt(abs(spenrg2-esub*esub))
c
      return
      end
c
c----------------------------------------------------------------
c
      subroutine mfhmlt(iter)
*
***********************************************************************
*                                                                     *
*        compute the mean field hamiltonian .                         *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
      include 'fields.inc'
      include 'derivs.inc'
      include 'levels.inc'
*
      REAL(8), ALLOCATABLE :: workenerg(:,:)
      logical trepeat
      data    trepeat/.false./
c
c     parameters for additional octupole kick in a few first iterations
      data octfac/0.7464D0/
      data octampl/0D0/   !6.0D0/               ! weight for additional octupole kick
      data itmino,itmaxo/3,43/
      REAL(8),PARAMETER :: espherbox=0D0 !20D0   ! to simulate a sppherical box
c
c------------------------------------------------------------------------
c
*
      if(iter.eq.0) then
       redpot = 1.0
       potold = 0.0
      else
       redpot = adypot
       potold = one-adypot
      endif

c
c     first case: just an external oscillator field, for testing
c
      if(iparam.eq.3) then
       amas  = npart(1)+npart(2)
       hboma = 41.0/amas**(edr)
       cfac  = 0.25*hboma*hboma/h2mav
       do ir=0,nr
         r     = (sr+ir)*dr
         r2    = r*r
         do iz=0,nz
           z     = (sz+iz-icm)*dz
           pot(iz,ir,1) = cfac*(r2+z*z)
           pot(iz,ir,2) = pot(iz,ir,1)
           bmass(iz,ir,1) = h2meff(1)
           bmass(iz,ir,2) = h2meff(2)
         enddo
       enddo
c
c     second case: external WS potential 
c
      elseif(iparam.EQ.-1) then
        if(trepeat) return
        call wspot()
        trepeat = .true.
        return
      ELSEIF(ioptsk==13) THEN
!        WRITE(7,'(a)') 'FY branch'
        ALLOCATE(workenerg(0:nzd,0:nrd))
c
c     prepare 'd2rho' and first derivatives on 'psip' and 'psin'
c
        call laplace_scalar(rho(0,0,1),d2rho(0,0,1),+1,0,1,one)
        call laplace_scalar(rho(0,0,2),d2rho(0,0,2),+1,0,1,one)
        call nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,one)
        call nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,one)
        call nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,one)
        call nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,one)
        do ir=0,nr
         r     = (sr+ir)*dr
         r2    = r*r
         do iz=0,nz
          rp      = rho(iz,ir,1)
          rc      = rp                            ! preliminary
          rn      = rho(iz,ir,2)
          tp      = tau(iz,ir,1)
          tn      = tau(iz,ir,2)
          d2p     = d2rho(iz,ir,1)
          d2n     = d2rho(iz,ir,2)
          ssp     = drhols(iz,ir,1)
          ssn     = drhols(iz,ir,2)
          sp      = 0D0                           ! preliminary
          sn      = 0D0                           ! preliminary
          d1pr    = psip(iz,ir,1)
          d1nr    = psip(iz,ir,2)
          d1pz    = psin(iz,ir,1)
          d1nz    = psin(iz,ir,2)
          CALL fyfunc(rp,rn,rc,d1pr,d1pz,d1nr,d1nz,d2p,d2n,tp,tn,
     &                  0D0,0D0,sp,sn,ssp,ssn,
     &                  potp,potn,bmasp,bmasn,bmaslp,bmasln,
     &                  wlspr,wlspz,wlsnr,wlsnz,energd)
          pot(iz,ir,1) = potp
          pot(iz,ir,2) = potn
          potls(iz,ir,1,1) = wlspr
          potls(iz,ir,1,2) = wlspz
          potls(iz,ir,2,1) = wlsnr
          potls(iz,ir,2,2) = wlsnz
          bmass(iz,ir,1) = h2meff(1)+bmasp
          bmass(iz,ir,2) = h2meff(2)+bmasn
          workenerg(iz,ir) = energd
         enddo
        enddo
        energyFy = f1int(workenerg)
        IF(itest>0) WRITE(7,*) '1. energyFy=',energyFy,epair
        CAll fy_pairpot()
        CALL fyfunc_grad()
c
*****************************************************
*
*    time odd hamiltonian (minimal Galilean)
*    S_q, C_q, A_q
*
****************************************************
*
        IF(iSTsym>1) THEN
          DO ir=0,nr
          DO iz=0,nz
!           calculate S_r 
            sp = rotjr(iz,ir,1)
            sn = rotjr(iz,ir,2)
            ss = sp+sn
            potS(iz,ir,1,1) = redpot*(-clsp*ss - clsm*(sp-sn))
     &                        +potold*potS(iz,ir,1,1)
            potS(iz,ir,2,1) = redpot*(-clsp*ss + clsm*(sp-sn))
     &                        +potold*potS(iz,ir,2,1)
!           calculate S_z 
            sp = rotjz(iz,ir,1)
            sn = rotjz(iz,ir,2)
            ss = sp+sn
            potS(iz,ir,1,2) = redpot*(-clsp*ss - clsm*(sp-sn))
     &                        +potold*potS(iz,ir,1,2)
            potS(iz,ir,2,2) = redpot*(-clsp*ss + clsm*(sp-sn))
     &                        +potold*potS(iz,ir,2,2)
!           calculate A_phi 
            rp = rotspden(iz,ir,1)
            rn = rotspden(iz,ir,2) 
            rs = rp+rn
            potA(iz,ir,1) = redpot*(-clsp*rs - clsm*(rp-rn))
     &                      + potold*potA(iz,ir,1)
            potA(iz,ir,2) = redpot*(-clsp*rs + clsm*(rp-rn))
     &                      + potold*potA(iz,ir,2)
!           contribution to energy
            CALL fyfunc_odd(rotjr(iz,ir,1),rotjr(iz,ir,2),
     &                    rotjz(iz,ir,1),rotjz(iz,ir,2),
     &                    spdenr(iz,ir,1),spdenr(iz,ir,2),
     &                    spdenz(iz,ir,1),spdenz(iz,ir,2),
     &                    workenerg(iz,ir))
          END DO
          END DO
!          energspin = energyFy
          energyFy = f1int(workenerg)
!          energspin = energyFy-energspin
!          IF(iter==itmax) 
!     &      WRITE(7,'(a,2(1pg13.5))') 'spin energy=',energspin,energyFy
        END IF

        DEALLOCATE(workenerg)
c
c     add coulomb potential
c
      if (chg.gt.0.0) then
       call coulmb(iter,0)
cPK       call wrt_2D_scalar(coulpt,'field.coul',1,'tst')
       do ir=0,nr
       do iz=0,nz
         wrk1(iz,ir) = crho(iz,ir)*coulpt(iz,ir)          ! ?? test, remove
         pot(iz,ir,1) = redpot*coulpt(iz,ir)+pot(iz,ir,1)
       enddo
       enddo
       ecoul = f1int(wrk1)
        if(mod(iter,iprint).eq.0 .and. itest.ge.1) then
          call testpr(coulpt,0,nzo,0,nro,'coulp')
        endif
      endif
!        WRITE(7,'(a)') 'FY: potentials ready'

      else
c
c   Now the main part for the Skyrme hamiltonian:
c   ---------------------------------------------
c
c     prepare 'd2rho' and first derivatives on 'psip' and 'psin'
c
      call laplace_scalar(rho(0,0,1),d2rho(0,0,1),+1,0,1,one)
      call laplace_scalar(rho(0,0,2),d2rho(0,0,2),+1,0,1,one)
cPK      call mulr(derr2a(0,0,0),rho(0,0,1),d2rho(0,0,1),one,0)
cPK      call mulr(derr2a(0,0,0),rho(0,0,2),d2rho(0,0,2),one,0)
cPK      call mulz(derz2(0,0,1),rho(0,0,1),d2rho(0,0,1),one,1)
cPK      call mulz(derz2(0,0,1),rho(0,0,2),d2rho(0,0,2),one,1)
      if(t4.ne.zero) then
         call nablaz_scalar(rho(0,0,1),psin(0,0,1), 1,1,one)
         call nablaz_scalar(rho(0,0,2),psin(0,0,2), 1,1,one)
         call nablar_scalar(rho(0,0,1),psip(0,0,1), 0,1,one)
         call nablar_scalar(rho(0,0,2),psip(0,0,2), 0,1,one)
cPK     call mulr(derr1(0,0,1),rho(0,0,1),psip(0,0,1),one,0)
cPK     call mulr(derr1(0,0,1),rho(0,0,2),psip(0,0,2),one,0)
cPK     call mulz(derz1(0,0,1),rho(0,0,1),psin(0,0,1),one,0)
cPK     call mulz(derz1(0,0,1),rho(0,0,2),psin(0,0,2),one,0)
      endif
coddtest
*        call nablar_scalar(curphi(0,0,1),rotjz(0,0,1),0,1,one)
*        call nablar_scalar(curphi(0,0,2),rotjz(0,0,2),0,1,one)
*        call nablaz_scalar(curphi(0,0,1),rotjr(0,0,1),1,1,one)
*        call nablaz_scalar(curphi(0,0,2),rotjr(0,0,2),1,1,one)
coddtest
c
c     compose combination of skyrme parameters
c
      t0a    = t0*(one+oh*x0)
      t0b    = t0*(oh+x0)
      t3a    = t3*(one+oh*x3)/6.0
      t3b    = t3*(oh+x3)/6.0
      t3as   = oh*power*t3a
      t3bs   = oh*power*t3b
      t1p2   = 0.25*(t1+t1*x1*0.5+t2+t2*x2*0.5)
      t1m2   = 0.25*(0.5*t1+t1*x1-0.5*t2-t2*x2)
      t13p2  = 0.125*(1.5*t1+t1*x1*3.0+0.5*t2+t2*x2)
      t13m2  = 0.125*(3.0*t1+t1*x1*1.5-t2-t2*x2*0.5)
      b4     = oh*t4
codd
        
      b6     = oh*oh*t0*x0
      b6p    = oh*oh*t0
      b7     = oh*oh*oh*t3*x3/3
      b7p    = oh*oh*oh*t3/3
      b8     = oh*oh*oh*oh*(3*t1*x1-t2*x2)
      b8p    = oh*oh*oh*oh*(3*t1+t2)

      if(isTSym.le.2) then
        b6  = zero
        b6p = zero
        b7  = zero
        b7p = zero
        b8  = zero
        b8p = zero
      elseif(isTSym.eq.3) then
        b8  = zero
        b8p = zero
      endif

! Switch additional potentials on after some iterations
       if(iter.ge.25) then
         fac = oh
       elseif(iter.ge.50) then
         fac = oh+oh
       else 
         fac = zero
       endif
       b6  = fac*b6
       b6p = fac*b6p
       b7  = fac*b7
       b7p = fac*b7p    
       b8  = fac*b8
       b8p = fac*b8p   

  

codd
      excfac = -1.418033*chgx               ! coulomb exchange in Slater appr.
c                                             for t_1,t_2 contr. to spin-orbit
      if(ifls.eq.1) then
        t1m20  = (t1-t2)/8.0
        t1p2x  = (x1*t1+x2*t2)/8.0
      else
        t1m20  = 0.0
        t1p2x  = 0.0
      endif
c
c     compose potentials
c

      do 10 iq=1,2            ! ?? note: double work, remove this loop
                              !          and make iq=1,2 explicitely.

      ic      = 3-iq
      do 10 ir=0,nr
      do 10 iz=0,nz
      rp      = rho(iz,ir,iq)
      rm      = rho(iz,ir,ic)
      rs      = max(rm+rp,smal)
      rspow   = (rm+rp)**power
c      rspow   = rs**power
      tp      = tau(iz,ir,iq)
      ts      = tau(iz,ir,ic)+tp
      d2p     = d2rho(iz,ir,iq)
      d2s     = d2rho(iz,ir,ic)+d2p
      ssp     = drhols(iz,ir,iq)
      sss     = drhols(iz,ir,ic)+ssp
c
codd
      rsigp    = spdenr(iz,ir,iq)
      rsigs    = spdenr(iz,ir,ic)+rsigp
      zsigp    = spdenz(iz,ir,iq)
      zsigs    = spdenz(iz,ir,ic)+zsigp

      sig2     = rsigs*rsigs+zsigs*zsigs
      sigq2    = rsigp*rsigp+zsigp*zsigp
      sigc2    = spdenz(iz,ir,ic)*spdenz(iz,ir,ic)
     &         + spdenr(iz,ir,ic)*spdenr(iz,ir,ic)


codd
codd
* terms added: -b7*alpha*rho^{alpha-1}*sigma^2+b7'*alpha*rho^{alpha-1}*Sum_q sigma_q^2
codd
      pot(iz,ir,iq) = redpot*
     &               ( t0a*rs-t0b*rp  +t1p2*ts-t1m2*tp
     &          -t13m2*d2s+t13p2*d2p -b4*sss-b4p*ssp
     &          +rspow*(t3a*rs-t3b*rp +t3as*rs
     &                  -t3bs*(rp*rp+rm*rm)/rs ) 
     &          +b7*power*(rspow/rs)*sig2 
     &          -b7p*power*(rspow/rs)*(sigq2+sigc2)  )
     &               + potold*pot(iz,ir,iq)
      if(iq.eq.1 .and. chgx.ne.zero)
     &   pot(iz,ir,1) = redpot*excfac*max(rp,smal)**edr + pot(iz,ir,1)
ctest wrk1(iz,ir) = -b4*sss-b4p*ssp *rho(iz,ir,iq)     ! sign befoore b4
c
      IF(ipair>6 .AND. iq==2) THEN
        addparp = (-pairfc(1)* paircr(iz,ir,1)* paircr(iz,ir,1)
     &             -pairfc(2)* paircr(iz,ir,2)* paircr(iz,ir,2))/rho0pr
        pot(iz,ir,1) = -redpot*addparp+pot(iz,ir,1)
        pot(iz,ir,2) = -redpot*addparp+pot(iz,ir,2)
!        IF(iz==1 .AND. IR==1) WRITE(7,'(a,1pg13.5)') 'addparp=',addparp
      END IF
      bmass(iz,ir,iq) = redpot*(h2meff(iq) +t1p2*rs-t1m2*rp)
     &                 + potold*bmass(iz,ir,iq)
c
      if(t4.ne.zero) then
       d1pr    = psip(iz,ir,iq)
       d1sr    = psip(iz,ir,ic)+d1pr
       d1pz    = psin(iz,ir,iq)
       d1sz    = psin(iz,ir,ic)+d1pz
       rlsrp   = rholsr(iz,ir,iq)
       rlszp   = rholsz(iz,ir,iq)
       potls(iz,ir,iq,1) = redpot*(b4*d1sr+b4p*d1pr+
     &                     t1m20*rlsrp-t1p2x*(rlsrp+rholsr(iz,ir,ic)))
     &                    + potold*potls(iz,ir,iq,1)
       potls(iz,ir,iq,2) = redpot*(b4*d1sz+b4p*d1pz+
     &                     t1m20*rlszp-t1p2x*(rlszp+rholsz(iz,ir,ic)))
     &                    + potold*potls(iz,ir,iq,2)
      endif
c

*****************************************************
*
*    time odd hamiltonian
*    S_q, C_q, A_q
*
****************************************************

      b5 = t1p2x
      b5p = t1m20
*  calculate S_r 
      sp = rotjr(iz,ir,iq)
      ss = rotjr(iz,ir,ic) + sp
      taup = spdenkr(iz,ir,iq) 
      taus = spdenkr(iz,ir,ic)+taup
      d2rsigp    = d2spdenr(iz,ir,iq)
      d2rsigs    = d2spdenr(iz,ir,ic)+d2rsigp

      potS(iz,ir,iq,2) = redpot*((-b4)*ss - b4p*sp+
     & b5*taus-b5p*taup
     & +2.0*(b6+b7*rspow+dglan-dglanp)*rsigs    !  add extra g_0 "dglan" as "dglan*rsigs"
c     & -2*(b6p*rsigp+b7p*rspow*rsigp) 
     & -2.0*(b6p+b7p*rspow-2D0*dglanp)*rsigp 
     & -b8*d2rsigs + b8p*d2rsigp)
     & +potold*potS(iz,ir,iq,2)

*  calculate S_z 
      sp = rotjz(iz,ir,iq)
      ss = rotjz(iz,ir,ic) + sp
      taup = spdenkz(iz,ir,iq)
      taus = spdenkz(iz,ir,ic)+taup
      d2zsigp    = d2spdenz(iz,ir,iq)
      d2zsigs    = d2spdenz(iz,ir,ic)+d2zsigp

      potS(iz,ir,iq,1) = redpot*((-b4)*ss - b4p*sp+
     & b5*taus-b5p*taup
     & +2.0*(b6*zsigs+b7*rspow*zsigs)+dglan*zsigs
c     & -2*(b6p*zsigp+b7p*rspow*zsigp)
     & -2.0*(b6p+b7p*rspow)*zsigp
     & -b8*d2zsigs + b8p*d2zsigp)
     & +potold*potS(iz,ir,iq,1)

*  calculate C_r 

      sp = spdenr(iz,ir,iq)
      ss = spdenr(iz,ir,ic)+sp

      potC(iz,ir,iq,2)=redpot*(b5*ss-b5p*sp)
     & + potold*potC(iz,ir,iq,2)

*  calculate C_z 
      sp = spdenz(iz,ir,iq)
      ss = spdenz(iz,ir,ic)+sp

      potC(iz,ir,iq,1)=redpot*(b5*ss-b5p*sp)
     & + potold*potC(iz,ir,iq,1)

*   calculate A_phi 
      cp = curphi(iz,ir,iq)
      cs = curphi(iz,iz,ic)+cp 
      rp = rotspden(iz,ir,iq)
      rs = rotspden(iz,ir,ic) + rp

      potA(iz,ir,iq)=redpot*(-2*t1p2*cs + 2*t1m2*cp 
     & - b4*rs - b4p*rp)
     & + potold*potA(iz,ir,iq)

   10 continue
ctest if(itest.eq.1 .and. t4.ne.zero .and. mod(iter,iprint).eq.0)
ctest&  write(7,'(a,g12.4)') ' spin-orbit energy=',f1int(wrk1)
c
c     add coulomb potential
c
      if (chg.gt.0.0) then
       call coulmb(iter,0)
cPK       call wrt_2D_scalar(coulpt,'field.coul',1,'tst')
       do ir=0,nr
       do iz=0,nz
         wrk1(iz,ir) = crho(iz,ir)*coulpt(iz,ir)          ! ?? test, remove
         pot(iz,ir,1) = redpot*coulpt(iz,ir)+pot(iz,ir,1)
       enddo
       enddo
       ecoul = f1int(wrk1)
        if(mod(iter,iprint).eq.0 .and. itest.ge.1) then
          call testpr(coulpt,0,nzo,0,nro,'coulp')
        endif
      endif
c
c     add steep wall at the boundaries to prevent particle loss
c
      IF(espherbox>0D0) THEN
        do ir=0,nr
        do iz=0,nz
          IF(SQRT(1D0*ir**2+iz**2).GE.nr) 
     &      pot(iz,ir,1:2) = espherbox+pot(iz,ir,1:2)
        enddo
        enddo
      END IF
c      if (1.eq.1) then
c      stepp = 5.0                  !  defines height of wall, may be changed
c      addpr = zero
c      do ir=nr-2,nr
c        addpr = stepp+addpr
c        stepp = stepp+stepp
c        addpt = redpot*addpr
c        do iz=0,nz
c          pot(iz,ir,1) = addpt+pot(iz,ir,1)
c          pot(iz,ir,2) = addpt+pot(iz,ir,2)
c        enddo
c      enddo
c      stepp = 5.0                  !  defines height of wall, may be changed
c      addpz = zero
c      do iz=nzm1,nz
c        addpz = stepp+addpz
c        addpt = redpot*addpz
c        do ir=0,nr
c          pot(iz,ir,1) = addpt+pot(iz,ir,1)
c          pot(iz,ir,2) = addpt+pot(iz,ir,2)
c        enddo
c      enddo
c      if(isym.eq.0)then
c      do iz=0,1
c        addpz = stepp+addpz
c        addpt = redpot*addpz
c        do ir=0,nr
c          pot(iz,ir,1) = addpt+pot(iz,ir,1)
c          pot(iz,ir,2) = addpt+pot(iz,ir,2)
c        enddo
c      enddo
c      endif
c
c      endif
c
c     add steep wall at r**2 + z**2 > R0**2 to 
c     simulate radial symmetric box
c    
c      do ir=0,nr
c         do iz=0,nz
c            pot(iz,ir,1) = redpot*potdamp(iz,ir) + pot(iz,ir,1)
c            pot(iz,ir,2) = redpot*potdamp(iz,ir) + pot(iz,ir,2)
c         enddo
c      enddo     
c
c
c     optionally add octupole kick to explore non-symmetric shapes
c
c      write(6,*)'MFHMLT:iter,beta3,econ,devh2=',iter,beta3,econ,devh2
      if(isym.eq.0 .and. octampl.gt.1D-10 .and. 
     &   abs(beta3).lt.4.0D-2) then
      if(iter.gt.itmino .and. iter.lt.itmaxo) then
c .and.
c     &   econ.gt.10.0D0*errenr .and. devh2.gt.10.0D0*errvar) then
        write(6,*) ' octupole kick active '
        do ir=0,nr
          r2 = (ir+sr)*dr
          r2 = r2*r2
          do iz=0,nz
            z = (iz+sz-icm)*dz - zcm
            addpt = (octampl*octfac)*atan(z*((z*z)-1.5D0*r2))
            pot(iz,ir,1) = addpt+pot(iz,ir,1)
            pot(iz,ir,2) = addpt+pot(iz,ir,2)
          enddo
        enddo
      endif
      endif
c
c     adjust potential at iz=0 in case of isym.eq.0
c
      if (isym.eq.0) then
         do ir=0,nr
            pot(0,ir,1)=0.0D0
            pot(0,ir,2)=0.0D0
            potls(0,ir,1,1)=0.0D0
            potls(0,ir,2,1)=0.0D0
            potls(0,ir,1,2)=0.0D0
            potls(0,ir,2,2)=0.0D0
         enddo
      endif

c
      endif     ! end of Skyrme branch
c
c     compute derivatives of 'bmass'
c
      do  ir=0,nr
      do  iz=0,nz
        wrk1(iz,ir) = bmass(iz,ir,1)- h2meff(1)
        wrk2(iz,ir) = bmass(iz,ir,2)- h2meff(2)
      enddo
      enddo

      call nablaz_scalar(wrk1,dzbmas(0,0,1), 1,1,one)
      call nablaz_scalar(wrk2,dzbmas(0,0,2), 1,1,one)
      call nablar_scalar(wrk1,drbmas(0,0,1), 0,1,one)
      call nablar_scalar(wrk2,drbmas(0,0,2), 0,1,one)
c
c     prints at test level
c
!      WRITE(7,'(a,2i5)') 'TEST FIELDS:',iter,iplot
      if(itest.ge.1 .and. mod(iter,iprint).eq.0) then
        write(7,'(a,i3,a)') ' potentials in mfhmlt (iter=',iter,'):'
        call testpr(pot,0,nzo,0,nro,'  pot')
        call testpr(bmass,0,nzo,0,nro,'  bmass')
        call testpr(drbmas,0,nzo,0,nro,' drbmas')
        if(t4.ne.zero .OR. ioptsk==13) then
            call testpr(potls(0,0,1,1),0,nzo,0,nro,'potls_r')
            call testpr(potls(0,0,1,2),0,nzo,0,nro,'potls_z')
        endif
      endif
c
c     prints
c
      if (mod(iter,iplot).eq.0.and.iter.gt.0) then
        write(18,'(3(a,i3),4(/a),a/a)')
     &  ' Z=',npart(1),' N=',npart(2),' ITER=',iter,
     &  'X:z [fm]',
     &  'Y:fields',
     &  'S: DIST(4.0) SCAL(0.8) THICKN(0.6) XMIN(0.0) ',
     &  '   z     U_proton   U_neutron    B_proton   B_neutron ',
     &  '   rho_p    rho_n',
     &  'H:   X   Y1,L   Y2,P'
        do iz=0,nz
          write(18,'(1x,f5.2,8(1pg13.5))')
     &    (iz+sz-icm)*dz,pot(iz,0,1),pot(iz,0,2),
     &    bmass(iz,0,1),bmass(iz,0,2),rho(iz,0,1),rho(iz,0,2),
     &    drhols(iz,0,:)
        enddo
        write(18,'(/4(/a),a/a)')
     &  'X:r [fm]',
     &  'Y:fields',
     &  'S: DIST(4.0) SCAL(0.8) THICKN(0.6) XMIN(0.0) ',
     &  '   r     U_proton   U_neutron    B_proton   B_neutron ',
     &  '   rho_p    rho_n',
     &  'H:   X   Y1,L'
        do ir=0,nr
          write(18,'(1x,f5.2,6g12.4)')
     &    (ir+sr)*dr,pot(0,ir,1),pot(0,ir,2),
     &    bmass(0,ir,1),bmass(0,ir,2),rho(0,ir,1),rho(0,ir,2)
        enddo
        if(ifstab.eq.1) then
          write(18,'(4(/a))')
     x      'X: r (fm) ',
     x      'Y: potadd (1/fm**3) ',
     x      'S: SCAL(0.9) THICKN(0.6) XMIN(0.0) ',
     x      'H: X   Y1,L   Y2,P  '
          write(18,'(1x,f6.3,2g12.4)')
     x          ((ir+sr)*dr,
     x           potadd(0,ir,1),potadd(0,ir,2),
     x           ir=0,nr)
          write(18,'(4(/a))')
     x      'X: z (fm) ',
     x      'Y: potadd (1/fm**3) ',
     x      'S: SCAL(0.9) THICKN(0.6) XMIN(0.0) ',
     x      'H: X   Y1,L   Y2,P  '
          write(18,'(1x,f6.3,2g12.4)')
     x          ((iz+sz-icm)*dz,
     x           potadd(iz,0,1),potadd(iz,0,2),
     x           iz=0,nz)
        endif
      endif
c
      if(mod(iter,iprint).eq.10000) then
        call testpr(pot,0,nz,0,nr,'pot')
        call testpr(bmass,0,nz,0,nr,'bmass')
        call testpr(drbmas,0,nz,0,nr,'drbmas')
        call testpr(dzbmas,0,nz,0,nr,'dzbmas')
        call testpr(rho,0,nz,0,nr,'rho')
        call testpr(tau,0,nz,0,nr,'tau')
        call testpr(d2rho,0,nz,0,nr,'d2rho')
      endif
c
      return
      end
c-----stapot-----------------------------------------------------------
c
      subroutine stapot(delr)
c
c     computes stabilizing potential for posistive energy states.
c
c     A deformation dependent "coordinate" is deduced from the density
c     of surface states, accumulated on 'potadd' as
c
c      effrad  =  max(0,-log(potadd/rhocri))
c
c      rhocri  =  rho0*rhocut
c
c     where 'rhocri' (recommended 0.08 or 0.07) sets the fractional
c     density at the outer surface.
c     The slope of 'effrad' is then scaled to 1 and the desired
c     slope 'stslop' of the additional potential is reconstructed.
c
c     The additional potential is switched with the particle energy as
c
c                 / 0                   for  epsilon < -wext
c      theta_w = {  (1-|epsilon|)/wext  for  -wext < epsilon < 0
c                 \ 1                   for  0 < epsilon
c
c     the energy-width 'wext' is choosen as
c
c      wext = epsilon_F/2
c
c     This switch is implemented in 'hpsi'.
c
c     The use of the additional potential requires an energy sorted
c     orthonormalization.
c
c
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      include 'observ.inc'
      include 'levels.inc'
      include 'wavefr.inc'
c
      dimension delr(0:nzd,0:nrd)   ! field for distance to surface
c
c------------------------------------------------------------------------------
c
      do 10 iq=1,2
c
c     normalize stabilizing density
c
      xnorm1 = one/f1int(rhosur(0,0,iq))
      do ir=0,nr
        r = (ir+sr)*dr
        do iz=0,nz
          delr(iz,ir) = xnorm1*rhosur(iz,ir,iq)
          z     = (iz+sz-icm)*dz
          wrk1(iz,ir) = (z*z+r*r)*delr(iz,ir)
        enddo
      enddo
      rad2ad = 1.666666667*f1int(wrk1)
c
c     prepare parameters
c
c      rho0   = 0.159
c      rhoprt = rho0*npart(iq)/anumb    ! average density for particles 'iq'
      rhoprt = 0.238732414*rad2ad**(-0.66666667)
      rhocri = rhocut*rhoprt            ! density defining surface
      rhotrd = rhoprt/3.0               ! density defining interior
c
c     compute unscaled effective radius
c
      do ir=0,nr
      do iz=0,nz
        rhoact = max(delr(iz,ir),smal)
        if(rho(iz,ir,iq).gt.rhotrd) rhoact = rhocri+rhocri
        delr(iz,ir) = max(zero,-log(rhoact/rhocri))
      enddo
      enddo
c
c     compute proper scalaing and slope
c
      iz1 = nz
      do iz=nz-1,0,-1
        if(delr(iz,0).le.zero) goto 39
        iz1 = iz
      enddo
   39 continue
      if(iz1.gt.nz-3) stop ' grid-box to small for surfac density'
      if(iz1.le.0) stop ' something is wrong with the surface density'
      slope  = (delr(iz1+2,0)-delr(iz1,0))/(dz+dz)
      slpfac = one/slope
      do ir=0,nr
      do iz=0,nz
        delr(iz,ir)      = slpfac*delr(iz,ir)
        potadd(iz,ir,iq) = stslop*delr(iz,ir)
      enddo
      enddo
      write(6,'(a,8g12.5)')
     & ' potadd:',(potadd(iz,0,1),iz=1,8)
c
   10 continue
c
      return
      end
c
c=====densit===========================================================
c
      subroutine densit(psiw,weight,ist,rhod,taud,divlsd,rlsr,rlsz)
*
***********************************************************************
*                                                                     *
*        accumulates the densities rhod,taud,divlsd, calculated       *
*        for the state in wave function 'psiw' with the occupation    *
*        probability 'weight'. This state is supposed to be related   *
*        to the global label 'ist', i.e.                              *
*          ist    =  absolute adress of the state, to be used in      *
*                    all fields running in the range (1:npsid)        *
*                                                                     *
*        Input parameters:                                            *
*          psiw        s.p. state with label ist                      *
*          weight      occupation probability of the state ist        *
*                      note: this does not include the multiplicity   *
*                            'wate', it is accounted for inside here  *
*          ist         label of state                                 *
*          rhod        density                 sum |psiw|^2           *
*          taud        kinetic density         sum |grad psiw|^2      *
*          divlsd      divergence ls-density   sum                    *
*          rlsr        ls-density r-component                         *
*          rlsz        ls-density z-component                         *
*        Note that the fields for the density do no carry the         *
*        isosipn label within this routine. The appropriate           *
*        placement needs to arranged when calling.                    *
*                                                                     *
*        M.B. 04.05.1995                                              *
*        revised P.-G.R. 19.08.95
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'wavefr.inc'
      include 'observ.inc'
*     include 'resden.inc'
*
      dimension psiw(0:nzd,0:nrd,nsd)
      dimension rhod(0:nzd,0:nrd)
      dimension taud(0:nzd,0:nrd)
      dimension divlsd(0:nzd,0:nrd)
      dimension rlsr(0:nzd,0:nrd)
      dimension rlsz(0:nzd,0:nrd)
*
*-----------------------------------------------------------
*
      fpair  = weight*wate(ist)

*
*---- some combinations of <j_z>
*
      mm     = mpsi(ist)
      xm2    = mm*mm
      xm2l   = (mm+1)*(mm+1)
      nparr  = 1+mod(mm,2)
      nparrc = 3-nparr
*
*---- calculate derivatives of psiw
*
      call nablaz_spin(psiw,psin,nparz(1,ist),1,one)
codd      call nablar_spin(psiw,psip,   mpsi(ist),1,one)
      call nablar_spin(psiw,psip,mpsi(ist),1,one)

cPK      call mulz(derz1(0,0,nparz(1,ist)),psiw(0,0,1),
cPK     &          psin(0,0,1),one,0)
cPK      call mulr(derr1(0,0,nparr),psiw(0,0,1),
cPK     &          psip(0,0,1),one,0)
cPK      call mulz(derz1(0,0,nparz(2,ist)),psiw(0,0,2),
cPK     &          psin(0,0,2),one,0)
cPK      call mulr(derr1(0,0,nparrc),psiw(0,0,2),
cPK     &          psip(0,0,2),one,0)
*
*---- calculate rhod and taud --------------------------------------------
*
*---- value at 'r=0'
*
      if(mm.eq.0) then
        do iz=0,nz
          rhod(iz,0) = fpair * psiw(iz,0,1)*psiw(iz,0,1)
     &                     + rhod(iz,0)
          addtau = fpair * ( psin(iz,0,1)*psin(iz,0,1)
     &                     + 2.0*psip(iz,0,2)*psip(iz,0,2) )
          wrk1(iz,0) = addtau                 ! save that for spin-orbit
          taud(iz,0) = addtau  + taud(iz,0)
        enddo
      elseif(mm.eq.1) then
        do iz=0,nz
           addtau = (fpair+fpair) * ( psip(iz,0,1)*psip(iz,0,1) )
           wrk1(iz,0) = addtau                 ! save that for spin-orbit
           taud(iz,0) = addtau + taud(iz,0)
        enddo
      endif
codd Zeitgespiegelter Fall
      if(mm.eq.(-1)) then
        do iz=0,nz
          rhod(iz,0) = fpair * psiw(iz,0,2)*psiw(iz,0,2)
     &                     + rhod(iz,0)
          addtau = fpair * ( psin(iz,0,2)*psin(iz,0,2)
     &                     + 2.0*psip(iz,0,1)*psip(iz,0,1) )
          wrk1(iz,0) = addtau                 ! save that for spin-orbit
          taud(iz,0) = addtau  + taud(iz,0)
        enddo
      elseif(mm.eq.(-2)) then
        do iz=0,nz
           addtau = (fpair+fpair) * ( psip(iz,0,2)*psip(iz,0,2) )
           wrk1(iz,0) = addtau                 ! save that for spin-orbit
           taud(iz,0) = addtau + taud(iz,0)
        enddo
      endif
codd


*
*---- all other 'r'
*
      do ir=1,nr
         r      = (sr+ir)*dr
         zentfu = xm2/(r*r)
         zentfl = xm2l/(r*r)
         do iz=0,nz
            rhoacu =  psiw(iz,ir,1)*psiw(iz,ir,1)
            rhoacl =  psiw(iz,ir,2)*psiw(iz,ir,2)
            rhod(iz,ir) = fpair*(rhoacu+rhoacl) + rhod(iz,ir)
            addtau = fpair * ( psin(iz,ir,1)*psin(iz,ir,1)
     &                        +psin(iz,ir,2)*psin(iz,ir,2)
     &                        +psip(iz,ir,1)*psip(iz,ir,1)
     &                        +psip(iz,ir,2)*psip(iz,ir,2)
     &                        +zentfu*rhoacu+zentfl*rhoacl)
            wrk1(iz,ir) = addtau          ! save that for spin-orbit
*           rtau(iz,ir) = addtau + rtau(iz,ir)
            taud(iz,ir) = addtau + taud(iz,ir)
        enddo
      enddo

*




*----- calculate spin-orbit density ----------------------------------
*
*     uses divlsd = sum_ist |sigma*nabla psiw|^2-|nabla psiw|^2
*
      if(t4.ne.zero .OR. ioptsk==13) then
      if(mm.eq.0) then
        do iz=0,nz
           sppsiu = psin(iz,0,1)+ 2.0*psip(iz,0,2)
           addls = fpair*(sppsiu*sppsiu)-wrk1(iz,0)
           wrk2(iz,0) = -addls
           divlsd(iz,0) = -addls + divlsd(iz,0)
*          rrhols(iz,0) = -addls + rrhols(iz,0)
            rlsr(iz,0) = fpair*(
     &               psiw(iz,ir,1)*(                  psin(iz,ir,2)))
     &                     + rlsr(iz,0)
            rlsz(iz,0) = fpair*(
     &              - psiw(iz,ir,1)*(2.0*psip(iz,ir,2)           ) )
     &                     + rlsz(iz,0)
        enddo
      elseif(mm.eq.1) then
        do iz=0,nz
           addls  = -wrk1(iz,0)
           wrk2(iz,0) = -addls
           divlsd(iz,0) = -addls + divlsd(iz,0)
*          rrhols(iz,0) = -addls + rrhols(iz,0)
        enddo
      else
        do iz=0,nz
           wrk2(iz,0) = 0.0
        enddo
      endif
codd VZ-Wechsel
      if(mm.eq.(-1)) then
        do iz=0,nz
           sppsiu = -psin(iz,0,2)+ 2.0*psip(iz,0,1)
           addls = fpair*(sppsiu*sppsiu)-wrk1(iz,0)
           wrk2(iz,0) = -addls
           divlsd(iz,0) = -addls + divlsd(iz,0)
*          rrhols(iz,0) = -addls + rrhols(iz,0)
            rlsr(iz,0) = fpair*(
     &               -psiw(iz,ir,2)*(                  psin(iz,ir,1)))
     &                     + rlsr(iz,0)
            rlsz(iz,0) = fpair*(
     &               psiw(iz,ir,2)*(2.0*psip(iz,ir,1)           ) )
     &                     + rlsz(iz,0)
        enddo
      elseif(mm.eq.(-2)) then
        do iz=0,nz
           addls  = -wrk1(iz,0)
           wrk2(iz,0) = -addls
           divlsd(iz,0) = -addls + divlsd(iz,0)
*          rrhols(iz,0) = -addls + rrhols(iz,0)
        enddo
      else
        do iz=0,nz
           wrk2(iz,0) = 0.0
        enddo
      endif


codd

      do ir=1,nr
         rinv  = one/((sr+ir)*dr)
         xmup  = mm*rinv
         xmlow = (mm+1)*rinv
         do iz=0,nz
            sppsiu = psin(iz,ir,1)+psip(iz,ir,2)+xmlow*psiw(iz,ir,2)
            sppsil = -psin(iz,ir,2)+psip(iz,ir,1)-xmup*psiw(iz,ir,1)
            addls  = fpair*(sppsiu*sppsiu+sppsil*sppsil)-wrk1(iz,ir)
            wrk2(iz,ir)   = - addls
            divlsd(iz,ir) = -addls + divlsd(iz,ir)
*           rrhols(iz,ir) = -addls + rrhols(iz,ir)
            rlsr(iz,ir) = fpair*(
     &               psiw(iz,ir,1)*(xmup*psiw(iz,ir,1)+psin(iz,ir,2))
     &              -psiw(iz,ir,2)*(psin(iz,ir,1)+xmlow*psiw(iz,ir,2)))
     &                     + rlsr(iz,ir)
            rlsz(iz,ir) = fpair*(
     &                psiw(iz,ir,2)*
     &                    (psip(iz,ir,1)-xmup*psiw(iz,ir,1) )
     &              - psiw(iz,ir,1)*
     &                    (psip(iz,ir,2)+xmlow*psiw(iz,ir,2) ) )
     &                     + rlsz(iz,ir)
        enddo
      enddo
      endif
*
      return
      end
*
! ! c-----densitodd--------------------------------------------------------

      subroutine densitodd(psiw,weight,ist,
     & rrotj,zrotj,rotsig,rsig,zsig,slarsig,slazsig,
     & curr,spdenkinr,spdenkinz)
*
***********************************************************************
*                                                                     *
*        accumulates the time-odd densities 
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
      include 'derivs.inc'
      include 'wavefr.inc'
      include 'observ.inc'
*     include 'resden.inc'
*
      dimension psiw(0:nzd,0:nrd,nsd)
      dimension rrotj(0:nzd,0:nrd)
      dimension zrotj(0:nzd,0:nrd)
      dimension rotsig(0:nzd,0:nrd)
      dimension rsig(0:nzd,0:nrd)
      dimension zsig(0:nzd,0:nrd)
      dimension slarsig(0:nzd,0:nrd)
      dimension slazsig(0:nzd,0:nrd)
      dimension curr(0:nzd,0:nrd)
      dimension spdenkinr(0:nzd,0:nrd)
      dimension spdenkinz(0:nzd,0:nrd)

      dimension psinn(0:nzd,0:nrd,nsd)
      dimension psipp(0:nzd,0:nrd,nsd)

*-----------------------------------------------------------
*
*  in case of calculating one the half spectrum, set 
*  time-odd densitys to zero
*
      if(isTSym.ge.1)then
       fpair  = weight*wate(ist)
      else
       fpair = zero
      endif
c      fpair = zero

*
*---- some combinations of <j_z>
*
      mm     = mpsi(ist)
      xm2    = mm*mm
      xm2l   = (mm+1)*(mm+1)
      nparr  = 1+mod(mm,2)
      nparrc = 3-nparr
*
*---- calculate derivatives of psiw
*
      call nablaz_spin(psiw,psin,nparz(1,ist),1,one)
      call nablar_spin(psiw,psip,   mpsi(ist),1,one)
      call nablaz_spin(psin,psinn,nparz(1,ist),1,one)
      call nablar_spin(psip,psipp,   mpsi(ist),1,one)
      if(ist.eq.1) then
c        call testpr(psin(0,0,1),0,10,0,10,'Dz oben')
c        call testpr(psin(0,0,2),0,10,0,10,'Dz unten')
      endif

***********************************************************************
*         testweise eingef�gt
*
*
***********************************************************************

***********************************************************************
*       rrotj r-comp of rot j
*       zrotj z-comp of rot j
***********************************************************************
*
*----- r = 0
*
*
c      fpair = (-1.0)*fpair
* only z-comp changes
      if(mm.eq.0) then
         do iz=0,nz
           zrotj(iz,0)=zrotj(iz,0) + 2*fpair*psip(iz,0,2)*psip(iz,0,2)
         enddo
      elseif(mm.eq.1) then
         do iz=0,nz
           zrotj(iz,0)=zrotj(iz,0) + 2*fpair*psip(iz,0,1)*psip(iz,0,1)
         enddo
      elseif(mm.eq.(-1)) then
         do iz=0,nz
           zrotj(iz,0)=zrotj(iz,0) - 2*fpair*psip(iz,0,1)*psip(iz,0,1)
         enddo
      elseif(mm.eq.(-2)) then
         do iz=0,nz
           zrotj(iz,0)=zrotj(iz,0) - 2*fpair*psip(iz,0,2)*psip(iz,0,2)
         enddo
      endif

*---- all other 'r'
*
      do ir=1,nr
         rinv  = one/((sr+ir)*dr)
         xmup  = 2*mm*rinv
         xmlow = 2*( mm+1)*rinv
         do iz=0,nz
		zrotj(iz,ir)=zrotj(iz,ir) +
     &                       fpair*(xmup*psiw(iz,ir,1)*psip(iz,ir,1)+
     &                       xmlow*psiw(iz,ir,2)*psip(iz,ir,2))
		rrotj(iz,ir)=rrotj(iz,ir) -
     &                       fpair*(xmup*psiw(iz,ir,1)*psin(iz,ir,1)+
     &                       xmlow*psiw(iz,ir,2)*psin(iz,ir,2))
        enddo
      enddo

c      fpair = (-1.0)*fpair
***********************************************************************
*                Current j
*		curr = phi-comp. Current
***********************************************************************
*---- r=0
* nothing has to be done

      do iz=0,nz
        curr(iz,0)=zero
      enddo
*---- all other r
      do ir=1,nr
         rinv  = one/((sr+ir)*dr)
         fup  = mm*rinv
         flow = (mm+one)*rinv
   	 do iz=0,nz
            add = fpair*(fup*psiw(iz,ir,1)*
     &              psiw(iz,ir,1) + flow*psiw(iz,ir,2)*psiw(iz,ir,2))
	    curr(iz,ir)= curr(iz,ir)+add
	 enddo
      enddo

***********************************************************************
*
*              sigma
*             rsig = r-comp of sigma
*             zsig = z-comp of sigma
*
***********************************************************************
       fpair2= (1.0)*fpair/1.0
      do ir=0,nr
	 do iz=0,nz
            rsig(iz,ir)=rsig(iz,ir)+fpair2*2*psiw(iz,ir,1)*psiw(iz,ir,2)
	    zsig(iz,ir)=zsig(iz,ir)+fpair2*(
     &       psiw(iz,ir,1)*psiw(iz,ir,1) - psiw(iz,ir,2)*psiw(iz,ir,2))
	 enddo
      enddo







***********************************************************************
*
*             rot sigma
*             rotsig = phi-comp of rot sigma
*
***********************************************************************

      do ir=0,nr
         do iz=0,nz
	     rotsig(iz,ir)= rotsig(iz,ir) + 2*fpair*(psiw(iz,ir,1)*
     &                     psin(iz,ir,2) + psiw(iz,ir,2)*psin(iz,ir,1)
     &                     -psiw(iz,ir,1)*psip(iz,ir,1)
     &                     + psiw(iz,ir,2)*psip(iz,ir,2))
	 enddo
      enddo

***********************************************************************
*
*             Laplace sigma
*             larsig = r-comp of Laplace sigma
*             lazsig = z-comp of Laplace sigma
*
***********************************************************************
*
*----- r = 0
*
*
* only z-comp changes

      if((mm.eq.0).or.(mm.eq.(-2))) then
         do iz=0,nz
            slazsig(iz,ir)=slazsig(iz,ir) -
     &    4*fpair*psip(iz,ir,2)*psip(iz,ir,2)
         enddo
      endif
      if((mm.eq.1).or.(mm.eq.(-1))) then
         do iz=0,nz
            slazsig(iz,ir)=slazsig(iz,ir) +
     &    4*fpair*psip(iz,ir,1)*psip(iz,ir,1)
         enddo
      endif

* now all other r
      do ir=1,nr
         rinv  = one/((sr+ir)*dr)
	 do iz=0,nz
	    slarsig(iz,ir)=slarsig(iz,ir)+fpair*(2*
     & (psiw(iz,ir,1)*psipp(iz,ir,2) +psiw(iz,ir,2)*psipp(iz,ir,1)
     & + 2*psip(iz,ir,1)*psip(iz,ir,2)) + rinv*2*
     &(psip(iz,ir,1)*psiw(iz,ir,2) + psip(iz,ir,2)*psiw(iz,ir,1))
     & - rinv*rinv*2*psiw(iz,ir,1)*psiw(iz,ir,2)
     & + 2*(psiw(iz,ir,1)*psinn(iz,ir,2) + psiw(iz,ir,2)*psinn(iz,ir,1)
     & + 2*psin(iz,ir,1)*psin(iz,ir,2)))

	    slazsig(iz,ir)=slazsig(iz,ir) + 2*fpair*(
     & psip(iz,ir,1)*psip(iz,ir,1) -psip(iz,ir,2)*psip(iz,ir,2)
     &+psiw(iz,ir,1)*psipp(iz,ir,1)-psiw(iz,ir,2)*psipp(iz,ir,2)
     &+rinv*(psiw(iz,ir,1)*psip(iz,ir,1)-psiw(iz,ir,2)*psip(iz,ir,2))
     &+psinn(iz,ir,1)*psiw(iz,ir,1)-psinn(iz,ir,2)*psiw(iz,ir,2)
     &+psin(iz,ir,1)*psin(iz,ir,1)-psin(iz,ir,2)*psin(iz,ir,2))

	 enddo
      enddo





*
**********************************************************************
*                kinetic spin density
*
*********************************************************************
* r = 0
      if((mm.eq.0).or.(mm.eq.(-1)))	 then
      ism = abs(mm)+1
      sigism = 3 -2*ism
      do iz=0,nz
* evtl spdenkinr = 0
          spdenkinz(iz,0)=spdenkinz(iz,0)+sigism*fpair*(
     &      psin(iz,0,ism)*psin(iz,0,ism))
      enddo
      endif 
      if((mm.eq.0).or.(mm.eq.(-2))) then
       do iz=0,nz
          spdenkinz(iz,0)=spdenkinz(iz,0)-2*fpair*
     &     psin(iz,0,2)*psin(iz,0,2)
       enddo
      elseif((mm.eq.1).or.(mm.eq.(-1))) then
       do iz=0,nz
          spdenkinz(iz,0)=spdenkinz(iz,0)+2*fpair*
     &     psin(iz,0,1)*psin(iz,0,1)
       enddo
      endif
* now all other r
      do ir=1,nr
         r      = (sr+ir)*dr
         f = mm*(mm+1)/(r*r)
         fu = mm/r
         fl = (mm+1)/r
  	 do iz=0,nz
	    spdenkinr(iz,ir)= spdenkinr(iz,ir)+2*fpair*(psip(iz,ir,1)*
     &         psip(iz,ir,2) + f*psiw(iz,ir,1)*psiw(iz,ir,2)+ 
     &         psin(iz,ir,2)*psin(iz,ir,1))
            spdenkinz(iz,ir)=spdenkinz(iz,ir)+fpair*(psip(iz,ir,1)*
     &         psip(iz,ir,1)+fu*fu*psiw(iz,ir,1)*psiw(iz,ir,1)
     &         + psin(iz,ir,1)*psin(iz,ir,1) - (
     &         psip(iz,ir,2)*psip(iz,ir,2)+fl*fl*psiw(iz,ir,2)*
     &         psiw(iz,ir,2) + psin(iz,ir,2)*psin(iz,ir,2) ) )
            
	 enddo
      enddo      

      return
      end

c-----pairobs-----------------------------------------------------------
c
      subroutine pairobs(sumuv,sumuv2,sumuve)
*
***********************************************************************
*                                                                     *
*        accumulates the a couple of measures of pairing:             *
*         sumuv    = the sum of proton and neutron u*v as a measure   *
*                    for pair current.                                *
*         sumuv2   = the sum of proton and neutron (u*v)**2 as        *
*                    a measure for the variance of pair current.      *
*         sumuve   = the sum of proton and neutron (u*v)*spenergy**2  *
*                    as a measure for the energiy span around the     *
*                    Fermi energy.                                    *
*        information on s.p. occupation and energies comes through    *
*        'levels.inc'.                                                *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
!      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
*
      REAL(8),INTENT(OUT) :: sumuv(kqm),sumuv2(kqm),sumuve(kqm)
*
*
      sumuv=0D0
      sumuv2=0D0
      sumuve=0D0
      DO jst=1,nst(1)+nst(2)

        iq = icharg(jst)
        IF(tpairspace) THEN
          fpair2 = wtpair(jst)*wate(jst)
        ELSE
          fpair2 = wate(jst)
        END IF
        uvact=SQRT(vsqf(jst)-vsqf(jst)**2)
        sumuv(iq)=sumuv(iq)+uvact*fpair2
        sumuv2(iq)=sumuv2(iq)+uvact**2*fpair2
        sumuve(iq)=sumuve(iq)+uvact*(spenrg(jst)-fer(iq))**2*fpair2

      END DO
c
      return
      end
c-----pairobs-----------------------------------------------------------
c
      REAL(8) FUNCTION sum_uv(iq)
*
***********************************************************************
*                                                                     *
*     The sum of proton (iq=1) or neutron (iq=2) u*v as a measure     *
*     for pair current.                                               *
*     Information on s.p. occupation and energies comes through       *
*     'levels.inc'.                                                   *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
!      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
!
      accuv=0D0
      DO jst=1,nst(1)+nst(2)
        IF(iq.NE.icharg(jst)) CYCLE
        IF(tpairspace) THEN
          fpair2 = wtpair(jst)*wate(jst)
        ELSE
          fpair2 = wate(jst)
        END IF
        uvact=SQRT(vsqf(jst)-vsqf(jst)**2)
        accuv=accuv+uvact*fpair2
      END DO
      sum_uv=accuv
c
      return
      end

c-----fdenti-----------------------------------------------------------
c
      subroutine fdenti(psiw,rhoout)
*
***********************************************************************
*                                                                     *
*        accumulates the total local density 'rhoout' for the set     *
*        of wave functions given in 'psiw'.  The states are weighted  *
*        by occupation probability 'vsqf' and phase space weights     *
*        'wtpair' and 'wate'. These three arrays are communicated     *
*        through 'levels.inc'.                                        *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
*
      REAL(8),INTENT(IN) ::  psiw(0:nzd,0:nrd,nsd,npsid) ! set of s.p. wave functions
      REAL(8),INTENT(OUT) :: rhoout(0:nzd,0:nrd,kqm)         ! resulting density
*
      REAL(8),PARAMETER :: lowadd=1D-9
*
*
*
      rhoout=0D0
      facrho=1D0
      DO jst=1,nst(1)+nst(2)

        iq = icharg(jst)
        IF(tpairspace) THEN
          fpair2 = facrho*vsqf(jst)*wtpair(jst)*wate(jst)
        ELSE
          fpair2 = facrho*vsqf(jst)*wate(jst)
        END IF
        if (fpair2.lt.lowadd) CYCLE
        do ir=0,nr
        do iz=0,nz
          rhoout(iz,ir,iq) = fpair2*(psiw(iz,ir,1,jst)*psiw(iz,ir,1,jst)
     &                          +psiw(iz,ir,2,jst)*psiw(iz,ir,2,jst))
     &                    + rhoout(iz,ir,iq)
        enddo
        enddo

      END DO
c
      return
      end


c-----sdenti-----------------------------------------------------------
c
      subroutine sdenti(psiw,rhoacc,addwgt,ist)
*
***********************************************************************
*                                                                     *
*        accumulates just the densities 'rho' in /common/ curren  *
*        for the state in wave function 'psiw' . this state is        *
*        supposed to be related to the global label 'ist', i.e.       *
*          ist    =  absolute adress of the state, to be used in      *
*                    all fields running in the range (1:npsid)        *
*        the fields have to be multiplied with (1.0-addrho) before    *
*        the m-loop in the calling routine.                           *
*                                                                     *
*        this routine is used in connection with the iteration of     *
*        quadrupole deformations.                                     *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
*
      dimension psiw(0:nzd,0:nrd,nsd)      ! wave function to be accumulated
      dimension rhoacc(0:nzd,0:nrd)        ! density on which one accumulates
*
      data lowadd/1.0D-9/
*
*
*
      if (addwgt.lt.lowadd) return
      fpair  = wate(ist)*addwgt
      do ir=0,nr
      do iz=0,nz
        rhoacc(iz,ir) = fpair * (psiw(iz,ir,1)*psiw(iz,ir,1)
     &                        +psiw(iz,ir,2)*psiw(iz,ir,2))
     &                  + rhoacc(iz,ir)
      enddo
      enddo
c
      return
      end

c-----cpairpot-----------------------------------------------------------
c
      subroutine cpairpot(psiw1,psiw2,rhoacc,addwgt,ist)
*
***********************************************************************
*                                                                     *
*        accumulates the just the densities 'rho' in /common/ curren  *
*        for the state in wave function 'psiw' . this state is        *
*        supposed to be related to the global label 'ist', i.e.       *
*          ist    =  absolute adress of the state, to be used in      *
*                    all fields running in the range (1:npsid)        *
*        the fields have to be multiplied with (1.0-addrho) before    *
*        the m-loop in the calling routine.                           *
*                                                                     *
*        this routine is used in connection with the iteration of     *
*        quadrupole deformations.                                     *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'dynpar.inc'
*
      dimension psiw1(0:nzd,0:nrd,nsd)      ! wave function to be accumulated
      dimension psiw2(0:nzd,0:nrd,nsd)      ! wave function partner to be accumulated
      dimension rhoacc(0:nzd,0:nrd)        ! density on which one accumulates
*
      data lowadd/1.0D-9/
*
*
*
c      if (addwgt.lt.lowadd) return
      fpair  = wate(ist)*addwgt
      do ir=0,nr
      do iz=0,nz
        rhoacc(iz,ir) = fpair * (psiw1(iz,ir,1)*psiw2(iz,ir,2)
     &                        -psiw1(iz,ir,2)*psiw2(iz,ir,1))
     &                  + rhoacc(iz,ir)
      enddo
      enddo
c
      return
      end
c
c-----damint-----------------------------------------------------------
*
*
      subroutine damint(rhoin, rhosw, delr)
*
      include 'paramr.inc'
      include 'dynpar.inc'
*
      double precision rhoin, rhosw, delr
      dimension rhoin(0:nzd,0:nrd)
      dimension delr(0:nzd,0:nrd)
*
      integer maxp
      parameter(maxp=nzd*nrd)
*
      double precision surf
      dimension surf(maxp,2)
      integer iz, ir, izz, irr, point, dir, p
      double precision frac, diag, maxdist, z, r, dist  
*
c     if(itest.gt.0)  write(6,'(a,g12.4)')
c    &   ' DAMINT called with rhosw=',rhosw
      diag = sqrt(2.0)
      maxdist = 1.0e10
      point = 1
      do iz=0,nz
        do ir=0,nr
          do 120 dir=1,4 
            if(dir.eq.1 .and. ir.ne.0 .and. iz.ne.0) then 
              izz=iz-1
              irr=ir-1
            elseif(dir.eq.2 .and. iz.ne.0) then 
              izz=iz-1
              irr=ir
            elseif(dir.eq.3 .and. ir.ne.nr .and. iz.ne.0) then 
              izz=iz-1
              irr=ir+1
            elseif(dir.eq.4 .and. ir.ne.nr) then 
              izz=iz
              irr=ir+1
            else 
              goto 120
            endif
            if((rhoin(iz,ir).lt.rhosw .and. rhoin(izz,irr).lt.rhosw) 
     &         .or.
     &         (rhoin(iz,ir).ge.rhosw .and. rhoin(izz,irr).ge.rhosw))
     &         goto120
            frac=(rhosw-rhoin(iz,ir))/(rhoin(izz,irr)-rhoin(iz,ir))
            if(iz.eq.izz) then        
              surf(point,1) = (iz+sz-icm) * dz
            elseif(ir.eq.irr) then
              surf(point,1) = (iz+sz-frac-icm) * dz
            else
              surf(point,1) = (iz+sz-frac/diag-icm) * dz
            endif
            if(ir.eq.irr) then        
              surf(point,2) = (ir+sr) * dr
            elseif(iz.eq.izz) then
              surf(point,2) = (ir+sr+frac) * dr
            elseif(irr.lt.ir) then
              surf(point,2) = (ir+sr-frac/diag) * dr
            else
              surf(point,2) = (ir+sr+frac/diag) * dr
            endif
            point = point + 1       
  120     continue
        enddo
      enddo
c
c     if(itest.gt.0) 
c    &   write(6,'(a,i6)') ' nr. of point in damint=',point
      do iz=0,nz
        z = (iz+sz-icm) * dz
        do ir=0,nr
          r = (ir+sr) * dr
          delr(iz,ir) = maxdist
          do p=1,point-1
            dist = (z-surf(p,1))*(z-surf(p,1)) +
     &             (r-surf(p,2))*(r-surf(p,2))
            if(delr(iz,ir).gt.dist)
     &         delr(iz,ir) = dist
          enddo
          delr(iz,ir) = sqrt(delr(iz,ir))
          if(rhoin(iz,ir).gt.rhosw) then
            delr(iz,ir) = -delr(iz,ir)
          endif
        enddo
      enddo
*
      return
      end
c-----damcst-----------------------------------------------------------
*
      subroutine damcst(iiter,drho0l)
*
c     computes 'distsf'=distance to rho/10 surface
c     and subsequently the damping for the constraint
c
c      iiter   = actual iteration number
c      drho0l  = density at which the surface line is defined (output)
c
      include "paramr.inc"
*
      include "dynpar.inc"
      include "fields.inc"
      include "constr.inc"
*     include "coulox.inc"
*     include "cntl.inc"
*
c      double precision dzmin0, dzmin, z, drho0, drho0l,
c     &                 dzmax, drmax,
c     &                 dcadd, dferm, dtmp(0:nzd,0:nrd)
c      integer iz, ir, iz2, ir2, iiter
      logical lch
*
* laenge der gitterdiagonalen
      dzmin0 = dz*dz * real((nz+1)*(nz+1)) +
     &         dr*dr * real((nr+1)*(nr+1))
*
* maximale dichte
      drho0 = d0
      do ir = 0, nr
      do iz = 0, nz
        drho0 = max(drho0, rho(iz,ir,1)+rho(iz,ir,2))
      enddo
      enddo
      drho0l = drho0/10.0d0
*     write(*,*) 'drho0l', drho0l
* hat sich die oberflaeche des Kerns veraendert?
*
      lch = .false.
      do ir = 0, nr
      do iz = 0, nz
        if( (rho(iz,ir,2)+rho(iz,ir,2).ge.drho0l) .neqv.
     &      lrho(iz,ir) ) then
          lrho(iz,ir) = .not. lrho(iz,ir)
          lch = .true.
        endif
      enddo
      enddo
c      if (iiter.eq.0) then
      if (iiter.gt.-1) then
        lch = .true.
      endif
*      write(*,*) 'lch   ', lch
*
*  entfernung zum kernrand
      if( lch ) then
        do ir = 0, nr
          do iz = 0, nz
            dzmin = dzmin0
            dzmax = dzmin0
            do ir2 = 0, nr
              do iz2 = 0, nz
                z = dz*dz*(iz2-iz)*(iz2-iz) +
     &              dr*dr*(ir2-ir)*(ir2-ir)
                if( lrho(iz2,ir2) ) then
                  dzmin = min(dzmin, z)
                else
                  dzmax = min(dzmax, z)
                endif
              enddo
            enddo
            if( lrho(iz,ir) ) then
              distsf(iz,ir) = -sqrt(dzmax)
            else
              distsf(iz,ir) = sqrt(dzmin)
            endif
          enddo
        enddo
      endif
*           call testpr(dr,0,100,0,100,' ')
*
*  konstruktion der gedaempften mpol-operatoren
c      if( iiter.eq.1 ) then
c        dcadd=1.0d0
c      else
c        dcadd=0.2d0
c      endif
*      write(*,*) 'dcadd', dcadd
*     call testpr(drho(0,0,xvmes),0,100,0,100,'rho')
      do ir = 0, nr
      do iz = 0, nz
        qfac(iz,ir) = qcfac(iz,ir) /
     &               (one+exp( (distsf(iz,ir)-dalp) / dgam ))
      enddo 
      enddo 
*
      return
      end
c
c-----delrmn-------------------------------------------------------------------
c
      subroutine delrmn(rhoin,rhosw,delr)
c
c     to compute the field of minimum distance to the surface line
c     defined by the "switch-density":  rho(r,z) = rhosw.
c     The definition is:
c            2
c     Delta_R (r,z) = Min  [(r-r')^2+(z-z')^2] with rho(r',z') = rhosw
c                     r'z'
c
c     The minimum is searched in two steps:
c     1. The surface line r_sw= r_sw(z) is established,
c     2. For each point on the grid, the minimum distance to
c        the points on r_sw is looked for. This search proceeds
c        in strips of constant z while restarting from the minimum
c        point found previously.
c
c     the routine works presently only for the case SYM=1
c
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'observ.inc'
c
      logical tifzmn                 ! flag for "cap" at lower z
      dimension rhoin(0:nzd,0:nrd)   ! input density
      dimension delr(0:nzd,0:nrd)    ! output field of min. distances
      dimension rsw(-1:nzd)           ! r values for surface line
      dimension zsw(-1:nzd)           ! z values for surface line
c
c
c------------------------------------------------------------------------------
c
c     the program assumes presently the same grid size in r and z
c
      if(dr.ne.dz) stop ' subroutine "DELRMN" requires DR=DZ'
c
      if(itest.ge.2) then
        write(7,'(a,g12.4)')
     &   ' DELRMN called with rhosw=',rhosw
        call testpr(rhoin,0,nz,0,nr,'rhoin')
      endif

c
c     Search surface line and store on 'rsw' with linear interpolation:
c
      do iz=0,nz
        do ir=nr,0,-1
          ira  = ir
          if(rhoin(iz,ir).gt.rhosw) goto 19
        enddo
        ira = -1
   19   continue
        if(ira.eq.-1) then
           rsw(iz) = -nr            ! signals that 'rhoin < rhosw' always
        else
           rsw(iz) = ira
     &           +(rhoin(iz,ira)-rhosw)/(rhoin(iz,ira)-rhoin(iz,ira+1))
           zsw(iz) = iz
        endif
      enddo
c
c     determine minimum and maximum relevant 'iz'
c
      do iz=0,nz
        izmin = iz
        if(rsw(iz).gt.zero) goto 39
      enddo
   39 continue
      tifzmn = izmin.gt.0
      do iz=nz,0,-1
        izmax = iz
        if(rsw(iz).gt.zero) goto 49
      enddo
   49 continue
c
c     determine 'z' values of the end points, reset 'rsw'
c
      if(izmin.gt.0) then
        sqr1 = sqrt(rsw(izmin))
        sqr2 = sqrt(rsw(izmin+1))
        izmin = izmin-1
        addz  = min(one,sqr2/(sqr1+sqr2))
        zsw(izmin) = izmin+addz
        rsw(izmin) = zero
        zsw(izmin-1) = zsw(izmin+1)      ! symm. point for cap-parabola
        rsw(izmin-1) = -rsw(izmin+1)
      else
        izmin = izmin-1
        zsw(-1) = izmin
        if(isym.eq.1)rsw(-1) = rsw(1)                 ! requires SYM=1
        if(isym.eq.0)rsw(-1) = rsw(nz)                ! requires SYM=0
      endif
      if(izmax.lt.nz) then
        sqr1 = sqrt(rsw(izmax))
        sqr2 = sqrt(rsw(izmax-1))
        izmax = izmax+1
        addz  = min(one,sqr1/(sqr2-sqr1))
        zsw(izmax) = izmax-1+addz
        rsw(izmax) = zero
        zsw(izmax+1) = zsw(izmax-1)      ! symm. point for cap-parabola
        rsw(izmax+1) = -rsw(izmax-1)
      endif
      if(itest.ge.2) then
        write(7,'(a,2i5)')
     &   ' izmin,izmax=',izmin,izmax
        write(7,'(a,20(/1x,5g12.4))')
     &   ' rsw:',(rsw(iz),iz=-1,nz)
        write(7,'(a,20(/1x,5g12.4))')
     &   ' zsw:',(zsw(iz),iz=-1,nz)
      endif
c
c     Scan all grid points and determine minimum distance to 'rsw'
c
      if(itest.ge.2) write(7,'(a)')
     &     '  iz  irmbasidirizinizen iz0   drsqmn   delr'
      do iz=0,nz
        mbase = 0              ! relative base point in 'iz' for search
        do ir=0,nr
          iz0     = min(max(izmin+1,iz+mbase),izmax-1)
          izm     = iz0-1
c                             first triplet determines direction of search
          drsq0   = (rsw(iz0)-ir)**2+(zsw(iz0)-iz)**2
          drsqm1  = (rsw(izm)-ir)**2+(zsw(izm)-iz)**2
          drsqp1  = (rsw(iz0+1)-ir)**2+(zsw(iz0+1)-iz)**2
          if(drsqp1.gt.drsq0) then
            if(drsqm1.lt.drsq0) then
              izini  = min(max(iz0-2,izmin),izmax)
              if(izini.ge.0) then
                idir   = -1
                izend  = izmin
                x      = drsqm1          ! order the "drsq.." for interpolation
                drsqm1 = drsqp1
                drsqp1 = x
              else
                idir   = 0
              endif
            else
              idir   = 0             ! signals end of search
            endif
          else
            idir   = +1
            izini  = min(max(iz0+2,izmin),izmax)
            izend  = izmax
          endif
c                             search in 'idir' until distance rises
          if(idir.ne.0) then
            do iz1=izini,izend,idir
              iz0    = iz1
              drsqm1 = drsq0
              drsq0  = drsqp1
              drsqp1 = (rsw(iz0)-ir)**2+(zsw(iz0)-iz)**2
              if(drsqp1.gt.drsq0) goto 29
            enddo
   29       continue            ! minimum embraced in 'drsqm1,drsq0,drsqp1'
          endif
c
c         check outer cap-discontinuity
c
          if(idir.gt.0 .and. iz0.eq.izmax) then
            if(drsqp1.lt.drsqm1) then  ! complement cap
              drsqm1 = drsq0
              drsq0  = drsqp1
              drsqp1 = (rsw(iz0+1)-ir)**2+(zsw(iz0+1)-iz)**2
            else
              drsqp1 = drsq0
              drsq0  = drsqm1
              drsqm1 = (rsw(iz0-3)-ir)**2+(zsw(iz0-3)-iz)**2
c              write(6,'(a,g12.4)') ' complement dist:',drsqm1
            endif
          endif
c
c         check inner cap-discontinuity
c
          if(tifzmn .and. idir.lt.0 .and. iz0.eq.izmin) then
            if(drsqp1.lt.drsqm1) then  ! complement cap
              drsqm1 = drsq0
              drsq0  = drsqp1
              drsqp1 = (rsw(iz0-1)-ir)**2+(zsw(iz0-1)-iz0)**2
            else
              drsqp1 = drsq0
              drsq0  = drsqm1
              drsqm1 = (rsw(iz0+3)-ir)**2+(zsw(iz0+3)-iz0)**2
            endif
          endif
c                             quadratic interpolation to minimum distance
          x      = drsqp1-drsqm1
          drsqmn = drsq0-x*x/(8.0*(drsqp1-drsq0-drsq0+drsqm1))
c                             save and privide proper sign
          delr(iz,ir) = sign(sqrt(abs(drsqmn)),rhosw-rhoin(iz,ir))
          if(itest.ge.2) write(7,'(7i4,2(1pg12.4))')
     &          iz,ir,mbase,idir,izini,izend,iz0,drsqmn,delr(iz,ir)
          mbase  = iz-iz0     ! save minimum position as base for next search
          mbase  = max(0,mbase)
c
        enddo
      enddo
c
      return
      end
c-----damnuc-----------------------------------------------------------
c
      subroutine damnuc(psiin,psiout,ist)
*
***********************************************************************
*                                                                     *
*        kinetic energy damping applied to the the wavefunction       *
*        'psiin', result is returned on 'psiout', 'ist' is the        *
*        absolute position of the state in the spectrum.              *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'derivs.inc'
      include 'levels.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
*
      dimension psiin(0:nzd,0:nrd,nsd),psiout(0:nzd,0:nrd,nsd)
*
c
c     the straightforward imaginary-time step in case of "e0inv=0"
c
      if(e0inv.eq.zero) then
       do is=1,2
       do ir=0,nr
       do iz=0,nz
        psiout(iz,ir,is) = x0inv*psiin(iz,ir,is)
       enddo
       enddo
       enddo
      else
c
c     now the full kinetic energy damping
c
      step   = x0inv/e0inv
      call damp_nucleon(rdnuc,zdnuc,
     &                  psiin,psiout,nparz(1,ist),mpsi(ist),step)
c
c      do is=1,2
c        mm     = mpsi(ist)+is-1
c        call mulz(zdnuc(0,0,nparz(is,ist)),psiin(0,0,is),
c     &            wrk1,one,0)
c        call mulr(rdnuc(0,0,mm),wrk1,psiout(0,0,is),step,0)

c
c       set boundary zeroes if parity requires it.
c
      do is=1,2
        mm = mpsi(ist)+is-1
        if (mm.gt.0) then              ! ABS(mm) ?
          do iz=0,nz
            psiout(iz,0,is)  =  zero
          enddo
        endif
c
        if (nparz(is,ist).eq.2) then
          do ir=0,nr
            psiout(0,ir,is)  =  zero
          enddo
        endif
c
      enddo
c
      endif
c
      return
      end

c-----schmit-----------------------------------------------------------
c
      subroutine schmit(psiw,istadd,nhi)
*
***********************************************************************
*                                                                     *
*        schmit - subroutine to orthonormalize the set of             *
*                 states psi                                          *
*                 new states are replaced in psi                      *
*                 variables are as follows:                           *
*                                                                     *
*        if sym is 1.0 then npar (taken from fiels npar(ist) at       *
*        common /levels/) is a good quantum number and only those     *
*        states having the same z-parity need be orthogonalized.      *
*        if sym is 0.0 all states need be orthogonalized.             *
*        the routine runs over one 'm'-block of the states.           *
*        the wave functions are contained in 'psiw' running from      *
*        'n' = 1 through 'nhi'. the corresponding global quantum      *
*        number is computed with 'istadd' as base point:              *
*                      ist = istadd + n                               *
*                                                                     *
*        the norm of each component separately is returned on         *
*        'cmpn' at common /levels/; this serves as a diagnostic       *
*        tool.                                                        *
*                                                                     *
*        the orthonormalization should be done in order of            *
*        increasing energy. the s.p. energies on 'spenrg' are         *
*        first sorted in energy order. this is done on a              *
*        a workspace 'eord'. The ordering is stored in the            *
*        pointer 'isort'. 'eord' and 'spenrg' from the                *
*        common block /levels/.                                       *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'levels.inc'
!      include 'wavefr.inc'
c
c
      dimension psiw(0:nzd,0:nrd,nsd,nhi)
      dimension eord(nhwd)                    ! workspace
      dimension isort(nhwd)
      data iford/1/
      LOGICAL, PARAMETER :: tsortwf = .FALSE.
      REAL(8), ALLOCATABLE :: esav(:),spesav(:),occsav(:)
      INTEGER, ALLOCATABLE :: iparsav(:)
*
***********************************************************************
*
      if(ifstab.eq.1. .and. iford.ne.1)
     &  stop ' IFSTAB=1 needs energy ordering in SCHMIT. set IFORD=1'
*
***********************************************************************
*                                                                     *
*     energy ordering on workspace 'eord' to pointer 'isort'          *
*                                                                     *
***********************************************************************
*
      if(nhi.gt.nhwd)
     &  stop ' in SCHMIT: #states "nhi" > supply "nhwd" '
      do in=1,nhi
        eord(in)   = spenrg(in+istadd)
        isort(in) = in+istadd
      enddo
      if(iford.eq.1) then
      do n=1,nhi
        emin   = 1D32
        imin   = n
        do i=n,nhi
          if(eord(i).lt.emin) then
            emin   = eord(i)
            imin   = i
          endif
        enddo
        isav         = isort(imin)
        eord(imin)   = eord(n)
        isort(imin)  = isort(n)
        eord(n)      = emin
        isort(n)     = isav
      enddo
      endif
!
! optionally sort wavefunctions 
!
!      IF(tsortwf) THEN
!          itest = 0
!          DO in=1,nhi
!             IF(isort(in).NE.in+istadd) itest = 1+itest
!          END DO
!          IF(itest>0) THEN
!          ALLOCATE(spesav(nhwd),esav(nhwd),occsav(nhwd),iparsav(nhwd))
!!         reorder
!          DO in=1,nhi
!            ist = isort(in)
!            inw = ist-istadd
!            iparsav(in) = ipar(ist)
!            spesav(in) = spenrg(ist)
!            esav(in) = eord(inw)
!            occsav(in) = vsqf(ist)
!            dpsiz(:,:,:,in) = psiw(:,:,:,inw)
!          END DO
!!         copy back
!          DO in=1,nhi
!            ist = in+istadd
!            ipar(ist) = iparsav(in)
!            spenrg(ist) = spesav(in)
!            eord(in) = esav(in)
!            vsqf(ist) = occsav(in)
!            psiw(:,:,:,in) = dpsiz(:,:,:,in) 
!          END DO
!          DEALLOCATE(spesav,occsav,iparsav)
!        END IF
!      END IF
*
***********************************************************************
*                                                                     *
*     first loop (do 2000) runs over the state to be orthogonalised   *
*     to all the previous states and then to be normalised.           *
*     second loop (do 1000) runs over the previous states.            *
*                                                                     *
***********************************************************************
*
c      ist    = istadd
      do 2000 in=1,nhi
c      ist     = ist+1
      ist = isort(in)
      inw = ist-istadd
      npari   = npar(ist)
      in1    = in-1
      if(in1 .lt. 1) go to 3000
*
c      jst    = istadd
      do 1000 jn=1,in1
c      jst    = jst+1
      jst = isort(jn)
      jnw = jst-istadd
*
*        test if parity is a good quantum number select those states  *
*
      nparj   = npar(jst)
      if(nparj.ne.npari .and. sym.gt.zero) go to 1000
*
*     compute overlap
*
      ovlp   = wfovlp(psiw(0,0,1,inw),psiw(0,0,1,jnw))
*
*     extract psi(jn)-component in psi(in)
*
      do is=1,nsd
      do ir=0,nr
      do iz=0,nz
        psiw(iz,ir,is,inw)=psiw(iz,ir,is,inw)
     &                     -ovlp*psiw(iz,ir,is,jnw)
      enddo
      enddo
      enddo
*
1000  continue
*
*     normalise state psi(in), and add up partial norm.
*
3000  continue
      wnor   = wfovlp(psiw(0,0,1,inw),psiw(0,0,1,inw))
      if(wnor.lt.smal) then
        write(7,'(a,4i5,g12.4)')
     &   ' SCHMIT: ist,inw,istadd,nhi,wnor=',ist,inw,istadd,nhi,wnor
        stop ' too small norm in SCHMIT'
      endif
      wnor   = one/sqrt(wnor)
*
      do is=1,nsd
      do ir=0,nr
      do iz=0,nz
        psiw(iz,ir,is,inw) = psiw(iz,ir,is,inw)*wnor
      enddo
      enddo
      enddo
*
      sum2 = zero
      do is=1,nsd
       cmpn(is,ist) = wfovpc(psiw(0,0,is,inw))
       sum2         = cmpn(is,ist) + sum2
      enddo
      spn(ist) = sqrt(sum2)
*
2000  continue
*
      return
      end

c
c-----projec-----------------------------------------------------------
c
      subroutine projec(psiw,nparps,istadd,nhi,offnor)
*
***********************************************************************
*                                                                     *
*        projec - subroutine to project out of 'psiw' the states      *
*                 starting at 'istadd' and going up to 'nhi'.         *
*                 'nparps' is the z-parity of 'psiw'.                 *
*                 'offnor' is output of the squared sum of projected  *
*                 elements.                                           *
*                 This routine makes sense only for pure HF states!   *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'levels.inc'
      include 'wavefr.inc'
c
c
      dimension psiw(0:nzd,0:nrd,nsd)
*
***********************************************************************
*                                                                     *
*     first loop (do 2000) runs over the states to be projected.      *
*                                                                     *
***********************************************************************
*
      sumnor = 0.0
      ist    = istadd
      do 2000 in=1,nhi
      ist     = ist+1
      npari   = npar(ist)
*        test if parity is a good quantum number select those states  *
      if(nparps.eq.npari .or. sym.eq.zero) then
*
*     compute overlap
*
        ovlp   = wfovlp(psi(0,0,1,ist),psiw(0,0,1))
        sumnor = ovlp*ovlp + sumnor
c        write(6,'(a,i5,g12.4)') ' ist,ovlp=',ist,ovlp   !testprint
*
*     extract psi(in)-component in psiw
*
         do is=1,nsd
         do ir=0,nr
         do iz=0,nz
           psiw(iz,ir,is)=psiw(iz,ir,is)-ovlp*psi(iz,ir,is,ist)
         enddo
         enddo
         enddo
*
      endif
*
2000  continue
      offnor = sumnor
*
      return
      end
c
c-----ondiag-----------------------------------------------------------
c
      subroutine ondiag(psih,jst,nparps,istadd,nhi)
*
***********************************************************************
*                                                                     *
*        ondiag - subroutine to force the diagonalization amongst     *
*                 the occupied states.                                *
*                 Input h*psi for state nr. 'istadd' on 'psih'.       *
*                 Computes all h-overlaps up to 'nhi' and performs    *
*                 a perturbative correction towards diagonality.      *
*                 'jst' is the nr. of the state associated to 'psih'  *
*                 in absolute addres mode.                            *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'dynpar.inc'
      include 'levels.inc'
      include 'wavefr.inc'
c
c
      dimension psih(0:nzd,0:nrd,nsd)
*
      data stepmx/0.1/                ! maximum admixture of other states
      stepac = 0.5*stepmx/pi
*
***********************************************************************
*                                                                     *
*     first loop (do 2000) runs over the states to be projected.      *
*                                                                     *
***********************************************************************
*
      ist    = istadd
      do 2000 in=1,nhi
      ist     = ist+1
      npari   = npar(ist)
*        test if parity is a good quantum number select those states  *
      if((nparps.eq.npari .or. sym.eq.zero).and. ist.ne.jst) then
*
*     compute overlap
*
        ovlp   = wfovlp(psi(0,0,1,ist),psih(0,0,1))
        corfac = stepac*atan(ovlp/((spenrg(jst)-spenrg(ist))*stepac))
*
*     extract psi(in)-component in psih
*
         do is=1,nsd
         do ir=0,nr
         do iz=0,nz
           psi(iz,ir,is,jst)=psi(iz,ir,is,jst)
     &                       +corfac*psi(iz,ir,is,ist)
         enddo
         enddo
         enddo
*
      endif
*
2000  continue
*
      return
      end
c
c-----vsqtr -----------------------------------------------------------
c
      subroutine vsqtr(vsqf,wate,nst,eord,npiq,iorder)
*
***********************************************************************
*                                                                     *
*        subroutine to dynamically compute the occupation numbers     *
*        by simple filling in energy order from below.                *
*        The energy order is returned on 'iorder'.                    *
*        In case of pairing the 1 st guess to the occupation          *
*        is done here.                                                *
*                                                                     *
*        Input:                                                       *
*                  wate(*)  -  spin and ml degeneracy factors for     *
*                              the orbits                             *
*                  eord(*)  -  s.p. energies on which ordering is     *
*                              performed                              *
*                  nst      -  number of single particle states of    *
*                              charge q                               *
*                  npiq     -  number of particles of charge q        *
*                                                                     *
*        Output:                                                      *
*                  vsqf(*)  -  resulting occupation numbers           *
*                              for the filling approximation          *
*                  iorder(*)   ordering of the states by energy       *
*                              note the original order is by q and    *
*                              then by ml                             *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
c
      dimension vsqf(nst),wate(nst),eord(nst),iorder(nst),iorderinv(nst)
c
c     initialize pointer field with trivial ordering
c
      do j=1,nst
        iorder(j) = j
      enddo
c
c     determine energy ordering and fill 'vsqf'
c
      rpcur=0.0
      do 1 j=1,nst
      emin=1.0e12
      jbm = j
      do 2 jb=j,nst
      if(eord(iorder(jb)).lt.emin) then
        emin=eord(iorder(jb))
        jbm=jb
      endif
 2    continue
      iosav      = iorder(j)
      iorder(j)  = iorder(jbm)
      iorder(jbm)= iosav
      iosav      = iorder(j)
c      write(7,'(a,i5,a,60(1x,20i5/))')
c     &  ' j=',j,':',(iorder(i),i=1,nst)
      ve2=min(one,(npiq-rpcur)/wate(iosav))
      ve2=max(ve2,zero)
      rpcur=rpcur+wate(iosav)*ve2
      vsqf(iosav)=ve2
 1    continue
      return
      end


c-----sortenerg -----------------------------------------------------------
c
      subroutine sortenerg(nst,eord,iorder)
*
***********************************************************************
*                                                                     *
*        The energy order is returned on 'iorder'.                    *
*                                                                     *
*        Input:                                                       *
*                  eord(*)  -  s.p. energies on which ordering is     *
*                              performed                              *
*                  nst      -  number of single particle states of    *
*                              charge q                               *
*                                                                     *
*        Output:                                                      *
*                  iorder(*)   ordering of the states by energy       *
*                              note the original order is by q and    *
*                              then by ml                             *
*                                           
*                                                                     *
*        for given m (xk respectively) we do the damped gradient      *
*        step with the effective hamiltonian for the upper            *
*        components alone (sequence: lowcmp,efhmlt,dampng,do-loop).   *                          *
***********************************************************************
*
      include 'paramr.inc'
c
      dimension eord(nst),iorder(nst)
c
c     initialize pointer field with trivial ordering
c
      do j=1,nst
        iorder(j) = j
      enddo
      rpcur=0.0
      do j=1,nst
      emin=1.0e12
      jbm = j
      do jb=j,nst
      if(eord(iorder(jb)).lt.emin) then
        emin=eord(iorder(jb))
        jbm=jb
      endif
      enddo
        iosav      = iorder(j)
        iorder(j)  = iorder(jbm)
        iorder(jbm)= iosav
      enddo

      return
      end


