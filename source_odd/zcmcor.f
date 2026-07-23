c
c-----zcmobs-----------------------------------------------------------
c
      function zcmobs(rhoini)
c
c     computes the centre-of-mass for the density distribution 
c     given in two components on 'rhoini'.
c
      include 'paramr.inc'
      include 'fields.inc'
      include 'dynpar.inc'
      dimension rhoini(0:nzd,0:nrd,kqm)
c
c      data three,seven,eight,twfour /3.0d0,7.0,8.0,24.0/
c
c----------------------------------------------------------------------
c    
c     compute actual center of mass
c
      if(isym.ne.1) then
c
        drezw  = one/12.0D0
        d1     = zero
        at1    = zero
        do iz=1,nzm1
          z     = ((iz)+sz-icm)*dz
          rhot  = rhoini(iz,0,1)+rhoini(iz,0,2)
          at1   = rhot     + at1
          d1    = rhot*z   + d1 
        enddo
        at1   = drezw * at1
        d1    = drezw * d1 
c
        at2   = zero
        d2    = zero
        z     = ((sz)-icm) * dz
        do ir=1,nr
          rhot  = (rhoini(0,ir,1)+rhoini(0,ir,2))*((ir)+sr)
          at2   = rhot     + at2
          d2    = rhot*z   + d2
        enddo
        at2   = oh * at2
        d2    = oh * d2
c
        at3   = zero
        d3    = zero
        do ir=1,nr
          dir   = (ir)+sr
          do iz=1,nzm1
            z     = ((iz)+sz-icm)*dz
            rhot  = (rhoini(iz,ir,1)+rhoini(iz,ir,2))*dir
            at3   = rhot     + at3
            d3    = rhot*z   + d3
          enddo
        enddo
c
c     the different parts are summed up:
c
        vol   = (sym+one)*dvol
        atnum = (at1+at2+at3)*vol
        d     = (d1 +d2 +d3)*vol
c        write(7,'(a,i5,2(1pg12.4))') 
c     &     'ZCMOBS: icm,sz,zcm=',icm,sz,d/atnum
        zcmobs = d/atnum
      else
        zcmobs = zero
      endif
c
      return
      end

c
c-----corzcm------------------------------------------------------------
c
      subroutine corzcm(zcm)
c
c     shifts the wavefunctions such that the center-of-mass tends
c     to zero.
c     zcm  = c.m. at input, corrected c.m. at output
c
      include 'paramr.inc'
      include 'derivs.inc'
      include 'levels.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'wavefr.inc'
c
      data addzcm/-0.5D0/
      logical ttest
      data ttest/.false./
c
c
c----------------------------------------------------------------------
c
      write(6,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
      write(7,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
      stepfc = zcm*addzcm
      do iq=1,2
      do ir=0,nr
      do iz=0,nz
         rhoti(iz,ir,iq) = zero
      enddo
      enddo
      enddo
c
c     loop over wavefunctions in blocks of same quantum numbers
c
      do im=1,nom
c
        nhi=mhi(im)
        nlo=mlo(im)
        nadd=nhi-nlo+1
        kpsic=mim(im)
        iq=miq(im)
c
c       shift by 'stepfc'
c
        if(ttest) then
        do jst=nlo,nhi
           call nablaz_spin(psi(0,0,1,jst),psi1,nparz(1,jst),1,one)
           do is=1,2    
cPK          call mulz(derz1(0,0,nparz(is,jst)),psi(0,0,is,jst),
cPK     &                wrk4,one,0)
            do ir=0,nr   
            do iz=0,nz   
              psi(iz,ir,is,jst) = psi(iz,ir,is,jst) 
     &                           - stepfc*psi1(iz,ir,is)
            enddo
            enddo
          enddo
        enddo
c
c      Schmidt orthonormmalization of the states
c
        call schmit(psi(0,0,1,nlo),nlo-1,nadd)
        endif
c
c       new density
c
        do jst=nlo,nhi
          call sdenti(psi(0,0,1,jst),rhoti(0,0,iq),vsqf(jst),jst)
        enddo
c
      enddo
c
      zcm = zcmobs(rhoti)
      write(6,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
      write(7,'(a,1pg12.4)')  ' CORZCM with zcm=',zcm
c
      return
      end



c-----shiftcm------------------------------------------------------------
c
      subroutine shiftcm(shift,jst)
c
c     shifts the wavefunctions such that the center-of-mass tends
c     to zero.
c     zcm  = c.m. at input, corrected c.m. at output
c
      include 'paramr.inc'
      include 'derivs.inc'
      include 'levels.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'wavefr.inc'
c
      data     addzcm/-0.5D0/
c
c
c----------------------------------------------------------------------
c
      call nablaz_spin(psi(0,0,1,jst),psi1,nparz(1,jst),1,1.0D0)
c
      do is=1,2   
c        call mulz(derz1(0,0,nparz(is,jst)),psi(0,0,is,jst),
c     &            wrk4,one,0)
        do ir=0,nr  
          psi(0,ir,is,jst) = 0.0D0
        do iz=1,nz ! iz=0 should be zero 
          psi(iz,ir,is,jst) = psi(iz,ir,is,jst)
     &                       + shift*psi1(iz,ir,is)
        enddo
        enddo
      enddo
c
      return
      end
