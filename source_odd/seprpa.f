c-----seprpa-----------------------------------------------------------
c
      subroutine seprpa
c
c     separable RPA
c     -------   presently a testing code  --------------------
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
      if(kcrank.le.0) return
c
      call inie2q(0,1)                           ! 2QP energies on 'e2quas'
c
c     c.m. test
c
      if(isym.eq.0) then
        call ini_cmz(presp)
        call hcomm0(presp,ptest)
        pp   = opovlp(presp,presp)
        ph0p = opovlp(presp,ptest)
        call hcommP(presp,qresp)
        php   = opovlp(presp,qresp)
        write(6,'(a,5(1pg12.4))')
     &  ' CM TEST: <P*P>,<PHP>,<PH0P>=',pp,php,ph0p
      endif
c
c     test of 'hcommP' with sum rule model
c
      call quadms(qresp,0,1)                    ! initialize Q
c
      call hcommQ(qresp,ptest,presp,0)
      qq    = opovlp(qresp,qresp)
      qhq   = opovlp(qresp,presp)
c
      call hcommP(presp,qresp)
      pp    = opovlp(presp,presp)
      php   = opovlp(presp,qresp)
c
      call hcomm0(presp,ptest)
      ph0p  = opovlp(presp,ptest)
c
      call quadms(qresp,0,1)                    ! initialize Q
c
      call hcomm0(qresp,presp)
      qh0q  = opovlp(qresp,presp)
c
      call hcomm0(presp,ptest)
      qh02q = opovlp(presp,ptest)
c
      write(6,'(a,5(1pg12.4))')
     &  ' sumrule estimate: <Q*Q>,<QHQ>,<P*P>,<PHP>,omega=',
     &   qq,qhq,pp,php,sqrt(php/qhq)
      write(6,'(a,5(1pg12.4))')
     &  '                   <QH0Q>,<PH0P>,omega0=',
     &   qh0q,ph0p,sqrt(qh02q/qh0q)
c
      return
      end
c-----hcomm0-----------------------------------------------------------
c
      subroutine hcomm0(opin,opout)
c
c     Commutator [H_0,A] of mean-field Hamiltonian with A-response matrix,
c     matrix elements of A are entered on 'opin',
c     result is returned on 'opout'.
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
      dimension opin(ktriac),opout(ktriac)
c
c-----------------------------------------------------------------------
c
      ipostr = 0
      do ist=1,nst(1)+nst(2)
        mist   = mpsi(ist)
        npist  = nparz(1,ist)
        if(ist.le.nst(1)) then
          jstmax = nst(1)
        else
          jstmax = nst(1)+nst(2)
        endif
        do jst=ist,jstmax
          if(mist.eq.mpsi(jst) .and. npist.eq.nparz(1,jst)) then
            ipostr = 1 + ipostr
            opout(ipostr) = opin(ipostr)*e2quas(ipostr)
          endif
        enddo
      enddo
      write(6,'(a,i6)') ' HCOMM=: ipostr=',ipostr
c
      return
      end




c-----ini_cmz-----------------------------------------------------------
c
      subroutine ini_cmz(popout)
c
c     prepare c.m. momemtum on 'popout'
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
      dimension popout(ktriac)
c
c      data intest/0/
c      data itstpr/1/
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
      if(isym.ne.0) stop ' c.m. test runs only without z-symmetry'
c
      sumwz  = zero
      sumez  = zero  
      sumdir = zero
      ipostr = 0
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         if(ist.le.nst(1)) then
            jstmax = nst(1)
            eferm  = fer(1)
         else
            jstmax = nst(1)+nst(2)
            eferm  = fer(2)
         endif
c
c     gradient operator Nabla_z on 'psip' 
c
         call nablaz_spin(psi(0,0,1,ist),psi1,nparz(1,ist),1,one)
         do is=1,nsd
            mm        = mist+is-1
            npist  = nparz(is,ist)
            nparr     = 1+mod(mm,2)
cPK         call mulz(derz1(0,0,npist),psi(0,0,is,ist),wrk4,one,0)
            do ir=0,nr
               do iz=0,nz
                  psip(iz,ir,is) = psi1(iz,ir,is)
               enddo
            enddo
         enddo
c                                       auxiliary variables
         npist  = nparz(1,ist)
         vi    = sqrtb(vsqf(ist))
         ui    = sqrtb(one-vi*vi)
         do jst=ist,jstmax
            if(mist.eq.mpsi(jst)) then
               vj    = sqrtb(vsqf(jst))
               uj    = sqrtb(one-vj*vj)
               pij    = wfovlp(psi(0,0,1,jst),psip)
               pij    = pij*(vj*ui-uj*vi)
               ipostr = 1 + ipostr
               popout(ipostr) = pij
            endif
         enddo
      enddo
      write(6,'(a,i6)') ' C.M. Z-MODE initialized, IPOSTR=',ipostr
c     
      return
      end
