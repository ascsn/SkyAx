c-----zpe_cm -----------------------------------------------------------
c
      subroutine zpe_cm
c
c     zero point energy for c.m. motion
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
c
      data intest/0/
      data itstpr/1/
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
      IF(ifzpe<1) THEN
        zpecmz = 0D0
        widcmz = 0D0
        widdz  = 0D0
        RETURN
      END IF
c
c     1. double loop to accumulate x-y-width:
c
      if(intest.eq.1) write(7,'(a)') ' Nabla_xy:'
      sumwxy = zero
      sumexy = zero
      sumdir = zero
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         if(ist.le.nst(1)) then
            jstmax = nst(1)
            jstmin = 1
            eferm  = fer(1)
CR            efermo = ferolp
         else
            jstmax = nst(1)+nst(2)
            jstmin = nst(1)+1
            eferm  = fer(2)
CR            efermo = feroln
         endif
c
c     gradient operator Nabla+ on 'psip' and Nabla- on 'psin'
c
         call nablar_spin(psi(0,0,1,ist),psi1,mpsi(ist),1,one)
         
         do is=1,nsd
            mm        = mist+is-1
            npist  = nparz(is,ist)
            nparr     = 1+mod(mm,2)
cPK            call mulr(derr1(0,0,nparr),psi(0,0,is,ist),
cPK     &                psi1(0,0,is),one,0)
            ir = 0
            if(mm.eq.1) then
               do iz=0,nz
                  psip(iz,ir,is) = 0D0
                  psin(iz,ir,is) = 2D0*psi1(iz,ir,is)
               enddo
            elseif(mm.eq.-1) then
               do iz=0,nz
                  psin(iz,ir,is) = 0D0
                  psip(iz,ir,is) = 2D0*psi1(iz,ir,is)
               enddo
            else
               do iz=0,nz
                  psip(iz,ir,is) = 0D0 !psi1(iz,ir,is)
                  psin(iz,ir,is) = 0D0 !psip(iz,ir,is)
               enddo
            endif
            do ir=1,nr
               r    = (ir+sr)*dr
               rinv = mm/r
               do iz=0,nz
                  psip(iz,ir,is) = psi1(iz,ir,is)-rinv*psi(iz,ir,is,ist)
                  psin(iz,ir,is) = psi1(iz,ir,is)+rinv*psi(iz,ir,is,ist)
               enddo
            enddo
         enddo
c                                       auxiliary variables
         npist  = nparz(1,ist)
         vi    = sqrtb(vsqf(ist))
         ui    = sqrtb(one-vi*vi)
         IF(mist<0) ui=-ui
         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-eferm
         equasi = sqrt(edifi*edifi+deltai*deltai)
c                                       the overlaps
         tstovl = wfovlp(psip,psip)
         sumdir = 2.0*vi*vi*(tstovl+wfovlp(psin,psin)) +sumdir
         if(intest.eq.1) write(7,'(a,i3,1pg12.4)') 
     &        ' test: ist, <Nabla**2>=',ist,tstovl
         do jst=jstmin,jstmax
            if((npist.eq.nparz(1,jst) .and.
     &           (mist+1.eq.mpsi(jst) .or. mist-1.eq.mpsi(jst)) )
     &           .or.
     &           (npist.eq.nparz(2,jst) .and.
     &           ((mist.eq.0 .and. mpsi(jst).eq.0)
     &           .or.(mist.eq.-1 .and. mpsi(jst).eq.-1)))      ) then
               vj    = sqrtb(vsqf(jst))
               uj    = sqrtb(one-vj*vj)
               IF(mpsi(jst)<0) uj=-uj
               deltaj = deltaf(jst)*wtpair(jst)
               edifj  = spenrg(jst)-eferm
               equasj = sqrt(edifj*edifj+deltaj*deltaj)
               if(istsym.le.0 .and. mist.eq.0 .and. mpsi(jst).eq.0) then
                  pij    = wfovoc(psi(0,0,1,jst),psin(0,0,2))
     &                 + wfovoc(psi(0,0,2,jst),psin(0,0,1))
               elseif(mist+1.eq.mpsi(jst)) then
                  pij    = wfovlp(psi(0,0,1,jst),psip)
               elseif(mist-1.eq.mpsi(jst)) then
                  pij    = wfovlp(psi(0,0,1,jst),psin)
               else
                  pij    = 0D0
               endif
               sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)
               if(intest.eq.1)
     &              write(7,'(a,3(2i2,1x),2(1pg12.4))')
     &              ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
     &              ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij,
     &              pij*(vj*ui-uj*vi)
               if(itstpr.ne.0) 
     &              write(6,'(a,2i4,2f8.2,4(1pg12.4))')
     &              ' i,j,equasi/j,phase,pij=',
     &              ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
               pij    = pij*(vj*ui-uj*vi)
               pij    = pij*pij
               sumexy = pij/(equasi+equasj) + sumexy
               sumwxy  = pij + sumwxy
            endif
         enddo
      enddo
      IF(istsym>0) THEN
        sumdir = sumdir/2D0
        sumexy = sumexy/2D0
        sumwxy = sumwxy/2D0
      END IF
      widcmx = sumwxy/2.0
      widdxy = sumdir/2.0
      h2mav0 = 0.5D0*(h2m(1)+h2m(2))
      zpecmx = h2mav0*widcmx/(npart(1)+npart(2))
      write(7,'(a,2(1pg12.4))') 
     &     ' cm_xy: ZPE,Inglis ZPE=',zpecmx,sumwxy/sumexy/2.0
      write(7,'(a,2(1pg12.4))') '  dir.:      ZPE,width=',
     &     h2mav0*sumdir/(npart(1)+npart(2)),sumdir
c
c     2. loop to accumulate c.m. width along z
c
      if(intest.eq.1) write(7,'(a)') ' Nabla_z:'
      sumwz  = zero
      sumez  = zero  
      sumdir = zero
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         if(ist.le.nst(1)) then
            jstmax = nst(1)
            jstmin = 1
            eferm  = fer(1)
CR            efermo = ferolp
         else
            jstmax = nst(1)+nst(2)
            jstmin = nst(1)+1
            eferm  = fer(2)
CR            efermo = feroln
         endif
c
c     gradient operator Nabla_z on 'psip' 
c
         call nablaz_spin(psi(0,0,1,ist),psi1,nparz(1,ist),1,one)

         do is=1,nsd
cPK            mm        = mist+is-1
cPK            npist  = nparz(is,ist)
cPK            nparr     = 1+mod(mm,2)
cPK            call mulz(derz1(0,0,npist),psi(0,0,is,ist),wrk4,one,0)
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
         IF(mist<0) ui=-ui
         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-eferm
         equasi = sqrt(edifi*edifi+deltai*deltai)
c                                       the overlaps
         tstovl = wfovlp(psip,psip)
         if(intest.ne.0) write(7,'(a,i3,1pg12.4)') 
     &        ' test: ist, <Delta_z**2>=',ist,tstovl
         sumdir = 2.0*vi*vi*tstovl +sumdir
         do jst=jstmin,jstmax
            if(npist.eq.nparz(2,jst) .and. mist.eq.mpsi(jst)) then
               vj    = sqrtb(vsqf(jst))
               uj    = sqrtb(one-vj*vj)
               IF(mpsi(jst)<0) uj=-uj
               deltaj = deltaf(jst)*wtpair(jst)
               edifj  = spenrg(jst)-eferm
               equasj = sqrt(edifj*edifj+deltaj*deltaj)
               pij    = wfovlp(psi(0,0,1,jst),psip)
               if(intest.eq.1)   write(7,'(a,3(2i2,1x),2(1pg12.4))')
     &              ' ist,jst,nparzi,nparzj,mpsii,mpsij,pij=',
     &              ist,jst,npist,nparz(1,jst),mist,mpsi(jst),pij*pij,
     &              (pij*(vj*ui-uj*vi))**2
               if(itstpr.ne.0) 
     &              write(6,'(a,2i4,2f8.2,4(1pg12.4))')
     &              ' i,j,equasi/j,phase,pij=',
     &              ist,jst,equasi,equasj,(vj*ui-uj*vi),pij
               sumdir = sumdir - 2.0*pij*pij*vi*vj*(vi*vj+ui*uj)
               pij    = pij*(vj*ui-uj*vi)
               pij    = pij*pij
               sumez  = pij/(equasi+equasj) + sumez
               sumwz  = pij + sumwz
            endif
         enddo
      enddo
      IF(istsym>0) THEN
        sumdir = sumdir/2D0
        sumez = sumez/2D0
        sumwz = sumwz/2D0
      END IF
      zpecmz = h2mav0*sumwz/(npart(1)+npart(2))
      widcmz = sumwz
      widdz  = sumdir
      write(7,'(a,2(1pg12.4))') 
     &     '  cm_z: ZPE,Inglis ZPE=',zpecmz,sumwz/sumez/2.0
      write(7,'(a,2(1pg12.4))') 
     &     '  dir.:      ZPE,width=',h2mav0*sumdir/(npart(1)+npart(2)),
     &     sumdir

c     
      return
      end


c-----trizpe -----------------------------------------------------------
c
      subroutine trizpe
c
c     zero point energy for triaxial vibrations,
c     deduced simply from quadrupole momentum.
c
      include 'paramr.inc'
      include 'constr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      include 'derivs.inc'
c
      data intest/0/
c      data itstpr/0/
c
c     buffered SQRT as internal statement function
c
      sqrtb(x) = sqrt(max(x,smal))
c
c-----------------------------------------------------------------------
c
c
c     double loop to accumulate r**2 width with 'm' stepping by +-2
c
      if(intest.eq.1) write(7,'(a)') ' Nabla_xy:'
      sumQ   = zero
      sumQe  = zero
      sumQi  = zero
      sumQi2 = zero
      sumQi3 = zero
      do ist=1,nst(1)+nst(2)
         mist   = mpsi(ist)
         if(ist.le.nst(1)) then
            jstmax = nst(1)
            jstmin = 1
            eferm  = fer(1)
CR            efermo = ferolp
         else
            jstmax = nst(1)+nst(2)
            jstmin = nst(1)+1
            eferm  = fer(2)
CR            efermo = feroln
         endif
c
c     transition operator r**2 on 'psip'
c
         do is=1,nsd
            do ir=0,nr
               r    = (ir+sr)*dr
               r2   = r*r
               do iz=0,nz
                  psip(iz,ir,is) = r2*psi(iz,ir,is,ist)
               enddo
            enddo
         enddo
c                                       auxiliary variables
         npist  = nparz(1,ist)
         mistup = mist+2
         mistdw = abs(mist-2)
         vi    = sqrtb(vsqf(ist))
         ui    = sqrtb(one-vi*vi)
         deltai = deltaf(ist)*wtpair(ist)
         edifi  = spenrg(ist)-eferm
         equasi = sqrt(edifi*edifi+deltai*deltai)
c                                       the overlaps
         do jst=jstmin,jstmax
            if(npist.eq.nparz(1,jst) .and.
     &           (mistup.eq.mpsi(jst) .or. mistdw.eq.mpsi(jst)) ) then
               vj    = sqrtb(vsqf(jst))
               uj    = sqrtb(one-vj*vj)
               deltaj = deltaf(jst)*wtpair(jst)
               edifj  = spenrg(jst)-eferm
               equasj = sqrt(edifj*edifj+deltaj*deltaj)
c??               pij    = wfovoc(psi(0,0,1,jst),psip(0,0,2))
c??     &                 + wfovoc(psi(0,0,2,jst),psip(0,0,1))
               pij    = wfovoc(psi(0,0,1,jst),psip(0,0,1))
     &                 + wfovoc(psi(0,0,2,jst),psip(0,0,2))
               volQ   = (pij*pij)*(ui*vj+uj*vi)**2
               eqainv  = 1.0/(equasi+equasj)
               sumQ   = volQ                        + sumQ
               sumQe  = volQ*(equasi+equasj)        + sumQe
               sumQi  = volQ*eqainv                 + sumQi
               sumQi2 = volQ*eqainv*eqainv          + sumQi2
               sumQi3 = volQ*eqainv*eqainv*eqainv   + sumQi3
            endif
         enddo
      enddo
      zpetri = 0.5*sumQi2/sumQi3
      omgtri = sqrt(sumQi/sumQi3)
      write(7,'(a,5(1pg12.4))') 
     &  ' trizpe: sumQ,sumQe,sumQi,sumQi2,sumQi3=',
     &            sumQ,sumQe,sumQi,sumQi2,sumQi3
      write(7,'(a,3(1pg12.4))') 
     &     '      ZPE,omega,GR=',zpetri,omgtri,sqrt(sumQe/sumQi)
c     
      return
      end
