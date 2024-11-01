!------------------------------------------------------------------------------
! MODULE: steputil
!------------------------------------------------------------------------------
! DESCRIPTION: 
!> @brief
!!This module contains the routines to compute the mean-field Hamiltonian, 
!!applying it to a wave function, the Gram-Schmidt orthonormalization and other 
!!helper functions for the damped gradient step.
!------------------------------------------------------------------------------
MODULE steputil
IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------------  
! DESCRIPTION: hpsi
!> @brief
!!Action of the mean field hamiltonian on the wave function \c psiin.
!> @details
!!\f[ \hat h_q = U_q(\vec r)-\nabla\cdot\left[B_q(\vec r)\nabla\right]
!!+\I\vec W_q\cdot(\vec\sigma\times\nabla) \f]
!!Result is returned on \c psiout. At the same time \c hpsi evaluates the 
!!single particle energy and the variances of the state \c ist. \c psiin 
!!is supposed to be related to the state \c ist in the global ordering.
!>
!> @param[in] psiin
!> REAL(8), array, takes the input wave function.
!> @param[out] psiout
!> REAL(8), array, returns the Hamiltonian applied to the input wave 
!!function.
!> @param[in] ist
!> INTEGER, absolute adress of the state, to be used in all fields.
!--------------------------------------------------------------------------- 

SUBROUTINE hpsi(psiin,psiout,ist)

USE paramr
USE constr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid, ONLY: nr,nz,sr,laplace_spin,nablar_spin,nablaz_spin
IMPLICIT NONE

REAL(8), INTENT(IN)                      :: psiin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)                     :: psiout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)                      :: ist

INTEGER,PARAMETER :: ifclean=1

REAL(8),ALLOCATABLE :: psi1(:,:,:),psi2(:,:,:),psi3(:,:,:)

REAL(8) :: ebin,eferm2,esub,r,rinv,sigis,spefac,spenrg2,xmm
INTEGER :: iq,ir,is,isi,mm,nparr

ALLOCATE(psi1(0:nz,0:nr,nsd),psi2(0:nz,0:nr,nsd),psi3(0:nz,0:nr,nsd))
psi1=0D0;psi2=0D0;psi3=0D0

iq      = icharg(ist)
ebin    = spenrg(ist)

!     first, zero the output wavfunction.
psiout(0:nz,0:nr,1:2) = zero

CALL laplace_spin(psiin(0:nz,0:nr,1:2),psi1(0:nz,0:nr,1:2),nparz(1,ist),mpsi(ist),1,one)
CALL nablar_spin(psiin(0:nz,0:nr,1:2),psi2(0:nz,0:nr,1:2),   mpsi(ist),1,one)
CALL nablaz_spin(psiin(0:nz,0:nr,1:2),psi3(0:nz,0:nr,1:2),nparz(1,ist),1,one)

DO  is=1,2
  isi       = 3-is
  mm        = mpsi(ist)+is-1
  nparr     = 1+MOD(mm,2)
  sigis     = MOD(one*(is),two)*two - one
  
!Determine spefac, the factor for stabilazation of positive energy states.  
  IF(tstab) THEN
    eferm2 = fer(iq)*oh
    IF(ebin > zero) THEN
      spefac = one
    ELSE IF(ebin > eferm2) THEN
      spefac = one-ebin/eferm2
    ELSE
      spefac = zero
    END IF
  ELSE
    spefac = zero
  END IF
!Calculate potential and effective mass part.
  IF(spefac <= zero) THEN
    psiout(0:nz,0:nr,is) =  &
        (pot(0:nz,0:nr,iq)-confield(0:nz,0:nr)-ebin)*psiin(0:nz,0:nr,is)  &
        - bmass(0:nz,0:nr,iq)*psi1(0:nz,0:nr,is) &
        - drbmas(0:nz,0:nr,iq)*psi2(0:nz,0:nr,is)  &
        - dzbmas(0:nz,0:nr,iq)*psi3(0:nz,0:nr,is) + psiout(0:nz,0:nr,is)
  ELSE
    psiout(0:nz,0:nr,is) = (pot(0:nz,0:nr,iq)+spefac*potadd(0:nz,0:nr,iq)  &
        -confield(0:nz,0:nr)-ebin)*psiin(0:nz,0:nr,is)  &
        - bmass(0:nz,0:nr,iq)*psi1(0:nz,0:nr,is) &
        - drbmas(0:nz,0:nr,iq)*psi2(0:nz,0:nr,is)  &
        - dzbmas(0:nz,0:nr,iq)*psi3(0:nz,0:nr,is) + psiout(0:nz,0:nr,is)
  END IF
  IF(t4 /= zero) THEN
    
!   the spin-orbit part. action of the s.o. force is computed simply from
!   W_ls |psi> ~ (nabla rho)*sigma sigma*nabla |psi> - (nabla rho) * nabla |psi>

!   r=0 requires l'Hospital for 1/r*psiin
    IF(mm == 1) THEN
      psiout(0:nz,0,is) =  &
          +sigis*potls(0:nz,0,iq,1)*psi2(0:nz,0,is)   & 
          + psiout(0:nz,0,is)   ! 1 = r   2 = z  
      psiout(0:nz,0,isi)= -sigis*potls(0:nz,0,iq,1)*psi3(0:nz,0,is)  &
          -potls(0:nz,0,iq,2)*(-sigis*psi2(0:nz,0,is) +psi2(0:nz,0,is) )  &
          + psiout(0:nz,0,isi)
    END IF
!   now all other r
    DO ir=1,nr
      r   = (sr+ir)*dr
      rinv = one/r
      xmm  = mm*rinv
      psiout(0:nz,ir,is) = +sigis*xmm*potls(0:nz,ir,iq,1)*psiin(0:nz,ir,is)  &
          + psiout(0:nz,ir,is)
      psiout(0:nz,ir,isi)= -sigis*potls(0:nz,ir,iq,1)*psi3(0:nz,ir,is)  &
          -potls(0:nz,ir,iq,2)*(-sigis*psi2(0:nz,ir,is) +xmm*psiin(0:nz,ir,is) )  &
          + psiout(0:nz,ir,isi)
    END DO
  END IF
  
END DO

!     compute expectation value and variance

esub        = wfovlp(psiin(0:nz,0:nr,1:2),psiout(0:nz,0:nr,1:2))
spenrg2     = wfovlp(psiout(0:nz,0:nr,1:2),psiout(0:nz,0:nr,1:2))
spevar(ist) = SQRT(ABS(spenrg2-esub*esub))
spenrg(ist) = esub+ebin

DEALLOCATE(psi1,psi2,psi3)
RETURN
END SUBROUTINE hpsi

!---------------------------------------------------------------------------  
! DESCRIPTION: mfhmlt
!> @brief
!!Compute the mean field hamiltonian.
!> @details
!!First the 
!!routine computes the mean-field potential, the inverse effective mass, the 
!!spin-orbit potential, and then the Coulomb potential. An octupule kick can 
!!be added to explore asymmetric shapes or an isovector dipole potential. The 
!!derivatives of the inverse mass are determined and finally prints and test 
!!plots are generated.
!!
!> @param[in] iter
!> INTEGER, takes the iteration number.
!--------------------------------------------------------------------------- 

SUBROUTINE mfhmlt(iter)

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE fyfunc
USE grid, ONLY: nr,nz,f1int,icm,tsym,nr_testp,nz_testp,sr,sz,&
                nablaz_scalar,nablar_scalar,laplace_scalar,test_print
USE coulomb
IMPLICIT NONE


INTEGER, INTENT(IN)                  :: iter
!LOGICAL :: trepeat=.false.

LOGICAL :: talpinv,tdens2                                     

!     parameters for additional octupole kick in a few first iterations
REAL(8),PARAMETER :: octkick=0.7464D0,octampl=6D0

INTEGER,PARAMETER :: itmino=3,itmaxo=43

LOGICAL,SAVE :: tfirstcall=.true.,tfirstround=.true.

REAL(8) :: accn,accp,addpt,adparp,alpeff,d1pr,&
           d1pz,d1sr,d1sz,d2p,d2s,facne,facpr,parfac,potn,potold,&
           potp,r,r2,redpot,rlsrp,rlszp,rm,rnr,rp,rpr,rs,rsfac,rsfac2,&
           rsinv,rspow,rspow7,ssp,sss,tp,ts,z
INTEGER :: ic,iq,ir,iz
INTEGER, SAVE :: iterold

REAL(8),ALLOCATABLE :: wrk1(:,:),wrk2(:,:)

ALLOCATE(wrk1(0:nz,0:nr),wrk2(0:nz,0:nr))
wrk1=0D0;wrk2=0D0  !;psip=0D0;psin=0D0

IF(tfirstcall .AND. tfirstround) THEN
  iterold = iter-1
END IF

IF(iterold >= iter) tfirstround = .false.

iterold = iter
tfirstcall = .false.
tpair = (ipair == 7 .OR. ipair == 11) .AND. tpairon
IF(ABS(b7) > 1D-10) THEN
  tdens2 = .true.
ELSE
  tdens2 = .false.
END IF
IF(ABS(alpinv) > 1D-10) THEN
  talpinv = .true.
ELSE
  talpinv = .false.
END IF                                             

IF(iter == 0) THEN
  redpot = 1.0
  potold = 0.0
ELSE
  redpot = adypot
  potold = one-adypot
END IF


!***************************************************************  
!
!   Main part for the Skyrme or Fayans hamiltonian:
!
!***************************************************************

  IF(tFayans) THEN
    CALL Fayans_potentials()  
    CALL fy_pairpot()
    CALL fyfunc_grad()
!    pot(0,0,:)=(pot(0,1,:)+pot(1,0,:))/2D0
!    IF(MOD(iter,iprint)==0 .AND. itest>0)  CALL test_print(pot(0,0,1),0,nz_testp,0,nr_testp,'potp0')
  ELSE
    CALL Skyrme_potentials()  
  END IF

!***************************************************************   
!     add Coulomb potential
!***************************************************************   
  IF (chg > 0.0) THEN
    CALL coulmb(iter)
    DO ir=0,nr
      DO iz=0,nz
        wrk1(iz,ir) = crho(iz,ir)*coulpt(iz,ir)          ! ?? test, remove
        pot(iz,ir,1) = redpot*coulpt(iz,ir)+pot(iz,ir,1)
        wrk2(iz,ir) = MAX(rho(iz,ir,1),smal)**1.33333333333D0
      END DO
    END DO
    ecoul  = f1int(wrk1(0:nz,0:nr))
    ecolex = f1int(wrk2(0:nz,0:nr))*excfac*0.75D0
    IF(MOD(iter,iprint) == 0 .AND. itest >= 1) THEN
      CALL test_print(coulpt(0,0),0,nz_testp,0,nr_testp,'coulpt')
    END IF
  END IF
!  WRITE(7,'(a,2(1pg13.5))') 'COULPT(0,0)=',coulpt(0,0),redpot
!  IF(MOD(iter,iprint)==0 .AND. itest>0)  CALL test_print(pot(0,0,1),0,nz_testp,0,nr_testp,'potp1')

!***************************************************************   
!     optionally add octupole kick to explore non-symmetric shapes
!***************************************************************  
  IF(.NOT.tsym .AND. octampl > 1D-10 .AND. ABS(beta3) < 4.0D-2) THEN

    IF(iter > itmino .AND. iter < itmaxo) THEN

      WRITE(6,*) ' octupole kick active '
      DO ir=0,nr
        r2 = (ir+sr)*dr
        r2 = r2*r2
        DO iz=0,nz
          z = (iz+sz-icm)*dz - zcm
          addpt = (octampl*octkick)*ATAN(z*((z*z)-1.5D0*r2))
          pot(iz,ir,1) = addpt+pot(iz,ir,1)
          pot(iz,ir,2) = addpt+pot(iz,ir,2)
        END DO
      END DO
    END IF
  END IF

!***************************************************************   
!     add isovector dipole if asked for
!***************************************************************    
  IF(.NOT.tsym .AND. dipin /= 0D0) THEN
    anumb  = npart(1)+npart(2)
    facpr = npart(2)/anumb * redpot * dipin
    facne = -npart(1)/anumb * redpot * dipin
    DO iz=0,nz
      z = (iz+sz-icm)*dz - zcm
      pot(iz,0:nr,1) = facpr*z+pot(iz,0:nr,1)
      pot(iz,0:nr,2) = facne*z+pot(iz,0:nr,2)
    END DO
  END IF
  
!***************************************************************   
!     adjust potential at iz=0 in case of tsym being false
!***************************************************************   
  IF (.NOT.tsym) THEN
    pot(0,0:nr,1:2)=0.0D0
    potls(0,0:nr,1:2,1:2)=0.0D0
  END IF
  
!***************************************************************   
!     compute derivatives of inverse mass
!***************************************************************   
  wrk1(0:nz,0:nr) = bmass(0:nz,0:nr,1)- h2meff(1)
  wrk2(0:nz,0:nr) = bmass(0:nz,0:nr,2)- h2meff(2)
  
!  CALL test_print(wrk1(0:nz,0:nr),0,nz_testp,0,nr_testp,'BMASS_p')
!  CALL test_print(wrk2(0:nz,0:nr),0,nz_testp,0,nr_testp,'BMASS_n')
  IF(SUM(ABS(wrk1(0:nz,0:nr)))+SUM(ABS(wrk2(0:nz,0:nr)))> 0D0) THEN
  CALL nablaz_scalar(wrk1(0:nz,0:nr),dzbmas(0:nz,0:nr,1), 1,1,one)
  CALL nablaz_scalar(wrk2(0:nz,0:nr),dzbmas(0:nz,0:nr,2), 1,1,one)
  CALL nablar_scalar(wrk1(0:nz,0:nr),drbmas(0:nz,0:nr,1), 0,1,one)
  CALL nablar_scalar(wrk2(0:nz,0:nr),drbmas(0:nz,0:nr,2), 0,1,one)
  ELSE
    dzbmas(0:nz,0:nr,:)=0D0
    drbmas(0:nz,0:nr,:)=0D0
  END IF
!  IF(MOD(iter,iprint)==0 .AND. itest>0)  CALL test_print(pot(0,0,1),0,nz_testp,0,nr_testp,'potp3')
  
!***************************************************************    
!     prints and plots
!***************************************************************  
  IF (ipri_fields>0) THEN
  IF (MOD(iter,ipri_fields) == 0.AND.iter > 0) THEN
    WRITE(18,'(3(a,i3),a,f10.5,4(/a),a)')  &
        ' Z=',npart(1),' N=',npart(2),' ITER=',iter,' BETA2=',beta2,  &
        'X:z [fm]', 'Y:fields',  &
        'S: DIST(4.0) SCAL(0.8) THICKN(0.6) XMIN(0.0) ',  &
        '   z     U_proton    U_neutron     B_proton    B_neutron ',  &
        '      rho_p       rho_n       int(r*rho_p)  int(r*rho_n)' 
    DO iz=0,nz
      accp=0D0
      accn=0D0
      DO ir=1,nr
        r   = ((ir)+sr)*dr
        accp = rho(iz,ir,1)*r + accp
        accn = rho(iz,ir,2)*r + accn
      END DO
      WRITE(18,'(1x,f6.2,8(1pg13.5))')  &
          (iz+sz-icm)*dz,pot(iz,0,1),pot(iz,0,2),  &
          bmass(iz,0,1),bmass(iz,0,2),rho(iz,0,1),rho(iz,0,2), accp*dr,accn*dr
    END DO
    WRITE(18,'(1x/1x/4(/a),a)') 'X:r [fm]',  &
        'Y:fields', 'S: DIST(4.0) SCAL(0.8) THICKN(0.6) XMIN(0.0) ',  &
        '   r     U_proton   U_neutron    B_proton   B_neutron ',  &
        '     rho_p      rho_n'
    DO ir=0,nr
      WRITE(18,'(1x,f6.2,6g12.4)') (ir+sr)*dr,pot(0,ir,1),pot(0,ir,2),  &
          bmass(0,ir,1),bmass(0,ir,2),rho(0,ir,1),rho(0,ir,2)
    END DO
    WRITE(18,'(1x/1x)')
    IF(tstab) THEN
      WRITE(18,'(4(/a))') 'X: r (fm) ',  &
          'Y: potadd (1/fm**3) ', 'S: SCAL(0.9) THICKN(0.6) XMIN(0.0) ',  &
          'H: X   Y1,L   Y2,P  '
      WRITE(18,'(1x,f6.3,2g12.4)') ((ir+sr)*dr,  &
          potadd(0,ir,1),potadd(0,ir,2), ir=0,nr)
      WRITE(18,'(4(/a))') 'X: z (fm) ',  &
          'Y: potadd (1/fm**3) ', 'S: SCAL(0.9) THICKN(0.6) XMIN(0.0) ',  &
          'H: X   Y1,L   Y2,P  '
      WRITE(18,'(1x,f6.3,2g12.4)') ((iz+sz-icm)*dz,  &
          potadd(iz,0,1),potadd(iz,0,2), iz=0,nz)
    END IF
  END IF
  END IF
  
!***************************************************************   
!     prints at test level
!***************************************************************   
  IF(itest >= 1 .AND. MOD(iter,iprint) == 0) THEN
    WRITE(7,'(a,i3,a)') ' potentials in mfhmlt (iter=',iter,'):'
    CALL test_print(pot(0:nz,0:nr,1),0,nz_testp,0,nr_testp,'  potp')
    CALL test_print(pot(0:nz,0:nr,2),0,nz_testp,0,nr_testp,'  potn')
    CALL test_print(bmass(0:nz,0:nr,1),0,nz_testp,0,nr_testp,'  bmassp')
    CALL test_print(bmass(0:nz,0:nr,2),0,nz_testp,0,nr_testp,'  bmassn')
    CALL test_print(drbmas(0:nz,0:nr,1),0,nz_testp,0,nr_testp,' drbmasp')
    CALL test_print(drbmas(0:nz,0:nr,2),0,nz_testp,0,nr_testp,' drbmasn')
    IF(t4 /= zero) THEN
      CALL test_print(potls(0:nz,0:nr,1,1),0,nz_testp,0,nr_testp,'potls_r')
      CALL test_print(potls(0:nz,0:nr,1,2),0,nz_testp,0,nr_testp,'potls_z')
    END IF
  END IF  

IF(MOD(iter,iprint) == 10000) THEN
  CALL test_print(pot(0:nz,0:nr,1),0,nz,0,nr,'potp')
  CALL test_print(bmass(0:nz,0:nr,1),0,nz,0,nr,'bmassp')
  CALL test_print(drbmas(0:nz,0:nr,1),0,nz,0,nr,'drbmasp')
  CALL test_print(dzbmas(0:nz,0:nr,1),0,nz,0,nr,'dzbmasp')
  CALL test_print(rho(0:nz,0:nr,1),0,nz,0,nr,'rhop')
  CALL test_print(tau(0:nz,0:nr,1),0,nz,0,nr,'taup')
  CALL test_print(d2rho(0:nz,0:nr,1),0,nz,0,nr,'d2rhop')
  CALL test_print(pot(0:nz,0:nr,2),0,nz,0,nr,'potn')
  CALL test_print(bmass(0:nz,0:nr,2),0,nz,0,nr,'bmassn')
  CALL test_print(drbmas(0:nz,0:nr,2),0,nz,0,nr,'drbmasn')
  CALL test_print(dzbmas(0:nz,0:nr,2),0,nz,0,nr,'dzbmasn')
  CALL test_print(rho(0:nz,0:nr,2),0,nz,0,nr,'rhon')
  CALL test_print(tau(0:nz,0:nr,2),0,nz,0,nr,'taun')
  CALL test_print(d2rho(0:nz,0:nr,2),0,nz,0,nr,'d2rhon')
END IF

DEALLOCATE(wrk1,wrk2)

RETURN
CONTAINS

SUBROUTINE Skyrme_potentials()

REAL(8),ALLOCATABLE :: psip(:,:,:),psin(:,:,:)

ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd))

!***************************************************************   
!     first and secondderivatives on 'psip' and 'psin'
!***************************************************************   
  CALL laplace_scalar(rho(0:nz,0:nr,1),d2rho(0:nz,0:nr,1),+1,0,1,one)
  CALL laplace_scalar(rho(0:nz,0:nr,2),d2rho(0:nz,0:nr,2),+1,0,1,one)

  IF(t4 /= zero) THEN
    CALL nablaz_scalar(rho(0:nz,0:nr,1),psin(0:nz,0:nr,1), 1,1,one)
    CALL nablaz_scalar(rho(0:nz,0:nr,2),psin(0:nz,0:nr,2), 1,1,one)
    CALL nablar_scalar(rho(0:nz,0:nr,1),psip(0:nz,0:nr,1), 0,1,one)
    CALL nablar_scalar(rho(0:nz,0:nr,2),psip(0:nz,0:nr,2), 0,1,one)
  END IF
  
  IF(tpair) parfac = 1D0/rho0pr         ! factor on pairing potential
  
!***************************************************************   
!     compose mean-field potentials
!***************************************************************   
  DO  iq=1,2            
    ic      = 3-iq
    DO  ir=0,nr
      DO  iz=0,nz
        rp      = rho(iz,ir,iq)
        rm      = rho(iz,ir,ic)
        rs      = MAX(rm+rp,smal)
        rspow   = (rm+rp)**power
        tp      = tau(iz,ir,iq)
        ts      = tau(iz,ir,ic)+tp
        d2p     = d2rho(iz,ir,iq)
        d2s     = d2rho(iz,ir,ic)+d2p
        ssp     = drhols(iz,ir,iq)
        sss     = drhols(iz,ir,ic)+ssp
        
        pot(iz,ir,iq) = redpot* ( t0a*rs-t0b*rp  +t1p2*ts-t1m2*tp  &
            -t13m2*d2s+t13p2*d2p -b4*sss-b4p*ssp) + potold*pot(iz,ir,iq)
        IF(talpinv.AND.(iq == 2)) THEN                         !new
          rpr = rm
          rnr = rp
          alpeff = alpinv+betinv
          rsinv = 1D0/(1D0+alpeff*rs)
          rsfac = (rs+betinv*(rs*rs))*rsinv*alpeff
          rsfac2 = (1D0+betinv*(rs+rs)+(alpeff*betinv)*rs*rs)*  &
              (rsinv*rsinv)*alpeff
          potp = rsfac*(t3a*rs-t3b*rpr)  &
              +0.5D0*rsfac2*(t3a*(rs*rs)-t3b*(rpr*rpr+rnr*rnr))
          potn = rsfac*(t3a*rs-t3b*rnr)  &
              +0.5D0*rsfac2*(t3a*(rs*rs)-t3b*(rpr*rpr+rnr*rnr))
          pot(iz,ir,1) = redpot*potp + pot(iz,ir,1)
          pot(iz,ir,2) = redpot*potn + pot(iz,ir,2)
        ELSE IF(.NOT.talpinv) THEN
          rspow   = rs**power
          pot(iz,ir,iq) = redpot*rspow*(t3a*rs-t3b*rp +t3as*rs  &
              -t3bs*(rp*rp+rm*rm)/rs ) + pot(iz,ir,iq)
        END IF
        IF(tdens2.AND.(iq == 2)) THEN
          rpr = rm
          rnr = rp
          rspow7 = rs**pow7
          potp =  rspow7*(b7a*rs - b7pa*rpr  &
              -b7pb*(rpr*rpr+rnr*rnr)/MAX(rs,xsmall) )
          potn =  rspow7*(b7a*rs - b7pa*rnr  &
              -b7pb*(rpr*rpr+rnr*rnr)/MAX(rs,xsmall) )
          pot(iz,ir,1) = redpot*potp + pot(iz,ir,1)
          pot(iz,ir,2) = redpot*potn + pot(iz,ir,2)
        END IF                                                  !new
        IF(iq == 1 .AND. chgx /= zero)  &
            pot(iz,ir,1) = redpot*excfac*MAX(rp,smal)**edr + pot(iz,ir,1)
        IF(tpair .AND. iq == 2) THEN
          adparp = redpot*(pairfc(1)*(paircr(iz,ir,1)*paircr(iz,ir,1))  &
              +pairfc(2)*(paircr(iz,ir,2)*paircr(iz,ir,2)))*parfac
          IF(tregpair) THEN
            pot(iz,ir,1) = adparp*regpairadd(iz,ir,1) + pot(iz,ir,1)
            pot(iz,ir,2) = adparp*regpairadd(iz,ir,2) + pot(iz,ir,2)
          ELSE
            pot(iz,ir,1) = adparp + pot(iz,ir,1)
            pot(iz,ir,2) = adparp + pot(iz,ir,2)
          END IF
        END IF
        
!***************************************************************   
!     compose inverse mass
!***************************************************************  
        bmass(iz,ir,iq) = redpot*(h2meff(iq) +t1p2*rs-t1m2*rp)  &
            + potold*bmass(iz,ir,iq)
            
!***************************************************************   
!     compose spin-orbit potential
!***************************************************************         
        IF(t4 /= zero) THEN
          d1pr    = psip(iz,ir,iq)
          d1sr    = psip(iz,ir,ic)+d1pr
          d1pz    = psin(iz,ir,iq)
          d1sz    = psin(iz,ir,ic)+d1pz
          rlsrp   = rholsr(iz,ir,iq)
          rlszp   = rholsz(iz,ir,iq)
          potls(iz,ir,iq,1) = redpot*(b4*d1sr+b4p*d1pr+  &
              t1m20*rlsrp-t1p2x*(rlsrp+rholsr(iz,ir,ic)))  &
              + potold*potls(iz,ir,iq,1)
          potls(iz,ir,iq,2) = redpot*(b4*d1sz+b4p*d1pz+  &
              t1m20*rlszp-t1p2x*(rlszp+rholsz(iz,ir,ic)))  &
              + potold*potls(iz,ir,iq,2)
        END IF
        
      END DO
    END DO
  END DO

DEALLOCATE(psip,psin)

END SUBROUTINE Skyrme_potentials

SUBROUTINE Fayans_potentials()

REAL(8),ALLOCATABLE :: psip(:,:,:),psin(:,:,:),workenerg(:,:)
REAL(8) :: energd

IF(ABS(c0tens)>0D0 .OR. ABS(c1tens)>0D0) STOP 'tensor spin-orbit not yet implemented'

ALLOCATE(psip(0:nz,0:nr,nsd),psin(0:nz,0:nr,nsd),workenerg(0:nz,0:nr))

!***************************************************************   
!     first and secondderivatives on 'psip' and 'psin'
!***************************************************************   
  CALL laplace_scalar(rho(0:nz,0:nr,1),d2rho(0:nz,0:nr,1),+1,0,1,one)
  CALL laplace_scalar(rho(0:nz,0:nr,2),d2rho(0:nz,0:nr,2),+1,0,1,one)

  CALL nablaz_scalar(rho(0:nz,0:nr,1),psin(0:nz,0:nr,1), 1,1,one)
  CALL nablaz_scalar(rho(0:nz,0:nr,2),psin(0:nz,0:nr,2), 1,1,one)
  CALL nablar_scalar(rho(0:nz,0:nr,1),psip(0:nz,0:nr,1), 0,1,one)
  CALL nablar_scalar(rho(0:nz,0:nr,2),psip(0:nz,0:nr,2), 0,1,one)
  
  IF(tpair) parfac = 1D0/rho0pr         ! factor on pairing potential
  
!***************************************************************   
!     compose mean-field potentials
!***************************************************************   
!  pot(0,0,:) =0D0
  iq=1
  ic=3-iq
  DO  ir=0,nr
    DO  iz=0,nz
      rp      = rho(iz,ir,iq)
      rm      = rho(iz,ir,ic)
      rs      = MAX(rm+rp,smal)
      rspow   = (rm+rp)**power
      tp      = tau(iz,ir,iq)
      ts      = tau(iz,ir,ic)+tp
      d2p     = d2rho(iz,ir,iq)
      d2s     = d2rho(iz,ir,ic)+d2p
      ssp     = drhols(iz,ir,iq)
      sss     = drhols(iz,ir,ic)+ssp
      
      CALL fyfunc_local(rho(iz,ir,1),rho(iz,ir,2),crho(iz,ir),  &
                psip(iz,ir,1),psin(iz,ir,1),psip(iz,ir,2),psin(iz,ir,2),  &
                d2rho(iz,ir,1),d2rho(iz,ir,2),  &
                tau(iz,ir,1),tau(iz,ir,2),  &
                0D0,0D0,drhols(iz,ir,1),drhols(iz,ir,2),  &
                pot(iz,ir,1),pot(iz,ir,2),bmass(iz,ir,1),bmass(iz,ir,2),  &
                potls(iz,ir,1,1),potls(iz,ir,1,2), &
                potls(iz,ir,2,1),potls(iz,ir,2,2), &
                energd)
      
      workenerg(iz,ir) = energd
    END DO
  END DO
  bmass(0:nz,0:nr,1) = bmass(0:nz,0:nr,1) + h2meff(1)
  bmass(0:nz,0:nr,2) = bmass(0:nz,0:nr,2) + h2meff(2)
  energyFy = f1int(workenerg)

  DEALLOCATE(psip,psin,workenerg)

END SUBROUTINE Fayans_potentials

END SUBROUTINE mfhmlt

!---------------------------------------------------------------------------  
! DESCRIPTION: stapot 
!> @brief
!!Computes stabilizing potential for positive-energy states.
!!The resulting potential is returned on array 'potadd' via 
!! module 'fields'.
!> @details
!!A deformation dependent "coordinate" is deduced from the density of surface 
!!states, accumulated on \c potadd as
!!\f[ \texttt{effrad}  =  \max(0,-\log(\texttt{potadd}/\texttt{rhocri})) \f]
!!\f[ \texttt{rhocri}  =  \texttt{rho0}\cdot\texttt{rhocut} \f]
!!Where \c rhocri (recommended 0.08 or 0.07) sets the fractional density at 
!!the outer surface. The slope of 'effrad' is then scaled to 1 and the 
!!desired slope \c stslop of the additional potential is reconstructed.
!!
!!The additional potential is switched with the particle energy as
!!\f[\theta_w=\begin{array}{lr} 
!!  0 & \mbox{for }\epsilon < -\texttt{wext}\\
!!  (1-|\epsilon|)/\texttt{wext} & \mbox{for }-\texttt{wext} < \epsilon < 0\\
!!  1 & \mbox{for }0 < \epsilon 
!!\end{array} \f]
!!The energy-width \c wext is choosen as \f$ \texttt{wext} = \epsilon_F/2 \f$
!>
!> @param[out] delrm
!> REAL(8), array, returns for a point (z,r) the closest distance
!! to the half-density surface. That is an auxiliary array for
!! computation of the stabilizing potential 'potadd'.
!--------------------------------------------------------------------------- 

SUBROUTINE stapot(delr)

USE paramr
USE constr
USE dynpar
USE fields
USE observ
USE wavefr
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(OUT)   :: delr(0:nz,0:nr)

REAL(8) :: rad2ad,rhoact,rhocri,rhoprt,rhotrd,slope,slpfac,xnorm1
INTEGER :: iq,ir,iz,iz1
REAL(8),ALLOCATABLE :: wrk1(:,:)

ALLOCATE(wrk1(0:nz,0:nr))
wrk1=0D0
DO  iq=1,2
  
!     normalize stabilizing density
  xnorm1 = one/f1int(rhosur(0:nz,0:nr,iq))
  DO ir=0,nr
    DO iz=0,nz
      delr(iz,ir) = xnorm1*rhosur(iz,ir,iq)
      wrk1(iz,ir) = (z(iz)**2+r(ir)**2)*delr(iz,ir)
    END DO
  END DO
  rad2ad = 1.666666667*f1int(wrk1(0:nz,0:nr))
  
!     prepare parameters

  rhoprt = 0.238732414*rad2ad**(-0.66666667)
  rhocri = rhocut*rhoprt            ! density defining surface
  rhotrd = rhoprt/3.0               ! density defining interior
  
!     compute unscaled effective radius
  DO ir=0,nr
    DO iz=0,nz
      rhoact = MAX(delr(iz,ir),smal)
      IF(rho(iz,ir,iq) > rhotrd) rhoact = rhocri+rhocri
      delr(iz,ir) = MAX(zero,-LOG(rhoact/rhocri))
    END DO
  END DO
  
!     compute proper scalaing and slope
  
  iz1 = nz
  DO iz=nz-1,0,-1
    IF(delr(iz,0) <= zero) EXIT
    iz1 = iz
  END DO

  IF(iz1 > nz-3) STOP ' grid-box to small for surfac density'
  IF(iz1 <= 0) STOP ' something is wrong with the surface density'
  slope  = (delr(iz1+2,0)-delr(iz1,0))/(dz+dz)
  slpfac = one/slope
  delr(0:nz,0:nr)      = slpfac*delr(0:nz,0:nr)
  potadd(0:nz,0:nr,iq) = stslop*delr(0:nz,0:nr)
  WRITE(6,'(a,8g12.5)') ' potadd:',(potadd(iz,0,1),iz=1,8)
  
END DO
DEALLOCATE(wrk1)
RETURN
END SUBROUTINE stapot

!---------------------------------------------------------------------------  
! DESCRIPTION: damnuc
!> @brief
!!Kinetic energy damping applied to the the wavefunction \c psiin, result is 
!!returned on \c psiout, \c ist is the absolute position of the state in the 
!!spectrum. 
!>
!> @param[in] psiin
!> REAL(8), array, takes the input wave function.
!> @param[out] psiout
!> REAL(8), array, returns the Hamiltonian applied to the input wave 
!!function.
!> @param[in] ist
!> INTEGER, absolute adress of the state.
!--------------------------------------------------------------------------- 

SUBROUTINE damnuc(psiin,psiout,ist)

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)                         :: psiin(0:nz,0:nr,nsd)
REAL(8), INTENT(OUT)                        :: psiout(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)                         :: ist

REAL(8) :: step
INTEGER :: mm, is

!     the straightforward imaginary-time step in case of "e0inv=0"

IF(e0inv == zero) THEN
  psiout(0:nz,0:nr,1:2) = x0inv*psiin(0:nz,0:nr,1:2)
ELSE
  
!     now the full kinetic energy damping
  step   = x0inv/e0inv
  CALL damp_nucleon(rdnuc(0:nr,0:nr,:),zdnuc(0:nz,0:nz,:),psiin(0:nz,0:nr,:),&
                    psiout(0:nz,0:nr,:),nparz(1,ist),mpsi(ist),step)

!       set boundary zeroes if parity requires it.
  DO is=1,2
    mm = mpsi(ist)+is-1
    IF (mm > 0) THEN
      psiout(0:nz,0,is)  =  zero
    END IF
    IF (nparz(is,ist) == 2) THEN
      psiout(0,0:nr,is)  =  zero
    END IF
  END DO
  
END IF

RETURN
END SUBROUTINE damnuc

!---------------------------------------------------------------------------  
! DESCRIPTION: schmit
!> @brief
!!subroutine to orthonormalize the set of states \c psiw new states are 
!!replaced in \c psiw.
!> @details
!!If \c tsym=.TRUE. then \c npar (taken from fields \c npar(ist)) is a good 
!!quantum number and only those states having the same z-parity need be 
!!orthogonalized. If tsym is .FALSE. all states need be orthogonalized. The 
!!routine runs over one 'm'-block of the states. The wave functions are 
!!contained in \c psiw running from \c n = 1 through \c nhi. The 
!!corresponding global quantum number is computed with \c istadd as base 
!!point: \c ist = \c istadd + \c n . The norm of each component separately 
!!is returned on \c cmpn ; this serves as a diagnostic tool.                                                        *
!!The orthonormalization should be done in order of increasing energy. The 
!!s.p. energies on \c spenrg are first sorted in energy order. This is done 
!!on a workspace \c eord. The ordering is stored in the pointer \c isort. 
!!\c eord and \c spenrg.
!>
!> @param[in,out] psiw
!> REAL(8), array, takes the input wave functions.
!> @param[in] istadd
!> INTEGER, takes the base point of the wave functions.
!> @param[in] nhi
!> INTEGER, takes the number of the wave functions to be orthonormalized.
!--------------------------------------------------------------------------- 

SUBROUTINE schmit(psiw,istadd,nhi)

USE paramr
USE dynpar
USE levels
USE grid
USE wavefr, ONLY: wfovlp,wfproj,wfovpc
!$ use omp_lib

IMPLICIT NONE

REAL(8), INTENT(IN OUT)     :: psiw(0:nz,0:nr,nsd,nhi)
INTEGER, INTENT(IN)         :: istadd
INTEGER, INTENT(IN)         :: nhi

INTEGER,PARAMETER :: iford=1

REAL(8) :: eord(nhwd)
INTEGER :: isort(nhwd)

REAL(8) :: emin,ovlp,sum2,wnor
INTEGER :: i,imin,in,in1,inw,is,isav,ist,jn,jnw,jst,n,npari,nparj,mythread

IF(tstab .AND. iford /= 1)  &
    STOP ' TSTAB=.TRUE. needs energy ordering in SCHMIT. set IFORD=1'

!**********************************************************************
!                                                                     *
!     energy ordering on workspace 'eord' to pointer 'isort'          *
!                                                                     *
!**********************************************************************
!  mythread=OMP_GET_THREAD_NUM()
IF(nhi > nhwd) STOP ' in SCHMIT: #states "nhi" > supply "nhwd" '
DO in=1,nhi
  eord(in)  = spenrg(in+istadd)
  isort(in) = in+istadd
END DO
IF(iford == 1) THEN
  DO n=1,nhi
    emin = 1.0E32
    imin   = n
    DO i=n,nhi
      IF(eord(i) < emin) THEN
        emin   = eord(i)
        imin   = i
      END IF
    END DO
    isav         = isort(imin)
    eord(imin)   = eord(n)
    isort(imin)  = isort(n)
    eord(n)      = emin
    isort(n)     = isav
  END DO
END IF
!IF(thfb) WRITE(6,*) 'SCHMIT 1: mythread=',mythread
!**********************************************************************
!                                                                     *
!     Perform Gram-Schmidt orthonormalization of ordered w.f.         *
!                                                                     *
!**********************************************************************
DO  in=1,nhi
  ist = isort(in)
  inw = ist-istadd
  npari   = npar(ist)
  in1    = in-1
!IF(thfb) WRITE(6,*) 'SCHMIT 2: mythread,ist=',mythread,ist,SUM(psiw(0:nz,0:nr,1:2,inw)**2)
  
  DO  jn=1,in1
    jst = isort(jn)
    jnw = jst-istadd
    
!        test if parity is a good quantum number select those states  *
    nparj   = npar(jst)
    IF(nparj /= npari .AND. tsym) CYCLE
    
!     compute overlap
    ovlp   = wfovlp(psiw(0:nz,0:nr,1:2,inw),psiw(0:nz,0:nr,1:2,jnw))
!IF(thfb) WRITE(6,*) 'SCHMIT 2: mythread,ist,jst=',mythread,ist,jst,ovlp
    
!     extract psi(jn)-component in psi(in)
    psiw(0:nz,0:nr,1:2,inw)=psiw(0:nz,0:nr,1:2,inw) -ovlp*psiw(0:nz,0:nr,1:2,jnw)
    
  END DO
  
!     normalise state psi(in), and add up partial norm.
  wnor   = wfovlp(psiw(0:nz,0:nr,1:2,inw),psiw(0:nz,0:nr,1:2,inw))
!IF(thfb) WRITE(6,*) 'SCHMIT 2: mythread,wnor=',mythread,wnor
  IF(wnor < smal) THEN
    WRITE(7,'(a,4i5,g12.4)')  &
        ' SCHMIT: ist,inw,istadd,nhi,wnor=',ist,inw,istadd,nhi,wnor
    STOP ' too small norm in SCHMIT'
  END IF
  wnor   = one/SQRT(wnor)
  
  psiw(0:nz,0:nr,1:2,inw) = psiw(0:nz,0:nr,1:2,inw)*wnor
  
  sum2 = zero
  DO is=1,nsd
    cmpn(is,ist) = wfovpc(psiw(0:nz,0:nr,is,inw))
    sum2         = cmpn(is,ist) + sum2
  END DO
  spn(ist) = SQRT(sum2)
  
END DO

RETURN
END SUBROUTINE schmit

!---------------------------------------------------------------------------  
! DESCRIPTION: projec
!> @brief
!!Subroutine to project out of \c psiw the states starting at \c istadd and 
!!going up to \c nhi. \c nparps is the z-parity of \c psiw. \c offnor is 
!!output of the squared sum of projected elements.
!> @details
!>
!> @param[in,out] psiw
!> REAL(8), array, takes the input wave functions.
!> @param[in] nparps
!> INTEGER, takes the z-parity.
!> @param[in] istadd
!> INTEGER, takes the base point of the wave functions.
!> @param[in] nhi
!> INTEGER, takes the highest label of the wave functions to be 
!!projected out.
!> @param[out] offnor
!> REAL(8), returns the squared sum of projected elements.
!--------------------------------------------------------------------------- 

SUBROUTINE projec(psiw,nparps,istadd,nhi,offnor)

USE paramr
USE dynpar
USE wavefr
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN OUT)                  :: psiw(0:nz,0:nr,nsd)
INTEGER, INTENT(IN)                      :: nparps
INTEGER, INTENT(IN)                      :: istadd
INTEGER, INTENT(IN)                      :: nhi
REAL(8), INTENT(OUT)                     :: offnor

REAL(8) :: ovlp,sumnor
INTEGER :: in,ist,npari

sumnor = 0.0D0
ist    = istadd
DO  in=1,nhi
  ist     = ist+1
  npari   = npar(ist)
!        test if parity is a good quantum number select those states  *
  IF(nparps == npari .OR. .NOT.tsym) THEN
    
!     compute overlap
    ovlp   = wfovlp(psi(0:nz,0:nr,1:2,ist),psiw(0:nz,0:nr,1:2))
    sumnor = ovlp*ovlp + sumnor

!     extract psi(in)-component in psiw
    psiw(0:nz,0:nr,1:2)=psiw(0:nz,0:nr,1:2)-ovlp*psi(0:nz,0:nr,1:2,ist)
  END IF
  
END DO
offnor = sumnor

RETURN
END SUBROUTINE projec

!---------------------------------------------------------------------------  
! DESCRIPTION: vsqtr
!> @brief
!!subroutine to dynamically compute the occupation numbers by simple filling 
!!in energy order from below. The energy order is returned in \c iorder. 
!> @details
!!In case of pairing the first guess to the occupation is done here.
!>
!> @param[out] vsqf
!> REAL(8), array, returns resulting occupation numbers for the filling 
!!approximation.
!> @param[in] wate
!> REAL(8), array, takes spin and ml degeneracy factors for the orbits. 
!> @param[in] nst
!> INTEGER, takes the number of single particle states of charge q.
!> @param[in] eord
!> REAL(8), array, takes s.p. energies on which ordering is performed.
!> @param[in] npiq
!> INTEGER, takes the number of particles of charge q.
!> @param[out] iorder
!> INTEGER, array, ordering of the states by energy note the original order 
!!is by q and then by ml. 
!--------------------------------------------------------------------------- 

SUBROUTINE vsqtr(vsqf,wate,nst,eord,npiq,iorder)

USE paramr
IMPLICIT NONE

REAL(8), INTENT(OUT)                     :: vsqf(nst)
REAL(8), INTENT(IN)                      :: wate(nst)
INTEGER, INTENT(IN)                      :: nst
REAL(8), INTENT(IN)                      :: eord(nst)
INTEGER, INTENT(IN OUT)                  :: npiq
INTEGER, INTENT(OUT)                     :: iorder(nst)

REAL(8) :: rpcur,ve2
INTEGER :: iosav,j

!     determine energy ordering and fill 'vsqf'
rpcur=0.0
CALL update_energy_order(nst,eord,iorder)
DO  j=1,nst
  iosav      = iorder(j)
  ve2=MIN(one,(npiq-rpcur)/wate(iosav))
  ve2=MAX(ve2,zero)
  rpcur=rpcur+wate(iosav)*ve2
  vsqf(iosav)=ve2
END DO

RETURN
END SUBROUTINE vsqtr

!---------------------------------------------------------------------------  
! DESCRIPTION: update_energy_order
!> @brief
!!The energy order is returned in \c iorder. 
!>
!> @param[in] nst
!> INTEGER, takes the number of single particle states of charge q.
!> @param[in] eord
!> REAL(8), array, takes s.p. energies on which ordering is performed.
!> @param[out] iorder
!> INTEGER, array, ordering of the states by energy note the original order 
!!is by q and then by ml. 
!--------------------------------------------------------------------------- 

SUBROUTINE update_energy_order(nst,eord,iorder)

USE paramr
IMPLICIT NONE

INTEGER, INTENT(IN)                      :: nst
REAL(8), INTENT(IN)                      :: eord(nst)
INTEGER, INTENT(OUT)                     :: iorder(nst)

REAL(8) :: emin
INTEGER :: iosav,j,jb,jbm

!     initialize pointer field with trivial ordering

DO j=1,nst
  iorder(j) = j
END DO


DO  j=1,nst
  emin=1.0E12
  jbm = j
  DO  jb=j,nst
    IF(eord(iorder(jb)) < emin) THEN
      emin=eord(iorder(jb))
      jbm=jb
    END IF
  END DO
  iosav      = iorder(j)
  iorder(j)  = iorder(jbm)
  iorder(jbm)= iosav
END DO

RETURN
END SUBROUTINE update_energy_order

!---------------------------------------------------------------------------  
! DESCRIPTION: shiftcm
!> @brief
!!Shifts the wavefunctions such that the center-of-mass tends to zero. zcm  
!!= c.m. at input, corrected c.m. at output
!>
!> @param[in] shift
!> REAL(8), takes c.m. correction.
!> @param[in] jst
!> INTEGER, takes the label of wave function to be shifted. 
!--------------------------------------------------------------------------- 

SUBROUTINE shiftcm(shift,jst)

USE paramr
USE dynpar
USE fields
USE wavefr
USE levels
USE grid
IMPLICIT NONE

REAL(8), INTENT(IN)   :: shift
INTEGER, INTENT(IN)   :: jst

!REAL(8),PARAMETER     :: addzcm=-0.5D0

REAL(8),ALLOCATABLE  :: psi1(:,:,:)
REAL(8),ALLOCATABLE,SAVE  :: expmat(:,:),forward(:,:)
INTEGER :: i1,i2,iz
REAL(8) :: workz(0:nz),sh2
REAL(8),SAVE :: shiftold=1D99
LOGICAL,PARAMETER :: tunitary=.FALSE. ! switch between 1. order and 4. order unitary cm shift
! External procedures defined in LAPACK
EXTERNAL DGETRF
EXTERNAL DGETRI

IF(tsym) STOP 'SHIFTCM only without reflection symmetry'
!WRITE(*,*) 'SHIFTCM: shift,tunitary=',shift,tunitary


ALLOCATE(psi1(0:nz,0:nr,nsd))

IF(tunitary) THEN

  IF(.NOT.ALLOCATED(expmat)) ALLOCATE(expmat(0:nz,0:nz),forward(0:nz,0:nz))

  IF(shift.NE.shiftold) THEN
    
    shiftold = shift
! compose translation operator
    forward=0D0
    DO iz=0,nz
      forward(iz,iz) = 1D0
    END DO
    sh2 = shift *0.5D0
    forward = forward + sh2*derz1(:,:,1)+0.5D0*sh2**2*derz2(:,:,1)
    expmat = forward - sh2*derz1(:,:,1)+0.5D0*sh2**2*derz2(:,:,1)
    CALL matr_inv(expmat,nz+1,nz+1)

    DO i1=0,nz
      DO i2=0,nz
        workz(i2) =  SUM(expmat(i1,:)*forward(:,i2))
!      expmat(i1,i2) =  SUM(forward(i1,:)*forward(:,i2))
      END DO
      expmat(i1,:)=workz(:)
    END DO
!    WRITE(*,*) 'EXPMAT generated'
  END IF

!  expmat=forward
!  DO is=1,2; DO i2=0,nr
!    DO iz=0,nz
!      workz(iz) = SUM(expmat(iz,:)*psi(:,i2,is,jst))
!    END DO
!    psi(:,i2,is,jst) = workz(:)
!  END DO; END DO
!  IF(jst==1) WRITE(*,*) 'EXPMAT multiplied'
  CALL mulz(expmat,psi(:,:,1,jst),psi1(:,:,1),1D0,0)
  CALL mulz(expmat,psi(:,:,2,jst),psi1(:,:,2),1D0,0)
  psi(:,:,:,jst) = psi1

ELSE

!  IF(jst==1) WRITE(*,*) 'SHIFTCM with shift=',shift
  psi1=0D0
  CALL nablaz_spin(psi(0:nz,0:nr,1:2,jst),psi1(0:nz,0:nr,1:2),nparz(1,jst),1,1.0D0)

  psi(0,0:nr,1:2,jst) = 0.0D0
  psi(1:nz,0:nr,1:2,jst) = psi(1:nz,0:nr,1:2,jst) + shift*psi1(1:nz,0:nr,1:2)

END IF


DEALLOCATE(psi1)
RETURN
END SUBROUTINE shiftcm


END MODULE steputil

