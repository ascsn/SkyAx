c-----radial_dens ------------ part of the skyax package ----------------------
c
      SUBROUTINE radial_dens(rhoin)
c
      include 'paramr.inc'
c
c     radial density for L=0,2,4 from 2D array 'rhoin'
c
      include 'dynpar.inc'

      REAL(8),INTENT(IN) :: rhoin(0:nzd,0:nrd)
c
      REAL(8),PARAMETER :: delrad=0.15
      REAL(8),PARAMETER :: y20fac=SQRT(5D0/(16D0*3.1415926D0))
      REAL(8),PARAMETER :: y40fac=0.846284375D0     !    24/(16*SQRT(PI))
      REAL(8),ALLOCATABLE :: rhorad(:)
      LOGICAL,PARAMETER :: ttest=.FALSE.
      
c
c     prepare parameters, arrays and output file
c
      radmax=MIN(dz*(nz-1),dr*(nr-1))
      nrad=INT(radmax/delrad)
      ALLOCATE(rhorad(0:nrad))

c
c     run along the radial grid
c
      WRITE(71,'(f6.2,3(1pg14.6))') 0D0,rhoin(0,0),0D0,0D0
      DO irad=1,nrad-1
        radact=irad*delrad
        radact2=radact*radact
        numz=2*NINT(radact/dz)+40
        dzact=radact/(numz)
        IF(ttest) WRITE(*,*) 'irad,numz,radact,dzact=',
     &        irad,numz,radact,dzact
        acc0=0D0
        acc2=0D0
        DO iz=0,numz
          zact=iz*dzact
          zact2=zact*zact
          ract2=(radact2-zact2)
          ract=SQRT(MAX(ract2,1D-200))
          IF(ttest) WRITE(*,*) 'iz,zact,ract=',iz,zact,ract
          rhoact=rho_intpol2D(rhoin,zact,ract)
          IF(iz==0 .OR. iz==numz) rhoact=rhoact*0.5D0
          acc0=acc0+rhoact
          acc2=acc2+rhoact*(2D0*zact*zact-ract*ract)
          acc2=acc2+rhoact*(2D0*zact*zact-ract*ract)
          acc4=acc4+rhoact*(zact2*zact2-3D0*zact2*ract2
     &                      +0.375D0*ract2*ract2)
        END DO
        acc0=acc0/numz !*dzact
        acc2=acc2*y20fac/(radact2*numz)!*dzact
        acc4=acc4*y40fac/(radact2**2*numz)!*dzact
        WRITE(71,'(f6.2,4(1pg14.6),i4)')radact,acc0,acc2,acc4,dzact,numz
      END DO

c
      DEALLOCATE(rhorad)
      RETURN

      END SUBROUTINE radial_dens

c-----intpol2D ------------ part of the skyax package ----------------------
c
      REAL(8) FUNCTION rho_intpol2D(rhoin,zint,rint)
c
      include 'paramr.inc'
c
c     Interpolates density from array 'rhoin' to points 'zint','rint'
c
      include 'fields.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
c
      REAL(8),INTENT(IN) :: rhoin(0:nzd,0:nrd),zint,rint
c
      REAL(8) :: rhoaux(4,4),zaux(4),raux(4),rhoauxr(4)
      LOGICAL,PARAMETER :: ttest=.FALSE.
      
c
      IF(isym==0) RETURN !STOP 'rho_intpol2D requires z-symmetry'

!          r     = dr*(ir+sr)
!            z      = ((iz)+sz-icm)*dz
!      WRITE(6,'(a)')'rhoin(0:4,0:4):'
!      DO ir=0,4
!        WRITE(6,'(i2,5(1pg13.5))') ir,rhoin(0:4,ir)
!      END DO
c
c     find closest grid points for cubic interpolation
c     and compose interpolation arrays
c
      nrlow=INT(rint/dr)
      IF(nrlow>nr-2) THEN
        rho_intpol2D=0D0
      END IF
      IF(nrlow<0) STOP 'lower r-grid point <0'
      nzlow=INT(zint/dz)
      IF(nzlow>nz-2) THEN
        rho_intpol2D=0D0
      END IF
      IF(nzlow<0) THEN
         WRITE(6,*) 'nzlow,nrlow,zint,rint,dz=',nzlow,nrlow,zint,rint,dz
         STOP 'lower z-grid point <0'
      END IF
        
         
      zaux=(/(nzlow-1)*dz,nzlow*dz,(nzlow+1)*dz,(nzlow+2)*dz/)
      raux=(/(nrlow-1)*dr,nrlow*dr,(nrlow+1)*dr,(nrlow+2)*dr/)
      IF(ttest) THEN
         WRITE(6,'(a,2i4)') 'nzlow,nrlow=',nzlow,nrlow
         WRITE(6,'(a,2(1pg13.5)/4(1pg13.5)/4(1pg13.5))')
     &         'zint,rint,aux:',zint,rint,zaux,raux
         flush(6)
      END IF
      DO iz=1,4
        izin=iz-2+nzlow
        DO ir=1,4
          irin=ir-2+nrlow
          IF(izin<0 .AND. irin<0) THEN
            rhoaux(iz,ir)=rhoin(1,1)
          ELSE IF(izin<0) THEN
            rhoaux(iz,ir)=rhoin(1,irin)
          ELSE IF(irin<0) THEN
            rhoaux(iz,ir)=rhoin(izin,1)
          ELSE
            rhoaux(iz,ir)=rhoin(izin,irin)
          END IF
      END DO
      END DO
      IF(ttest) THEN
         WRITE(6,'(a/4(i2,4(1pg13.5)/))')
     &       'rhoaux:',(ir,rhoaux(:,ir),ir=1,4)
         flush(6)
      END IF

c
c     interpolate along z
c
      DO ir=1,4
        CALL polint(zaux,rhoaux(:,ir),4,zint,rhoauxr(ir),errorz)
      END DO

c
c     interpolate along r
c
      IF(ttest) WRITE(6,'(a,4(1pg13.5))') 'rhoauxr:',rhoauxr
      flush(6)
      CALL polint(raux,rhoauxr,4,rint,fin_inter,errorr)

      rho_intpol2D=fin_inter

      RETURN

      END FUNCTION rho_intpol2D


c-----radial_dens ------------ part of the skyax package ----------------------
c
      SUBROUTINE test_intpol2D(rhoin)
c
      include 'paramr.inc'
c
c     test 2D interpolation
c
      include 'dynpar.inc'

      REAL(8),INTENT(IN) :: rhoin(0:nzd,0:nrd)
c
      REAL(8),PARAMETER :: delrad=0.3

      WRITE(*,*) 'test 2D interpolation:'
      DO i=1,10
        radact=i*delrad
        WRITE(6,*) 'now radact=',radact
        rhointz0=rho_intpol2D(rhoin,0D0,radact)
        rhointr0=rho_intpol2D(rhoin,radact,0D0)
        rhointal=rho_intpol2D(rhoin,radact,radact)
        WRITE(6,'(f6.2,3(1pg13.5))') radact,rhointz0,rhointr0,rhointal
      END DO

      RETURN

      END SUBROUTINE test_intpol2D
