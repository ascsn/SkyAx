
c-----fourf2D ------------ part of the skyax package ----------------------
c
      subroutine fourf2D(cmwid,qr,fop2D,fon2D)
c
      include 'paramr.inc'
c
c     Computes the 2D proton and neutron formfactors
c
c     Input variables are:
c      cmwid : parameter for center of mass correction
c      qr    : array of radial momenta
c     Grid parameters comein through 'dynpar.inc'.
c     The densities and the grid parameter are entered via common.
c     output paramaters:
c      fop2D  :   proton formfactor
c      fon2D  :   neutron formfactor
c
      include 'fields.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
c
      REAL(8),INTENT(IN) :: cmwid,qr(0:nr)
      REAL(8),INTENT(OUT) :: fop2D(0:nz,0:nr),fon2D(0:nz,0:nr)
c
c
c----------------------------------------------------------------
c
      IF(isym=0) STOP 'FOURF2D only for z symmetry'
c
c     outer loop over q grid
c 
      zmaxpi=dz*nz/pi
      DO iqr=0,nr; DO iqz=0,nz
        qz=zmaxpi/(0.5D0+iqz)
c
c     inner loop over coordinate space grid
c
        do ir=0,nr
          r     = dr*(ir+sr)
          besfacr=bessel_j0(qr(iqr)*r)
          do iz=0,nz
            z      = ((iz)+sz-icm)*dz
            besfac=COS(qz*z)*besfacr
            wrk1(iz,ir) = rho(iz,ir,1)*besfac
            wrk2(iz,ir) = rho(iz,ir,2)*besfac
          enddo
        enddo
        sump  = f1int(wrk1)
        sumn  = f1int(wrk2)
c
c     c.m. correction
c
        qa     = qz**2+qr(iqr)**2
        exphf  = exp(cmwid*qa)
c
c     finally assembling the formfactors
c
        fop2D(iqz,iqr)    = exphf*sump
        fon2D(iqz,iqr)    = exphf*sumn
c
      END DO;END DO
c
      return
      end
c-----fourb2D ------------ part of the skyax package ----------------------
c
      subroutine fourb2D(fop2D,fon2D,)
c
      include 'paramr.inc'
c
c     Computes the 2D proton and neutron formfactors
c
c     Input variables are:
c      cmwid : parameter for center of mass correction
c      qr    : array of radial momenta
c     Grid parameters comein through 'dynpar.inc'.
c     The densities and the grid parameter are entered via common.
c     output paramaters:
c      fop2D  :   proton formfactor
c      fon2D  :   neutron formfactor
c
      include 'fields.inc'
      include 'dynpar.inc'
      include 'wavefr.inc'
c
      REAL(8),INTENT(IN) :: cmwid,qr(0:nr)
      REAL(8),INTENT(OUT) :: fop2D(0:nz,0:nr),fon2D(0:nz,0:nr)
c
c
c----------------------------------------------------------------
c
      IF(isym=0) STOP 'FOURF2D only for z symmetry'
c
c     outer loop over q grid
c 
      zmaxpi=dz*nz/pi
      DO iqr=0,nr; DO iqz=0,nz
        qz=zmaxpi/(0.5D0+iqz)
c
c     inner loop over coordinate space grid
c
        do ir=0,nr
          r     = dr*(ir+sr)
          besfacr=bessel_j0(qr(iqr)*r)
          do iz=0,nz
            z      = ((iz)+sz-icm)*dz
            besfac=COS(qz*z)*besfacr
            wrk1(iz,ir) = rho(iz,ir,1)*besfac
            wrk2(iz,ir) = rho(iz,ir,2)*besfac
          enddo
        enddo
        sump  = f1int(wrk1)
        sumn  = f1int(wrk2)
c
c     c.m. correction
c
        qa     = qz**2+qr(iqr)**2
        exphf  = exp(cmwid*qa)
c
c     finally assembling the formfactors
c
        fop2D(iqz,iqr)    = exphf*sump
        fon2D(iqz,iqr)    = exphf*sumn
c
      END DO;END DO
c
      return
      end

c-----prep_four2D------------ part of the skyax package ----------------------
c
      subroutine prep_four2D(qr,rmax,nr)
c
      include 'paramr.inc'
c
c     prepares grid of q_r values for Bessel transformation
c
c     rmax  = maximum r value on r-grid
c     qr    = array of q_r
c     nr    = nr. of r-grid points
c
      REAL(8),INTENT(OUT) :: qr(0:nr)
      REAL(8),INTENT(IN) :: rmax
      INTEGER,INTENT(IN) :: nr
c
      REAL(8) :: beszero(0:100)        ! zeroes of J_0(x)
c

c
c     define zeroes of Bessel function
c
      res(  0)=  2.40482555769577D+00
      res(  1)=  5.52007811028632D+00
      res(  2)=  8.65372791291101D+00
      res(  3)=  1.17915344390143D+01
      res(  4)=  1.49309177084878D+01
      res(  5)=  1.80710639679109D+01
      res(  6)=  2.12116366298793D+01
      res(  7)=  2.43524715307493D+01
      res(  8)=  2.74934791320403D+01
      res(  9)=  3.06346064684320D+01
      res( 10)=  3.37758202135736D+01
      res( 11)=  3.69170983536640D+01
      res( 12)=  4.00584257646282D+01
      res( 13)=  4.31997917131767D+01
      res( 14)=  4.63411883716618D+01
      res( 15)=  4.94826098973978D+01
      res( 16)=  5.26240518411150D+01
      res( 17)=  5.57655107550200D+01
      res( 18)=  5.89069839260809D+01
      res( 19)=  6.20484691902272D+01
      res( 20)=  6.51899648002069D+01
      res( 21)=  6.83314693298568D+01
      res( 22)=  7.14729816035937D+01
      res( 23)=  7.46145006437018D+01
      res( 24)=  7.77560256303881D+01
      res( 25)=  8.08975558711376D+01
      res( 26)=  8.40390907769382D+01
      res( 27)=  8.71806298436412D+01
      res( 28)=  9.03221726372105D+01
      res( 29)=  9.34637187819448D+01
      res( 30)=  9.66052679509963D+01
      res( 31)=  9.97468198586806D+01
      res( 32)=  1.02888374254195D+02
      res( 33)=  1.06029930916452D+02
      res( 34)=  1.09171489649805D+02
      res( 35)=  1.12313050280495D+02
      res( 36)=  1.15454612653667D+02
      res( 37)=  1.18596176630873D+02
      res( 38)=  1.21737742087951D+02
      res( 39)=  1.24879308913233D+02
      res( 40)=  1.28020877006008D+02
      res( 41)=  1.31162446275214D+02
      res( 42)=  1.34304016638305D+02
      res( 43)=  1.37445588020284D+02
      res( 44)=  1.40587160352854D+02
      res( 45)=  1.43728733573690D+02
      res( 46)=  1.46870307625797D+02
      res( 47)=  1.50011882456955D+02
      res( 48)=  1.53153458019228D+02
      res( 49)=  1.56295034268534D+02
      res( 50)=  1.59436611164263D+02
      res( 51)=  1.62578188668947D+02
      res( 52)=  1.65719766747955D+02
      res( 53)=  1.68861345369236D+02
      res( 54)=  1.72002924503078D+02
      res( 55)=  1.75144504121903D+02
      res( 56)=  1.78286084200074D+02
      res( 57)=  1.81427664713731D+02
      res( 58)=  1.84569245640639D+02
      res( 59)=  1.87710826960049D+02
      res( 60)=  1.90852408652582D+02
      res( 61)=  1.93993990700109D+02
      res( 62)=  1.97135573085661D+02
      res( 63)=  2.00277155793332D+02
      res( 64)=  2.03418738808199D+02
      res( 65)=  2.06560322116244D+02
      res( 66)=  2.09701905704294D+02
      res( 67)=  2.12843489559949D+02
      res( 68)=  2.15985073671534D+02
      res( 69)=  2.19126658028041D+02
      res( 70)=  2.22268242619084D+02
      res( 71)=  2.25409827434859D+02
      res( 72)=  2.28551412466099D+02
      res( 73)=  2.31692997704039D+02
      res( 74)=  2.34834583140383D+02
      res( 75)=  2.37976168767276D+02
      res( 76)=  2.41117754577268D+02
      res( 77)=  2.44259340563296D+02
      res( 78)=  2.47400926718653D+02
      res( 79)=  2.50542513036970D+02
      res( 80)=  2.53684099512193D+02
      res( 81)=  2.56825686138564D+02
      res( 82)=  2.59967272910604D+02
      res( 83)=  2.63108859823095D+02
      res( 84)=  2.66250446871066D+02
      res( 85)=  2.69392034049776D+02
      res( 86)=  2.72533621354705D+02
      res( 87)=  2.75675208781537D+02
      res( 88)=  2.78816796326153D+02
      res( 89)=  2.81958383984615D+02
      res( 90)=  2.85099971753160D+02
      res( 91)=  2.88241559628188D+02
      res( 92)=  2.91383147606255D+02
      res( 93)=  2.94524735684065D+02
      res( 94)=  2.97666323858459D+02
      res( 95)=  3.00807912126411D+02
      res( 96)=  3.03949500485021D+02
      res( 97)=  3.07091088931505D+02
      res( 98)=  3.10232677463195D+02
      res( 99)=  3.13374266077528D+02
      res(100)=  3.16515854772043D+02


c
c     deduces q-grid
c
      DO n=0,nr
        qr(n)=rmax/res(n)
      END DO

c
      RETURN
      END
