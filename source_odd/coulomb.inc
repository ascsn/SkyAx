c=======================================================================
c
c     common block "coulomb"
c
c     the common fields of the coulomb solver
c
c       cbgrho :  background densities
c       cbgfld :  fields of the background densities
c       cbgmom :  multipole moments of the background densities
c
c       cfldbuf:  buffer               
c       cdmpmas:  damping mass for coulomb iteration
c
c=======================================================================
      common  /coulomb/ cbgrho(0:nzd,0:nrd,0:KANGMOL,0:KANGMOM),
     &                  cbgfld(0:nzd,0:nrd,0:KANGMOL,0:KANGMOM),
     &                  cbgmom(0:KANGMOL),
     &                  cfldbuf(0:nzd,0:nrd),
     &                  cdmpmas(0:nzd,0:nrd)
c     &                  nqrun,nhrun,ndrun,norun