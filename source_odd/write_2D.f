c-----writ2D-----------------------------------------------------------
c
      subroutine writ2D()
*
***********************************************************************
*                                                                     *
*     save wavefunctions and associated quantities for restart        *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
*
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
c      include 'derivs.inc'
c
c-------------------------------------------------------------------------
c
c     generate output unit on standard name for 2D code
c
      open(55,status='unknown',form='unformatted',
     &        file='restart2D.data')
c
c     write basic information onf grid and levels
c
      WRITE(55) dr,dz,nr,nz
      WRITE(55) anumb,nst,npart,iokk,numsta,lst
c
c     write single-particle information , remind doubling of states
c
      DO iq = 1,2  
        WRITE(55) (mpsi(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (npar(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (nrpar(ist),ist=1,nst(1)+nst(2))
        WRITE(55) ((nparz(is,ist),is=1,2),ist=1,nst(1)+nst(2))
        WRITE(55) (nqr(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (nqz(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (icharg(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (iorder(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (spenrg(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (spn(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (spevar(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (wate(ist),ist=1,nst(1)+nst(2))
        WRITE(55) (watep(ist),ist=1,nst(1)+nst(2))
      END DO
c
c     write wavefunctions
c
      DO ist=1,nst(1)+nst(2)
        WRITE(55) (((psi(iz,ir,is,ist),iz=0,nz),ir=0,nr),is=1,2)
      END DO
      CLOSE(55)
      write(6,'(a)') ' 2D restart file written'
c
c     all done
c
      RETURN
      END       
c-----read2D-----------------------------------------------------------
c
      subroutine read2D()
*
***********************************************************************
*                                                                     *
*     save wavefunctions and associated quantities for restart        *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
*
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
      include 'wavefr.inc'
      dimension nparti(2)
c      include 'derivs.inc'
c
c-------------------------------------------------------------------------
c
c     generate output unit on standard name for 2D code
c
      open(55,status='old',form='unformatted',
     &        file='restart2D.data')
c
c     read basic information onf grid and levels and compare
c
      READ(55) dri,dzi,nri,nzi
      IF(dri.NE.dr .OR. dzi.NE.dz .OR. nri.NE.nr .OR. nzi.NE.nz) THEN
        write(7,'(a,2f8.3,2i5)')
     &   ' grid from input:   dr,dz,nr,nz=',dr,dz,nr,nz,
     &   ' grid from restart: dr,dz,nr,nz=',dri,dzi,nri,nzi
        STOP ' READ2D: grid does not match at restart'
      ENDIF
      READ(55) anumb,nst,nparti,iokk,numsta,lst
      IF(npart(1).NE.nparti(1) .OR. npart(2).NE.nparti(2)) THEN
        write(7,'(a,2i5)')
     &   ' particle number from input:   nprot,nneut=',npart,
     &   ' particle number from restart: nprot,nneut=',nparti
        STOP ' READ2D: particle number does not match '
      ENDIF
c
c     read single-particle information , remind doubling of states
c
      DO iq = 1,2  
        READ(55) (mpsi(ist),ist=1,nst(1)+nst(2))
        READ(55) (npar(ist),ist=1,nst(1)+nst(2))
        READ(55) (nrpar(ist),ist=1,nst(1)+nst(2))
        READ(55) ((nparz(is,ist),is=1,2),ist=1,nst(1)+nst(2))
        READ(55) (nqr(ist),ist=1,nst(1)+nst(2))
        READ(55) (nqz(ist),ist=1,nst(1)+nst(2))
        READ(55) (icharg(ist),ist=1,nst(1)+nst(2))
        READ(55) (iorder(ist),ist=1,nst(1)+nst(2))
        READ(55) (spenrg(ist),ist=1,nst(1)+nst(2))
        READ(55) (spn(ist),ist=1,nst(1)+nst(2))
        READ(55) (spevar(ist),ist=1,nst(1)+nst(2))
        READ(55) (wate(ist),ist=1,nst(1)+nst(2))
        READ(55) (watep(ist),ist=1,nst(1)+nst(2))
      END DO
c
c     read wavefunctions
c
      DO ist=1,nst(1)+nst(2)
        READ(55) (((psi(iz,ir,is,ist),iz=0,nz),ir=0,nr),is=1,2)
      END DO
      CLOSE(55)
      write(6,'(a)') ' 2D restart file read'
c
c     all done
c
      RETURN
      END       


c-----pri_2D -----------------------------------------------------------
c
      subroutine pri_2D(outnam)
*
***********************************************************************
*                                                                     *
*     print 2D info of densities                                      *
*                                                                     *
***********************************************************************
*
      include 'paramr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
*
      character*3  outnam,outnumber
c
      data icallno/0/
c
c----------------------------------------------------------------------
c
      icallno = 1+icallno
      open(77,status='scratch')
      write(77,'(i4)') icallno+1000
      rewind(77)
      read(77,'(1x,a3)') outnumber
      close(77)
      open(77,file='fields_2D.'//outnumber//'.'//outnam)
      write(77,'(a,2i5,a,f6.0//a,3f7.3//a,1pg14.6,1pg12.4)')
     &   '# Nr. proton/neutron=',npart,
     &   ' total Nr. Nucleons=',anumb,
     &   '# deformations beta2,beta3,beta4=',beta2,beta3,beta4,
     &   '# energy, rms-radius=',energy,rms
c
       do ir=0,nr
         r   = (sr+ir)*dr
         do iz=0,nz
           z     = (sz+iz-icm)*dz
           write(77,'(2i5,2f7.2,3(1pg14.6))')
     &        ir,iz,r,z,
     &        rho(iz,ir,1),rho(iz,ir,2),rho(iz,ir,1)+rho(iz,ir,2)
         enddo
         write(77,'(1x)') 
       enddo
c
       close(77)
c
       return
       end          


c======================================================================
c
c     wrt_2D_spin(fldin,filename,label,incridx,outnam)
c
c     print 2D info of a spinor field (wavefunction, or 
c     isospin-dependent fields)
c
c=======================================================================
      subroutine wrt_2D_spin(fldin,filename,incridx,outnam)
c
c----------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
c
c----------------------------------------------------------------------
c
      double precision fldin(0:nzd,0:nrd,2)
      character*3      outnam,outnumber
      character*10     filename
      character*80     label
      integer          incridx
      integer          nfld,ir,iz
      double precision r,z
      data             icallno/0/
c
c----------------------------------------------------------------------
c
c     generate a string containing the index of the file
c----------------------------------------------------------------------
      if (incridx.eq.1) then    ! increase index if incridx==1
         icallno = 1 + icallno
      endif
c
      open(77,status='scratch')
      write(77,'(i4)') icallno + 1000
      rewind(77)
      read(77,'(1x,a3)') outnumber
      close(77)
c
c     open file and write header
c----------------------------------------------------------------------
      open(77,file=filename//'.'//outnumber//'.'//outnam)
      write(77,'(a,2i5,a,f6.0/a,3f7.3/a,1pg14.6,1pg12.4)')
     &   '# Nr. proton/neutron=',npart,
     &   ' total Nr. Nucleons=',anumb,
     &   '# deformations beta2,beta3,beta4=',beta2,beta3,beta4,
     &   '# energy, rms-radius=',energy,rms
c
c     write the field (could either be a spinor or a isospinor)
c----------------------------------------------------------------------
      do ir=0,nr
         r   = (sr+ir)*dr
         do iz=0,nz
            z     = (sz+iz-icm)*dz
            write(77,'(2i5,2f7.2,2(1pg16.5))')
     &           ir,iz,r,z,
     &           fldin(iz,ir,1),fldin(iz,ir,2)
         enddo
         write(77,'(1x)') 
      enddo
c     
      close(77)
c
c----------------------------------------------------------------------
c
      return
      end          


c======================================================================
c
c     wrt_2D_scalar(fldin,filename,label,incridx,outnam)
c
c     print 2D info of a scalar field (like coulomb density)  
c
c=======================================================================
      subroutine wrt_2D_scalar(fldin,filename,incridx,outnam)
c
c----------------------------------------------------------------------
c
      include 'paramr.inc'
      include 'observ.inc'
      include 'dynpar.inc'
      include 'fields.inc'
      include 'levels.inc'
c
c----------------------------------------------------------------------
c
      double precision fldin(0:nzd,0:nrd)
      character*3      outnam,outnumber
      character*10     filename
      character*80     label
      integer          incridx
      integer          nfld,ir,iz
      double precision r,z
      data             icallno/0/
c
c----------------------------------------------------------------------
c
c     generate a string containing the index of the file
c----------------------------------------------------------------------
      if (incridx.eq.1) then    ! increase index if incridx==1
         icallno = 1 + icallno
      endif
c
      open(77,status='scratch')
      write(77,'(i4)') icallno + 1000
      rewind(77)
      read(77,'(1x,a3)') outnumber
      close(77)
c
c     open file and write header
c----------------------------------------------------------------------
      open(77,file=filename//'.'//outnumber//'.'//outnam)
      write(77,'(a,2i5,a,f6.0/a,3f7.3/a,1pg14.6,1pg12.4)')
     &   '# Nr. proton/neutron=',npart,
     &   ' total Nr. Nucleons=',anumb,
     &   '# deformations beta2,beta3,beta4=',beta2,beta3,beta4,
     &   '# energy, rms-radius=',energy,rms
c
c     write the scalar field
c----------------------------------------------------------------------
      do ir=0,nr
         r   = (sr+ir)*dr
         do iz=0,nz
            z     = (sz+iz-icm)*dz
            write(77,'(2i5,2f7.2,1pg16.5)')
     &           ir,iz,r,z,
     &           fldin(iz,ir)
         enddo
         write(77,'(1x)') 
      enddo
c     
      close(77)
c
c----------------------------------------------------------------------
c
      return
      end          

