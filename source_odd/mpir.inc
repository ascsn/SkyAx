      include 'mpif.h'
      integer mpierror,mpirank,mpisize
      integer mpichunks(0:255,2)
      common /MPI/ mpierror,mpirank,mpisize,mpichunks
