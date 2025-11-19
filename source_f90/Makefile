q/SHELL 	       = /bin/bash
#
# ----- targets--------------------------------------------------------------
#
.PHONY : all info params gf gf_d ifort params info gf_omp
#
# ----- define the compiler being used---------------------------------------
# TODO: may need to change ifort to ifx
gf gf_omp gf_d  : export COMPILER_SKYAX = gfortran
ifort ifort_omp : export COMPILER_SKYAX = ifort

# ----- define HDF5 flagser being used---------------------------------------
# I have $(HDF5_DIR) set in as an environment variable. On my Ubuntu machine 
# it points to /usr/lib/x86_64-linux-gnu/hdf5/serial/
HDF5_INCLUDE = $(HDF5_DIR)/include
HDF5_LIB = $(HDF5_DIR)/lib

# Linker flags to link with the HDF5 library
H5_LDFLAGS = -L$(HDF5_LIB) -lhdf5_fortran -lhdf5

# Compiler flags
H5_FFLAGS = -I$(HDF5_INCLUDE)

export H5_LDFLAGS
export H5_FFLAGS
#
# ----- define compiler flags and linker flags-------------------------------
#
#gf        : export COMPILERFLAGS_SKYAX = -O3 -msse4.2 -mavx
#gf_omp    : export COMPILERFLAGS_SKYAX = -O3 -fopenmp  -msse4.2 -mavx
# NOTE: the -fPIC flag was added because of how HDF5 was built on my Ubuntu 22.x machine (DL, 11/18/25).
# For omp code I probably need that flag as well; I don't know how this will port to other
# machines!
# gf        : export COMPILERFLAGS_SKYAX = -Ofast -march=native
gf        : export COMPILERFLAGS_SKYAX = -Ofast -march=native -fPIC
gf_omp    : export COMPILERFLAGS_SKYAX = -Ofast -march=native -fopenmp
gf_d      : export COMPILERFLAGS_SKYAX = -g -fbacktrace -Wall -fbounds-check -finit-real=NaN
ifort     : export COMPILERFLAGS_SKYAX = -O3 -xhost -mcmodel=large
ifort_omp : export COMPILERFLAGS_SKYAX = -O3 -xhost -qopenmp -mcmodel=large
#
# ----- executables----------------------------------------------------------
#
gf        : export EXEC_SKYAX = skyax.gfortran.seq
gf_omp    : export EXEC_SKYAX = skyax.gfortran.omp
gf_d      : export EXEC_SKYAX = skyax.gfortran.dbg
ifort     : export EXEC_SKYAX = skyax.ifort.seq
ifort_omp : export EXEC_SKYAX = skyax.ifort.omp
#
# ----- libraries--------------------------------------------------------------
# TODO: change here for lapack library
gf gf_d gf_omp ifort ifort_omp: export LIBS_SKYAX = -llapack -lblas
#gf gf_d ifort : export LIBS_SKYAX = -L$(/usr/lib) -lblas -L$(/usr/lib) -llapack 
#
# ----- restart make, clean or print info------------------------------------
#

# make all targets
all: 	info params

# clean all objects and temporary files
clean:  
	@rm -f *.o
	@rm -f *.mod
#
clean-exec: clean
	@rm -fr ../skyax.*.*
#
info: 
	@echo "------------------------------------------------------"
	@echo "    SkyAx  - the axial Skyrme-Hartree-Fock Code       "
	@echo "             >>> version 1.0.0 <<<                    "         
	@echo "------------------------------------------------------"
	@rm -f compile.log
#
params:
	@echo "  syntax: "
	@echo "     make gf        compile and link SkyAx with gfortran"
	@echo "     make ifort     compile and link SkyAx using ifort"
	@echo "     make gf_omp    compile and link SkyAx with gfortran using OpenMP"
	@echo "     make ifort_omp compile and link SkyAx using ifort using OpenMP"
	@echo
	@echo "     make clean       remove objects and temporary files"
	@echo "     make clean-exec  remove objects and temporary files and executables"
	@echo 
	@echo "     make all      print version information and this help screen"
	@echo 


gf gf_d gf_omp ifort ifort_omp : 
	@$(MAKE) -f $(firstword $(MAKEFILE_LIST)) makeit

#
# ***** after restart********************************************************
#

LINKER_SKYAX         = $(COMPILER_SKYAX)
LINKERFLAGS_SKYAX    = $(COMPILERFLAGS_SKYAX)

SPATH           = .

MODULES        =  m_paramr.o m_dynpar.o m_levels.o m_grid.o m_wavefr.o \
                  m_skyax_hdf5.o \
                  m_fields.o m_constr.o m_observ.o m_coulomb.o \
                  m_staticstep.o m_pairing.o m_util.o m_init.o \
                  m_hdiag.o m_steputil.o user_interface.o \
                  skyax.o m_fyfunc_2D.o
# MODULES        =  m_paramr.o m_dynpar.o m_levels.o m_grid.o m_wavefr.o \
#                   m_fields.o m_constr.o m_observ.o m_coulomb.o \
#                   m_staticstep.o m_pairing.o m_util.o m_init.o \
# 	 	  m_hdiag.o m_steputil.o user_interface.o \
#                   skyax.o m_fyfunc_2D.o

#
#-----------------------------------------------------------------------------
makeit : ../$(EXEC_SKYAX)

m_skyax_hdf5.o: Makefile

m_paramr.o : Makefile

m_grid.o : m_paramr.o Makefile

m_dynpar.o : m_paramr.o forces.data Makefile

m_levels.o : m_paramr.o Makefile

m_fields.o : m_paramr.o m_grid.o Makefile

m_wavefr.o : m_paramr.o m_grid.o Makefile

m_constr.o : m_paramr.o m_dynpar.o m_levels.o m_grid.o m_fields.o Makefile

m_hdiag.o : m_grid.o m_levels.o m_paramr.o m_wavefr.o Makefile

# m_observ.o : m_paramr.o m_grid.o m_constr.o m_dynpar.o m_wavefr.o \
#              m_levels.o m_fields.o Makefile
m_observ.o : m_skyax_hdf5.o m_paramr.o m_grid.o m_constr.o m_dynpar.o \
             m_wavefr.o m_levels.o m_fields.o Makefile

m_coulomb.o : m_paramr.o m_fields.o m_wavefr.o m_levels.o m_grid.o \
              m_constr.o m_dynpar.o m_observ.o Makefile

m_steputil.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
               m_levels.o m_observ.o m_paramr.o m_wavefr.o Makefile

m_fyfunc_2D.o : m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
             m_observ.o m_paramr.o m_wavefr.o Makefile

m_util.o : m_dynpar.o m_fields.o m_grid.o m_levels.o m_observ.o \
           m_paramr.o m_wavefr.o Makefile

m_pairing.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
              m_levels.o m_observ.o m_paramr.o m_wavefr.o m_fyfunc_2D.o Makefile

m_init.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
           m_levels.o m_observ.o m_paramr.o m_wavefr.o m_pairing.o \
           m_util.o m_steputil.o Makefile

m_staticstep.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
                 m_grid.o m_levels.o m_observ.o m_paramr.o m_wavefr.o m_hdiag.o \
		         m_pairing.o m_util.o m_steputil.o Makefile

# user_interface.o : Makefile
user_interface.o : m_fields.o m_skyax_hdf5.o Makefile

# skyax.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
#           m_grid.o m_init.f90 m_levels.o m_observ.o m_paramr.o \
#           m_wavefr.o m_staticstep.o m_util.o Makefile
skyax.o : m_skyax_hdf5.o m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
          m_grid.o m_init.f90 m_levels.o m_observ.o m_paramr.o \
          m_wavefr.o m_staticstep.o m_util.o Makefile

# executables
../$(EXEC_SKYAX): $(MODULES)
	 @echo -n "linking " $@ "... "
	 $(LINKER_SKYAX) $(LINKERFLAGS_SKYAX) $(H5_LDFLAGS) -o $@ $(MODULES) $(LIBS_SKYAX)
	 @echo "OK"

.SUFFIXES:  .f90 .o

%.o : %.mod

%.o: %.f90
	@echo -n "compiling "$<" ... "
	$(COMPILER_SKYAX) $(COMPILERFLAGS_SKYAX) $(H5_FFLAGS) -c $<
	@echo "OK"

