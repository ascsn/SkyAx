q/SHELL 	       = /bin/bash
#
# ----- targets--------------------------------------------------------------
#
.PHONY : all info params gf gf_d ifort params info gf_omp
#
# ----- define the compiler being used---------------------------------------
#
gf gf_omp gf_d  : export COMPILER_SKYAX = gfortran
ifort ifort_omp ifort_dbg ifort_domp: export COMPILER_SKYAX = ifx
#
# ----- define compiler flags and linker flags-------------------------------
#
#gf        : export COMPILERFLAGS_SKYAX = -O3 -msse4.2 -mavx
#gf_omp    : export COMPILERFLAGS_SKYAX = -O3 -fopenmp  -msse4.2 -mavx
gf        : export COMPILERFLAGS_SKYAX = -Ofast -march=native -mcmodel=large
gf_omp    : export COMPILERFLAGS_SKYAX = -Ofast -march=native -fopenmp
gf_d      : export COMPILERFLAGS_SKYAX = -g -fbacktrace -Wall -fbounds-check -finit-real=NaN -ffpe-trap=invalid,zero,overflow
ifort     : export COMPILERFLAGS_SKYAX = -O3 -xhost -mcmodel=large
ifort_dbg : export COMPILERFLAGS_SKYAX = -g -traceback -CB -mcmodel=large -check uninit
ifort_omp : export COMPILERFLAGS_SKYAX = -O3 -xhost -qopenmp -mcmodel=large
ifort_domp : export COMPILERFLAGS_SKYAX = -g -traceback -CB -qopenmp -mcmodel=large -check uninit
#
# ----- executables----------------------------------------------------------
#
gf        : export EXEC_SKYAX = skyax.gfortran.seq
gf_omp    : export EXEC_SKYAX = skyax.gfortran.omp
gf_d      : export EXEC_SKYAX = skyax.gfortran.dbg
ifort     : export EXEC_SKYAX = skyax.ifort.seq
ifort_dbg : export EXEC_SKYAX = skyax.ifort.dbg
ifort_omp : export EXEC_SKYAX = skyax.ifort.omp
ifort_domp : export EXEC_SKYAX = skyax.ifort.domp
#
# ----- libraries--------------------------------------------------------------
#
# Uncomment below for general lapack and blas libraries
# gf gf_d gf_omp ifort ifort_omp ifort_dbg: export LIBS_SKYAX = -llapack -lblas
# Uncomment below for ICER, MSU.
gf gf_d gf_omp ifort ifort_omp ifort_dbg: export LIBS_SKYAX = -llapack -lblas
# gf gf_d gf_omp ifort ifort_omp ifort_dbg: export LIBS_SKYAX = -lmkl
# ifort ifort_omp ifort_dbg: export LIBS_SKYAX = -qmkl
#gf gf_d gf_omp: export LIBS_SKYAX = -mkl
# gf gf_d ifort : export LIBS_SKYAX = -lblas -llapack 
#gf gf_d ifort : export LIBS_SKYAX = -L$(y/usr/lib) -lblas -L$(/usr/lib) -llapack 
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


gf gf_d gf_omp ifort ifort_omp ifort_dbg ifort_domp: 
	@$(MAKE) -f $(firstword $(MAKEFILE_LIST)) makeit

#
# ***** after restart********************************************************
#

LINKER_SKYAX         = $(COMPILER_SKYAX)
LINKERFLAGS_SKYAX    = $(COMPILERFLAGS_SKYAX)

SPATH           = .

MODULES        =  m_paramr.o m_dynpar.o m_levels.o m_grid.o m_wavefr.o \
                  m_fields.o m_constr.o m_observ.o m_coulomb.o \
                  m_staticstep.o m_pairing.o m_util.o m_init.o \
	 	  m_hdiag.o m_steputil.o user_interface.o \
                  skyax.o m_fyfunc_2D.o hfb-natural.o m_blocking.o

#
#-----------------------------------------------------------------------------
makeit : ../$(EXEC_SKYAX)

m_paramr.o : Makefile

m_grid.o : m_paramr.o Makefile

m_dynpar.o : m_paramr.o forces.data Makefile

m_levels.o : m_paramr.o Makefile

m_fields.o : m_paramr.o m_grid.o Makefile

m_wavefr.o : m_paramr.o m_grid.o Makefile

m_constr.o : m_paramr.o m_dynpar.o m_levels.o m_grid.o m_fields.o Makefile

m_hdiag.o : m_grid.o m_levels.o m_paramr.o m_wavefr.o Makefile

m_observ.o : m_paramr.o m_grid.o m_constr.o m_dynpar.o m_wavefr.o \
             m_levels.o m_fields.o Makefile

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

hfb-natural.o : m_pairing.o m_dynpar.o m_fields.o m_grid.o \
              m_levels.o m_observ.o m_paramr.o m_wavefr.o m_fyfunc_2D.o Makefile

m_staticstep.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
                 m_grid.o m_levels.o m_observ.o m_paramr.o m_wavefr.o m_hdiag.o \
		         m_pairing.o m_util.o m_steputil.o m_init.o Makefile

m_init.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
           m_levels.o m_observ.o m_paramr.o m_wavefr.o m_pairing.o \
           m_util.o m_steputil.o m_blocking.o Makefile

m_blocking.o : m_dynpar.o m_grid.o m_levels.o m_paramr.o m_wavefr.o m_pairing.o \
               m_staticstep.o m_util.o Makefile

user_interface.o : Makefile

skyax.o : m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
          m_grid.o m_init.f90 m_levels.o m_observ.o m_paramr.o \
          m_wavefr.o m_staticstep.o m_util.o m_blocking.o Makefile

# executables
../$(EXEC_SKYAX): $(MODULES)
	 @echo -n "linking " $@ "... "
	 $(LINKER_SKYAX) $(LINKERFLAGS_SKYAX) -o $@  $(MODULES) $(LIBS_SKYAX)
	 @echo "OK"

.SUFFIXES:  .f90 .o

%.o : %.mod

%.o: %.f90
	@echo -n "compiling "$<" ... "
	$(COMPILER_SKYAX) $(COMPILERFLAGS_SKYAX)  -c $<
	@echo "OK"

