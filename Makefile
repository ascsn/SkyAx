q/SHELL = /bin/bash

# =============================================================================
#  SkyAx — Axial Skyrme-Hartree-Fock Code  |  Build System
# =============================================================================

.PHONY: all info params clean clean-exec \
        gf gf_omp gf_d \
        ifort ifort_omp ifort_dbg ifort_domp \
        ifx ifx_omp \
        makeit

# -----------------------------------------------------------------------------
#  Compiler selection
# -----------------------------------------------------------------------------

gf gf_omp gf_d                         : COMPILER_SKYAX = gfortran
ifort ifort_omp ifort_dbg ifort_domp   : COMPILER_SKYAX = ifort
ifx ifx_omp                            : COMPILER_SKYAX = ifx

export COMPILER_SKYAX

# -----------------------------------------------------------------------------
#  Compiler flags
# -----------------------------------------------------------------------------

gf         : COMPILERFLAGS_SKYAX = -Ofast -march=native -mcmodel=large
gf_omp     : COMPILERFLAGS_SKYAX = -Ofast -march=native -fopenmp
gf_d       : COMPILERFLAGS_SKYAX = -g -fbacktrace -Wall -fbounds-check \
                                    -finit-real=NaN -ffpe-trap=invalid,zero,overflow

ifort      : COMPILERFLAGS_SKYAX = -O3 -xhost -mcmodel=large
ifort_omp  : COMPILERFLAGS_SKYAX = -O3 -xhost -qopenmp -mcmodel=large
ifort_dbg  : COMPILERFLAGS_SKYAX = -g -traceback -CB -mcmodel=large -check uninit
ifort_domp : COMPILERFLAGS_SKYAX = -g -traceback -CB -qopenmp -mcmodel=large -check uninit

ifx        : COMPILERFLAGS_SKYAX = -O3 -xhost -mcmodel=large
ifx_omp    : COMPILERFLAGS_SKYAX = -Ofast -march=native -fopenmp

export COMPILERFLAGS_SKYAX

# -----------------------------------------------------------------------------
#  Executable names
# -----------------------------------------------------------------------------

gf         : EXEC_SKYAX = skyax.gfortran.seq
gf_omp     : EXEC_SKYAX = skyax.gfortran.omp
gf_d       : EXEC_SKYAX = skyax.gfortran.dbg

ifort      : EXEC_SKYAX = skyax.ifort.seq
ifort_omp  : EXEC_SKYAX = skyax.ifort.omp
ifort_dbg  : EXEC_SKYAX = skyax.ifort.dbg
ifort_domp : EXEC_SKYAX = skyax.ifort.domp

ifx        : EXEC_SKYAX = skyax.ifx.seq
ifx_omp    : EXEC_SKYAX = skyax.ifx.omp

export EXEC_SKYAX

# -----------------------------------------------------------------------------
#  Libraries
#  Uncomment the appropriate line for your environment:
#    - Generic LAPACK/BLAS (most systems / ICER, MSU):
#        ALL_TARGETS_LIBS = -llapack -lblas
#    - MKL via ifort/ifx (Intel compilers):
#        INTEL_TARGETS_LIBS = -qmkl
# -----------------------------------------------------------------------------

# Generic: uncomment for gfortran + system LAPACK/BLAS
# gf gf_d gf_omp : LIBS_SKYAX = -llapack -lblas

# Intel MKL (ifort / ifx)
ifort ifort_omp ifort_dbg ifort_domp ifx ifx_omp : LIBS_SKYAX = -qmkl

export LIBS_SKYAX

# -----------------------------------------------------------------------------
#  Top-level targets
# -----------------------------------------------------------------------------

all: info params

info:
	@echo "======================================================"
	@echo "    SkyAx  —  Axial Skyrme-Hartree-Fock Code         "
	@echo "              >>> version 1.0.0 <<<                   "
	@echo "======================================================"
	@rm -f compile.log

params:
	@echo "  Syntax:"
	@echo "    make gf          gfortran, sequential"
	@echo "    make gf_omp      gfortran + OpenMP"
	@echo "    make gf_d        gfortran, debug"
	@echo "    make ifort       ifort, sequential"
	@echo "    make ifort_omp   ifort + OpenMP"
	@echo "    make ifort_dbg   ifort, debug"
	@echo "    make ifort_domp  ifort, debug + OpenMP"
	@echo "    make ifx         ifx, sequential"
	@echo "    make ifx_omp     ifx + OpenMP"
	@echo
	@echo "    make clean       remove objects and module files"
	@echo "    make clean-exec  remove objects, modules, and executables"
	@echo
	@echo "    make all         print version info and this help screen"
	@echo

clean:
	@rm -f *.o *.mod
	@echo "  Cleaned objects and module files."

clean-exec: clean
	@rm -f ../skyax.*.*
	@echo "  Cleaned executables."

# Restart make with the proper environment exported from the target above
gf gf_d gf_omp ifort ifort_omp ifort_dbg ifort_domp ifx ifx_omp:
	@$(MAKE) --no-print-directory -f $(firstword $(MAKEFILE_LIST)) makeit

# =============================================================================
#  After restart — build rules
# =============================================================================

LINKER_SKYAX      = $(COMPILER_SKYAX)
LINKERFLAGS_SKYAX = $(COMPILERFLAGS_SKYAX)

MODULES = \
    m_paramr.o      \
    m_dynpar.o      \
    m_levels.o      \
    m_grid.o        \
    m_wavefr.o      \
    m_fields.o      \
    m_constr.o      \
    m_observ.o      \
    m_coulomb.o     \
    m_staticstep.o  \
    m_pairing.o     \
    m_util.o        \
    m_init.o        \
    m_hdiag.o       \
    m_steputil.o    \
    user_interface.o \
    skyax.o         \
    m_fyfunc_2D.o   \
    hfb-natural.o   \
    m_blocking.o

# -----------------------------------------------------------------------------
#  Module dependencies
# -----------------------------------------------------------------------------

makeit: ../$(EXEC_SKYAX)

m_paramr.o      : Makefile

m_grid.o        : m_paramr.o Makefile

m_dynpar.o      : m_paramr.o forces.data Makefile

m_levels.o      : m_paramr.o Makefile

m_fields.o      : m_paramr.o m_grid.o Makefile

m_wavefr.o      : m_paramr.o m_grid.o m_levels.o Makefile

m_constr.o      : m_paramr.o m_dynpar.o m_levels.o m_grid.o m_fields.o Makefile

m_hdiag.o       : m_grid.o m_levels.o m_paramr.o m_wavefr.o Makefile

m_observ.o      : m_paramr.o m_grid.o m_constr.o m_dynpar.o m_wavefr.o \
                  m_levels.o m_fields.o Makefile

m_coulomb.o     : m_paramr.o m_fields.o m_wavefr.o m_levels.o m_grid.o \
                  m_constr.o m_dynpar.o m_observ.o Makefile

m_steputil.o    : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
                  m_levels.o m_observ.o m_paramr.o m_wavefr.o lagpar.f90 Makefile

m_fyfunc_2D.o   : m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
                  m_observ.o m_paramr.o m_wavefr.o Makefile

m_util.o        : m_dynpar.o m_fields.o m_grid.o m_levels.o m_observ.o \
                  m_paramr.o m_wavefr.o Makefile

m_pairing.o     : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
                  m_levels.o m_observ.o m_paramr.o m_wavefr.o m_fyfunc_2D.o Makefile

hfb-natural.o   : m_pairing.o m_dynpar.o m_fields.o m_grid.o \
                  m_levels.o m_observ.o m_paramr.o m_wavefr.o m_fyfunc_2D.o Makefile

m_staticstep.o  : m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
                  m_grid.o m_levels.o m_observ.o m_paramr.o m_wavefr.o \
                  m_hdiag.o m_pairing.o m_util.o m_steputil.o m_init.o Makefile

m_init.o        : m_constr.o m_coulomb.o m_dynpar.o m_fields.o m_grid.o \
                  m_levels.o m_observ.o m_paramr.o m_wavefr.o m_pairing.o \
                  m_util.o m_steputil.o Makefile

m_blocking.o    : m_dynpar.o m_grid.o m_levels.o m_paramr.o m_wavefr.o \
                  m_pairing.o m_staticstep.o m_util.o m_init.o Makefile

user_interface.o : Makefile

skyax.o         : m_constr.o m_coulomb.o m_dynpar.o m_fields.o \
                  m_grid.o m_init.o m_levels.o m_observ.o m_paramr.o \
                  m_wavefr.o m_staticstep.o m_util.o m_blocking.o Makefile

# -----------------------------------------------------------------------------
#  Link
# -----------------------------------------------------------------------------

../$(EXEC_SKYAX): $(MODULES)
	@echo -n "Linking $@ ... "
	@$(LINKER_SKYAX) $(LINKERFLAGS_SKYAX) -o $@ $(MODULES) $(LIBS_SKYAX)
	@echo "OK"

# -----------------------------------------------------------------------------
#  Compile rule  (.f90 → .o)
# -----------------------------------------------------------------------------

.SUFFIXES:
.SUFFIXES: .f90 .o

%.o: %.f90
	@echo -n "Compiling $< ... "
	@$(COMPILER_SKYAX) $(COMPILERFLAGS_SKYAX) -c $<
	@echo "OK"
