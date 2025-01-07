#!/bin/bash --login
#SBATCH --job-name=job_101Sn_hfb_1.0d0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH --time=1-23:59:59
#SBATCH --mem=64G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/101Sn/errs/error_job_101Sn_hfb_1.0d0.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/101Sn/errs/output_job_101Sn_hfb_1.0d0.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo "&main outnam='job_101Sn_hfb_1.0d0',title='test 101Sn',iprint=4000,itmax=6000," > for005_3
echo " iterbcs=1000, iterhfb=2000, pairenhance=4D0, iternat=10, modnat=1, itest=0 " >> for005_3
echo "/ " >> for005_3
echo "&force " >> for005_3
echo " name='FyDrHFB', tFayans=.TRUE., pairstab=-0.3D0,HFBwidth=1.0d0," >> for005_3
echo "h2m(2)=20.72126,h2m(1)=20.74982," >> for005_3
echo "avp=-9.809519866,h1vp=0.6174527519,h2vp=0.2142736685," >> for005_3
echo "avm=23.41886322,h1vm=-0.2861245253,h2vm=22.80909400," >> for005_3
echo "asp=0.5654841300,hsp=0D0,hgradp=0.4579587200,c0nabj=0.1983338500,c1nabj=4.4008347000E-02," >> for005_3
echo "fxi=-4.455644500,h1xi=4.243982800,hgradxi=3.112633900,delfxi=0D0," >> for005_3
echo "ipair=8,ex=1,zpe=1" >> for005_3
echo "pair_cutoff=-15D0,tcutoff_uvspace=.TRUE." >> for005_3
echo "isovectormodel=0 " >> for005_3
echo "/ " >> for005_3
echo "&grids tsym=.T.,nz=31,nr=31,dr=0.7D0,dz=0.7D0" >> for005_3
echo "/ " >> for005_3
echo "&static nprot=50,nneut=51,e0inv=110D0,x0inv=0.1D0,addrho=0.15D0,addocc=0.3D0," >> for005_3
echo "deocc(1)=3.1,deocc(2)=3.1,errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0," >> for005_3
echo "tvaryx0=F,inodiag=99999,irepc=100 " >> for005_3
echo "/ " >> for005_3
echo "&constraint" >> for005_3
echo "/ " >> for005_3
echo "&thermal " >> for005_3
echo "/ " >> for005_3

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
# cd ../source_f90/
# make clean
# make ifort_omp
# cd ../slurmjobs
../skyax.ifort.omp < for005/for005_3 > outs/job_101Sn_hfb_1.0d0_out.log

