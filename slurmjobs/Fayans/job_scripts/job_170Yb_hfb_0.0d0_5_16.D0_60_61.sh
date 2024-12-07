#!/bin/bash --login
#SBATCH --job-name=job_170Yb_hfb_0.0d0_5_16.D0_60_61
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --propagate=CORE
#SBATCH --time=6-23:59:59
##SBATCH --mem=64G
##SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/Fayans/errs/error_job_170Yb_hfb_0.0d0_5_16.D0_60_61.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/Fayans/outs/output_job_170Yb_hfb_0.0d0_5_16.D0_60_61.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo "&main outnam='job_170Yb_hfb_0.0d0_5_16.D0_60_61',title='test 170Yb',iprint=4000,itmax=4000,itest=0, tpri_2d=F," > for005_9
echo " iterbcs=500, iterhfb=1000,pairenhance=1D0,iternat=10,modnat=1, " >> for005_9
echo "/ " >> for005_9
echo "&force " >> for005_9
echo " name='none', tFayans=.TRUE., ipair=8, pairstab=-0.3D0, ex=1, zpe=1, HFBwidth=0.0d0," >> for005_9
echo " h2m(2)=20.72126, h2m(1)=20.74982, " >> for005_9
echo " avp=-9.809519866,h1vp=0.6174527519,h2vp=0.2142736685," >> for005_9
echo " avm=23.41886322,h1vm=-0.2861245253,h2vm=22.80909400," >> for005_9
echo " asp=0.5654841300,hsp=0.0D0,hgradp=0.4579587200,c0nabj=0.1983338500,c1nabj=4.4008347000E-02," >> for005_9
echo " fxi=-5,h1xi=4.243982800,hgradxi=3.112633900,delfxi=0D0," >> for005_9
echo " pair_cutoff=-60,tcutoff_uvspace=.TRUE." >> for005_9
echo " isovectormodel=0 " >> for005_9

echo "/ " >> for005_9
echo "&grids tsym=.T.,nz=31,nr=31,dr=0.7D0,dz=0.7D0" >> for005_9
echo "/ " >> for005_9
echo "&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.1D0,addrho=0.15D0,addocc=0.3D0," >> for005_9
echo "deocc(1)=16.D0,deocc(2)=16.D0,errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0," >> for005_9
echo "tvaryx0=F,inodiag=9999,irepc=100 " >> for005_9
echo "/ " >> for005_9
echo "&constraint" >> for005_9
echo "/ " >> for005_9
echo "&thermal " >> for005_9
echo "/ " >> for005_9

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../../skyax.ifort.omp < for005_9 > job_170Yb_hfb_0.0d0_5_16.D0_60_61_out.log

