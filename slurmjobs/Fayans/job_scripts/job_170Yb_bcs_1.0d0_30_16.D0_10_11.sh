#!/bin/bash --login
#SBATCH --job-name=job_170Yb_bcs_1.0d0_30_16.D0_10_11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH --time=1-23:59:59
#SBATCH --mem=64G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/Fayans/errs/error_job_170Yb_bcs_1.0d0_30_16.D0_10_11.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/Fayans/outs/output_job_170Yb_bcs_1.0d0_30_16.D0_10_11.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo "&main outnam='job_170Yb_bcs_1.0d0_30_16.D0_10_11',title='test 170Yb',iprint=1000,itmax=4000," > for005_1
echo " iterbcs=500, iterhfb=1000,pairenhance=1D0,iternat=10,modnat=1, " >> for005_1
echo "/ " >> for005_1
echo "&force " >> for005_1
echo " name='none', tFayans=.TRUE., ipair=7, pairstab=-0D0, ex=1, zpe=1, HFBwidth=1.0d0," >> for005_1
echo " h2m(2)=20.72126, h2m(1)=20.74982, " >> for005_1
echo " avp=-10.05903239,h1vp=0.6079313192,h2vp=0.2776243076," >> for005_1
echo " avm=47.20760687,h1vm=-0.3746093562,h2vm=54.42781346," >> for005_1
echo " asp=0.5283255600,hsp=0.0,hgradp=0.3274157900,c0nabj=0.1788187900,c1nabj=6.8951346000E-02," >> for005_1
echo " fxi=-30,h1xi=4.443955000,hgradxi=2.648612100,delfxi=-0.1881396400," >> for005_1
echo " pair_cutoff=-10,tcutoff_uvspace=.TRUE." >> for005_1
echo " isovectormodel=2 " >> for005_1

echo "/ " >> for005_1
echo "&grids tsym=.T.,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0," >> for005_1
echo "/ " >> for005_1
echo "&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.1D0,addrho=0.15D0,addocc=0.3D0," >> for005_1
echo "deocc(1)=16.D0,deocc(2)=16.D0,errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0," >> for005_1
echo "tvaryx0=F,inodiag=99999,irepc=100 " >> for005_1
echo "/ " >> for005_1
echo "&constraint" >> for005_1
echo "/ " >> for005_1
echo "&thermal " >> for005_1
echo "/ " >> for005_1

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../../skyax.ifort.omp < for005_1 > job_170Yb_bcs_1.0d0_30_16.D0_10_11_out.log

