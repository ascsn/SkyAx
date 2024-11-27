#!/bin/bash --login
#SBATCH --job-name=job_170Yb_hfb_1.0d0_1758_13.D0_25_26
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH --time=12:59:59
#SBATCH --mem-per-cpu=4G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/error_job_170Yb_hfb_1.0d0_1758_13.D0_25_26.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/output_job_170Yb_hfb_1.0d0_1758_13.D0_25_26.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo "&main outnam='job_170Yb_hfb_1.0d0_1758_13.D0_25_26',title='test 170Yb',iprint=1000,itmax=1000" > for005_5
echo "/ " >> for005_5
echo "&force " >> for005_5
echo " name='none',pairstab=-0D0,HFBwidth=1.0d0," >> for005_5
echo " t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921," >> for005_5
echo " x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530," >> for005_5
echo " t4=124.633,b4p=34.1116989,power=0.3," >> for005_5
echo " h2m=20.749811,20.721249," >> for005_5
echo " chg=1D0,ex=1,riftm=0D0,zpe=1," >> for005_5
echo " ipair=8,v0prot=1758,v0neut=1520,rho0pr=0.2011591," >> for005_5
echo " pair_cutoff=-25,tregpair=F" >> for005_5
echo "/ " >> for005_5
echo "&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0," >> for005_5
echo "/ " >> for005_5
echo "&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0," >> for005_5
echo "deocc(1)=13.D0,deocc(2)=13.D0,errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0," >> for005_5
echo "tvaryx0=T,inodiag=9999 " >> for005_5
echo "/ " >> for005_5
echo "&constraint" >> for005_5
echo "/ " >> for005_5
echo "&thermal " >> for005_5
echo "/ " >> for005_5

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
# cd ../source_f90/
# make clean
# make ifort_omp
# cd ../slurmjobs
../skyax.ifort.omp < for005_5 > job_170Yb_hfb_1.0d0_1758_13.D0_25_26_output.log

