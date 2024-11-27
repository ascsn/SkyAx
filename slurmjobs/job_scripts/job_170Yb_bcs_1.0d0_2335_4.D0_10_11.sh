#!/bin/bash --login
#SBATCH --job-name=job_170Yb_bcs_1.0d0_2335_4.D0_10_11
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --propagate=CORE
#SBATCH --time=1-23:59:59
#SBATCH --mem-per-cpu=4G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/error_job_170Yb_bcs_1.0d0_2335_4.D0_10_11.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/output_job_170Yb_bcs_1.0d0_2335_4.D0_10_11.out

module purge
module load powertools
module load OpenBLAS/0.3.23-GCC-12.3.0
module load intel/2023a
# module load imkl/2023.1.0

# Prepare input file
echo "&main outnam='job_170Yb_bcs_1.0d0_2335_4.D0_10_11',title='test 170Yb',iprint=1000,itmax=10000," > for005_6
echo " iterbcs=5000, iterhfb=10000 " >> for005_6
echo "/ " >> for005_6
echo "&force " >> for005_6
echo " name='none',pairstab=-0D0,HFBwidth=1.0d0," >> for005_6
echo " t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921," >> for005_6
echo " x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530," >> for005_6
echo " t4=124.633,b4p=34.1116989,power=0.3," >> for005_6
echo " h2m=20.749811,20.721249," >> for005_6
echo " chg=1D0,ex=1,riftm=0D0,zpe=1," >> for005_6
echo " ipair=7,v0prot=2335,v0neut=2081,rho0pr=0.2011591," >> for005_6
echo " pair_cutoff=-10,tregpair=F" >> for005_6
echo "/ " >> for005_6
echo "&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0," >> for005_6
echo "/ " >> for005_6
echo "&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0," >> for005_6
echo "deocc(1)=4.D0,deocc(2)=4.D0,errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0," >> for005_6
echo "tvaryx0=T,inodiag=9999 " >> for005_6
echo "/ " >> for005_6
echo "&constraint" >> for005_6
echo "/ " >> for005_6
echo "&thermal " >> for005_6
echo "/ " >> for005_6

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
# cd ../source_f90/
# make clean
# make ifort_omp
# cd ../slurmjobs
../skyax.ifort.seq < for005_6 > job_170Yb_bcs_1.0d0_2335_4.D0_10_11_out.log

