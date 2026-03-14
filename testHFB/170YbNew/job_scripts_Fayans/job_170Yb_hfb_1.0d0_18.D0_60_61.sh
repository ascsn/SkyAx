#!/bin/bash
#SBATCH --job-name=job_170Yb_hfb_1.0d0_18.D0_60_61
#SBATCH --output=output_job_170Yb_hfb_1.0d0_18.D0_60_61.out
#SBATCH --error=error_job_170Yb_hfb_1.0d0_18.D0_60_61.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH -A frib-nodes
#SBATCH --time=1-23:59:59
#SBATCH --mem=64G

# Load any necessary modules here
module purge
module load powertools
module load intel/2023a
# module load intel
# module load mkl

# Prepare input file (for005)
cat > for005_49 << 'EOF'
&main outnam='job_170Yb_hfb_1.0d0_18.D0_60_61',title='test 170Yb',iprint=1000,itmax=6000,iterbcs=1000,iterhfb=2000,
      pairenhance=4D0,iternat=10,modnat=1,itest=0,
/
&force
 tFayans=.TRUE.,
 name='none',pairstab=-0.3D0,HFBwidth=1.0d0,
 avp=-11.37022366,h1vp=0.5846452878,h2vp=0.5291677813,
 avm=48.75238475,h1vm=-0.7554109405,h2vm=65.87341867,
 asp=0.5823215900,hsp=0.0,hgradp=0.4560273859,
 asm=0.0,hsm=65.87341867,hgradm=0.0,hcoul=0.0,
 fxi=-8.332310968,h1xi=5.382384194,hgradxi=12.62794586,delfxi=0.0,
 c0nabj=0.1999503224,c1nabj=-3.7238150000E-02,c0tau=0.0,c1tau=0.0,
 h2m=20.749811,20.721249,
 tcutoff_uvspace=.TRUE.,isovectormodel=2,
 chg=1D0,ex=1,zpe=1,tcoulfold=.TRUE.,
 ipair=8,pair_cutoff=-60D0,tregpair=F
/
&grids tsym=T,nz=25,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,
/
&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.05D0,addrho=0.15D0,addocc=0.15D0,
deocc(1)=18.D0,deocc(2)=18.D0,errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0,
tvaryx0=F,inodiag=9999, mmax(1)=50, mmax(2)=50/
&constraint
/
&thermal /
EOF

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../../skyax.ifort.omp < for005/for005_49 > job_170Yb_hfb_1.0d0_18.D0_60_61_output.log

echo "Job job_170Yb_hfb_1.0d0_18.D0_60_61 completed"
