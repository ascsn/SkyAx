import os
import subprocess
# Setting Epair = 15.5
edf = f"Fayans"
hfbwidths = ["0.0d0","1.0d0"]

ipair = [7, 8]



# Directory to store job scripts
job_script_dir = "job_scripts"
os.makedirs(job_script_dir, exist_ok=True)

def create_and_submit_job_script(job_name, hfbwidth, ipair_val, protNum, neutNum, inpNum):
    # SLURM job script content
    job_script_content = f"""#!/bin/bash --login
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH --time=1-23:59:59
#SBATCH --mem=64G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/101Sn/errs/error_{job_name}.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/101Sn/errs/output_{job_name}.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo \"&main outnam='{job_name}',title='test {nucleus}',iprint=4000,itmax=6000,\" > for005_{inpNum}
echo \" iterbcs=1000, iterhfb=2000, pairenhance=4D0, iternat=10, modnat=1, itest=0 \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&force \" >> for005_{inpNum}
echo \" name='FyDrHFB', tFayans=.TRUE., pairstab=-0.3D0,HFBwidth={hfbwidth},\" >> for005_{inpNum}
echo \"h2m(2)=20.72126,h2m(1)=20.74982,\" >> for005_{inpNum}
echo \"avp=-9.809519866,h1vp=0.6174527519,h2vp=0.2142736685,\" >> for005_{inpNum}
echo \"avm=23.41886322,h1vm=-0.2861245253,h2vm=22.80909400,\" >> for005_{inpNum}
echo \"asp=0.5654841300,hsp=0D0,hgradp=0.4579587200,c0nabj=0.1983338500,c1nabj=4.4008347000E-02,\" >> for005_{inpNum}
echo \"fxi=-4.455644500,h1xi=4.243982800,hgradxi=3.112633900,delfxi=0D0,\" >> for005_{inpNum}
echo \"ipair={ipair_val},ex=1,zpe=1\" >> for005_{inpNum}
echo \"pair_cutoff=-15D0,tcutoff_uvspace=.TRUE.\" >> for005_{inpNum}
echo \"isovectormodel=0 \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&grids tsym=.T.,nz=31,nr=31,dr=0.7D0,dz=0.7D0\" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&static nprot={protNum},nneut={neutNum},e0inv=110D0,x0inv=0.1D0,addrho=0.15D0,addocc=0.3D0,\" >> for005_{inpNum}
echo \"deocc(1)=3.1,deocc(2)=3.1,errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0,\" >> for005_{inpNum}
echo \"tvaryx0=F,inodiag=99999,irepc=100 \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&constraint\" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&thermal \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
# cd ../source_f90/
# make clean
# make ifort_omp
# cd ../slurmjobs
../skyax.ifort.omp < for005/for005_{inpNum} > outs/{job_name}_out.log

"""
    subprocess.run(["cp", f"for005_{inpNum}", f"for005/for005_{inpNum}"])
    # Save the script to a file
    job_script_path = os.path.join(job_script_dir, f"{job_name}.sh")
    with open(job_script_path, "w") as job_script_file:
        job_script_file.write(job_script_content)
    
    # Submit the job script
    subprocess.run(["sbatch", job_script_path])

# Loop through parameter combinations and submit jobs
inpNum = 0
protNum = 50
neutNum = 51
for j in ipair:
    cond = "hfb" if j == 8 else "bcs"
    combined = zip(hfbwidths)
    aNum = protNum + neutNum
    nucleus = f"{aNum}Sn"

    for hfbwidth in hfbwidths:
        job_name = f"job_{nucleus}_{cond}_{hfbwidth}"
        create_and_submit_job_script(job_name, hfbwidth,j, protNum,neutNum, inpNum)
        inpNum +=1

print("All jobs submitted!")
