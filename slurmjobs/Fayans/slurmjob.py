import os
import subprocess
# Setting Epair = 15.5
nucleus = f"170Yb"
edf = f"FyDr"
# hfbwidths = ["0.0d0","0.0d0", "0.0d0", "0.0d0", "0.0d0"]
hfbwidths = ["1.0d0","1.0d0", "1.0d0", "1.0d0", "1.0d0"]


# pairStrengthsBCS = ["4.7", "4.59", "4.473", "4.394", "4.15"] # HFBWIDTH = 0
# pairStrengthsBCS = ["7", "6.75", "6.5", "6.6", "6.62"] # HFBWIDTH = 1

pairStrengthsHFB = ["8", "8", "6.5", "6.6", "6.6"] # HFBWIDTH = 1

ipair = [8]
# deoccList = ["2.D0","4.D0","6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]   
# deoccList = ["12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]    
deoccList = ["16.D0"]
runList = [[5,6], [10,11], [20,21], [40,41], [60,61]]
# runList = [[40,41], [60,61]]


# Directory to store job scripts
job_script_dir = "job_scripts"
os.makedirs(job_script_dir, exist_ok=True)

def create_and_submit_job_script(job_name, hfbwidth, pairStrength, deocc_val, ipair_val, run_start, run_end, inpNum):
    # SLURM job script content
    job_script_content = f"""#!/bin/bash --login
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --propagate=CORE
#SBATCH --time=1-23:59:59
##SBATCH --mem=64G
##SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/Fayans/errs/error_{job_name}.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/Fayans/outs/output_{job_name}.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo \"&main outnam='{job_name}',title='test {nucleus}',iprint=4000,itmax=4000,itest=0, tpri_2d=F,\" > for005_{inpNum}
echo \" iterbcs=500, iterhfb=1000,pairenhance=1D0,iternat=10,modnat=1, \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&force \" >> for005_{inpNum}
echo \" name='none', tFayans=.TRUE., ipair={ipair_val}, pairstab=-0.3D0, ex=1, zpe=1, HFBwidth={hfbwidth},\" >> for005_{inpNum}
echo \" h2m(2)=20.72126, h2m(1)=20.74982, \" >> for005_{inpNum}
echo \" avp=-9.809519866,h1vp=0.6174527519,h2vp=0.2142736685,\" >> for005_{inpNum}
echo \" avm=23.41886322,h1vm=-0.2861245253,h2vm=22.80909400,\" >> for005_{inpNum}
echo \" asp=0.5654841300,hsp=0.0D0,hgradp=0.4579587200,c0nabj=0.1983338500,c1nabj=4.4008347000E-02,\" >> for005_{inpNum}
echo \" fxi=-{pairStrength},h1xi=4.243982800,hgradxi=3.112633900,delfxi=0D0,\" >> for005_{inpNum}
echo \" pair_cutoff=-{run_start},tcutoff_uvspace=.TRUE.\" >> for005_{inpNum}
echo \" isovectormodel=0 \" >> for005_{inpNum}

echo \"/ \" >> for005_{inpNum}
echo \"&grids tsym=.T.,nz=31,nr=31,dr=0.7D0,dz=0.7D0\" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.1D0,addrho=0.15D0,addocc=0.3D0,\" >> for005_{inpNum}
echo \"deocc(1)={deocc_val},deocc(2)={deocc_val},errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0,\" >> for005_{inpNum}
echo \"tvaryx0=F,inodiag=9999,irepc=100 \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&constraint\" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&thermal \" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../../skyax.ifort.omp < for005_{inpNum} > {job_name}_out.log

"""

    # Save the script to a file
    job_script_path = os.path.join(job_script_dir, f"{job_name}.sh")
    with open(job_script_path, "w") as job_script_file:
        job_script_file.write(job_script_content)
    
    # Submit the job script
    subprocess.run(["sbatch", job_script_path])

# Loop through parameter combinations and submit jobs
inpNum = 0
for j in ipair:
    for val in deoccList:
        cond = "hfb" if j == 8 else "bcs"
        pair_strength_list = pairStrengthsHFB if j == 8 else pairStrengthsBCS
        combined = zip(hfbwidths, pair_strength_list, runList)

        for hfbwidth, pair_strength, run_val in combined:
            run_start, run_end = run_val
            job_name = f"job_{nucleus}_{cond}_{hfbwidth}_{pair_strength}_{val}_{run_start}_{run_end}"
            create_and_submit_job_script(job_name, hfbwidth, pair_strength, val, j, run_start, run_end,inpNum)
            inpNum +=1

print("All jobs submitted!")
