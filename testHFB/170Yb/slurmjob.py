#!/usr/bin/env python3
"""
SLURM job script for SkyAx calculations.
Uses Fayans functional for 170Yb nucleus.
"""
import os
import subprocess
import sys

# Setting Epair = 15.5
nucleus = "170Yb"
edf = "FyDr"
# hfbwidths = ["0.0d0", "0.0d0", "0.0d0", "0.0d0", "0.0d0"]
hfbwidths = ["1.0d0", "1.0d0", "1.0d0", "1.0d0", "1.0d0"]

ipair = [8]
deoccList = ["2.D0", "4.D0", "6.D0", "8.D0", "10.D0", "12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]
runList = [[5, 6], [10, 11], [20, 21], [40, 41], [60, 61]]

# Directory to store job scripts
job_script_dir = "job_scripts_Fayans"
os.makedirs(job_script_dir, exist_ok=True)


def create_and_submit_job_script(job_name, hfbwidth, deocc_val, ipair_val, run_start, run_end, inpNum):
    """Create and submit a SLURM job script with Fayans functional."""
    # SLURM job script content
    job_script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=output_{job_name}.out
#SBATCH --error=error_{job_name}.err
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
cat > for005_{inpNum} << 'EOF'
&main outnam='{job_name}',title='test {nucleus}',iprint=1000,itmax=6000,iterbcs=1000,iterhfb=2000,
      pairenhance=4D0,iternat=10,modnat=1,itest=0,
/
&force
 tFayans=.TRUE.,
 name='none',pairstab=-0.3D0,HFBwidth={hfbwidth},
 avp=-11.37022366,h1vp=0.5846452878,h2vp=0.5291677813,
 avm=48.75238475,h1vm=-0.7554109405,h2vm=65.87341867,
 asp=0.5823215900,hsp=0.0,hgradp=0.4560273859,
 asm=0.0,hsm=65.87341867,hgradm=0.0,hcoul=0.0,
 fxi=-8.332310968,h1xi=5.382384194,hgradxi=12.62794586,delfxi=0.0,
 c0nabj=0.1999503224,c1nabj=-3.7238150000E-02,c0tau=0.0,c1tau=0.0,
 h2m=20.749811,20.721249,
 tcutoff_uvspace=.TRUE.,isovectormodel=2,
 chg=1D0,ex=1,zpe=1,tcoulfold=.TRUE.,
 ipair={ipair_val},pair_cutoff=-{run_start}D0,tregpair=F
/
&grids tsym=T,nz=25,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,
/
&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.05D0,addrho=0.15D0,addocc=0.15D0,
deocc(1)={deocc_val},deocc(2)={deocc_val},errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0,
tvaryx0=F,inodiag=9999, mmax(1)=50, mmax(2)=50/
&constraint
/
&thermal /
EOF

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../../skyax.ifort.omp < for005/for005_{inpNum} > {job_name}_output.log

echo "Job {job_name} completed"
"""

    # Save the script to a file
    job_script_path = os.path.join(job_script_dir, f"{job_name}.sh")
    with open(job_script_path, "w") as job_script_file:
        job_script_file.write(job_script_content)
    
    # Submit the job script
    try:
        result = subprocess.run(["sbatch", job_script_path], 
                              capture_output=True, 
                              text=True, 
                              check=True)
        print(f"Submitted: {job_name} - {result.stdout.strip()}")
    except subprocess.CalledProcessError as e:
        print(f"Error submitting {job_name}: {e.stderr}")
    except FileNotFoundError:
        print(f"Error: sbatch command not found. Are you on a SLURM cluster?")
        sys.exit(1)


def main():
    """Generate and submit all job configurations."""
    print(f"Submitting SLURM jobs for {nucleus} ({edf})...")
    
    job_count = 0
    
    # Loop through parameter combinations and submit jobs
    for j in ipair:
        for val in deoccList:
            cond = "hfb" if j == 8 else "bcs"
            combined = zip(hfbwidths, runList)

            for hfbwidth, run_val in combined:
                run_start, run_end = run_val
                job_name = f"job_{nucleus}_{cond}_{hfbwidth}_{val}_{run_start}_{run_end}"
                
                create_and_submit_job_script(job_name, hfbwidth, val, j, run_start, run_end, job_count)
                job_count += 1

    print(f"\nAll {job_count} jobs submitted!")
    print(f"Check job status with: squeue -u $USER")
    print(f"Job scripts saved in: {job_script_dir}/")


if __name__ == "__main__":
    main()
