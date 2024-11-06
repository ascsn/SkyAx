import os
import subprocess

nucleus = f"120Sn"
edf = f"SVbas"
hfbwidths = ["1.0d0", "1.0d0", "1.0d0", "1.0d0", "1.0d0", "1.0d0"]
pairStrengthsBCS = ["537", "478", "451", "438", "433", "432"]
# pairStrengthsHFB = ["628", "462", "408", "374", "328.3", "389.2"] # HFBWIDTH = 0
pairStrengthsHFB = ["1734", "1470", "1328", "1284", "1270", "1268"] # HFBWIDTH = 1.0
# pairStrengthsHFB = ["1736", "1470", "1328", "1284", "1270", "1268"]
ipair = [8]
deoccList = ["6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]   
# deoccList = ["16.D0"]    
runList = [[10,11], [20,21], [40,41], [60,61],[80,81],[100,101]]

# Directory to store job scripts
job_script_dir = "job_scripts"
os.makedirs(job_script_dir, exist_ok=True)

def create_and_submit_job_script(job_name, hfbwidth, pairStrength, deocc_val, ipair_val, run_start, run_end, inpNum):
    # SLURM job script content
    job_script_content = f"""#!/bin/bash --login
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH --time=12:59:59
#SBATCH --mem-per-cpu=4G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/error_{job_name}.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/output_{job_name}.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo \"&main outnam='{job_name}',title='test {nucleus}',iprint=1000,itmax=1000\" > for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&force \" >> for005_{inpNum}
echo \" name='none',pairstab=-0D0,HFBwidth={hfbwidth},\" >> for005_{inpNum}
echo \" t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921,\" >> for005_{inpNum}
echo \" x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530,\" >> for005_{inpNum}
echo \" t4=124.633,b4p=34.1116989,power=0.3,\" >> for005_{inpNum}
echo \" h2m=20.749811,20.721249,\" >> for005_{inpNum}
echo \" chg=1D0,ex=1,riftm=0D0,zpe=1,\" >> for005_{inpNum}
echo \" ipair={ipair_val},v0prot=10,v0neut={pairStrength},rho0pr=0.2011591,\" >> for005_{inpNum}
echo \" pair_cutoff=-{run_start},tregpair=F\" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\" >> for005_{inpNum}
echo \"/ \" >> for005_{inpNum}
echo \"&static nprot=50,nneut=70,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0,\" >> for005_{inpNum}
echo \"deocc(1)={deocc_val},deocc(2)={deocc_val},errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0,\" >> for005_{inpNum}
echo \"tvaryx0=T,inodiag=9999 \" >> for005_{inpNum}
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
../skyax.ifort.omp < for005_{inpNum} > {job_name}_output.log

# Post-process output
module load Conda
conda activate slpy
python3 process_results.py --output {job_name}_output.log --deocc {deocc_val} --run_val {run_start} --pair_strength {pairStrength} --ipair {ipair_val}"
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
            job_name = f"job_{cond}_{hfbwidth}_{pair_strength}_{val}_{run_start}_{run_end}"
            create_and_submit_job_script(job_name, hfbwidth, pair_strength, val, j, run_start, run_end,inpNum)
            inpNum +=1

print("All jobs submitted!")
