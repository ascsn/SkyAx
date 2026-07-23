import os
import subprocess

nucleus = "120Sn"
edf = "SVbas"
hfbwidths = ["0.0d0"]
pairStrengthsBCS = ["537", "478", "451", "438", "433", "432"]
pairStrengthsHFB = ["500", "460", "450", "438", "433", "432"]
ipair = [8]
deoccList = ["6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]      
runList = [[10, 11], [20,21], [40,41], [60,61],[80,81],[100,101]]

# Directory to store job scripts
job_script_dir = "job_scripts"
os.makedirs(job_script_dir, exist_ok=True)

def create_and_submit_job_script(job_name, hfbwidth, pairStrength, deocc_val, ipair_val, run_start, run_end):
    # SLURM job script content
    job_script_content = f"""#!/bin/bash
#SBATCH --job-name={job_name}
#SBATCH --output=output_{job_name}.out
#SBATCH --error=error_{job_name}.err
#SBATCH --ntasks=1
#SBATCH --time=00:30:00
#SBATCH --mem=4G

module load your_module_name # Load any necessary modules here

# Prepare input file
echo "&main outnam='tst',title='test {nucleus}',iprint=1000,itmax=1000\n/" > for005
echo "&force name='none',pairstab=-0D0,HFBwidth={hfbwidth}," >> for005
echo " t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921," >> for005
echo " x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530," >> for005
echo " t4=124.633,b4p=34.1116989,power=0.3," >> for005
echo " h2m=20.749811,20.721249," >> for005
echo " chg=1D0,ex=1,riftm=0D0,zpe=1," >> for005
echo " ipair={ipair_val},v0prot=10,v0neut={pairStrength},rho0pr=0.2011591," >> for005
echo " pair_cutoff=-{run_start},tregpair=F\n/" >> for005
echo "&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\n/" >> for005
echo "&static nprot=50,nneut=70,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0," >> for005
echo "deocc(1)={deocc_val},deocc(2)={deocc_val},errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0,\n/" >> for005
echo "&constraint\n/" >> for005
echo "&thermal /\n" >> for005

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../skyax.ifort.seq < for005 > {job_name}_output.log

# Post-process output
python3 process_results.py --output {job_name}_output.log --deocc {deocc_val} --run_val {run_start} --pair_strength {pairStrength} --ipair {ipair_val}
"""

    # Save the script to a file
    job_script_path = os.path.join(job_script_dir, f"{job_name}.sh")
    with open(job_script_path, "w") as job_script_file:
        job_script_file.write(job_script_content)
    
    # Submit the job script
    subprocess.run(["sbatch", job_script_path])

# Loop through parameter combinations and submit jobs
for j in ipair:
    for val in deoccList:
        cond = "hfb" if j == 8 else "bcs"
        pair_strength_list = pairStrengthsHFB if j == 8 else pairStrengthsBCS
        combined = zip(hfbwidths, pair_strength_list, runList)

        for hfbwidth, pair_strength, run_val in combined:
            run_start, run_end = run_val
            job_name = f"job_{cond}_{hfbwidth}_{pair_strength}_{val}_{run_start}_{run_end}"
            create_and_submit_job_script(job_name, hfbwidth, pair_strength, val, j, run_start, run_end)

print("All jobs submitted!")
