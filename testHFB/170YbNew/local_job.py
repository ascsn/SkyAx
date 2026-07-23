#!/usr/bin/env python3
"""
Local job script for SkyAx calculations (no SLURM).
Converts from slurmjob.py to run locally on your machine.
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

# Directory to store job scripts and input files
job_script_dir = "job_scripts"
os.makedirs(job_script_dir, exist_ok=True)
os.makedirs("for005", exist_ok=True)
os.makedirs("outs", exist_ok=True)


def create_for005_file(job_name, hfbwidth, deocc_val, ipair_val, run_start, inpNum):
    """Create the for005 input file with specified parameters."""
    lines = [
        f"&main outnam='{job_name}',title='test {nucleus}',iprint=1000,itmax=6000,iterbcs=1000,iterhfb=2000,",
        "      pairenhance=4D0,iternat=10,modnat=1,itest=0,",
        "/",
        "&force",
        " tFayans=.TRUE.,",
        f" name='none',pairstab=-0.3D0,HFBwidth={hfbwidth},",
        " avp=-11.37022366,h1vp=0.5846452878,h2vp=0.5291677813,",
        " avm=48.75238475,h1vm=-0.7554109405,h2vm=65.87341867,",
        " asp=0.5823215900,hsp=0.0,hgradp=0.4560273859,",
        " asm=0.0,hsm=65.87341867,hgradm=0.0,hcoul=0.0,",
        " fxi=-8.332310968,h1xi=5.382384194,hgradxi=12.62794586,delfxi=0.0,",
        " c0nabj=0.1999503224,c1nabj=-3.7238150000E-02,c0tau=0.0,c1tau=0.0,",
        " h2m=20.749811,20.721249,",
        " tcutoff_uvspace=.TRUE.,isovectormodel=2,",
        " chg=1D0,ex=1,zpe=1,tcoulfold=.TRUE.,",
        f" ipair={ipair_val},pair_cutoff=-{run_start}D0,tregpair=F",
        "/",
        "&grids tsym=T,nz=25,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,",
        "/",
        "&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.05D0,addrho=0.15D0,addocc=0.15D0,",
        f"deocc(1)={deocc_val},deocc(2)={deocc_val},errenr=1D-10,errvar=0.05,colerr=1D-9,itpini=0,",
        "tvaryx0=F,inodiag=9999, mmax(1)=50, mmax(2)=50/",
        "&constraint",
        "/",
        "&thermal /",
    ]

    for005_filename = f"for005_{inpNum}"
    for005_path = os.path.join("for005", for005_filename)
    with open(for005_path, "w") as f:
        f.write("\n".join(lines))

    return for005_path


def create_job_script(job_name, for005_path, inpNum):
    """Create a shell script to run a single job locally."""
    script_content = f"""#!/bin/bash
# Job script for {job_name}

echo "Running job: {job_name}"
echo "Input file: {for005_path}"

# Make sure skyax executable exists
if [ ! -f "../../skyax.ifort.omp" ]; then
    echo "Error: ../../skyax.ifort.omp not found. Please build it first."
    exit 1
fi

# Set environment variables for OpenMP
export OMP_STACKSIZE=4500m
ulimit -s unlimited

# Run the calculation
../../skyax.ifort.omp < {for005_path} > outs/{job_name}_out.log

echo "Job {job_name} completed. Output in outs/{job_name}_out.log"
"""

    script_path = os.path.join(job_script_dir, f"{job_name}.sh")
    with open(script_path, "w") as f:
        f.write(script_content)
    os.chmod(script_path, 0o755)

    return script_path


def main():
    """Generate all job configurations."""
    print(f"Generating jobs for {nucleus} ({edf})...")
    print(f"Job scripts: {job_script_dir}/")
    print(f"Input files: for005/")
    print(f"Output logs: outs/")

    inpNum = 0
    job_scripts = []

    for j in ipair:
        for val in deoccList:
            cond = "hfb" if j == 8 else "bcs"
            combined = zip(hfbwidths, runList)

            for hfbwidth, run_val in combined:
                run_start, run_end = run_val
                job_name = f"job_{nucleus}_{cond}_{hfbwidth}_{val}_{run_start}_{run_end}"

                # Create for005 input file
                for005_path = create_for005_file(job_name, hfbwidth, val, j, run_start, inpNum)

                # Create shell script
                script_path = create_job_script(job_name, for005_path, inpNum)
                job_scripts.append((job_name, script_path))

                print(f"  Created: {job_name}")
                inpNum += 1

    print(f"\nGenerated {len(job_scripts)} job configurations!")
    print(f"\nTo run all jobs sequentially:")
    print(f"  for script in {job_script_dir}/*.sh; do bash \"$script\"; done")
    print(f"\nOr run a single job:")
    if job_scripts:
        print(f"  bash {job_scripts[0][1]}")
    print(f"\nTo run jobs in parallel (using GNU parallel or xargs):")
    print(f"  find {job_script_dir} -name '*.sh' | xargs -P 4 -I {{}} bash {{}}")


if __name__ == "__main__":
    main()
