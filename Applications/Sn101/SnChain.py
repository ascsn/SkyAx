""" Calculates all properties for a given range of Z and N
Script written for HPC at MSU. Creates jobs for all nuclei and runs them.
For collecting the results (ensurf and other files) use transfer.py script.
Author: Sudhanva Lalit
Last updated: 2025-01-30
"""

import os
import subprocess
import numpy as np

# Directory to store job scripts
job_script_dir = "job_scripts"
os.makedirs(job_script_dir, exist_ok=True)


def write_for005(edf, Z, N, hfbwidth, pairType, betaInp, index):
    """Writes the for005 file for a given nucleus and EDF"""

    if edf == "FyDr":
        with open(f"for005_{index}", "w") as f:
            f.write(
                f"&main outnam='{index}',title='Fy(Dr,HFB)',iprint=1000,itmax=8000,iterbcs=300,iterhfb=600,\n"
            )
            f.write("      pairenhance=1.0D0,iternat=0,modnat=99999,itest=0,\n")
            f.write("/\n")
            f.write("&force\n")
            f.write(" tFayans=.TRUE.,\n")
            f.write(" name='none',pairstab=-0.3D0,HFBwidth=0.0d0,\n")
            f.write("  avp=-9.809519866,h1vp=0.6174527519,h2vp=0.2142736685,\n")
            f.write("  avm=23.41886322,h1vm=-0.2861245253,h2vm=22.80909400,\n")
            f.write("  asp=0.5654841300,hsp=0.0,hgradp=0.4579587200,\n")
            f.write("  asm=0.0,hsm=22.80909400,hgradm=0.0,hcoul=0.0,\n")
            f.write(
                "  fxi=-4.455644500,h1xi=4.243982800,hgradxi=3.112633900,delfxi=0.0,\n"
            )
            f.write(
                "  c0nabj=0.1983338500,c1nabj=4.4008347000E-02,c0tau=0.0,c1tau=0.0,\n"
            )
            f.write(" h2m=20.749811,20.721249,\n")
            f.write("  tcutoff_uvspace=.TRUE.,isovectormodel=2,\n")
            f.write(" chg=1D0,ex=1,zpe=1,tcoulfold=.TRUE.,\n")
            f.write(f" ipair={pairType},pair_cutoff=-15D0,tregpair=F\n")
            f.write("/\n")
            f.write("&grids tsym=T,nz=25,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\n")
            f.write("/ \n")
            f.write(
                f"&static nprot={Z},nneut={N},e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.15D0,\n"
            )
            f.write(
                f"deocc(1)=3.2,deocc(2)=3.2,betain={betaInp},errenr=1D-12,errvar=0.01,colerr=1D-9,itpini=0,\n"
            )
            f.write("tvaryx0=F,inodiag=9999/\n")
            f.write("&constraint\n")
            f.write("/\n")
            f.write("&thermal\n")
            f.write("/\n")
    elif edf == "UNE1":  # Old UNEDF1 EDF
        with open(f"for005_{index}", "w") as f:
            f.write(
                f"&main outnam='{index}', ipri_fields=0, iprint=2000, itmax=2000,\n"
            )
            f.write("/ \n")
            f.write("&force name='UNEDF1', ipair=5, \n")
            f.write("/ \n")
            f.write("&grids tsym=T,nz=25,nr=25,dz=0.7,dr=0.7 \n")
            f.write("/\n")
            f.write(
                f"&static nprot={Z},nneut={N},e0inv=100D0,x0inv=0.3,tvaryx0=F,addocc=0.3,addrho=0.15,\n"
            )
            f.write(
                "deocc(1)=2.1D0,deocc(2)=2.1D0,errenr=1D-10,errvar=0.01,colerr=1D-9,itpini=0,\n"
            )
            f.write(f"tvaryx0=F,inodiag=9999,betain={betaInp}/\n")
            f.write("/ \n")
            f.write("&constraint \n")
            f.write("/ \n")
            f.write("&thermal\n")
            f.write("/ \n")

    elif edf == "UNE1HFB":  # New UNEDF1 EDF
        with open(f"for005_{index}", "w") as f:
            f.write(
                f"&main outnam='{index}', ipri_fields=0, iprint=2000, itmax=8000,\n"
            )
            f.write(
                f" iterbcs=300,iterhfb=600,pairenhance=1.0D0,iternat=0,modnat=99999,itest=0, \n"
            )
            f.write("/ \n")
            f.write("&force \n")
            f.write(" t0=-2078.328,t1=239.401,t2=1574.243,t3=14263.646,t4=76.736,  \n")
            f.write(
                " x0=0.54304371E-01,x1=-5.0781477,x2=-1.3665705,x3=-0.16116885,b4p=71.31652, \n"
            )
            f.write(" power=0.270018,h2m=20.73553,20.73553,\n")
            f.write(" chg=1D0,ex=1,riftm=0D0,zpe=-2, \n")
            f.write(f" ipair={pairType},v0prot=967.427,v0neut=884.339,rho0pr=0.16, \n")
            f.write(" pair_cutoff=-5D0,HFBwidth=0.0D0,tregpair=F \n")
            f.write("/\n")
            f.write("&grids tsym=T,nz=25,nr=25,dz=0.7,dr=0.7\n")
            f.write("/\n")
            f.write(
                f"&static nprot={Z},nneut={N},e0inv=100D0,x0inv=0.3,tvaryx0=F,addocc=0.3,addrho=0.15,\n"
            )
            f.write(
                "deocc(1)=2.1D0,deocc(2)=2.1D0,errenr=1D-12,errvar=0.01,colerr=1D-9,itpini=0,\n"
            )
            f.write(f"tvaryx0=F,inodiag=9999,betain={betaInp}/\n")
            f.write("/ \n")
            f.write("&constraint \n")
            f.write("/ \n")
            f.write("&thermal\n")
            f.write("/ \n")

    elif edf == "FyIVP":
        with open(f"for005_{index}", "w") as f:
            f.write(
                f"&main outnam='{index}',title='Fy(IVP,s05)',iprint=1000,itmax=8000,iterbcs=300,iterhfb=600,\n"
            )
            f.write("      pairenhance=1.0D0,iternat=0,modnat=9999,itest=0,\n")
            f.write("/ \n")
            f.write("&force \n")
            f.write(" tFayans=.TRUE.,\n")
            f.write(f" name='none',pairstab=-0D0,HFBwidth={hfbwidth},\n")
            f.write("  avp=-9.744306796,h1vp=0.6164610864,h2vp=0.2103424521,\n")
            f.write("  avm=40.51973521,h1vm=-0.5882145708,h2vm=51.73369965,\n")
            f.write("  asp=0.5333352800,hsp=0.0,hgradp=0.3458154900,\n")
            f.write("  asm=0.0,hsm=51.73369965,hgradm=0.0,\n")
            f.write(
                "  fxi=-4.599406100,h1xi=4.605883700,hgradxi=2.329780300,delfxi=-0.1685570700,\n"
            )
            f.write("  pairstab=-0.05,\n")
            f.write(
                "  c0nabj=0.1828568700,c1nabj=6.6576914000E-02,c0tau=0.0,c1tau=0.0,\n"
            )
            f.write(" h2m=20.749811,20.721249,\n")
            f.write("  tcutoff_uvspace=.TRUE.,isovectormodel=2,\n")
            f.write(" chg=1D0,ex=1,zpe=1,tcoulfold=.TRUE.,\n")
            f.write(f" ipair={pairType},pair_cutoff=-15D0,tregpair=F\n")
            f.write("/ \n")
            f.write("&grids tsym=T,nz=25,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\n")
            f.write("/ \n")
            f.write(
                f"&static nprot={Z},nneut={N},e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.15D0,\n"
            )
            f.write(
                "deocc(1)=3,deocc(2)=3,errenr=1D-12,errvar=0.01,colerr=1D-9,itpini=0,\n"
            )
            f.write(f"tvaryx0=F,inodiag=9999,betain={betaInp},\n")
            f.write(f"/ \n")
            f.write("&constraint\n")
            f.write("/\n")
            f.write("&thermal\n")
            f.write("/\n")
    else:
        print("Invalid EDF")


def create_and_submit_job_script(job_name, index):
    job_script_content = f"""#!/bin/bash --login
#SBATCH --job-name={job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --propagate=CORE
#SBATCH --time=1-23:59:59
#SBATCH --mem=64G
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAxOdd/101Sn/slurmjobs/errs/error_{job_name}.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAxOdd/101Sn/slurmjobs/errs/output_{job_name}.out

module purge
module load powertools
module load OpenBLAS/0.3.23-GCC-12.3.0
module load intel/2023a

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
../../SkyAx/wrapper.exe < for005_{index} > /dev/null

"""
    # Save the job script to a file
    job_script_path = os.path.join(job_script_dir, f"{job_name}.sh")
    with open(job_script_path, "w") as job_script_file:
        job_script_file.write(job_script_content)
    # Submit the job script
    subprocess.run(["sbatch", job_script_path])


def main():
    # List of nuclei
    nucleusList = ["Cd", "In", "Sn", "Sb", "Te"]
    # Given EDF
    edf = "UNE1HFB"
    ipair = [8]  # for BCS and HFB
    # HFB width for all EDFs
    hfbwidth = "0.0d0"
    ZList = np.arange(48, 53, 1)
    AList = np.arange(94, 110, 1)
    # Configurations for a given nucleus
    betaInVals = ["-0.2D0", "0.0D0", "0.2D0"]

    # Make directories for functionals
    os.makedirs(edf, exist_ok=True)

    for j in ipair:
        cond = "bcs" if j == 7 else "hfb"
        print(f"cond = {j}, hfbwidth = {hfbwidth}")
        count = 0
        for i, Z in enumerate(ZList):
            NList = AList - Z
            for N in NList:
                nucleus = f"{N+Z}{nucleusList[i]}"
                for betaInp in betaInVals:
                    # Write the for005 file
                    write_for005(edf, Z, N, hfbwidth, j, betaInp, str(count).zfill(3))
                    # Write and execute job scripts
                    job_name = f"{nucleus}_{cond}_{betaInp}"
                    create_and_submit_job_script(job_name, str(count).zfill(3))
                    count += 1

            print(f"{cond} Done!")


if __name__ == "__main__":
    main()
