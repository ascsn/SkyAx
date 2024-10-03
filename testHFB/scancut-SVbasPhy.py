import os
import subprocess
import re

nucleus = "120Sn"
edf = "SVbas"

# Set HFB widths
hfbwidths = ["0.0d0", "1.0d0", "1.5d0", "2.0d0", "3.0d0"]

# BCS pairing strengths
pairStrengthsBCS = ["530", "1620", "5920", "19210", "41820"]

# HFB pairing strengths
pairStrengthsHFB = ["530", "1610", "136", "1500", "130", "1500"]

ipair = [7]  # 7 for BCS 8 for HFB

# Set hfbwidth
hfbwidth = "0.0d0"

def preprocess_line(line):
    # Add a space before negative signs that are not part of the scientific notation
    line = re.sub(r"(\d)(-)", r"\1 -", line)
    return line

def read_nmass_ncharge(file_path):
    with open(file_path, "r") as f:
        for line in f:
            if "# Nmass" in line:
                parts = line.split()
                return f"Nmass: {parts[2]} Ncharge: {parts[4]}"
    return ""


def read_energy(file_path):
    energy_info = ""
    with open(file_path, "r") as f:
        lines = f.readlines()
        for idx, line in enumerate(lines):
            if "epairp" in line:
                next_line = lines[idx + 1]
                parts = preprocess_line(next_line).split()
                energy_info = f"energy: {parts[12]} epairp: {parts[13]} epairn: {parts[14]} rmsc: {parts[18]}"
                break
    return energy_info


def read_gaps(file_path):
    proton_gap = neutron_gap = ""
    with open(file_path, "r") as f:
        for line in f:
            if "proton:" in line:
                proton_gap = f"proton gap: {line.split()[3]}"
            elif "neutron:" in line:
                neutron_gap = f"neutron gap: {line.split()[3]}"
    return f"{proton_gap} {neutron_gap}"

for j in ipair:

    if j == 7:
        cond = "bcs"
        pairStrength = "129.8"
        combined = zip(hfbwidths, pairStrengthsBCS)
    else:
        cond = "hfb"
        pairStrength = "1500"
        combined = list(zip(hfbwidths, pairStrengthsHFB))

    # for hfbwidth, pairStrength in combined:
    print(f"cond = {cond}, hfbwidth={hfbwidth}, pairStrength = {pairStrength}")
    fileName = f"collect_cat_{edf}_{nucleus}_{hfbwidth}_{cond}.dat"
    for i in range(100, 101):
        with open("for005", "w") as f:
            f.write(
                "&main outnam='tst',title='test 120Sn',iprint=1000,itmax=1000\n"
            )
            f.write("/\n")
            f.write("&force\n")
            f.write(f" name='none',pairstab=-0D0,HFBwidth={hfbwidth},\n")
            f.write(
                " t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921,\n"
            )
            f.write(
                " x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530,\n"
            )
            f.write(" t4=124.633,b4p=34.1116989,power=0.3,\n")
            f.write(" h2m=20.749811,20.721249,\n")
            f.write(" chg=1D0,ex=1,riftm=0D0,zpe=1,\n")
            f.write(f" ipair=7,v0prot=10,v0neut={pairStrength},rho0pr=0.2011591,\n")
            f.write(f" pair_cutoff=-{i},tregpair=F\n")
            f.write("/\n")
            # f.write("&grids tsym=F,nz=81,nr=41,dr=0.5D0,dz=0.5D0,idrtyp=0,\n")
            f.write("&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\n")
            f.write("/\n")
            f.write(
                "&static nprot=50,nneut=70,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0,\n"
            )
            f.write(
                f"deocc(1)=13.0,deocc(2)=13.0,errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0,\n"
            )
            f.write("tvaryx0=T,inodiag=9999/\n")
            f.write("&constraint\n")        
            f.write("/\n")
            f.write("&thermal /\n")

        os.system("ulimit -s unlimited")
        os.system("export OMP_STACKSIZE=4500m")
        subprocess.run(
            ["../skyax.gfortran.omp"], 
            stdin=open("for005"),
            stdout=subprocess.DEVNULL,
        )

        with open(fileName, "a") as f:
            f.write(
                f"cutoff=-{i},hfbwidth={hfbwidth},pair.strength={pairStrength}\n"
            )

        # Extract required information from output files using readlines
        nmass_ncharge = read_nmass_ncharge("ensurf.tst")
        energy_info = read_energy("ensurf.tst")
        gap_info = read_gaps("out_details.tst")

        # Combine and write output to file
        with open(fileName, "a") as f:
            f.write(f"{nmass_ncharge} {energy_info} {gap_info}\n")

    print(f"{cond} Done!")
os.system(f"mv collect_cat_*.dat collect_cat/")
