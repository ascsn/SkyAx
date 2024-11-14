import os
import subprocess
import pandas as pd
import re
from io import StringIO  # For testing with string data

nucleus = "120Sn"
edf = "SVbas"

# Set HFB widths
hfbwidths = ["1.0d0", "1.0d0", "1.0d0", "1.0d0", "1.0d0", "1.0d0"]

# BCS pairing strengths
# pairStrengthsBCS = ["142.5", "1860", "5920", "19210", "41820"]
# pairStrengthsBCS = ["530", "1620", "136", "1500",  "139", "1500"]
# pairStrengthsBCS = ["432"]#, "630.5", "630.5", "630.5", "630.5", "630.5 "]
# pairStrengthsBCS = ["630.5", "2297","5984", "12360", "36100"] #Cutoff 100
pairStrengthsBCS = ["537", "478", "451", "438", "433", "432"]

# HFB pairing strengths
# pairStrengthsHFB = ["530", "1610", "136", "1500", "130", "1500"]#, "10815", "30301"]
pairStrengthsHFB = ["728", "562", "508", "474", "428", "489"]  # HFBWIDTH = 1.0
# pairStrengthsHFB =["537","478","451","438","433","432"]

ipair = [8]  # 7 for BCS 8 for HFB

# Set hfbwidth
# hfbwidth = "0.0d0"


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
                energy_info = f"energy: {parts[12]} epairp: {parts[13]} epairn: {parts[14]} rmsc: {parts[18]} beta2: {parts[3]}"
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


def find_max_energy(file_path):

    # Read the entire file into memory
    with open(file_path, "r") as file:
        content = file.read()

    # Use regex to find all relevant tables
    # The regex assumes that each table starts with 'Nr' and ends with a blank line
    tables = re.findall(r"Nr.*?(?=\n\s*\n)", content, re.DOTALL)

    # Initialize a list to hold DataFrames for each table
    dataframes = []

    # Process each table
    for table in tables:
        # Split the table into lines and clean up
        lines = table.strip().split("\n")

        # Load the relevant data into a pandas DataFrame
        try:
            df = pd.read_csv(StringIO("\n".join(lines)), delim_whitespace=True)
            dataframes.append(df)
        except pd.errors.EmptyDataError:
            print("Empty DataFrame encountered for a table. Skipping...")

    # Display extracted DataFrames and find maximum energy for J_z = 0.5 in each
    for i, df in enumerate(dataframes):
        # print(f"Table {i + 1}:")
        is_numeric = pd.to_numeric(df.iloc[-1], errors="coerce").notnull().all()
        if is_numeric == False:
            df.drop(df.tail(1).index, inplace=True)  # drop last n rows
        # Filter for J_z = 0.5 and find the maximum energy
        # Somehow the column heads are messed up, so needs to be fixed (later)
        if "energy" in df.columns:
            df["Iso"] = df["Iso"].astype(float)
            filtered_df = df[df["Iso"] == 0.5]
            if not filtered_df.empty:
                max_energy = filtered_df["delta"].max()
                # print(f"The maximum energy corresponding to J_z = 0.5 in Table {i + 1} is: {max_energy}\n")
            else:
                print(f"No data for J_z = 0.5 found in Table {i + 1}.\n")
        if i == 4:
            return f"Max_SP_energy: {max_energy}"


# deoccList = ["6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]
# runList = [[10, 11],[20,21], [40, 41], [60, 61], [80, 81], [100, 101]]
deoccList = ["16.D0", "16.D0", "16.D0", "16.D0", "16.D0", "16.D0"]
runList = [[10, 11], [20, 21], [40, 41], [60, 61], [80, 81], [100, 101]]

for j in ipair:
    for val in deoccList:
        if j == 7:
            cond = "bcs"
            pairStrength = "1610"
            combined = zip(hfbwidths, pairStrengthsBCS, runList)
        else:
            cond = "hfb"
            pairStrength = "1500"
            combined = list(zip(hfbwidths, pairStrengthsHFB, runList))

        print("Val =", val, "j=", j)
        for hfbwidth, pairStrength, runVal in combined:
            print(
                f"cond = {cond}, hfbwidth={hfbwidth}, pairStrength = {pairStrength}, deocc={val}"
            )
            fileName = f"test_collect_cat_{edf}_{nucleus}_{hfbwidth}_{cond}.dat"
            a, b = runVal
            for i in range(a, b):
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
                    f.write(
                        f" ipair={j},v0prot=10,v0neut={pairStrength},rho0pr=0.2011591,\n"
                    )
                    f.write(f" pair_cutoff=-{i},tregpair=F\n")
                    f.write("/\n")
                    # f.write("&grids tsym=F,nz=81,nr=41,dr=0.5D0,dz=0.5D0,idrtyp=0,\n")
                    f.write("&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\n")
                    f.write("/\n")
                    f.write(
                        "&static nprot=50,nneut=70,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0,\n"
                    )
                    f.write(
                        f"deocc(1)={val},deocc(2)={val},errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0,\n"
                    )
                    f.write("tvaryx0=T,inodiag=9999/\n")
                    f.write("&constraint\n")
                    f.write("/\n")
                    f.write("&thermal /\n")

                os.system("ulimit -s unlimited")
                os.system("export OMP_STACKSIZE=4500m")
                subprocess.run(
                    ["../skyax.ifort.omp"],
                    stdin=open("for005"),
                    stdout=subprocess.DEVNULL,
                )

                with open(fileName, "a") as f:
                    f.write(f"deocc(1,2)= {val}\n")
                    f.write(
                        f"cutoff=-{i},hfbwidth={hfbwidth},pair.strength={pairStrength}\n"
                    )

                # Extract required information from output files using readlines
                nmass_ncharge = read_nmass_ncharge("ensurf.tst")
                energy_info = read_energy("ensurf.tst")
                gap_info = read_gaps("out_details.tst")
                max_energy = find_max_energy("out_details.tst")

                # Combine and write output to file
                with open(fileName, "a") as f:
                    f.write(f"{nmass_ncharge} {energy_info} {gap_info} {max_energy}\n")

    print(f"{cond} Done!")
os.system(f"mv test_collect_cat_*.dat collect_cat_120Sn/")

# def main():
#     data = read_energy("ensurf.tst")
#     print("data", data)

# if __name__=="__main__":
#     main()
