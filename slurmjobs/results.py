import pandas as pd
import re
import glob
import os
import numpy as np

BASE_PATH=os.getcwd()

pairStrengthsBCS = ["674.6", "2335", "2031", "1976", "1958"]
pairStrengthsHFB = ["674.6","1862", "1674", "1563", "1533"] # HFBWIDTH = 1.0
# deoccList = ["6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]   
runList = [[5,6], [10,11], [25,26], [40,41], [60,61]]

# Define file path pattern to match all job files
file_pattern = 'job_170Yb_*out.log'

# Verify files matched by glob
files = glob.glob(file_pattern)

columns = [
    "bet1", "bet1p", "bet1n", "bet2", "bet2p", "bet2n", "bet3", "bet3p", "bet3n", 
    "bet4", "bet4p", "bet4n", "energy", "epairp", "epairn", "ferp", "fern", "gapp", "gapn",  "iter", 
    "rmsc", "rms", "rmsp", "rmsn", "N_boundp", "N_boundn", "last_value"
]
with open('ensurf.job', 'r') as f:
    readLines = f.readlines()
    for line in readLines:
        # Extract test label
        test_match = re.search(r"#test\s+(\w+)", line)
        if test_match:
            test_label = test_match.group(1)
    
        if "# Nmass" in line:
            match = re.search(r"# Nmass=\s*(\d+), Ncharge=\s*(\d+)", line)
            if match:
                nmass = int(match.group(1))
                ncharge = int(match.group(2))
        

data = pd.read_csv('ensurf.job', sep='\\s+', comment='#', names=columns, skiprows=8)

def extract_last_energy(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Ensure the file has at least one line
    if len(lines) > 0:
        # If there are two or more lines, extract the energies from both
        energies = []
        
        # Check if there is a second-to-last line and extract energy
        if len(lines) >= 3:
            target_line_1 = lines[-3].strip()
            match_1 = re.search(r"(-?\d+\.\d+)$", target_line_1)
            if match_1:
                energies.append(float(match_1.group(1)))
        
        # Extract energy from the last line
        target_line_2 = lines[-2].strip()
        match_2 = re.search(r"(-?\d+\.\d+)$", target_line_2)
        if match_2:
            energies.append(float(match_2.group(1)))
        
        # Return the list of energies (1 or 2 energies depending on the file)
        if energies:
            return energies
        else:
            print("No matching energy value found.")
            return None
    else:
        print("File is empty.")
        return None
    

def extract_second_energy(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    energies = []  # List to store extracted energies
    for line in lines:
        # Match the line containing the energies
        match = re.search(r'iter=\s*\d+,\s*energies:\s*from rearr. vs. from densities=\s*([0-9.-]+)\s+([-+]?[0-9]*\.?[0-9]+)', line)
        if match:
            # Append the second energy found to the list
            energies.append(float(match.group(2)))  # Extract the second energy

    # Return the second energy if it exists
    if len(energies) > 1:
        return energies[1]  # Return the second occurrence
    else:
        return None  # Or handle the case where it does not exist


# Function to extract energy from the last three lines of a file
def extract_energy(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
        # print(re.search(r"(-?\d+\.\d+)$", lines[-3].strip()).group(1))
        # Extract the energy value from the last relevant line
        match = re.search(r"energies:\s*from rearr\.\s*vs\.\s*from densities=\s*[\d.]+\s*(-?\d+\.\d+)", line)
        # match = re.search(r"(-?\d+\.\d+)$", lines[-3].strip())  # Matches the last floating number
        return float(match.group(1)) if match else None
        

# List to collect matched data
matched_data = []
file_dict = dict()
data['energy_rounded'] = np.round(data['energy'], 4)

# Process each file
for file_path in files:
    # Extract energy from the file
    # energy = extract_energy(file_path)
    energy = extract_second_energy(file_path)
    # energy = extract_last_energy(file_path)
    # print("Energies", energy, file_path)
    
    # Find matches in the main data and create dictionary
    match_rows = data[data["energy_rounded"] == energy]
    # match_rows = data[data['energy'].round(4) == energy]
    
    if not match_rows.empty:
        file_dict[file_path] = match_rows
        matched_data.append((file_path, match_rows))

# pattern = r"job_hfb_0\.0d0_(\d+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
# pattern = r"job__1\.0d0_([\d.]+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
# pattern = r"job_hfb_1\.0d0_([\d.]+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
# pattern = r"job_(\d+[A-Za-z]+)_(bcs|hfb)_(\d+\.\d+d0)_(\d+)_([\d.]+)\.D0_(\d+)_(\d+)_out\.log"
pattern = r"job_(\d+[A-Za-z]+)_(bcs|hfb)_(\d+\.\d+d0)_(\d+(?:\.\d+)?)_([\d]+)\.D0_(\d+)_(\d+)_out\.log"



finalDF = pd.DataFrame()

# Process each file and its corresponding DataFrame row
try:
    for filename, df in file_dict.items():
        # Extract values from the filename
        file_match = re.search(pattern, filename)
        print("file_match", file_match)
        dfdict = dict()
        if file_match:
            print("groups", file_match.groups())
            Nucleus, pairType, hfbwidth, pair_strength, deocc, cutoff, cutoffHigh = file_match.groups()
            print("pair", pair_strength, deocc, cutoff)
            dfdict["Cutoff"] = cutoff
            dfdict["Deocc"] = str(deocc).zfill(2)
            dfdict["Pair Strength"] = pair_strength
            dfdict["Pair Type"] = pairType
            dfdict["hfb width"] = hfbwidth
            
            # Extract specific values from the DataFrame
            energy = df['energy'].iloc[0]
            rmsc = df['rmsc'].iloc[0]
            rmsn = df['rmsn'].iloc[0]
            N_boundn = df['N_boundn'].iloc[0]
            epairp = df["epairp"].iloc[0]
            epairn = df["epairn"].iloc[0]
            gapp = df["gapp"].iloc[0]
            gapn = df["gapn"].iloc[0]
            beta2 = df["bet2"].iloc[0]
            beta2p = df["bet2p"].iloc[0]
            beta2n = df["bet2n"].iloc[0]

            # print stuff
            print("Energy:",energy)
                    
            dfdict["Energy"] = energy
            dfdict["epairp"] = epairp
            dfdict["epairn"] = epairn
            dfdict["gap_p"] = gapp
            dfdict["gap_n"] = gapn
            dfdict["beta2"] = beta2
            dfdict["beta2p"] = beta2p
            dfdict["beta2n"] = beta2n
            dfdict["RMSC"] = rmsc
            dfdict["RMSN"] = rmsn
            dfdict["N_boundn"] = N_boundn

            dfupdate = pd.DataFrame([dfdict])
            finalDF = pd.concat([finalDF, dfupdate], ignore_index=True)

    finalDF = finalDF.sort_values(by=['Cutoff', 'Pair Type', 'Deocc'], ascending=[True, True, True])

    with open("summary_170Yb.txt", 'w') as f:
        f.write(f"Nucleus: {test_label}\n")
        f.write(f"Nmass: {nmass} \t Ncharge: {ncharge}\n")
        f.write(finalDF.to_string(header=True, index=False))

except KeyError:
    print("Data file empty")

# print(finalDF)
        
