import pandas as pd
import re
import glob
import os
import numpy as np

BASE_PATH=os.getcwd()

pairStrengthsHFB = ["628", "462", "408", "374", "328.3", "389.2"] #Final values [20] - 462
deoccList = ["6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]   
runList = [[10,11], [20,21], [40,41], [60,61],[80,81],[100,101]]

# Define file path pattern to match all job files
file_pattern = 'job_hfb_*.log'

# Verify files matched by glob
files = glob.glob(file_pattern)

columns = [
    "bet1", "bet1p", "bet1n", "bet2", "bet2p", "bet2n", "bet3", "bet3p", "bet3n", 
    "bet4", "bet4p", "bet4n", "energy", "epairp", "epairn", "ferp", "fern", "iter", 
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
        

data = pd.read_csv('ensurf.job', delim_whitespace=True, comment='#', names=columns, skiprows=8)

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
        # Extract the energy value from the last relevant line
        # match = re.search(r"energies:\s*from rearr\.\s*vs\.\s*from densities=\s*[\d.]+\s*(-?\d+\.\d+)", line)
        match = re.search(r"(-?\d+\.\d+)$", lines[-3].strip())  # Matches the last floating number
        return float(match.group(1)) if match else None
        

# List to collect matched data
matched_data = []
file_dict = dict()
data['energy_rounded'] = np.round(data['energy'], 4)

# Process each file
for file_path in files:
    # Extract energy from the file
    energy = extract_energy(file_path)
    energy1 = extract_second_energy(file_path)
    
    # Find matches in the main data and create dictionary
    match_rows = data[data["energy_rounded"] == energy]
    # match_rows = data[data['energy'].round(4) == energy]
    
    if not match_rows.empty:
        file_dict[file_path] = match_rows
        matched_data.append((file_path, match_rows))

# pattern = r"job_hfb_0\.0d0_(\d+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
pattern = r"job_hfb_0\.0d0_([\d.]+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"

finalDF = pd.DataFrame()

# Process each file and its corresponding DataFrame row
for filename, df in file_dict.items():
    # Extract values from the filename
    file_match = re.search(pattern, filename)
    dfdict = dict()
    if file_match:
        pair_strength, deocc, cutoff = file_match.groups()
        dfdict["Cutoff"] = cutoff
        dfdict["Deocc"] = str(deocc).zfill(2)
        dfdict["Pair Strength"] = pair_strength
        
        # Extract specific values from the DataFrame
        energy = df['energy'].iloc[0]
        rmsc = df['rmsc'].iloc[0]
        rmsn = df['rmsn'].iloc[0]
        N_boundn = df['N_boundn'].iloc[0]
        epairn = df["epairn"].iloc[0]
                
        dfdict["Energy"] = energy
        dfdict["epairn"] = epairn
        dfdict["RMSC"] = rmsc
        dfdict["RMSN"] = rmsn
        dfdict["N_boundn"] = N_boundn

        dfupdate = pd.DataFrame([dfdict])
        finalDF = pd.concat([finalDF, dfupdate], ignore_index=True)

finalDF = finalDF.sort_values(['Cutoff', 'Deocc'], ascending=[True, True])

with open("summary.txt", 'w') as f:
    f.write(f"Nucleus: {test_label}\n")
    f.write(f"Nmass: {nmass} \t Ncharge: {ncharge}\n")
    f.write(finalDF.to_string(header=True, index=False))

# print(finalDF)
        
