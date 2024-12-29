import pandas as pd
import re
import glob
import os
import numpy as np

BASE_PATH=os.getcwd()

# Define file path pattern to match all job files
file_pattern = 'job_170Yb_*out.log'

# Verify files matched by glob
files = glob.glob(file_pattern)

columns = ["bet2", "bet2p", "bet2n", "bet4", "bet4p", "bet4n", "energy", "epairp", "epairn", "ferp", "fern", "gap_p","gap_n", "iter", "rmsc", "rms", "rmsp", "rmsn", "N_boundp", "N_boundn", "last_value"]
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

def preprocess_file(input_file, output_file):
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for line in infile:
            parts = line.split()
            if len(parts) >= 20:  # Ensure enough columns in each row
                # Merge the first three columns into one
                first_column = " ".join(parts[:3])
                rest = parts[3:]
                new_line = f"{first_column} {' '.join(rest)}\n"
                outfile.write(new_line)

# Preprocess the file
# preprocess_file("ensurf.job", "preprocessed_file.txt")

# # Reload the preprocessed file
# data = pd.read_csv(
#     "preprocessed_file.txt",
#     delim_whitespace=True,
#     header=None,
# )

# data.columns = columns

data = pd.read_csv('ensurf.job', sep='\\s+', comment='#', names=columns, skiprows=8)

def extract_energy(file_path):
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
    elif len(energies) == 1:
        return energies[0]
    else:
        return None  # Or handle the case where it does not exist

# List to collect matched data
matched_data = []
file_dict = dict()
data['energy_rounded'] = np.round(data['energy'].astype(float), 4)

print("data energy:", data['energy'], data['energy_rounded'])

# Process each file
for file_path in files:
    # Extract energy from the file
    energy = extract_energy(file_path)
    print("Energies", energy, file_path)
    
    # Find matches in the main data and create dictionary
    # match_rows = data[data["energy_rounded"] == energy]
    match_rows = data[data['energy'].round(4) == energy]
    
    if not match_rows.empty:
        file_dict[file_path] = match_rows
        matched_data.append((file_path, match_rows))

# pattern = r"job_hfb_0\.0d0_(\d+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
# pattern = r"job__1\.0d0_([\d.]+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
# pattern = r"job_hfb_1\.0d0_([\d.]+)_([\d.]+)\.D0_(\d+)_\d+_output\.log"
# pattern = r"job_(\d+[A-Za-z]+)_(bcs|hfb)_(\d+\.\d+d0)_(\d+)_([\d.]+)\.D0_(\d+)_(\d+)_out\.log"
pattern = r"job_(\d+[A-Za-z]+)_(bcs|hfb)_(\d+\.\d+d0)_(\d+(?:\.\d+)?)_([\d]+)\.D0_(\d+)_(\d+)_out\.log"

# print("Matched Data", matched_data)
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
            # print("pair", Nucleus, pairType, hfbwidth, pair_strength, deocc, cutoff)
            dfdict["Cutoff"] = cutoff.zfill(2)
            dfdict["Deocc"] = str(deocc).zfill(2)
            dfdict["Pair_Strength"] = pair_strength
            dfdict["Pair_Type"] = pairType
            dfdict["hfb_width"] = hfbwidth
            
            # Extract specific values from the DataFrame
            energy = df['energy'].iloc[0]
            rmsc = df['rmsc'].iloc[0]
            rmsn = df['rmsn'].iloc[0]
            N_boundn = df['N_boundn'].iloc[0]
            epairp = df["epairp"].iloc[0]
            epairn = df["epairn"].iloc[0]
            gapp = df["gap_p"].iloc[0]
            gapn = df["gap_n"].iloc[0]
            beta2 = df["bet2"].iloc[0]
            beta2p = df["bet2p"].iloc[0]
            beta2n = df["bet2n"].iloc[0]

            # print stuff
            print("Energy:", energy, rmsc, rmsn, N_boundn)
                    
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

    finalDF = finalDF.sort_values(by=['Cutoff', 'Pair_Type', 'Deocc'], ascending=[True, True, True])

    with open("summary_170Yb_Fayans_hfb_0.txt", 'w') as f:
        f.write(f"Nucleus: {test_label}\n")
        f.write(f"Nmass: {nmass} \t Ncharge: {ncharge}\n")
        f.write(finalDF.to_string(header=True, index=False))

except KeyError:
    print("Data file empty")

# print(finalDF)
        
