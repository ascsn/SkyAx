import pandas as pd
import re
import glob
import os
import numpy as np

BASE_PATH=os.getcwd()

# pairStrengthsBCS = ["674","2335", "2031", "1976", "1958"]
# pairStrengthsHFB = ["674.6","1862", "1674", "1563", "1533"] # HFBWIDTH = 1.0
# deoccList = ["6.D0", "8.D0", "10.D0","12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]   
deoccList = ["2.D0", "4.D0", "6.D0", "8.D0", "10.D0", "12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]
runList = [[5, 6], [10,11], [20,21], [40,41], [60,61]]

# Define file path pattern to match all job files (try both patterns)
file_pattern1 = 'logFayans/job_170Yb_*_output.log'
file_pattern2 = 'logFayans/job_170Yb_*.log'

# Verify files matched by glob
files1 = glob.glob(file_pattern1)
files2 = glob.glob(file_pattern2)
# Combine and remove duplicates
files = list(set(files1 + files2))
print(f"Found {len(files)} log files to process")
if len(files) == 0:
    print("ERROR: No log files found matching patterns:")
    print(f"  - {file_pattern1}")
    print(f"  - {file_pattern2}")
    print(f"Current directory: {os.getcwd()}")
    import sys
    sys.exit(1)

columns = [
    "bet2", "bet2p", "bet2n", "bet4", "bet4p", "bet4n", "energy", "epairp", "epairn", "ferp", "fern", "gapp", "gapn",  "iter", 
    "rmsc", "rms", "rmsp", "rmsn", "N_boundp", "N_boundn", "last_value"
]

# Check if ensurf.job exists
if not os.path.exists('ensurf.job'):
    print("ERROR: ensurf.job file not found in current directory")
    print(f"Current directory: {os.getcwd()}")
    import sys
    sys.exit(1)

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

# Read and clean the ensurf.job file
# Skip header lines and filter out malformed lines
with open('ensurf.job', 'r') as f:
    lines = f.readlines()

# Keep only valid data lines (skip header, comments, and malformed lines)
clean_lines = []
for i, line in enumerate(lines):
    # Skip header lines (first 9 lines)
    if i < 9:
        continue
    # Skip empty lines or lines with "part.nr."
    if not line.strip() or 'part.nr.' in line:
        continue
    # Split and check if we have the expected number of fields (21)
    fields = line.split()
    if len(fields) == len(columns):
        clean_lines.append(line)
    elif len(fields) > len(columns):
        # Handle concatenated lines - split into separate entries
        print(f"Warning: Line {i+1} has {len(fields)} fields (expected {len(columns)}), attempting to split...")
        # Try to split into chunks of 21 fields
        for j in range(0, len(fields), len(columns)):
            if j + len(columns) <= len(fields):
                chunk = ' '.join(fields[j:j+len(columns)]) + '\n'
                clean_lines.append(chunk)

# Write cleaned data to a temporary object and read with pandas
from io import StringIO
clean_data = ''.join(clean_lines)
data = pd.read_csv(StringIO(clean_data), sep=r'\s+', names=columns, engine='python')
print(f"Loaded {len(data)} entries from ensurf.job")
print(f"Sample energy values: {data['energy'].head().tolist()}")

# Function to extract energy and rms_C from the last three lines of a file
def extract_energy_and_rmsc(file_path):
    """
    Extract energy and rms_C from the last three lines of the log file.
    Expected format:
     iter= 6000, energies: from rearr. vs. from densities=        0.0000    -1378.0953
    rms_C^2:   28.741       28.226      0.72970     -0.16771      3.29670E-02  5.51782E-02 -6.35803E-02 -7.14287E-02   28.740    
      rms_C:   5.3610       5.3610    
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()
        
    # Get last several lines to search more thoroughly
    if len(lines) < 3:
        return None, None
    
    # Search last 10 lines for more robustness
    search_lines = lines[-10:] if len(lines) >= 10 else lines
    energy = None
    rmsc = None
    
    for line in search_lines:
        # Extract energy from the line with "energies: from rearr. vs. from densities="
        if 'energies: from rearr. vs. from densities=' in line:
            # Match the second number after "energies: from rearr. vs. from densities="
            match = re.search(r'energies:\s*from rearr\.\s*vs\.\s*from densities=\s*[\d.eE+-]+\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', line)
            if match:
                energy = float(match.group(1))
        
        # Extract rms_C from the line with "rms_C:" (but not "rms_C^2:")
        if 'rms_C:' in line and 'rms_C^2:' not in line:
            # Match the first number after "rms_C:"
            match = re.search(r'rms_C:\s+([-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)', line)
            if match:
                rmsc = float(match.group(1))
    
    return energy, rmsc
        

# List to collect matched data
matched_data = []
file_dict = dict()
skipped_files = []
print(f"\nEnergy values in ensurf.job: {len(data)} entries")
print(f"Energy range: {data['energy'].min():.2f} to {data['energy'].max():.2f} MeV\n")

# Process each file
for i, file_path in enumerate(files):
    # Extract energy and rms_C from the file
    energy, rmsc = extract_energy_and_rmsc(file_path)
    
    if energy is None:
        print(f"Warning: Could not extract energy from {file_path}")
        skipped_files.append((file_path, "No energy found"))
        continue
    
    print(f"File: {file_path}")
    rmsc_str = f"{rmsc:.4f}" if rmsc is not None else "N/A"
    print(f"  Extracted - Energy: {energy:.4f}, rms_C: {rmsc_str}")
    
    # Find matches using tolerance-based approach
    # First match on energy with tolerance of 0.1 MeV (increased for better matching)
    energy_tolerance = 0.1
    energy_matches = data[np.abs(data['energy'] - energy) < energy_tolerance]
    
    if energy_matches.empty:
        print(f"  No match found in ensurf.job (energy tolerance: {energy_tolerance})")
        print(f"     Available energies range: {data['energy'].min():.2f} to {data['energy'].max():.2f}")
        skipped_files.append((file_path, f"No energy match (E={energy:.2f})"))
        continue
    
    # If rms_C is available, use it as additional matching criterion with tolerance
    match_rows = energy_matches
    if rmsc is not None and len(energy_matches) > 1:
        rmsc_tolerance = 0.1  # Increased tolerance
        rmsc_matches = energy_matches[np.abs(energy_matches['rmsc'] - rmsc) < rmsc_tolerance]
        if not rmsc_matches.empty:
            match_rows = rmsc_matches
            print(f"  Matched with ensurf.job entry (energy and rms_C match)")
        else:
            # Use best energy match if no rms_C match found
            match_rows = energy_matches.iloc[[0]]
            print(f"  Matched with ensurf.job entry (energy match only)")
    else:
        # Single match or no rmsc to compare
        match_rows = energy_matches.iloc[[0]]
        print(f"  Matched with ensurf.job entry")
    
    file_dict[file_path] = match_rows
    matched_data.append((file_path, match_rows))

# Pattern to match both formats:
# - job_170Yb_hfb_1.0d0_10.D0_10_11_output.log
# - job_170Yb_hfb_1.0d0_10.D0_10_11.log (without _output)
pattern = r"job_(\d+[A-Za-z]+)_(bcs|hfb)_(\d+\.\d+d0)_([\d.]+)\.D0_(\d+)_(\d+)(?:_output)?\.log"


finalDF = pd.DataFrame()

# Process each file and its corresponding DataFrame row
if not file_dict:
    print("No files were successfully matched with ensurf.job entries.")
    print("Please check that:")
    print("  1. Log files contain the expected last 3 lines with energy and rms_C")
    print("  2. Energy values are within tolerance range of ensurf.job entries")
else:
    try:
        for filename, df in file_dict.items():
            # Extract values from the filename
            file_match = re.search(pattern, filename)
            # print("file_match", file_match)
            
            if not file_match:
                print(f"Warning: Filename {filename} does not match expected pattern, skipping...")
                skipped_files.append((filename, "Pattern mismatch"))
                continue
                
            dfdict = dict()
            # print("groups", file_match.groups())
            # Updated to match 6 groups: nucleus, pairType, hfbwidth, deocc, cutoff, cutoffHigh
            Nucleus, pairType, hfbwidth, deocc, cutoff, cutoffHigh = file_match.groups()
            print("deocc, cutoff:", deocc, cutoff)
            dfdict["Cutoff"] = cutoff
            dfdict["Deocc"] = str(deocc).zfill(2)
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
            # print("Energy:",energy)
                    
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

        if not finalDF.empty:
            # Convert Cutoff to int for proper numerical sorting
            finalDF['Cutoff'] = finalDF['Cutoff'].astype(int)
            finalDF = finalDF.sort_values(by=['Cutoff', 'Pair Type', 'Deocc'], ascending=[True, True, True])

            with open("summary_170Yb_Fayans.txt", 'w') as f:
                f.write(f"Nucleus: {test_label}\n")
                f.write(f"Nmass: {nmass} \t Ncharge: {ncharge}\n")
                f.write(finalDF.to_string(header=True, index=False))
            
            print(f"\nSummary written to summary_170Yb_Fayans.txt")
            print(f"Total entries: {len(finalDF)}")
            print(f"Successfully processed: {len(file_dict)} files")
            print(f"Skipped: {len(skipped_files)} files")
            
            if skipped_files:
                print("\n--- Skipped Files ---")
                for fname, reason in skipped_files:
                    print(f"  {os.path.basename(fname)}: {reason}")
        else:
            print("No valid entries to write to summary file.")

    except KeyError as e:
        print(f"Error processing data: {e}")
        print("This might indicate missing columns in the ensurf.job file.")

# print(finalDF)
        
