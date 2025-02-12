import os
import re
import pandas as pd

def process_data(file_path):
    # Define column headers based on your data
    columns = [
        "bet2", "bet2p", "bet2n", "bet4", "bet4p", "bet4n", "energy", "epairp",
        "epairn", "iter", "rms", "rmsp", "rmsn", "jz_p", "jz_n", "j"
    ]
    rows = []  # List to store the extracted rows as dictionaries
    nmass = None
    ncharge = None
    min_energy_row = None
    block_rows = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

        for i, line in enumerate(lines):
            # Check for Nmass and Ncharge line and extract their values
            if re.search(r"# Nmass=", line):
                # Process the previous block
                if nmass is not None:
                    if nmass % 2 == 0 and min_energy_row:
                        row_data = dict(zip(columns, min_energy_row.split()))
                        row_data["Nmass"] = nmass
                        row_data["Ncharge"] = ncharge
                        rows.append(row_data)
                    elif nmass % 2 != 0 and block_rows:
                        min_row = min(block_rows, key=lambda x: float(x.split()[6]))
                        row_data = dict(zip(columns, min_row.split()))
                        row_data["Nmass"] = nmass
                        row_data["Ncharge"] = ncharge
                        rows.append(row_data)

                # Reset variables for the new block
                nmass = int(re.search(r"Nmass=\s+(\d+)", line).group(1))
                ncharge = int(re.search(r"Ncharge=\s+(\d+)", line).group(1))
                min_energy_row = None
                block_rows = []

            # Check for the exact #  bet2 line and get the energy row
            elif re.match(r"#\s+bet2\s+bet2p\s+bet2n", line):
                if i + 1 < len(lines):  # Ensure the next line exists
                    energy_line = lines[i + 1].strip()
                    if energy_line:
                        min_energy_row = energy_line  # Store the entire row

            # Check for # blocking lines and collect rows
            elif re.match(r"# blocking", line):
                if i + 1 < len(lines):  # Ensure the next line exists
                    block_line = lines[i + 1].strip()
                    if block_line:
                        block_rows.append(block_line)  # Store the entire row

        # Handle the last block
        if nmass is not None:
            if nmass % 2 == 0 and min_energy_row:
                row_data = dict(zip(columns, min_energy_row.split()))
                row_data["Nmass"] = nmass
                row_data["Ncharge"] = ncharge
                rows.append(row_data)
            elif nmass % 2 != 0 and block_rows:
                min_row = min(block_rows, key=lambda x: float(x.split()[6]))
                row_data = dict(zip(columns, min_row.split()))
                row_data["Nmass"] = nmass
                row_data["Ncharge"] = ncharge
                rows.append(row_data)

    # Convert to a DataFrame
    df = pd.DataFrame(rows)

    # Ensure all numeric values are in scientific notation
    for col in df.columns:
        if col not in {"iter", "Nmass", "Ncharge"}:  # Skip integer columns
            df[col] = pd.to_numeric(df[col], errors="coerce").map("{:.6E}".format)

    return df

def find_minimum_energy_in_folder(folder_path):
    # Collect DataFrames for all `ensurf` files
    all_dfs = []
    for filename in os.listdir(folder_path):
        if filename.startswith("ensurf.") and filename.split('.')[-1].isdigit():
            file_path = os.path.join(folder_path, filename)
            print(f"Processing file: {file_path}")
            df = process_data(file_path)
            all_dfs.append(df)

    # Combine all DataFrames into one
    if all_dfs:
        combined_df = pd.concat(all_dfs, ignore_index=True)

        # Find the row with the minimum energy
        combined_df["energy"] = pd.to_numeric(combined_df["energy"], errors="coerce")  # Convert energy to numeric
        min_energy_row = combined_df.loc[combined_df["energy"].idxmin()]  # Row with minimum energy
        return min_energy_row
    else:
        return None

# Path to the parent folder containing subfolders
parent_folder = "FyIVP"

# Iterate through subfolders and find minimum energy in each
results = []
for subfolder in os.listdir(parent_folder):
    subfolder_path = os.path.join(parent_folder, subfolder)
    if os.path.isdir(subfolder_path):
        print(f"Processing folder: {subfolder_path}")
        min_energy_row = find_minimum_energy_in_folder(subfolder_path)
        if min_energy_row is not None:
            results.append(min_energy_row)

# Combine results and save to a file
if results:
    results_df = pd.DataFrame(results)
    results_df = results_df.sort_values(by=["Ncharge", "Nmass"], ascending=[True, True],ignore_index=True)
    results_df.to_csv(f"minimum_energies_{parent_folder}.dat", sep="\t", index=False, header=True)
    print("Results saved to 'minimum_energies_across_folders.dat'.")
