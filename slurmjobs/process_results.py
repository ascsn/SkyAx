import argparse
import os
import re
import pandas as pd
from io import StringIO

def preprocess_line(line):
    line = re.sub(r"(\d)(-)", r"\1 -", line)
    return line

def read_nmass_ncharge(file_path):
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            for line in f:
                if "# Nmass" in line:
                    parts = line.split()
                    return f"Nmass: {parts[2]} Ncharge: {parts[4]}"
    return "Nmass and Ncharge not found"

def read_energy(file_path):
    energy_info = ""
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                if "epairp" in line:
                    next_line = lines[idx + 1]
                    parts = preprocess_line(next_line).split()
                    energy_info = f"energy: {parts[12]} epairp: {parts[13]} epairn: {parts[14]} rmsc: {parts[18]} beta2: {parts[3]}"
                    break
    return energy_info or "Energy information not found"

def read_gaps(file_path):
    proton_gap = neutron_gap = ""
    if os.path.exists(file_path):
        with open(file_path, "r") as f:
            for line in f:
                if "proton:" in line:
                    proton_gap = f"proton gap: {line.split()[3]}"
                elif "neutron:" in line:
                    neutron_gap = f"neutron gap: {line.split()[3]}"
    return f"{proton_gap} {neutron_gap}" or "Gap information not found"

def find_max_energy(file_path):
    max_energy = None
    if os.path.exists(file_path):
        with open(file_path, 'r') as file:
            content = file.read()
            tables = re.findall(r'Nr.*?(?=\n\s*\n)', content, re.DOTALL)
            dataframes = []

            for table in tables:
                lines = table.strip().split('\n')
                try:
                    df = pd.read_csv(StringIO('\n'.join(lines)), delim_whitespace=True)
                    dataframes.append(df)
                except pd.errors.EmptyDataError:
                    continue

            for df in dataframes:
                is_numeric = pd.to_numeric(df.iloc[-1], errors='coerce').notnull().all()
                if not is_numeric:
                    df.drop(df.tail(1).index,inplace=True)

                if 'energy' in df.columns:
                    df['Iso'] = df['Iso'].astype(float)
                    filtered_df = df[df['Iso'] == 0.5]
                    if not filtered_df.empty:
                        max_energy = filtered_df['delta'].max()
    return f"Max_SP_energy: {max_energy}" if max_energy else "Max SP energy not found"

def main(output_file, deocc, run_val, pair_strength, ipair):
    # Define output paths
    ensurf_file = "ensurf.tst"
    details_file = "out_details.tst"
    summary_file = "summary_results.txt"

    # Read data
    nmass_ncharge = read_nmass_ncharge(ensurf_file)
    energy_info = read_energy(ensurf_file)
    gap_info = read_gaps(details_file)
    max_energy = find_max_energy(details_file)

    # Collect results
    results = f"Deocc: {deocc}, Run: {run_val}, Pair Strength: {pair_strength}, Ipair: {ipair}\n"
    results += f"{nmass_ncharge}\n{energy_info}\n{gap_info}\n{max_energy}\n"

    # Append to summary file
    with open(summary_file, "a") as f:
        f.write(results + "\n")
    print("Results processed and saved to", summary_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process Fortran job results.")
    parser.add_argument("--output", required=True, help="The output file name.")
    parser.add_argument("--deocc", required=True, help="Deoccupation value.")
    parser.add_argument("--run_val", required=True, help="Run value.")
    parser.add_argument("--pair_strength", required=True, help="Pair strength value.")
    parser.add_argument("--ipair", required=True, help="Ipair value.")
    args = parser.parse_args()
    
    main(args.output, args.deocc, args.run_val, args.pair_strength, args.ipair)
