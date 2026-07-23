import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import io # Needed to handle string data like a file

# --- 1. Setup ---
# Define File Paths 
script_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(script_dir)
folders = ['x0inv0.1', 'x0inv0.2', 'x0inv0.3']
files = ['conver.ZR', 'conver.FR'] # Using uppercase as per your file names
data_paths = {
    folder_name: os.path.join(parent_dir, folder_name) for folder_name in folders
}

# Define the full list of column names for the data section
COL_NAMES = [
    'iter', 'energy', 'variance', 'difference', 'qtotal', 
    'devoc2', 'devmf2', 'rmsp', 'x0inv', 'epair_1', 'epair_2' # epair(1:2) split into two
]

# --- 2. Load Data Function (UPDATED to read ALL columns) ---

def load_full_convergence_data(filepath):
    """
    Loads ALL data columns, skipping header lines and filtering out 'itermcount'.
    """
    HEADER_LINES_TO_SKIP = 6
    
    try:
        with open(filepath, 'r') as f:
            all_lines = f.readlines()
        
        data_lines = all_lines[HEADER_LINES_TO_SKIP:]
        cleaned_lines = [line for line in data_lines if 'itermcount' not in line]
        
        if not cleaned_lines:
            # print(f"Warning: No data found after filtering in {filepath}")
            return None

        data_string = "".join(cleaned_lines)
        
        # NOTE: usecols is removed to read all 11 columns
        df = pd.read_csv(
            io.StringIO(data_string),
            sep='\\s+',
            header=None,
            names=COL_NAMES
        )
        return df
        
    except FileNotFoundError:
        # print(f"Error: Data file not found at {filepath}.")
        return None
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None

# --- 3. Load the Data ---

data_zr = {}
data_fr = {}

for folder_name, full_folder_path in data_paths.items():
    for file in files:
        filepath = os.path.join(full_folder_path, file)
        
        # Check for existence and try .upper() only if not found (safer approach)
        if not os.path.exists(filepath):
             filepath_upper = os.path.join(full_folder_path, file.upper())
             if os.path.exists(filepath_upper):
                 filepath = filepath_upper
             # else, let the load function fail if it truly doesn't exist
        
        df = load_full_convergence_data(filepath)

        if df is not None:
            # Use the folder name for the label and store the data in the correct dict
            if file == 'conver.ZR': 
                data_zr[folder_name] = df
            elif file == 'conver.FR':
                data_fr[folder_name] = df


# --- 4. Plotting Function (Generalized for any column) ---

def plot_convergence_variable(variable_name, zr_data, fr_data):
    """
    Creates the two-panel plot (ZR and FR) for a given variable.
    """
    
    # Check if data for the variable exists in the dataframes
    if not zr_data and not fr_data:
        print(f"Skipping plot for {variable_name}: No data loaded.")
        return

    # Set up the plot figure
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10, 8), sharex=True, dpi=300)
    plt.subplots_adjust(hspace=0.3)
    
    y_label = f'{variable_name.title()} ({("MeV" if variable_name in ["energy", "variance", "epair_1", "epair_2"] else "fm")})'

    # Panel 1: conver.ZR data
    ax1 = axes[0]
    if zr_data:
        for label, df in zr_data.items():
            ax1.plot(df['iter'][2000:], df[variable_name][2000:], label=label, marker='o', markersize=3) # linestyle Added for distinction
    
    ax1.set_title(f'Convergence of {variable_name.title()} (Zero Range)', fontsize=14)
    ax1.set_ylabel(y_label, fontsize=12)
    ax1.set_yscale('log') if variable_name in ['difference', 'variance'] else None
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(loc='best')

    # Panel 2: conver.FR data
    ax2 = axes[1]
    if fr_data:
        for label, df in fr_data.items():
            ax2.plot(df['iter'][2000:], df[variable_name][2000:], label=label, marker='o', markersize=3) # Added linestyle for distinction

    ax2.set_title(f'Convergence of {variable_name.title()} (Finite Range)', fontsize=14)
    ax2.set_xlabel('Iteration', fontsize=12)
    ax2.set_ylabel(y_label, fontsize=12)
    ax2.set_yscale('log') if variable_name in ['difference', 'variance'] else None
    # ax2.set_xscale('log') # Removed log scale unless specifically requested for this variable
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(loc='best')

    # Save the plot with the requested naming convention
    plot_filename = f'comparison_{variable_name}_2000.png'
    plt.savefig(plot_filename)
    plt.close(fig) # Close the figure to free up memory
    print(f"Plot saved as {plot_filename}")

# --- 5. Execution Loop ---

# Variables to be plotted against 'iter'
# VARIABLES_TO_PLOT = ['energy', 'variance', 'difference','rmsp', 'epair_1', 'epair_2']
VARIABLES_TO_PLOT = ['qtotal', 'devoc2', 'devmf2']

for var in VARIABLES_TO_PLOT:
    plot_convergence_variable(var, data_zr, data_fr)

print("All convergence plots have been successfully generated.")