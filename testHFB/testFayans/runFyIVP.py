import os
import subprocess
import shutil
from pathlib import Path
import re
import sys

# Configuration
# Path to the executable
current_dir = Path(__file__).resolve().parent
print("Current directory:", current_dir)
grandparent_dir = current_dir.parent.parent
print("Grandparent directory:", grandparent_dir)
exe_path = grandparent_dir / "skyax.ifort.omp"
EXECUTABLE = str(exe_path)
print(f"Executable path: {EXECUTABLE}")
print(f"Executable exists: {os.path.exists(EXECUTABLE)}")
BASE_OUTPUT_DIR = current_dir / "testNuclei"

# Nuclear data: {nucleus_name: {"Z": proton_number, "N": neutron_number}}
NUCLEAR_DATA = {
    "Sn120": {"Z": 50, "N": 70},
    "Sn132": {"Z": 50, "N": 82}, 
    "Gd146": {"Z": 64, "N": 82},
    "Pb208": {"Z": 82, "N": 126},
    "Yb158": {"Z": 70, "N": 88}
}

NUCLEI = list(NUCLEAR_DATA.keys())
X0INV_VALUES = [0.1, 0.2, 0.3]
INPUT_FILES = ["for005.FyIVP3FR_nat", "for005.FyIVP3ZR_nat","for005.FyIVP3FR_nonat", "for005.FyIVP3ZR_nonat"]

def modify_input_file(input_file, x0inv_value, z_value, n_value, nucleus_name):
    """Read input file and modify x0inv, nprot, nneut parameters, and title"""
    with open(input_file, 'r') as f:
        content = f.read()
    
    # Replace the title with the current nucleus
    content = re.sub(
        r"title='test [^']*'",
        f"title='test {nucleus_name}'",
        content
    )
    
    # Replace x0inv value
    content = re.sub(
        r'x0inv\s*=\s*[\d.]+D?\d*',
        f'x0inv={x0inv_value}D0',
        content
    )
    
    # Replace nprot (Z) and nneut (N) values on the same line
    content = re.sub(
        r'nprot\s*=\s*\d+\s*,\s*nneut\s*=\s*\d+',
        f'nprot={z_value},nneut={n_value}',
        content
    )
    
    return content

def run_skyax(input_content, output_dir, input_filename):
    """Run skyax_mkl.exe with input file redirection and save output"""
    os.makedirs(output_dir, exist_ok=True)
    
    # Get current working directory to restore later
    original_cwd = os.getcwd()
    
    try:
        # Change to output directory so output files are generated there
        os.chdir(output_dir)
        
        # Write modified input content to the input file
        with open(input_filename, 'w') as f:
            f.write(input_content)
        
        # Run skyax with shell redirection - ensure executable path is properly quoted
        command = f'"{EXECUTABLE}" < "{input_filename}"'
        print(f"    Running command: {command}")
        
        process = subprocess.Popen(
            command,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            cwd=output_dir
        )
        stdout, stderr = process.communicate()
        
        # Save stdout and stderr
        with open("stdout.txt", 'w') as f:
            f.write(stdout)
        if stderr:
            with open("stderr.txt", 'w') as f:
                f.write(stderr)
        
        # List all files generated in output directory
        output_files = [f for f in os.listdir('.') if os.path.isfile(f) and f not in ['stdout.txt', 'stderr.txt', input_filename]]
        if output_files:
            print(f"    Generated files: {', '.join(output_files)}")
        
        return process.returncode == 0
    except Exception as e:
        print(f"Error running {EXECUTABLE}: {e}")
        return False
    finally:
        # Restore original working directory
        os.chdir(original_cwd)

def run_plotting(nucleus_dir, nucleus_name):
    """Run plotConver.py to generate plots for a nucleus"""
    plots_dir = os.path.join(nucleus_dir, "plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    # Get current working directory to restore later
    original_cwd = os.getcwd()
    
    try:
        # Copy plotConver.py to the plots directory
        plotconver_source = os.path.join(os.path.dirname(__file__), "plotConver.py")
        plotconver_dest = os.path.join(plots_dir, "plotConver.py")
        shutil.copy2(plotconver_source, plotconver_dest)
        
        # Change to plots directory
        os.chdir(plots_dir)
        
        # Run plotConver.py
        process = subprocess.Popen(
            [sys.executable, "plotConver.py"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        stdout, stderr = process.communicate()
        
        # Save output
        with open("plot_stdout.txt", 'w') as f:
            f.write(stdout)
        if stderr:
            with open("plot_stderr.txt", 'w') as f:
                f.write(stderr)
        
        if process.returncode == 0:
            print(f"Plots generated successfully for {nucleus_name}")
            # List generated plot files
            plot_files = [f for f in os.listdir('.') if f.endswith('.png')]
            if plot_files:
                print(f"    Generated plots: {', '.join(plot_files)}")
        else:
            print(f"Warning: Plot generation failed for {nucleus_name}")
            
        return process.returncode == 0
    except Exception as e:
        print(f"    Error generating plots for {nucleus_name}: {e}")
        return False
    finally:
        # Restore original working directory
        os.chdir(original_cwd)

def main():
    for nucleus in NUCLEI:
        nucleus_dir = os.path.join(BASE_OUTPUT_DIR, nucleus)
        os.makedirs(nucleus_dir, exist_ok=True)
        
        # Get nuclear parameters
        z_val = NUCLEAR_DATA[nucleus]["Z"]
        n_val = NUCLEAR_DATA[nucleus]["N"]
        
        print(f"\n Processing {nucleus} (Z={z_val}, N={n_val})...")
        
        # Run simulations for all x0inv values
        for x0inv in X0INV_VALUES:
            x0inv_dir = os.path.join(nucleus_dir, f"x0inv{x0inv}")
            os.makedirs(x0inv_dir, exist_ok=True)
            
            print(f"\n  📁 Processing x0inv={x0inv}...")
            
            # Process input files and organize by nat/nonat
            for input_file in INPUT_FILES:
                if not os.path.exists(input_file):
                    print(f"    Warning: {input_file} not found")
                    continue
                
                # Determine subdirectory based on file extension
                if input_file.endswith("_nat"):
                    subdir_name = "nat"
                elif input_file.endswith("_nonat"):
                    subdir_name = "nonat"
                else:
                    # Fallback for files without _nat or _nonat extension
                    subdir_name = "other"
                
                subdir = os.path.join(x0inv_dir, subdir_name)
                os.makedirs(subdir, exist_ok=True)
                
                # Determine file type (FR or ZR) for progress reporting
                file_type = "FR" if "FR" in input_file else "ZR"
                
                modified_input = modify_input_file(input_file, x0inv, z_val, n_val, nucleus)
                success = run_skyax(modified_input, subdir, os.path.basename(input_file))
                
                status = "✓" if success else "✗"
                print(f"    {status} {nucleus} {file_type} {subdir_name} x0inv={x0inv}")
        
        # Generate plots for this nucleus
        print(f"\n Generating plots for {nucleus}...")
        run_plotting(nucleus_dir, nucleus)

if __name__ == "__main__":
    main()