#!/bin/bash
#SBATCH --job-name=my_array_job
#SBATCH --array=1-10             # This submits 10 jobs (1 through 10)
#SBATCH --time=00:10:00          # Specify the time limit for each job
#SBATCH --ntasks=1               # Number of tasks per job
#SBATCH --cpus-per-task=1        # Number of CPU cores per task
#SBATCH --mem=4G                 # Memory per job
#SBATCH --output=output_%A_%a.out    # Output file for each job
#SBATCH --error=error_%A_%a.err      # Error file for each job

# Load any necessary modules or environments
module load my_fortran_module

# Run the executable
./my_fortran_program

