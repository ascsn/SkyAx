#!/bin/bash
# Job script for job_170Yb_hfb_1.0d0_14.D0_10_11

echo "Running job: job_170Yb_hfb_1.0d0_14.D0_10_11"
echo "Input file: for005/for005_36"

# Make sure skyax executable exists
if [ ! -f "../../skyax.ifort.omp" ]; then
    echo "Error: ../../skyax.ifort.omp not found. Please build it first."
    exit 1
fi

# Set environment variables for OpenMP
export OMP_STACKSIZE=4500m
ulimit -s unlimited

# Run the calculation
../../skyax.ifort.omp < for005/for005_36 > outs/job_170Yb_hfb_1.0d0_14.D0_10_11_out.log

echo "Job job_170Yb_hfb_1.0d0_14.D0_10_11 completed. Output in outs/job_170Yb_hfb_1.0d0_14.D0_10_11_out.log"
