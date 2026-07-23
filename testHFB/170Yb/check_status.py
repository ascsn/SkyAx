#!/usr/bin/env python3
"""
Check status of completed jobs and create script for remaining ones.
"""
import os
import subprocess
from pathlib import Path

# Parameters from local_job.py
nucleus = "170Yb"
edf = "FyDr"
hfbwidths = ["1.0d0", "1.0d0", "1.0d0", "1.0d0", "1.0d0"]
ipair = [8]
deoccList = ["2.D0", "4.D0", "6.D0", "8.D0", "10.D0", "12.D0", "13.D0", "14.D0", "16.D0", "18.D0"]
runList = [[5, 6], [10, 11], [20, 21], [40, 41], [60, 61]]

job_script_dir = "job_scripts"
outs_dir = "outs"

def get_expected_jobs():
    """Generate list of expected job names based on parameters."""
    jobs = []
    job_count = 0
    for j in ipair:
        for val in deoccList:
            cond = "hfb" if j == 8 else "bcs"
            combined = zip(hfbwidths, runList)
            for hfbwidth, run_val in combined:
                run_start, run_end = run_val
                job_name = f"job_{nucleus}_{cond}_{hfbwidth}_{val}_{run_start}_{run_end}"
                jobs.append(job_name)
                job_count += 1
    return jobs

def get_completed_jobs():
    """Get list of jobs with output files."""
    completed = []
    if os.path.isdir(outs_dir):
        for file in os.listdir(outs_dir):
            if file.endswith("_out.log"):
                job_name = file.replace("_out.log", "")
                completed.append(job_name)
    return sorted(completed)

def main():
    expected = get_expected_jobs()
    completed = get_completed_jobs()
    
    print("="*80)
    print(f"Job Status Report for {nucleus} ({edf})")
    print("="*80)
    print(f"\nExpected jobs: {len(expected)}")
    print(f"Completed jobs: {len(completed)}")
    print(f"Remaining jobs: {len(expected) - len(completed)}")
    
    remaining = [job for job in expected if job not in completed]
    
    if remaining:
        print(f"\n❌ INCOMPLETE - {len(remaining)} jobs still need to run:")
        for job in remaining:
            print(f"  - {job}")
    else:
        print("\n✓ ALL JOBS COMPLETED!")
    
    print("\n" + "="*80)
    print("Completed jobs by deocc value:")
    print("="*80)
    
    # Group by deocc value
    by_deocc = {}
    for job in completed:
        parts = job.split("_")
        deocc_val = parts[4]  # Extract deocc value
        if deocc_val not in by_deocc:
            by_deocc[deocc_val] = 0
        by_deocc[deocc_val] += 1
    
    for deocc in deoccList:
        count = by_deocc.get(deocc, 0)
        total = len([j for j in expected if f"_{deocc}_" in j])
        status = "✓" if count == total else "✗"
        print(f"  deocc={deocc}: {count}/{total} {status}")
    
    # Generate script to run remaining jobs
    if remaining:
        script_content = f"""#!/bin/bash
# Script to run remaining {len(remaining)} jobs for {nucleus}

set -e  # Exit on error

remaining_jobs=(
"""
        for job in remaining:
            script_content += f'    "{job}"\n'
        
        script_content += """)

echo "Running ${#remaining_jobs[@]} remaining jobs..."
echo "=================="

for job_name in "${remaining_jobs[@]}"; do
    echo "Starting: $job_name"
    script="${job_name}.sh"
    
    if [ -f "${script_content.split(chr(10))[0].replace(chr(39), '')}job_scripts/$script" ]; then
        bash "job_scripts/$script"
        echo "Completed: $job_name"
    else
        echo "ERROR: Script not found: job_scripts/$script"
    fi
done

echo "=================="
echo "All remaining jobs completed!"
"""
        
        script_path = "run_remaining_jobs.sh"
        with open(script_path, "w") as f:
            f.write(f"""#!/bin/bash
# Script to run remaining {len(remaining)} jobs for {nucleus}

export OMP_STACKSIZE=4500m
ulimit -s unlimited

echo "Running {len(remaining)} remaining jobs..."
echo "=================="

remaining_jobs=(
""")
            for job in remaining:
                f.write(f'    "{job}"\n')
            f.write(""")

for job_name in "${remaining_jobs[@]}"; do
    echo "Starting: $job_name"
    bash "job_scripts/${job_name}.sh"
    echo "Completed: $job_name"
    echo ""
done

echo "=================="
echo "All remaining jobs completed!"
""")
        os.chmod(script_path, 0o755)
        print(f"\nCreated: {script_path}")
        print(f"Run with: bash {script_path}")
    else:
        print("\n✓ No script needed - all jobs are complete!")

if __name__ == "__main__":
    main()
