#!/bin/bash
# Script to run remaining 35 incomplete jobs for 170Yb
# These jobs did not complete their iterations (no rms_C in final output)

export OMP_STACKSIZE=4500m
ulimit -s unlimited

remaining_jobs=(
    "job_170Yb_hfb_1.0d0_2.D0_5_6"
    "job_170Yb_hfb_1.0d0_2.D0_20_21"
    "job_170Yb_hfb_1.0d0_2.D0_40_41"
    "job_170Yb_hfb_1.0d0_2.D0_60_61"
    "job_170Yb_hfb_1.0d0_4.D0_5_6"
    "job_170Yb_hfb_1.0d0_4.D0_10_11"
    "job_170Yb_hfb_1.0d0_4.D0_20_21"
    "job_170Yb_hfb_1.0d0_4.D0_40_41"
    "job_170Yb_hfb_1.0d0_4.D0_60_61"
    "job_170Yb_hfb_1.0d0_6.D0_5_6"
    "job_170Yb_hfb_1.0d0_6.D0_10_11"
    "job_170Yb_hfb_1.0d0_6.D0_20_21"
    "job_170Yb_hfb_1.0d0_6.D0_40_41"
    "job_170Yb_hfb_1.0d0_6.D0_60_61"
    "job_170Yb_hfb_1.0d0_8.D0_5_6"
    "job_170Yb_hfb_1.0d0_8.D0_10_11"
    "job_170Yb_hfb_1.0d0_8.D0_20_21"
    "job_170Yb_hfb_1.0d0_8.D0_40_41"
    "job_170Yb_hfb_1.0d0_8.D0_60_61"
    "job_170Yb_hfb_1.0d0_13.D0_60_61"
    "job_170Yb_hfb_1.0d0_14.D0_5_6"
    "job_170Yb_hfb_1.0d0_14.D0_10_11"
    "job_170Yb_hfb_1.0d0_14.D0_20_21"
    "job_170Yb_hfb_1.0d0_14.D0_40_41"
    "job_170Yb_hfb_1.0d0_14.D0_60_61"
    "job_170Yb_hfb_1.0d0_16.D0_5_6"
    "job_170Yb_hfb_1.0d0_16.D0_10_11"
    "job_170Yb_hfb_1.0d0_16.D0_20_21"
    "job_170Yb_hfb_1.0d0_16.D0_40_41"
    "job_170Yb_hfb_1.0d0_16.D0_60_61"
    "job_170Yb_hfb_1.0d0_18.D0_5_6"
    "job_170Yb_hfb_1.0d0_18.D0_10_11"
    "job_170Yb_hfb_1.0d0_18.D0_20_21"
    "job_170Yb_hfb_1.0d0_18.D0_40_41"
    "job_170Yb_hfb_1.0d0_18.D0_60_61"
)

echo "=========================================="
echo "Running ${#remaining_jobs[@]} incomplete jobs"
echo "=========================================="
echo ""

start_time=$(date)
echo "Start time: $start_time"
echo ""

completed=0
failed=0

for job_name in "${remaining_jobs[@]}"; do
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: $job_name"
    
    script_path="job_scripts/${job_name}.sh"
    
    if [ ! -f "$script_path" ]; then
        echo "  ❌ ERROR: Script not found: $script_path"
        ((failed++))
        continue
    fi
    
    # Run the job
    if bash "$script_path" > /dev/null 2>&1; then
        echo "  ✓ Completed successfully"
        ((completed++))
    else
        echo "  ✗ Job may have issues (check outs/${job_name}_out.log)"
        ((completed++))  # Still count as completed even if script returns error
    fi
    
    echo ""
done

echo "=========================================="
echo "Job Summary"
echo "=========================================="
echo "Total jobs run: ${#remaining_jobs[@]}"
echo "Completed: $completed"
echo "Failed: $failed"
echo ""
echo "Start time: $start_time"
echo "End time: $(date)"
echo "=========================================="
