#!/bin/bash --login

# Parameters
nucleus="120Sn"
edf="SVbas"
hfbwidths=("0.0d0" "0.0d0" "0.0d0" "0.0d0" "0.0d0" "0.0d0")
pairStrengthsBCS=("537" "478" "451" "438" "433" "432")
pairStrengthsHFB=("537" "478" "451" "438" "433" "432")
ipair=(8)
deoccList=("18.D0")
runList=("10 11" "20 21" "40 41" "60 61" "80 81" "100 101")

# Directory to store job scripts
job_script_dir="job_scripts"
mkdir -p $job_script_dir

create_and_submit_job_script() {
    job_name=$1
    hfbwidth=$2
    pairStrength=$3
    deocc_val=$4
    ipair_val=$5
    run_start=$6
    run_end=$7

    # SLURM job script content
    job_script_content="#!/bin/bash --login

#SBATCH --job-name=${job_name}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --propagate=CORE
#SBATCH --time=05:59:59
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-100
#SBATCH -A ptg
#SBATCH -e /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/error_${job_name}.err
#SBATCH -o /mnt/home/lalitsud/Witek/SkyAx/slurmjobs/errs/output_${job_name}.out

module purge
module load powertools
module load intel/2023b
module load imkl/2023.2.0

# Prepare input file
echo \"&main outnam='tst',title='test ${nucleus}',iprint=1000,itmax=1000\" > for005
echo \"/ \" >> for005
echo \"&force \" >> for005
echo \" name='none',pairstab=-0D0,HFBwidth=${hfbwidth},\" >> for005
echo \" t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921,\" >> for005
echo \" x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530,\" >> for005
echo \" t4=124.633,b4p=34.1116989,power=0.3,\" >> for005
echo \" h2m=20.749811,20.721249,\" >> for005
echo \" chg=1D0,ex=1,riftm=0D0,zpe=1,\" >> for005
echo \" ipair=${ipair_val},v0prot=10,v0neut=${pairStrength},rho0pr=0.2011591,\" >> for005
echo \" pair_cutoff=-${run_start},tregpair=F\" >> for005
echo \"/ \" >> for005
echo \"&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0,\" >> for005
echo \"/ \" >> for005
echo \"&static nprot=50,nneut=70,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0,\" >> for005
echo \"deocc(1)=${deocc_val},deocc(2)=${deocc_val},errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0,\" >> for005
echo \"/ \" >> for005
echo \"&constraint\" >> for005
echo \"/ \" >> for005
echo \"&thermal \" >> for005
echo \"/ \" >> for005

# Run the program
ulimit -s unlimited
export OMP_STACKSIZE=4500m
cd ../source_f90/
make clean
make ifort_omp
cd ../slurmjobs
../skyax.ifort.omp < for005 > ${job_name}_output.log

# Post-process output
source ~/anaconda3/bin/activate
conda activate slpy
python3 process_results.py --output ${job_name}_output.log --deocc ${deocc_val} --run_val ${run_start} --pair_strength ${pairStrength} --ipair ${ipair_val}"

    # Save the script to a file
    job_script_path="${job_script_dir}/${job_name}.sh"
    echo "$job_script_content" > $job_script_path

    # Submit the job script
    sbatch $job_script_path
}

# Loop through parameter combinations and submit jobs
for j in "${ipair[@]}"; do
    for val in "${deoccList[@]}"; do
        if [[ $j -eq 8 ]]; then
            cond="hfb"
            pair_strength_list=("${pairStrengthsHFB[@]}")
        else
            cond="bcs"
            pair_strength_list=("${pairStrengthsBCS[@]}")
        fi

        for hfbwidth in "${hfbwidths[@]}"; do
            for pair_strength in "${pair_strength_list[@]}"; do
                for run_val in "${runList[@]}"; do
                    run_start=$(echo $run_val | cut -d ' ' -f 1)
                    run_end=$(echo $run_val | cut -d ' ' -f 2)
                    job_name="job_${cond}_${hfbwidth}_${pair_strength}_${val}_${run_start}_${run_end}"
                    create_and_submit_job_script $job_name $hfbwidth $pair_strength $val $j $run_start $run_end
                done
            done
        done
    done
done

echo "All jobs submitted!"
