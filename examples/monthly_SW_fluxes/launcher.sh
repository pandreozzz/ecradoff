#!/bin/bash
#
#
# Calculate one month offline fluxes
#

# Change this to the path of ecradoff if needed
ecradoff_path=$(realpath "$0"/../..)
outputs_root="$0/outputs"
python=$(realpath "$ecradoff_path/.venv/bin/python3")


year=$1
month=$2
skipecrad=${3:-0}

yearmonth=$(printf "%04d-%02d" "$year" "$month")
grid="1p0"
expname="vXX"

# Sets the time shift rule for the aerosol fields
# For years before 2000, the shifted setup uses
# aerosols from year+30, otherwise year-30
if (( year <= 2000 ));
then
  timeshift=30 #Plus 30 years for shifted aerosols
else
  timeshift=-30
fi
echo "timeshift: $timeshift years"

# 3-hourly data
timefreq=3

# SLURM job settings
array_chunk_size=35
jobmem_gb=64 #64
time_per_job="01:30:00"
num_cpus=32 #32



yshift_str=$(printf "%+d" "$timeshift")
outputs_dir="${outputs_root}/${expname}_${timefreq}hourly_${yshift_str}yshift_${grid}"
mkdir -p "$outputs_dir"

# Experiments to run (defined in the `offline_ecrad_driver.py`)
declare -a explist_all=("ifsdef" "ndlut" "ndonlyshift" "ndlutshift")
declare -a explist
for exp in "${explist_all[@]}"; do
  basepath="$outputs_dir/ecrad_output_${exp}_${yearmonth}"
  merged_path=${basepath}_merged.nc
  if [ ! -f "$merged_path" ]; then
    explist+=("$exp")
  else
    echo "Merged output for experiment $exp already exists."
  fi
done

if [ ${#explist[@]} -eq 0 ]; then
  echo "All merged outputs already exist. Exiting."
  exit 0
fi


driver_script="./offline_ecrad_driver.py"
driver_script=$(realpath "$driver_script")

datestr=$(date +%Y%m%d_%H%M%S)

workdir="./workdir/$(basename $outputs_dir)_${yearmonth}_${datestr}"
echo "Logs under $workdir"
mkdir -p "$workdir"

ln -s "$ecradoff_path" "$workdir/ecradoff"
ln -s "$driver_script" "$workdir/$(basename "$driver_script")"

wrapper_driver="$workdir/wrapper_driver.sh"

cat << EOF > "$wrapper_driver"
#!/bin/bash
module load intel
this_day=\$((\$SLURM_ARRAY_TASK_ID))
echo "Running job for day \$this_day"
$python $driver_script $yearmonth \$this_day $grid $expname $timeshift $timefreq $outputs_dir
EOF

ini_day=1
end_day=$(cal $month $year | awk 'NF {DAYS = $NF}; END {print DAYS}')

# Is SLURM available?
if ! command -v sbatch &> /dev/null; then
  echo "SLURM not found. Simple local run"
  SLURM_FOUND=0
else
  SLURM_FOUND=1
fi

if [ "$skipecrad" -eq 0 ]; then
  if [ "$SLURM_FOUND" -eq 1 ]; then
    jobid=$(sbatch --parsable --array="${ini_day}-${end_day}%${array_chunk_size}" --wait\
      -t $time_per_job -c $num_cpus \
      -n 1 --chdir="$workdir" \
      --mem="${jobmem_gb}GB" \
      --output="job_output_%a.txt" \
      --error="job_output_%a.txt" \
      "$wrapper_driver")

    driver_job_name=${wrapper_driver##*/}
    nfailed=$(sacct -j "$jobid" -n -o state,JobName | grep -Ev 'COMPLETED' \
      | grep ${driver_job_name:0:8} | wc -l)
    if [ "$nfailed" -gt 0 ]; then
      echo "Calculations failed for $nfailed jobs."
      exit 1
    fi
  else
    for day in $(seq $ini_day $end_day); do
      SLURM_ARRAY_TASK_ID=$day bash "$wrapper_driver"
    done
  fi
fi


jobids=()
for exp in "${explist[@]}"; do
  basepath="$outputs_dir/ecrad_output_${exp}_${yearmonth}"
  merged_path=${basepath}_merged.nc
  echo "Merging output files with base path $basepath"
  shopt -s nullglob
  files=($basepath-*.nc)
  # Check if any files were found
  if [ ${#files[@]} -eq 0 ]; then
    echo "No output files found for merging. Exiting."
    exit 1
  else
    echo "Found ${#files[@]} files for merging."
  fi
  shopt -u nullglob

  cmd="module load cdo && cdo --history -O mergetime ${files[@]} ${merged_path} && rm ${files[@]}"

  echo "Merging outputs for experiment $exp"
  if [ "$SLURM_FOUND" -eq 1 ]; then
    #srun -n 1 -c 1 --mem="${jobmem_gb}GB" \
    #bash -c "$cmd"
    jobid=$(sbatch --parsable \
      -t "00:30:00" -c 1 -n 1 \
      --mem="${jobmem_gb}GB" \
      --wrap="$cmd")
    echo "Submitted merging job with ID $jobid for experiment $exp"
    jobids+=($jobid)

    while true; do
      sleep 5
      remaining=0
      for jobid in "${jobids[@]}"; do
        jobstate=$(sacct -j "$jobid" --format=State --noheader | head -n 1 | tr -d ' ')

        case "$jobstate" in
          COMPLETED)
            ;;
          FAILED|CANCELLED|TIMEOUT|NODE_FAIL*)
            echo "Merging job $jobid failed ($jobstate)"
            ;;
          ""|PENDING|RUNNING)
            ;;
          *)
            echo "Job $jobid has unknown state ($jobstate)"
            remaining=$((remaining + 1))
        esac
        if [ "$jobstate" == "" ] || [ "$jobstate" == "PENDING" ] || [ "$jobstate" == "RUNNING" ]; then
          remaining=$((remaining + 1))
        fi
      done
      if [ "$remaining" -eq 0 ]; then
        echo "All merging jobs completed."
        break
      fi
      echo "Waiting for $remaining merging jobs to finish..."
    done
  else
    bash -c "$cmd"
  fi
done