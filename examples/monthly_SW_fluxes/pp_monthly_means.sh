#!/bin/bash
#
#
# Merge all times in a monthly mean file

outfile_dir="${1%/}"

nprocs=16
jobmem_gb=64
time_per_job="01:30:00"

module load cdo

# Look for all _merged.nc files in the directory
merged_files=("$outfile_dir"/*_merged.nc)
# Sort the files by name
IFS=$'\n' sorted_files=($(sort <<<"${merged_files[*]}"))
unset IFS


merger_function() {
  local exp_type=$1
  local files=("$outfile_dir"/ecrad_output_"$exp_type"_????-??_merged.nc)
  IFS=$'\n' files=($(sort <<<"${files[*]}"))
  unset IFS
  local output_file="$outfile_dir/ecrad_output_${exp_type}_monthly_means.nc"

  if [ ${#files[@]} -eq 0 ]; then
    echo "No files found for experiment type $exp_type. Skipping."
    return
  fi

  echo "Merging ${#files[@]} files for experiment type $exp_type into $output_file"
  local runner_script
  runner_script=$(mktemp "$outfile_dir/.${exp_type}_merge_XXXXXX.sh")

  {
    printf '#!/bin/bash\n'
    printf 'module load cdo\n'
    printf 'output_file=%q\n' "$output_file"
    # Run cdo jobs in batches of nprocs
    printf 'failed=0\n'
    local file
    local batch_count=0
    local batch_inputs=()
    local batch_outputs=()
    for file in "${files[@]}"; do
      local tmp="${file%.nc}_monthly.nc"
      batch_inputs+=("$file")
      batch_outputs+=("$tmp")
      batch_count=$((batch_count + 1))
      if [ "$batch_count" -ge "$nprocs" ]; then
        # emit parallel launch of this batch
        local i
        for (( i=0; i<${#batch_inputs[@]}; i++ )); do
          printf 'cdo -setday,15 -monmean %q %q &\n' "${batch_inputs[$i]}" "${batch_outputs[$i]}"
        done
        printf 'wait\n'
        for (( i=0; i<${#batch_outputs[@]}; i++ )); do
          printf '[ -f %q ] || { echo "ERROR: not created: %q" >&2; failed=1; }\n' "${batch_outputs[$i]}" "${batch_outputs[$i]}"
        done
        batch_count=0
        batch_inputs=()
        batch_outputs=()
      fi
    done
    # remaining files (last partial batch)
    if [ "${#batch_inputs[@]}" -gt 0 ]; then
      local i
      for (( i=0; i<${#batch_inputs[@]}; i++ )); do
        printf 'cdo -setday,15 -monmean %q %q &\n' "${batch_inputs[$i]}" "${batch_outputs[$i]}"
      done
      printf 'wait\n'
      for (( i=0; i<${#batch_outputs[@]}; i++ )); do
        printf '[ -f %q ] || { echo "ERROR: not created: %q" >&2; failed=1; }\n' "${batch_outputs[$i]}" "${batch_outputs[$i]}"
      done
    fi
    printf 'if [ "$failed" -ne 0 ]; then echo "One or more cdo jobs failed" >&2; exit 1; fi\n'
    printf 'tmp_files=(\n'
    for file in "${files[@]}"; do
      printf '  %q\n' "${file%.nc}_monthly.nc"
    done
    printf ')\n'
    printf 'for f in "${tmp_files[@]}"; do [ -f "$f" ] || { echo "Missing expected output: $f" >&2; exit 1; }; done\n'
    printf 'cdo -O mergetime "${tmp_files[@]}" "$output_file"\n'
    printf 'rm -f "${tmp_files[@]}"\n'
  } > "$runner_script"

  chmod +x "$runner_script"

  srun -t "$time_per_job" -c ${nprocs} -n 1 --mem="${jobmem_gb}GB" bash "$runner_script"
  rm -f "$runner_script"
}

# File names are of the form ecrad_output_{exptype}_YYYY-MM_merged.nc
current_exp_type=""
ncurrent_exp_type=0
for f in "${sorted_files[@]}"; do
  # Extract the experiment type from the file name
  this_exp_type=$(basename "$f" | cut -d'_' -f3)
  this_yearmonth=$(basename "$f" | cut -d'_' -f4 | cut -d'-' -f1-2)
  if [ "$this_exp_type" != "$current_exp_type" ]; then
    # If we have a new experiment type, merge the previous one (if any)
    if [ -n "$current_exp_type" ]; then
      echo "Found ${ncurrent_exp_type} files for experiment type $current_exp_type"
      merger_function "$current_exp_type"
    fi
    current_exp_type="$this_exp_type"
    ncurrent_exp_type=1
    echo "$current_exp_type $this_yearmonth"
  else
    echo "$current_exp_type $this_yearmonth"
    ncurrent_exp_type=$((ncurrent_exp_type + 1))
  fi
done

# Merge the last experiment type
if [ -n "$current_exp_type" ]; then
  echo "Found ${ncurrent_exp_type} files for experiment type $current_exp_type"
  merger_function "$current_exp_type"
fi