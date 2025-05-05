#!/bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem 0
#SBATCH --time=2-00:00:00
#SBATCH --output=make.out
#SBATCH --error=make.err

# Define an array of directories
directories=(
  "/scratch/ibp5092/simStateSpace"
)

# SIF path
sif_path="/scratch/ibp5092/sif/docs.sif"

# Iterate over directories
for dir in "${directories[@]}"; do
  if [ -d "$dir" ]; then
    echo "Processing: $dir"
    cd "$dir" || exit
    apptainer exec "$sif_path" make all
    apptainer exec "$sif_path" make auto
  fi
done

echo "done"

exit
