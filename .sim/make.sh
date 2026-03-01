#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=2-00:00:00
#SBATCH --output=make.out
#SBATCH --error=make.err

# Define an array of directories
directories=(
  "/home/ibp5092/Downloads/simStateSpace"
)

# SIF path
sif_path="/home/ibp5092/work/docs.sif"

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
