#! /bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem 0
#SBATCH --time=2-00:00:00
#SBATCH --output=make.out
#SBATCH --error=make.err

cd /scratch/ibp5092/simStateSpace || exit
apptainer exec /scratch/ibp5092/sif/docs.sif make all
apptainer exec /scratch/ibp5092/sif/docs.sif make auto
