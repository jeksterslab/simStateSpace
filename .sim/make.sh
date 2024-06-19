#! /bin/bash
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem 0
#SBATCH --time=2-00:00:00
#SBATCH --output=make-simStateSpace.out
#SBATCH --error=make-simStateSpace.err

cd /scratch/ibp5092/simStateSpace
apptainer exec /scratch/ibp5092/sif/docs-mplus.sif make all
apptainer exec /scratch/ibp5092/sif/docs-mplus.sif make auto
