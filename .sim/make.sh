#! /bin/bash
#SBATCH --job-name=make
#SBATCH --partition=basic
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=04:00:00
#SBATCH --output=make.out
#SBATCH --error=make.err

###############################################################################
# Usage (on the submit node, from .sim/):
#   sbatch make.sh                         # build ~/Downloads/simStateSpace
#   sbatch --account=open make.sh          # charge an account
#   PROJECT=/path/to/clone sbatch make.sh  # build another clone
#
# Description:
#   make all, then make auto (PUB) in a fresh clone of the private repository.
#   make all is serial (compile, vignettes, pkgdown, coverage), so one CPU
#   and one thread per library keep the CPU busy; more CPUs would sit idle.
#   Sizes (2026-10-17, docs.sif, job 56531470): make all to pkgdown took
#   14 min, CPU 94%, peak 1.66 GB (coverage and make auto not yet run).
#   After a run, check with: seff <jobid>, and lower --mem to the peak
#   plus 10% (or raise it if the job ran out of memory).
###############################################################################

set -euo pipefail

# Define project variables
PROJECT=${PROJECT:-${HOME}/Downloads/simStateSpace}
SIF=${SIF:-${HOME}/work/sif/docs.sif}
# ~/work links outside the home folder; bind it so that git in the container
# finds the SSH keys (make auto clones and pushes)
BIND=${BIND:-$(readlink -f "${HOME}/work")}

# One thread per library, one make job
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MAKEFLAGS=-j1

if [ ! -d "${PROJECT}" ]; then
  echo "Not found: ${PROJECT}" >&2
  exit 1
fi

echo "Processing: ${PROJECT}"
cd "${PROJECT}"
apptainer exec --bind "${BIND}" "${SIF}" make all
apptainer exec --bind "${BIND}" "${SIF}" make auto

echo "done"
