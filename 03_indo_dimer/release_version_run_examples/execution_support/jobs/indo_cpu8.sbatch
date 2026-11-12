#!/bin/bash
#SBATCH --job-name=AHPC_indo_cpu8
#SBATCH --partition=noci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --chdir=/scratch/alpha/AHPC
#SBATCH --output=/scratch/alpha/AHPC/logs/indo_cpu8-%j.log
set -euo pipefail
srun --cpu-bind=threads bash /scratch/alpha/AHPC/scripts/run_case.sh 03_indo_dimer cpu8 8 false false
