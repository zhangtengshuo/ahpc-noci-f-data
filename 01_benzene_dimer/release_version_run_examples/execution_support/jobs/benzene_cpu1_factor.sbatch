#!/bin/bash
#SBATCH --job-name=AHPC_benzene_cpu1_factor
#SBATCH --partition=noci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/scratch/alpha/AHPC
#SBATCH --output=/scratch/alpha/AHPC/logs/benzene_cpu1_factor-%j.log
set -euo pipefail
srun --cpu-bind=threads bash /scratch/alpha/AHPC/scripts/run_case.sh 01_benzene_dimer cpu1_factor 1 false false
