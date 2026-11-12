#!/bin/bash
#SBATCH --job-name=AHPC_indo_gpu8
#SBATCH --partition=noci
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --chdir=/scratch/alpha/AHPC
#SBATCH --output=/scratch/alpha/AHPC/logs/indo_gpu8-%j.log
#SBATCH --gres=gpu:v100:8
set -euo pipefail
srun --cpu-bind=threads bash /scratch/alpha/AHPC/scripts/run_case.sh 03_indo_dimer gpu8 8 true false
