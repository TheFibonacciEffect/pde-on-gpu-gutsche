#!/bin/bash -l
#SBATCH --job-name="my_gpu_run"
#SBATCH --output=my_gpu_run.%j.o
#SBATCH --error=my_gpu_run.%j.e
#SBATCH --time=02:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account class04
#sbatch submit.sh

module load daint-gpu
module load Julia/1.7.2-CrayGNU-21.09-cuda
srun julia -O3 --check-bounds=no --project=../.. ./PorousConvection_2D_xpu_daint.jl