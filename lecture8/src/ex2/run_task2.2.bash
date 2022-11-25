#!/bin/bash -l
#SBATCH --job-name="Diff2D_xpu"
#SBATCH --output=Diff2D_xpu.%j.o
#SBATCH --error=Diff2D_xpu.%j.e
#SBATCH --time=00:05:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account class04

module load daint-gpu
module load Julia/1.7.2-CrayGNU-21.09-cuda
srun julia -O3 --check-bounds=no --project=../..  l8_diffusion_2D_perf_xpu.jl true
srun julia -O3 --check-bounds=no --project=../..  l8_diffusion_2D_perf_xpu.jl false