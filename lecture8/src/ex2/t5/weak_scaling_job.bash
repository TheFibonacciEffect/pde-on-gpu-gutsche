#!/bin/bash -l
#SBATCH --job-name="weak_scaling"
#SBATCH --output=weak_scaling.%j.o
#SBATCH --error=weak_scaling.%j.e
#SBATCH --time=00:00:20
#SBATCH --ntasks-per-node=1
#SBATCH --partition=normal
#SBATCH --constraint=gpu
#SBATCH --account class04

module load daint-gpu
module load Julia/1.7.2-CrayGNU-21.09-cuda

export MPICH_RDMA_ENABLED_CUDA=1
export IGG_CUDAAWARE_MPI=1

srun -n4 bash -c 'LD_PRELOAD="/usr/lib64/libcuda.so:/usr/local/cuda/lib64/libcudart.so" julia -O3 --check-bounds=no --project=../../.. weak_scaling.jl'
srun julia --project=../../.. "plot result.jl"
