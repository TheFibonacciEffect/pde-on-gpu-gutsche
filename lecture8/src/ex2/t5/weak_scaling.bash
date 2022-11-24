#!/bin/bash -l
# write for loop
for i in 1 2 4
do
    # run the job
    sbatch --nodes=$i weak_scaling_job.bash $i
done

# julia --project=../../.. "plot result.jl"