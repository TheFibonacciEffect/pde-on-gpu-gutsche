# write for loop
for i in 1 2 4 8 16 32 64
do
    # run the job
    sbatch  --nodes=$i weak_scaling_job.bash
done