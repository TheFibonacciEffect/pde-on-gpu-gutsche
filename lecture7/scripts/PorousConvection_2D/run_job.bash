# clean everything
git fetch
git reset --hard
git clean -x -d -f -n
sleep 4
git clean -x -d -f #-n for what if

git pull
sbatch ./runPorousConvection.bash
sleep 5
squeue -u class202
