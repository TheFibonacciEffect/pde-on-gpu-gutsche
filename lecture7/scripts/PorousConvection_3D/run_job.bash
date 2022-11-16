read -p "Do you want to clean everything (y/n)? " answer
case ${answer:0:1} in
    y|Y )
        echo "cleaning everything"
        # clean everything
        git fetch
        git reset --hard
        git clean -x -d -f -n
        sleep 4
        git clean -x -d -f #-n for what if
        git pull
    ;;
    * )
        echo "not cleaning"
    ;;
esac


sbatch ./runPorousConvection3D.bash
sleep 5
squeue -u class202
