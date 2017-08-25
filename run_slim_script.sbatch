#!/bin/bash
#SBATCH -p long
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH --mem-per-cpu=4G
#SBATCH --ntasks-per-core=1

# change to directory you were in when you submit the job
cd $SLURM_SUBMIT_DIR

: ${PARAMS?Must define PARAMS}
if [ -z "$PARAMS" ]
then
    echo "Must define PARAMS (is empty)."
    exit 1
fi

SEED=$(printf "%06d" $RANDOM); 

echo "Task ID $SLURM_ARRAY_TASK_ID - seed $SEED"
ALL_PARAMS="$PARAMS -s $SEED"
echo "Parameters for SLiM run:
echo "  $ALL_PARAMS"

/usr/bin/time --format='elapsed: %E / kernel: %S / user: %U / mem: %M' \
     /projects/haldane/SLiM/bin/slim $ALL_PARAMS

echo "Done!"

