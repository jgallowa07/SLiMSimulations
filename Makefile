run:
	mkdir -p ${PARAMS%.E}; sbatch -o"./${PARAMS%.E}/slurm-run-%j.out" run_slim_script.sbatch;	
