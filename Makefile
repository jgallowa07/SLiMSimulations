run:
	mkdir -p ./Output/${PARAMS%.E}; sbatch -o"./Output/${PARAMS%.E}/slurm-run-%j.out" run_slim_script.sbatch;	


