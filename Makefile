run:
	mkdir -p ./Output1/$${PARAMS%.E};
	sbatch -o"./Output1/$${PARAMS%.E}/slurm-run-%j.out" run_slim_script.sbatch;	


