snakemake all --jobname "{jobid}" --jobs 100 \
		--keep-going \
		--latency-wait 30 \
		--rerun-incomplete \
		--use-conda \
		--use-singularity \
		--singularity-args "--bind $SCRATCH" \
		--printshellcmds \
		--cluster-config config_cluster.yaml \
		--cluster-status slurm_status.py \
		--cluster "sbatch --output {cluster.output} --time {cluster.time} --mem {cluster.mem} --ntasks {cluster.ntasks} --cpus-per-task {cluster.cpus} --account {cluster.account} --parsable" \
		> snakemake.log 2>&1 &
