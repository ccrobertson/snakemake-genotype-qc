snakemake all --jobname "{jobid}" --jobs 100 \
		--keep-going \
		--latency-wait 30 \
		--rerun-incomplete \
		--use-conda \
		--use-singularity \
		--singularity-args "--bind $SCRATCH" \
		--printshellcmds \
		> snakemake.log 2>&1 &
