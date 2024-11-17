### Set up the configs
This is a pipeline to perform preliminary QC on array-based genotyping. It required input files to be array-based genotypes aligned to hg19 and in plink v1.9 binary format. 
Input file paths should be specified in a config file called "config.yaml". An example config file and example input files are provided for reference.

### Set up the software environment
This pipeline runs in a containerized environment using singularity. Running it requires Snakemake version 7 (it is not compatible with Snakemake version 8) and singularity to be available in your work environment. 
If you will be running with slurm, slurm job submission parameters should be defined in a file called "config_cluster.yaml" (see config_cluster_example.yaml as an example).


### Set up the reference data
The pipeline looks for a 1000G plink file
to use for ancestry inference. To get this set up, run this from the pipeline directory.
You may want to do this in an interactive session since it will take a bit of time. The ref files are about 3.9GB.
```
bash download_king_reference_files.sh
```


### Run snakemake pipeline
```
# Load singularity and snakemake7 (current on great lakes as of 2024-11-16)
module load openjdk/18.0.1.1
module load singularity/4.1.3
module load snakemake/7.32.4

# This will submit each job to the slurm scheduler
bash run.sh

# This will run jobs on the login node (not recommended if you have a large data set)
bash run_no_cluster.sh
```





