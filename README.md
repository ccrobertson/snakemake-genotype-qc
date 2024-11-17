### Set up the configs
This is a pipeline to perform preliminary QC on array-based genotyping. It takes a config file "config.yaml" which specifies input files. The config.yaml file included in this repos is just an example.

### Set up the software environment
Running this pipeline requires Snakemake7 or earlier with slurm job submission parameters defined in "config_cluster.yaml".

I plan to create a docker image for this pipeline, but for now you'll have to make the dependencies available in your local work environment.
Running this pipeline requires a few software programs which can be downloaded as binary files and placed in one of the directories referenced in your PATH:
* plink1.9 (https://www.cog-genomics.org/plink/)
* plink2 (https://www.cog-genomics.org/plink/2.0/)
* king (https://www.kingrelatedness.com/)

In addition to the plink and king programs, you also need to have a conda environment called "renv" with R and the R package "e1071" installed. You can build this with the following commands:
```
mamba create -n renv -c conda-forge r-base
mamba install -n renv -c R r-e1071
```

### Set up the reference data
The pipeline also looks for a 1000G plink file
to use for ancestry inference. To get this set up, run this from the pipeline directory.
You may want to do this in an interactive session since it will take a bit of time. The ref files are about 3.9GB.
```
bash download_king_reference_files.sh
```

### Run snakemake pipeline
```
module load openjdk/18.0.1.1
module load singularity/4.1.3
module load snakemake/7.32.4

bash run.sh
```

### TO DO
* create docker environment
* add example data set for quick run
* add generalized plotting code




