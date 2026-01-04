#!/bin/bash
#SBATCH --partition=componc_cpu
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=17.1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=preskaa@mskcc.org
#SBATCH --output=slurm%j_snkmk.out

### example slurm submission script ###

## set directories
pipeline_dir=$HOME/SVChordinator
outdir=/data1/shahs3/users/preskaa/SarcAtlas/data/APS046_ont_ill_sv_analysis/TCDO-SAR-061_svchordinator/
config_yaml=config/config.yml
workflow_profile=${pipeline_dir}/workflow/profiles/config.v8+.cluster-generic
snakefile=${pipeline_dir}/workflow/Snakefile
## switch to the right conda environment
source /home/preskaa/miniforge3/bin/activate snakemake

mkdir -p ${outdir}

### run the pipeline
echo "Current working directory: $(pwd)"
echo "Snakefile path: ${snakefile}"

# this will run the process, for example you could switch the below to a python command:
# python script.py
cd ${pipeline_dir}

uv run snakemake \
  --snakefile ${snakefile} \
  --workflow-profile ${workflow_profile} \
  --configfile ${config_yaml}\
  --conda-prefix /data1/shahs3/users/preskaa/conda \
  --singularity-prefix /data1/shahs3/users/preskaa/singularity \
  --singularity-args "--bind /data1/shahs3" \
  #--dry-run
