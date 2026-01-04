# Snakemake 9.x Execution Guide

This guide provides updated commands for running SVChordinator with Snakemake 9.x.

## Prerequisites

### Install Snakemake 9.x
```bash
# Create a new conda environment (recommended)
conda create -n snakemake9 -c conda-forge -c bioconda snakemake=9.14.5
conda activate snakemake9
```

### Install Executor Plugins

Choose ONE of the following based on your preferred execution method:

**Option 1: SLURM-specific executor (recommended for SLURM clusters)**
```bash
pip install snakemake-executor-plugin-slurm
```

**Option 2: Cluster-generic executor (more flexible, similar to old --cluster)**
```bash
pip install snakemake-executor-plugin-cluster-generic
```

## Execution Methods

### Method 1: Using SLURM Executor Plugin (Recommended)

This is the native SLURM integration for Snakemake 9.

```bash
snakemake \
  --profile workflow/profiles \
  --config-file config/config.yml \
  --executor slurm \
  --jobs 300 \
  --software-deployment-method apptainer \
  --default-resources \
    partition=componc_cpu \
    mem_mb=20000 \
    time=180
```

Or use the pre-configured profile:
```bash
snakemake \
  --profile workflow/profiles \
  --configfile config/config.yml \
  --profile-name config.v8+.yaml
```

### Method 2: Using Cluster-Generic Executor

This method is most similar to the old `--cluster` approach:

```bash
snakemake \
  --profile workflow/profiles \
  --configfile config/config.yml \
  --profile-name config.v8+.cluster-generic.yaml
```

Or manually:
```bash
snakemake \
  --configfile config/config.yml \
  --executor cluster-generic \
  --cluster-generic-submit-cmd "sbatch --partition={resources.partition} --cpus-per-task={threads} --mem={resources.mem_mb} --time={resources.time}" \
  --jobs 300 \
  --software-deployment-method apptainer \
  --default-resources partition=componc_cpu mem_mb=20000 time=180
```

### Method 3: Local Execution (for testing)

For running locally without cluster submission:

```bash
snakemake \
  --configfile config/config.yml \
  --cores 10 \
  --software-deployment-method apptainer
```

## Key Changes from Snakemake 7.x

### Command-line Options

| Old (v7.x)              | New (v9.x)                              |
|-------------------------|------------------------------------------|
| `--use-singularity`     | `--software-deployment-method apptainer` or `--sdm apptainer` |
| `--use-conda`           | `--software-deployment-method conda` or `--sdm conda` |
| `--cluster "sbatch ..."` | `--executor cluster-generic --cluster-generic-submit-cmd "sbatch ..."` |
| `--restart-times 1`     | `--retries 1` |
| No change needed        | Container/singularity directives in workflow still work! |

### Profile Configuration

- Old profiles: `config.yaml`
- New versioned profiles: `config.v8+.yaml` (for Snakemake 8+)
- The workflow now requires `min_version("9.0.0")`

## Example Workflow Runs

### Dry Run (Test Configuration)
```bash
snakemake \
  --configfile config/config.yml \
  --executor slurm \
  --dry-run \
  --printshellcmds
```

### Full Production Run with SLURM
```bash
snakemake \
  --profile workflow/profiles \
  --configfile config/config.yml \
  --profile-name config.v8+.yaml \
  --rerun-incomplete \
  --keep-going
```

### Unlock Working Directory (if needed)
```bash
snakemake --unlock
```

### Generate DAG Visualization
```bash
snakemake \
  --configfile config/config.yml \
  --dag | dot -Tpdf > dag.pdf
```

## Troubleshooting

### Error: "unknown option --use-singularity"
- Solution: Use `--software-deployment-method apptainer` instead

### Error: "executor 'slurm' not found"
- Solution: Install the executor plugin: `pip install snakemake-executor-plugin-slurm`

### Error: "Container runtime not available"
- Solution: Ensure Apptainer/Singularity is installed and in your PATH

### Performance Issues
- Adjust `--jobs`, `--max-jobs-per-second`, and `--latency-wait` parameters
- Consider using `--scheduler greedy` for better resource utilization

## Additional Resources

- [Snakemake 9 Documentation](https://snakemake.readthedocs.io/)
- [Migration Guide](https://snakemake.readthedocs.io/en/stable/getting_started/migration.html)
- [Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/)
