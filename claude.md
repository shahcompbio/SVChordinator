# SVChordinator - Claude Development Guide

## Project Overview

SVChordinator is a Snakemake workflow for merging, annotating, and visualizing structural variants (SVs) from multiple callers. It creates "chord" visualizations in circos plots.

**Key Features:**
- Multi-caller SV merging using MINDA
- Support for both ONT (Oxford Nanopore) and Illumina (ILL) technologies
- Annotation with gene databases (OncoKB, gene annotations)
- Optional genotyping for somatic ONT SVs
- Circos visualization of SVs

## Technology Stack

- **Workflow Manager:** Snakemake 9.x (minimum version 9.0.0)
- **Package Manager:** UV (recommended) or Conda
- **Containerization:** Apptainer/Singularity
- **Cluster Execution:** SLURM (via snakemake-executor-plugin-slurm or cluster-generic)

## Project Structure

```
SVChordinator/
├── config/
│   ├── config.yml              # Main configuration file
│   ├── samples.tsv             # Input samples (path, caller, technology)
│   ├── example_samples.tsv     # Example format
│   └── run_full_workflow.sh    # SLURM submission script
├── workflow/
│   ├── Snakefile               # Main workflow entry point
│   ├── rules/
│   │   ├── common.smk          # Common functions and helpers
│   │   ├── sv_merge.smk        # MINDA merging step
│   │   ├── genotype.smk        # Genotyping with Sniffles2
│   │   ├── annotate_sv.smk     # Gene annotation
│   │   └── visualize.smk       # Circos plotting
│   └── profiles/
│       ├── config.v8+/         # Default SLURM executor profile
│       │   └── config.yaml
│       └── config.v8+.cluster-generic/  # Cluster-generic executor profile
│           └── config.yaml
├── README.md                   # User documentation
├── SNAKEMAKE_V9_EXECUTION.md   # Snakemake 9 migration guide
└── CHANGELOG.md                # Version history
```

## Supported SV Callers

### ONT (Oxford Nanopore) Callers
- nanomonsv
- SAVANA
- Severus
- cuteSV
- Sniffles2

### Illumina Callers
- SvABA
- manta
- GRIPSS

## Configuration Files

### 1. config/config.yml

Main configuration controlling workflow behavior:

```yaml
samples: config/samples.tsv        # Path to samples TSV
out_dir: /path/to/output          # Output directory
sample_name: sample               # Sample identifier
reference: "hg38"                 # hg19 or hg38
filter_bed: /path/to/regions.bed  # BED file for filtering
sv_type: "somatic"                # "germline" or "somatic"
min_callers: 2                    # Minimum callers for consensus

genotype:
  activate: False                 # Enable genotyping (ONT somatic only)
  tumor_bam: /path/to/tumor.bam
  normal_bam: /path/to/normal.bam

annotate:
  activate: True
  oncokb: /path/to/oncokb.tsv
  gene_annotation: /path/to/genes.txt

visualize:
  activate: True
```

### 2. config/samples.tsv

Tab-separated file listing input VCFs:

```
/path/to/nanomonsv.vcf	nanomonsv	ONT
/path/to/savana.vcf	SAVANA	ONT
/path/to/svaba.vcf	SvABA	ILL
/path/to/manta.vcf	manta	ILL
```

**Format:** `vcf_path<TAB>caller<TAB>technology`

### 3. workflow/profiles/

Snakemake profiles define execution parameters. Profiles must be **directories** containing a `config.yaml` file:

- `config.v8+/` - SLURM executor (native Snakemake 9 SLURM integration)
- `config.v8+.cluster-generic/` - Cluster-generic executor (flexible, similar to old --cluster)

**Important:** Profile paths should reference the directory, not the YAML file:
- ✅ Correct: `--workflow-profile workflow/profiles/config.v8+`
- ❌ Wrong: `--workflow-profile workflow/profiles/config.v8+.yaml`

## Workflow Steps

1. **SV Merging** ([workflow/rules/sv_merge.smk](workflow/rules/sv_merge.smk))
   - Uses MINDA to merge SVs across callers
   - Filters by region (if filter_bed provided)
   - Requires minimum number of callers (min_callers)

2. **Genotyping** ([workflow/rules/genotype.smk](workflow/rules/genotype.smk))
   - Optional step for somatic ONT SVs
   - Uses Sniffles2 to extract read support from normal sample
   - Only works for ONT-only, somatic SV calls

3. **Annotation** ([workflow/rules/annotate_sv.smk](workflow/rules/annotate_sv.smk))
   - Annotates SVs with gene information
   - Uses OncoKB cancer gene list
   - Adds gene annotations from reference

4. **Visualization** ([workflow/rules/visualize.smk](workflow/rules/visualize.smk))
   - Creates circos plots showing SVs as "chords"
   - Visualizes breakpoints and connections

## Running the Workflow

### Local Execution (Testing)

```bash
cd /home/preskaa/SVChordinator
uv run snakemake \
  --snakefile workflow/Snakefile \
  --configfile config/config.yml \
  --cores 4 \
  --software-deployment-method apptainer \
  --dry-run
```

### SLURM Execution

Using the submission script:

```bash
bash config/run_full_workflow.sh
```

Or manually with profile:

```bash
uv run snakemake \
  --snakefile workflow/Snakefile \
  --workflow-profile workflow/profiles/config.v8+ \
  --configfile config/config.yml \
  --conda-prefix /path/to/conda \
  --singularity-prefix /path/to/singularity \
  --singularity-args "--bind /data1/shahs3"
```

### Common Snakemake Commands

```bash
# Dry run to check workflow
uv run snakemake --dry-run --printshellcmds

# Unlock directory if stuck
uv run snakemake --unlock

# Clean specific outputs
uv run snakemake --delete-all-output

# Generate workflow DAG
uv run snakemake --dag | dot -Tpdf > dag.pdf

# Run with retries
uv run snakemake --retries 3 --keep-going
```

## Development Guidelines

### Adding New SV Callers

1. Update caller list in documentation
2. Add caller-specific reformatting rules in [workflow/rules/common.smk](workflow/rules/common.smk)
3. Ensure VCF format compatibility with MINDA
4. Update [config/example_samples.tsv](config/example_samples.tsv)

### Modifying Rules

- **Common functions:** [workflow/rules/common.smk](workflow/rules/common.smk)
- **Resource requirements:** Adjust in rule definitions or profile config
- **Container images:** Specified in `container:` directives within rules

### Profile Customization

Edit profile configs to adjust:
- Default resource allocations (mem, time, partition)
- Job limits and scheduling
- Executor-specific parameters

**Example:** Increase default memory in [workflow/profiles/config.v8+/config.yaml](workflow/profiles/config.v8+/config.yaml)

### Rule Execution Order

The workflow includes a rule order declaration:
```python
ruleorder: filter_minda > reformat_minda_no_genotype
```

This resolves ambiguity when genotyping is disabled.

## Important Constraints

1. **Genotyping limitations:**
   - Only works for somatic SV calls
   - Only works with ONT-only samples (no Illumina)
   - Requires tumor and normal BAM files

2. **Reference genomes:**
   - Must be `hg19` or `hg38`
   - Ensure annotations match reference version

3. **Minimum callers:**
   - `min_callers` determines consensus threshold
   - Set appropriately based on available callers

## Troubleshooting

### Profile Directory Error
```
NotADirectoryError: [Errno 20] Not a directory: '.../config.v8+.yaml'
```
**Solution:** Reference the profile directory, not the YAML file directly

### Container Runtime Issues
```
Container runtime not available
```
**Solution:** Ensure Apptainer/Singularity is installed and in PATH

### Missing Executor Plugin
```
executor 'slurm' not found
```
**Solution:** `pip install snakemake-executor-plugin-slurm`

### Genotyping Failures
- Check that SVs are somatic (not germline)
- Verify only ONT callers are used (no Illumina)
- Confirm BAM files are accessible and indexed

## Common File Paths (User: preskaa)

- **Pipeline:** `/home/preskaa/SVChordinator`
- **Conda envs:** `/data1/shahs3/users/preskaa/conda`
- **Singularity cache:** `/data1/shahs3/users/preskaa/singularity`
- **Reference data:** `/data1/shahs3/reference/ref-sarcoma/`
- **Output directory:** `/data1/shahs3/users/preskaa/ThreeByThreeSarcoma/data/`

## Key Migration Notes (Snakemake 7 → 9)

| Old (v7.x)              | New (v9.x)                              |
|-------------------------|------------------------------------------|
| `--use-singularity`     | `--software-deployment-method apptainer` |
| `--use-conda`           | `--software-deployment-method conda`     |
| `--cluster "sbatch ..."` | `--executor cluster-generic --cluster-generic-submit-cmd "sbatch ..."` |
| `--restart-times 1`     | `--retries 1`                            |

See [SNAKEMAKE_V9_EXECUTION.md](SNAKEMAKE_V9_EXECUTION.md) for detailed migration guide.

## Useful References

- [Snakemake 9 Documentation](https://snakemake.readthedocs.io/)
- [MINDA SV Merger](https://github.com/mkirsche/Minda)
- [Sniffles2 Genotyper](https://github.com/fritzsedlazeck/Sniffles)
