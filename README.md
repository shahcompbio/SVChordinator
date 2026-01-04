# SVChordinator

Merge, annotate, and visualize structural variants (as "chords" in a circos plot).

## Installation

### Quick Start with UV (Recommended)

```bash
# Install uv (if not already installed)
curl -LsSf https://astral.sh/uv/install.sh | sh

# Clone the repository
git clone https://github.com/shahcompbio/SVChordinator.git
cd SVChordinator

# Create environment and install all dependencies
uv sync

# Verify installation
uv run snakemake --version  # Should show 9.x.x
```

## Supported SV Callers

Currently supports the following callers:

- `SvABA`
- `manta`
- `GRIPSS`
- `nanomonsv`
- `SAVANA`
- `Severus`
- `cuteSV`
- `Sniffles2`

## Configuration

Input TSV file with paths to VCF must follow format of the `config/example_samples.tsv` file,
in which callers are labeled with the names of the tools
