# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Updated minimum Snakemake version requirement from 7.31.1 to 9.0.0
- Migrated cluster profiles to Snakemake 9.x standards
- Added new execution profiles for SLURM and cluster-generic executors

### Added
- Created Snakemake 9.x compatible profile configurations
- Added comprehensive execution guide (SNAKEMAKE_V9_EXECUTION.md)
- Added versioned profile configurations (config.v8+.yaml)

## [v0.0.7-beta] - 2025-04-08

### Fixed
- Fixed bug in MINDA integration

## [v0.0.6-beta] - 2025-03-14

### Added
- Multimodal SV calling with ONT and Illumina (#1)
- Introduced germline calling functionality
- Support for both somatic and germline SV detection
- Example VCF table included

### Changed
- Updated README.md with multimodal calling documentation
- Configuration file edits for multimodal support

### Fixed
- Fixes for weird strand annotations
- Option to reduce minimum number of callers for consensus calls

## [v0.0.5-beta] - 2025-01-28

### Changed
- Cleaned up log output
- Updated config with reference genome option (hg19/hg38)

### Fixed
- Bug fixes for hg37/hg19 reference genome support
- Bug fix for reading hg37 VCFs
- Bug fix for ONT genotyping

## [v0.0.4-beta.1] - 2025-01-27

### Fixed
- Bug fixes for hg37 reference genome compatibility
- Fixed issues with variant reading for hg37/hg19

## [v0.0.4-beta] - 2024-12-02

### Added
- Updates to capture strand information from individual callers
- Support for viola-sv for Illumina variant parsing

### Changed
- Docker file updates for viola-sv integration
- Enhanced strand annotation capture across multiple callers

### Fixed
- Small config file edits and corrections

## [v0.0.3-beta] - 2024-11-25

### Fixed
- Translocation bug fix
- Fix to strand annotation for reconstruction plots

### Changed
- Allow for retries in workflow execution
- Branch cleanup

## [v0.0.2-beta] - 2024-11-22

### Added
- SV type annotation functionality
- Retry capability for failed jobs

### Fixed
- Variant annotation fixes across multiple rules
- Improved handling of complex SV types

### Changed
- Enhanced annotation pipeline for structural variants

## [v0.0.1-beta] - 2024-11-15

### Added
- Initial release of SVChordinator
- Core SV merging functionality using MINDA
- Genotyping with Sniffles2 for somatic ONT SVs
- Gene annotation pipeline
- OncoKB integration for cancer gene annotations
- Circos plot visualization
- Filtering based on read support from normal samples
- Support for multiple SV callers:
  - ONT: nanomonsv, SAVANA, Severus, cuteSV, Sniffles2
  - Illumina: SvABA, manta, GRIPSS

### Features
- VCF merging across multiple SV callers
- Optional genotyping to extract read support
- Gene and oncogene annotation
- Configurable minimum caller support threshold
- SLURM cluster execution support
- Containerized workflow with Docker/Singularity

## [v0.0.0] - 2024-11-09

### Added
- Initial project structure
- Basic Snakemake workflow skeleton
- Initial commit and repository setup

---

## Release Notes

### Beta Release Information

All releases are currently in beta. The software is under active development and APIs may change between releases.

### Compatibility Notes

- **v0.0.7-beta and earlier**: Snakemake 7.31.1
- **Unreleased/Next**: Snakemake 9.0.0+

### Migration Guide

When upgrading to the next release (Snakemake 9.x compatible):
1. Install Snakemake 9.x: `conda install -c conda-forge -c bioconda snakemake=9.14.5`
2. Install required executor plugin: `pip install snakemake-executor-plugin-slurm`
3. Update execution commands (see SNAKEMAKE_V9_EXECUTION.md)
4. Use new profile configurations in `workflow/profiles/config.v8+.yaml`

[Unreleased]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.7-beta...HEAD
[v0.0.7-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.6-beta...v0.0.7-beta
[v0.0.6-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.5-beta...v0.0.6-beta
[v0.0.5-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.4-beta.1...v0.0.5-beta
[v0.0.4-beta.1]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.4-beta...v0.0.4-beta.1
[v0.0.4-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.3-beta...v0.0.4-beta
[v0.0.3-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.2-beta...v0.0.3-beta
[v0.0.2-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.1-beta...v0.0.2-beta
[v0.0.1-beta]: https://github.com/shahcompbio/SVChordinator/compare/v0.0.0...v0.0.1-beta
[v0.0.0]: https://github.com/shahcompbio/SVChordinator/releases/tag/v0.0.0
