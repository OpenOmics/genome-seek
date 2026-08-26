# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.9.4] - 2026-08-26

### Fixed

- Prevented Purple driver analysis failures by filtering somatic variants with missing reference or alternate allele depth (`AD`) values or missing read depth (`DP`).
- Fixed ClairS-TO VCF sample naming so output files use the tumor sample identifier instead of the generic `SAMPLE` name.

See [#27](https://github.com/OpenOmics/genome-seek/pull/27).

## [0.9.3] - 2025-10-27

### Fixed

- Updated PAVE processing of SAGE tumor-only calls to read only passing variants, preventing filtered variants from being included in downstream processing.

## [0.9.2] - 2025-07-22

### Added

- Added CPU and GPU execution paths for the DeepVariant `call_variants` step.
- Added GPU-accelerated DeepVariant execution when the `--use-gpus` option is enabled.

### Changed

- Updated DeepVariant output handling to support the sharded TFRecord files produced by GPU execution.
- Added separate resource allocations for CPU and GPU DeepVariant jobs.

## [0.9.1] - 2025-03-17

### Added

- Added HMF Tools SAGE somatic variant calling for tumor-normal and tumor-only samples.
- Added ClairS-TO somatic variant calling for tumor-only samples.
- Added DeepSomatic somatic variant calling for tumor-normal samples, with CPU and GPU execution support.
- Added the `--use-gpus` option for GPU acceleration of supported workflow steps.
- Added PAVE annotation and panel-of-normals filtering for SAGE tumor-only calls.
- Added AMD and Intel processor detection for selecting a compatible BWA-MEM2 binary.
- Added Zenodo citation and DOI information to the project documentation.

### Changed

- Replaced Strelka with DeepSomatic in the tumor-normal somatic variant workflow.
- Updated the somatic caller ensemble:
  - Tumor-normal samples use Mutect2, Octopus, SAGE, MuSE, and DeepSomatic.
  - Tumor-only samples use Mutect2, Octopus, SAGE, and ClairS-TO.
- Required tumor-only somatic variants to be supported by at least two callers.
- Updated DeepSomatic to version 1.8.0.
- Restricted DeepVariant and DeepSomatic WES processing to the configured capture regions.
- Improved WES BED normalization to produce valid numeric scores and strand annotations.
- Removed an additional Octopus genotype-level filter that reduced recall.
- Removed the Mutect2 contamination-table filter after it was found to reduce recall.
- Updated ClairS-TO resource allocations and region handling.
- Adjusted Octopus resource allocations for observed memory requirements.

### Fixed

- Fixed DeepVariant and DeepSomatic WES calling with normalized BED input.
- Fixed PAVE temporary-directory creation for SAGE tumor-only processing.
- Fixed BWA-MEM2 execution on AMD Zen 4 processors by selecting the compatible AVX2 binary.
- Disabled Snakemake HTML report generation after an upstream Vega dependency change caused report creation to fail.

### Removed

- Removed Strelka from the active somatic variant caller ensemble.

## [0.9.0] - 2024-06-28

### Added

- Added an optional GATK4 germline short-variant workflow enabled with `--gatk-germline`.
- Added GATK HaplotypeCaller processing scattered across genomic regions.
- Added cohort joint genotyping with GenomicsDBImport and GenotypeGVCFs.
- Added SNP and indel Variant Quality Score Recalibration.
- Added population-based genotype refinement and final cohort VCF generation.
- Added support and documentation for the UGE execution backend.

### Changed

- DeepVariant remains the default germline caller; enabling `--gatk-germline` runs the GATK4 best-practices workflow in addition to DeepVariant.
- Enforced use of Snakemake versions earlier than 8.0 because of incompatible upstream changes.
- Expanded germline workflow resource configuration for SLURM and UGE.

## [0.8.0] - 2024-05-07

### Added

- Added UGE/SGE cluster execution through the `--mode uge` option.
- Added scheduler-specific resource configurations for local, SLURM, and UGE execution.
- Added support for a custom reference resource bundle through `--resource-bundle`.
- Added automatic rewriting of bundled reference paths when a custom resource bundle is supplied.
- Added generation of a matching FastQ Screen configuration for custom resource bundles.
- Added reporting of detected Snakemake and Singularity versions.

### Changed

- Split DeepVariant into separate `make_examples`, `call_variants`, and `postprocess_variants` jobs for more efficient resource use.
- Moved the original single-step DeepVariant rule to the deprecated workflow rules.
- Updated memory calculations for UGE/SGE, where memory is allocated per CPU.
- Reorganized cluster configurations by execution backend.
- Improved Singularity cache handling outside scheduled jobs.

## [0.7.1] - 2024-02-06

### Fixed

- Moved FastQC intermediate processing to local scratch or the configured temporary directory to avoid locking and performance issues on GPFS filesystems.
- Added local scratch allocations for raw-read and aligned-BAM FastQC jobs.

## [0.7.0] - 2024-02-02

### Added

- Added whole-exome sequencing support through `--wes-mode`.
- Added support for custom exome capture regions through `--wes-bed`.
- Added preparation, padding, sorting, compression, and indexing of exome capture BED files.
- Added WES-aware DeepVariant, GLnexus, GATK, Qualimap, Manta, MuSE, Strelka, Octopus, and Mutect2 processing.
- Added CNVkit copy-number analysis for paired tumor-normal WES samples.
- Added Sequenza purity, ploidy, and copy-number analysis for paired tumor-normal WES samples.
- Added WES pipeline documentation and CI dry-run coverage.

### Changed

- Restricted WES variant calling and downstream processing to configured capture regions where supported.
- Selected WES-specific DeepVariant and GLnexus models when running in exome mode.
- Disabled WGS-specific HMF Tools CNV processing for WES datasets.
- Updated the project description from a WGS pipeline to a combined WGS and WES pipeline.

### Fixed

- Fixed parsing of SRA-style read identifiers when extracting flowcell and lane information.
- Ensured normalized VCFs are compressed and indexed before interval-based GATK selection.

## [0.6.0] - 2023-12-06

_First public release of genome-seek._

### Added

- Added comprehensive clinical whole-genome sequencing workflows for germline and somatic small variants, structural variants, copy-number variants, HLA typing, variant annotation, and quality control.
- Added filtered Manta somatic structural-variant output using read support, allele fraction, coverage, and contig realignment.
- Added Purple-derived somatic MAF files and cohort-level maftools summaries and plots.
- Added population-frequency filtering of variants supplied to Purple.

### Changed

- Reworked somatic VCF normalization and merging across callers.
- Replaced caller intersection with prioritized caller-aware variant merging.
- Updated Purple to consume the merged somatic callset rather than separate caller-specific VCFs.
- Updated Strelka output to include allelic-depth annotations required by downstream processing.
- Standardized tumor-only Manta output naming with tumor-normal output naming.
- Increased and tuned resource allocations for somatic calling, annotation, and CNV analysis.

### Fixed

- Sorted scattered Octopus VCF chunks before concatenation.
- Improved filtering of Octopus calls to retain passing somatic variants.
- Removed unnecessary splitting of normalized somatic VCFs into separate tumor and normal files.

## [0.1.0-beta] - 2022-02-18

### Added

- Added the initial recommended scaffold for building a Snakemake pipeline.

[0.9.4]: https://github.com/OpenOmics/genome-seek/compare/v0.9.3...v0.9.4
[0.9.3]: https://github.com/OpenOmics/genome-seek/compare/v0.9.2...v0.9.3
[0.9.2]: https://github.com/OpenOmics/genome-seek/compare/v0.9.1...v0.9.2
[0.9.1]: https://github.com/OpenOmics/genome-seek/compare/v0.9.0...v0.9.1
[0.9.0]: https://github.com/OpenOmics/genome-seek/compare/v0.8.0...v0.9.0
[0.8.0]: https://github.com/OpenOmics/genome-seek/compare/v0.7.1...v0.8.0
[0.7.1]: https://github.com/OpenOmics/genome-seek/compare/v0.7.0...v0.7.1
[0.7.0]: https://github.com/OpenOmics/genome-seek/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/OpenOmics/genome-seek/releases/tag/v0.6.0
[0.1.0-beta]: https://github.com/OpenOmics/genome-seek/commits/v0.1.0-beta
