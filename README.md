<div align="center">
   
  <h1>genome-seek 🔬</h1>
  
  **_Whole Genome and Exome Clinical Sequencing Pipeline._**

  [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12586652.svg)](https://doi.org/10.5281/zenodo.12586652) [![GitHub release (latest SemVer including pre-releases)](https://img.shields.io/github/v/release/OpenOmics/genome-seek?color=blue&include_prereleases)](https://github.com/OpenOmics/genome-seek/releases) [![Docker Pulls](https://img.shields.io/docker/pulls/skchronicles/genome-seek)](https://hub.docker.com/repository/docker/skchronicles/genome-seek)<br>[![tests](https://github.com/OpenOmics/genome-seek/workflows/tests/badge.svg)](https://github.com/OpenOmics/genome-seek/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/genome-seek/workflows/docs/badge.svg)](https://github.com/OpenOmics/genome-seek/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/genome-seek?color=brightgreen)](https://github.com/OpenOmics/genome-seek/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/genome-seek)](https://github.com/OpenOmics/genome-seek/blob/main/LICENSE) 
  
  <i>
    This is the home of the pipeline, genome-seek. Its long-term goals: to accurately call germline and somatic variants, to infer SVs & CNVs, and to boldly annotate variants like no pipeline before!
  </i>
</div>

## Overview

Welcome to genome-seek! Before getting started, we highly recommend reading through [genome-seek's documentation](https://openomics.github.io/genome-seek/).

The **`./genome-seek`** pipeline is composed of several interrelated sub-commands to set up and run the pipeline across different systems. Each of the available sub-commands performs different functions: 

 * [<code>genome-seek <b>run</b></code>](https://openomics.github.io/genome-seek/usage/run/): Run the genome-seek pipeline with your input files.
 * [<code>genome-seek <b>unlock</b></code>](https://openomics.github.io/genome-seek/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>genome-seek <b>cache</b></code>](https://openomics.github.io/genome-seek/usage/cache/): Cache software containers locally.

**genome-seek** is a comprehensive clinical WGS and WES pipeline that is focused on speed. Each tool in the pipeline was benchmarked and selected due to its low run times without sacrificing accuracy or precision. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster (recommended). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM or UGE (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/genome-seek/usage/run/) section of each available sub-command.

For more information about issues or troubleshooting a problem, please check out our [FAQ](https://openomics.github.io/genome-seek/faq/questions/) before [opening an issue on Github](https://github.com/OpenOmics/genome-seek/issues).

## Dependencies

**Requires:** `singularity>=3.5`  `snakemake<=7.32.3`

At the current moment, the pipeline only has two dependencies: snakemake and singularity/apptainer. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and ([singularity](https://singularity.lbl.gov/all-releases) or [apptainer](https://apptainer.org/docs/admin/main/installation.html)) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. The second dependency, i.e singularity/apptainer, handles downloading/distribution of the pipeline's software dependencies. To guarantee the highest level of reproducibility, each step of the pipeline relies on versioned images from [DockerHub](https://hub.docker.com/repositories/skchronicles). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity/apptainer will be the only two dependencies in the future.

## Installation

Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/genome-seek.git
# Change your working directory
cd genome-seek/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./genome-seek -h
```

## Contribute

This site is a living document, created for and by members like you. genome-seek is maintained by the members of NCBR and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull requests to our [GitHub repository](https://github.com/OpenOmics/genome-seek).

## Cite

If you use this software, please cite it as below:  

<details>
  <summary><b><i>@BibText</i></b></summary>
 
```text
@software{Kuhn_OpenOmics_genome-seek_2024,
  author       = {Kuhn, Skyler and
                  Lack, Justin and
                  Talsania, Keyur},
  title        = {OpenOmics/genome-seek: v0.9.0},
  month        = jun,
  year         = 2024,
  publisher    = {Zenodo},
  version      = {v0.9.0},
  doi          = {10.5281/zenodo.12586652},
  url          = {https://doi.org/10.5281/zenodo.12586652}
}
```

</details>

<details>
  <summary><b><i>@APA</i></b></summary>

```text
Kuhn, S., Lack, J., & Talsania, K. (2024). OpenOmics/genome-seek: v0.9.0 (v0.9.0). Zenodo. https://doi.org/10.5281/zenodo.12586652
```

</details>

For more citation style options, please visit the pipeline's [Zenodo page](https://doi.org/10.5281/zenodo.12586652).

## References

<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
