<div align="center">

  <h1 style="font-size: 250%">genome-seek 🔬</h1>

  <b><i>Whole Genome Clinical Sequencing Pipeline</i></b><br> 
  <a href="https://github.com/OpenOmics/genome-seek/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/genome-seek/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/genome-seek/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/genome-seek/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/genome-seek/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/genome-seek?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/genome-seek/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/genome-seek">
  </a>

  <p>
    This is the home of the pipeline, genome-seek. Its long-term goals: to accurately call germline and somatic variants, to infer SVs & CNVs, and to boldly annotate variants like no pipeline before!
  </p>

</div>  


## Overview
Welcome to genome-seek's documentation! This guide is the main source of documentation for users that are getting started with the [whole genome sequencing pipeline](https://github.com/OpenOmics/genome-seek/). 

The **`./genome-seek`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions: 

 * [<code>genome-seek <b>run</b></code>](usage/run.md): Run the genome-seek pipeline with your input files.
 * [<code>genome-seek <b>unlock</b></code>](usage/unlock.md): Unlocks a previous runs output directory.
 * [<code>genome-seek <b>cache</b></code>](usage/cache.md): Cache software containers locally.

**genome-seek** is a comprehensive clinical WGS pipeline that is focused on speed. Each tool in the pipeline was benchmarked and selected due to its low run times without sacrificing accuracy or precision. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance, on-premise using a cluster, or on the cloud (feature coming soon!). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM, or run on AWS using Tibanna (feature coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/genome-seek/issues).

## Contribute 

This site is a living document, created for and by members like you. genome-seek is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/genome-seek).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
