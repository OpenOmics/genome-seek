<div align="center">

  <h1 style="font-size: 250%">genome-seek 🔬</h1>

  <b><i>Whole Genome and Exome Clinical Sequencing Pipeline</i></b><br> 
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

Welcome to genome-seek's documentation! This guide is the main source of documentation for users who are getting started with the OpenOmics [genome-seek pipeline](https://github.com/OpenOmics/genome-seek/). 

The **`./genome-seek`** pipeline is composed of several interrelated sub-commands to set up and run the pipeline across different systems. Each of the available sub-commands performs different functions: 

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">genome-seek <b>run</b></code>](usage/run.md)   
    Run the genome-seek pipeline with your input files.

!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">genome-seek <b>unlock</b></code>](usage/unlock.md)  
    Unlocks a previous runs output directory.

</section>

<section align="center" markdown="1" style="display: flex; flex-wrap: row wrap; justify-content: space-around;">


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">genome-seek <b>install</b></code>](usage/install.md)  
    Download remote reference files locally.


!!! inline custom-grid-button ""

    [<code style="font-size: 1em;">genome-seek <b>cache</b></code>](usage/cache.md)  
    Cache remote software containers locally.  

</section>

**genome-seek** is a comprehensive clinical WGS and WES pipeline that is focused on speed. Each tool in the pipeline was benchmarked and selected due to its low run times without sacrificing accuracy or precision. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a set of FastQ files and can be run locally on a compute instance or on-premise using a cluster (recommended). A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM or UGE (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/genome-seek/usage/run/) section of each available sub-command.

For more information about issues or troubleshooting a problem, please check out our [FAQ](https://openomics.github.io/genome-seek/faq/questions/) before [opening an issue on Github](https://github.com/OpenOmics/genome-seek/issues).

## Contribute

This site is a living document, created for and by members like you. genome-seek is maintained by the members of NCBR and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull requests to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/genome-seek).

## References

<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
