[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

# deluxpore
**deluxpore** is a bioinformatic pipeline designed to demultiplex [Oxford Nanopore](https://nanoporetech.com/) reads that have previously been multiplexed using [Illumina Dual-Index](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-dual-indexes.html) identifiers. 


## Table of contents

* [Installation and dependencies](#install)
* [Quick usage examples](#quick-usage-examples)
* [How it works](#how-it-works)
    * [Find matching adapter sets](#find-matching-adapter-sets)
    * [Trim adapters from read ends](#trim-adapters-from-read-ends)
    * [Split reads with internal adapters](#split-reads-with-internal-adapters)
    * [Discard reads with internal adapters](#discard-reads-with-internal-adapters)
    * [Barcode demultiplexing](#barcode-demultiplexing)
    * [Barcode demultiplexing with Albacore](#barcode-demultiplexing-with-albacore)
    * [Output](#output)
    * [Verbose output](#verbose-output)
* [Full usage](#full-usage)
* [Acknowledgements](#acknowledgements)
* [License](#license)

<a name="install"></a>
## Installation and dependencies
This pipeline was built using [Nextflow](https://www.nextflow.io/). The whole workflow runs through a built-in conda environment. You should:
- Install [Nextflow version >=23.04.0](https://www.nextflow.io/docs/latest/getstarted.html#installation)
- Install [Conda]
The complete workflow can be run on a personal computer, as [Oxford Nanopore](https://nanoporetech.com/) reads are tipically reported as **fastq.gz chunk files**. The smart [Nextflow wildcard](https://www.nextflow.io/docs/latest/working-with-files.html) path matcher allows for easy paralellization. The pipelie can also be run on HPC clusters: [Nextflow](https://www.nextflow.io/) offers multiple [executor](https://www.nextflow.io/docs/latest/executor.html) options; however, this pipeline is only prepared for local and [Slurm](https://slurm.schedmd.com/documentation.html) profiles. 
