[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

# deluxpore
**deluxpore** is a bioinformatic pipeline designed to demultiplex [Oxford Nanopore](https://nanoporetech.com/) reads that have previously been multiplexed using [Illumina Dual-Index](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-dual-indexes.html) identifiers. 


## Table of contents

* [Installation and dependencies](#install)
* [Quick usage](#quick-usage)
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
- Install [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html)

The complete workflow can be run on a personal computer, as [Oxford Nanopore](https://nanoporetech.com/) reads are tipically reported as **fastq.gz chunk files**. This combined with the smart [Nextflow wildcard](https://www.nextflow.io/docs/latest/working-with-files.html) path matcher, allows for easy paralellization. The pipelie can also be run on HPC clusters: [Nextflow](https://www.nextflow.io/) offers multiple [executor](https://www.nextflow.io/docs/latest/executor.html) options; however, this pipeline is only prepared for local and [Slurm](https://slurm.schedmd.com/documentation.html) profiles. 

<a name="quick-usage"></a>
## Quick usage
Simple steps to run deluxpore:
  1. Create `experimental_design.csv` file, where first row is the final desired sample ID, second row is the i5 Illumina index, and third row is the i7 Illumina index:
     
     ```csv
     sample1,i501,i701
     sample2,i501,i702
     sample3,i501,i703
     ```
  2. Create `params_file.my_experiment.json`

Running deluxpore in **local profile**:
```bash
# Run BAT on each input genome, saving all results to the same folder
CAT bins -b ${genome_name}.fna -d ${path_to_CAT_database} -t ${path_to_CAT_tax_folder} -o BAT_results/${genome_name}

# Optional: to check what taxa were assigned, you can add names to them
CAT add_names -i BAT_results/${genome_name}.bin2classification.txt -o BAT_results/${genome_name}.name.txt -t ${path_to_CAT_tax_folder}
```
