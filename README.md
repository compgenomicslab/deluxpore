[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

# deluxpore
**deluxpore** is a bioinformatic pipeline designed to demultiplex [Oxford Nanopore](https://nanoporetech.com/) reads that have previously been multiplexed using [Illumina Dual-Index](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-dual-indexes.html) identifiers. 


## Table of contents

* [Installation and Dependencies](#install)
* [Simple Usage](#simple-usage)
* [Full Usage](#full-usage)
* [Acknowledgements](#acknowledgements)

<a name="install"></a>
## Installation and Dependencies
This pipeline was built using [Nextflow](https://www.nextflow.io/). The whole workflow runs through a built-in conda environment. You should:
- Install [Conda](https://docs.conda.io/projects/conda/en/stable/user-guide/getting-started.html) or [Mamba](https://mamba.readthedocs.io/) (recommended)
- Install [Nextflow version >=23.04.0](https://www.nextflow.io/docs/latest/getstarted.html#installation)

No further installation is required — Nextflow downloads the pipeline automatically on first run. 

### Conda Environment Options

**Option 1: Let Nextflow handle the environment (simplest)**
```bash
nextflow run ktlina/deluxpore -profile local,conda -params-file params.json
```

The conda environment is built on the first run and cached for future use.

**Option 2: Pre-build the environment (faster)**
```bash
# Create environment once
mamba env create -f https://raw.githubusercontent.com/ktlina/deluxpore/main/envs/deluxpore.yml -n deluxpore

# Run with pre-built environment
nextflow run ktlina/deluxpore -profile local,conda --conda_env /path/to/envs/deluxpore -params-file params.json
```

### Execution Profiles 
[Oxford Nanopore](https://nanoporetech.com/) sequencing runs natively output reads as multiple fastq.gz chunk files, enabling parallelization without additional preprocessing. Deluxpore leverages this structure alongside Nextflow's wildcard path matching to distribute demultiplexing across available computational resources, scaling efficiently from personal workstations to HPC clusters. The pipeline can also be run on HPC clusters: [Nextflow](https://www.nextflow.io/) offers multiple [executor](https://www.nextflow.io/docs/latest/executor.html) options; however, this pipeline is only prepared for `local` and `slurm` profiles. 

| Profile | Description |
|---------|-------------|
| `local` | Run on local machine |
| `slurm` | Run on HPC cluster with SLURM scheduler |
| `conda` | Enable conda environment management |

> [!NOTE]
> SLURM queue names defined in the configuration files (e.g., 'fast', 'medium') are specific to our institutional HPC system. Users should modify these values in the SLURM configuration to match their local cluster queue names.

<a name="simple-usage"></a>
## Simple Usage
Simple steps to run deluxpore:
  1. Create `experimental_design.tsv` file, where first row is the final desired sample ID, second row is the i5 Illumina index, and third row is the i7 Illumina index:
     
     ```tsv
     sample1  i501  i701
     sample2  i502  i702
     sample3  i503  i703
     ```
     > [!NOTE]
     > Note on sample names: Although all characters are allowed in sample names, it is encouraged to avoid special characters particularly ".", to prevent processing issues.
     
   2. Copy the example params file and edit with your paths:
      ```bash
      cp examples/params_file.json my_params.json
      ```
3. Run:
   ```bash
   nextflow run ktlina/deluxpore -profile local,conda -params-file my_params.json
   ```
Code to run deluxpore on test dataset:
   ```bash
   git clone https://github.com/ktlina/deluxpore.git
   cd deluxpore

   nextflow run main.nf -profile local,conda -params-file test/test_params.json
   ```

<a name="quick-usage"></a>
## Full Usage
```angular2html
nextflow run -latest ktlina/deluxpore/main.nf --help

    =========================================
     D E L U X P O R E   P I P E L I N E
    =========================================
    
    Usage:
      nextflow run ktlina/deluxpore -profile local,conda -params-file params.json
    
    Required parameters:
      --projectName          Name for your project
      --readsDir             Path to directory containing Nanopore reads
      --readsFileExtension   Glob pattern to match input read files (e.g., *.fastq.gz, *.fq.gz, *.fastq, *.fq)
                             Each matched file is processed as a separate chunk in parallel.
                             Examples:
                                "*.fastq.gz"     - Process all gzipped fastq files as separate chunks
                                "sample1.fq.gz"  - Process a single file
                                "batch_*.fq"     - Process all files matching the pattern
      --experimentalDesign   Path to sample-to-index mapping file (TSV)
      --outDir               Output directory
      --libraryIndexSeqs     Illumina index kit used for multiplexing
                             Accepted values: NEBNext, NEXTERA

    Optional parameters:
      --trimandfilterNanopore  Enable Nanopore read trimming/filtering [default: true]
      --nanoQscore             Minimum quality score [default: 20]
      --nanoLength             Minimum read length [default: 100]
      --trimmIlluminaIndexes   Trim Illumina index sequences [default: false]

    Resource limits:
      --max_cpus             Maximum CPUs to use [default: auto-detected]
      --max_memory           Maximum memory to use [default: 16 GB]

    Other:
      --conda_env            Path to pre-built conda environment [default: null]
      --publishIntermediate  Publish intermediate files [default: false]
      
      --version              Show pipeline version    
      --help                 Show this help message

    Examples:
      # Using NEBNext indexes
      nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs nebnext -params-file params.json

      # Using Nextera indexes
      nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs nextera -params-file params.json

```

<a name="acknowledgements"></a>
## Acknowledgements
The original demultiplexing approach was conceived and prototyped by [Claudia Sanchis López](https://github.com/compgenomicslab/demultiplex-ont-illumina). 
