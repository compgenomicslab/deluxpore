[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

# deluxpore
**deluxpore** is a bioinformatic pipeline designed to demultiplex [Oxford Nanopore](https://nanoporetech.com/) reads that have previously been multiplexed using [Illumina Dual-Index](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-dual-indexes.html) identifiers. 


## Table of contents

* [Installation and Dependencies](#install)
* [Simple Usage](#simple-usage)
* [Full Usage](#full-usage)
* [Custom Index Sequences](#custom-indexes)
* [Ambiguous Read Assignments](#ambiguous-reads)
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
                         Accepted values: NEBNext, NEXTERA, custom
                         When set to 'custom', also provide:
                           --customCompleteIndexes  Path to complete index sequences FASTA (adapter + barcode). (See assets folder for reference formatting). 
                           --customUniqueIndexes    Path to unique barcode-only index sequences FASTA

Optional parameters:
  --trimandfilterNanopore  Enable Nanopore read trimming/filtering [default: false]
  --nanoQscore             Minimum quality score [default: 20]
  --nanoLength             Minimum read length [default: 100]
  --trimmIlluminaIndexes   Trim Illumina index sequences [default: false]

Resource limits:
  --max_cpus             Maximum CPUs to use [default: auto-detected]
  --max_memory           Maximum memory to use [default: 16 GB]

Other:
  --conda_env            Path to pre-built conda environment [default: null]
  --help                 Show this help message

Examples:
  # Using NEBNext indexes
  nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs NEBNext -params-file params.json

  # Using Nextera indexes
  nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs NEXTERA -params-file params.json

  # Using custom index sequences
  nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs custom \
    --customCompleteIndexes /path/to/complete_indexes.fna \
    --customUniqueIndexes /path/to/unique_indexes.fna \
    -params-file params.json
```

<a name="custom-indexes"></a>
## Custom Index Sequences

If you used an index kit other than NEBNext or NEXTERA, you can provide your own index sequence files by setting `libraryIndexSeqs` to `custom` and supplying two FASTA files:

| Parameter | Description |
|-----------|-------------|
| `customCompleteIndexes` | FASTA file containing the **complete** index sequences, i.e. the full adapter + barcode sequence used to build the BLAST mapping database |
| `customUniqueIndexes` | FASTA file containing the **unique barcode-only** sequences (8 bp UDI barcodes) used for Levenshtein distance matching |

The sequence IDs in these files must match the index names used in the `experimentalDesign` TSV (e.g. `i501`, `i701`). Reverse complement sequences are computed automatically from `customUniqueIndexes` — you do not need to provide them separately.

> [!NOTE]
> Built-in kits (NEBNext, NEXTERA) include pre-built `*_rc.fna` files in `assets/`, but these are not read by the pipeline — reverse complements are always derived on the fly from the unique index sequences. The same applies to custom kits.

Example params file for custom indexes: `examples/params_file_custom_indexes.json`

<a name="ambiguous-reads"></a>
## Ambiguous Read Assignments

During demultiplexing, some reads cannot be unambiguously assigned to a sample. This can happen when:

- **`tie_both_valid`** — A read's detected barcodes match two different valid sample combinations with equal edit distance. The read is excluded from all sample files to avoid misassignment.
- **`single_barcode_multi_sample`** — Only one barcode (i5 or i7) was detected in the read, but that barcode is shared by more than one sample in the experimental design.

After each run, deluxpore writes a per-chunk report to:
```
{outDir}/ambiguous_reads_report/ambiguous_reads.{chunkID}.tsv
```

The TSV has four columns:

| Column | Description |
|--------|-------------|
| `read_id` | Nanopore read identifier |
| `ambiguity_type` | `tie_both_valid` or `single_barcode_multi_sample` |
| `barcode_info` | The barcode combination(s) involved |
| `possible_samples` | Pipe-separated list of sample names the read could belong to |

Use this report to identify which samples are affected by barcode collisions and verify whether the ambiguous reads are consistent with your plate layout.

<a name="acknowledgements"></a>
## Acknowledgements
The original demultiplexing approach was conceived and prototyped by [Claudia Sanchis López](https://github.com/compgenomicslab/demultiplex-ont-illumina). 
