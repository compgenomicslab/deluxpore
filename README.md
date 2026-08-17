[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

# deluxpore
**deluxpore** is a bioinformatic pipeline designed to demultiplex [Oxford Nanopore](https://nanoporetech.com/) reads that have previously been multiplexed using [Illumina Dual-Index](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-dual-indexes.html) identifiers. 


## Table of contents

* [Installation and Dependencies](#install)
* [Simple Usage](#simple-usage)
* [Full Usage](#full-usage)
* [Custom Index Sequences](#custom-indexes)
* [Post-demultiplexing Trimming and Chimera Removal](#trimming)
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
> SLURM queue names defined in the configuration files (e.g., 'compute', 'bigmem') are specific to our institutional HPC system. Users should modify these values in the SLURM configuration to match their local cluster queue names.

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
nextflow run -latest ./deluxpore/main.nf --help

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
                           --customCompleteIndexes  Path to complete index sequences FASTA (adapter + barcode)
                           --customUniqueIndexes    Path to unique barcode-only index sequences FASTA

Optional parameters:
  --trimandfilterNanopore  Enable Nanopore read trimming/filtering [default: true]
  --nanoQscore             Minimum quality score [default: 20]
  --nanoLength             Minimum read length [default: 100]
  --trimmIlluminaIndexes   Trim Illumina adapter sequences from demultiplexed reads.
                           Trimming runs per-sample after demultiplexing; also enables
                           chimera detection (see --removeChimeras) [default: false]
  --removeChimeras         Split reads at confident internal adapter occurrences
                           (chimera detection) instead of leaving them uncorrected [default: false]
  --removeChimerasCoverage Minimum fraction of the adapter template an internal alignment
                           must cover to be treated as a genuine chimeric junction rather
                           than noise; 0.7 sits in the valley between coincidental short
                           matches (~0.15-0.20) and genuine chimeras (~0.90-1.0) [default: 0.7]

Resource limits:
  --max_cpus             Maximum CPUs to use [default: auto-detected]
  --max_memory           Maximum memory to use [default: 16 GB]

Other:
  --conda_env                Path to pre-built conda environment [default: null]
  --publishIntermediate      Publish intermediate files [default: false]
  --rcCollisionWithholdDist  RC collision exclusion threshold. Reads whose best
                             RC collision Levenshtein distance is <= this value are
                             excluded from their assigned sample FASTA and reported
                             as rc_collision_excluded in ambiguous_reads.tsv.
                             0 (default) excludes only exact RC mirrors (dist=0) --
                             an 8bp barcode matching a different, wrong-slot index
                             with zero edits is not plausible ONT noise and is
                             treated as a genuine collision. Set to -1 to disable
                             exclusion entirely (flag in TSV only); set to 1 to
                             also withhold near-exact (dist<=1) collisions.
                             [default: 0]

  --version              Show pipeline version
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

<a name="trimming"></a>
## Post-demultiplexing Trimming and Chimera Removal

When `--trimmIlluminaIndexes` is enabled, deluxpore runs a per-sample trimming and (optionally) chimera-splitting step **after** demultiplexing. This design is intentional: barcode assignment uses the full, untrimmed reads (giving BLAST the most signal), and trimming is applied only once reads are correctly sorted into per-sample files.

### What gets trimmed

For each demultiplexed FASTA, the pipeline runs `blastn` against the same complete-adapter BLAST database used during demultiplexing, then applies terminal trimming:

- **Start-terminal** hits (adapter alignment ending within the first 73 bp) — the adapter and everything before it are removed.
- **End-terminal** hits (adapter alignment starting within the last 73 bp) — the adapter and everything after it are removed.

### Chimera detection (`--removeChimeras`)

Some ONT reads contain internal Illumina adapter sequences from accidental ligation during library preparation (chimeric reads). With `--removeChimeras true`, internal BLAST hits whose alignment covers at least `--removeChimerasCoverage` of the adapter template are treated as genuine chimeric junctions: the read is split at each junction into separate fragments.

Fragment IDs get a `_frag1`, `_frag2`, … suffix. A merged chimera report across all samples is written to:
```
{outDir}/ambiguous_reads_report/chimera_reads.tsv
```

| Column | Description |
|--------|-------------|
| `original_read_id` | Read ID before splitting |
| `num_fragments` | Number of fragments produced (1 = detected but not split) |
| `fragment_id` | Fragment ID (`_frag1`, `_frag2`, …), or same as `original_read_id` when not split |
| `fragment_start` | Start position in the original read (0-based) |
| `fragment_end` | End position in the original read (0-based exclusive) |
| `junction_coverage` | Fraction of the adapter template covered by the best chimeric junction hit |
| `chimera_split` | `yes` if the read was split; `no` if detected but `--removeChimeras` was false |

> [!NOTE]
> Setting `--removeChimerasCoverage` too low risks splitting reads at coincidental short adapter matches. The default of 0.7 was chosen based on the bimodal distribution of internal alignment coverage in ONT UCE libraries, where noise clusters below ~0.2 and genuine chimeras cluster above ~0.9.

<a name="ambiguous-reads"></a>
## Ambiguous Read Assignments

During demultiplexing, some reads cannot be unambiguously assigned to a sample, or are assigned correctly but flagged for inspection. Six situations are reported:

- **`tie_both_valid`** — A read's detected barcodes match two different valid sample combinations with equal edit distance. Unresolvable; always excluded.
- **`single_barcode_multi_sample`** — Only one barcode (i5 or i7) was detected in the read, but that barcode is shared by more than one sample in the experimental design. Unresolvable; always excluded.
- **`invalid_index_pair`** — Both i5 and i7 were detected confidently, but the combination doesn't match any sample in your experimental design (e.g. index hopping, a chimeric read). Always excluded.
- **`no_barcode_match`** — bin/04 extracted *something* from this read, but nothing came within `MAX_BARCODE_MATCH_DIST` of any catalog barcode in either slot (see [Barcode matching threshold](#barcode-threshold) below). `barcode_info` notes the closest distance actually found, so you can tell this apart from a read where no adapter was detected at all. Always excluded, unless rescued (see [Dual-confirmation rescue](#rescue) below).
- **`rescued_barcode_match`** — A `no_barcode_match` read whose two individually-too-loose candidates (one per slot) turned out to agree on a real, valid sample pair. Always included — see [Dual-confirmation rescue](#rescue) below.
- **`rc_collision`** — The barcode extracted from one adapter slot (i5 or i7) is a near-exact reverse complement of a barcode used by a different sample. See [RC collision handling](#rc-collision) below — whether it's kept or excluded depends on whether a second barcode corroborates the assignment.

After each run, deluxpore writes a merged report to:
```
{outDir}/ambiguous_reads_report/ambiguous_reads.tsv
```

The TSV has five columns:

| Column | Description |
|--------|-------------|
| `read_id` | Nanopore read identifier |
| `ambiguity_type` | `tie_both_valid`, `single_barcode_multi_sample`, `invalid_index_pair`, `no_barcode_match`, `rescued_barcode_match`, or `rc_collision` |
| `barcode_info` | The slot extracted and the colliding index with its Levenshtein distance (for `rc_collision`); the mismatched i5+i7 combination (for `invalid_index_pair`); the closest distance found (for `no_barcode_match`); the two rescued candidates and their distances (for `rescued_barcode_match`); etc. |
| `possible_samples` | Sample the read was assigned to (or would have been assigned to, for excluded reads) |
| `decision` | `included` (kept in its sample FASTA) or `excluded` (withheld from all sample FASTAs) |

Use this report to identify which samples are affected by barcode collisions and verify whether the ambiguous reads are consistent with your plate layout.

`{outDir}/ambiguous_reads_report/index_assignment_summary.tsv` gives the same picture as one-line-per-chunk counts (`both`/`i5_only`/`i7_only`/`unassigned`/`no_barcode_match`/`total_reads`) — `total_reads` there is the *entire* population entering barcode matching (i.e. everything that survived Nanopore quality/length filtering), so every read is accounted for in one place: `total_reads = both + i5_only + i7_only + unassigned + no_barcode_match`.

<a name="barcode-threshold"></a>
### Barcode matching threshold

| Parameter | Default | Description |
|---|---|---|
| `--maxBarcodeMatchDist` | `2` | Maximum Levenshtein distance allowed between an extracted barcode and a catalog barcode for it to count as a match. Absolute edit distance, not relative to barcode length — adjust per index kit. |

<a name="rescue"></a>
### Dual-confirmation rescue

A single barcode match beyond `--maxBarcodeMatchDist` is too unreliable to trust on its own. But if a read's i5 *and* i7 **both** land on a real sample pair at that looser distance, that agreement is strong evidence the read is genuine — two independent, error-prone signals landing on the same specific real combination by chance is far less likely than either one doing so alone. This rescue exists to recover those reads instead of discarding them outright. It never invents a pair that isn't in your experimental design; it can only promote a read when both slots agree on one that's already real.

| Parameter | Default | Description |
|---|---|---|
| `--rescueBarcodeMatch` | `false` | If a read fails `--maxBarcodeMatchDist` in both slots, retry with `--rescueMaxDist`; promote it to a match only if both slots agree on a real, valid sample pair. |
| `--rescueMaxDist` | `maxBarcodeMatchDist + 1` | Distance threshold used for the rescue retry. Only applies when `--rescueBarcodeMatch true`. |

<a name="rc-collision"></a>
### RC collision handling

Some index kits contain pairs whose i5 and i7 barcodes are reverse complements of each other (e.g. sample A's i5 = RC of sample B's i7). This means a read from sample A sequenced in reverse-complement orientation presents barcodes that are indistinguishable from sample B.

An `rc_collision` is only ever reported when two conditions both hold:
1. The candidate's edit distance to the colliding, wrong-slot barcode is <= 1 (looser matches are far more likely to be ordinary ONT basecalling noise on an 8bp barcode than a genuine collision).
2. That colliding index also has its own genuine BLAST-supported extraction on the *same read* — i.e. there's real alignment evidence for it, not just an 8bp string coincidence discovered by comparing the extracted barcode against the full catalog. A read whose only "evidence" for the colliding index is a low-significance stray BLAST hit (or none at all) never gets flagged in the first place.

Once flagged, `decision` is determined as follows:

- **Dual-confirmed (always `included`)** — if the read's i5 *and* i7 independently confirm a valid, real sample pair, the collision is noted (`barcode_info` includes `dual_confirmed`) but the read is always kept. The slot itself (i5 vs i7) is established by BLAST-aligning ~60bp of adapter backbone, not the 8bp barcode, so an 8bp cross-slot coincidence doesn't undermine a match that a second, independently-extracted barcode already corroborates.
- **Single-barcode assignments** — no second barcode to corroborate, so `--rcCollisionWithholdDist` controls the outcome:
  - **`0` (default)** — `excluded` for exact RC mirrors (dist=0); `included` (flagged only) for dist=1, since those are plausibly just sequencing error.
  - **`-1`** — always `included`; flag but never exclude.
  - **`1`** — `excluded` for both dist=0 and dist=1.

<a name="acknowledgements"></a>
## Acknowledgements
The original demultiplexing approach was conceived and prototyped by [Claudia Sanchis López](https://github.com/compgenomicslab/demultiplex-ont-illumina). 
