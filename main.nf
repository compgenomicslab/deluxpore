nextflow.enable.dsl=2
/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   Define the default parameters
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

params.bin = "scripts"
params.customCompleteIndexes = null
params.customUniqueIndexes   = null

/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   Help message
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

def helpMessage() {
    log.info """
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
      --conda_env            Path to pre-built conda environment [default: null]
      --publishIntermediate  Publish intermediate files [default: false]

      --version              Show pipeline version
      --help                 Show this help message

    Examples:
      # Using NEBNext indexes
      nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs NEBNext -params-file params.json

      # Using Nextera indexes
      nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs NEXTERA -params-file params.json

      # Using custom index sequences
      nextflow run ktlina/deluxpore -profile local,conda --libraryIndexSeqs custom \\
        --customCompleteIndexes /path/to/complete_indexes.fna \\
        --customUniqueIndexes /path/to/unique_indexes.fna \\
        -params-file params.json

    """.stripIndent()
}

// Show version if --version is specified
if (params.version) {
    log.info "deluxpore version: ${workflow.manifest.version}"
    exit 0
}

// Show help message if --help is specified
if (params.help) {
    helpMessage()
    exit 0
}

log.info """\
    ================================================
           D E L U X P O R E  P I P E L I N E
    ================================================
    Project Name        : ${params.projectName}
    Reads Dir           : ${params.readsDir}
    Experimental Design : ${params.experimentalDesign}
    Output Dir          : ${params.outDir}
    Work dir            : ${workDir}
    Config Profile      : ${workflow.profile}
    Start Time          : ${workflow.start}
    """
    .stripIndent()


/*
 *  ~~~~~~~~~~~~~~~~~~
 *        Steps
 *  ~~~~~~~~~~~~~~~~~~
*/

include { generateIndexFiles } from './modules/00-generate_index_files'

include { removeNanoporeIndexes } from './modules/00-trim_and_filter'
include { filterNanoporeReads }   from './modules/00-trim_and_filter'

include { transFastqtoFasta } from './modules/01-transform_to_fasta'

include { createDB }      from './modules/02-reads2database'
include { mapReads2DB }   from './modules/02-reads2database'

include { extractUniqQueryIndex } from './modules/04-extract_uniq_query_index'

include { calcLevDistance } from './modules/05-calc_lev_distance'

include { parseBestDemulti } from './modules/07-parse_best_and_demultiplex'

include { concatenateSamples }         from './modules/08-concat_sample_fna_files'
include { concatenateSummaries }       from './modules/08-concat_sample_fna_files'
include { concatenateAmbiguousReport } from './modules/08-concat_sample_fna_files'
include { concatenateChimeraReports }  from './modules/08-concat_sample_fna_files'

include { trimAndRemoveChimeras } from './modules/09-trim_and_remove_chimeras'


/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow {

    // 0) Generate index files
    runIndexFilesInput = Channel.fromPath("${params.experimentalDesign}", type: 'file')

    if (params.libraryIndexSeqs.toLowerCase() == "custom") {
        if (!params.customCompleteIndexes || !params.customUniqueIndexes) {
            error "When libraryIndexSeqs is 'custom', you must provide --customCompleteIndexes and --customUniqueIndexes"
        }
        indexCompleteFile = file(params.customCompleteIndexes, checkIfExists: true)
        indexUniqueFile   = file(params.customUniqueIndexes,   checkIfExists: true)
    } else {
        indexCompleteFile = file("${projectDir}/assets/${params.libraryIndexSeqs}.complete_indexes.fna", checkIfExists: true)
        indexUniqueFile   = file("${projectDir}/assets/${params.libraryIndexSeqs}.unique_indexes.fna",   checkIfExists: true)
    }

    generateIndexFilesInput = runIndexFilesInput.map { expDesign ->
        tuple(expDesign, indexCompleteFile, indexUniqueFile)
    }
    runIndexFilesOutput = generateIndexFiles(generateIndexFilesInput)

    read_ch = Channel.fromPath("${params.readsDir}/${params.readsFileExtension}", type: 'file')
    read_ch = read_ch.map { file ->
        def extension = params.readsFileExtension
            .replace("*", "")
        def chunkID = file.name.replace(extension, "")
        return [chunkID, file]
    }

    if (params.trimandfilterNanopore) {
        // 1) Remove Nanopore indexes and filter by quality/length
        removeNanoporeIndexesOutput = removeNanoporeIndexes(read_ch)
        filterNanoporeReadsOutput   = filterNanoporeReads(removeNanoporeIndexesOutput)
        transFastqtoFastaOutput     = transFastqtoFasta(filterNanoporeReadsOutput)
    } else {
        transFastqtoFastaOutput = transFastqtoFasta(read_ch)
    }

    // 2) Map reads to Illumina index database
    createDBInput  = runIndexFilesOutput.map { tuple -> return tuple[0] }
    // .first() converts the single-item queue channel to a value channel so it
    // can be consumed by both mapReads2DB (step 4) and trimAndRemoveChimeras (step 9).
    createDBOutput = createDB(createDBInput).first()

    mapReads2DBInput  = transFastqtoFastaOutput.combine(createDBOutput)
    mapReads2DBOutput = mapReads2DB(mapReads2DBInput)

    // 3) Extract unique query indexes
    extractUniqQueryIndexInput  = transFastqtoFastaOutput.join(mapReads2DBOutput)
    extractUniqQueryIndexOutput = extractUniqQueryIndex(extractUniqQueryIndexInput)

    // 4) Calculate Levenshtein distance; also detects RC collision candidates
    calcLevDistanceInput  = extractUniqQueryIndexOutput.combine(runIndexFilesOutput.map { [it] })
    calcLevDistanceOutput = calcLevDistance(calcLevDistanceInput)

    // Sample assignment always runs on untrimmed reads.
    // Illumina-index trimming and chimera splitting happen after demultiplexing (step 09).
    parseBestDemultiInput = transFastqtoFastaOutput
        .join(calcLevDistanceOutput)
        .combine(runIndexFilesInput)

    // 7) Parse distance matrix, extract best distance values per read and demultiplex
    parseBestDemultiOutput = parseBestDemulti(parseBestDemultiInput)

    // 8) Concatenate per-chunk sample files into final per-sample files
    allSampleFiles = parseBestDemultiOutput
        .map { chunkID, sampleFilesList, jsonFile, tsvReport, ambiguousFastas ->
            return sampleFilesList
        }
        .collect()
        .flatten()
        .filter { file -> file.name.endsWith('.fna') }
        .map { file ->
            def sampleName = file.name.split('\\.')[0]
            return [sampleName, file]
        }
        .groupTuple()

    concatenatedSamples = concatenateSamples(allSampleFiles)

    // 9) Per-sample Illumina-index trimming and chimera splitting
    //    merge all per-sample chimera reports into one file in ambiguous_reads_report/
    if (params.trimmIlluminaIndexes || params.removeChimeras) {
        trimInput = concatenatedSamples
            .combine(createDBOutput)
            .combine(Channel.value(indexCompleteFile))
        trimOutput = trimAndRemoveChimeras(trimInput)
        concatenateChimeraReports(trimOutput.chimeraReport.collect())
    }

    // 10) Merge per-chunk ambiguous FASTA files into one file per ambiguity type
    //     rc_collision reads are reported in ambiguous_reads.tsv only — no FASTA needed
    allAmbiguousFastas = parseBestDemultiOutput
        .map { chunkID, sampleFilesList, jsonFile, tsvReport, ambiguousFastas -> ambiguousFastas }
        .collect()
        .flatten()
        .filter { file -> !file.name.startsWith('rc_collision') }
        .map { file ->
            def type = file.name.split('\\.')[0]
            return [type, file]
        }
        .groupTuple()

    concatenateSummaries(allAmbiguousFastas)

    // 10) Merge per-chunk ambiguous read TSV reports into a single report
    allTsvReports = parseBestDemultiOutput
        .map { chunkID, sampleFilesList, jsonFile, tsvReport, ambiguousFastas -> tsvReport }
        .collect()

    concatenateAmbiguousReport(allTsvReports)

}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
