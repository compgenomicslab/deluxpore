nextflow.enable.dsl=2
/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   Define the default parameters 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

params.bin = "scripts"

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

//Import required modules. Repeated modules will be loaded with corresponding aliases

include { generateIndexFiles } from './modules/00-generate_index_files'

include { removeNanoporeIndexes } from './modules/00-trim_and_filter'
include { filterNanoporeReads } from './modules/00-trim_and_filter'

include { transFastqtoFasta } from './modules/01-transform_to_fasta'

include { createDB } from './modules/02-reads2database'
include { mapReads2DB } from './modules/02-reads2database'

include { parseBlastOut } from './modules/03-parse_blast_output'

include { extractUniqQueryIndex } from './modules/04-extract_uniq_query_index'

include { calcLevDistance } from './modules/05-calc_lev_distance'

include { removeIlluminaIndexes } from './modules/06-remove_illumina_indexes'

include { parseBestDemulti } from './modules/07-parse_best_and_demultiplex'

include { concatenateSamples } from './modules/08-concat_sample_fna_files'




/*
    ~~~~~~~~~~~~~~~~~~
     Run workflow
    ~~~~~~~~~~~~~~~~~~
*/
workflow {
    
    // 0) Generate index files one per demultiplexing experiment
    runIndexFilesInput = Channel.fromPath("${params.experimentalDesign}", type: 'file')
    runIndexFilesOutput = generateIndexFiles(runIndexFilesInput)

    read_ch = Channel.fromPath("${params.readsDir}/${params.readsFileExtension}", type: 'file')
    // read_ch = read_ch.map { file ->
    //     def chunkID = file.name.replaceAll("\\.fastq\\.gz", "")
    //     return [chunkID, file]
    // }
    read_ch = read_ch.map { file ->
        def extension = params.readsFileExtension
            .replace("*", "")  // Remove glob asterisk: "*.fastq.gz" → ".fastq.gz"
        def chunkID = file.name.replace(extension, "")
    return [chunkID, file]
    }

    if (params.trimandfilterNanopore){
        // 1) Remove Nanopore indexes from read files
        removeNanoporeIndexesOutput = removeNanoporeIndexes(read_ch)
        filterNanoporeReadsOutput = filterNanoporeReads(removeNanoporeIndexesOutput)

        // 3) Transform read files to fasta format
        transFastqtoFastaOutput = transFastqtoFasta(filterNanoporeReadsOutput)
    } else {
        transFastqtoFastaOutput = transFastqtoFasta(read_ch)
    }

    // 4) Map fastas to illumina index sequences
    createDBInput = runIndexFilesOutput.map { tuple ->
        return tuple[0]}
    createDBOutput = createDB(createDBInput)

    // Map reads making sure to use fasta formated reads
    mapReads2DBInput = transFastqtoFastaOutput.combine(createDBOutput)
    mapReads2DBOutput = mapReads2DB(mapReads2DBInput)
    // mapReads2DBOutput.view()

    // // 5) Parse BLAST output to extract complete query index sequences
    // parseBlastOutputInput = transFastqtoFastaOutput.join(mapReads2DBOutput)
    // parseBlastOutOutput = parseBlastOut(parseBlastOutputInput)

    // 6) Extract unique query indexes from query complete fasta sequence files based on mapping to subject index sequences and their fixed position
    extractUniqQueryIndexInput = transFastqtoFastaOutput.join(mapReads2DBOutput)
    extractUniqQueryIndexOutput = extractUniqQueryIndex(extractUniqQueryIndexInput)

    // 6) Calculate Levenshtein distance between read unique index sequences (and RC) and illumina unique index sequences
    calcLevDistanceInput = extractUniqQueryIndexOutput.combine(runIndexFilesOutput.map { [it] })
    calcLevDistanceOutput = calcLevDistance(calcLevDistanceInput)

    if (params.trimmIlluminaIndexes) {
        // 5) Parse BLAST output - only needed for trimming
        parseBlastOutputInput = transFastqtoFastaOutput.join(mapReads2DBOutput)
        parseBlastOutOutput = parseBlastOut(parseBlastOutputInput)
        
        // Remove Illumina index from Nanopore clean fasta sequences
        removeIlluminaIndexesInput = transFastqtoFastaOutput
            .join(parseBlastOutOutput)
        removeIlluminaIndexesOutput = removeIlluminaIndexes(removeIlluminaIndexesInput)

        // Combine output to distance and index sequences
        parseBestDemultiInput = removeIlluminaIndexesOutput
            .join(calcLevDistanceOutput)
            .combine(runIndexFilesInput)
    } else {
        // Instead use complete Nanopore fasta sequences
        parseBestDemultiInput = transFastqtoFastaOutput
            .join(calcLevDistanceOutput)
            .combine(runIndexFilesInput)
    }

    // 7) Parse distance matrix, extract best distance values per read and demultiplex
    parseBestDemultiOutput = parseBestDemulti(parseBestDemultiInput)

    // 8) Join chunk files by sample name and concatenate into final sample files
    allSampleFiles = parseBestDemultiOutput
        .map { chunkID, sampleFilesList, jsonFile ->
            return sampleFilesList
        }
        .collect()  // Wait for all chunks to complete
        .flatten()  // Flatten all files into single channel
        .filter { file ->
            file.name.endsWith('.fna')
        }
        .map { file ->
            def sampleName = file.name.split('\\.')[0]
            return [sampleName, file]
        }
        .groupTuple()

    concatenatedSamples = concatenateSamples(allSampleFiles)
    

}

workflow.onComplete {
    println "Pipeline completed at: ${workflow.complete}"
    println "Time to complete workflow execution: ${workflow.duration}"
    println "Execution status: ${workflow.success ? 'Succesful' : 'Failed' }"
}

workflow.onError {
    println "Oops... Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}
