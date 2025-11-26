nextflow.enable.dsl=2
/*
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *   Define the default parameters 
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */

params.general_env = "/home/acatalina/miniforge3/envs/metag_long"
params.bin = "scripts"


log.info """\
    ================================================
     N X F - D E M U L T I P L E X   P I P E L I N E
    ================================================
    Project Name   : ${params.projectName}
    Reads Dir      : ${params.readsDir}
    Experimental Design: ${params.experimentalDesign}
    Output Dir     : ${params.outDir}
    Config Profile : ${workflow.profile}
    Start Time     : ${workflow.start}
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

    read_ch = Channel.fromPath("${params.readsDir}/${params.readsFileExtension}", type: 'file')
    read_ch = read_ch.map { file ->
        def chunkID = file.name.replaceAll("\\.fastq\\.gz", "")
        return [chunkID, file]
    }

    // 0) Generate index files one per demultiplexing experiment
    runIndexFilesInput = Channel.fromPath("${params.experimentalDesign}", type: 'file')
    runIndexFilesOutput = generateIndexFiles(runIndexFilesInput)

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

    // 5) Parse BLAST output to extract unique query index sequences
    parseBlastOutputInput = transFastqtoFastaOutput.join(mapReads2DBOutput)
    parseBlastOutOutput = parseBlastOut(parseBlastOutputInput)

    // 6) Extract unique query indexes from query complete fasta sequence files based on mapping to subject index sequences and their fixed position
    extractUniqQueryIndexInput = transFastqtoFastaOutput.join(mapReads2DBOutput)
    extractUniqQueryIndexOutput = extractUniqQueryIndex(extractUniqQueryIndexInput)

    // 6) Calculate Levenshtein distance between read unique index sequences (and RC) and illumina unique index sequences
    calcLevDistanceInput = extractUniqQueryIndexOutput.combine(runIndexFilesOutput.map { [it] })
    calcLevDistanceOutput = calcLevDistance(calcLevDistanceInput)

    // 6.1) Remove Illumina index sequences 
    if (params.trimmIlluminaIndexes) {
        // Remove Ilumina indcex from Nanopore clean fasta sequences
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