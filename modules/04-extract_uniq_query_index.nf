process extractUniqQueryIndex {    
    label 'fast'

    conda "${params.general_env}"

    tag { "${params.projectName}.rExtractUniqueQuery.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFileFasta), path(blastOut)

    output:
    tuple val(chunkID), path ("04-extract_uniq_query_indexes/${chunkID}.unique_query_indexes.fna")

    script:
    """

    mkdir -p "04-extract_uniq_query_indexes"
    
    04-extract_unique_query_indexes.py \
    -i "${blastOut}" \
    -r "${readFileFasta}" \
    -o "04-extract_uniq_query_indexes/${chunkID}.unique_query_indexes.fna"

    """
}

