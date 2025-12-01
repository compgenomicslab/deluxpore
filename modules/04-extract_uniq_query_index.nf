process extractUniqQueryIndex {    
    label 'fast'

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rExtractUniqueQuery.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false'
    
    input:
    tuple val(chunkID), path(readFileFasta), path(blastOut)

    output:
    tuple val(chunkID), path ("04-extract_uniq_query_indexes/${chunkID}.unique_query_indexes.fna")

    script:
    """

    mkdir -p "04-extract_uniq_query_indexes"
    
    04-extract_unique_query_indexes.py \
    -i "${blastOut}" \
    -ik "${params.libraryIndexSeqs}" \
    -r "${readFileFasta}" \
    -o "04-extract_uniq_query_indexes/${chunkID}.unique_query_indexes.fna"

    """
}

