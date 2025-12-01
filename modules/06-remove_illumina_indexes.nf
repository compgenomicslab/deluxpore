process removeIlluminaIndexes {    
    label 'medium'

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rRemoveIllIndexes.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'true'


    input:
    tuple val(chunkID), path(readFileFasta), path(completeIndexSeqs)

    output:
    tuple val(chunkID), path ("06-trimm_illumina_indexes_per_chunk/*.fna")

    script:
    """

    mkdir -p "06-trimm_illumina_indexes_per_chunk"
    
    06-remove_illumina_indexes.py \
    -i  ${readFileFasta} \
    -ir ${completeIndexSeqs} \
    -o "06-trimm_illumina_indexes_per_chunk"

    """
}