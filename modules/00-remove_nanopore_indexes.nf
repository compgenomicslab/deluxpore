process removeNanoporeIndexes {
    conda "${params.general_env}"
    
    label 'medium'

    tag { "${params.projectName}.rRemoveNanoIndexes.${chunkID}" }

    publishDir "${params.outDir}/00-trimm_nanopore_indexes_per_chunk", mode: 'copy', overwrite: 'true'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFile)

    output:
    tuple val(chunkID), path("${chunkID}.nano_trimmed.fastq.gz")

    script:
    """
    mkdir -p "00-trimm_nanopore_indexes_per_chunk"

    porechop \
    -i ${readFile} \
    -o "${chunkID}.nano_trimmed.fastq.gz" \
    --verbosity 2 
    """
}