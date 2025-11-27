process removeNanoporeIndexes {
    conda params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    label 'filter_trim'

    tag { "${params.projectName}.rRemoveNanoIndexes.${chunkID}" }

    publishDir "${params.outDir}/00-trim_and_filter_nanopore", mode: 'copy', overwrite: 'true'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFile)

    output:
    tuple val(chunkID), path("${chunkID}.nano_trimmed.fastq.gz")

    script:
    """
    mkdir -p "00-trim_and_filter_nanopore"

    porechop \
    -i ${readFile} \
    -o "${chunkID}.nano_trimmed.fastq.gz" \
    --verbosity 2 -t ${task.cpus}
    """
}

process filterNanoporeReads {
    conda params.conda_env ?: "${projectDir}/envs/deluxpore.yml"
    
    label 'filter_trim'

    tag { "${params.projectName}.rFilterNanoReads.${chunkID}" }

    publishDir "${params.outDir}/00-trim_and_filter_nanopore", mode: 'copy', overwrite: 'true'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(trimmedFile)

    output:
    tuple val(chunkID), path("${chunkID}.nano_trimmed_filtered.fastq.gz")

    script:
    """
    mkdir -p "00-trim_and_filter_nanopore"

    chopper -t ${task.cpus} -q ${params.nanoQscore} -l ${params.nanoLength} -i ${trimmedFile} | gzip > "${chunkID}.nano_trimmed_filtered.fastq.gz"

    """
}