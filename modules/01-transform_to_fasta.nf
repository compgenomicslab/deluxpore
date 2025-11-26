process transFastqtoFasta {
    conda "${params.general_env}"
    
    label 'fast'

    tag { "${params.projectName}.rTransToFasta.${chunkID}" }

    publishDir "${params.outDir}/01-fastq2fasta", mode: 'copy', overwrite: 'false'


    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFile)

    output:
    tuple val(chunkID), path("${chunkID}.nano_trimmed_filtered.fasta")

    script:
    """
    mkdir -p "01-fastq2fasta"
    zcat ${readFile} | seqtk seq -a - > "${chunkID}.nano_trimmed_filtered.fasta"
    """
}