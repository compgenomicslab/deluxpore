process transFastqtoFasta {
    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"
    
    label 'fast'

    tag { "${params.projectName}.rTransToFasta.${chunkID}" }

    publishDir "${params.outDir}/01-fastq2fasta", mode: 'copy', overwrite: 'false'
    
    input:
    tuple val(chunkID), path(readFile)

    output:
    tuple val(chunkID), path("${chunkID}.nano_trimmed_filtered.fasta")

    script:
    """
    mkdir -p "01-fastq2fasta"
    if file "${readFile}" | grep -q "compressed"; then
        zcat ${readFile} | seqtk seq -a - > "${chunkID}.nano_trimmed_filtered.fasta"
    else
        cat ${readFile} | seqtk seq -a - > "${chunkID}.nano_trimmed_filtered.fasta"
    fi
    """
}