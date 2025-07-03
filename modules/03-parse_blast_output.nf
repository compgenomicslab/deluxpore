process parseBlastOut {    
    label 'fast'

    conda "${params.general_env}"
    
    tag { "${params.projectName}.rParseBlastOut.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFileFasta), path(blastOut)

    output:
    tuple val(chunkID), path ("03-parse_reads2indexes/${chunkID}.complete_index.fna")

    script:
    """

    mkdir -p "03-parse_reads2indexes"
    
    03-parse_blastn_output.py \
    -i "${blastOut}" \
    -r "${readFileFasta}" \
    -o "03-parse_reads2indexes/${chunkID}.complete_index.fna"

    """
}