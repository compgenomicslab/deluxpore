process parseBlastOut {    
    label 'medium'

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"
    
    tag { "${params.projectName}.rParseBlastOut.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false', enabled: params.publishIntermediate

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