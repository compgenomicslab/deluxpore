process parseBestDemulti {    
    label 'fast'

    conda "${params.general_env}"

    tag { "${params.projectName}.rparseBestDemulti.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false', pattern:"*.json"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFileFasta), path(levDistMat), path(expDesFile)

    output:
    tuple val(chunkID), path ("06-parse_best_and_demultiplex_per_chunk/*.fna"), path ("06-parse_best_and_demultiplex_per_chunk/*.json")

    script:
    """

    mkdir -p "06-parse_best_and_demultiplex_per_chunk"
    
    06-parse_best_and_demultiplex.py \
    -i  ${readFileFasta} \
    -id ${levDistMat} \
    -ed ${expDesFile} \
    -o "06-parse_best_and_demultiplex_per_chunk"

    """
}