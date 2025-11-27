process parseBestDemulti {
    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rparseBestDemulti.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'true', pattern:"07-parse_best_and_demultiplex_per_chunk/*.json"

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFileFasta), path(levDistMat), path(expDesFile)

    output:
    tuple val(chunkID), path ("07-parse_best_and_demultiplex_per_chunk/*.fna"), path ("07-parse_best_and_demultiplex_per_chunk/*.json")

    script:
    """

    mkdir -p "07-parse_best_and_demultiplex_per_chunk"
    
    07-parse_best_and_demultiplex.py \
    -i  ${readFileFasta} \
    -id ${levDistMat} \
    -ed ${expDesFile} \
    -o "07-parse_best_and_demultiplex_per_chunk"

    """
}