process parseBestDemulti {
    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rparseBestDemulti.${chunkID}" }

    publishDir "${params.outDir}/demultiplex_assignments.per_chunk", 
        mode: 'copy', overwrite: 'true', 
        pattern:"07-parse_best_and_demultiplex_per_chunk/*.json",
        saveAs: { it.replace("07-parse_best_and_demultiplex_per_chunk/", "") }

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