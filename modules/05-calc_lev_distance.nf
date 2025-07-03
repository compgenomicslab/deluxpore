process calcLevDistance {    
    label 'not_so_fast'

    conda "${params.general_env}"

    tag { "${params.projectName}.rCalcLevDistance.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'true'

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readUniqIndex), path(projectIndexFiles)

    output:
    tuple val(chunkID), path ("05-calc_lev_distance/${chunkID}.distance_matrix.tsv")

    script:
    """

    mkdir -p "05-calc_lev_distance"
    
    05-calculate_uniq_index_text_similarity.py \
    -ru "${readUniqIndex}" \
    -iu "${projectIndexFiles[1]}" \
    -iur "${projectIndexFiles[2]}" \
    -o "05-calc_lev_distance/${chunkID}.distance_matrix.tsv"

    """
}
