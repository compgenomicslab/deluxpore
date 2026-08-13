process calcLevDistance {    

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rCalcLevDistance.${chunkID}" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'true', enabled: params.publishIntermediate
    

    input:
    tuple val(chunkID), path(readUniqIndex), path(projectIndexFiles)

    output:
    tuple val(chunkID), path ("05-calc_lev_distance/${chunkID}.distance_matrix.tsv"), path ("05-calc_lev_distance/${chunkID}.rc_collision_events.tsv")

    script:
    """

    mkdir -p "05-calc_lev_distance"

    05-calculate_uniq_index_text_similarity.py \
    -ru "${readUniqIndex}" \
    -iu "${projectIndexFiles[1]}" \
    -iur "${projectIndexFiles[2]}" \
    -o "05-calc_lev_distance/${chunkID}.distance_matrix.tsv" \
    -rco "05-calc_lev_distance/${chunkID}.rc_collision_events.tsv"

    """
}
