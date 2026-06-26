process generateIndexFiles {    
    label 'fast'

    conda params.conda_env ?: "${projectDir}/envs/deluxpore.yml"
    
    tag { "${params.projectName}.rGenIndex" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false', enabled: params.publishIntermediate

    input:
    tuple path(experimentalDesignInput), path(completeIndexes), path(uniqueIndexes)

    output:
    path ("indexes/*")

    script:
    """

    mkdir -p "indexes"
    00-generate_index_files.py \
    -i ${experimentalDesignInput} \
    -ci ${completeIndexes} \
    -ui ${uniqueIndexes} \
    -o "indexes"
    """
}
