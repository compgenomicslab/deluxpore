process generateIndexFiles {    
    label 'fast'

    conda params.conda_env ?: "${projectDir}/envs/deluxpore.yml"
    
    tag { "${params.projectName}.rGenIndex" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false'

    input:
    file (experimentalDesignInput)

    output:
    path ("indexes/*")

    script:
    """

    mkdir -p "indexes"
    00-generate_index_files.py \
    -i ${experimentalDesignInput} \
    -ik ${params.libraryIndexSeqs} \
    -o "indexes"
    """
}
