process generateIndexFiles {    
    label 'fast'

    conda params.conda_env ?: "${projectDir}/envs/deluxpore.yml"
    
    tag { "${params.projectName}.rGenIndex" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false'

    when:
    !file("${params.outDir}/indexes").exists()

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

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
