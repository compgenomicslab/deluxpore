process generateIndexFiles {    
    label 'fast'

    conda "${params.general_env}"
    
    tag { "${params.projectName}.rGenIndex" }

    publishDir "${params.outDir}/", mode: 'copy', overwrite: 'false'

    when:
    !file("${params.outDir}/00-indexes").exists()

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    file (experimentalDesignInput)

    output:
    path ("00-indexes/*")

    script:
    """
    mkdir -p "00-indexes"
    00-generate_index_files.py \
    -i ${experimentalDesignInput} \
    -o "00-indexes"
    """
}
