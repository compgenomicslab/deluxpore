process concatenateSamples {
    label 'fast'
    
    tag { "${params.projectName}.rconcatenateSamples.${sampleName}" }

    publishDir "${params.outDir}/07-final_demultiplex", mode: 'copy'
    
    input:
    tuple val(sampleName), path(sampleFiles)
    
    output:
    tuple val(sampleName), path("${sampleName}.trimmed.fna")
    
    script:
    """
    cat ${sampleFiles.join(' ')} > ${sampleName}.trimmed.fna
    """
}