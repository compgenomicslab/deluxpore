process concatenateSamples {
    label 'fast'
    
    tag { "${params.projectName}.rconcatenateSamples.${sampleName}" }

    publishDir "${params.outDir}/demultiplexed_samples", mode: 'copy'
    
    input:
    tuple val(sampleName), path(sampleFiles)
    
    output:
    tuple val(sampleName), path("${sampleName}*fna")
    
    script:
    def suffix = params.trimmIlluminaIndexes ? ".trimmed.fna" : ".fna"
    """
    cat ${sampleFiles.join(' ')} > ${sampleName}${suffix}
    """
}