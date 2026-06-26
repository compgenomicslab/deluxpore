process concatenateSamples {
    label 'fast'

    tag { "${params.projectName}.rconcatenateSamples.${sampleName}" }

    publishDir "${params.outDir}/demultiplexed_samples", mode: 'copy'

    input:
    tuple val(sampleName), path(sampleFiles)

    output:
    tuple val(sampleName), path("${sampleName}*.fna")

    script:
    def suffix = params.trimmIlluminaIndexes ? ".trimmed.fna" : ".fna"
    """
    cat ${sampleFiles.join(' ')} > ${sampleName}${suffix}
    """
}

process concatenateSummaries {
    label 'fast'

    tag { "${params.projectName}.rconcatenateSummaries.${ambiguityType}" }

    publishDir "${params.outDir}/ambiguous_reads_report", mode: 'copy'

    input:
    tuple val(ambiguityType), path(chunkFastas)

    output:
    tuple val(ambiguityType), path("${ambiguityType}.fna")

    script:
    """
    cat ${chunkFastas.join(' ')} > ${ambiguityType}.fna
    """
}

process concatenateAmbiguousReport {
    label 'fast'

    tag { "${params.projectName}.rconcatenateAmbiguousReport" }

    publishDir "${params.outDir}/ambiguous_reads_report", mode: 'copy'

    input:
    path(tsvFiles)

    output:
    path("ambiguous_reads.tsv")

    script:
    """
    files=(\$(ls ambiguous_reads.*.tsv | sort -V))
    head -1 "\${files[0]}" > ambiguous_reads.tsv
    for f in "\${files[@]}"; do tail -n +2 "\$f"; done >> ambiguous_reads.tsv
    """
}
