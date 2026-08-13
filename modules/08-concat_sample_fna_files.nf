process concatenateSamples {
    label 'fast'

    tag { "${params.projectName}.rconcatenateSamples.${sampleName}" }

    publishDir "${params.outDir}/demultiplexed_samples", mode: 'copy'

    input:
    tuple val(sampleName), path(sampleFiles)

    output:
    tuple val(sampleName), path("${sampleName}${params.trimmIlluminaIndexes ? '.trimmed' : ''}.fna")

    script:
    def outFile = "${sampleName}${params.trimmIlluminaIndexes ? '.trimmed' : ''}.fna"
    """
    cat ${sampleFiles.join(' ')} > ${outFile}
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

process concatenateChimeraReports {
    label 'fast'

    tag { "${params.projectName}.rconcatenateChimeraReports" }

    publishDir "${params.outDir}/ambiguous_reads_report", mode: 'copy'

    input:
    path(tsvFiles)

    output:
    path("chimera_reads.tsv")

    script:
    """
    files=(\$(ls *.chimera_report.tsv | sort -V))
    head -1 "\${files[0]}" > chimera_reads.tsv
    for f in "\${files[@]}"; do tail -n +2 "\$f"; done >> chimera_reads.tsv
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

process concatenateIndexAssignmentSummary {
    label 'fast'

    tag { "${params.projectName}.rconcatenateIndexAssignmentSummary" }

    publishDir "${params.outDir}/ambiguous_reads_report", mode: 'copy'

    input:
    path(tsvFiles)

    output:
    path("index_assignment_summary.tsv")

    script:
    """
    files=(\$(ls index_assignment_summary.*.tsv | sort -V))
    head -1 "\${files[0]}" > index_assignment_summary.tsv
    for f in "\${files[@]}"; do tail -n +2 "\$f"; done >> index_assignment_summary.tsv

    awk -F'\\t' '{both+=\$2;i5+=\$3;i7+=\$4;un+=\$5;tot+=\$6} \
        END{print "TOTAL\\t"both"\\t"i5"\\t"i7"\\t"un"\\t"tot}' \
        <(tail -n +2 index_assignment_summary.tsv) >> index_assignment_summary.tsv
    """
}
