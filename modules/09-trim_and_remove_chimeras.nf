process trimAndRemoveChimeras {
    label 'medium'

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rTrimChimeras.${sampleName}" }

    publishDir "${params.outDir}/demultiplexed_samples", mode: 'copy', overwrite: 'true',
        saveAs: { filename -> filename.endsWith('.fna') ? filename : null }

    input:
    tuple val(sampleName), path(sampleFasta), path(blastDB), path(completeIndexesFna)

    output:
    tuple val(sampleName), path("${sampleName}.trimmed.fna"), emit: trimmedFna
    path("${sampleName}.chimera_report.tsv"),                 emit: chimeraReport

    script:
    def chimeraFlags = params.removeChimeras ?
        "--remove_chimeras true --min_chimera_coverage ${params.removeChimerasCoverage}" :
        "--remove_chimeras false"
    """
    blastn -task blastn \\
        -query ${sampleFasta} \\
        -db ${blastDB}/db \\
        -perc_identity 90 \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" \\
        -out ${sampleName}.vs_index.out \\
        -num_threads ${task.cpus}

    09-trim_and_remove_chimeras.py \\
        --fasta_reads ${sampleFasta} \\
        --blast_output ${sampleName}.vs_index.out \\
        --complete_indexes_fna ${completeIndexesFna} \\
        ${chimeraFlags} \\
        --output ${sampleName}.trimmed.fna \\
        --report ${sampleName}.chimera_report.tsv
    """
}
