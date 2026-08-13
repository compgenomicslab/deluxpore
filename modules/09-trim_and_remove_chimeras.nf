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

    # --output is written to a temp name and moved into place afterwards:
    # ${sampleFasta} can be identically named ${sampleName}.trimmed.fna (when
    # trimmIlluminaIndexes concatenates samples under that name already), and
    # it is staged in as a symlink to the upstream work dir. Writing directly
    # to that name while 09-trim_and_remove_chimeras.py is still lazily
    # reading --fasta_reads (Bio.SeqIO.index) would truncate the file through
    # the symlink mid-read, corrupting it and clobbering the upstream cache.
    09-trim_and_remove_chimeras.py \\
        --fasta_reads ${sampleFasta} \\
        --blast_output ${sampleName}.vs_index.out \\
        --complete_indexes_fna ${completeIndexesFna} \\
        ${chimeraFlags} \\
        --output ${sampleName}.trimmed.fna.tmp \\
        --report ${sampleName}.chimera_report.tsv

    mv ${sampleName}.trimmed.fna.tmp ${sampleName}.trimmed.fna
    """
}
