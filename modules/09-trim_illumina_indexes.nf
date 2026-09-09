process trimIlluminaIndexes {
    label 'medium'

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rTrimIllumina.${sampleName}" }

    publishDir "${params.outDir}/demultiplexed_samples", mode: 'copy', overwrite: 'true',
        saveAs: { filename -> filename.endsWith('.fna') ? filename : null }

    input:
    tuple val(sampleName), path(sampleFasta), path(blastDB)

    output:
    tuple val(sampleName), path("${sampleName}.trimmed.fna")

    script:
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
    # concatenateSamples names it that way because trimmIlluminaIndexes is on)
    # and is staged in as a symlink to the upstream work dir. Writing directly
    # to that name while 09-trim_illumina_indexes.py is still lazily reading
    # --fasta_reads (Bio.SeqIO.index) would truncate the file through the
    # symlink mid-read, corrupting it.
    09-trim_illumina_indexes.py \\
        --fasta_reads ${sampleFasta} \\
        --blast_output ${sampleName}.vs_index.out \\
        --min_fragment_length ${params.nanoLength} \\
        --output ${sampleName}.trimmed.fna.tmp

    mv ${sampleName}.trimmed.fna.tmp ${sampleName}.trimmed.fna
    """
}
