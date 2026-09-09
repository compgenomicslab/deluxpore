process detectChimeras {
    label 'medium'

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rDetectChimeras.${chunkID}" }

    publishDir "${params.outDir}/01b-detect_chimeras/", mode: 'copy', enabled: params.publishIntermediate

    input:
    tuple val(chunkID), path(readFasta), path(blastDB), path(completeIndexesFna)

    output:
    tuple val(chunkID), path("${chunkID}.chimera_checked.fna"), emit: reads
    path("${chunkID}.chimera_report.tsv"),                      emit: chimeraReport

    script:
    def removeFlag = params.removeChimeras ? "true" : "false"
    """
    blastn -task blastn \\
        -query ${readFasta} \\
        -db ${blastDB}/db \\
        -perc_identity 90 \\
        -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" \\
        -out ${chunkID}.vs_index.out \\
        -num_threads ${task.cpus}

    01b-detect_chimeras.py \\
        --fasta_reads ${readFasta} \\
        --blast_output ${chunkID}.vs_index.out \\
        --complete_indexes_fna ${completeIndexesFna} \\
        --remove_chimeras ${removeFlag} \\
        --min_chimera_coverage ${params.removeChimerasCoverage} \\
        --min_fragment_length ${params.nanoLength} \\
        --output ${chunkID}.chimera_checked.fna \\
        --report ${chunkID}.chimera_report.tsv
    """
}
