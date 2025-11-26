process createDB {    
    label 'fast'
    
    conda "${params.general_env}"

    tag { "${params.projectName}.rcreateDB" }

    publishDir "${params.outDir}/02-trim_reads_2_complete_indexes/", mode: 'copy'

    when:
    !file("${params.outDir}/02-trim_reads_2_complete_indexes/${params.projectName}_db").exists()

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    path (indexFile)

    output:
    path ("${params.projectName}_db/")

    script:
    """
    mkdir -p "${params.projectName}_db"
    makeblastdb -in "${indexFile}" -dbtype nucl -out "db"
    mv db* "${params.projectName}_db/"
    """
}

process mapReads2DB {    

    conda "${params.general_env}"

    tag { "${params.projectName}.rReads2DB.${chunkID}" }

    publishDir "${params.outDir}/02-trim_reads_2_complete_indexes/", mode: 'copy', overwrite: 'false'

    when:
    !file("${params.outDir}/02-trim_reads_2_complete_indexes/").exists()

    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    maxRetries 2

    input:
    tuple val(chunkID), path(readFileFasta), path(blastDB)

    output:
    tuple val(chunkID), path ("${chunkID}.nano_trimmed_filt_2_comp_index.out")

    script:
    """
    echo "Using conda environment: ${params.general_env}"

    blastn -task blastn \
        -query "${readFileFasta}" \
        -db ${blastDB}/db \
        -perc_identity 90 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" \
        -out "${chunkID}.nano_trimmed_filt_2_comp_index.out" \
        -num_threads 8

    """
}