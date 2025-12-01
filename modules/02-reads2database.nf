process createDB {    
    label 'fast'
    
    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rcreateDB" }

    publishDir "${params.outDir}/02-trim_reads_2_complete_indexes/", mode: 'copy'

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

    params.conda_env ?: "${projectDir}/envs/deluxpore.yml"

    tag { "${params.projectName}.rReads2DB.${chunkID}" }

    publishDir "${params.outDir}/02-trim_reads_2_complete_indexes/", mode: 'copy', overwrite: 'false'

    // when:
    // !file("${params.outDir}/02-trim_reads_2_complete_indexes/").exists()

    input:
    tuple val(chunkID), path(readFileFasta), path(blastDB)

    output:
    tuple val(chunkID), path ("${chunkID}.nano_trimmed_filt_2_comp_index.out")

    script:
    """
    blastn -task blastn \
        -query "${readFileFasta}" \
        -db ${blastDB}/db \
        -perc_identity 90 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand" \
        -out "${chunkID}.nano_trimmed_filt_2_comp_index.out" \
        -num_threads ${task.cpus}

    """
}