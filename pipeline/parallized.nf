#!/usr/bin/env nextflow

// Project folders
BIN = "$workflow.projectDir/../bin"

Channel
    .fromPath("${params.output}/*.fna")
    .set { fasta_files_ch }
Channel
    .fromPath("${params.output}/interleaved_reads.fq")
    .set { fastq_file_ch }


/*
 * The ParallelKMC process performs k-mer counting on each FASTA file from a specified
 * directory using the KMC tool. Files named 'unplaced.scaf.fna' are skipped.
 * Intermediate directories created during the process are cleaned up after k-mer counting.
 */
process ParallelKMC {
    input:
    file fastaFile from fasta_files_ch
    script:
    """
    name=\$(basename $fastaFile .fna)
    genomes_g_2="\${name}.kmer"
    mkdir -p ${params.output}/\${name}
    if [[ \${name} != "unplaced.scaf.fna" ]]; then
        "${BIN}/kmc/kmc" -k28 -ci1 -hp -fm $fastaFile ${params.output}/\${genomes_g_2} ${params.output}/\${name}
    fi
    rm -rf  ${params.output}/\${name}
    """
}
