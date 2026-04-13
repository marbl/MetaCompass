#!/usr/bin/env nextflow

// -----------------------------------------------------------------------------
// Default parameters
// -----------------------------------------------------------------------------
params.forward      = params.forward ?: ''
params.reverse      = params.reverse ?: ''
params.unpaired     = params.unpaired ?: ''
params.output       = params.output ?: ''
params.help         = params.help ?: false
params.filter_refs  = params.filter_refs ?: false

// -----------------------------------------------------------------------------
// Helper values
// -----------------------------------------------------------------------------
final String REFERENCE_SELECTION_DIR = params.output ?
    "${params.output}/reference_selection" :
    "reference_selection"

// Build read argument string for standalone use
String reads = ''
List<String> readInputs = []

if (params.forward && params.reverse) {
    readInputs << params.forward
    readInputs << params.reverse
}

if (params.unpaired) {
    readInputs << params.unpaired
}

reads = readInputs.join(' ')

// -----------------------------------------------------------------------------
// Usage
// -----------------------------------------------------------------------------
def usage = { int status ->
    log.info """
    Usage:
      nextflow run ref_selection.nf \\
        [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] |
         [--forward /path/to/forwardReads --reverse /path/to/reverseReads] |
         [--unpaired /path/to/unpairedReads]] \\
        --output /path/to/outputDir

    Required:

      --forward                  Path to forward read
      --reverse                  Path to reverse read
      --unpaired                 Path to unpaired read(s)
      --output                   Path to output folder
      --reference_db             Path to reference database
      --ms                       K-mer size
      --clean                    K-mer clean percentage
      --match                    K-mer match percentage
      --readlen                  Read length (default set elsewhere)
      --depth_of_coverage        Depth of coverage
      --breadth_of_coverage      Breadth of coverage
      --percent_markers_covered  Percentage of marker coverage
      --threads                  Number of threads to use

    Optional:

      --help                     Print help message
    """.stripIndent().trim()

    exit status
}

if (params.help) {
    usage(0)
}

// -----------------------------------------------------------------------------
// Processes
// -----------------------------------------------------------------------------

// Keeps only reads that map to marker genes
process filter_reads {
    publishDir "${REFERENCE_SELECTION_DIR}", mode: 'copy'

    input:
    val reads

    output:
    path "sorted_mapped_reads.fq"

    script:
    """
    minimap2 -ax sr -t ${params.threads} ${params.reference_db}/marker_index/marker_clustered.fna ${reads} | \
    samtools view -bS -F 4 - | \
    samtools bam2fq - | \
    seqkit sort - > sorted_mapped_reads.fq
    """
}

// Map reads to each marker gene
process map_to_gene {
    publishDir "${REFERENCE_SELECTION_DIR}", mode: 'copy'

    input:
    path reads
    val gene

    output:
    path "${gene}_marker_cov_per_genome.txt"

    script:
    """
    minimap2 -ax sr -t ${params.threads} ${params.reference_db}/marker_index/${gene}/${gene}_clustered.fna ${reads} | \
    samtools view -b > ${gene}.match.bam

    samtools sort ${gene}.match.bam -o ${gene}.match.sorted.bam

    # Calculate coverage
    bedtools genomecov -ibam ${gene}.match.sorted.bam \
        -max ${params.depth_of_coverage} > ${gene}_genomeCov.txt

    # Select markers with good breadth of coverage at the given depth threshold
    echo -e "seq\tcoverage" > ${gene}_marker_cov.txt
    awk -v doc=${params.depth_of_coverage} -v boc=${params.breadth_of_coverage} \
        '{if (\$2 == doc && \$5 >= boc){print \$1"\\t"\$5}}' \
        ${gene}_genomeCov.txt >> ${gene}_marker_cov.txt

    # Pull out the genomes corresponding to all covered marker genes
    python3 ${projectDir}/scripts/marker_coverage_per_candidate_genome.py \
        ${gene} \
        ${gene}_marker_cov.txt \
        ${gene}_marker_cov_per_genome.txt \
        ${params.reference_db}/marker_index
    """
}

// Combine results from all genes to select genomes with enough coverage
process select_genomes {
    publishDir "${REFERENCE_SELECTION_DIR}", mode: 'copy'

    input:
    path "*_marker_cov_per_genome.txt"

    output:
    path "reference_candidates.txt", emit: candidates
    path "ref_genome_marker_gene_coverage.tsv"

    script:
    """
    cat *_marker_cov_per_genome.txt > marker_cov_per_genome.txt

    python3 ${projectDir}/scripts/identify_candidate_genomes.py \
        ${params.percent_markers_covered} \
        . \
        ${params.reference_db}
    """
}
