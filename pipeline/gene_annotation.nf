#!/usr/bin/env nextflow

process annotateContigGenes {

    publishDir "${workflow.outputDir}/output", mode: 'copy'

    input:
    path refctgs

    output:
    path "contig_gene_annotations.tsv"
    path "all_contigs.gff3"
    path "missing_gff_files.txt"
    path "unparsed_contig_headers.txt"

    script:
    """
    GFF_BASE="${params.reference_db}/collect_refs/ncbi_dataset/data"

    python ${projectDir}/scripts/annotate_contig_genes.py \\
        --gff-base "\$GFF_BASE" \\
        --output-tsv contig_gene_annotations.tsv \\
        --output-gff all_contigs.gff3 \\
        --refctgs ${refctgs}
    """
}