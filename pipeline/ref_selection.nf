#!/usr/bin/env nextflow

// required
reads = ""
params.forward = "" // initialize, so that params.forward != null when modifying "reads"
params.reverse = ""
params.unpaired = ""
params.output = ""
params.workdir = ""
params.help = false
params.filter_refs = false
def usage(status) {
    log.info "Usage: \nnextflow run ref_selection.nf \n\
            [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] | \n\
            [--forward /path/to/forwardReads --reverse /phat/to/reverseReads] | \n\
            [--unpaired /path/to/unpairedReads]] \n\
            -output-dir /path/to/outputDir"

    log.info ""
    log.info "Required:"
    log.info ""

    log.info " --forward\n                  Path to forward read."
    log.info " --reverse\n                  Path to reverse read."
    log.info " --unpaired\n                 Paths to unpaired read."
    log.info " -output-dir\n                   Path to output folder."
    log.info " --log\n                      Path to log file."
    log.info " --workdir\n                  Path to working directory. (MetaCompass Folder)"
    log.info " --reference\n                Path to reference file. Default selected by kmer-mask."
    log.info " --ms\n                       Kmer size."
    log.info " --clean\n                    Kmer clean percentage."
    log.info " --match\n                    Kmer match percentage."
    log.info " --readlen\n                  Read length. Default set to 105."
    log.info " --depth_of_coverage          Deapth of coverge. Default set to 5."
    log.info " --breadth_of_coverage        Breadth of coverage. Default set to 0.9."
    log.info " --percent_marker_covered\n   Percentage of marker coverage."
    log.info " --threads\n                  Number of threads to use. Default set to 12."

    log.info ""
    log.info "Optional:"
    log.info ""
    
    log.info " --help\n                     Print help message."

    exit status
}

if (params.help){
    usage(0)
}

if( params.forward != "" && params.reverse != "" ){

    reads += "--forward $params.forward --reverse $params.reverse "
}

if( params.unpaired != "" ){
    reads += "--unpaired $params.unpaired"
}

OUTPUT = "${workflow.outputDir}/reference_selection"

// keeps only reads that map to marker genes
process filter_reads { 
    publishDir {
        path: file("${workflow.outputDir}/reference_selection"),
        mode: 'copy'
    }

    input:
        val reads

    output:
        path "sorted_mapped_reads.fq"

    script:
    """
        minimap2 -ax sr -t ${params.threads} ${params.reference_db}/marker_index/marker_clustered.fna $reads | \
        samtools view -bS -F 4 - | \
        samtools bam2fq - | \
        seqkit sort - > sorted_mapped_reads.fq
    """
}

// Map reads to each marker gene
process map_to_gene{
    publishDir { 
       path: file("${workflow.outputDir}/reference_selection"), 
       mode:'copy'
    }

    input:
       path reads
       val gene

    output:
       file "${gene}_marker_cov_per_genome.txt"

    script:
    """
        minimap2 -ax sr -t ${params.threads} ${params.reference_db}/marker_index/${gene}/${gene}_clustered.fna ${reads} | \
        samtools view -b > ${gene}.match.bam
        samtools sort ${gene}.match.bam -o ${gene}.match.sorted.bam

        # capture marker gene lengths
       # tail -n+2 ${params.reference_db}/marker_index/${gene}/${gene}_genome.tsv | \
       # cut -f 1,4 > ${gene}_genome.txt

        # calculate coverage
        bedtools genomecov -ibam ${gene}.match.sorted.bam \
        -max ${params.depth_of_coverage} > ${gene}_genomeCov.txt
        
        # select markers with good breadth of coverage at given depth of coverage
        echo -e "seq\tcoverage" > ${gene}_marker_cov.txt
        awk -v doc=${params.depth_of_coverage} \
            -v boc=${params.breadth_of_coverage} \
            '{if (\$2 == doc && \$5 >= boc){print \$1"\t"\$5}}' ${gene}_genomeCov.txt >> \
            ${gene}_marker_cov.txt

        # pull out the genomes corresponding to all marker genes covered
        python3 ${projectDir}/scripts/marker_coverage_per_candidate_genome.py ${gene} \
        ${gene}_marker_cov.txt \
        ${gene}_marker_cov_per_genome.txt \
        ${params.reference_db}/marker_index

        # output should now be in ${gene}_marker_cov_per_genome.txt
    """
}

// Combine results from all the genes to select genomes with enough coverage
process select_genomes {
    publishDir { 
       path: file("${workflow.outputDir}/reference_selection"), 
       mode:'copy'
    }

    input:
       path "*_marker_cov_per_genome.txt"
       
    output:
       path "reference_candidates.txt", emit: candidates
       path "ref_genome_marker_gene_coverage.tsv"

    script:
    """
       cat *_marker_cov_per_genome.txt > marker_cov_per_genome.txt
       python3 ${projectDir}/scripts/identify_candidate_genomes.py ${params.percent_markers_covered} . ${params.reference_db}
    """
}

