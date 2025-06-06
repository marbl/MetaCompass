#!/usr/bin/env nextflow

params.output = ""
params.workdir = ""
params.threads = 8
params.help = false

def usage(status) {
    log.info "Usage: \nnextflow run de_novo.nf"
    log.info ""
    log.info "Required:"
    log.info ""
    log.info "  --workdir               Path to working directory"
    log.info "  --output                Path to output folder"
    log.info ""
    log.info "Optional:"
    log.info ""
    log.info "  --threads               Number of threads for MEGAHIT (default: 8)"
    log.info "  --help                  Print help message"
    exit status
}

if (params.help) {
    usage(0)
}

if (!file(params.output).isDirectory() || 
    !file(params.workdir).isDirectory()) {
    println "ERROR: Please check required inputs."
    usage(1)
}

// Process to find unmapped reads
process find_unmapped {
    cache 'lenient'
    
    output:
    path 'paths.txt' into unmapped_paths
    
    script:
    """
    mkdir -p ${params.output}/denovo_assembly
    find ${params.output}/reference_assembly -name "unmapped.fq" > paths.txt
    if [ ! -s paths.txt ]; then
        echo "No unmapped files found"
        exit 1
    fi
    """
}

// Process for denovo assembly
process denovo_assembly {
    cache 'lenient'
    publishDir "${params.output}/denovo_assembly", mode: 'copy', overwrite: true
    
    input:
    path 'paths.txt' from unmapped_paths
    
    output:
    path "megahit_output_denovo"
    
    script:
    """
    # Create output directory
    mkdir -p ${params.output}/denovo_assembly

    # Get the unmapped file path
    unmapped_path=\$(cat paths.txt)
    
    # Create temporary directory for read splitting
    mkdir -p tmp_reads
    
    # Copy unmapped reads
    cp \$unmapped_path tmp_reads/unpaired_reads.fastq
    
    # Run MEGAHIT
    megahit \\
        -r tmp_reads/unpaired_reads.fastq \\
        -o megahit_output_denovo \\
        -t ${params.threads}
    
    # Ensure output directory exists and is accessible
    if [ ! -d megahit_output_denovo ]; then
        echo "MEGAHIT failed to create output directory"
        exit 1
    fi
    
    # Clean up temporary files
    rm -rf tmp_reads
    """
}
