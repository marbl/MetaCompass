#!/usr/bin/env nextflow

// required
params.reference_db = ""
params.output = ""

def usage(status) {
    log.info "Usage: \nnextflow run initialize.nf --reference_db /path/to/reference_db --output /path/to/outputDir"

    log.info "Required:"

    log.info " --output           Path to output folder."
    log.info " --reference_db     Path to reference database."

    log.info "Optional:"

    log.info " --help    Print usage information."
    
    exit status
}

if (params.help){
    usage(0)
}

// invalid reference or output directory
// TODO: Allow users to run if reference_db not provided
if(!file(params.output).isDirectory()){
    log.info "Invalid output directory"
    log.info "Output:   \t$params.output"
    
    usage(1)
}

if(!file(params.reference_db).isDirectory()){
    log.info "Invalid reference directory"
    log.info "Reference database:\t$params.reference_db"

    usage(1)
}


SCRIPTS = "$workflow.projectDir/../scripts"
LOG = file("$params.output/initialize.log")
LOG.append("Parameters: \n$params\n\n")

process check_dependencies{
    cache 'lenient'

    output:
    val done

    exec:
	done = true

    script:
    '''
    printf %s "-----Running check_dependencies-----\n$(date)\n\n" >> !{LOG}

    cmd="!{SCRIPTS}/check_exist.sh -b bowtie2,samtools,python3,datasets,bedtools,jellyfish"
    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd
    '''
}

workflow {
    check_dependencies
}

