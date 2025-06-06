#!/usr/bin/env nextflow

// required
reads = ""
params.forward = ""
params.reverse = ""
params.unpaired = ""
params.contigs = ""
params.metacarvel_path = ""
params.output = ""
params.workdir = ""
params.help = false

def usage(status) {
    log.info "Usage: \nnextflow run scaffolding.nf"

    log.info ""
    log.info "Required:"
    log.info ""

    log.info "  --forward               Path to forward reads."
    log.info "  --reverse               Path to reverse reads."
    log.info "  --unpaired              Path to unpaired reads."
    log.info "  --contigs               Path to contigs file."
    log.info "  --metacarvel_path       Path to metacarvel <run.py>."
    log.info "  --workdir               Path to working directory. (MetaCompass Folder)"
    log.info "  --log                   Path to log file."
    log.info "  --output                Path to output folder"

    log.info ""
    log.info "Optional:"
    log.info ""
    
    log.info "  --help                  Print help message"

    exit status
}

if (params.help){
    usage(0)
}

// TODO: Validate required set of reads given
if( params.forward != "" && params.reverse != "" ){

    reads += "--forward $params.forward --reverse $params.reverse "
}

if( params.unpaired != "" ){
    reads += "--unpaired $params.unpaired"
}

if( reads == "" ){
    println "ERROR: No reads given."
    usage(1)
}

if( !file(params.contigs).isFile() ){
    println "ERROR: No contig file given."
    usage(1)
}

if( !file(params.metacarvel_path).exists() ){
    println "ERROR: Need path to MetaCarvel executable."
    usage(1)
}

if( !file(params.output).isDirectory() || \
    !file(params.workdir).isDirectory() ){
    println "ERROR: Please check required inputs."
    usage(1)
}

if ( ! (params.output != "" && file(params.output).isDirectory() ) ) {
    println "ERROR: Output folder does not exist"
    usage(1)
}

if( ! file("$params.output/scaffolding").exists() && ! file("$params.output/scaffolding").mkdir()){
    println "ERROR: Cannot create output directory for reference assembly"
    usage(1)
}

WORKDIR = "$params.workdir"
OUTPUT = "$params.output/scaffolding"
LOG = file("$params.output/scaffolding/scaffolding.log")
LOG.append("Parameters: \n$params\n\n")

process scaffolding {
    cache 'lenient'

    shell:
    '''
    printf %s "-----Running MetaCarvel for Scaffolding-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/metacarvel"
    if [[ ! -d $out ]]; then
        mkdir $out
    fi
    log=!{OUTPUT}"/metacarvel.log"

    cmd="bash !{WORKDIR}/scripts/scaffolding.sh    !{reads} \
                                                   --contigs !{params.contigs} \
                                                   --metacarvel_path !{params.metacarvel_path} \
                                                   --out $out \
                                                   --log $log"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a $log
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running MetaCarvel successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: MetaCarvel process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}
