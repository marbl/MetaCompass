#!/usr/bin/env nextflow
// required

reads = ""
params.forward = ""
params.reverse = ""
params.unpaired = ""
paired = false
unpaired = false
params.output = ""
params.workdir = ""
params.help = false
params.minlen = "1"
params.threads  = "12"
params.log = ""

def usage(status) {
    log.info "Usage: \nnextflow run denovo_assembly.nf"

    log.info ""
    log.info "Required:"
    log.info ""

    log.info "  --forward\n     Path to forward paired-end read."
    log.info "  --reverse\n     Path to reverse paried-end read."
    log.info "  --unpaired\n    Path to unpaired read fasta file(s)."
    log.info "  --output\n      Path to output folder."

    log.info ""
    log.info "Optional:"
    log.info ""
    
    log.info "  --log\n         Path to log file." // log is at default in the denovo assembly folder
    log.info "  --workdir\n     Path to working directory. (MetaCompass Folder)" //might need it in the future, not used in this one
    log.info "  --threads\n     Number of threads to use. Default set to 12."
    log.info "  --help\n        Print help message."

    exit status
}

if (params.help){
    usage(0)
}

if( !file(params.output).isDirectory() ){
    println "ERROR: Invalid output path."
    usage(1)
}

if( params.workdir != "" && !file(params.workdir).isDirectory() ) {
    println "ERROR: Cannot find workdir."
    usage(1)
}

if( !file("$params.output/denovo_assembly").exists() && !file("$params.output/denovo_assembly").mkdir()){
    println "ERROR: Cannot create $params.output/denovo_assembly output directory!"
    usage(1)
}

// check input reads
if (params.forward != "" && params.reverse != ""){
    if(!file(params.forward).isFile() ||
       !file(params.reverse).isFile()){
        println "ERROR: Incorrect filepath to paired reads!"
        usage(1)
    }

    reads += "--forward " + file(params.forward) + " --reverse " + file(params.reverse)
    paired = true
}

if (params.unpaired != ""){
    if(!file(params.unpaired).isFile()){
        println "ERROR: Incorrect filepath to unpaired reads!"
        usage(1)
    }

    reads += " --unpaired " + file(params.unpaired)
    unpaired = true
}

if (paired == false && unpaired == false){
    println "ERROR: Incorrect combination of reads given!"
    usage(1)
}

WORKDIR = "$params.workdir"
OUTPUT = "$params.output/denovo_assembly"
LOG = file"$params.output/denovo_assembly/denovo_assembly.log"
LOG.append("Parameters: \n$params\n\n")

process megahit{
    cache 'lenient'

    shell:
    '''
    printf %s "-----Running MEGAHIT-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/megahit"
    log="!{OUTPUT}/megahit.log"
   
    cmd="!{WORKDIR}/scripts/denovo_assembly.sh  !{reads} \
                                                --minlen !{params.minlen} \
                                                --threads !{params.threads} \
                                                --output $out \
                                                --log $log"
   
    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee --append $log
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running MEGAHIT successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: MEGAHIT process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}