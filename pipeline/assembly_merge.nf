#!/usr/bin/env nextflow

// required
params.denovo = ""
params.pilon = ""
params.output = ""
params.workdir = ""
params.help = false


def usage(status) {
    log.info "Usage: \nnextflow run assembly_merge.nf"

    log.info ""
    log.info "Required:"
    log.info ""

    log.info "  --denovo        Path to config file of denovo assembly"
    log.info "  --pilon         Path to folder containing config files of reference-guided assembly"
    log.info "  --output        Path to output folder"
    log.info "                  Do not include / at the end of the folder"
    log.info "  --workdir       Path to working directory. (MetaCompass Folder)"
    log.info "                  Do not include / at the end of the folder"

    log.info ""
    log.info "Optional:"
    log.info ""

    log.info "  --help          Print help message"

    exit status
}

if (params.help){
    usage(0)
}

if ( ! (params.workdir != "" && file(params.workdir).isDirectory() ) ){
    println "Path to metacompass workdir: $params.workdir"
    println "ERROR: Metacompass workdir does not exist"
    usage(1)
}

MERGE_SCRIPT = "$params.workdir/scripts/merge.py"
if ( ! file(MERGE_SCRIPT).isFile() ){
    println "Path to workdir: $params.workdir"
    println "ERROR: merge.py cannot be found in <workdir>/scripts"
    usage(1)
}

if ( ! (params.denovo != "" && file(params.denovo).isFile() ) ){
    println "Path to denovo config file: $params.denovo"
    println "ERROR: Denovo config file does not exist"
    usage(1)
}

if ( ! (params.pilon != "" && file(params.pilon).isDirectory() ) ){
    println "Path to pilon config folder: $params.pilon"
    println "ERROR: Folder containing pilon config files does not exist"
    usage(1)
}

if ( ! (params.output != "" && file(params.output).isDirectory() ) ) {
    println "ERROR: Output folder does not exist"
    usage(1)
}

if( ! file("$params.output/assembly_merge").mkdir()){
    println "ERROR: Cannot create output directory for assembly_merge"
    usage(1)
}

OUTPUT = "$params.output/assembly_merge"
LOG = file("$params.output/assembly_merge/assembly_merge.log")
LOG.append("Parameters: \n$params\n\n")

process merge {
    cache 'lenient'

    shell:
    '''
    printf %s "-----Running assembly merge-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/contigs.merge.fasta"
    touch $out

    cmd="python3 !{MERGE_SCRIPT} !{params.denovo} \
                                 !{params.pilon} \
                                 $out"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a !{LOG}
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running merge process successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: Merge process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}