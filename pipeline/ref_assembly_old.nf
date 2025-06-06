#!/usr/bin/env nextflow

// required
params.forward = ""
params.reverse = ""
params.unpaired = ""
params.aligned_reads = ""
params.references = ""
params.pilon = ""
params.output = ""
params.workdir = ""
params.help = false

// optional
// params.threads = "12"
// params.flank = "5"
// params.mindepth = "3"
// params.memory = "20"

def usage(status) {
    log.info "Usage: \nnextflow run ref_assembly.nf"

    log.info ""
    log.info "Required:"
    log.info ""

    log.info "  --aligned_reads      Path to 'align_reads' output folder from reference_culling."
    log.info "  --references         Path to folder containing references."
    log.info "                       Format is that of downloading references using NCBI 'datasets' tool."
    log.info "  --pilon              Path to pilon executable. Default: pilon-1.23.jar"
    log.info "  --valet              Path to valet"
    log.info "  --workdir            Path to working directory. (MetaCompass Folder)"
    log.info "                       Do not include / at the end of the folder"
    log.info "  --output             Path to output folder"
    log.info "                       Do not include / at the end of the folder"

    log.info ""
    log.info "Optional:"
    log.info ""

    log.info "  --help               Print help message"
    log.info "  --run_valet          Whether or not to run VALET."
    log.info "  --flank              Flank arg for pilon (default: params.flank or 5)"
    log.info "  --threads            Threads arg for pilon (default: params.threads or 12)"
    log.info "  --mindepth           Mindepth arg for pilon (default: params.mindepth or 3)"
    log.info "  --memory             Max memory (in GB) assigned for pilon call (default: params.memory or 20)"

    exit status
}

if (params.help){
    usage(0)
}

if ( ! (params.aligned_reads != "" && file(params.aligned_reads).isDirectory() )){
    println "ERROR: Alignments does not exist"
    usage(1)
}

if ( ! (params.references != "" && file(params.references).isDirectory() ) ){
    println "ERROR: Reference folder does not exist"
    usage(1)
}

if ( ! (params.pilon != "" && file(params.pilon).exists() ) ){
    println "ERROR: Pilon jar file does not exist"
    usage(1)
}

if ( ! (params.output != "" && file(params.output).isDirectory() ) ) {
    println "ERROR: Output folder for MetaCompass does not exist"
    usage(1)
}

// if( ! file("$params.output/reference_assembly").exists() && ! file("$params.output/reference_assembly").mkdir() ){
//     println "ERROR: Cannot create output directory for reference assembly"
//     usage(1)
// }

WORKDIR = "$params.workdir"
OUTPUT = "$params.output/reference_assembly"
LOG = file("$params.output/reference_assembly/reference_assembly.log")
LOG.append("Parameters: \n$params\n\n")

process valet {
    cache 'lenient'

    output:
    env VALET_REFS into valet

    shell:
    '''
    printf %s "-----Running VALET-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/VALET"
    if [[ ! -d $out ]]; then
        mkdir $out
    fi
    log=!{OUTPUT}"/valet.log"

    cmd="bash !{WORKDIR}/scripts/valet.sh    --reads !{params.aligned_reads} \
                                             --references !{params.references} \
                                             --valet !{params.valet} \
                                             --run_valet !{params.run_valet} \
                                             --out $out \
                                             --log $log"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a $log
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running VALET successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: VALET process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi

    VALET_REFS=$out
    '''
}

// TODO:
// Add script that uses valet output to modify references.
process modify_refs {
    cache 'lenient'

    input:
    env VALET_REFS from valet

    output:
    env MODIFIED_REFS into modify_refs

    shell:
    '''
    # currently does not use valet output
    MODIFIED_REFS=$VALET_REFS
    '''

}

/*
This step will call pilon.sh. Currently assume we call pilon right after ref culling
Will need to fix this later to add VALLEY process before PILON process
*/
process pilon {
    cache 'lenient'

    input:
    env MODIFIED_REFS from modify_refs

    shell:
    '''
    printf %s "-----Running pilon-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/pilon"
    if [[ ! -d $out ]]; then
        mkdir $out
    fi
    log=$out"/pilon.log"

    cmd="bash !{WORKDIR}/scripts/pilon.sh   --alignment $MODIFIED_REFS \
                                            --reference !{params.references}\
                                            --pilon-path !{params.pilon} \
                                            --flank !{params.flank} \
                                            --threads !{params.threads} \
                                            --mindepth !{params.mindepth} \
                                            --memory !{params.memory} \
                                            --out $out \
                                            --log $log"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a $log
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running pilon successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: Pilon process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}
