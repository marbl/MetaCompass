#!/usr/bin/env nextflow

reads = ""
params.forward = ""
params.reverse = ""
params.unpaired = ""
gzip_flag = false
paired = false
unpaired = false
forward_gz = ""
reverse_gz = ""
unpaired_gz = ""
params.reference_db = ""
params.output = ""
params.skip_rs = false
params.skip_rc = false
params.clean_uf = false
params.de_novo = 1



def usage(status) {
    log.info "Usage: \nnextflow run metacompass.nf \n\
            [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] | \n\
            [--forward /path/to/forwardReads --reverse /path/to/reverseReads] | \n\
            [--unpaired /path/to/unpairedReads]] \n\
            --output /path/to/outputDir"
    
    log.info ""
    log.info "Required:"
    log.info ""

    log.info " --forward        Path to forward paired-end read."
    log.info " --reverse        Path to reverse paried-end read."
    log.info " --unpaired       Path to unpaired read fasta file(s)."
    log.info " --output         Path to output folder."
    
    log.info ""
    log.info "Optional:"
    log.info ""
    
    log.info " --reference_db   Path to reference marker gene database."
    log.info " --ref_sel        Reference selection method. Default set to 'tax'."
    log.info " --ref_pick       Pick reference selection method. Default set to 'breadth'."
    log.info " --readlen        Read length for filtering. Default set to 200."
    log.info " --mincov         Minimum coverage. Default set to 1."
    log.info " --minctglen      Minimum contig length. Default set to 1."
    log.info " --run_valet      Whether or not run VALET during reference-guided assembly."
    log.info " --de_novo        Set to 0 to skip de_novo assembly. Default set to 1."
    log.info " --skip_rs        Set to true to skip ref selection process. Default set to false."
    log.info " --skip_rc        Set to true to skip ref culling process. Default set to false."
    log.info " --tracks         Tracks. Default set to false."
    log.info " --threads        Number of threads to use. Default set to 12."
    log.info " --memory         Amount of memory to use."
    log.info " --clean_uf       Remove unnecessary files while running the job. Default is false."
    log.info " --executor       Executor to use."
    log.info " --help           Print help message."

    exit status
}

if (params.help){
    usage(0)
}

// TODO: implement resuming option with 
//       previously created output folder


// if the path is invalid, throw an error
if (params.output == "" || !file(params.output).mkdirs() ){
    println "ERROR: Need proper output directory path!"
    usage(1)
}

// check if reference file exists
if (params.reference_db != "" && !file(params.reference_db).isDirectory() ){
    println "ERROR: Reference genome files not found!"
    usage(1)
}

String getExtension(String path){
    file = file(path)
    String Extension = file.getExtension()
    return Extension
}


// check input reads
if (params.forward != "" && params.reverse != ""){
        if(!file(params.forward).isFile() ||
       !file(params.reverse).isFile()){
        println "ERROR: Incorrect filepath to paired reads!"
        usage(1)
    }
    String forwardExtension = getExtension(params.forward)
    String reverseExtension = getExtension(params.reverse)

    if(forwardExtension == "gz" && reverseExtension == "gz"){

        // save gzip reads to file
        gzip_flag = true
        forward_gz = params.output + "/forward.fastq"
        reverse_gz = params.output + "/reverse.fastq"

        println "Forward GZ: ${forward_gz}"
        println "Reverse GZ: ${reverse_gz}"

        // Forward lock

        forward_lock_path=  params.output+"_forward.lock"
        forward_lock_file = new File(forward_lock_path)

        // check if there is already a lock
        while (file(forward_lock_path).exists()) {
                // Wait for a specified time before checking again
                println "Waiting for ${forward_lock_path} to get unlocked."
                sleep(5000) // Wait for 30 seconds
        }

        println "Acquiring a new lock: ${forward_lock_path}"
        if (forward_lock_file.createNewFile()) {

            try {
                // Lock acquired, proceed with unzipping
                Channel
                    .fromPath(params.forward)
                    .splitFastq( by: 1000 )
                    .collectFile( name: forward_gz, sort: false)
                    .subscribe onComplete: { println "Done reading gzipped forward reads." }

            } finally {
                // Ensure that the lock file is deleted after unzipping
                forward_lock_file.delete()
                println "Deleting the acquired lock: ${forward_lock_path}"
            }
        }

        // Reverse lock

        reverse_lock_path=  params.output+"_reverse.lock"
        reverse_lock_file = new File(reverse_lock_path)

        // check if there is already a lock
        while (file(reverse_lock_path).exists()) {
                // Wait for a specified time before checking again
                println "Waiting for ${reverse_lock_path} to get unlocked."
                sleep(30000) // Wait for 30 seconds
        }
        println "Acquiring a new lock: ${reverse_lock_path}"
        if (reverse_lock_file.createNewFile()) {
            try {
                // Lock acquired, proceed with unzipping
                Channel
                    .fromPath(params.reverse)
                    .splitFastq( by: 1000 )
                    .collectFile( name: reverse_gz, sort: false)
                    .subscribe onComplete: { println "Done reading gzipped reverse reads." }

            } finally {
                // Ensure that the lock file is deleted after unzipping
                reverse_lock_file.delete()
                println "Deleting the acquired lock: ${reverse_lock_path}"
            }
        }

        reads += "--forward " + file(forward_gz) + " --reverse " + file(reverse_gz)

    }else{

        reads += "--forward " + file(params.forward) + " --reverse " + file(params.reverse)
    }

    paired = true
}
if (params.unpaired != ""){
    if(!file(params.unpaired).isFile()){
        println "ERROR: Incorrect filepath to unpaired reads!"
        usage(1)
    }
    String unpairedExtension = getExtension(params.unpaired)

    if(unpairedExtension == "gz"){
        gzip_flag = true
        // save gzip reads to file
        unpaired_gz = params.output + "/unpaired.fastq"

        Channel
            .fromPath(params.unpaired)
            .splitFastq( by: 1000 )
            .collectFile( name: unpaired_gz, sort: false)
            .subscribe onComplete: { println "Done reading gzipped unpaired reads." }

        reads += " --unpaired " + file(unpaired_gz)

    }else{

        reads += " --unpaired " + file(params.unpaired)
    }

    unpaired = true
}
if (paired == false && unpaired == false){
    println "ERROR: Incorrect combination of reads given!"
    usage(1)
}

// set necessary paths
WORKFLOWS = "$workflow.projectDir/pipeline"
WORKDIR = "$workflow.projectDir"
OUTPUT = file(params.output)
LOG = file("$OUTPUT/metacompass.log")
LOG.append("Parameters: \n$params\n\n")

if (!file("$OUTPUT/timeline/").mkdirs()){
    println "ERROR: Could not create directory for resource usage!"
    usage(1)
}

process initialize {
    cache 'lenient'

    output:
    stdout initialize

    shell:
    '''
    printf %s "-----Initializing MetaCompass-----\n$(date)\n\n" >> !{LOG}
    # for gzipped files, need to wait until fully unzipped by nextflow
    while [[ !{gzip_flag} == true ]]
    do

        if [[ "!{paired}" == true ]]; then
            if [[ -f "!{forward_gz}" && -f "!{reverse_gz}" ]]; then

                break
            fi
        elif [[ "!{unpaired}" == true ]]; then
            if [[ -f "!{unpaired_gz}" ]]; then
                break
            fi
        fi
        sleep 5

    done

    cmd="nextflow run !{WORKFLOWS}/initialize.nf    --reference_db !{params.reference_db} \
                                                    --output !{OUTPUT} \
                                                    -with-timeline !{OUTPUT}/timeline/initialize.html"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a !{LOG}
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finished Initialization successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: Initialization process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}

process ref_selection {
    cache 'lenient'

    input:
    val init from initialize

    output:
    env REFS into ref_selection

    shell:
    '''
    if [ !{params.skip_rs} == false ]
    then
        printf %s "-----Starting Reference Selection-----\n$(date)\n\n" >> !{LOG}

        cmd="nextflow run !{WORKFLOWS}/ref_selection.nf    !{reads} \
                                                        --reference_db !{params.reference_db} \
                                                        --filter_refs !{params.filter_refs} \
                                                        --ms !{params.ms} \
                                                        --clean !{params.clean} \
                                                        --match !{params.match} \
                                                        --readlen !{params.readlen} \
                                                        --depth_of_coverage !{params.depth_of_coverage} \
                                                        --breadth_of_coverage !{params.breadth_of_coverage} \
                                                        --percent_markers_covered !{params.percent_markers_covered} \
                                                        --threads !{params.threads} \
                                                        --output !{OUTPUT} \
                                                        --workdir !{WORKDIR} \
                                                        -with-timeline !{OUTPUT}/timeline/ref_selection.html"

        printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
        $cmd 2>&1 | tee -a !{LOG}
        exitcode=${PIPESTATUS[0]}

        if [[ $exitcode = "0" ]]; then
            printf %s "-----Finished Reference Selection successfully-----\n$(date)\n\n" >> !{LOG}
        else
            printf %s "-----ERROR: Reference Selection process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
            exit $exitcode
        fi

        REFS=!{OUTPUT}/reference_selection/cluster_refs/reference_candidates.txt

    else
        printf %s "-----Skipping Reference Selection-----\n$(date)\n\n" >> !{LOG}
        REFS=!{OUTPUT}/reference_selection/cluster_refs/reference_candidates.txt

    fi

    export REFS="!{OUTPUT}/reference_selection/cluster_refs/reference_candidates.txt"
    '''
}

process ref_culling {
    cache 'lenient'

    input:
    env REFS from ref_selection

    output:
    env REFS into culled_refs
    
    shell:
    '''
    reference_culling_dir="!{params.output}/reference_culling"
    mkdir -p \${reference_culling_dir}
    reference_culling_log="!{params.output}/reference_culling/reference_culling.log"
    touch \$reference_culling_log
    align_reads_out_dir="!{params.output}/reference_assembly"
    mkdir -p \${align_reads_out_dir}

    if [ !{params.skip_rc} == false ]
    then
        printf %s "-----Starting Reference Culling-----\n$(date)\n\n" >> !{LOG}

        num_refs=$(wc -l < $REFS)

        if [[ $num_refs > 0 ]]; then

            cmd="nextflow run !{WORKFLOWS}/ref_culling.nf !{reads} \
                                 --reference_db !{params.reference_db} \
                                 --ref_candidates ${REFS} \
                                 --ms !{params.ref_culling_ms} \
                                 --num_matches !{params.num_matches}  \
                                 --threads !{params.threads}  \
                                 --output !{params.output} \
                                 --workdir !{WORKDIR}  \
                                 --cull_stop !{params.stop} \
                                 -with-timeline !{OUTPUT}/timeline/ref_culling.html"

            printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
            $cmd 2>&1 | tee -a !{LOG}
            exitcode=${PIPESTATUS[0]}

            if [[ $exitcode = "0" ]]; then
                printf %s "-----Finished Reference Culling successfully-----\n$(date)\n\n" >> !{LOG}
            else
                printf %s "-----ERROR: Reference Culling process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
                exit $exitcode
            fi

        else

            printf %s "No references found. Skipping reference culling step and sending all reads to denovo assembly.\n" >> !{LOG}
            printf %s "-----Finished Reference Culling successfully-----\n$(date)\n\n" >> !{LOG}
        fi



    else
        printf %s "-----Skipping Reference Culling-----\n$(date)\n\n" >> !{LOG}
        MIN_REFS=!{OUTPUT}/reference_culling/cull_candidates/min_reference_candidates.txt
        ALL_REFS=!{OUTPUT}/reference_culling/collect_refs/ncbi_dataset/data/
        REF_INFO=!{OUTPUT}/reference_culling/collect_refs/ncbi_dataset/data/assembly_data_report.jsonl
    fi


    '''
}

process denovo {
    cache 'lenient'

    input:
    env REFS from culled_refs

    shell:
    '''
    if [ !{params.de_novo} == 0 ]; then
        printf %s "-----Skipping Denovo Assembly as requested-----\n$(date)\n\n" >> !{LOG}
        exit 0
    fi
    printf %s "-----Starting Denovo Assembly-----\n$(date)\n\n" >> !{LOG}
    
    denovo_dir="!{params.output}/denovo_assembly"
    mkdir -p ${denovo_dir}
    
    num_refs=$(wc -l < $REFS)
    
    if [[ $num_refs > 0 ]]; then
        printf %s "Found $num_refs reference(s). Processing unmapped reads for denovo assembly.\n" >> !{LOG}
    else
        printf %s "No references found. All reads will be used for denovo assembly.\n" >> !{LOG}
    fi

    cmd="nextflow run !{WORKFLOWS}/de_novo.nf \
        --workdir !{WORKDIR} \
        --output !{params.output} \
        --threads !{params.threads} \
        -with-timeline !{OUTPUT}/timeline/denovo.html"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a !{LOG}
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finished Denovo Assembly successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: Denovo Assembly process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}

// process report_generation {
//     cache 'lenient'
//
//     input:
//     env ALL_REFS from report_gen
//     env MIN_REFS from ref_culling_min_refs
//     env REF_INFO from ref_culling_ref_info
//     env ALIGNED_READS from ref_assembly_align
//     env UNALIGNED_READS from ref_assembly_unalign
//     env REFGUIDED_CONTIGS from ref_assembly
//     // env DENOVO_CONTIGS from de_novo_assembly
//     // env MERGED_CONTIGS from assembly_merge
//
//     output:
//     env REFGUIDED_CONTIGS into report_generation
//
//     shell:
//     '''
//     printf %s "-----Starting Report Generation-----\n$(date)\n\n" >> !{LOG}
//     printf %s "-----Checking-----\n$REFGUIDED_CONTIGS\n\n" >> !{LOG}
//
//     if [[ -f "!{forward_gz}" && -f "!{reverse_gz}" ]]; then
//         rm !{forward_gz}
//         rm !{reverse_gz}
//     fi
//     if [[ -f "!{unpaired_gz}" ]]; then
//         rm !{unpaired_gz}
//     fi
//     printf %s "-----Checking-----\n$REFGUIDED_CONTIGS\n\n" >> !{LOG}
//
//     cmd="python3 !{WORKDIR}/scripts/output_folder_creation.py -o !{params.output} -mr $MIN_REFS -rc $REFGUIDED_CONTIGS -ar $ALIGNED_READS -ur $UNALIGNED_READS"
//     printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
//     $cmd 2>&1 | tee -a !{LOG}
//
//     if [ !{params.clean_uf} != false ]
//     then
//             cmd="python3 !{WORKDIR}/scripts/clean_output.py -o !{params.output} -pc"
//             printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
//             $cmd 2>&1 | tee -a !{LOG}
//     fi
//     '''
// }


// report_generation
//     .subscribe onComplete{ file(forward_gz).delete(); file(reverse_gz).delete(); file(unpaired_gz).delete() }
