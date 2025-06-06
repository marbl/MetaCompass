#!/usr/bin/env nextflow

// required
reads = ""
params.forward = ""
params.reverse = ""
params.unpaired = ""
paired = false
unpaired = false
params.reference_db = ""
params.ref_filter = true
params.filter_kmer_size = ""
params.culling_kmer_size = ""
params.depth_of_cov = ""
params.breadth_of_cov = ""
params.references = ""
params.references_info = ""
params.min_ref_list = ""
params.aligned_reads = ""
params.ref_contigs = ""
params.denovo_contigs = ""
params.merged_contigs = ""
params.frc_path = ""
params.libsize_avg = 165
params.libsize_sdv = 92.8
params.orientation_a = "+"
params.orientation_b = "-"
params.mp_size_z_cutoff = 2.58
params.mp_edge_z_cutoff = 2.58
params.doc_logp_cutoff = 7
params.output = ""
params.help = false

def usage(status) {
    log.info "Usage: \nnextflow run report_generation.nf \n\
            [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] | \n\
            [--forward /path/to/forwardReads --reverse /phat/to/reverseReads] | \n\
            [--unpaired /path/to/unpairedReads]] \n\
            --reference_db /path/to/reference_db --ref_filter [True | False] --filter_kmer_size n \n\
            --culling_kmer_size c --depth_of_cov d --breadth_of_cov b --references /path/to/references \n\
            --references_info assembly_data_report.jsonl --aligned_reads /path/to/align_reads --min_ref_list min_reference_candidates.txt \n\
            --ref_contigs /path/to/pilon --denovo_contigs contigs.fa  --merged_contigs contigs.merge.fasta \n\
            --output /path/to/outputDir"
    
    log.info ""
    log.info "Required:"
    log.info ""

    log.info " --forward                Path to forward paired-end read."
    log.info " --reverse                Path to reverse paried-end read."
    log.info " --unpaired               Path to unpaired read fasta file(s)."
    log.info " --reference_db           Path to reference marker gene database."
    log.info " --ref_filter             Whether reference filtering was use. True/False."
    log.info " --filter_kmer_size       Kmer size used for reference filtering."
    log.info " --culling_kmer_size      Kmer size used for culling references."
    log.info " --depth_of_cov           Depth of coverage cutoff for clustering references."
    log.info " --breadth_of_cov         Breadth of coverage cutoff for clustering references."
    log.info " --references             Path to downloaded references file. Directory structure should be same as 'datasets download assembly' command."
    log.info " --references_info        Path to 'assembly_data_report.jsonl' file downloaded in Reference-Selection step with the datasets tool."
    log.info " --run_ref_info           Whether or not to run reference_info step."
    log.info " --min_ref_list           List of references used in Reference-Guided Assembly."
    log.info " --aligned_reads          Path to align_reads output folder from reference culling process."
    log.info " --ref_contigs            Path to folder containing folders named by each reference used, which contain contig files for that particular assembly. From reference-guided assembly."
    log.info " --denovo_contigs         Path to contigs file generated from Denovo Assembly."
    log.info " --merged_conigs          Path to contigs file generated from Assembly Merge."
    log.info " --run_frc                Whether or not to run FRC curve generation step."
    log.info " --frc_path               Path to FRC executable."
    log.info " --libsize_avg            Average library size of the reads."
    log.info " --libsize_sdv            Standard deviation of library size of reads."
    log.info " --orientation_a          Orientation of forward reads."
    log.info " --orientation_b          Orientation of reverse reads."
    log.info " --mp_size_z_cutoff       Z score for mate pair size cutoff."
    log.info " --mp_edge_z_cutoff       Z score for mate pair edge curoff."
    log.info " --doc_logp_cutoff        Log(p) cutoff for depth of coverage."
    log.info " --output                 Path to output folder."
    
    log.info ""
    log.info "Optional:"
    log.info ""
    

    log.info " --workdir                All paths will be relative to this directory."
    log.info " --help                   Print help message."

    exit status
}
println "Params: $params"
if (params.help){
    usage(0)
}

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

if ( ! (params.reference_db != "" && file(params.reference_db, isDirectory: true)) ){
    println "ERROR: Reference db path does not exist!"
    usage(1)
}

if ( params.filter_kmer_size == ""){
    println "ERROR: filter_kmer_size not given!"
    usage(1)
}

if ( params.culling_kmer_size == ""){
    println "ERROR: culling_kmer_size not given!"
    usage(1)
}

if ( params.depth_of_cov == ""){
    println "ERROR: depth_of_cov not given!"
    usage(1)
}

if ( params.breadth_of_cov == ""){
    println "ERROR: breadth_of_cov not given!"
    usage(1)
}

if ( ! (params.references != "" && file(params.references, isDirectory: true)) ){
    println "ERROR: references path not valid directory!"
    usage(1)
}

if ( ! (params.references_info != "" && file(params.references_info, checkIfExists: true)) ){
    println "ERROR: references_info file does not regular file!"
    usage(1)
}

if ( ! (params.min_ref_list != "" && file(params.min_ref_list, checkIfExists: true)) ){
    println "ERROR: min_ref_list is not a regular file!"
    usage(1)
}

if ( ! (params.aligned_reads != "" && file(params.aligned_reads, isDirectory: true)) ){
    println "ERROR: aligned_reads is not a valid directory!"
    usage(1)
}

if ( ! (params.ref_contigs != "" && file(params.ref_contigs, isDirectory: true)) ){
    println "ERROR: ref_contigs path not valid directory!"
    usage(1)
}

if ( ! (params.denovo_contigs != "" && file(params.denovo_contigs, checkIfExists: true)) ) {
    println "ERROR: denovo_contigs is not a regular file!"
    usage(1)
}

if ( ! (params.merged_contigs != "" && file(params.merged_contigs, checkIfExists: true)) ) {
    println "ERROR: merged_contigs is not a regular file!"
    usage(1)
}

if ( ! (params.frc_path != "" && file(params.frc_path, checkIfExists: true)) ) {
    println "ERROR: frc_path is not valid path!"
    usage(1)
}

if ( ! (params.output != "" && file(params.output, checkIfExists: true)) ) {
    println "ERROR: Output folder does not exist!"
    usage(1)
}

if( ! file("$params.output/report_generation").mkdir()){
    println "ERROR: Cannot create output directory for report generation!"
    usage(1)
}

WORKDIR = "$params.workdir"
OUTPUT = "$params.output/report_generation"
LOG = file("$params.output/report_generation/report_generation.log")
LOG.append("Parameters: \n$params\n\n")

process reference_info {
    cache 'lenient'

    output:
    env REF_INFO into reference_info

    shell:
    '''
    printf %s "-----Running Reference Info-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/ref_info"
    mkdir $out
    log=!{OUTPUT}"/ref_info.log"

    cmd="bash !{WORKDIR}/scripts/ref_info.sh    !{reads} \
                                                --references !{params.references} \
                                                --depth !{params.depth_of_cov} \
                                                --run !{params.run_ref_info} \
                                                --output $out \
                                                --log $log"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a $log
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running Reference Info successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: Reference Info process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi

    REF_INFO=$out/ref_metrics.tsv
    '''
}

process frc_data_gen {
    cache 'lenient'

    output:
    env REFGUIDED_FRC_INFO into ref_frc_info
    env DENOVO_FRC_INFO into denovo_frc_info
    env MERGED_FRC_INFO into merged_frc_info

    shell:
    '''
    printf %s "-----Running FRC Data Generation-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/frc_data_gen"
    mkdir $out
    log=!{OUTPUT}"/frc_data_gen.log"

    if [[ !{params.run_frc} == "true" ]]; then

        cmd="bash !{WORKDIR}/scripts/frc_data_gen.sh    !{reads} \
                                                        --split_reads !{params.aligned_reads} \
                                                        --refguided_contigs !{params.ref_contigs} \
                                                        --denovo_contigs !{params.denovo_contigs} \
                                                        --merged_contigs !{params.merged_contigs} \
                                                        --frc !{params.frc_path} \
                                                        --libsize_avg !{params.libsize_avg} \
                                                        --libsize_sdv !{params.libsize_sdv} \
                                                        --orientation_a !{params.orientation_a} \
                                                        --orientation_b !{params.orientation_b} \
                                                        --mp_size_z_cutoff !{params.mp_size_z_cutoff} \
                                                        --mp_edge_z_cutoff !{params.mp_edge_z_cutoff} \
                                                        --doc_logp_cutoff !{params.doc_logp_cutoff} \
                                                        --out $out \
                                                        --log $log"

        printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
        $cmd 2>&1 | tee -a $log
        exitcode=${PIPESTATUS[0]}

        if [[ $exitcode = "0" ]]; then
            printf %s "-----Finish running FRC Data Generation successfully-----\n$(date)\n\n" >> !{LOG}
        else
            printf %s "-----ERROR: FRC Data Generation process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
            exit $exitcode
        fi
    else

        mkdir -p $out/ref_guided
        mkdir -p $out/denovo
        mkdir -p $out/merged

        touch $out/denovo/denovo.tsv
        touch $out/merged/merged.tsv

        printf %s "-----FRC Data Generation will be skipped-----\n$(date)\n\n" >> !{LOG}
        printf %s "-----Finish running FRC Data Generation successfully-----\n$(date)\n\n" >> !{LOG}
    fi

    REFGUIDED_FRC_INFO=$out/ref_guided
    DENOVO_FRC_INFO=$out/denovo/denovo.tsv
    MERGED_FRC_INFO=$out/merged/merged.tsv
    '''
}

process frc_curve_gen {
    cache 'lenient'

    input:
    env REFGUIDED_FRC_INFO from ref_frc_info
    env DENOVO_FRC_INFO from denovo_frc_info
    env MERGED_FRC_INFO from merged_frc_info

    output:
    env FRC_CURVES into frc_curves

    shell:
    '''
    printf %s "-----Running FRC Curve Generation-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/frc_curve_gen"
    mkdir $out
    log=!{OUTPUT}"/frc_curve_gen.log"

    if [[ !{params.run_frc} == "true" ]]; then

        frc_info=""

        if [[ -d $REFGUIDED_FRC_INFO ]]; then

            frc_info="--ref_info ${REFGUIDED_FRC_INFO}"

            ref_out="${out}/refguided"
            mkdir $ref_out
            for i in $(ls $REFGUIDED_FRC_INFO); do
                curr_ref_out="${ref_out}/${i}"
                mkdir $curr_ref_out
            done
        fi
        if [[ -f $DENOVO_FRC_INFO ]]; then

            denovo_out="${out}/denovo"
            mkdir $denovo_out
            frc_info="${frc_info} --denovo_info ${DENOVO_FRC_INFO}"
        fi
        if [[ -f $MERGED_FRC_INFO ]]; then

            merged_out="${out}/merged"
            mkdir $merged_out
            frc_info="${frc_info} --merged_info ${MERGED_FRC_INFO}"
        fi

        cmd="python !{WORKDIR}/scripts/frc_plot.py    $frc_info \
                                                    --out $out \
                                                    --log $log"

        printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
        $cmd 2>&1 | tee -a $log
        exitcode=${PIPESTATUS[0]}

        if [[ $exitcode = "0" ]]; then
            printf %s "-----Finish running FRC Curve Generation successfully-----\n$(date)\n\n" >> !{LOG}
        else
            printf %s "-----ERROR: FRC Curve Generation process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
            exit $exitcode
        fi
    else
        printf %s "-----FRC Curve Generation will be skipped-----\n$(date)\n\n" >> !{LOG}
        printf %s "-----Finish running FRC Curve Generation successfully-----\n$(date)\n\n" >> !{LOG}
    fi
    
    FRC_CURVES=$out
    '''
}

process report_metrics {
    cache 'lenient'

    input:
    env REF_INFO from reference_info
    env FRC_CURVES from frc_curves

    output:
    stdout report_metrics

    shell:
    '''
    printf %s "-----Running Report Metrics-----\n$(date)\n\n" >> !{LOG}

    out="!{OUTPUT}/report_metrics"
    mkdir $out
    log=$out"/report_metrics.log"
    out_json=$out"/data.json"
    report_template=!{WORKDIR}/scripts/report

    # pre-process number of reads
    if [[ -f "!{params.forward}" ]]; then
        tmp_num=$(wc -l < !{params.forward})
        n_forward=$(($tmp_num/4))
        forward=!{params.forward}
    else
        n_forward="NA"
        forward="NA"
    fi

    if [[ -f "!{params.reverse}" ]]; then
        tmp_num=$(wc -l < !{params.reverse})
        n_reverse=$(($tmp_num/4))
        reverse=!{params.reverse}
    else
        n_reverse="NA"
        reverse="NA"
    fi

    if [[ -f "!{params.unpaired}" ]]; then
        tmp_num=$(wc -l < !{params.unpaired})
        n_unpaired=$(($tmp_num/4))
        unpaired=!{params.unpaired}
    else
        n_unpaired="NA"
        unpaired="NA"
    fi

    cmd="python !{WORKDIR}/scripts/report_generation.py --forward $forward \
                                                        --reverse $reverse \
                                                        --unpaired $unpaired \
                                                        --forward_reads $n_forward \
                                                        --reverse_reads $n_reverse \
                                                        --unpaired_reads $n_unpaired \
                                                        --ref_db !{params.reference_db} \
                                                        --ref_filter !{params.ref_filter} \
                                                        --filter_kmer_size !{params.filter_kmer_size} \
                                                        --culling_kmer_size !{params.culling_kmer_size} \
                                                        --depth_of_cov !{params.depth_of_cov} \
                                                        --breadth_of_cov !{params.breadth_of_cov} \
                                                        --references_info !{params.references_info} \
                                                        --references_metrics $REF_INFO \
                                                        --min_ref_list !{params.min_ref_list} \
                                                        --ref_contigs !{params.ref_contigs} \
                                                        --denovo_contigs !{params.denovo_contigs} \
                                                        --merged_contigs !{params.merged_contigs} \
                                                        --frc_curves $FRC_CURVES \
                                                        --templates !{WORKDIR}/scripts/report \
                                                        --out_json $out_json \
                                                        --out !{OUTPUT}"

    printf "Command Executed: \n\n$cmd\n\n" >> !{LOG}
    $cmd 2>&1 | tee -a $log
    exitcode=${PIPESTATUS[0]}

    if [[ $exitcode = "0" ]]; then
        printf %s "-----Finish running Report Metrics successfully-----\n$(date)\n\n" >> !{LOG}
    else
        printf %s "-----ERROR: Report Metrics process exit with exit code $exitcode-----\n$(date)\n\n" >> !{LOG}
        exit $exitcode
    fi
    '''
}