#!/usr/bin/env nextflow

// required parameters
reads = ""
params.forward = ""
params.reverse = ""
params.unpaired = ""
params.ref_candidates = ""
params.ref_genomes = ""
params.output = ""
params.threads = 16
params.help = false
python_cull_input = ""
params.reference_db = ""  // Define a default value or leave it empty

// Workflow variables
REFS = ""
WORKDIR = "$params.workdir"

// Project folders
BIN = "$workflow.projectDir/../bin"
PIPELINE="$workflow.projectDir/../pipeline"
SCRIPTS = "$workflow.projectDir/../scripts"

// Output variables
OUTPUT = "$params.output/reference_culling"
CULLING_OUTPUT_PATH="${params.output}/reference_culling"


def usage(status) {
    log.info "Usage: \nnextflow run ref_culling.nf \n\
            [[--forward /path/to/forwardReads --reverse /path/to/reverseReads --unpaired /path/to/unpairedReads] | \n\
            [--forward /path/to/forwardReads --reverse /phat/to/reverseReads] | \n\
            [--unpaired /path/to/unpairedReads]] \n\
            --ref_candidates /path/to/reference_candidates.txt \n\
            --ref_genomes /path/to/reference_genome_files \n\
            --ms kmer size \n\
            --num_matches \n\
            --output /path/to/outputDir"

    log.info ""
    log.info "Required:"
    log.info ""

    log.info " --forward            Path to forward read."
    log.info " --reverse            Path to reverse read."
    log.info " --unpaired           Paths to unpaired read."
    log.info " --ref_candidates || --ref_genomes    path to a txt file containing reference candidate accession ids (Reference files will be downloaded if this argument is given). Or a path to reference genomes in fasta format (directory must contain only relevant files)."
    log.info " --ms                 Kmer size used when collecting kmers for comparison from reads and references."
    log.info " --num_matches        Number matches between read kmers and ref kmers. Threshold for assigning reads to a reference (Default: 10)."
    log.info " --output             Path to output folder."
    log.info " --log                Path to log file."
    log.info " --workdir            Path to working directory. (MetaCompass Folder)"

    log.info ""
    log.info "Optional:"
    log.info ""

    log.info " --ref_genomes        Path to reference genome files. If genomes are not available, they will be downloaded from ncbi."
    log.info " --threads            Number of threads to use. Default set to 12."
    log.info " --help               Print help message."

    exit status
}
println "The Reference DB is set to: ${params.reference_db}"
println "Output: ${params.output}"
println "Workdir: ${params.workdir}"
if (params.help){
    usage(0)
}

// TODO: Validate required set of reads given
if( params.forward != "" && params.reverse != "" ){
    reads += "--forward $params.forward --reverse $params.reverse "
    python_cull_input += "$params.forward,$params.reverse"
}

if( params.unpaired != "" ){
    reads += "--unpaired $params.unpaired"
    python_cull_input += "$params.unpaired"
}

if( !file(params.output).isDirectory() || \
    !file(params.workdir).isDirectory() ){

    println "ERROR: Please check required inputs."
    usage(1)
}


LOG = file("$params.output/reference_culling/reference_culling.log")
LOG.append("Parameters: \n$params\n\n")

println "Debug: params.reference_db is ${params.reference_db}"
if (params.reference_db) {
    params.downloadedList = "${params.reference_db}/downloaded.txt"

    println "The Reference DB is set to: ${params.reference_db}/downloaded.txt"
    println "Debug: params.downloadedList is ${params.downloadedList}"
} else {
    println "Warning: params.reference_db is not set"
}
if( params.ref_genomes != "" && file(params.ref_genomes, checkIfExists: true)){
    REFS = "$params.ref_genomes"
}

// Channel to count the number of lines in ref_candidates
ref_candidates_count_ch = Channel.fromPath(params.ref_candidates).splitText().count()

/*
 * This process collects references, searches for matching files in the
 * downloaded directories based on reference candidates, and then outputs
 * the paths of these matching files to 'matching_files.txt'.
 */
process collect_refs{
    cache 'lenient'
    publishDir "${CULLING_OUTPUT_PATH}", mode: 'copy'

    output:
    file('matching_files.txt') into matching_files_ch

    shell:
    '''
    #bash "!{SCRIPTS}/download.sh" T "!{params.reference_db}" "!{params.output}/reference_culling/downloaded.txt" "!{LOG}" "!{params.ref_candidates}" "!{params.output}/reference_culling"
    bash "!{SCRIPTS}/download.sh" T "!{params.reference_db}" "!{params.downloadedList}" "!{LOG}" "!{params.ref_candidates}" "!{params.output}/reference_culling"
    if [[ $? -ne 0 ]]; then
        echo "Error downloading references."
        exit 1
    fi

    # Initialize an empty array to store matching files
    matching_files_array=()

    # Read ref_candidates file line-by-line into an array
    ref_candidates_array=()
    while read -r line; do
        ref_candidates_array+=("\${line}")
    done < !{params.ref_candidates}

    # Search for matching files in the downloaded directories
    for x in "\${ref_candidates_array[@]}"; do
        # Assuming that files are downloaded to a path like downloaded_path/${x}/
        #files=($(find !{params.output}/reference_culling/ncbi_dataset/data/\${x}/ -type f -name "GCA_*_genomic.fna"))
        files=($(find !{params.reference_db}/collect_refs/ncbi_dataset/data/\${x}/ -type f -name "GCA_*_genomic.fna"))
        # Add them to the matching_files_array
        for file_path in "\${files[@]}"; do
            matching_files_array+=("\${file_path}")
        done
    done

    # Write the matching files to a text file
    printf "%s\n" "\${matching_files_array[@]}" > matching_files.txt
    '''
}

/*
 * This process utilizes the 'Skani Triangle' tool to analyze the files
 * from 'matching_files.txt' and outputs the results to 'skani_matrix_stool.txt'.
 */
process SkaniTriangle {
    executor 'local'
    publishDir "${CULLING_OUTPUT_PATH}", mode: 'copy'
    input:
    file('matching_files.txt') from matching_files_ch
    input:
    val(ref_candidates_count) from ref_candidates_count_ch

    output:
    file("skani_matrix_stool.txt") into skani_output_ch

    script:
    """
    # Read matching_files.txt into an array
    readarray -t matching_files < matching_files.txt

    if [ ${ref_candidates_count} -gt 1 ]; then
        # Run Skani Triangle
        skani triangle -l matching_files.txt -t ${params.threads} > skani_matrix_stool.txt
    else
        # If ref_candidates_count is not greater than 1, create an empty skani_matrix_stool.txt
        touch skani_matrix_stool.txt
    fi
    """
}

/*
 * This process runs a Python script to perform clustering based on the
 * 'skani_matrix_stool.txt' input. It produces two files: 'clusters_name.txt'
 * and 'clusters.txt' that represent the clustering results.
 */
process Cluster {

    publishDir "${CULLING_OUTPUT_PATH}", mode: 'copy'
    input:
    file("skani_matrix_stool.txt") from skani_output_ch

    output:
    file("clusters_name.txt") into clusters_name_ch
    file("clusters.txt") into clusters_ch
    script:
    """
    python3 "${SCRIPTS}/cluster.py"  skani_matrix_stool.txt 5 "${CULLING_OUTPUT_PATH}"
    ln -s "${CULLING_OUTPUT_PATH}/clusters_name.txt" clusters_name.txt
    ln -s "${CULLING_OUTPUT_PATH}/clusters.txt" clusters.txt
    """
}

/*
 * This process concatenates FASTA files from clusters. The process generates
 * a 'cluster_content.csv' file, providing a summary of the concatenated clusters.
 */
process ConcatFasta {
    input:
    file("clusters_name.txt") from clusters_name_ch
    file("clusters.txt") from clusters_ch
    output:
    file("cluster_content.csv") into concat_done_ch
    val("trigger") into trigger_ch  // Trigger signal
    script:
    """
    sh "${SCRIPTS}/concat_fasta.sh" "${CULLING_OUTPUT_PATH}"
    ln -s "${CULLING_OUTPUT_PATH}/cluster_content.csv" cluster_content.csv
    """

}

println "params.forward: ${params.forward}"

/*
 * This process uses the BBmap tool to reformat input read files. If paired-end
 * reads are provided, they are interleaved; if single-end reads are provided,
 * they are simply copied. The output is the 'interleaved_reads.fq' file.
 * Additionally, the process runs the KMC tool for k-mer counting.
 */

process run_BBmap {
    publishDir "${CULLING_OUTPUT_PATH}", mode: 'copy'

    output:
    file("interleaved_reads.fq") into interleaved_ch

    script:
    """
    if [ ! -z "${params.forward}" ] && [ ! -z "${params.reverse}" ]; then

        "${BIN}/bbmap/reformat.sh" in1=${params.forward} in2=${params.reverse} out=interleaved_reads.fq
   #     stat -c "%s" interleaved_reads.fq >/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/kmc/scripts/a.txt
    elif [ ! -z "${params.unpaired}" ]; then

        cp ${params.unpaired}  interleaved_reads.fq
    else
        echo "No input files provided."
        exit 1
    fi
    "${BIN}/kmc/kmc" -k28 -ci1 -hp -fq interleaved_reads.fq  "${CULLING_OUTPUT_PATH}/reads.base_2.kmers" "${CULLING_OUTPUT_PATH}"
    """
}

Channel.create().set { startAlignReadsChannel }

/*
 * This process executes an external Nextflow script 'parallized.nf',
 * which performs tasks related to the KMC tool (k-mer counting).
 * The process serves as a bridge or trigger for parallelized tasks in
 * the external script.
 */
process run_KMC {
    executor 'local'
    input:
    file interleaved_file from interleaved_ch
    val trigger from trigger_ch

    output:
    val true into startAlignReadsChannel

    script:
    """
    nextflow run "${PIPELINE}/parallized.nf" --output "${CULLING_OUTPUT_PATH}" --threads ${params.threads}
    """

}

/*
 * This process runs the 'align_reads' Python module. It aligns interleaved reads
 * with reference genome sketches, producing an alignment output. The module
 * operates based on cluster lists, read sketches, reference genome sketches,
 * and other provided parameters. Results are logged to 'process.log'.
 */
process run_align_reads {
    executor 'local'
    input:
    val canStart from startAlignReadsChannel

    script:
    """
    output_path=${params.output}
    align_reads_output="\${output_path}/reference_assembly"

    mkdir -p \${align_reads_output}

    cluster_list="${CULLING_OUTPUT_PATH}/clusters.txt"
    all_reads_sketches="${CULLING_OUTPUT_PATH}/reads.base_2.kmers"
    ref_genome_sketches="${CULLING_OUTPUT_PATH}"
    interleaved_reads="${CULLING_OUTPUT_PATH}/interleaved_reads.fq"
    metacompass_scripts="${SCRIPTS}"
    log_file="\${align_reads_output}/process.log"


    cd ${SCRIPTS}

    python -u -m align_reads \
        -ir \${interleaved_reads} \
        -cl \${cluster_list} \
        -rs \${ref_genome_sketches} \
        -as \${all_reads_sketches}  \
        -o "\${align_reads_output}" \
        -debug \
        -mcl 2000 \
        -t ${params.threads}  \
        > \${log_file}

    cd ..
    """
}
