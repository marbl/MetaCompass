#!/usr/bin/env nextflow

// -----------------------------------------------------------------------------
// Default parameters
// -----------------------------------------------------------------------------
params.output    = params.output ?: ''
params.ref_cache = params.ref_cache ?: params.output

// -----------------------------------------------------------------------------
// Helper values
// -----------------------------------------------------------------------------
final String REFERENCE_CULLING_DIR = params.output ?
    "${params.output}/reference_culling" :
    "reference_culling"

// -----------------------------------------------------------------------------
// Download or locate reference genomes for candidate references
// MOUMI: update: added ref_cache 
// -----------------------------------------------------------------------------
process collect_refs {
    publishDir "${REFERENCE_CULLING_DIR}", mode: 'copy'

    input:
    path ref_candidates

    output:
    path "matching_files.txt"

    script:
    """
    bash "${projectDir}/scripts/download.sh" \
        T \
        "${params.reference_db}" \
        "${params.ref_cache}" \
        "download.log" \
        "${ref_candidates}" \
        .

    if [[ \$? -ne 0 ]]; then
        echo "Error downloading references."
        exit 1
    fi

    # Initialize an array to store matching FASTA files
    matching_files_array=()

    # Read candidate reference IDs into an array
    ref_candidates_array=()
    while read -r line; do
        ref_candidates_array+=("\${line}")
    done < "${ref_candidates}"

    # Search for matching files in both:
    # 1) shared reference DB
    # 2) local writable ref_cache
    for fl in "\${ref_candidates_array[@]}"; do
        fl=\$(echo "\$fl" | tr -d '[:space:]')
        [[ -z "\$fl" ]] && continue

        for root in \
            "${params.reference_db}/collect_refs/ncbi_dataset/data" \
            "${params.ref_cache}/collect_refs/ncbi_dataset/data"
        do
            if [[ -d "\${root}/\${fl}" ]]; then
                files=(\$(find "\${root}/\${fl}" -type f -name "GCA_*_genomic.fna"))
                for file_path in "\${files[@]}"; do
                    matching_files_array+=("\${file_path}")
                done
            fi
        done
    done

    # Remove duplicates and write results
    printf "%s\n" "\${matching_files_array[@]}" | awk '!seen[\$0]++' > matching_files.txt
    """
}

// -----------------------------------------------------------------------------
// Run skani triangle on collected reference FASTA files
// -----------------------------------------------------------------------------
process SkaniTriangle {
    publishDir "${REFERENCE_CULLING_DIR}", mode: 'copy'

    input:
    path matching_files
    val ref_candidates_count

    output:
    path "skani_matrix_stool.txt"

    script:
    """
    if [[ ${ref_candidates_count} -gt 1 ]]; then
        skani triangle -l "${matching_files}" -t ${params.threads} > skani_matrix_stool.txt
    else
        touch skani_matrix_stool.txt
    fi
    """
}

// -----------------------------------------------------------------------------
// Cluster references based on skani matrix
// -----------------------------------------------------------------------------
process Cluster {
    publishDir "${REFERENCE_CULLING_DIR}", mode: 'copy'

    input:
    path skani_matrix

    output:
    path "clusters_name.txt", emit: clusters_name
    path "clusters.csv", emit: clusters

    script:
    """
    python3 "${projectDir}/scripts/cluster.py" "${skani_matrix}" 5 .
    """
}

// -----------------------------------------------------------------------------
// Concatenate FASTA files within each cluster
// -----------------------------------------------------------------------------
process ConcatFasta {
    publishDir "${REFERENCE_CULLING_DIR}", mode: 'copy'

    input:
    path clusters

    output:
    path "cluster_content.csv"
    path "cluster_*.fna", emit: cluster_files

    script:
    """
    bash "${projectDir}/scripts/concat_fasta.sh" . "${clusters}"
    """
}

// -----------------------------------------------------------------------------
// Build k-mer index from input reads
// -----------------------------------------------------------------------------
process IndexReads {
    publishDir "${REFERENCE_CULLING_DIR}", mode: 'copy'

    input:
    val readFiles

    output:
    path "reads.base_2.kmers.*"

    script:
    """
    echo "${readFiles}" > inputs.txt
    kmc -k28 -ci1 -hp -fq @inputs.txt reads.base_2.kmers .
    """
}

// -----------------------------------------------------------------------------
// Build k-mer index for each reference cluster FASTA
// -----------------------------------------------------------------------------
process ClusterIndex {
    publishDir "${REFERENCE_CULLING_DIR}", mode: 'copy'

    input:
    path cluster_file

    output:
    path "*.kmc_*"

    script:
    """
    name=\$(basename "${cluster_file}" .fna)
    mkdir "d_\${name}"
    kmc -k28 -ci -hp -fm "${cluster_file}" "\${name}" "d_\${name}"
    """
}