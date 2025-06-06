#!/usr/bin/env nextflow

process collect_refs {
    publishDir {
        path: file("$params.output/reference_culling"), 
        mode: 'copy'
    } 	

    input: 
        path ref_candidates	

    output: 
        path "matching_files.txt"

    script:
    """
   bash "${projectDir}/scripts/download.sh" T "${params.reference_db}" "${params.reference_db}/downloaded.txt" "download.log" "${ref_candidates}" .
    if [[ \$? -ne 0 ]]; then
        echo "Error downloading references."
        exit 1
    fi

    # Initialize an empty array to store matching files
    matching_files_array=()

    # Read ref_candidates file line-by-line into an array
    ref_candidates_array=()
    while read -r line; do
        ref_candidates_array+=("\${line}")
    done < ${ref_candidates}

    # Search for matching files in the downloaded directories
    for fl in "\${ref_candidates_array[@]}"; do
        # Assuming that files are downloaded to a path like downloaded_path/\${fl}/
        files=(\$(find ${params.reference_db}/collect_refs/ncbi_dataset/data/\${fl}/ -type f -name "GCA_*_genomic.fna"))
        # Add them to the matching_files_array
        for file_path in "\${files[@]}"; do
            matching_files_array+=("\${file_path}")
        done
    done

    # Write the matching files to a text file
    printf "%s\n" "\${matching_files_array[@]}" > matching_files.txt
    """
}

/*
 * This process utilizes the 'Skani Triangle' tool to analyze the files
 * from 'matching_files.txt' and outputs the results to 'skani_matrix_stool.txt'
.
 */
process SkaniTriangle { 
    publishDir {
        path: file("$params.output/reference_culling"), 
        mode: 'copy'
    } 	

    input: 
    path matching_files
    val ref_candidates_count

    output:
    path "skani_matrix_stool.txt"

    script:
    """
    # Read matching_files.txt into an array
    readarray -t matching_files < $matching_files

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
    publishDir {
        path: file("$params.output/reference_culling"), 
        mode: 'copy'
    } 	

    input:
    path skani_matrix

    output:
    path "clusters_name.txt", emit: clusters_name
    path "clusters.csv", emit: clusters

    script:
    """
    python3 "${projectDir}/scripts/cluster.py"  skani_matrix_stool.txt 5 .
    """
}

/*
 * This process concatenates FASTA files from clusters. The process generates
 * a 'cluster_content.csv' file, providing a summary of the concatenated clusters.
 */
process ConcatFasta {
    publishDir {
        path: file("$params.output/reference_culling"), 
        mode: 'copy'
    } 	

    input:
    path clusters

    output:
    path "cluster_content.csv"
    path "cluster_*.fna", emit: cluster_files
    //val("trigger") into trigger_ch  // Trigger signal

    script:
    """
    bash "${projectDir}/scripts/concat_fasta.sh" . ${clusters}
    """
}

/*
* This process creates a k-mer index from the input reads
*/
process IndexReads {
    publishDir {
        path: file("$params.output/reference_culling"), 
        mode: 'copy'
    } 	

    input:
    val readFiles 

    output:
    path "reads.base_2.kmers.*"

    script:
    """
    echo "$readFiles" > inputs.txt
    kmc -k28 -ci1 -hp -fq @inputs.txt reads.base_2.kmers .
    """
}

/*
* This process computes the k-mer index of each cluster
*/
process ClusterIndex {
    publishDir {
        path: file("$params.output/reference_culling"), 
        mode: 'copy'
    } 	

    input:
    path cluster_file

    output:
    path "*.kmc_*"

    script:
    """
    name=\$(basename ${cluster_file} .fna)
    mkdir d_\${name}
    kmc -k28 -ci -hp -fm ${cluster_file} \${name} d_\${name}
    
    """
}
