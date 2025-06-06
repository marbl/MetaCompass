#!/bin/bash

# Get the directory path from command line arguments
CLUSTER_DIR="${1}"
CLUSTERS="${2}"
# Move to the directory where the fna files are

# counter variable for naming files
counter=1

while IFS= read -r line
do
    # remove comma from the line and replace it with space to create a space-separated list
    files=$(echo $line| tr ',' ' ')
    echo processing... ${files}
    cluster_name="cluster_${counter}"
    # concatenate files and create a new fasta file for the cluster
    cat $files > "${CLUSTER_DIR}/cluster_${counter}.fna"
    cluster_content=${line}
    echo "${cluster_name},${cluster_content}" >> "${CLUSTER_DIR}/cluster_content.csv"
    ((counter++))
done < ${CLUSTERS}

