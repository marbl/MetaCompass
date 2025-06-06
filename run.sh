#!/bin/bash

#forward_read="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/hmp_samples/posterior_fornix/filtered_files/SRR628270_clean_R1.fastq.gz"
#reverse_read="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/hmp_samples/posterior_fornix/filtered_files/SRR628270_clean_R2.fastq.gz"

forward_read='/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/hmp_samples/tongue_dorsum/filtered_files/SRR514250_clean_R1.fastq.gz'
reverse_read='/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/hmp_samples/tongue_dorsum/filtered_files/SRR514250_clean_R2.fastq.gz'

output_folder="/fs/cbcb-scratch/ayyangar/metacompass_new_version/SRR514250_again"


mkdir -p $output_folder

 NXF_VER=21.10.6 NXF_VER=21.10.6 NXF_OPTS="-Dleveldb.mmap=false -Xmx500g"  nextflow  -log "${output_folder}/nextflow.log" nextflow run metacompass.nf \
      --forward "$forward_read" \
      --reverse "$reverse_read" \
      --output "$output_folder"
      --trace_file_name "$output_folder/trace.txt"
      -with-dag "$output_folder/${read}_dag.png"
      -with-timeline "$output_folder/timeline.html"

