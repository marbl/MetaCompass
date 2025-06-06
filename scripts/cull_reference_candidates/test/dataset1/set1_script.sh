#!/bin/bash

# `-m 2`: only choose a reference candidate if it shares at least 2 kmers with the reads
# `-c set1_ref_kmer_folders/set1_ref_folder1`: folder of reference kmer files.
# Can input more reference kmer folders using `-c path_to_ref_kmer_folder`
# `-o .` : output folder. Result will be in output_folder/min_reference_candidates.txt
# ../../refcull -m 1 \
#     -r set1_read_kmer.txt \
#     -c set1_ref_kmer_folders/set1_ref_folder1 \
#     -c set1_ref_kmer_folders/set1_ref_folder2 \
#     -o .

# for cpp implementation
../../refcull -m 2 \
    -r set1_read_kmer.fasta \
    -c set1_ref_kmer_folders/set1_ref_folder_combined \
    -o .

# for python implementation
# -r arg need a directory containing the read_kmer files
# read_kmer files can be split or kept in one file
python ../../../ref_culling.py -m 2 \
    -r read_kmers \
    -c set1_ref_kmer_folders/set1_ref_folder_combined \
    -o .