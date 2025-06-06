#!/bin/bash

# This is assuming you have already done a step to split the input file, for example:
# split -l 5000 /fs/cbcb-lab/mpop/Fraunhofer_metagenomics/marker_gene_ref_pkg/ref_pkg_generation_decontaminated_full/RefSeq_V2_db/data/prokaryotes.txt /fs/cbcb-lab/mpop/Fraunhofer_metagenomics/marker_gene_ref_pkg/ref_pkg_generation_decontaminated_full/RefSeq_V2_db/data/prokaryotes_

database_dir=${1}

for file in ${database_dir}/RefSeq_V2_db/data/prokaryotes_*; do 
	echo $file;
	basefile=$(basename $file)
	echo $basefile
	TMP=`echo $file | sed 's/^.*prokaryotes_//g'`; 
	echo $TMP

	# Can we replace sbatch here for generality purposes?
	sbatch build_ref.sh ${database_dir} $basefile $TMP 

	# #re-run failed jobs (if needed)
	# # EXAMPLE: Check if TMP suffix is within 'br' through 'ck'
    # if [[ "$TMP" =~ ^(br|bs|bt|bu|bv|bw|bx|by|bz|ca|cb|cc|cd|ce|cf|cg|ch|ci|cj|ck)$ ]]; then
	# 	echo $TMP
	# 	sbatch build_ref.sh ${database_dir} $basefile $TMP 
    # fi

done