#!/bin/bash

# This is the all-in-one setup script for running a relatively SMALL input file
# Must activate metacompass conda environment before running
# Example usage:
# ./setup_ref_db.sh /fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2 prokaryotes.txt

# database_dir=/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2
database_dir=${1}
accession_file=${2}
accession_path=${database_dir}/RefSeq_V2_db/data/${accession_file}

## build_ref.sh -------------------------------------------------------------------------------------

cd ${database_dir}

export PATH=$PATH:${database_dir}/fetchMG

yes | python ${database_dir}/get_sequences_for_refpkg_updated.py  ${database_dir} ${accession_path}

# Remove uncleaned ncbi-datasets files
if [ -d "${database_dir}/RefSeq_V2_db/marker_index/ncbi_dataset" ]; then
    rm -r "${database_dir}/RefSeq_V2_db/marker_index/ncbi_dataset"
    rm "${database_dir}/RefSeq_V2_db/marker_index/ncbi_dataset.zip"
    rm "${database_dir}/RefSeq_V2_db/marker_index/README.md"
    echo "Extraneous ncbi datasets files removed."
fi

## generate_indices.sh -------------------------------------------------------------------------------------


# Make these available in conda environment
# module load cdhit
# PATH=$PATH:/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/env/bin/meryl-r2013/meryl
# module load bowtie2/2.3.4

WRKPTH="${database_dir}/RefSeq_V2_db/marker_index"

for file in $( ls $WRKPTH ); do  
#for file in $( ls $WRKPTH | head -n 1 ); do	
	echo $file; 

	echo "Clustering"

	cd-hit-est -i ${WRKPTH}/${file}/${file}.fna \
		-c 0.99 -n 10 -d 0 -M 128000 -T 12 \
		-o ${WRKPTH}/${file}/${file}_clustered

	mv ${WRKPTH}/${file}/${file}_clustered \
		${WRKPTH}/${file}/${file}_clustered.fna

	echo "Building bowtie2 index"

	bowtie2-build ${WRKPTH}/${file}/${file}_clustered.fna \
  		${WRKPTH}/${file}/${file}_clustered

  	echo "Building kmer-mask database"

  	# Build database
  	meryl -B -C -m 28 -s ${WRKPTH}/${file}/${file}_clustered.fna \
  		-o ${WRKPTH}/${file}/${file}_clustered.ms28

  	rm ${WRKPTH}/${file}/${file}_clustered*.fastaidx
  	
  	#Get stats
  	meryl -Dc -s ${WRKPTH}/${file}/${file}_clustered.ms28
done




## org_clusters.sh -------------------------------------------------------------------------------------

WRKPTH="${database_dir}/RefSeq_V2_db/marker_index/"

for file in $( ls $WRKPTH ); do  
	
	echo $file; 

	# python /fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2/org_clusters.py $file
	python ${database_dir}/org_clusters.py $file $WRKPTH
done

## get_taxinfo.sh -------------------------------------------------------------------------------------

cd ${database_dir}

#echo "Getting accession info"
tail -n+2 RefSeq_V2_db/marker_index/*/*_genome.tsv | cut -f 3 | sort -u >> ./RefSeq_V2_db/data/used_accessions.txt

#echo "Getting taxonomy info"
python assign_taxa_label_mc.py ./RefSeq_V2_db/data/used_accessions.txt ./RefSeq_V2_db/data/${accession_file} > ./RefSeq_V2_db/data/taxonomy.info

#python code was written with v2-- so when I used python3, I got error AttributeError: 'itertools._grouper' object has no attribute 'next'
echo "Getting genes per acc data"

# python get_genes_per_acc_mc.py ${database_dir}/RefSeq_V2_db/marker_index ${database_dir}/RefSeq_V2_db/taxinfo.tsv
python get_genes_per_acc_mc.py ${database_dir} ${accession_path}
# python get_genes_per_acc_mc.py /fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2

##  -------------------------------------------------------------------------------------

## remove extraneous marker gene
if [ -d "${database_dir}/RefSeq_V2_db/marker_index/RpoC_COG0086" ]; then
    rm -r "${database_dir}/RefSeq_V2_db/marker_index/RpoC_COG0086"
    echo "Extraneous marker gene (RpoC_COG0086) removed."
else
    echo "Extraneous marker gene did not exist."
fi

## do this right after the first step -> NEED TO ALSO CHECK FOR AND REMOVE LEFTOVER NCBI DATASET FILES (README, ncbi_dataset directory, and zip folder)

echo "Finished reference database setup!"