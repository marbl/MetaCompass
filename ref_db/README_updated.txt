Building Marker Gene Reference Package - Updated 2023

# First we have to download bacterial and archeal reference genomes and extract and format marker gene data 
#     Based off of https://github.com/shahnidhi/TIPP_reference_package/blob/master/src/get_sequences.py
#     Relies on FetchMG_COGid2gene_name.tsv from TIPP_reference_package
#     Relies on list of prokaryotes from https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt
# The list of prokaryotes (prokaryotes.txt) can be replaced with a custom accession list to make a custom database if desired.
#     If an alternate list is used with a different file name, change usage of "prokaryotes" below and in "submit_build_ref.sh" 
# 
# This first step takes a long time, and as a result should be split up among multiple jobs.
# To do this, we split prokaryotes.txt into several smaller files, run the build_ref scripts, and then combine the outputs.

# This process relies on some packages in the MetaCompass conda environment
# ** Make sure MetaCompass Conda environment is activated for these steps!!! **

# Set repository path (example parent path shown)

repository_path=/fs/cbcb-scratch/zbowen/metacompass/ref_db

# Download the accession text file (prokayotes.txt) and place in /RefSeq_V2_db/data/
cd ${repository_path}/RefSeq_V2_db/data
wget https://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt

# Split the accession text file:

split -l 5000 ${repository_path}/RefSeq_V2_db/data/prokaryotes.txt ${repository_path}/RefSeq_V2_db/data/prokaryotes_

# Submit several jobs that will take these input files as inputs

cd $repository_path
./submit_build_ref.sh ${repository_path}

# --------------------------------------------------------------------------------------------------------------------------------
# Once ALL jobs complete successfully, combine the outputs:

cd $repository_path
./combineOutputs.sh ${repository_path}

# (Optional) Remove split directories now that everything is combined:

rm -r marker_index_[a-z][a-z]

# Run the rest of the reference database setup:

cd $repository_path
./setup_ref_db.sh ${repository_path} prokaryotes.txt

