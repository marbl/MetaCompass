#!/bin/bash
#SBATCH -J build_refDB # Job name
#SBATCH -o build_refDB.o%j # Name of output file
#SBATCH -e build_refDB.e%j # Name of error file
#SBATCH --mail-user=ujjwal.097@gmail.com # Email for job info
#SBATCH --mail-type=BEGIN,END,FAIL # Get email for begin, end, and fail
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --account=cbcb
#SBATCH --mem=36gb

# must activate metacompass conda environment before running
# database_dir=/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2
database_dir=${1}
accession_file=${2}
accession_path=${database_dir}/RefSeq_V2_db/data/${accession_file}

## build_ref.sh -------------------------------------------------------------------------------------

version_tag=${3} # optional input argument
# version_tag='' # optional input argument

cd ${database_dir}

export PATH=$PATH:${database_dir}/fetchMG

# logfile=get_sequences_log_${version_tag}.txt
# yes | python ${database_dir}/get_sequences_for_refpkg_updated.py  ${database_dir} ${accession_path} ${version_tag} > ${logfile}
yes | python ${database_dir}/get_sequences_for_refpkg_updated.py  ${database_dir} ${accession_path} ${version_tag}
