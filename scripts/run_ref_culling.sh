#!/bin/bash
#SBATCH -J nextflow_2
#SBATCH -o nextflow_2.o
#SBATCH -e nextflow_2.e
#SBATCH --mail-user=tluan@terpmail.umd.edu
#SBATCH --mail-type=all
#SBATCH --time=18:00:00
#SBATCH --qos=highmem
#SBATCH --partition=cbcb
#SBATCH --mem=64gb
#SBATCH --account=cbcb
#SBATCH --ntasks=8

#rm -rf /fs/nexus-scratch/tluan/skani/SRR628270_new
#mkdir /fs/nexus-scratch/tluan/skani/SRR628270_new
#source /fs/cbcb-scratch/tluan/conda_2/bin/activate /fs/cbcb-scratch/tluan/conda_2
#conda activate metacompass

ref_db="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/marker_gene_ref_pkg/ref_pkg_generation_V3"
forward="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/kmc/SRR628270/forward.fastq"
reverse="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/kmc/SRR628270/reverse.fastq"
ref_cans="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/MetaCompass_paper_data/kmc/SRR628270/reference_candidates.txt"
output_dir="/fs/cbcb-scratch/ayyangar/metacompass/metacompass_new/results/output"
culling_out_dir="${output_dir}/culling_output"
mkdir -p ${culling_out_dir}
align_reads_out_dir="${output_dir}/align_reads"
mkdir -p ${align_reads_out_dir}

work_dir="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/marker_gene_ref_pkg/ref_pkg_generation_V2/ref_download/metacompass-develop-nextflow_Jun_26/metacompass-develop-nextflow"

nextflow run ref_culling.nf --reference_db ${ref_db} --forward ${forward} --reverse ${reverse}  --ref_candidates ${ref_cans} --ms 16 --num_matches 2  --threads 8 --output ${output_dir} --workdir ${work_dir} --cull_stop 0.01