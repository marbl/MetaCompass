#!/bin/bash
#SBATCH -J reference_selection # Job name
#SBATCH -o reference_selection.o%j # Name of output file
#SBATCH -e reference_selection.e%j # Name of error file
#SBATCH --mail-user=meiselj@umiacs.umd.edu # Email for job info
#SBATCH --mail-type=all # Get email for begin, end, and fail
#SBATCH --time=100:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --qos=large
#SBATCH --mem=128gb

#----------------------------------------#
# define variables and files
#----------------------------------------#

refDir="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/RefSeq_V1_db/"
outdir="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/ref_selection_test/"
file1="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/hmp_stool/SRS011302.denovo_duplicates_marked.trimmed.1.fastq"
file2="/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/hmp_stool/SRS011302.denovo_duplicates_marked.trimmed.2.fastq"
declare -a files=($file1 $file2)
paired="true"
threads=12
readLength=105
out="SRS011302"
depthOfCoverage=5
breadthOfCoverage=0.9
percent_markers_covered=75

# For each marker gene in the reference database
 for gene in $(ls ${refDir}marker_index); do
 	
 	echo $gene
	
	#----------------------------------------#
	# run kmer-mask to filter references
	#----------------------------------------#

	## IF PAIRED END READS, NEED TO PROCESS THROUGH BOTH FORWARD AND \
	## REVERSE READ AND CONCATENATE OUTPUT
	## SYNTAX PROBABLY NEEDS TO BE FIXED
	for file in ${files[@]}
    do
		/cbcb/sw/RedHat-7-x86_64/users/meiselj/src/meryl-r2013/meryl/kmer-mask \
			-mdb ${refDir}marker_index/${gene}/${gene}_clustered.ms28 \
 			-1 $file \
 			-ms 28 \
 			-clean 0.0 \
 			-match 0.01 \
 			-nomasking \
 			-t $threads \
  			-l $readLength \
  			-o ${outdir}${out}_${gene}$RANDOM
	done
	
	# Keep only files containing masked sequences (found in the reference)
  	cat ${outdir}${out}_${gene}*.match.fastq > ${outdir}${out}_${gene}.keep.fastq
  	
	# Remove files we don't need
	rm *clean*fastq
	rm *murky*fastq
	rm *mixed*fastq
	rm *match*fastq
	
	#----------------------------------------#
	# run bowtie2 to map to references
	#----------------------------------------#
	module load bowtie2/2.3.0
	module load samtools/1.5
	
	# Align reads to marker gene reference
	bowtie2 --threads $threads \
	    -x ${refDir}marker_index/${gene}/${gene}_clustered \
	    -U ${outdir}${out}_${gene}.keep.fastq | \
		samtools view -bS > ${out}_${gene}.match.bam

	# Sort alignment file
	samtools sort ${out}_${gene}.match.bam -o ${out}_${gene}.match.sorted.bam
	rm ${out}_${gene}.match.bam
	
	#----------------------------------------#
	# use bedtools to compute coverage
	#----------------------------------------#
	
	module load bedtools/2.26.0
	
	# Make a temp file with the marker gene lengths
	tail -n+2 ${refDir}marker_index/${gene}/${gene}_genome.tsv | cut -f 1,4 > ${outdir}${gene}_genome.txt
	
	# Caclulate coverage
	# -max flag says that anything equal to or above the specified coverage is output into a single line
	bedtools genomecov -ibam ${outdir}${out}_${gene}.match.sorted.bam \
	  -g ${outdir}${gene}_genome.txt \
	  -max $depthOfCoverage > ${outdir}${out}_${gene}_genomeCov.txt
	
	# Remove temp file with the marker gene lengths
	rm ${outdir}${gene}_genome.txt
	
	# Extract markers with good breadth of coverage at specified depth
	# First column is the marker gene name, 
	# Second column is depth of coverage and 
	# fifth column is the percentage of the gene covered at the specified depth (breadth)
	echo -e "seq\tcoverage" > ${outdir}${out}_${gene}_marker_cov.txt
	awk -v doc=$depthOfCoverage \
		-v boc=$breadthOfCoverage \
		'{if($2 == doc && $5 >= boc){print $1"\t"$5}}' ${out}_${gene}_genomeCov.txt >> \
		${outdir}${out}_${gene}_marker_cov.txt

	# Markers are representative of a cluster of sequences
	# For each marker that is covered, pull out the corresponding sequences in the cluster
	python marker_coverage_per_candidate_genome.py $gene \
		${outdir}${out}_${gene}_marker_cov.txt \
		${outdir}${out}_${gene}_marker_cov_per_genome.txt
	
	# Remove unnecessary files
	rm ${outdir}${out}_${gene}.keep.*.fastq
	rm ${out}_${gene}.match.sorted.bam
	rm ${out}_${gene}_genomeCov.txt
	rm ${outdir}${out}_${gene}_marker_cov.txt
	
done

# Concatenate coverage for all genes
cat *_marker_cov_per_genome.txt >> marker_cov_per_genome.txt
rm *_marker_cov_per_genome.txt 

# Output a table with the reference genomes and how many marker genes are covered
# Output a list of candidate reference genomes
python identify_candidate_genomes.py $percent_markers_covered
rm marker_cov_per_genome.txt

