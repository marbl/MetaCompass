This gets the data needed to make the FRC curves. You can see an example of how to run it in run_FRC_All.sh

Output is a TSV file
Each row is one contig. Columns are:

Contig_Name - The name of this contig
Contig_Len - The length of this contig in BP
Lo_DOC - Number of low read coverage features found on this contig (LOW_COV_PE)
Hi_DOC - Number of many high read coverage features found on this contig (HIGH_COV_PE)
Lo_MP_DOC - Number of low fragment coverage features found on this contig (LOW_NORM_COV_PE)
Hi_MP_DOC - Nomber of high fragment coverage features found on this contig (HIGH_NORM_COV_PE)
Short_MP_DOC - Number of shortened mate-pair features found on this contig (COMPR_PE)
Long_MP_DOC - Number of long mate-pair features found on this contig (STRETCH_PE)
Singleton_DOC -  Number of singleton (reads with no mate) features found on this contig (HIGH_SINGLE_PE)
Misoriented_DOC - Number of misoriented mate-pair features found on this contig (HIGH_OUTIE_PE)

I listed the features from GAGE that these are roughly equivalent to, though I did NOT use the same methods to decide which features are significant, so they are not exactly the same. Instead, I used the methods that VALET was already using. So the resulting curves may look different. Also, HIGH_SPAN_PE is not included yet.

Here's an example of how to run it, including the preproscessing steps of creating and properly sorting the BED file from the SAM file. If you already have everything you need, you can just run the command at the bottom.

#files
sam="/fs/cbcb-scratch/jhauzel/metacompass/final_tests/metacompass_out_sha/frc/all.sam"

bam="/fs/cbcb-scratch/jhauzel/metacompass/final_tests/metacompass_out_sha/frc/all.bam"

bed="/fs/cbcb-scratch/jhauzel/metacompass/final_tests/metacompass_out_sha/frc/all.bed"

fasta="/fs/cbcb-scratch/jhauzel/metacompass/final_tests/metacompass_out_sha/assembly_merge/contigs.merge.fasta"

samtools view -Sbh $sam > $bam;
bedtools bamtobed -i $bam > $bed;
samtools faidx $fasta
env LC_COLLATE=C sort -k 1,1 -k 4,4 -o $bed $bed


fai="${fasta}.fai"
out="FRC_Data_All.tsv"

#data properties
libsize_avg="165"
libsize_sdv="92.8"
#read_len="101"
orientation_a="+"
orientation_b="-"

#statistical parameters
mp_size_z_cutoff="2.58"
mp_edge_z_cutoff="2.58"
doc_logp_cutoff="7"


./get_FRC_data $bed $fai $out $libsize_avg $libsize_sdv $orientation_a $orientation_b $mp_size_z_cutoff $mp_edge_z_cutoff $doc_logp_cutoff


