#!/bin/bash

path=$1 #metacompass path
out=$2 #metacompass out directory
assembly=$3 #metacompass assembly directory
contigs=$4 #$dir/intermediate_files/assembly_output/contigs.fasta
minlen=$5 
#echo "parameters:"
#echo "path $1 $path"
#echo "out $2 $out" 
#echo "assembly $3 $assembly"
#echo "contigs: $4 $contigs"
#echo "minlen $5 $minlen"

samtools faidx $contigs

if [ ! -d $out/$assembly/contigs ]
then
	mkdir -p $out/$assembly/contigs
	$path/bin/extractContigs $contigs $out/$assembly/contigs >$out/$assembly/extractContigs.log

fi
declare -i count=1
for contig in $(ls $out/$assembly/contigs|sed 's/.fasta//g');do
if [ $count -eq 1 ]
then
	#echo "python $path/bin/assembly_stats.py $out/$assembly/contigs/${contig}.fasta $minlen" #> $out/metacompass_assembly_pergenome_stats.tsv	
	python $path/bin/assembly_stats.py $out/$assembly/contigs/${contig}.fasta $minlen #> $out/metacompass_assembly_pergenome_stats.tsv
else
	#echo "python $path/bin/assembly_stats.py $out/$assembly/contigs/${contig}.fasta $minlen |grep -v File"
	python $path/bin/assembly_stats.py $out/$assembly/contigs/${contig}.fasta $minlen |grep -v sample_name #>> $out/metacompass_assembly_pergenome_stats.tsv
fi
count=$((count +1))
done

