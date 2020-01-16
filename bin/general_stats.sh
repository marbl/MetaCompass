#!/bin/bash

#   Calculation of ...

#   This code creates an output.txt file with all of the statistics

#   Usage: python general_stats.py <> <> <> <>

all_refs_ids=$1
assembled_refs_ids=$2
bowtie2reads=$3 #'%s/%s.%s.bowtie2map.log
bowtie2contigs=$4 #'%s/%s.%s.pilon.map.log
bowtie2sam=$5 #'%s/%s.%s.pilon.map.log only reads used by builcontig
output1=$6
output2=$7
        
selected_refs=$(wc -l $all_refs_ids|cut -f1 -d ' ')
assembled_refs=$(wc -l $assembled_refs_ids|cut -f1 -d ' ')
total_reads=$(head -n1 $bowtie2reads |cut -f1 -d ' ')
mapped_reads=$(grep "aligned 0 times" $bowtie2reads |tr -d ' '|cut -f1 -d '(')
bowtie2=$(grep overall $bowtie2reads |cut -f1 -d ' ')
pilon=$(grep overall $bowtie2contigs |cut -f1 -d ' ')
#printf "Selected References\tAssembled References\tMapped reads %\tMapped contigs %paired-end\tMapped contigs %single-end\n"
printf "# Selected References\t# Assembled References\t# reads\t# Mapped reads\tMapped reads%%\tMapped contigs to paired-end%%\tMapped contigs to single-end%%\n">$output1
printf "%d\t" ${selected_refs} ${assembled_refs} ${total_reads} ${mapped_reads}>>$output1
printf "%s\t" ${bowtie2} ${pilon}>>$output1
printf "\n">>$output1


assembled_reads=$(cat $bowtie2sam |cut -f3 |sort|uniq -c| awk '{print $2"\t"$1}')
printf "Ref ID\t#Mapped reads\tMapped reads%%\tAssembled reads%%\n">$output2

IFS=$'\n' 
for i in $assembled_reads;do
	id=$(echo $i|cut -f1)
	reads=$(echo $i|cut -f2)
	echo $id $reads ${total_reads} ${mapped_reads} | awk '{print $1"\t"$2"\t"$2/$3"\t"$2/$4}'>>$output2
done