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
unmapped_reads=$(grep "aligned 0 times" $bowtie2reads |awk '{print $1}')
mapped_reads=$(($total_reads-$unmapped_reads))
bowtie2=$(grep overall $bowtie2reads |cut -f1 -d '%')
pilon_paired=$(grep 'were paired' $bowtie2contigs|awk '{print $1}')
pilon_unpaired=$(grep 'were unpaired' $bowtie2contigs|awk '{print $1}')
pilon1=$(grep overall $bowtie2contigs |cut -f1 -d '%'|head -n1)
pilon2=$(grep overall $bowtie2contigs |cut -f1 -d '%'|tail -n1)
if [ -z $pilon_paired ];then
	pilon1=0
fi
if [ -z $pilon_unpaired ];then
	pilon2=0
fi

printf "sample_name\tnum_selected_refs\tnum_assembled_refs\tnum_reads\tnum_mapped_reads\tmapped_reads%%\tmapped_contigs_paired%%\tmapped_contigs_unpaired%%\n">$output1

printf "contigs.pilon.fasta\t">>$output1
printf "%d\t" ${selected_refs} ${assembled_refs} ${total_reads} ${mapped_reads}>>$output1
printf "%s\t" ${bowtie2} ${pilon1} ${pilon2}>>$output1
printf "\n">>$output1


assembled_reads=$(cat $bowtie2sam |cut -f3 |sort|uniq -c| awk '{print $2"\t"$1}')
printf "ref_ID\tnum_mapped_reads\tmapped_reads%%\tassembled_reads%%\n">$output2


IFS=$'\n' 
for i in $assembled_reads;do
	id=$(echo $i|cut -f1)
	reads=$(echo $i|cut -f2)
	echo $id $reads ${total_reads} ${mapped_reads} | awk '{printf $1"\t"$2"\t";if($3>0)printf $2/$3"\t";else printf "0\t";if($4>0)printf $2/$4"\n";else printf "0\n";}'>>$output2
done
