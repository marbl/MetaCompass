######################################################################
#  PUBLIC DOMAIN NOTICE
#
#  This software is "United States Government Work" under the terms of the United
#  States Copyright Act. It was written as part of the authors' official duties
#  for the United States Government and thus cannot be copyrighted. This software
#  is freely available to the public for use without a copyright
#  notice. Restrictions cannot be placed on its present or future use.
#
#  Although all reasonable efforts have been taken to ensure the accuracy and
#  reliability of the software and associated data, the National Human Genome
#  Research Institute (NHGRI), National Institutes of Health (NIH) and the
#  U.S. Government do not and cannot warrant the performance or results that may
#  be obtained by using this software or data. NHGRI, NIH and the U.S. Government
#  disclaim all warranties as to performance, merchantability or fitness for any
#  particular purpose.
#
#  Please cite the authors in any work or product based on this material.
######################################################################
#
# This script can be run as sh extractSRARuns.sh <jobid> where jobid is an index in a list of folders containing SRA runs for a sample. 
# All SRA files in a folder are merged and filtered by DACC protocol
# http://hmpdacc.org/doc/ReadProcessing_SOP.pdf
#
# SRA tools must be in your path

jobid=$SGE_TASK_ID
if [ x$jobid = x -o x$jobid = xundefined -o x$jobid = x0 ]; then
jobid=$1
fi

if test x$jobid = x; then
  echo Error: I need SGE_TASK_ID set, or a job index on the command line
  exit 1
fi

# now figure out which contig we are
NUM_JOBS=`ls -lhad SRS*  | wc -l |awk '{print $1}'`
if [ $jobid -le 0 ]; then
   echo "Invalid job id, must be 1 or greater"
   exit
fi

if [ $jobid -gt $NUM_JOBS ]; then
   echo "Invalid job id, max is $NUM_JOBS"
   exit
fi

sample=`ls -lhad SRS* |awk '{print $NF}' |sort |head -n $jobid |tail -n 1`

echo "Processing $sample"
   if [ ! -e $sample/$sample.trimmed.1.fastq ]; then
      cd $sample
      for file in `ls *.sra`; do
         run=`echo $file |sed s/.sra//g`
         fastq-dump --split-3 -E -A $run -D $run/ --defline-seq '@$sn/$ri' --defline-qual '+$sn/$ri' -O . $file
      done
      cat SRR*_1.fastq > ${sample}_1.fastq
      cat SRR*_2.fastq > ${sample}_2.fastq
      rm -f SRR*.fastq
#      rm -f SRR*.sra

      # now convert to bam
      java -Xmx256g -jar ./picard.jar FastqToSam F1=${sample}_1.fastq F2=${sample}_2.fastq O=$sample.duplicates.bam SM=$sample SORT_ORDER=coordinate MAX_RECORDS_IN_RAM=150000000
      # remove duplicates 
      java -Xmx256g -jar ./picard.jar EstimateLibraryComplexity I=$sample.duplicates.bam O=$sample.summary.txt MAX_RECORDS_IN_RAM=150000000
      # trim by qv
      perl ./trimBWAstyle.usingBam.pl -f $sample.bam -q2
      rm -f SRS*.bam
      cd -
   else
      echo "Already done"
   fi
