'''

   Calculation of [Total # of Contigs], [Total Length], [Total # of trimmed Contigs], [Trimmed Length], [GC content],
   [Min Contig Size [bp]], [Median Contig Size [bp]], [Mean Contig Size [bp]], [Max Contig Size [bp]],
   [N50[bp] [# of Contigs],
   [Total # of Contigs]

   This code creates an output.txt file with all of the statistics

   Usage: python assembly_stats.py <assembly.fasta> <minimum contig size>

   Author: Nicolas Schmelling

'''

from __future__ import division
from Bio import SeqIO
import sys

def assembly_stats(contigsMultifasta, mini):
    
    contigsLength = []
    trimmedLength = []
    total = 0
    sum = 0
    thres = mini - 1
    GC_count = 0
    GC_cont  = 0           
    
    seqs = open(contigsMultifasta, 'r')

    # Create lists for Total Length and Trimmed Length
    for seq_record in SeqIO.parse(open(contigsMultifasta), 'fasta'):
        contigsLength.append(len(seq_record.seq))
        total += len(seq_record.seq)
        # Min Contig Length Threshold
        if len(seq_record.seq) > thres: 
            sum += len(seq_record.seq)
            trimmedLength.append(len(seq_record.seq))
        
    # Calculating GC content
    for seq_record in seqs.read():
        if seq_record.startswith('>'):
            continue
        else:
            if 'G' in seq_record:
                GC_count += 1
            if 'C' in seq_record:
                GC_count += 1
    GC_cont = float((GC_count/total)*100)

    # Sorting the Trimmed Contigs from Large to Small
    trimmedLength.sort()
    trimmedLength.reverse()
    

    # Calculating Mean Contig Size
    meancontig = int(sum/len(trimmedLength))

    # Calculating Median Contig Size
    median = []
    medcon = []

    for con in trimmedLength:
        medcon.append(con)
        if len(medcon) > len(trimmedLength)/2:
            median.append(con)
            break

    # Theoretic NXX and NGXX Sizes
    teoN50 = sum / 2.0 
    # Checking N50 [bp] [# of Contigs]  
    totalSum = 0
    contigcount=0
    N50 = 0
    N50con = 0
    n1m = 0
    s1m = 0
    n2m = 0
    s2m = 0
    n4m = 0
    s4m = 0
    n10m = 0
    s10m = 0
    for con in trimmedLength:
        totalSum += con
        contigcount += 1
        if totalSum >= teoN50 and N50 == 0:
            N50 = con
            N50con = contigcount
        if totalSum >= 1000000 and n1m == 0:
    	    n1m = contigcount
    	    s1m = con
        if totalSum >= 2000000 and n2m == 0:
    	    n2m = contigcount
    	    s2m = con
        if totalSum >= 4000000 and n4m == 0:
    	    n4m = contigcount
    	    s4m = con
        if totalSum >= 10000000 and n10m == 0:
    	    n10m = contigcount
    	    s10m = con        
    #print '# of trimmed contigs: ' + str(len(trimmedLength))
    #print 'trimmed size [bp]: ' + str(sum)
    print ("File\t# Contigs\tTotal Size(Kbp)\tMin Size\tMax Size(Kbp)\tAverage Size\tMedian Size\tN50\t# N50 contigs\tSize at 1Mbp (Kbp)\tNumber @ 1Mbp\tSize at 2Mbp (Kbp)\tNumber @ 2Mbp\tSize at 4Mbp (Kbp)\tNumber @ 4Mbp\tSize at 10Mbp (Kbp)\tNumber @ 10Mbp\tGC content [%]")
    print ("%s\t%d" % (contigsMultifasta,len(contigsLength)) + '\t' \
    + "%d" % (total) + '\t' \
    + "%d\t%d\t%.2f\t%.2f" % (min(trimmedLength),max(trimmedLength),meancontig,median[0]) + '\t' \
    + "%d\t%d\t%.2f\t%d\t%.2f\t%d" % (N50,N50con,s1m/1000,n1m,s2m/1000,n2m) + '\t' \
    + "%.2f\t%d%.2f\t%d" % (s4m/1000,n4m,s10m/1000,n10m) + '\t' \
    + "%.2f" % (GC_cont) ) 
if __name__ == "__main__":
    contigsMultifasta = sys.argv[1]
    mini = int(sys.argv[2])
    
    assembly_stats(contigsMultifasta, mini)
