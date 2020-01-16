'''

   Calculation of [Total # of Contigs], [Total Length], [Total # of trimmed Contigs], [Trimmed Length], [GC content],
   [Min Contig Size [bp]], [Median Contig Size [bp]], [Mean Contig Size [bp]], [Max Contig Size [bp]],
   [N50[bp] [# of Contigs],
   [Total # of Contigs]

   This code creates an output.txt file with all of the statistics

   Usage: python assembly_stats.py <assembly.fasta> <minimum contig size>

  modified code from Nicolas Schmelling

'''

from __future__ import division
from Bio import SeqIO
from statistics import median
import sys
import os.path 

def median(trimmedLength):
    n = len(trimmedLength)
    if n < 1:
        return None
    if n % 2 == 1:
        return trimmedLength[n//2]
    else:
        return sum(trimmedLength[n//2-1:n//2+1])/2.0    
        
def assembly_stats(contigsMultifasta, minContigSize):
    
    #contigsLength = []
    trimmedLength = []
    sum_trimmed = 0
    thres = minContigSize - 1
    GC_count = 0
    GC_cont  = 0           
    
    seqs = open(contigsMultifasta, 'r')

    # Create lists for Total Length and Trimmed Length and Calculates GC content
    for seq_record in SeqIO.parse(open(contigsMultifasta), 'fasta'):                  
        #contigsLength.append(len(seq_record.seq))
        # Min Contig Length Threshold
        if len(seq_record.seq) > thres:
            sum_trimmed += len(seq_record.seq)
            trimmedLength.append(len(seq_record.seq))
            GC_count+=seq_record.seq.count("C")
            GC_count+=seq_record.seq.count("G")
            
    GC_cont = float((GC_count/sum_trimmed)*100)
                    
    # Sorting the Trimmed Contigs from Large to Small
    trimmedLength.sort()
    trimmedLength.reverse()
    

    # Calculating Mean Contig Size
    meancontig = int(sum_trimmed/len(trimmedLength))
    # Calculating Median Contig Size
    mediancontig=median(trimmedLength)
    # Checking N50 [bp] [# of Contigs]  
    totalSum = 0
    contigcount=0
    N50 = 0
    N50con = 0
    n1m = 0
    s1m = 0.0
    n2m = 0
    s2m = 0.0
    n4m = 0
    s4m = 0.0
    n10m = 0
    s10m = 0.0
    #contig len and contigtrimmedlen if cutoff >1
    for contiglen in trimmedLength:
        #print (contig)
        totalSum += contiglen
        contigcount += 1
        if totalSum >= sum_trimmed/2.0 and N50 == 0:
            N50con = contigcount
            N50 = contiglen
        if totalSum >= 1000000 and n1m == 0:
    	    n1m = contigcount
    	    s1m = contiglen
        if totalSum >= 2000000 and n2m == 0:
    	    n2m = contigcount
    	    s2m = contiglen
        if totalSum >= 4000000 and n4m == 0:
    	    n4m = contigcount
    	    s4m = contiglen
        if totalSum >= 10000000 and n10m == 0:
    	    n10m = contigcount
    	    s10m = contiglen      
 #   print ('N50 Length [bp]: ' + str(N50))
 #   print ('# N50 Contig: ' + str(N50con))
 #   print ('Size   @1Mbp [bp]: ' + str(s1m))
 #   print ('Number @1Mbp: ' + str(n1m))  
 #   print ('Size   @2Mbp [bp]: ' + str(s2m))
 #   print ('Number @2Mbp: ' + str(n2m))
 #   print ('Size   @4Mbp [bp]: ' + str(s4m))
 #   print ('Number @4Mbp: ' + str(n4m))
 #   print ('Size   @10Mbp [bp]: ' + str(s10m))
 #   print ('Number @10Mbp: ' + str(n10m)) 
 #   print '# of trimmed contigs: ' + str(len(trimmedLength))
 #   print 'trimmed size [bp]: ' + str( sum_trimmed)
 #   print ('trimmed total size [bp]: ' + str(sum_trimmed))
 #   print ('trimmed total size [bp] divided by 2: ' + str(sum_trimmed/2.0))
    
    print ("File\t# Contigs\tTotal Size(Kbp)\tMin Size\tMax Size(Kbp)\tAverage Size\tMedian Size\tN50\t# N50 contigs\tSize at 1Mbp (Kbp)\tNumber @ 1Mbp\tSize at 2Mbp (Kbp)\tNumber @ 2Mbp\tSize at 4Mbp (Kbp)\tNumber @ 4Mbp\tSize at 10Mbp (Kbp)\tNumber @ 10Mbp\tGC content [%]")
    print ("%s\t%d" % (os.path.basename(contigsMultifasta.replace(".fasta", "")),len(trimmedLength)) + '\t' \
    + "%d" % (sum_trimmed) + '\t' \
    + "%d\t%d\t%.2f\t%.2f" % (min(trimmedLength),max(trimmedLength),meancontig,mediancontig) + '\t' \
    + "%d\t%d\t%.2f\t%d\t%.2f\t%d" % (N50,N50con,s1m/1000,n1m,s2m/1000,n2m) + '\t' \
    + "%.2f\t%d\t%.2f\t%d" % (s4m/1000,n4m,s10m/1000,n10m) + '\t' \
    + "%.2f" % (GC_cont) ) 

if __name__ == "__main__":
    contigsMultifasta = sys.argv[1]
    minContigSize = int(sys.argv[2])
    assembly_stats(contigsMultifasta, minContigSize)
