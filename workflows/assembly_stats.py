'''

   Calculation of [Total # of Contigs], [Total Length], [Total # of trimmed Contigs], [Trimmed Length], [GC content],
   [Min Contig Size [bp]], [Median Contig Size [bp]], [Mean Contig Size [bp]], [Max Contig Size [bp]],
   [N50[bp] [# of Contigs]], [NG50[bp] [# of Contigs]], [N90 [bp] [# of Contigs]], [NG90 [bp] [# of Contigs]],
   [Total # of Contigs > Average Gene Size]

   This code creates an output.txt file with all of the statistics

   Usage: python assembly_stats.py <assembly.fasta> <minimum contig size> <estimated genome size> <average gene size>

   Author: Nicolas Schmelling

'''

from __future__ import division
from Bio import SeqIO
import sys

def assembly_stats(contigsMultifasta, mini, est_genome_size, average_gene_size):
    
    contigsLength = []
    trimmedLength = []
    total = 0
    sum = 0
    thres = mini - 1
    GC_count = 0
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
    
    # Theoretic NXX and NGXX Sizes
    teoN50 = sum / 2.0 
    teoNG50 = est_genome_size / 2.0
    teoN90 = sum * 0.9
    teoNG90 = est_genome_size * 0.9

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

    # Checking N50 [bp] [# of Contigs]
    testSum = 0
    N50 = 0
    N50con = 0
    for con in trimmedLength:
        testSum += con
        N50con += 1
        if teoN50 < testSum:
            N50 = con
            break

    # Checking NG50 [bp] [# of Contigs]
    testSum = 0
    NG50 = 0
    NG50con = 0
    for con in trimmedLength:
        if sum < (est_genome_size/2):
            break
        testSum += con
        NG50con += 1
        if teoNG50 < testSum:
            NG50 = con
            break
        
    # Checking N90 [bp] [# of Contigs]
    testSum = 0
    N90 = 0
    N90con = 0
    for con in trimmedLength:
        testSum += con
        N90con += 1
        if teoN90 < testSum:
            N90 = con
            break

    # Checking NG90 [bp] [# of Contigs]
    testSum = 0
    NG90 = 0
    NG90con = 0
    for con in trimmedLength:
        if sum < est_genome_size/2:
            break
        testSum += con
        NG90con += 1
        if teoNG90 < testSum:
            NG90 = con
            break
            
    Xkb = 0
    for con in trimmedLength:
        if con > average_gene_size:
            Xkb += 1
            
    print 'Filename: ' + str(contigsMultifasta)        
    print '# of contigs: ' + str(len(contigsLength))
    print 'total length [bp]: ' + str(total)
    print '# of trimmed contigs: ' + str(len(trimmedLength))
    print 'trimmed length [bp]: ' + str(sum)
    print 'GC content [%]: ' + str(GC_cont)
    print 'min. contig [bp]: ' + str(min(trimmedLength))
    print 'median contig [bp]: ' + str(median[0])
    print 'mean contig size [bp]: ' + str(meancontig)
    print 'max. contig [bp]: ' + str(max(trimmedLength))
    print 'N50 Length [bp]: ' + str(N50)
    print '# N50 Contig: ' + str(N50con)
    print 'NG50 Length [bp]: ' + str(NG50)
    print '# NG50 Contig: ' + str(NG50con)
    print 'N90 Length [bp]: ' + str(N90)
    print '# N90 Contig: ' + str(N90con)
    print 'NG90 Length [bp]: ' + str(NG90)
    print '# NG90 Contig: ' + str(NG90con)
    print '# of contigs > ' + str(average_gene_size) + 'bp: ' + str(Xkb)
    
    
    out = open(contigsMultifasta + '_stats.txt', 'w')
    out.write('Filename: ' + str(contigsMultifasta) + '\n')
    out.write('# of contigs: ' + str(len(contigsLength)) + '\n')
    out.write('total length [bp]: ' + str(total) + '\n')
    out.write('# of trimmed contigs: ' + str(len(trimmedLength)) + '\n')
    out.write('trimmed length [bp]: ' + str(sum) + '\n')
    out.write('GC content [%]: ' + str(GC_cont) + '\n')
    out.write('min. contig [bp]: ' + str(min(trimmedLength)) + '\n')
    out.write('median contig [bp]: ' + str(median[0]) + '\n')
    out.write('mean contig size [bp]: ' + str(meancontig) + '\n')
    out.write('max. contig [bp]: ' + str(max(trimmedLength)) + '\n')
    out.write('N50 Length[bp]: ' + str(N50) + '\n')
    out.write('# N50 Contig: ' + str(N50con) + '\n')
    out.write('NG50 Length[bp]: ' + str(NG50) + '\n')
    out.write('# NG50 Contig: ' + str(NG50con) + '\n')
    out.write('N90 Length[bp]: ' + str(N90) + '\n')
    out.write('# N90 Contig: ' + str(N90con) + '\n')
    out.write('NG90 Length[bp]: ' + str(NG90) + '\n')
    out.write('# NG90 Contig: ' + str(NG90con) + '\n')
    out.write('# of contigs > ' + str(average_gene_size) + 'bp: ' + str(Xkb) + '\n')
    out.close()

if __name__ == "__main__":
    contigsMultifasta = sys.argv[1]
    mini = int(sys.argv[2])
    est_genome_size = int(sys.argv[3])
    average_gene_size = int(sys.argv[4])
    
    assembly_stats(contigsMultifasta, mini, est_genome_size, average_gene_size)
