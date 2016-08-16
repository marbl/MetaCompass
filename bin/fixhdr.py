import os,sys,string

f1 = open(sys.argv[1],'r')
outf = open(sys.argv[1]+".fna",'w')

ctgs = f1.read().split(">")[1:]
for ctg in ctgs:
    hdr,fasta = ctg.split("\n",1)
    hdr = hdr.rsplit("_",1)[0]
    outf.write(">"+hdr+"\n")
    width = 80
    i = 0
    while i+width < len(fasta):
        outf.write(fasta[i:i+width]+"\n")
        i+=width
    if len(fasta[i:]) > 0:
        outf.write(fasta[i:]+"\n")
 

f1.close()
outf.close()
