import os,sys,string

f1 = open(sys.argv[1],'r')
outf = open(sys.argv[1]+".fna",'w')
hdr = f1.readline()
fasta = f1.read()
hdr = hdr.rsplit("_",1)[0]

outf.write(hdr+"\n")
outf.write(fasta)
f1.close()
