import os,sys,string

r1 = sys.argv[1]
g1 = sys.argv[2]
outf = open(sys.argv[3],'w')
g1f = open(sys.argv[2],'r')
mfilter = float(sys.argv[4])
#print(mfilter)
if mfilter < 1.0:
    os.system("mash sketch -r -p 64 -k 20 -s 100000 %s"%(r1))
    os.system("mash dist -i -p 64 %s.msh %s -t > %s.mash.out"%(r1,g1,r1))
    os.system("awk '{if($2<%f) print $1}' %s.mash.out > %s.mash.out.ids"%(mfilter,r1,r1))
else:
    os.system("cp %s %s.mash.out.ids"%(g1.replace(".fna",".ids"),r1))

idsf = open("%s.mash.out.ids"%(r1),'r')
idsf.readline()

ids = []
for line in idsf.readlines():
    ids.append(line.replace("\n","").split("\t")[0])

data = g1f.read().split(">")[1:]
seq_dict = {}
for seq in data:
    hdr,fasta = seq.split("\n",1)
    seq_dict[hdr.replace("\n","").split(" ")[0]] = fasta

for id in ids:
    outf.write(">%s\n"%(id))
    outf.write(seq_dict[id])

outf.close()
g1f.close()
#1. open up mash output
#2. keep ids that have distance <0.3
#3. filter reffile to remove genomes with distance >0.2
