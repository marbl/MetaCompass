import os,sys,string
import operator


pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])
print ("%s" % (pathbin))
r1 = sys.argv[1]
g1 = sys.argv[2]
results = sys.argv[3]
#print ("%s" % (outf))
#g1f = open(sys.argv[2],'r')
mfilter = float(sys.argv[4])
#shell:"/usr/bin/time %/bin/mash-Linux64-v2.0/mash sketch -p 12 -k 21 -s 1000 -i {input.g1} -o %/%.0.assembly.out/sketch; /usr/bin/time %/bin/mash-Linux64-v2.0/mash screen -w -p 12 {output.mashfile} {input.r1} >mash_screen.tab; cut -f4 mash_screen.tab >mc.refseq.filt.ids %bin/extractseq "
#if mfilter < 1.0:
cmd="%s/bin/mash-Linux64-v2.0/mash sketch -p 64 -k 21 -s 1000 -i %s"%(pathbin,g1)
print ("%s" % (cmd))
os.system(cmd)

cmd="%s/bin/mash-Linux64-v2.0/mash screen -w -p 64 %s.msh %s > %s.mash.tab"%(pathbin,g1,r1,r1)
print ("%s" % (cmd))
os.system(cmd)
#mfliter=0000000001
cmd="awk '{if($4<%s) print $5}' %s.mash.tab > %s.mash.ids"%(mfilter,r1,r1)
print ("%s" % (cmd))
os.system(cmd)
cmd="%s/bin/extractSeq %s %s.mash.ids > %s"%(pathbin,g1,r1,results)
#cmd="/bin/extractSeq {1} {2}.mash.ids  outf".format((pathbin,g1,r1))
print ("%s" % (cmd))
os.system(cmd)

#1. run mash screen
#2. keep ids that have identity >.7
#3. filter reffile to remove genomes with distance >.7
