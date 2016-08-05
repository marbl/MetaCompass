import os,sys,string,subprocess
import argparse

parser = argparse.ArgumentParser(description='snakemake and metacompass params')
parser.add_argument("-t",'--threads', type=int,help='num threads',default=1, nargs='?')
parser.add_argument("-i",'--iterations', type=int, help='num iterations',default=1, nargs='?')
parser.add_argument("-u",'--unlock', type=bool, help='unlock',default=False, nargs='?')
parser.add_argument("-r",'--ref', help='reference genome',default="",nargs='?')
parser.add_argument("-s",'--snakefile', help='snakemake file',default="",nargs='?',required=1)
#parser.add_argument("-I",'--SampleIds', help='sample ids',default="", nargs='*',required=1)
parser.add_argument("-S",'--Samples', help='',default="", nargs='?',required=1)
parser.add_argument("-c",'--config', help='config (json) file',default="",nargs='?',required=1)

#parser.add_argument("-F",'--force', help='force snakemake to rerun',type=bool, default="",nargs='?',required=1)



args = parser.parse_args()
unlock = args.unlock
threads = args.threads
iterations = args.iterations
ref = args.ref
snakefile = args.snakefile
config = args.config
samples = args.Samples.replace(" ","")


print(threads,iterations,ref,unlock)


#1. ensure required files are present
if not os.path.exists(snakefile):
    print("snakefile %s not found!"%(snakefile))
    sys.exit(1)

if not os.path.exists(config):
    print("configfile %s not found!"%(config))
    sys.exit(1)


#for reads in samples, check!
#if not os.path.exists
print("confirming sample file exists..")
print(samples)
if not os.path.exists(samples):
    print("sample file (-S) %s not found!"%(samples))
    sys.exit(1)
else:
    print("[OK]")
#2. check for snakemake, bowtie2
print("checking for dependencies (Bowtie2, Snakemake)")

print("Bowtie2--->",end="")
ret = subprocess.run("which bowtie2",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
if ret.returncode == 0:
    print("[OK]")
else:
    print("[FAIL]")
    sys.exit()

print("Snakemake--->",end="")
ret = subprocess.run("which snakemake",stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
if ret.returncode == 0:
    print("[OK]")
else:
    print("[FAIL]")
    sys.exit()

#3. process sample file
#read through sample file and get sample paths & ids
#SE or PE,interleaved or non-interleaved
allsamples = []
samplesf = open(samples,'r')
for line in samplesf.readlines():
    allsamples.append(line.replace("\n",""))

#CURRENTLY only single fastq file is supported
#s1.fq
#s2.fq

#3. for all samples, all iterations, go!
## process one sample at a time, so that we can track input/output easily and run in parallel if we want (qsub for each)


i = 0
while i < iterations:
    for s1 in allsamples:
        s1id = s1.split(".")[0]
        if unlock:
            ret = subprocess.run("snakemake -r --verbose --config ref=%s.0.assembly.out/mc.refseq.fna --snakefile %s --configfile %s --unlock"%(s1id,snakefile,config),shell=True)
        if i == 0:
            #%s.%s.assembly.out
            ret = subprocess.run("snakemake -r --verbose --snakefile %s --configfile %s --config sample=%s reads=%s ref=%s.%d.assembly.out/mc.refseq.fna iter=%d"%(snakefile,config,s1id,s1,s1id,i,i),shell=True)
        else:
            ret = subprocess.run("snakemake -r --verbose --snakefile %s --configfile %s --config sample=%s reads=%s ref=%s.%d.assembly.out/contigs.fasta iter=%d"%(snakefile,config,s1id,s1,s1id,i-i,i),shell=True)
        #this was working, just broke it to handle iterations etc. 
        #ret = subprocess.run("snakemake --snakefile %s --configfile --C final_asm contig.iter%d.fasta"%(snakefile,config),shell=True)
        if ret.returncode != 0:
            print("snakemake command failed; exiting..")
            sys.exit(1)
        i+=1






