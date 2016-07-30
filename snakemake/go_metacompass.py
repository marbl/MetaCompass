import os,sys,string,subprocess
import argparse

parser = argparse.ArgumentParser(description='snakemake and metacompass params')
parser.add_argument("-t",'--threads', type=int,help='num threads',default=1, nargs='?')
parser.add_argument("-i",'--iterations', type=int, help='num iterations',default=1, nargs='?')
parser.add_argument("-r",'--ref', help='reference genome',default="",nargs='?')
parser.add_argument("-s",'--snakefile', help='snakemake file',default="",nargs='?',required=1)
parser.add_argument("-S",'--Samples', help='sample ids',default="", nargs='*',required=1)
parser.add_argument("-c",'--config', help='config (json) file',default="",nargs='?',required=1)



args = parser.parse_args()
threads = args.threads
iterations = args.iterations
ref = args.ref
snakefile = args.snakefile
config = args.config
samples = args.Samples


print(threads,iterations,ref)


#1. ensure required files are present
if not os.path.exists(snakefile):
    print("snakefile %s not found!"%(snakefile))
    sys.exit(1)

if not os.path.exists(config):
    print("configfile %s not found!"%(config))
    sys.exit(1)

#for reads in samples, check!
#if not os.path.exists

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


#3. for all samples, all iterations, go!
## process one sample at a time, so that we can track input/output easily and run in parallel if we want (qsub for each)
i = 0
while i < iterations:
    for s1 in samples:
        if i == 0:
            ret = subprocess.run("snakemake --snakefile %s --configfile --C final_asm=contig.iter%d.fasta S1=%s"%(snakefile,config,i,s1),shell=True)
        else:
            ret = subprocess.run("snakemake --snakefile %s --configfile --C final_asm=contig.iter%d.fasta ref=./%s.assembly.out/contig.iter%d.fasta S1=%s"%(snakefile,config,i,i,s1,s1),shell=True)
        #this was working, just broke it to handle iterations etc. 
        #ret = subprocess.run("snakemake --snakefile %s --configfile --C final_asm contig.iter%d.fasta"%(snakefile,config),shell=True)
        if ret.returncode != 0:
            print("snakemake command failed; exiting..")
            sys.exit(1)
        i+=1






