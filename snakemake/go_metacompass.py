import os,sys,string,subprocess,signal
#psutil
import argparse

parser = argparse.ArgumentParser(description='snakemake and metacompass params')

group1 = parser.add_argument_group('required')
group1.add_argument("-s",'--snakefile', help='snakemake file',default="",nargs='?',required=1,type=str)
group1.add_argument("-c",'--config', help='config (json) file',default="",nargs='?',required=1,type=str)
group1.add_argument("-S",'--Samples', help='Provide file with fq reads (1 file per line)',default="", nargs='?',required=1,type=str)

group5 = parser.add_argument_group("metacompass")
group5.add_argument("-i",'--iterations', type=int, help='num iterations',default=1, nargs='?')
group5.add_argument("-r",'--ref', help='reference genomes',default="NA",nargs='?')

group2 = parser.add_argument_group('output')
group2.add_argument("-o",'--outdir', help='output directory? (cwd default)',default=".", nargs='?',type=str,required=0)
group2.add_argument("-d",'--sampleid', help='sample id (fq prefix is default)',default="NA", nargs='?',type=str,required=0)
group2.add_argument("-v",'--verbose', help='verbose',type=bool, default=False,nargs='?',required=0)

group3 = parser.add_argument_group('performance')
group3.add_argument("-t",'--threads', type=int,help='num threads',default=1, nargs='?')
group3.add_argument("-q",'--qsub', help='',default="", nargs='?',required=0)

group4 = parser.add_argument_group('snakemake')
group4.add_argument("-F",'--force', help='force snakemake to rerun',type=bool, default="False",nargs='?',required=0)
group4.add_argument("-u",'--unlock', type=bool, help='unlock',default=False, nargs='?')

args = parser.parse_args()
unlock = args.unlock
threads = int(args.threads)
iterations = args.iterations
ref = args.ref
snakefile = args.snakefile
config = args.config
samples = args.Samples.replace(" ","")
sampleid = args.sampleid
qsub = args.qsub
force = args.force
verbose = args.verbose
outdir = args.outdir
if not os.path.exists(outdir):
    prefix = os.getcwd
else:
    prefix = outdir

if force:
   if os.path.exists("Sample1.fasta"):
      os.system("rm Sample1.fasta")
   if os.path.exists("Sample1.marker.match.1.fastq"):
      os.system("rm Sample1.marker.match.1.fastq")
   os.system("rm -rf ./Sample1.*.assembly.out/")
#print(qsub)
#print(threads,iterations,ref,unlock)


#1. ensure required files are present
if not os.path.exists(snakefile):
    print("snakefile %s not found!"%(snakefile))
    sys.exit(1)

if not os.path.exists(config):
    print("configfile %s not found!"%(config))
    sys.exit(1)


if ref != "NA":
    print("confirming file containing reference genomes exists..")
    if not os.path.exists(samples):
        print("reference genome file %s not found!"%(ref))
        sys.exit(1)
    else:
        print("[OK]")
    
#for reads in samples, check!
#if not os.path.exists
print("confirming sample file exists..")
if not os.path.exists(samples):
    print("sample file (-S) %s not found!"%(samples))
    sys.exit(1)
else:
    print("[OK]")
#2. check for snakemake, bowtie2
print("checking for dependencies (Bowtie2, Blast, kmermask, Snakemake)")

if len(qsub) > 0:

    print("qsub--->",end="")
    ret = subprocess.call("which qsub",shell=True)
    if ret == 0:
        print("[OK]")
        qsub += " --jobs 4000"
    else:
        print("[FAIL]")
        qsub=""


print("Bowtie2--->",end="")
ret = subprocess.call("which bowtie2",shell=True)
if ret == 0:
    print("[OK]")
else:
    print("[FAIL]")
    sys.exit()

print("Blast+--->",end="")
ret = subprocess.call("which blastn",shell=True)
if ret == 0:
    print("[OK]")
else:
    print("[FAIL]")
    sys.exit()

print("kmer-mask--->",end="")
ret = subprocess.call("which kmer-mask",shell=True)
if ret == 0:
    print("[OK]")
else:
    print("[FAIL]")
    sys.exit()

print("Snakemake--->",end="")
ret = subprocess.call("which snakemake",shell=True)
if ret == 0:
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
        s1id = s1.split(os.sep)[-1].split(".")[0]
        if sampleid != "NA":
            s1id = sampleid

        #print(s1id)
            
        if unlock:
            ret = subprocess.call("snakemake -r --verbose --config ref=%s.0.assembly.out/mc.refseq.fna --snakefile %s --configfile %s --unlock"%(s1id,snakefile,config),shell=True)
        if i == 0:
            ret = 0
            if ref != "NA":
                cmd = "snakemake --cores %d --configfile %s --ignore-incomplete --config prefix=%s sample=%s reads=%s ref=%s iter=%d --snakefile %s"%(threads,config,prefix,s1id,s1,ref,i,snakefile)
            else:        
                cmd = "snakemake --cores %d --configfile %s --config prefix=%s sample=%s reads=%s ref=%s.%d.assembly.out/mc.refseq.fna iter=%d --snakefile %s"%(threads,config,prefix,s1id,s1,s1id,i,i,snakefile)

            if verbose:
                cmd += " --verbose"
            #if force:
            #    cmd += " -F"
            if len(qsub) > 0:
                cmd += " --cluster %s"%(qsub)
            else:
                try:
                    ret = subprocess.Popen(cmd,shell=True)
                    ret.communicate()
                except KeyboardInterrupt:
                    os.killpg(ret.pid,signal.SIGKILL)  
                except :
                    os.killpg(ret.pid,signal.SIGKILL)  
                    

        else:
            ret = 0
            if ref != "NA":
                cmd = "snakemake --cores %d --configfile %s --ignore-incomplete --config prefix=%s sample=%s reads=%s ref=%s.%d.assembly.out/contigs.fasta iter=%d --snakefile %s"%(threads,config,prefix,s1id,s1,s1id,i-1,i,snakefile)
            else:
                cmd = "snakemake --cores %d --configfile %s --config prefix=%s sample=%s reads=%s ref=%s.%d.assembly.out/contigs.fasta iter=%d --snakefile %s"%(threads,config,prefix,s1id,s1,s1id,i-1,i,snakefile)


            if verbose:
                cmd += " --verbose"
            #if force:
            #    cmd += " -F"
            if len(qsub) > 0:
                cmd += " --cluster %s"%(qsub)
                
            else:
                try:
                    ret = subprocess.Popen(cmd,shell=True)
                    ret.communicate()
                except KeyboardInterrupt:
                    os.killpg(ret.pid,signal.SIGKILL)  
                except:
                    os.killpg(ret.pid,signal.SIGKILL)  

        if ret.returncode != 0:
            print("snakemake command failed; exiting..")
            sys.exit(1)
        i+=1






