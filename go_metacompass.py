#!/usr/bin/env python3

import os,sys,string,subprocess,signal
#psutil
import argparse

mcdir = sys.path[0]
parser = argparse.ArgumentParser(description='snakemake and metacompass params')

group1 = parser.add_argument_group('required')
group1.add_argument("-s",'--snakefile', help='metacompass rules file',default="",nargs='?',required=1,type=str)
group1.add_argument("-c",'--config', help='config (json) file, set read length etc',default="",nargs='?',required=1,type=str)
group1.add_argument("-S",'--Samples', help='Provide file with fq reads (1 file per line)',default="", nargs='?',required=0,type=str)
group1.add_argument("-P",'--paired', help='Provide comma separated list of paired reads (r1.1.fq,r1.2.fq)',default="", nargs='?',required=0,type=str)
group1.add_argument("-U",'--unpaired', help='Provide comma separated list of unpaired reads (r1.fq,r2.fq,r3.fq)',default="", nargs='?',required=0,type=str)

group5 = parser.add_argument_group("metacompass")
group5.add_argument("-i",'--iterations', type=int, help='num iterations',default=1, nargs='?')
group5.add_argument("-r",'--ref', help='reference genomes',default="NA",nargs='?')
group5.add_argument("-p",'--pickref', help='coverage',default="breadth",nargs='?')

group2 = parser.add_argument_group('output')
group2.add_argument("-o",'--outdir', help='output directory? (cwd default)',default="./", nargs='?',type=str,required=1)
group2.add_argument("-d",'--sampleid', help='sample id (fq prefix is default)',default="NA", nargs='?',type=str,required=0)
group2.add_argument("-v",'--verbose', help='verbose',default=False,required=0,action='store_true')

group3 = parser.add_argument_group('performance')
group3.add_argument("-t",'--threads', type=int,help='num threads',default=1, nargs='?')
group3.add_argument("-q",'--qsub', help='',default="", nargs='?',required=0)

group4 = parser.add_argument_group('snakemake')
group4.add_argument("-F",'--force', help='force snakemake to rerun',default=False,required=0,action='store_true')
group4.add_argument("-m",'--retry', help='retry/resume a failed run',default=False,required=0,action='store_true')
group4.add_argument("-u",'--unlock',help='unlock snakemake locks',default=False, required=0,action='store_true')

args = parser.parse_args()
unlock = args.unlock
threads = int(args.threads)
iterations = args.iterations
ref = args.ref
snakefile = args.snakefile
config = args.config
samples = args.Samples.replace(" ","")
unpaired = args.unpaired
paired = args.paired
sampleid = args.sampleid
qsub = args.qsub
force = args.force
verbose = args.verbose
outdir = args.outdir
retry = args.retry
pickref = args.pickref
prefix = ""

if not os.path.exists(outdir):
    os.makedirs(outdir)
    prefix = outdir
else:
    if os.path.exists(outdir) and not force:
        print("ERROR: specified output directory %s exists! please remove first, or run with --force"%(outdir))
        sys.exit(1)
    elif os.path.exists(outdir) and force:
        prefix = outdir

#1. ensure required files are present
if not os.path.exists(snakefile):
    print("ERROR: snakefile %s not found!"%(snakefile))
    sys.exit(1)

if not os.path.exists(config):
    print("ERROR: configfile %s not found!"%(config))
    sys.exit(1)


if ref != "NA":
    print("confirming file containing reference genomes exists..")
    if not os.path.exists(ref):
        print("ERROR: reference genome file %s not found!"%(ref))
        sys.exit(1)
    else:
        os.system("cp %s %s/%s"%(ref,prefix,ref.split(os.sep)[-1]))
        print("[OK]")
    
#for reads in samples, check!
#if not os.path.exists
#print("confirming sample file exists..")
#if "," not in samples and not os.path.exists(samples):
#    print("ERROR: sample file (-S) %s not found!"%(samples))
#    sys.exit(1)
#else:
#    print("[OK]")
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
sout = open("out.txt",'w')
ret = subprocess.call("bowtie2 --version",stdout=sout,shell=True)
if ret == 0:# and "2.2.9" in sout:
    #print(stdout)
    #sys.exit(1)
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
#3 paths: sample file, -P, or -U
#-S or -P & -U
allsamples = []
if samples != "" and (paired != "" or unpaired != ""):
    print("ERROR: Must specific -S or -P/-U; please correct and relaunch")
    sys.exit()

if samples != "" and (paired == "" and unpaired == ""):
    samplesf = open(samples,'r')
    for line in samplesf.readlines():
        allsamples.append(line.replace("\n",""))


elif samples == "" and paired != "" and unpaired == "":
    #only paired end
    if not "," in paired:
        print("ERROR: --paired reads need to be provided as -P r1.1.fq,r1.2.fq")
        sys.exit()
    elif not os.path.exists(paired.split(",")[0]) or not os.path.exists(paired.split(",")[1]):
        print("ERROR: could not locate --paired fq files %s,%s"%(paired.split(",")[0],paired.split(",")[1]))
        sys.exit()
    allsamples = paired.split(",")


elif samples == "" and paired == "" and unpaired != "":
    #only unpaired
    allfiles = []
    if "," in unpaired:
        allfiles = unpaired.split(",")
    for ufile in allfiles:
        if not os.path.exists(ufile):
            print("ERROR: could not locate --unpaired file %s"%(ufile))
            sys.exit()

    allsamples = allfiles


elif samples == "" and paired != "" and unpaired != "":
    #only paired end
    if not "," in paired:
        print("ERROR: --paired reads need to be provided as -P r1.1.fq,r1.2.fq")
        sys.exit()
    elif not os.path.exists(paired.split(",")[0]) or not os.path.exists(paired.split(",")[1]):
        print("ERROR: could not locate --paired fq files %s,%s"%(paired.split(",")[0],paired.split(",")[1]))
        sys.exit()

    allsamples = paired.split(",")

    #only unpaired
    allfiles = []
    if "," in unpaired:
        allfiles = unpaired.split(",")
    else:
        allfiles = [unpaired]
    for ufile in allfiles:
        if not os.path.exists(ufile):
            print("ERROR: could not locate --unpaired file %s"%(ufile))
            sys.exit()

    allsamples.extend(allfiles)


#allsamples = []
#if "," not in samples:
#    samplesf = open(samples,'r')
#    for line in samplesf.readlines():
#        allsamples.append(line.replace("\n",""))
#else:
#    allsamples = [samples]
 

#CURRENTLY only single fastq file is supported
#s1.fq
#s2.fq

#3. for all samples, all iterations, go!
## process one sample at a time, so that we can track input/output easily and run in parallel if we want (qsub for each)


i = 0
isok = False
while i < iterations:
    for s1 in allsamples[0:1]:
        #s1a = ""
        #s1b = ""
        #if "," in s1: 
        #    s1a,s1b = s1.split(",")
        #else:
        #    s1a = s1
        #if not os.path.exists("%s"%(s1a)):# or not os.path.exists("%s"%(s1b)):
        #    print("ERROR: Cannot locate file %s within input file %s. Please verify and restart"%(s1a,args.Samples))
        #    sys.exit(1)
        s1id = s1.split(os.sep)[-1].split(".")[0]
        if sampleid != "NA":
            s1id = sampleid

        
        if force:
            if os.path.exists("%s.fasta"%(s1id)):
                os.system("rm %s.fasta"%(s1id))
            if os.path.exists("%s.marker.match.1.fastq"%(s1id)):
                os.system("rm %s.marker.match.1.fastq"%(s1id))
            if os.path.exists("%s.marker.match.2.fastq"%(s1id)):
                os.system("rm %s.marker.match.2.fastq"%(s1id))
            os.system("rm -rf ./%s.*.assembly.out/"%(s1id))
        elif os.path.exists("%s/%s.0.assembly.out/run.ok"%(prefix,s1id)):
            #run finished ok, don't allow to clobber
            print("ERROR: Output dir (%s/%s.0.assembly.out) exists and contains a previous, successful run. Please specify alternate output directory or force run with --force"%(prefix,s1id))
            sys.exit(1)
        elif retry and os.path.exists("%s/%s.0.assembly.out/run.fail"%(prefix,s1id)):
            #run finished ok, don't allow to clobber
            print("Output dir (%s/%s.0.assembly.out) exists and contains a previous, failed run. Attempting to resume failed run.."%(prefix,s1id))
        elif not retry and os.path.exists("%s/%s.0.assembly.out/run.fail"%(prefix,s1id)):
            print("ERROR: Output dir (%s/%s.0.assembly.out) exists and contains a previous, failed run. If you'd like to retry/resume this run, specify: --retry"%(prefix,s1id))
            sys.exit(1)
        if unlock:
            ret = subprocess.call("snakemake -r --verbose --config ref=%s.0.assembly.out/mc.refseq.fna --snakefile %s --configfile %s --unlock"%(s1id,snakefile,config),shell=True)
        if i == 0:
            ret = 0
            if ref != "NA":
                cmd = "snakemake --cores %d -a --configfile %s --config prefix=%s sample=%s pickref=breadth ref=%s mcdir=%s iter=%d "%(threads,config,prefix,s1id,ref,mcdir,i)
            else:        
                cmd = "snakemake --cores %d -a --configfile %s --config prefix=%s sample=%s pickref=breadth ref=%s/%s.%d.assembly.out/mc.refseq.fna mcdir=%s iter=%d "%(threads,config,prefix,s1id,prefix,s1id,i,mcdir,i)

            cmd += " reads="
            for fqfile in allsamples:
                cmd += fqfile+","
             
            if paired != "":
                cmd += " r1=%s r2=%s"%(paired.split(",")[0],paired.split(",")[1])
 
            cmd += " --snakefile %s"%(snakefile)

            if verbose:
                cmd += " --verbose"
            if retry:
                cmd += " --rerun-incomplete"
            else:
                cmd += " --ignore-incomplete"
            #cmd += " --mcdir %s"%(mcdir)

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
                cmd = "snakemake --cores %d -a --configfile %s --config prefix=%s sample=%s ref=%s/%s.%d.assembly.out/contigs.pilon.fasta mcdir=%s iter=%d pickref=%s "%(threads,config,prefix,s1id,prefix,s1id,i-1,mcdir,i,pickref)
            else:
                cmd = "snakemake --cores %d -a --configfile %s --config prefix=%s sample=%s ref=%s/%s.%d.assembly.out/contigs.pilon.fasta mcdir=%s iter=%d pickref=%s "%(threads,config,prefix,s1id,prefix,s1id,i-1,mcdir,i,pickref)

            cmd += " reads="
            for fqfile in allsamples:
                cmd += fqfile+","
             
            if paired != "":
                cmd += " r1=%s r2=%s"%(paired.split(",")[0],paired.split(",")[1])
 
            cmd += " --snakefile %s"%(snakefile)
            
            if verbose:
                cmd += " --verbose"
            if retry:
                cmd += " --rerun-incomplete"
            else:
                cmd += " --ignore-incomplete"
            #if force:
            #    cmd += " -F"
            #cmd += " --mcdir %s"%(mcdir)
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
            print("ERROR: snakemake command failed; exiting..")
            os.system("touch %s/%s.0.assembly.out/run.fail"%(prefix,s1id))
            sys.exit(1)
        i+=1


if os.path.exists("%s/%s.%d.assembly.out/contigs.fasta"%(prefix,s1id,i-1)):
    os.system("touch %s/%s.0.assembly.out/run.ok"%(prefix,s1id))
    print("MetaCompass finished succesfully!")
else:
    os.system("touch %s/%s.0.assembly.out/run.fail"%(prefix,s1id))
    print("MetaCompass run failed. See Log files for more info")

