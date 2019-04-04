#!/usr/bin/env python3

import os,sys,string,subprocess,signal,shutil
#psutil
import argparse

mcdir = sys.path[0]
parser = argparse.ArgumentParser(description='snakemake and metacompass params')

group1 = parser.add_argument_group('required')
#group1.add_argument("-s",'--snakefile', help='metacompass rules file',default="",nargs='?',required=1,type=str)
group1.add_argument("-c",'--config', help='config (json) file, set read length etc',default="",nargs='?',required=0,type=str)
group1.add_argument("-S",'--Samples', help='Provide file with fq reads (1 file per line)',default="", nargs='?',required=0,type=str)
group1.add_argument("-P",'--paired', help='Provide comma separated list of paired reads (r1.1.fq,r1.2.fq)',default="", nargs='?',required=0,type=str)
group1.add_argument("-U",'--unpaired', help='Provide comma separated list of unpaired reads (r1.fq,r2.fq,r3.fq)',default="", nargs='?',required=0,type=str)

group5 = parser.add_argument_group("metacompass")
group5.add_argument("-d",'--db', help='marker gene database directory',default="", nargs='?',type=str)
group5.add_argument("-i",'--iterations', type=int, help='num iterations',default=1, nargs='?')
group5.add_argument("-r",'--ref', help='reference genomes',default="NA",nargs='?')
group5.add_argument("-p",'--pickref', help='depth or breadth',default="breadth",nargs='?')
group5.add_argument("-m",'--mincov', help='min coverage to assemble',default="3",nargs='?',type=int)
group5.add_argument("-g",'--minctglen', help='min contig length',default="300",nargs='?',type=int)
group5.add_argument("-l",'--readlen', help='max read length',default="100",nargs='?',type=int)
group5.add_argument("-f",'--filter',help='filter recruited genomes with mash (experimental)',default=False,required=0, type=float)
group2 = parser.add_argument_group('output')
group2.add_argument("-b",'--clobber', help='clobber output directory (if exists?)',default=False,required=0,action='store_true')
group2.add_argument("-o",'--outdir', help='output directory? (cwd default)',default="./", nargs='?',type=str,required=1)
group2.add_argument("-e",'--sampleid', help='sample id (fq prefix is default)',default="NA", nargs='?',type=str,required=0)
group2.add_argument("-v",'--verbose', help='verbose',default=False,required=0,action='store_true')
group2.add_argument("-k",'--keepoutput', help='keep all output generated (default is to delete all but final fasta files)',default=False,required=0,action='store_true')

group3 = parser.add_argument_group('performance')
group3.add_argument("-t",'--threads', type=int,help='num threads',default=1, nargs='?')
group3.add_argument("-q",'--qsub', help='',default="", nargs='?',required=0)

group4 = parser.add_argument_group('snakemake')
group4.add_argument("-F",'--Force', help='force snakemake to rerun',default=False,required=0,action='store_true')
group4.add_argument("-u",'--unlock',help='unlock snakemake locks',default=False, required=0,action='store_true')


args = parser.parse_args()
minctglen = args.minctglen
db = str(args.db)
if db != "" and not os.path.isdir(db):
    print("provided marker gene database directory %s does not exist; try again"%(db))
    sys.exit(1)

mincov = args.mincov
readlen=args.readlen
clobber = args.clobber
unlock = args.unlock
threads = args.threads
iterations = args.iterations
ref = args.ref
mfilter = args.filter
keepoutput = args.keepoutput

#if args.filter:
#why todd?
#    #empirically determined on datasets with known truth, right way to do this is with contains operation
#    mfilter = 0.26

#snakefile = args.snakefile
config = args.config
samples = args.Samples.replace(" ","")
unpaired = args.unpaired
paired = args.paired
sampleid = args.sampleid
qsub = args.qsub
force = args.Force
verbose = args.verbose
outdir = args.outdir
pickref = args.pickref
prefix = "."
retry = False

if not os.path.exists(outdir):
    os.makedirs(outdir)
    prefix = outdir
else:
    if os.path.exists(outdir) and not clobber:
        print("ERROR: specified output directory %s exists! please remove first, or run with --clobber"%(outdir))
        sys.exit(1)
    elif os.path.exists(outdir) and force:
        os.system("rm -rf %s/*"%(outdir))
        os.system("mkdir %s"%(outdir))
        prefix = outdir
    elif os.path.exists(outdir):
        prefix = outdir
#1. ensure required files are present
#if not os.path.exists(snakefile):
#    print("ERROR: snakefile %s not found!"%(snakefile))
#    sys.exit(1)

if config == "":
    config = "%s/snakemake/config.json"%(mcdir)    
elif not os.path.exists(config):
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
print("checking for dependencies (Bowtie2, Blast, kmermask, Snakemake, etc)")

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

#print("bedtools--->",end="")
#ret = subprocess.call("which bedtools",shell=True)
#if ret == 0:
#    print("[OK]")
#else:
#    print("[FAIL]")
#    sys.exit()

if mfilter < 1.0:
    print("mash--->",end="")
    ret = subprocess.call("which mash",shell=True)
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
        allsamples.append(line.strip())


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
    else:
        allfiles.append(unpaired)
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

    if len(allfiles) != 0:
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
            #ret = subprocess.call("snakemake -r --verbose --config ref=%s.0.assembly.out/mc.refseq.fna --snakefile %s/snakemake/metacompass.iter0.unpaired.py --configfile %s --unlock"%(s1id,mcdir,config),shell=True)
            cmd_ret="snakemake -r --verbose --reason --unlock --cores %d -a --configfile %s --config prefix=%s sample=%s pickref=breadth reference=%s mcdir=%s iter=%d length=%d mincov=%d minlen=%d mfilter=%f nthreads=%d ref=%s.0.assembly.out/mc.refseq.fna"%(threads,config,prefix,s1id,ref,mcdir,i,readlen,mincov,minctglen,mfilter,threads,s1id)
            cmd_ret += " reads="
            for fqfile in allsamples:
                cmd_ret += str(fqfile)+","
                cmd_ret = cmd_ret[:-1]
            #todo:fix to work in all cases, add r1,r2,ru
            if paired != "":
                cmd_ret += " r1=%s r2=%s"%(paired.split(",")[0],paired.split(",")[1])
            if unpaired != "":
                cmd_ret += " ru=%s"%(unpaired)
            if unpaired != "" and paired =="":
                cmd_ret +=" --snakefile %s/snakemake/metacompass.iter0.unpaired.py"%(mcdir)
            elif paired != "" and unpaired =="":
                cmd_ret +=" --snakefile %s/snakemake/metacompass.iter0.paired.py"%(mcdir)
            elif paired != "" and unpaired !="":
                cmd_ret +=" --snakefile %s/snakemake/metacompass.iter0.py"%(mcdir)
#               elif samples =="":
#
            ret = subprocess.call(cmd_ret,shell=True)
#            ret = subprocess.call("snakemake -r --verbose --reason --cores %d -a --configfile %s --config prefix=%s sample=%s pickref=breadth reference=%s mcdir=%s iter=%d length=%d mincov=%d minlen=%d mfilter=%f nthreads=%d ref=%s.0.assembly.out/mc.refseq.fna --snakefile %s/snakemake/metacompass.iter0.unpaired.py reads=%s --unlock"%(threads,config,prefix,s1id,ref,mcdir,i,readlen,mincov,minctglen,mfilter,threads,s1id,mcdir,allsamples),shell=True)
# if paired != "":cmd += " r1=%s r2=%s"%(paired.split(",")[0],paired.split(",")[1])
        if i == 0:
            ret = 0
            #todo: fix to work with diff types of reads?
            if ref != "NA":
                cmd = "snakemake --verbose --reason --cores %d -a --configfile %s --config prefix=%s sample=%s pickref=breadth reference=%s mcdir=%s iter=%d length=%d mincov=%d minlen=%d mfilter=%f nthreads=%d"%(threads,config,prefix,s1id,ref,mcdir,i,readlen,mincov,minctglen,mfilter,threads)

            else:        
                cmd = "snakemake --verbose --reason --cores %d -a --configfile %s --config prefix=%s sample=%s pickref=breadth reference=%s/%s.%d.assembly.out/mc.refseq.fna mcdir=%s iter=%d length=%d mincov=%d minlen=%d mfilter=%f nthreads=%d "%(threads,config,prefix,s1id,prefix,s1id,i,mcdir,i,readlen,mincov,minctglen,mfilter,threads)

            cmd += " reads="
            for fqfile in allsamples:
                cmd += str(fqfile)+","
                #cmd = cmd[:-1]
            if paired != "":
                cmd += " r1=%s r2=%s"%(paired.split(",")[0],paired.split(",")[1])
            if unpaired != "":
                cmd += " ru=%s"%(unpaired)
  
            if ref != "NA":
                if unpaired != "" and paired =="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.ref.unpaired.py"%(mcdir)
                elif paired != "" and unpaired =="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.ref.paired.py"%(mcdir)
                elif paired != "" and unpaired !="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.ref.py"%(mcdir)
                #todo: fix to work with diff types of reads..
                elif samples =="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.ref.py"%(mcdir)
                
            else: 
                if unpaired != "" and paired =="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.unpaired.py"%(mcdir)
                elif paired != "" and unpaired =="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.paired.py"%(mcdir)
                elif paired != "" and unpaired !="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.py"%(mcdir)
                #todo: fix to work with diff types of reads..
                elif samples =="":
                    cmd += " --snakefile %s/snakemake/metacompass.iter0.py"%(mcdir)

            if verbose:
#iredunpaired reads -U?
                cmd += " --verbose"
            if retry:
                cmd += " --rerun-incomplete"
            else:
                cmd += " --ignore-incomplete"

            cmd += " --prioritize join_contigs"

            if len(qsub) > 0:
                cmd += " --cluster %s"%(qsub)
            else:
                try:
                    ret = subprocess.Popen(cmd,shell=True)
                    ret.communicate()
                except KeyboardInterrupt:
                    os.killpg(ret.pid,signal.SIGKILL)  
                except :
                    ret.returncode = 0
                    break                        
        else:
            ret = 0
            if ref != "NA":
                cmd = "snakemake --cores %d -a --configfile %s --config prefix=%s sample=%s reference=%s/%s.%d.assembly.out/contigs.pilon.fasta mcdir=%s iter=%d pickref=%s length=%d mincov=%d minlen=%d nthreads=%d"%(threads,config,prefix,s1id,prefix,s1id,i-1,mcdir,i,pickref,readlen,mincov,minctglen,threads)
            else:
                cmd = "snakemake --cores %d -a --configfile %s --config prefix=%s sample=%s reference=%s/%s.%d.assembly.out/contigs.pilon.fasta mcdir=%s iter=%d pickref=%s length=%d mincov=%d minlen=%d nthreads=%d"%(threads,config,prefix,s1id,prefix,s1id,i-1,mcdir,i,pickref,readlen,mincov,minctglen,threads)

            cmd += " reads="
            for fqfile in allsamples:
                cmd += fqfile+","
            cmd = cmd[:-1]         
            if paired != "":
                cmd += " r1=%s r2=%s"%(paired.split(",")[0],paired.split(",")[1])
 
            
            cmd += " --snakefile %s/snakemake/metacompass.py"%(mcdir)
            
            if verbose:
                cmd += " --verbose"
            if retry:
                cmd += " --rerun-incomplete"
            else:
                cmd += " --ignore-incomplete"
            cmd += " --prioritize pilon_contigs"
            if len(qsub) > 0:
                cmd += " --cluster %s"%(qsub)
                
            else:
                try:
                    ret = subprocess.Popen(cmd,shell=True)
                    ret.communicate()
                except KeyboardInterrupt:
                    os.killpg(ret.pid,signal.SIGKILL)  
                except:
                    #command finished but ran in background and timed out? h
                    ret.returncode = 0
                    break
                    #os.killpg(ret.pid,signal.SIGKILL)  
        if ret.returncode != 0:
            print("ERROR: snakemake command failed; exiting..")
            os.system("touch %s/%s.0.assembly.out/run.fail"%(prefix,s1id))
            sys.exit(1)
        i+=1


if os.path.exists("%s/%s.%d.assembly.out/contigs.final.fasta"%(prefix,s1id,i-1)):
    #cleanup output
    if not os.path.exists("%s/metacompass_output"%(prefix)):
        os.mkdir("%s/metacompass_output"%(prefix))
    if not os.path.exists("%s/metacompass_logs"%(prefix)):
        os.mkdir("%s/metacompass_logs"%(prefix))
    os.system("cp %s/%s.0.assembly.out/contigs.final.fasta %s/metacompass_output/metacompass.final.ctg.fa"%(prefix,s1id,prefix))
    #os.system("cp %s/metacompass_output/metacompass.final.ctg.fa %s/."%(prefix,prefix))
    os.system("mv %s/metacompass*tsv %s/metacompass_output/."%(prefix,prefix))
    os.system("mv %s/*.log %s/metacompass_logs/."%(prefix,prefix))
    os.system("mv %s/%s.0.assembly.out/*.log %s/metacompass_logs/."%(prefix,s1id,prefix))
    os.system("mv %s/%s.0.assembly.out/coverage.txt %s/metacompass_output/metacompass.genomes_coverage.txt"%(prefix,s1id,prefix))
    
    #os.system("cp %s/%s.0.assembly.out/contigs.pilon.fasta.fixed %s/metacompass_output/metacompass.only.ctg.fa"%(prefix,s1id,prefix))
    if mfilter < 1.0:
        if os.path.exists("%s/%s.merged.fq.mash.out.ids"%(prefix,s1id)):
            os.system("cp %s/%s.merged.fq.mash.out.ids  %s/metacompass_output/metacompass.recruited.ids"%(prefix,s1id,prefix))
            os.system("cp %s/%s.0.assembly.out/mc.refseq.filt.fna  %s/metacompass_output/metacompass.recruited.fa"%(prefix,s1id,prefix))
    else:
        if os.path.exists("%s/%s.%d.assembly.out/mc.refseq.ids"%(prefix,s1id,i-1)):
            os.system("cp %s/%s.0.assembly.out/mc.refseq.ids  %s/metacompass_output/metacompass.recruited.ids"%(prefix,s1id,prefix))
            os.system("cp %s/%s.0.assembly.out/mc.refseq.fna  %s/metacompass_output/metacompass.recruited.fa"%(prefix,s1id,prefix))

###create tsv output
#    os.system("grep '>' %s/%s.0.assembly.out/contigs.fasta| tr -d '>'|sed 's/_[0-9]\+ / /g' >%s/.tmp"%(prefix,s1id,prefix))
#    os.system("grep '>' %s/tr -d '>'|sed 's/_[0-9]\+ / /g' >%s/.tmp"%(prefix,s1id,prefix))


    if not keepoutput:
        print("Cleaning up files..")
        shutil.rmtree("%s/%s.0.assembly.out/"%(prefix,s1id))
        os.system("rm %s/*.fq "%(prefix))
        os.system("rm %s/*.fastq "%(prefix))
        os.system("rm %s/*.fasta "%(prefix))
        #os.system("rm %s/*.out "%(prefix))
        if os.path.exists("%s/%s.%d.assembly.out/mc.refseq.ids"%(prefix,s1id,i-1)):
            os.system("rm %s/*.ids "%(prefix))
        if os.path.exists("%s/%s.merged.fq.mash.out.ids"%(prefix,s1id)):
            os.system("rm %s/%s.merged.fq.mash.out.ids "%(prefix,s1id))
        #os.system("mv %s/*.log %s/metacompass_logs/."%(prefix,prefix))
    else:
        #if referece selection
        if os.stat("%s/%s.0.assembly.out/mc.refseq.fna"%(prefix,s1id)).st_size!= 0:
            os.makedirs("%s/%s.0.assembly.out/reference_selection_output"%(prefix,s1id), exist_ok=True)
            os.system("mv %s/%s.0.assembly.out/mc.* %s/%s.0.assembly.out/reference_selection_output"%(prefix,s1id,prefix,s1id))
#SRS016585.merged.fq.mash.
            os.system("mv %s/%s.0.assembly.out/*merged.fq.mash* %s/%s.0.assembly.out/reference_selection_output"%(prefix,s1id,prefix,s1id))
        else:
            os.system("rm %s/%s.0.assembly.out/mc.refseq.fna"%(prefix,s1id))
        os.makedirs("%s/%s.0.assembly.out/pilon_output"%(prefix,s1id), exist_ok=True)
        os.makedirs("%s/%s.0.assembly.out/assembly_output"%(prefix,s1id), exist_ok=True)
        ########moving mapping after best strata    
        os.system("mv %s/%s.0.assembly.out/contigs.fasta %s/%s.0.assembly.out/assembly_output"%(prefix,s1id,prefix,s1id))
        os.system("mv %s/%s.0.assembly.out/%s.sam %s/%s.0.assembly.out/assembly_output"%(prefix,s1id,s1id,prefix,s1id))
########testingstart
        os.system("mv %s/%s.0.assembly.out/%s.sam.all %s/%s.0.assembly.out/assembly_output"%(prefix,s1id,s1id,prefix,s1id))
########testing end
        if os.path.exists("%s/%s.0.assembly.out/selected_maps.sam"%(prefix,s1id)):
            os.system("mv %s/%s.0.assembly.out/selected_maps.sam %s/%s.0.assembly.out/assembly_output"%(prefix,s1id,prefix,s1id))
        os.system("mv %s/%s.0.assembly.out/*buildcontigs* %s/%s.0.assembly.out/assembly_output"%(prefix,s1id,prefix,s1id))
        os.system("mv %s/%s.0.assembly.out/sorted*.bam* %s/%s.0.assembly.out/pilon_output"%(prefix,s1id,prefix,s1id))
        os.system("mv %s/%s.0.assembly.out/*.pilon* %s/%s.0.assembly.out/pilon_output/"%(prefix,s1id,prefix,s1id))
 #if unmmaped reads#
        if os.stat("%s/%s.0.assembly.out/%s.mc.sam.unmapped.1.fq"%(prefix,s1id,s1id)).st_size!= 0 or os.stat("%s/%s.0.assembly.out/%s.mc.sam.unmapped.2.fq"%(prefix,s1id,s1id)).st_size!= 0 :
            if os.path.exists("%s/%s.0.assembly.out/unmapped_reads"%(prefix,s1id)):
                os.system("rm -rf %s/%s.0.assembly.out/unmapped_reads"%(prefix,s1id))
            os.mkdir("%s/%s.0.assembly.out/unmapped_reads"%(prefix,s1id))
            os.system("mv %s/%s.0.assembly.out/*sam.unmapped* %s/%s.0.assembly.out/unmapped_reads"%(prefix,s1id,prefix,s1id))
            os.system("mv %s/%s.0.assembly.out/*.megahit %s/%s.0.assembly.out/megahit_output"%(prefix,s1id,prefix,s1id))
        else:
            os.system("rm -rf %s/%s.0.assembly.out/*.megahit"%(prefix,s1id))
            os.system("rm %s/%s.0.assembly.out/*mc.sam.unmapped*"%(prefix,s1id))
        #provide mapped reads too? the bam is in pilon_output 
        os.mkdir("%s/%s.0.assembly.out/mapped_reads"%(prefix,s1id))
        os.system("mv %s/%s.0.assembly.out/*.mc*.sam* %s/%s.0.assembly.out/mapped_reads"%(prefix,s1id,prefix,s1id))
        os.system("rm %s/%s.0.assembly.out/*index "%(prefix,s1id))
        os.system("rm %s/%s.0.assembly.out/*.bt2 "%(prefix,s1id))
#####only for testing ???????
        #os.system("rm %s/%s.0.assembly.out/*.sam "%(prefix,s1id))
        #os.system("rm %s/%s.0.assembly.out/*.bam* "%(prefix,s1id))
        #os.system("mv %s/%s.0.assembly.out/*.sam "%(prefix,s1id))
        ###???os.system("mv %s/%s.0.assembly.out/*.sam %s/%s.0.assembly.out/mapped_reads"%(prefix,s1id,prefix,s1id))
        
        
        os.system("rm %s/%s.0.assembly.out/*.fasta "%(prefix,s1id))
        #os.system("rm %s/%s.0.assembly.out/*.log "%(prefix,s1id))
        os.system("rm %s/*.f*q "%(prefix))
        os.system("rm %s/*.fasta "%(prefix))
        if os.path.exists("%s/%s.0.assembly.out"%(prefix,s1id)):
            if os.path.exists("%s/intermediate_files"%(prefix)):
                os.system("rm -rf %s/intermediate_files"%(prefix))
            shutil.move("%s/%s.0.assembly.out"%(prefix,s1id) , "%s/intermediate_files"%(prefix))
        if os.path.exists("%s/%s.merged.fq.mash.out.ids"%(prefix,s1id)):
            os.system("rm %s/%s.merged.fq.mash.out.ids "%(prefix,s1id))
    #os.system("touch %s/%s.0.assembly.out/run.ok"%(prefix,s1id))
   # os.system("samtools faidx metacompass.final.ctg.fa")
       #|cut -f1-2> %s/.tmp2"%(prefix))
    print("MetaCompass finished succesfully!")
else:
    os.system("touch %s/%s.0.assembly.out/run.fail"%(prefix,s1id))
    print("MetaCompass run failed. See Log files for more info")

