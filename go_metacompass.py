#!/usr/bin/env python3
import os,sys,string,subprocess,signal,shutil,argparse
from subprocess import PIPE, run
#psutil

print("MetaCompass metagenome assembler version 2.0.0 by Victoria Cepeda (vcepeda@cs.umd.edu)\n")

mcdir = sys.path[0]
#1 READ/PARSE COMMAND LINE ARGUMENTS
#1.1 PROCESS COMMAND LINE ARGUMENTS

parser = argparse.ArgumentParser(description='snakemake and metacompass params')

group1 = parser.add_argument_group('required')
#group1.add_argument("-s",'--snakefile', help='metacompass rules file',default="",nargs='?',required=1,type=str)
group1.add_argument("-c",'--config', help='config (json) file, set read length etc',default="",nargs='?',required=0,type=str)
#group1.add_argument("-S",'--Samples', help='Provide file with fq reads (1 file per line)',default="", nargs='?',required=0,type=str)
group1.add_argument("-1",'--forward', help='Provide comma separated list of forward paired-end reads',default="", nargs='?',required=0,type=str)
group1.add_argument("-2",'--reverse', help='Provide comma separated list of reverse paired-end reads',default="", nargs='?',required=0,type=str)
group1.add_argument("-U",'--unpaired', help='Provide comma separated list of unpaired reads (r1.fq,r2.fq,r3.fq)',default="", nargs='?',required=0,type=str)

group5 = parser.add_argument_group("metacompass")
group5.add_argument("-r",'--ref', help='reference genomes',default="NA",nargs='?')
group5.add_argument("-s",'--refsel', help='reference selection [tax/all]',default="tax",nargs='?')
group5.add_argument("-p",'--pickref', help='depth or breadth',default="breadth",nargs='?')
group5.add_argument("-m",'--mincov', help='min coverage to assemble',default="3",nargs='?',type=int)
group5.add_argument("-g",'--minctglen', help='min contig length',default="300",nargs='?',type=int)
group5.add_argument("-l",'--readlen', help='max read length',default="100",nargs='?',type=int)

group2 = parser.add_argument_group('output')
group2.add_argument("-b",'--clobber', help='clobber output directory (if exists?)',default=False,required=0,action='store_true')
group2.add_argument("-o",'--outdir', help='output directory? (cwd default)',default="metacompass_output", nargs='?',type=str,required=1)
group2.add_argument("-k",'--keepoutput', help='keep all output generated (default is to delete all but final fasta files)',default=False,required=0,action='store_true')

group3 = parser.add_argument_group('performance')
group3.add_argument("-t",'--threads', type=int,help='num threads',default=1, nargs='?')
group3.add_argument("-y",'--memory', type=int,help='memory',default=8, nargs='?',required=1)

group4 = parser.add_argument_group('snakemake')
group4.add_argument('--Force', help='force snakemake to rerun',default=False,required=0,action='store_true')
group4.add_argument('--unlock',help='unlock snakemake locks',default=False, required=0,action='store_true')
group4.add_argument('--nolock', help='remove stale locks',default=False,required=0,action='store_true')
group4.add_argument('--verbose', help='verbose',default=False,required=0,action='store_true')
group4.add_argument('--reason', help='reason',default=False,required=0,action='store_true')
#--dryrun, -n :Do not execute anything.
group4.add_argument('--dryrun', help='dryrun',default=False,required=0,action='store_true')
#-cores

args = parser.parse_args()

#1.2 PROCESS COMMAND LINE ARGUMENTS

minctglen = args.minctglen
mincov = args.mincov
readlen=args.readlen
refsel=args.refsel
threads = args.threads
memory = args.memory
ref = args.ref
keepoutput = args.keepoutput
#snakefile = args.snakefile
config = args.config
#samples = args.Samples.replace(" ","")
unpaired = args.unpaired
fpaired = args.forward
rpaired = args.reverse
outdir = args.outdir
pickref = args.pickref
clobber = args.clobber
unlock = args.unlock
nolock = args.nolock
force = args.Force
verbose = args.verbose
reason = args.reason
dryrun = args.dryrun

#4. Check for existence output file:        
if os.path.exists(outdir):
    if not clobber and not force:#test this
        if os.path.exists("%s/run.ok"%(outdir)):
            print("ERROR: Output dir (%s) exists and contains a previous, successful run. Please specify alternate output directory or force run with --force"%(outdir))
            sys.exit(1)
        elif os.path.exists("%s/run.fail"%(outdir)):
            print("ERROR: Output dir (%s) exists and contains a previous, failed run. Please specify alternate output directory or force run with --force"%(outdir))
            sys.exit(1)
        else:
            print("ERROR: Output dir (%s) exists. Please specify alternate output directory or run with --clobber or --unlock"%(outdir))#force
            sys.exit(1)
    elif force:#test this
        os.system("rm -rf %s/*"%(outdir))
        os.system("mkdir %s"%(outdir))
        
else:
    os.makedirs(outdir)

optsfile="%s/opts.txt"%(outdir)
outfile = open(optsfile,'a')
print ( "minctglen: %s" % (minctglen), file= outfile)       
print ( "mincov: %s" % (mincov), file= outfile)       
print ( "readlen: %s" % (readlen), file= outfile)       
print ( "refsel: %s" % (refsel), file= outfile)       
print ( "threads: %s" % (threads), file= outfile)       
print ( "memory: %s" % (memory), file= outfile)       
print ( "ref: %s" % (ref), file= outfile)       
print ( "keepoutput: %s" % (keepoutput), file= outfile)       
print ( "unpaired: %s" % (unpaired), file= outfile)       
print ( "fpaired: %s" % (fpaired), file= outfile)       
print ( "rpaired: %s" % (rpaired), file= outfile)       
print ( "outdir: %s" % (outdir), file= outfile)       
print ( "pickref: %s" % (pickref), file= outfile)       
print ( "clobber: %s" % (clobber), file= outfile)       
print ( "unlock: %s" % (unlock), file= outfile)       
print ( "nolock: %s" % (nolock), file= outfile)       
print ( "force: %s" % (force), file= outfile)       
print ( "verbose: %s" % (verbose), file= outfile)       
print ( "reason: %s" % (reason), file= outfile)       
print ( "dryrun: %s" % (dryrun), file= outfile)       

#cmd="snakemake--prioritize join_contigs "
cmd="snakemake -T --printshellcmds "

#--printshellcmds Print out the shell commands that will be executed.
if verbose:
    cmd += " --verbose"
if unlock:
    cmd += " --unlock"
if nolock:
    cmd += " --nolock"    
if force:
    cmd += " --Force"    
if reason:
    cmd += " --reason"
if dryrun:
    cmd += " --dryrun"
    
#1.3 Checking config file
print("confirming config file exists...", file= outfile)

if config == "":
    config = "%s/snakemake/config.json"%(mcdir)
if not os.path.exists(config):
    print("ERROR: configfile %s not found!"%(config))
    sys.exit(1)
print ("config: %s" % (config), file= outfile)       

#1.4 Checking reference genomes
#fix
if ref != "NA":
    print("confirming file containing reference genomes exists...",file=outfile)
    if not os.path.exists(ref):
        print("ERROR: reference genome file %s not found!"%(ref))
        sys.exit(1)
    else:
    #   os.system("cp %s %s/%s"%(ref,outdir,ref.split(os.sep)[-1]))
        print("Reference genome file: %s" % (ref), file= outfile)
        
    
#2. Check for dependencies:
print("checking for assembly dependencies (Snakemake,Bowtie2,Samtools)",file=outfile)

#2.1 Check for assembly dependencies:
command = ['which', 'snakemake']
result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
print("Snakemake",file=outfile)
print(result.stdout,file=outfile)
if result.returncode != 0:
    print("Snakemake not found")
    sys.exit(1)

command = ['which', 'bowtie2']
result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
print("Bowtie2:",file=outfile)
print(result.stdout,file=outfile)
if result.returncode != 0:
    print("Bowtie2 not found. Bowtie2 v>=2.2.9 required")
    sys.exit(1)
#else:
    #version=str(os.system("bowtie2 --version|head -n1|cut -f3 -d ' '|cut -f1-3 -d '.'|head -n1"))
    #version=version.split("\n")[0]
    #required=str("2.2.9")
    #print ("version %s required %s" % (version, required))
    #if version < required:#"2.2.9":
    #    print("Bowtie2 2.2.9 or higher required")
    #    sys.exit(1)
command = ['which', 'samtools']
result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
print("Samtools:",file=outfile)
print(result.stdout,file=outfile)
if result.returncode != 0:
    print("Samtools not found")
    sys.exit(1)

#2.2 Check for reference selection dependencies:
print("checking for reference selection dependencies (Blast, kmer-mask, mash)",file=outfile)

if ref == "NA":
    command = ['which', 'blastn']
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    print("BLAST:",file=outfile)
    print(result.stdout,file=outfile)
    if result.returncode != 0:
        print("Blast not found")
        sys.exit(1)

    command = ['which', 'kmer-mask']
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    print("Kmer-mask:",file=outfile)
    print(result.stdout,file=outfile)
    if result.returncode != 0:
        print("Kmer-mask not found")
        sys.exit(1)
        
    command = ['which', 'mash']
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    print("Mash:",file=outfile)
    print(result.stdout,file=outfile)
    if result.returncode != 0:
        print("Mash not found")
        sys.exit(1)

allsamples=[]
if fpaired != "":
    for file in fpaired.split(","):
        if not os.path.exists(file):
            print("ERROR: could not locate file %s"%(file))
            sys.exit()
        else:
            allsamples.append(file)
    
if rpaired != "":
    for file in rpaired.split(","):
        if not os.path.exists(file):
            print("ERROR: could not locate file %s"%(file))
            sys.exit()
        else:
            allsamples.append(file)
paired=""
if fpaired != "" and rpaired != "":
   paired="yes"
    
#-U single end
if unpaired != "":
    for file in unpaired.split(","):
        if not os.path.exists(file):
            print("ERROR: could not locate file %s"%(file))
            sys.exit(1)
        else:
            allsamples.append(file)

##########################################################################
#reference_selection vs assembly only

cmd += " --cores %d -a --configfile %s --config outdir=%s pickref=%s mcdir=%s length=%d mincov=%d minlen=%d nthreads=%d memory=%d refsel=%s"%(threads,config,outdir,pickref,mcdir,readlen,mincov,minctglen,threads,memory,refsel)

cmd += " reads="
print("ALL READS found:",file=outfile)
for fqfile in allsamples:
    cmd += str(fqfile)+","
    print("%s"%(fqfile),file=outfile)
cmd = cmd[:-1]#remove last comma
    
if fpaired != "" and rpaired !="":
    cmd += " r1=%s r2=%s"%(fpaired,rpaired)
if unpaired != "":
    cmd += " ru=%s"%(unpaired)
  
if ref != "NA":
    print("REFERENCE genome file provided. Reference Selection step will be skipped.")
    
    cmd += " reference=%s" %(ref)
    if unpaired != "" and paired =="":
        cmd += " --snakefile %s/snakemake/metacompass.ref.unpaired.py"%(mcdir)
    elif paired != "" and unpaired =="":
        cmd += " --snakefile %s/snakemake/metacompass.ref.paired.py"%(mcdir)
    elif paired != "" and unpaired !="":
        cmd += " --snakefile %s/snakemake/metacompass.ref.py"%(mcdir)
    #elif samples =="":
    #    cmd += " --snakefile %s/snakemake/metacompass.ref.py"%(mcdir)
else: 
    cmd += " reference=reference_selection/mc.refseq.fna"
    if unpaired != "" and paired =="":
        cmd += " --snakefile %s/snakemake/metacompass.unpaired.py"%(mcdir)
    elif paired != "" and unpaired =="":
        cmd += " --snakefile %s/snakemake/metacompass.paired.py"%(mcdir)
    elif paired != "" and unpaired !="":
        cmd += " --snakefile %s/snakemake/metacompass.py"%(mcdir)

#print("Snakemake command:")
#print("%s\n"%(cmd))
print("Snakemake command:",file=outfile)
print("%s\n"%(cmd),file=outfile)
#RUN SNAKEMAKE!!    
try:
    ret = subprocess.Popen(cmd,shell=True)
    ret.communicate()
except KeyboardInterrupt:
#    print('Interrupted')
#    print("ERROR: snakemake command failed; exiting..")
#    os.system("touch %s/run.fail"%(outdir))
    try:
        sys.exit(1)
    except SystemExit:
#        print("ERROR: snakemake command failed; exiting..")
#        os.system("touch %s/run.fail"%(outdir))
        os._exit(1)
        #os.killpg(ret.pid,signal.SIGKILL)
        
    #os.killpg(ret.pid,signal.SIGKILL)
    #print("ERROR: SIGKILL")
    #sys.exit(1)  
except:
    ret.returncode = 1
if ret.returncode != 0:
    print("ERROR: snakemake command failed; exiting..")
    os.system("touch %s/run.fail"%(outdir))
    try:
        sys.exit(1)
    except SystemExit:
        #print("ERROR: snakemake command failed; exiting..")
        #os.system("touch %s/run.fail"%(outdir))
        os._exit(1)
else:
    if dryrun:
        sys.exit(0)
    #5 CLEANING output files
    if os.path.exists("%s/assembly/coverage.txt"%(outdir)):
        os.system("mv %s/assembly/coverage.txt %s/metacompass_output/metacompass.genomes_coverage.txt"%(outdir,outdir))
    print("Cleaning up files..")
    if os.path.exists("%s/intermediate_files"%(outdir)):
        os.system("rm -rf %s/intermediate_files"%(outdir))
    #reference_Selection
    if os.path.exists("%s/reference_selection"%(outdir)):
        os.system("rm %s/reference_selection/*.fastq "%(outdir))
        os.system("rm %s/reference_selection/mc.blastn* "%(outdir))
        os.system("rm %s/reference_selection/contigs_clusters "%(outdir))
        os.system("rm %s/reference_selection/*msh* "%(outdir))
        os.system("mv %s/reference_selection/*.log %s/logs "%(outdir,outdir))    
    #assembly
    if os.path.exists("%s/assembly/mc.index.1.bt2"%(outdir)):
        os.system("rm %s/assembly/mc.index* "%(outdir))
    if os.path.exists("%s/assembly/mc.sam"%(outdir)):
        os.system("mv %s/assembly/mc.sam %s/mapped_reads/"%(outdir,outdir))
    #error_correction
    if os.path.exists("%s/error_correction/mc.index.1.bt2"%(outdir)):
        os.system("rm %s/error_correction/mc.index* "%(outdir))
    if os.path.exists("%s/error_correction/mc.sam"%(outdir)):
        os.system("rm %s/error_correction/mc.sam "%(outdir))
    if os.path.exists("%s/error_correction/mc.sam.bam"%(outdir)):
        os.system("rm %s/error_correction/mc.sam.bam "%(outdir))
    if os.path.exists("%s/error_correction/mc_unpaired.sam"%(outdir)):
        os.system("rm %s/error_correction/mc_unpaired.sam "%(outdir))
    if os.path.exists("%s/error_correction/mc_unpaired.sam.bam"%(outdir)):
        os.system("rm %s/error_correction/mc_unpaired.sam.bam "%(outdir))
    if os.path.exists("%s/error_correction/mc.sam.unmapped.1.fq"%(outdir)):
        os.system("mv %s/error_correction/mc.sam.unmapped.1.fq %s/unmapped_reads/"%(outdir,outdir))
    if os.path.exists("%s/error_correction/mc.sam.unmapped.2.fq"%(outdir)):
        os.system("mv %s/error_correction/mc.sam.unmapped.2.fq %s/unmapped_reads/"%(outdir,outdir))
    if os.path.exists("%s/error_correction/mc.sam.unmapped.u.fq"%(outdir)):
        os.system("mv %s/error_correction/mc.sam.unmapped.u.fq %s/unmapped_reads/"%(outdir,outdir)) 
    #simplified output  
    if not keepoutput:
        os.system("rm -rf %s/intermediate_files"%(outdir))
        os.system("rm -rf %s/error_correction"%(outdir))
        if os.path.exists("%s/reference_selection"%(outdir)):
            os.system("rm -rf %s/reference_selection"%(outdir))
        os.system("rm -rf %s/assembly"%(outdir))
        if os.path.exists("%s/unmapped_reads"%(outdir)):
            os.system("rm -rf %s/unmapped_reads"%(outdir))
        if os.path.exists("%s/mapped_reads"%(outdir)):
            os.system("rm -rf %s/mapped_reads"%(outdir))

