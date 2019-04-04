#!/usr/bin/env python3

import os,sys,string,subprocess,signal,shutil
#psutil
import argparse

mcdir = sys.path[0]
parser = argparse.ArgumentParser(description='snakemake and metacompass params')

group1 = parser.add_argument_group('required')
group1.add_argument("-c",'--config', help='config (json) file, set read length etc',default="",nargs='?',required=0,type=str)
group1.add_argument("-S",'--Samples', help='Provide file with fq reads (1 file per line)',default="", nargs='?',required=0,type=str)
group1.add_argument("-P",'--paired', help='Provide comma separated list of paired reads (r1.1.fq,r1.2.fq)',default="", nargs='?',required=0,type=str)
group1.add_argument("-U",'--unpaired', help='Provide comma separated list of unpaired reads (r1.fq,r2.fq,r3.fq)',default="", nargs='?',required=0,type=str)

group5 = parser.add_argument_group("metacompass")
group5.add_argument("-d",'--db', help='marker gene database directory',default="", nargs='?',type=str)
#remove or keep?
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
mincov = args.mincov
readlen=args.readlen
clobber = args.clobber
unlock = args.unlock
threads = args.threads
iterations = args.iterations
ref = args.ref
mfilter = args.filter
keepoutput = args.keepoutput

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
prefix = outdir
i=1


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
        
        
for s1 in allsamples[0:1]:  
    s1id = s1.split(os.sep)[-1].split(".")[0]
    if sampleid != "NA":
        s1id = sampleid
        
print("%s %s %d"%(prefix,s1id,i-1))
        
if os.path.exists("%s/%s.%d.assembly.out"%(prefix,s1id,i-1)):
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
        if os.stat("%s/%s.0.assembly.out/%s.mc.sam.unmapped.1.fq"%(prefix,s1id,s1id)).st_size!= 0:
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