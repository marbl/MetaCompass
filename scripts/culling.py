import subprocess
import sys
import argparse
import os
import shutil
import datetime
import json

p_id = os.getpid()

# prints to log file
def logMsg(file, message):
    now = datetime.datetime.now()
    if args.verbose:
        print(message)
    file.write(now.strftime("%Y-%m-%d %H:%M:%S") + ": " + message + "\n")
    file.flush() # make sure things get written out as we go

KMCPATH = "/fs/cbcb-lab/mpop/mpop/devel/bin"  # TODO: use a configuration file to make sure we can find KMC

parser = argparse.ArgumentParser(description="Reference culling")
parser.add_argument('--input', '-i', required=True,
                    metavar='input',
                    help="Metacompass input files")  # input file
parser.add_argument('--reads', '-r', required=True,
                    metavar='reads',
                    help="Metacompass folder containing the forward and reverse fastq files")  # input file
parser.add_argument('--kmer_size', '-ks', required=True,
                    metavar='kmer size',
                    help="Kmer length to be used for analysis")
parser.add_argument('--db', '-d', required=True,
                    metavar='db',
                    help="Directory containing the genomes. Must contain a dataset_catalog.json file.")  # mutation size
parser.add_argument('--outfile', '-o', required=True,
                    metavar='outfolder',
                    help="Output folder")
parser.add_argument('--stop', '-s', default=0.1,
                    metavar='stop', type=float,
                    help="Stopping condition (genomes with fewer specified fraction of k-mers in intersection are discarded")
parser.add_argument('--logfile', '-l', metavar='logfile',
                    help="Log file")
parser.add_argument('--force', '-f', action="store_true",
                    help="Force recreation of k-mer indices")
parser.add_argument('--verbose', '-v', action="store_true",
                    help="Prints more information to stdout. Requires logfile to be set.")

args = parser.parse_args(sys.argv[1:])  # parse command line
print(args)
#Get Kmer Length
kmer_length = "-k" + str(args.kmer_size)

print(kmer_length)
temp = open(str(p_id) + "_reads.txt", "w")
input_files = str(args.input).split(",")
for i in input_files:
    temp.write(i+'\n')
temp.close()

if not os.path.exists(args.reads+ "/"+str(p_id)+"_tmp/"):
    os.mkdir(args.reads+ "/"+str(p_id)+"_tmp/")
tmpdir = args.reads+ "/"+str(p_id)+"_tmp/"

# Create logfile if necessary
logf = None
if args.logfile:
    try:
        logf = open(args.logfile, "w")
    except OSError as err:
        print("Cannot open logfile: {e}".format(e=err))
        exit(1)
    logMsg(logf, "START")

# Create output file
outf = None
try:
    outf = open(args.outfile + '/min_reference_candidates.txt', "w")
except OSError as err:
    print("Cannot open outfile: {e}".format(e=err))
    exit(1)


kmer_overlap_outf=None
try:
    kmer_overlap_outf = open(args.outfile + '/kmer_statistic.csv', "w")
    kmer_overlap_outf.write("Iteration"+','+"Best_Genome"+','+"Best_Kmer_Overlapped"+','+"Best_Genome_Total_Kmer"+','+"Best_Overlaping_Ratio"+"\n")
except OSError as err:
    print("Cannot open outfile: {e}".format(e=err))
    exit(1)


## create kmer files
### Can be run in parallel with the genome-level k-mer computations ###
if args.force or not os.path.isfile(tmpdir + "/reads.base.kmers.kmc_suf"): # need to create file
    if args.logfile:
        logMsg(logf, "Counting k-mers for reads")
    proc = None
    try:
        proc = subprocess.run([KMCPATH + "/kmc", kmer_length, "-ci1", "-hp", "-fq", "@" + str(p_id) + "_reads.txt", tmpdir + "/reads.base.kmers", tmpdir])
    except subprocess.CalledProcessError as err:
        print("Could not run kmc: {e}".format(e=err))
        exit(1)
    except Exception as err:
        print("Could not run kmc: {e}".format(e=err))
        exit(1)

try:
    shutil.copy(tmpdir + "/reads.base.kmers.kmc_suf", tmpdir + "/reads.cur.kmers.kmc_suf")
    shutil.copy(tmpdir + "/reads.base.kmers.kmc_pre", tmpdir + "/reads.cur.kmers.kmc_pre")
except OSError as err:
    print("Could not copy files: {e}".format(e=err))
    exit(1)

genomes = {} # genome information
ngenomes = 0

## create kmer lists for genomes
### This loop should be parallelized ###
genomesF = None
try:
    genomesF = open(args.db + "/dataset_catalog.json", "r")
except OSError as err:
    print("Could not open genomes file: {e}".format(e=err))
    exit(1)

genomeInfo = json.load(genomesF)

#print(genomeInfo['assemblies'])

for g in genomeInfo['assemblies']:
    acc = None
    try:
        acc = g['accession']
    except KeyError:
        continue  # skip genomes without accession

    if args.logfile:
        logMsg(logf, "Handling accession: " + acc)

    for f in g['files']:
        if f["fileType"] == "GENOMIC_NUCLEOTIDE_FASTA":
            temp_path = f["filePath"].split("/")[1]
            if temp_path != "unplaced.scaf.fna":
                if args.logfile:
                    logMsg(logf, "Found file " + f["filePath"])
                genomes[acc] = [acc, args.db + "/" + f["filePath"], tmpdir + "/" + acc + ".kmers",  0]
                break  # only take the first instance

ngenomes = len(genomes)

for g in genomes.keys():
    prefix = genomes[g][0]
    name = genomes[g][1]

    if args.force or not os.path.isfile(tmpdir + "/" + prefix + ".kmers.kmc_suf"):
        if args.logfile:
            logMsg(logf, "Counting k-mers for " + name)
        try:
            proc = subprocess.run([KMCPATH + "/kmc", kmer_length, "-ci1", "-hp", "-fm", name, genomes[g][2], tmpdir])
        except subprocess.CalledProcessError as err:
            print("Could not run kmc: {e}".format(e=err))
            exit(1)

    if args.logfile:
        logMsg(logf, "Retrieving number of k-mers in " + name)

    n = None # number of k-mers
    try:
        proc = subprocess.Popen([KMCPATH + "/kmc_tools", "info", tmpdir + "/" + prefix + ".kmers"],
                                stdout=subprocess.PIPE)

        for line in proc.stdout:
            fields = line.decode('ascii').strip().split(':')
            if fields[0].strip() == "total k-mers":
                n = int(fields[1].strip())
    except subprocess.CalledProcessError as err:
        print("Could not run kmc_tools: {e}".format(e=err))
        exit(1)

    if not n:
        print("kmc_tools did not complete correctly")
        exit(1)

    genomes[prefix][3] = n  # update number of k-mers
##Now loop through a number of iterations until intersection between all genomes and reads is below cutoff

iteration = 0
while ngenomes > 0:
    best_genome=None
    best_kmer_overlapped=0
    best_total_kmer=0
    best_overlaping=0
    iteration += 1
    if args.logfile:
        logMsg(logf, "Starting iteration {i} with {n} genomes".format(i=iteration, n=ngenomes))

    # compute intersection between each genome and reference
    overlaps = {}
    bestprefix = ""
    bestovl = 0
    toProcess = []
    for k in genomes.keys():
        toProcess.append(k)

    currGen = 0; # current genome count
#### This loop should be parallelized
    for prefix in toProcess:  # for each genome
        currGen += 1
        if args.logfile:
            logMsg(logf, str(iteration) + ": Finding read intersection with genome {g}/{n}: ".format(g=currGen, n=ngenomes) + prefix)
        n = -2 # number of k-mers in intersection
        try:
            proc = subprocess.run([KMCPATH + "/kmc_tools", "-hp", "simple", tmpdir + "/" + "reads.cur.kmers",
                                    "-ci1", genomes[prefix][2], "-ci1", "intersect",
                                    tmpdir + "/" + prefix + ".intersect"])

            ## TODO: Need code to directly read number of k-mers from .kmc_pre files
            proc = subprocess.Popen([KMCPATH + "/kmc_tools", "info",  tmpdir + "/" + prefix + ".intersect"],
                                    stdout=subprocess.PIPE)
            
            for line in proc.stdout:
                fields = line.decode('ascii').strip().split(':')
                if fields[0].strip() == "total k-mers":
                    n = int(fields[1].strip())

        except subprocess.CalledProcessError as err:
            print("Could not run kmc_tools: {e}".format(e=err))
            exit(1)

        if n < -1:
            print("kmc_tools did not complete correctly")
            exit(1)

        overlaps[prefix] = n
        if args.logfile:
            logMsg(logf, "Found {n} k-mers".format(n=n))
        print(float(overlaps[prefix]), float(genomes[prefix][3]), "Check Here")
        if float(overlaps[prefix]) / float(genomes[prefix][3]) < args.stop: # genome has too few k-mers selected
            if args.logfile:
                logMsg(logf, "Genome {p} has too few kmers {n} out of {m}".format(p=prefix,
                                                                                  n=overlaps[prefix],
                                                                                  m = genomes[prefix][3]))
            del genomes[prefix]
            del overlaps[prefix]
            ngenomes -= 1
        elif overlaps[prefix] > bestovl: # update the best
            ## IMPORTANT: Currently we use the genome that most overlaps the reads.
            ## An alternative implementation would select the genome that is best covered
            ## by the reads.
            ## TODO: Must try both approaches and compare results.
            overlap_p = float(overlaps[prefix]) / float(genomes[prefix][3])
            if args.logfile:
                logMsg(logf, "Genome {p} now is best with {k} kmers".format(p=prefix, k=overlaps[prefix]))
                
            bestovl = overlaps[prefix]
            bestprefix = prefix
            best_genome=prefix
            best_kmer_overlapped=overlaps[prefix]
            best_total_kmer=float(genomes[prefix][3])
            best_overlaping=overlap_p
            bestovl = overlaps[prefix]
            bestprefix = prefix

    if bestprefix != "": # if bestprefix is blank, we're done
        if args.logfile:
            logMsg(logf, "Genome selected is " + bestprefix)
            kmer_overlap_outf.write(str(iteration)+','+str(best_genome)+','+str(best_kmer_overlapped)+','+str(best_total_kmer)+','+str(best_overlaping)+"\n")
        outf.write(bestprefix + "\n")


        #subtract best genome's k-mers from the read kmers
        proc = subprocess.run([KMCPATH + "/kmc_tools", "-hp", "simple", tmpdir + "/" + "reads.cur.kmers", "-ci1",
                                 genomes[bestprefix][2], "-ci1", "kmers_subtract",
                                 tmpdir + "/" + "reads.tmp.kmers"])
        ngenomes -= 1
        del genomes[bestprefix]

        # update the current k-mer set
        shutil.copy(tmpdir + "/reads.tmp.kmers.kmc_suf", tmpdir + "/reads.cur.kmers.kmc_suf")
        shutil.copy(tmpdir + "/reads.tmp.kmers.kmc_pre", tmpdir + "/reads.cur.kmers.kmc_pre")


if args.logfile:
    logMsg(logf, "END")

os.system("rm -r " + tmpdir)
os.system("rm " + str(p_id) + "_reads.txt")
