#modified from TIPP_reference_package/src/get_genes_per_acc_mc.py
import os
import sys
from glob import glob
from itertools import groupby

refdir = sys.argv[1]
prokaryote_file = sys.argv[2]
# refdir = '/fs/cbcb-scratch/zbowen/metacompass/ref_pkg_generation_V2'

def fasta_iter(fasta_name):
    fh = open(fasta_name)
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        header = next(header)[1:].strip()
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq
    fh.close()

dirname = refdir + '/RefSeq_V2_db/marker_index'
files = glob(dirname+'/*/*faa')

count = {}
for f in files:
	fiter = fasta_iter(f)
	seen = {}
	acc = {}
	for ff in fiter:
		accname = '_'.join(ff[0].split('_')[0:2])
		seen[ff[0]] =1
		if accname in acc:
			print ("multiple genes of acc in file", accname, f)
		acc[accname] = 1
	

	fiter = fasta_iter(f.replace('faa','fna'))
	for ff in fiter:
		if ff[0] not in seen:
			print ("print not found in faa", f)
	for accname in acc:
		if accname in count:
			count[accname] += 1
		else:
			count[accname] = 1

## read in taxonomy 
taxa = {}
with open (refdir + '/RefSeq_V2_db/data/taxonomy.info') as f:
	for line in f:
		val = line.strip().split(',')
		taxa[val[0]] = line.strip()

refseq = {}
with open(prokaryote_file) as f:
	for line in f:
		if line.startswith('#'):
			continue
		val =line.strip().split('\t')
		print("Line: " + line)
		print(val)
		tokeep = [val[1], val[0], val[15]]
		refseq[val[18]] = tokeep

fw = open(refdir + '/RefSeq_V2_db/taxinfo.tsv', 'w')
fw.write('#accession\tnum_genes\trefseq_taxid\tspecies_taxid\torganism_name\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tassembly_level\n')
successes=0
fails=0
for accname in count:
	if accname in taxa: # 2023 EDIT SOME accname WERE NOT PRESENT
		taxaval = taxa[accname].split(',')
		taxaval = [str(x) if x != '' else "NA" for x in taxaval ]
		refseqval = refseq[accname]
		arr2print = [accname, str(count[accname]), str(refseqval[0])] + [str(refseqval[1])]+ taxaval[2:] + [str(refseqval[2])]
		fw.write('\t'.join(arr2print)+ '\n')
		successes+=1
	else:
		print(accname + " was missing!!")
		fails+=1

print("Script completed with " + str(successes) + " successes and " + str(fails) + " failures.")
# Usage python get_genes_per_acc.py  output_folder genes_per_acc.txt
