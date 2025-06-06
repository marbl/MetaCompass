#!/usr/bin/env python

import sys

def main():
	
	# Read in marker gene id, file with breadth of coverage for each genome and output file name
	gene = sys.argv[1]
	coverage_file = sys.argv[2]
	out_file = sys.argv[3]
	marker_index = sys.argv[4]

	# Make a dictionary of representative sequences (key) and the other sequences that are included in the cluster (val)
	markerGeneClusterDict=dict()
	with open( marker_index + '/' + gene + '/' + gene + '_clustered.clusters' ) as f:
		for num, line in enumerate(f):
			# The first line in the file is a header and can be skipped
			if num == 0:
				continue
			val = line.strip().split('\t')
			# Each line in the file contains the marker sequence id (seq) and the representative sequence of the cluster (repSeq)
			seq = val[0]
			repSeq = val[2]
			# For each representative seq, store the other sequences within the same cluster in the dictionary
			markerGeneClusterDict.setdefault(repSeq,[]).append(seq)

	# Make a dictionary of marker genes (key) and which genome they belong to (value)
	markerGeneGenomeDict=dict()
	with open( marker_index + '/' + gene + '/' + gene + '_genome.tsv' ) as f:
		for num, line in enumerate(f):
			# First line in the file is a header that should be ignored
			if num == 0:
				continue
			val = line.strip().split('\t')
			seq = val[0]
			genome = val[2]
			markerGeneGenomeDict[seq]=genome

	# For each marker gene sequence that is covered, write out to a file indicating 
	# the marker gene seq, the reference genome it comes from, and that it is covered (indicated by 1)
	f1 = open(out_file, 'w')

	with open(coverage_file) as f:
		for num, line in enumerate(f):
			if num==0:
				continue
			val = line.strip().split('\t')
			seq = val[0]
			cov = "1"

			# Get list of marker genes within cluster that share this sequence
			matches = markerGeneClusterDict[seq]
			for item in matches:
				genome = markerGeneGenomeDict[item]
				f1.write(str(genome)+"\t"+str(gene)+"\t"+cov+"\n")

	f1.close()
	
if __name__ == '__main__':
	main()
