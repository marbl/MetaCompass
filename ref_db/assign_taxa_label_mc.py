#modified from TIPP_reference_package/src/assign_taxa_label.py
import os
import sys

from ete3 import NCBITaxa
ncbi = NCBITaxa()

acc_list = {}
with open(sys.argv[1]) as f:
	for line in f:
		val = line.strip()
		acc_list[val] = 1

levels = {'superkingdom':1, 'phylum':1, 'class':1, 'order':1, 'family':1, 'genus':1, 'species':1}
with open(sys.argv[2]) as f:
	for line in f:
		if line.startswith('#'):
			continue
		val = line.strip().split('\t')
		print("Line: " + line)
		print(val)
		accno = val[18]
		if accno not in acc_list:
			continue
		taxid = val[1]
		final_map = {}
		for leveliter in levels:
			final_map[leveliter] = ""
		try:
			lineage = ncbi.get_lineage(taxid)
			rank_map = ncbi.get_rank(lineage)
			name_map = ncbi.get_taxid_translator(lineage)
			
			for taxid_iter in rank_map:
				if rank_map[taxid_iter] in levels:
					final_map[rank_map[taxid_iter]] = name_map[taxid_iter]
		except:
			print ('#'+accno + ',' + taxid +': not found')

		print (accno + ',' + taxid +',' + final_map['superkingdom'] + ',' + final_map['phylum'] + ',' + final_map['class'] + ',' + final_map['order']+ ',' + final_map['family']+ ',' + final_map['genus']+ ',' + final_map['species'])

#Usage python assign_taxa_label.py  accession_list.txt assembly_summary_refseq.txt > taxonomy.info




