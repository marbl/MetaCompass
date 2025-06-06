#!/bin/bash
maps=${1}

declare -A id2Sp
while IFS=$'\t' read accession	num_gene refseq_taxid	species_taxidorganism_name	Kingdom	Phylum	Class	Order	Family	Genus	Species	assembly_level
do
    id2Sp[${accession}]=${Species}
done < $maps

echo ${id2Sp[GCA_900082345.1]}
