#!/usr/bin/env python

import sys
import pandas as pd
import os


def main():
    try:

        # Determine percent of marker genes that are required to be present for a genome to be considered
        percent_markers = int(sys.argv[1])
        input_dir = sys.argv[2]
        references = sys.argv[3]

        if os.stat(
                input_dir + '/marker_cov_per_genome.txt').st_size == 0:  # check if input_dir/marker_cov_per_genome.txt is empty
            open(input_dir + '/ref_genome_marker_gene_coverage.tsv', 'a').close()  # create empty output file
            open(input_dir + '/reference_candidates.txt', 'a').close()  # create empty output file
            return 1

        # Read in table with coverages of each marker gene
        covTable = pd.read_table(input_dir + '/marker_cov_per_genome.txt', header=0,
                                 names=['genome', 'gene', 'coverage'])

        # Convert to a table where rows are genomes and columns are marker genes
        # If a marker gene is covered in a particular genome, the table value is 1
        covTableOrg = covTable.pivot_table(index='genome', columns='gene', values='coverage')

        # Add a column containing how many genes are covered per genome in a sample
        tmp = pd.DataFrame(covTableOrg.count(axis=1), columns=['num_covered_genes'])
        covTableOrg2 = pd.merge(tmp, covTableOrg, left_index=True, right_index=True)

        # Add taxonomic information about each genome (how many marker genes it contains in total and what the assembly level is)
        genomes = pd.read_table(references + '/taxinfo.tsv')
        genomes = genomes.rename(columns={'#accession': 'genome'})

        covTableResults = pd.merge(genomes[['genome', 'num_genes', 'assembly_level']], covTableOrg2, left_on='genome',
                                   right_index=True)

        # Calculate the percentage of marker genes in the genome that are covered in the sample
        covTableResults['percentage_covered_genes'] = covTableResults['num_covered_genes'] / covTableResults[
            'num_genes'] * 100

        # Output the table where rows are genomes and columns are marker genes
        # 1= marker is covered in sample, 0=marker is not
        covTableResults = covTableResults.fillna(0)
        covTableResults.to_csv(input_dir + '/ref_genome_marker_gene_coverage.tsv', sep='\t', index=False, header=True)

        # Output a list of genomes with specified percentage of genes covered
        covTableResultsSub = covTableResults[covTableResults.percentage_covered_genes >= percent_markers]
        covTableResultsSub.genome.to_csv(input_dir + '/reference_candidates.txt', index=False, header=False)
    except Exception as e:
        print("Failed due to error", str(e))
        print('percent_markers: ' + percent_markers)
        print('input_dir: ' + input_dir)
        print('references: ' + references)
        sys.exit(1)


if __name__ == '__main__':
    main()
