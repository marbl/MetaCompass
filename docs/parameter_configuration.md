# Metacompass Parameter Configuration


## Required Input
- **forward**: Path to the file with forward reads, providing sequence information from one end of the DNA fragment in a specific direction.
- **reverse**: Path to the file with reverse reads, representing sequence information from the opposite end of the DNA fragment in the reverse direction.
- **output**: The output directory or file path where the results will be stored.
- **threads**: The number of threads to use during the assembly.

## Reference Selection
- **reference_db**: The path to the reference database used for reference selection.
- **filter_refs**: A boolean flag indicating whether to filter the reference database based on specified criteria.
- **ms**: The minimum score threshold for filtering references during reference selection.
- **clean**: K-mers clean threshold used during reference selection.
- **match**: K-mers match threshold used during reference selection.
- **masking**: The masking option used during reference selection.

## Cluster Reference Selection
- **depth_of_coverage**: Minimum depth of coverage required for a reference during cluster reference selection.
- **breadth_of_coverage**: Minimum breadth of coverage required for a reference during cluster reference selection.
- **percent_markers_covered**: Minimum percentage of markers that need to be covered by reads in a reference sequence for it to be considered for further analysis.

## Reference Culling


- **mcl(ref_culling.nf,align_reads)**:the minimal max contig length of the assembly assembled based on the current reference, for the algorithm to continue assembly based on the next reference in the current cluster. For example, when the mcl is set to 2,000 base pairs (bp), this means that for the assembly process to continue with the next reference in the current cluster, the assembly based on the current reference must include at least one contig exceeding 2000 bp. This condition ensures that the assembly process proceeds only if the length of the contig surpasses the specified minimal length. The default is 10,000
- **boc(ref_culling.nf,align_reads)**:the minimal breadth of coverage of the current reference, for the algorithm to continue assembly based on the next reference in the current cluster. The default is 5, meaning 5%.
- **Hierarchical_clustering_dissimarity(ref_culling.nf)**:In reference culling, hierarchical clustering was conducted using the Average Nucleotide Identity (ANI) of the references. This clustering was executed via the hierarchy.fcluster method, employing an average linkage with a similarity cutoff set at 95%. In this case the Hierarchical_clustering_dissimarity is set to 5 which corresponds to 5% dissimilarity. It is advised not to change this parameter.



## Reference-Guided Assembly (Pilon)
- **mindepth**: Minimum depth threshold used during reference-guided assembly with Pilon.
- **ref_sel**: The reference selection method used during reference-guided assembly with Pilon.
- **ref_pick**: The reference picking method used during reference-guided assembly with Pilon.
- **mincov**: Minimum coverage threshold used during reference-guided assembly with Pilon.
- **minctglen**: Minimum contig length threshold used during reference-guided assembly with Pilon.
- **de_novo**: The number of de novo iterations to perform during reference-guided assembly with Pilon.
- **tracks**: A boolean flag indicating whether to generate assembly tracks during reference-guided assembly with Pilon.

- **memory**: The memory allocation (in GB) for Pilon and reference-guided assembly.

## Other Parameters
- **readlen**: Expected read length of the sequencing data (default: 200).

## Skip steps
- **de_novo** (0/1, default=1): Skip de novo assembly when set to 0. 
- **skip_rs** (true/false, default=false): Skip reference selection.
- **skip_rc** (true/false, default=false): Skip reference culling.

## Nextflow Parameters
- **trace_file_name** : Contains information about each process executed in your pipeline script, including: submission time, start time, completion time, cpu and memory used. 
- **with-timeline** : Generates an html file where you can view how each nextflow process was executed


Note: Some parameters are dependent on the condition of other parameters (e.g., not used if certain flags are set to false), so the grouping may overlap in some cases.
