#!/usr/bin/env python

## This script takes two inputs
## The first is a .csv file containing one line for each cluster
## and the comma separated records represent paths to the genome 
## sequences in the cluster
##
## the second parameter is the location where the genomes should be combined

# import sys
# import csv
# from typing import List
# from pathlib import Path
# from shutil import copyfileobj

# cluster_in_file = sys.argv[1]
# concat_out_file_path = sys.argv[2]

# with open(cluster_in_file, "r") as cluster_file:
#     cluster_list: List[List] = list(csv.reader(cluster_file))

#     all_files = []
#     for cluster_refs in cluster_list:
#         all_files.extend(cluster_refs)

#     with open(concat_out_file_path, "wb") as outf:
#         for i in range(0, len(all_files)):
#             with open(Path(all_files[i]).resolve().as_posix(), 'rb') as inf:
#                 copyfileobj(inf, outf)
#         outf.close()


# MOUMI: This version of the script adds more informative headers to the combined output file. 
# Each header will include the cluster index, genome ID (if it can be extracted from the file path), and contig ID. 
#!/usr/bin/env python

import sys
import csv
import re
from pathlib import Path
from typing import List


def extract_genome_id(file_path: Path) -> str:
    """
    Try to recover genome accession like GCA_000247755.2 or GCF_xxx
    from the path.
    """
    acc_pattern = re.compile(r'^(GC[AF]_\d+\.\d+)$')

    # Search parent directory names first
    for parent in file_path.parents:
        if acc_pattern.match(parent.name):
            return parent.name

    # Fallback: try file name
    m = re.search(r'(GC[AF]_\d+\.\d+)', file_path.name)
    if m:
        return m.group(1)

    return "UNKNOWN_GENOME"


cluster_in_file = sys.argv[1]
concat_out_file_path = sys.argv[2]

with open(cluster_in_file, "r") as cluster_file:
    cluster_list: List[List[str]] = list(csv.reader(cluster_file))

with open(concat_out_file_path, "w") as outf:
    for cluster_idx, cluster_refs in enumerate(cluster_list, start=1):
        for ref_path_str in cluster_refs:
            ref_path = Path(ref_path_str).resolve()
            genome_id = extract_genome_id(ref_path)

            with open(ref_path, "r") as inf:
                for line in inf:
                    if line.startswith(">"):
                        original_header = line[1:].strip()
                        contig_id = original_header.split()[0]
                        new_header = f">cluster_{cluster_idx}|{genome_id}|{contig_id}\n"
                        outf.write(new_header)
                    else:
                        outf.write(line)
