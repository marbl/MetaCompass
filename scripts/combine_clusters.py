#!/usr/bin/env python

## This script takes two inputs
## The first is a .csv file containing one line for each cluster
## and the comma separated records represent paths to the genome 
## sequences in the cluster
##
## the second parameter is the location where the genomes should be combined

import sys
import csv
from typing import List
from pathlib import Path
from shutil import copyfileobj

cluster_in_file = sys.argv[1]
concat_out_file_path = sys.argv[2]

with open(cluster_in_file, "r") as cluster_file:
    cluster_list: List[List] = list(csv.reader(cluster_file))

    all_files = []
    for cluster_refs in cluster_list:
        all_files.extend(cluster_refs)

    with open(concat_out_file_path, "wb") as outf:
        for i in range(0, len(all_files)):
            with open(Path(all_files[i]).resolve().as_posix(), 'rb') as inf:
                copyfileobj(inf, outf)
        outf.close()

