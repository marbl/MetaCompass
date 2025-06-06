import re
import numpy as np
import sys
from scipy.cluster import hierarchy
import scipy
from collections import defaultdict
from pathlib import Path

sys.setrecursionlimit(100000)

file = sys.argv[1]
path = sys.argv[3]

if Path(file).stat().st_size == 0:  # the file is empty then that means we have at max 1 reference candidates

    ref_culling_path = Path(path)
    matching_files = ref_culling_path / "matching_files.txt"
    clusters_txt = ref_culling_path / "clusters.csv"

    matched_genome = matching_files.read_text()
    clusters_txt.write_text(matched_genome)

    # create empty clusters_name.txt
    clusters_name = ref_culling_path / "clusters_name.txt"
    matched_genome_name = matched_genome.split("/")[-1]
    with open(clusters_name, "w") as f:
        modified_items = re.match("(GCA_[0-9]+\.[0-9]+)", matched_genome_name).group(1)
        f.write(str('\t'.join(modified_items) + '\n'))




else:
    if 'mash' in file:
        print("ANI matrix obtained from Mash detected.")
    if 'fastani' in file:
        print("ANI matrix obtained from FastANI detected.")

    counter = 0
    items = 0
    labels = []
    condensed = []
    matrix = []
    delim = '\t'
    for line in open(file, 'r'):
        if counter == 0:
            spl = line.split(delim)
            if len(spl) > 2:
                items = len(spl)
            else:
                items = int(line.split(delim)[-1])
            matrix = [[] for x in range(items)]
            counter += 1
            continue
        if delim in line:
            spl = line.split(delim)
        else:
            spl = line.split()
        labels.append(spl[0])
        endpoints = range(1, counter)
        print(counter)
        for i in endpoints:
            if 'mash' in file:
                matrix[i - 1].append(100 - 100 * float(spl[i]))
            elif 'fastani' in file:
                matrix[i - 1].append(float(spl[i]))
            else:
                print(i)
                print(spl)
                if float(spl[i]) <= 1:
                    matrix[i - 1].append(float(spl[i]) * 100)
                else:
                    matrix[i - 1].append(float(spl[i]))
        counter += 1

    for vec in matrix:
        for score in vec:
            condensed.append(100 - score)

    Z = hierarchy.linkage(condensed, 'average')
    # Z = hierarchy.linkage(condensed, method="ward")
    if len(sys.argv) > 2:
        vmax = float(sys.argv[2])

    clusters = hierarchy.fcluster(Z, vmax, criterion='distance')

    cluster_dict = defaultdict(list)
    for i, cluster_id in enumerate(clusters):
        cluster_dict[cluster_id].append(labels[i])

    for cluster_id, items in cluster_dict.items():
        print(f"Cluster {cluster_id}:")
        for item in items:
            print(f"\tItem: {item}")

    with open(path + "/clusters.csv", "w") as f:
        for cluster_id, items in cluster_dict.items():
            f.write(','.join(items) + '\n')

    with open(path + "/clusters_name.txt", "w") as f:
        for cluster_id, items in cluster_dict.items():
            items = [item.split('/')[-1] for item in items]
            modified_items = [re.match("(GCA_[0-9]+\.[0-9]+)", item).group(1) for item in items]
            print(f"Cluster {cluster_id}:")
            for item in modified_items:
                print(f"\tItem: {item}")
            f.write(str('\t'.join(modified_items) + '\n'))
