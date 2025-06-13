#!/usr/bin/env python

################################################################
###
### This script should be run in the reference assembly process
### It takes as input the mapped_genomes.txt file, and retrieves
### the list of genomes that were  used by the assembly, and then
### transfers the output of pilon to an output folder.
###
################################################################

import shutil
import os

inFName = "mapped_genomes.txt"

try: 
    inf = open(inFName, "r")
except OSError as err:
    print("Cannot open input: ", inFName, " ", err)

for name in inf: 
    name = name.strip()
    baseNm = os.path.splitext(name.strip())[0]
    outNm = baseNm + ".refctgs.fna"
    try:
        shutil.copy(name + "/pilon/contigs.pilon.fasta", outNm)
    except OSError as err:
        print("Cannot copy file ", name + "/pilon/contigs.pilon.fasta", " into ", outNm)
