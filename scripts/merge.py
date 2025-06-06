#!/usr/bin/env python3

import sys
from pathlib import Path

headers = set()

def append_to_file(file_to_read, file_to_append):
    with open(file_to_append, "a") as filetoappend:
        with open(file_to_read, "r") as filetoread:
            lines = filetoread.readlines()
            for line in lines:
                if line[0] != '>': #sequence line
                    filetoappend.write(line)
                else: #Header line
                    header = line.rstrip()[1:] #header without trailing whitespace
                    count = 1
                    temp_header = header

                    #append _1, _2, or _3 ... if header was seen before
                    while temp_header in headers:
                        temp_header = header
                        temp_header += ("_" + str(count)) 
                        count += 1
    
                    headers.add(temp_header)
                    new_header = ">" + temp_header + "\n"
                    filetoappend.write(new_header)
        filetoappend.write('\n')

# Concat config files. Make sure headers are unique
def main():
    denovo_file = sys.argv[1] #denovo assembly config file
    ref_folder = sys.argv[2] #folder contains reference assembly config files
    output_file = sys.argv[3]

    #Create/overwrite output file
    with open(output_file, "w") as _:
        
        #Read denovo config file and append to output file
        append_to_file(denovo_file, output_file)
        
        #Read ref_assembly config files and append to output file
        files = sorted(Path(ref_folder).rglob("contigs.pilon.fasta"))
        for file in files:
            append_to_file(file, output_file)

if __name__ == '__main__':
	exit_code = main()
	if exit_code != 0:
		sys.exit(exit_code)
