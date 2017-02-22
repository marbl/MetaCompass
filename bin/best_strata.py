#!/usr/bin/env python3

import sys

input_name=sys.argv[1]
output_name=sys.argv[2]

with open (input_name, 'rt') as infile,\
   open (output_name,'wt') as outfile:

   #skip header
   line = infile.readline()
   while line[0]=="@":
      last_pos = infile.tell()
      line = infile.readline()

   infile.seek(last_pos)

   #extract best-hits
   for line in infile:

      items = line.split('\t')

      if "XS" in items[12]: #multiple mapping --> test best-strata
         XS=int(items[12].split(':')[2],10)
         AS=int(items[11].split(':')[2],10)
         if AS>=XS:
            outfile.write(line)
      else: #unique mapping --> save
         outfile.write(line)
      
#print('end')


































