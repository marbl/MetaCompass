#!/usr/bin/env python

import sys

input_name=sys.argv[1]
output_name=sys.argv[2]

with open (input_name, 'rt') as infile,\
   open (output_name,'wt') as outfile:
   
   line = infile.readline()
   items = line.split('\t')
   oldquery = items[0]
   oldscore = items[-1].strip()
   infile.seek(0)

   for line in infile:
      items=line.split('\t')
      query=items[0]
      score=items[-1].strip()
      extract=False
      update=True

      if oldquery == query:  
         if oldscore == score: 
            extract=True
         else:
            update=False
      else:
         extract=True
   
      if extract and line:
         outfile.write(line)

      #keep the best score for same read
      if update:
         oldquery=query
         oldscore=score
