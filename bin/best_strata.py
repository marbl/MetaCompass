#!/usr/bin/env python

import sys

input_name=sys.argv[1]
output_name=sys.argv[2]

with open (input_name, 'rt') as infile,\
   open (output_name,'wt') as outfile:

   line = infile.readline()

   while line[0]=="@":
      last_pos = infile.tell()
      line = infile.readline()

   items = line.split('\t')
   oldquery = items[0]

   for score in items[12:]:
      if score[:2]=="NM":
         oldNM = score.split(':')[2]
         break
   infile.seek(last_pos)

   #print('start parsing...')
   for line in infile:
      items=line.split('\t')
      query=items[0]

      for score in items[12:]:
         if score[:2]=="NM":
            NM = score.split(':')[2]
            break

      extract=False
      update=True

      if oldquery == query:  
         if oldNM == NM: 
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
         oldNM=NM

   #print('end')
