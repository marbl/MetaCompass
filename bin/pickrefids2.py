#!/usr/bin/env python

#############################################
#
# Program: Pick reference genomes based on
#          metaphyler output
#
# Author: Victoria Cepeda, Mathieu Almeida, Bo Liu,
#
# Thu Sep 28 07:55:00 EDT 2017
#
#############################################

import sys
import os
import operator

#============================
#calculate median coverage
#============================
def median(L):
   if sum(L)>0:
         
         srt = sorted(L)
         n = len(L)
         mid = n//2
         
         if n % 2:
            return float(srt[mid])
         else:
            med = (srt[mid] + srt[mid-1]) / 2.0 
            return float(med)
   
   else:
      return 0.0

#============================
#calculte average coverage
#============================
def average(x):
   sum_x = sum(l)
   len_x = len(l)
   return( sum_x / float(len_x))
 
#good
#mejor filtrar antes este archivo
#da marker[nc...]=largo, en otras palabras
#========================================
#initalize marker coverage arrays
#========================================
def init_marker_data():
   pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])
   marker2length={}
   mlength=0
   #with open(pathbin + "/src/metaphyler/markers/markers.refseq.dna", "rt") as markerfile:
   with open(pathbin + "/src/metaphyler/markers/markers.length", "rt") as markerfile:
      id=""
      for line in markerfile:
         items = line.split("\t")
         #print ("%s" % items)
         id=items[0]
         #print ("%s" % id)
         mlength=int(items[1])
         #print ("%d" % mlength)
         ref = "_".join(id.split('_')[0:2])
         #print ("%s" % ref)
##marker2length[ref] = 0
         if ref in marker2length:
            marker2length[ref] += mlength
         else:
            marker2length[ref] = mlength
#print ("%s" % marker2length)
   return marker2length

#?
#===========================================================
#extract blast information from the read to marker alignment
#===========================================================
def extract_marker_data(blast_filename):
   marker2coverage={}
   markertotalcount={}
   with open(blast_filename, "rt") as blastfile:
      for line in blastfile:
         items = line.split("\t")
         marker = items[1]
         ref = "_".join(marker.split('_')[0:2])
         alglength = int(items[3])
         if ref in marker2coverage:
            marker2coverage[ref]+=alglength
         else:
            marker2coverage[ref]=alglength
         if ref in markertotalcount:
            markertotalcount[ref] += 1
         else:
            markertotalcount[ref] = 1
#print ("%s  %s" % (marker2coverage,markertotalcount))
   return marker2coverage,markertotalcount

#nope...
#============================================================
#generate reference genome list data
#marker2coverage
#markertotalcount
#covthreshold ???
#============================================================
def create_refid_file(marker2coverage, marker2lenght, markertotalcount,covthreshold):
   total_ref_cov={}   
   #list_species=[]
   #list_ref_median_cov={}
   for ref in marker2coverage:#is ref not marker
         total_ref_cov[ref] = float(marker2coverage[ref]/marker2lenght[ref])
         #total_ref_cov[ref] = float(markertotalcount*(marker2coverage[ref]/marker2lenght[ref]))

   #calculate median
   #for ref in total_ref_cov:
   #   list_ref_median_cov[ref] = float(median(total_ref_cov[ref]))
   #calculate avg
   #write result
   #for ref, median_cov in sorted(list_ref_median_cov.items(), key=operator.itemgetter(1), reverse=True):
   for ref, median_cov in sorted(total_ref_cov.items(), key=operator.itemgetter(1), reverse=True):
       #if median_cov >= covthreshold:
       print("%s\t%.5f" % (ref, median_cov))
#list_species.append(species)

#===========================
#MAIN
#===========================
def main():
   """Main program function
   """
   blast_filename = sys.argv[1]
   #coverage per marker gene or based on the all marker genes in a genome
   covthreshold = float(sys.argv[2])
   #what is marker2coverage and markertotalcount?
   #marker2length={}
   marker2length = init_marker_data()
   marker2coverage, markertotalcount = extract_marker_data(blast_filename)
   create_refid_file(marker2coverage, marker2length, markertotalcount, covthreshold)
          

if __name__ == '__main__':
    main()
 
