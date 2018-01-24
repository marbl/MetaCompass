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
#calculate median of list
#============================
def median2(L):
   median = 0
   sortedlist = sorted(L)
   lengthofthelist = len(sortedlist)
   centerofthelist = lengthofthelist / 2
   if len(sortedlist) % 2 == 0:
        temp = 0.0
        medianparties = []
        medianparties = sortedlist[centerofthelist -1 : centerofthelist +1 ]
        for value in medianparties:
            temp += value
            median = temp / 2
        return median
   else:
        return sortedlist[centerofthelist]
        
        
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
#calculte average of list
#============================
def average(L):
    return( sum(L)/ float(len(L)))
 
#good
#mejor filtrar antes este archivo
#da marker[nc...]=largo, en otras palabras
#========================================
#initalize marker coverage arrays
#========================================
def init_marker_data(marker2coverage):
   pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])
   marker2length={}
   mlength=0
   with open(pathbin + "/src/metaphyler/markers/markers.length", "rt") as markerfile:
      id=""
      for line in markerfile:
         items = line.split("\t")
         #print ("%s" % items)
         marker=items[0]
         #print ("%s" % id)
         mlength=int(items[1])
         #print ("%d" % mlength)
         #ref = "_".join(marker.split('_')[0:2])
         #print ("%s" % ref)
         if marker in marker2coverage:
            if marker in marker2length:
               marker2length[marker] += mlength
            else:
               marker2length[marker] = mlength
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
         #ref = "_".join(marker.split('_')[0:2])
         alglength = int(items[3])
         if marker in marker2coverage:
            marker2coverage[marker]+=alglength
         else:
            marker2coverage[marker]=alglength
         if marker in markertotalcount:
            markertotalcount[marker] += 1
         else:
            markertotalcount[marker] = 1
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
   total_marker_cov={}
   list_ref_median_cov={}
   list_ref_avg_cov={}
   #calculate coverage per gene
   for marker in marker2coverage:#is ref not marker
       total_marker_cov[marker] = float(marker2coverage[marker]/marker2lenght[marker])
   prev_ref=" "
   aux_list=[]
   #calculate median coverage per genome
   for marker in total_marker_cov:
       ref = "_".join(marker.split('_')[0:2])
       if ref == prev_ref:
          aux_list.append(total_marker_cov[marker])
       else:
          list_ref_median_cov[prev_ref] = median(aux_list)
          del aux_list
          aux_list = []
          aux_list.append(total_marker_cov[marker])
       prev_ref=ref
   #add last value
   list_ref_median_cov[ref] = average(aux_list)
   #for ref, median_cov in sorted(total_ref_cov.items(), key=operator.itemgetter(1), reverse=True):

   for ref, median_cov in sorted(list_ref_median_cov.items(), key=operator.itemgetter(1), reverse=True):
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
   marker2coverage, markertotalcount = extract_marker_data(blast_filename)
   marker2length = init_marker_data(marker2coverage)
   create_refid_file(marker2coverage, marker2length, markertotalcount, covthreshold)
          

if __name__ == '__main__':
    main()
 
