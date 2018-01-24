#!/usr/bin/env python

#############################################
#
# Program: Pick reference genomes based on
#          metaphyler output
#
# Author: Mathieu Almeida, Bo Liu, Victoria Cepeda
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
 
#============================
#extract metaphyler taxa data
#============================
def extract_reference_data():
   seq2tid={}
   tid2seqs={}
   tid2name={}
   tid2sp={}
   tid2seq={}
   pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])   

   with open(pathbin + "/refseq/tid2par.tab", "rt") as reffile:
     for line in reffile:
        
        seqid, taxonid, species, genus, name = line.split("\t")
        
        seq2tid[seqid] = taxonid.strip()
        tid2name[taxonid] = name.strip()
        tid2sp[taxonid] = species.strip()
        if not tid2seq.get(taxonid):
           tid2seq[taxonid] = []
        tid2seq[taxonid].append(seqid)
                   
   return seq2tid, tid2name, tid2sp, tid2seq


#========================================
#initalize marker coverage arrays
#========================================
def init_marker_data(marker2coverage,markertotalcount):
   #add blast list marker2coverage
   pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])
   #marker2coverage={}
   #markertotalcount={}
   with open(pathbin + "/src/metaphyler/markers/markers.refseq.dna", "rt") as markerfile:
      header=""
      seq=""
      flag=0
      for line in markerfile:
         if line[0]==">":
           #first line
           if header!="" and flag == 0:
                  ref = "_".join(header.split('_')[0:2]
                  print "%s" % [0]
                  marker2coverage[header] = [0] * (len(seq)-1)
              #markertotalcount[header] = 0
           header= line.strip()[1:]
           if marker2coverage.has_key(header):
               flag=0
           else:
               flag=1
               header=""
           seq=""
         else:
            if flag ==0:
                seq+= line.strip()
      #last line
      if flag ==0:
         ref = "_".join(header.split('_')[0:2])
         marker2coverage[header] = [0] * (len(seq)-1)            
      #markertotalcount[header]=0
   return marker2coverage


#===========================================================
#extract blast information from the read to marker alignment
#===========================================================
def extract_marker_data(blast_filename, marker2coverage, markertotalcount ):
  
   with open(blast_filename, "rt") as blastfile:
      for line in blastfile:
         items = line.split("\t")
         marker = items[1]
         mstart = int(items[8])
         mend = int(items[9])
         if int(mstart) < int(mend):
            pos=int(mstart)
            while pos < int(mend):
               if !marker2coverage.has_key(marker):
             
               marker2coverage[marker][pos-1] += 1
               #markertotalcount[marker]=0
               pos += 1
         else:
            pos=int(mend)
            while pos < int(mstart):
               marker2coverage[marker][pos-1] += 1
               #markertotalcount[marker]=0
               pos += 1
         markertotalcount[marker] += 1
         #markertotalcount[marker]=0
   #return marker2coverage, markertotalcount
   

#============================================================
#generate reference genome list data for metaphyler
#============================================================
def create_refid_file(marker2coverage, markertotalcount, seq2tid, tid2name, tid2sp, tid2seq, covthreshold, nbsamespeciesthreshold):
   total_ref_cov={}   
   #list_species=[]
   list_ref_median_cov={}

   for marker in marker2coverage:
      ref = "_".join(marker.split('_')[:2])
      if not total_ref_cov.get(ref):
         total_ref_cov[ref] = marker2coverage[marker]
      else:
         total_ref_cov[ref] = total_ref_cov.get(ref) + marker2coverage[marker]

   #calculate median
   for ref in total_ref_cov:
      list_ref_median_cov[ref] = float(median(total_ref_cov[ref]))

   #write result
   for ref, median_cov in sorted(list_ref_median_cov.items(), key=operator.itemgetter(1), reverse=True):
      taxa = seq2tid[ref]
      name = tid2name[taxa]
      #species = tid2sp[taxa]
      
      if median_cov >= 0.0:#covthreshold:
         #species_count = list_species.count(species) 
         #if species_count < nbsamespeciesthreshold:
            print("%s\t%s\t%.1f\t%s" % (ref,taxa, median_cov,name))
            #print("%s\t%.1f\t%s" % (ref,median_cov,name))
            #print ("%s\t%.1f" % (ref,median_cov))
            #for refid in tid2seq[taxa]:
            #   print("%s\t%.1f\t%s" % (refid,median_cov,name))
               #print("%s\t%s\t%.1f\t%s" % (refid,taxa, median_cov,name)) 
            #list_species.append(species)

#===========================
#MAIN
#===========================
def main():
   """Main program function
   """
   blast_filename = sys.argv[1]
   covthreshold = float(sys.argv[2])
   nbsamespeciesthreshold=5
   marker2coverage={}
   markertotalcount={}
   seq2tid, tid2name, tid2sp,tid2seq = extract_reference_data()
   extract_marker_data(blast_filename, marker2coverage, markertotalcount )
   init_marker_data(marker2coverage)
   create_refid_file(marker2coverage, markertotalcount, seq2tid, tid2name, tid2sp, tid2seq, covthreshold, nbsamespeciesthreshold)
          

if __name__ == '__main__':
    main()
 
