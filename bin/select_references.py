#!/usr/bin/env python

#############################################
#
# Program: Pick reference genomes based on
#          marker genes 
#
# Author: Victoria Cepeda
#
# Fri Apr 05 18:00:00 EDT 2019
#
# Description:uses "/refseq/markers/markers.length",
# "/refseq/markers/markers.numcogs", and blastn  
# output to calculate depth of coverage and 
# output the fasta file of candidate references
#############################################

import sys
import os
import operator


pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])


#========================================
#initalize marker coverage dictionaries
#========================================
# input:
#marker2alignlen: markers found in blastn output
# outputs:
#marker2length: dic of marker gene and [sumation of length per marker]
#marker2numcog: dic of genomes and [num cogs per genome]

def process_marker_db(marker2alignlen):
   marker2length={}
   marker2numcog={}
   mlength=0
   with open(pathbin + "/refseq/markers/markers.length", "rt") as markerfile:
   
      for line in markerfile:
         items = line.split("\t")
         marker=items[0]
         mlength=int(items[1])
         #ref = "_".join(marker.split('_')[0:2])
         #print ("%s" % ref)
         if marker in marker2alignlen:
            if marker in marker2length:
               marker2length[marker] += mlength
            else:
               marker2length[marker] = mlength
   with open(pathbin + "/refseq/markers/markers.numcogs", "rt") as markerfile:
      for line in markerfile:
          items = line.split("\t")
          marker=items[0]
          numcogs=int(items[1])
          marker2numcog[marker] = numcogs

   return marker2length,marker2numcog

#===========================================================
#extract blast information from the read to marker alignment
#===========================================================
#add alignment length per gene to genome
#sum number of marker genes per genome
# outputs:
#marker2alignlen is dictionary of marker and sumation of reads alglength per marker
#markertotalcount is dictionary of marker and sumation of reads per marker
def extract_blast_hits(blast_filename):
   marker2alignlen={}
   markertotalcount={}
   with open(blast_filename, "rt") as blastfile:
      for line in blastfile:
         items = line.split("\t")
         marker = items[1]
         #ref = "_".join(marker.split('_')[0:2])
         alglength = int(items[3])
         if marker in marker2alignlen:
            marker2alignlen[marker]+=alglength
         else:
            marker2alignlen[marker]=alglength
         if marker in markertotalcount:
            markertotalcount[marker] += 1
         else:
            markertotalcount[marker] = 1
   print ("extract_blast_hits:marker2alignlen,markertotalcount %s  %s" % (marker2alignlen,markertotalcount))
   return marker2alignlen,markertotalcount

#============================================================
#generate reference genome list data
#marker2alignlen
#markertotalcount
#marker2lenght
#covthreshold from stdin
#   refid_filename = outdir+"/"+prefix+".refseq.ids"#sys.argv[1]
#   assembly_id_filename = outdir+"/"+prefix+".assembly.ids"#sys.argv[1]
#   markercoverage_filename = outdir+"/"+prefix+ ".markercoverage.txt"#sys.argv[1]
#return total_marker_cov,myset
#total_marker_cov has dic of marker and summation marker2alignlen[marker]/marker2lenght[marker])
#============================================================
def marker_cov(marker2alignlen, marker2lenght,covthreshold):
   marker_cov={}
   #calculate coverage per marker gene
   #only output markers with minimum selected coverage=covthreshold of 1
   for marker in marker2alignlen:#is ref not marker
      coverage= float(marker2alignlen[marker]/marker2lenght[marker])  
      print ("marker, cov: %s %f" % (marker,coverage))
      if coverage >= covthreshold:
         marker_cov[marker] = coverage  
   print ("marker_cov: %s" % marker_cov)
   return marker_cov
#refid_filename = outdir+"/"+prefix+".refseq.ids"#sys.argv[1]
#assembly_id_filename = outdir+"/"+prefix+".assembly.ids"#sys.argv[1]
#assembly_markercoverage_filename or assembly_markercoverage_filename = outdir+"/"+prefix+ ".assemblymarkercoverage.txt"#sys.argv[1]   
#cog_filename = outdir+"/"+prefix+ ".cogmarkercoverage.txt"#sys.argv[1]   
#marker2numcog

def create_aux_dics2(markertotalcount, covthreshold,marker_cov):   
   #from highest to lowest coverage print marker and totalcoverage per marker
   #in addition, add refid(NC_) to ids list; add assemblyid(GCF) to ass list 
   cumulative_cov={}#assembly_acc cum cov
   cog_num={}#dict of list of cog acc per assembly_acc
   acc_ref={}#dict of list of acc ref per assembly_acc
   
   for marker, cov in sorted(marker_cov.items(), key=operator.itemgetter(1), reverse=True):
       ref_acc = "_".join(marker.split('_')[0:2])
       assembly_acc= "_".join(marker.split('_')[2:4])
       #print ("%s" % assembly_acc)
       cog_acc ="_".join(marker.split('_')[4:5])       
       if assembly_acc in cumulative_cov:
            cumulative_cov[assembly_acc]+=cov  
       else:
            cumulative_cov[assembly_acc]=cov
       if assembly_acc in cog_num:
            if cog_acc not in cog_num[assembly_acc]:
               cog_num[assembly_acc].append(cog_acc)
       else:
            cog_num[assembly_acc]=[cog_acc]
       if assembly_acc in acc_ref:
           if ref_acc not in acc_ref[assembly_acc]:
               acc_ref[assembly_acc].append(ref_acc)
       else:
            acc_ref[assembly_acc]=[ref_acc]
   print ("%s" % assembly_acc)
   return cumulative_cov,cog_num,acc_ref


def cog_stats2(markertotalcount,refid_filename, assembly_id_filename, assembly_markercoverage_filename, cog_filename,marker2numcog,cumulative_cov,cog_list,acc_ref):
   #cumulative_cov={}#assembly_acc cum cov
   #cog_list={}#assembly_acc cog numbers
   acc_list=[]#assembly_acc list
   #print marker and coverae
   outfile = open(cog_filename,'w')
   for acc, cog in cog_list.items():
       print ( "acc,cog, len(cog), cumulative_cov[acc]: %s\t%s\t%d\t%d" % (acc, cog,len(cog),cumulative_cov[acc]), file= outfile)       
#       if len(cog) >= 10 and cumulative_cov[acc]>=50:#at least 10 cogs with deep of coverage 1
       
       if len(cog) >= 1 and cumulative_cov[acc]>=len(cog):#at least 10 cogs with deep of coverage 1
           acc_list.append(acc)
           print ( "acc, cog,len(cog)):%s\t%s\t%d" % (acc, cog,len(cog)))       
   outfile.close()
   print (cog_list)
   print (acc_list)
   

   outfile = open(assembly_markercoverage_filename,'w')         
   for acc, cov in sorted(cumulative_cov.items(), key=operator.itemgetter(1), reverse=True): 
       if acc in acc_list:
           print ( "%s\t%.5f" % (acc, cov), file= outfile)
   outfile.close()
#   refid_filename = outdir+"/"+prefix+".refseq.ids"#sys.argv[1]
   outfile = open(refid_filename,'w')
   for acc, ref in acc_ref.items():
       if acc in acc_list:
           print ( "%s\t%s" % (acc, ref), file =outfile)
   outfile.close()
   #   assembly_id_filename = outdir+"/"+prefix+".assembly.ids"#sys.argv[1]
   outfile = open(assembly_id_filename,'w')
   #   ###arreglar
   for acc in acc_list:
       print ("%s" % (acc), file =outfile)
   outfile.close()

#input:ref2hash dictionary a
def calculate_mash(ref2hash,cutoff):    
    refids=[]
    if len(ref2hash)==1:
        for ref, hashn in ref2hash.items():
            return [ref]
    max_value=max(ref2hash.values())
    min_value=min(ref2hash.values())
    avg_value=sum(ref2hash.values())/len(ref2hash)

    #if min_value*1.05>=max_value:#all values are conserved or similar
    if (min_value*cutoff)*1.0>=max_value:#all values are conserved or similar
        for refid, hashes in sorted(ref2hash.items(), key=operator.itemgetter(1),reverse=True):
            if hashes == max_value:
                return [refid]
    refids.append(max(ref2hash.items(), key=operator.itemgetter(1))[0])        
    for refid, hashes in sorted(ref2hash.items(), key=operator.itemgetter(1),reverse=True):
        #if refid in refids:#skip first
        #    break
        #if hashes <= max_value: is sorted, it has to be <=
        if hashes*cutoff*1.0<max_value:# first will be spikked
            max_value=hashes
            refids.append(refid)
    return refids        
def highest_mash(ref2hash):    
    refids=[]
    if len(ref2hash)==1:
        print ("GFCUNIQ len(ref2hash)==1" )
        print ("ref2hash.keys:%s" % ref2hash.keys())
        print ("ref2hash.values:%s" % ref2hash.values())
        for ref, hashn in ref2hash.items():
            return [ref]
    max_value=max(ref2hash.values())
    min_value=min(ref2hash.values())
    avg_value=sum(ref2hash.values())/len(ref2hash)

    #return first highest mash
    for refid, hashes in sorted(ref2hash.items(), key=operator.itemgetter(1),reverse=True):
            if hashes == max_value:
                return [refid]

#select assembly tax <x%[x<=10] of variation in # of shared-hashes. Ties are resolved randomly.
#mash screen fields are [identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]:
def parse_mash_screen(screen_bacgeno,tax_dic,tax_dic_old,filt_assembly_id_filename,filt_assembly_fna_filename,filt_tax_file,cutoff):
   ref2tax={} #filtered results, 1 per tax
   ref2hash={}#mash screen results
   with open(screen_bacgeno, "rt") as mashfile:
        for line in mashfile:
         items = line.split("\t")
         identity=items[0]
         sharedhashes=items[1]
         hashes=int(sharedhashes.split("/")[0])
         medianmultiplicity=int(items[2])
         pvalue=float(items[3])
         queryid=items[4].split(".fasta")[0]
         ref2hash[queryid]=hashes
   temp2hash={}
   for tax,refs in tax_dic_old.items():
       for refid in refs:#refs is list of assembly IDS assigned to same tax
           if refid in ref2hash:#mash screen can have zero matches and need to be chacked here
               if refid not in temp2hash:
                   temp2hash[refid]=ref2hash[refid]
               else:
                   temp2hash[refid].append(ref2hash[refid])
       if temp2hash:
           ref2tax[tax]=highest_mash(temp2hash)
           #ref2tax[tax]=calculate_mash(ref2hash,cutoff)#option to keep some hashes per tax    
           temp2hash.clear()       

   #write filt_tax_file,filt_assembly_id_filename
   outfile = open(filt_tax_file,'w')
   for tax, ref in ref2tax.items():
       print ( "%s\t%s" % (tax, ref), file =outfile)
   outfile.close()

   outfile = open(filt_assembly_id_filename,'w')
   for tax, ref in ref2tax.items():
       for i in ref:
           print ("%s" % (i), file =outfile)
   outfile.close()
   

def main():
#----------------------------------------#
# read command line options
#----------------------------------------#
   query = ""
   blast = "blastn"
   outdir = ""
   prefix = "mc"
   nump = 0
   query  = sys.argv[1]
   reads  = sys.argv[2]
   outdir = sys.argv[3]
   nump   = sys.argv[4]
   covthreshold = float(sys.argv[5])
#----------------------------------------#
   ref = pathbin +"/refseq/markers/markers.refseq.dna"
   
   param = "-word_size 28"
# run blast
   #cmd = "%s %s -num_threads %s -evalue 1e-10 -perc_identity 97 -outfmt 6 -max_target_seqs 100 -query %s -db %s > %s/%s.%s.all"%(blast,param,nump,query,ref,outdir,prefix,blast)
   cmd = "%s %s -num_threads %s -evalue 1e-20 -perc_identity 97 -outfmt 6 -query %s -db %s > %s/%s.%s.all"%(blast,param,nump,query,ref,outdir,prefix,blast)

   print ("%s" % (cmd))
   os.system(cmd)
# filter blast best hit
   #cmd = "%s/bin/best-hits.py %s/%s.%s.all %s/%s.%s"%(pathbin,outdir,prefix,blast,outdir,prefix,blast)
   #print ("%s" % (cmd))
   #os.system(cmd)
    
#if (scalar @ARGV == 3) {
#{input} {output.out} {threads} {params.mincov} {params.readlen}

   blast_filename = outdir+"/"+prefix+"."+blast+".all"#sys.argv[1]
   #blast_filename = outdir+"/"+prefix+"."+blast#sys.argv[1]
   
   refid_filename = outdir+"/"+prefix+".refseq.assembly.ids"#sys.argv[1]
   assembly_id_filename = outdir+"/"+prefix+".assembly.ids"#sys.argv[1]
   #markercoverage_filename = outdir+"/"+prefix+ ".markercoverage.txt"#sys.argv[1]
   assembly_markercoverage_filename = outdir+"/"+prefix+ ".assemblymarkercoverage.txt"#sys.argv[1]   
   cog_filename = outdir+"/"+prefix+ ".cogmarkercoverage.txt"#sys.argv[1]   
   
   #coverage per marker gene or based on the all marker genes in a genome
   #extract blast info
   marker2alignlen, markertotalcount = extract_blast_hits(blast_filename)
   
   #extract gene info
   marker2length, marker2numcog = process_marker_db(marker2alignlen)

   marker_coverage={}

   #find total marker coverage =marker2alignlen[marker]/marker2lenght[marker])
   marker_coverage=marker_cov(marker2alignlen, marker2length,covthreshold)
   cumulative_cov={}#assembly_acc cum cov
   cog_num={}#dict of list of cog acc per assembly_acc
   acc_ref={}
   #find total marker coverage =marker2alignlen[marker]/marker2lenght[marker])
   cumulative_cov,cog_num,acc_ref=create_aux_dics2(markertotalcount, covthreshold,marker_coverage)
   #get statistics from cogs 
   cog_stats2(markertotalcount,refid_filename, assembly_id_filename, assembly_markercoverage_filename, cog_filename,marker2numcog,cumulative_cov,cog_num,acc_ref)

   print ("# Extract reference genome sequences")
   #search ids from same assembly project
   cmd ="for file in $(cat %s/%s.assembly.ids);do cat  %s/refseq/genomes/${file}/${file}.fasta >> %s/%s.refseq.fna;done"%(outdir,prefix,pathbin,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   cmd ="grep '>'  %s/%s.refseq.fna|tr -d '>'|cut -f1 -d ' '> %s/%s.refseq.ids" %(outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   refseq_tax = pathbin + "/refseq/refseq_tax"
   
   tax_file = outdir + "/"+ prefix + ".tax"
   cmd="grep -f %s/%s.assembly.ids %s |sort -k2,2 >>%s/%s.tax" %(outdir,prefix,refseq_tax,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   tax_file = outdir + "/"+ prefix + ".tax"
   print ("%s" % (tax_file))
#store taxes 
   tax_dic_old={}
   tax_dic={}
   prev_tax=""
   with open(tax_file, "rt") as taxfile:
       for line in taxfile:
           items = line.split("\t")
           ref_acc=items[0]
           tax=int(items[1])
           if tax in tax_dic_old:
               tax_dic_old[tax].append(ref_acc)
           else:
               tax_dic_old[tax] = [ref_acc]
   filt_tax_file=outdir + "/"+ prefix + ".assemblypertax.tax"
   outfile = open(filt_tax_file,'w')
###############to speed up for now
   for tax, ref in tax_dic_old.items():
       print ( "%s\t%s" % (tax, ref), file =outfile)
   outfile.close()                     
           
   with open(tax_file, "rt") as taxfile:
       for line in taxfile:
           items = line.split("\t")
           ref_acc=items[0]
           tax=int(items[1])
           if ref_acc in tax_dic:
               tax_dic[ref_acc].append(tax)
           else:
               tax_dic[ref_acc] = [tax]

   cmd ="for file in $(cat %s/%s.assembly.ids);do echo  %s/refseq/genomes/${file}/${file}.fasta.msh >> %s/%s.msh.list;done"%(outdir,prefix,pathbin,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)  
   
   cmd="mash paste -l %s/%s.assembly %s/%s.msh.list"%(outdir,prefix,outdir,prefix)
   #cmd="%s/bin/scratch/mash-Linux64-v2.1/mash paste -l %s/%s.assembly %s/%s.msh.list"%(pathbin,outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd) 
   #run mash screen
   screen=outdir+"/mc.mashscreen.tab"
   cmd="mash screen -p %s %s/%s.assembly.msh  %s >%s"%(nump,outdir,prefix,reads,screen)
   #cmd="%s/bin/scratch/mash-Linux64-v2.1/mash screen -p %s %s/%s.assembly.msh  %s >%s"%(pathbin,nump,outdir,prefix,reads,screen)
   print ("%s" % (cmd))
   os.system(cmd)
   #parse mash screen, sort from highest to lowest # shared-hashes
   filt_assembly_id_filename = outdir+"/"+prefix+".filt.assembly.ids"
   filt_assembly_fna_filename = outdir+"/"+prefix+".refseq.filt.fna"
   
   filt_tax_file=outdir + "/"+ prefix + ".filt.tax"
   cutoff=1.1
   parse_mash_screen(screen,tax_dic,tax_dic_old,filt_assembly_id_filename,filt_assembly_fna_filename,filt_tax_file,cutoff)#cut -f1
   cmd ="for file in $(cat %s);do cat %s/refseq/genomes/${file}/${file}.fasta >> %s;done"%(filt_assembly_id_filename,pathbin,filt_assembly_fna_filename)
   print ("%s" % (cmd))
   os.system(cmd)
   ###new addition to map acc to assemblies
   cmd ="for file in $(cat %s);do acc=$(cut -f1 %s/refseq/genomes/${file}/${file}.fasta.fai); echo $file'\t'$acc >> %s;done"%(filt_assembly_id_filename,pathbin,filt_assembly_id_filename+".map")
   print ("%s" % (cmd))
   os.system(cmd)
   
   cmd ="for file in $(cat %s);do acc=$(cut -f1 %s/refseq/genomes/${file}/${file}.fasta.fai); for j in $acc;do echo $file'\t'$j >> %s;done;done"%(filt_assembly_id_filename,pathbin,filt_assembly_id_filename+".map2acc")
   print ("%s" % (cmd))
   os.system(cmd)
   
   
   cmd ="grep '>'  %s/%s.refseq.filt.fna|tr -d '>'|cut -f1 -d ' '> %s/%s.refseq.filt.ids" %(outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)

if __name__ == '__main__':
    main() 
