#!/usr/bin/env python

#############################################
#
# Program: Pick reference genomes based on
#          marker genes 
#
# Author: Victoria Cepeda-Espinoza
#
# Fri Apr 05 18:00:00 EDT 2019
#
# Description:uses "/refseq/markers/markers.length",
# "/refseq/markers/markers.numcogs", 
#"/refseq/markers/markers.clusters",
#   cdhitclusters_file = pathbin +"/refseq/markers/contigs_clusters"

# and blastn  
# output to calculate depth of coverage and 
# output the fasta file of candidate references
#############################################

import sys
import os
import operator
import time
import datetime
from datetime import timedelta

pathbin="/".join(os.path.realpath(__file__).split('/')[:-2])
        
def main():
#----------------------------------------#
# read command line options
#----------------------------------------#
   query = ""
   outdir = ""
   prefix = "mc"
   nump = 0
   selection=sys.argv[1]
   query  = sys.argv[2]
   reads  = sys.argv[3].replace(',', ' ')
   outdir = sys.argv[4]
   nump   = sys.argv[5]
   cogcov = float(sys.argv[6])#cog coverage
   print ("COGCOV:%f"% (cogcov))
   cutoff = float(sys.argv[7])#0.8identity mash or ANI
#if (scalar @ARGV <6) {
#{selection.type} {input} {output.out} {threads} {params.mincov} {params.readlen}

#----------------------------------------#
# define variables and files
#----------------------------------------#
   ref = pathbin +"/refseq/markers/markers.clusters"
   print ("%s" % (pathbin))
#----------------------------------------#
# run blast
#----------------------------------------#
   cmd = "/usr/bin/time blastn -word_size 28 -num_threads %s -evalue 1e-20 -perc_identity 97 -outfmt 6 -query %s -db %s > %s/%s.blastn.all"%(nump,query,ref,outdir,prefix)
   #modifie blastdb with representative marker genes
   #create cluster files with more than 1 gene.
   #map blast output to full list of marker genes
   print ("%s" % (cmd))
   os.system(cmd)
   
   print ("# Extract clusters")
   blast_filename = outdir+"/"+prefix+".blastn.all"
   blast_names = outdir+"/"+prefix+".blastn.all_markers"
   
   cdhitclusters_file = pathbin +"/refseq/markers/contigs_clusters"   
   #search ids from same assembly project
   cmd ="cut -f2 %s|sort|uniq >%s; %s/bin/extractSeq  %s %s >%s/contigs_clusters"%(blast_filename,blast_names,pathbin,cdhitclusters_file,blast_names,outdir)
   print ("%s" % (cmd))
   os.system(cmd)
   cmd="/usr/bin/time %s/bin/breadth %s %s %s >%s/breadth.log"%(pathbin,cogcov,outdir,pathbin,outdir)
   print ("%s" % (cmd))
   os.system(cmd)
   assembly_id_filename = outdir+"/"+prefix+".assembly.ids"
   assembly_id_filename2 = outdir+"/"+prefix+".assembly.cov_cogs"
   
   print ("# Extract reference genome sequences")
   
   cmd ="cut -f1 %s |sort|uniq> %s"%(assembly_id_filename2,assembly_id_filename)
   print ("%s" % (cmd))
   os.system(cmd)
   
   #search ids from same assembly project
   cmd ="for file in $(cat %s/%s.assembly.ids);do cat  %s/refseq/genomes/${file}/${file}.fasta >> %s/%s.refseq.fna;done"%(outdir,prefix,pathbin,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   cmd ="grep '>'  %s/%s.refseq.fna|tr -d '>'|cut -f1 -d ' '> %s/%s.refseq.ids" %(outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)

####mash option
   if selection =="all":
       sys.exit(0)#Success!!
   refseq_tax = pathbin + "/refseq/refseq_acc_tax"
   print ("%s" % refseq_tax)
   tax_file = outdir + "/"+ prefix + ".tax"
   #3rd bottleneck...
   cmd="grep -f %s/%s.assembly.ids %s |sort -k2,2 >>%s/%s.tax" %(outdir,prefix,refseq_tax,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   print ("%s" % (tax_file))
#store taxes 
   tax_dic={}
   prev_tax=""
   with open(tax_file, "rt") as taxfile:
       for line in taxfile:
           items = line.split("\t")
           ref=items[0]
           acc=items[1]
           tax=int(items[2])
           if tax in tax_dic:
               tax_dic[tax].append((ref,acc))
           else:
               tax_dic[tax] = [(ref,acc)]
   filt_tax_file=outdir + "/"+ prefix + ".assemblypertax.tax"
   ###############to speed up for now
   outfile = open(filt_tax_file,'w')
   for tax, ref in tax_dic.items():
       print ( "%s\t%s" % (tax, ref), file =outfile)
   outfile.close()
   #clus_dic={}
   #read_clusters(tax_dic,clus_dic)
   
   cmd ="for file in $(cat %s/%s.assembly.ids);do echo  %s/refseq/genomes/${file}/${file}.msh >> %s/%s.msh.list;done"%(outdir,prefix,pathbin,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)  
   
   cmd="mash paste -l %s/%s.assembly %s/%s.msh.list 2>>%s/%s.mash.log"%(outdir,prefix,outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd) 
   #run mash screen
   screen=outdir+"/mc.mashscreenw.tab"
   cmd="mash screen -w -i 0.8 -p %s %s/%s.assembly.msh  %s >%s 2>>%s/%s.mash.log"%(nump,outdir,prefix,reads,screen,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   allhashes=outdir+"/mc.mashscreen.tab"
   
   cmd="mash screen -i 0.9 -p %s %s/%s.assembly.msh  %s >%s 2>>%s/%s.mash.log"%(nump,outdir,prefix,reads,allhashes,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   #parse mash screen, sort from highest to lowest # shared-hashes

  #write filt_tax_file,filt_assembly_id_filename
  #need assemblyperacc
   cmd ="mv %s/%s.refseq.ids %s/%s.refseq.all.ids" %(outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   cmd ="mv %s/%s.refseq.fna %s/%s.refseq.all.fna" %(outdir,prefix,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)
   cutoff=0.9
   filt_assembly_id_filename = outdir+"/"+prefix+".filt.assembly.ids"
   filt_assembly_fna_filename = outdir+"/"+prefix+".refseq.fna"
   filt_tax_file=outdir + "/"+ prefix + ".filt.tax"
   refid_filename = outdir+"/"+prefix+".refseq.ids"
   refseq_ani = pathbin + "/refseq/refseq_ani"
   
   cmd="/usr/bin/time %s/bin/processmash %s %s %s %s >%s/processmash.log"%(pathbin,refseq_ani,allhashes,tax_file,refid_filename,outdir)
   print ("%s" % (cmd))
   os.system(cmd)
   cmd="%s/bin/extractSeq  %s/%s.refseq.all.fna %s >%s/%s.refseq.fna"%(pathbin,outdir,prefix,refid_filename,outdir,prefix)
   print ("%s" % (cmd))
   os.system(cmd)

#///////////////   mash -w
#   filt_assembly_id_filename = outdir+"/"+prefix+".filt.assembly2.ids"
#   filt_assembly_fna_filename = outdir+"/"+prefix+".refseq2.fna"
#   filt_tax_file=outdir + "/"+ prefix + ".filt2.tax"
#   refid_filename = outdir+"/"+prefix+".refseq2.ids"
#   cmd="/usr/bin/time %s/bin/processmash %s %s %s >%s/processmashw.log"%(pathbin,refseq_ani,screen,refid_filename,outdir)
#   print ("%s" % (cmd))
#   os.system(cmd)
#   cmd="%s/bin/extractSeq %s/%s.refseq.all.fna %s >%s/%s.refseq2.fna"%(pathbin,outdir,prefix,refid_filename,outdir,prefix)
#   print ("%s" % (cmd))
#   os.system(cmd)
#///////////////   

if __name__ == '__main__':
    main() 
