#!/usr/bin/env python

"""Check the metagenomic assembly by comparing with reference genome(s) """ 

import argparse
import os
from os import listdir
from os.path import isfile, join
import sys
import subprocess
import time

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def isdir(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isdir(path):
        if os.path.isfile(path):
            msg = "{0} is a file.".format(path)
            raise argparse.ArgumentTypeError(msg)
        else:
            msg = "{0} does not exist.".format(path)        
    
    return path


def config_parameters():
    """Extract program options
    """
    parser = argparse.ArgumentParser(description=__doc__,
                            usage="{0} -h [options] [arg]".format(sys.argv[0]))
    parser.add_argument('-c', '--contig_file', required=True, 
                        dest='contig_file', type=isfile, 
                        help='metagenomic contig file')                          
    parser.add_argument('-f', '--ref_folder', dest='ref_folder',
                        required=True, 
                        type=isdir, default='', help='Folder with reference .fasta files')  
    parser.add_argument('-r', dest='result_dir', type=isdir,
                        default=os.curdir + os.sep, help='Path to result '
                        'directory.')                   
    parser.add_argument('-o', dest='output_name', type=str,
                        default="assembly_summary.txt", help='Summary output table name')
    parser.add_argument('-a', dest='num_threads', type=int,
                        default=1, help='#CPU for blast and marker search')                 
               
    return parser.parse_args()    

#=========================================================================================
#extract genome list to analyse
#=========================================================================================
def extract_list_ref(ref_folder):

   list_ref = []
   list_ref = [f for f in listdir(ref_folder) if os.path.isfile(join(ref_folder, f))]

   return list_ref

#=========================================================================================
#generate nucmer report file
#=========================================================================================
def compare_contigs_to_ref(contig_file, ref_folder, list_ref, result_dir):

   os.system("touch %s/ready.done" % result_dir)
   for ref in list_ref:
      complete_ref = ref_folder + "/" + ".".join(ref.split('.')[:-1]) + ".fasta" 
      genome_name = ".".join(ref.split('.')[:-1])
      if not os.path.isfile("%s/%s.done" % (result_dir, genome_name)): 
         #nucmer and gene prediction comparison
         cmd = "nucmer -l 100 {0} {1}/{2}.mc.fasta -p {1}/{2} &>{1}/{2}.log; "\
            "delta-filter -i 95 -l 100 {1}/{2}.delta > {1}/{2}.delta.filt; "\
            "dnadiff -d {1}/{2}.delta.filt -p {1}/{2} &>>{1}/{2}.log 2>&1; "\
            "prodigal -i {1}/{2}.mc.fasta -d {1}/{2}.mc.fna -a {1}/{2}.mc.faa -q > /dev/null; "\
            "sleep 1s; "\
            "fetchMG.pl -m extraction {1}/{2}.mc.faa -o {1}/{2}.fetchMG > /dev/null; "\
            "cat {1}/{2}.fetchMG/*.fna > {1}/{2}.mc.marker.fna; "\
            "touch {1}/{2}.done".format(complete_ref,result_dir,genome_name)
         
         os.system(cmd)

#=========================================================================================
#create quality report file using the genome, gene and marker files
#=========================================================================================
def create_report_file(list_ref, result_dir, output_name):

   with open (output_name, 'wt') as output_file:
      output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("genome_name", 
                        "#contig", "size_theo", "size_obs", "perc_recovery", 
                        "perc_identity", "avg_alignment_len", "#snps", 
                        "#gsnps","#breakpoints","#gene","#complete_gene","#both_side_incomplete_gene",
                        "#marker","#complete_marker","both_side_incomplete_marker"))
      
      for ref in list_ref:
         genome_name = ".".join(ref.split('.')[:-1])
         report_file_name = result_dir + "/" + genome_name + ".report"
         gene = dict.fromkeys(["total","complete","lack_both"],0)
         marker = dict.fromkeys(["total","complete","lack_both"],0)
         genome = dict.fromkeys(["size_theo","contigs","size_obs",
       	       	       	   "perc_identity","perc_recovery", "avg_alignment_len",
                           "snps","gsnps","breakpoints"],0)

         
         #gene status
         with open("%s/%s.mc.fna"%(result_dir,genome_name)) as gene_file:
            for line in gene_file:
               if line[0]==">":
                  gene["total"] += 1
                  if "partial=00" in line:
                     gene["complete"] += 1
                  elif "partial=11" in line:
                     gene["lack_both"] += 1
         
         #marker status
         with open("%s/%s.mc.marker.fna"%(result_dir,genome_name)) as marker_file:
            for	line in	marker_file:
       	       if line[0]==">":
       	       	  marker["total"] += 1
       	       	  if "partial=00" in line:
       	       	     marker["complete"] += 1
       	       	  elif "partial=11" in line:
       	       	     marker["lack_both"] += 1
         
         with open(report_file_name,'rt') as report_file:
            for line in report_file:
               if "AlignedSeqs" in line:
                  genome["nb_contig"] = line.split()[2].strip().split('(')[0]
               elif "TotalBases" in line:
                  genome["size_theo"] = line.split()[1].strip()
               elif "AlignedBases" in line:
                  genome["size_obs"] = line.split()[1].strip().split('(')[0]
               elif ("AvgLength" in line) and (genome["avg_alignment_len"] == 0):
                  genome["avg_alignment_len"] = line.split()[1].strip()
               elif "TotalSNPs" in line:
                  genome["snps"] = line.split()[2].strip()
               elif "TotalGSNPs" in line:
                  genome["gsnps"] = line.split()[2].strip()
               elif "Breakpoints" in line:
                  genome["breakpoints"] = line.split()[2].strip()

            genome["perc_recovery"] = (float(int(genome["size_obs"])) / float(int(genome["size_theo"]))) * 100.0
            genome["perc_identity"] = (float(int(genome["size_obs"]) - int(genome["snps"])) / float(int(genome["size_obs"]))) * 100.0

            output_file.write("%s\t%s\t%s\t%s\t%.3f\t%.3f\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n" % (
                              genome_name, genome["nb_contig"], genome["size_theo"], genome["size_obs"], 
                              genome["perc_recovery"], genome["perc_identity"],
                              genome["avg_alignment_len"], genome["snps"], genome["gsnps"], genome["breakpoints"],
                              gene["total"],gene["complete"],gene["lack_both"],
                              marker["total"],marker["complete"],marker["lack_both"]))

   os.system("rm {0}/*.delta {0}/*.1delta {0}/*.qdiff {0}/*.rdiff {0}/*.snps"\
             " {0}/*.unqry {0}/*.mcoords {0}/*.mdelta {0}/*.delta.filt"\
             " {0}/*.1coords".format(result_dir))


#=========================================================================================
#create fasta file for each genome reconstructed by metacompass
#=========================================================================================
def create_contig_files_shared(contig_file, ref_folder, list_ref, result_dir,contig_list):
   genome_contig_list = {}

   #create ref db that contain all reference
   for genome in list_ref:
      genome_name = genome.replace(".fasta","")
      os.system("sed 's/>/>{0} /g' {1}/{0}.fasta >> {2}/reference.fasta".format(
                                           genome_name,ref_folder, result_dir))

   print('align sequences...')
   os.system("makeblastdb -in {0}/reference.fasta -out {0}/reference -dbtype nucl > /dev/null".format(result_dir))
   os.system("blastn -db {0}/reference -query {1} -word_size 100 -evalue 1e-2"\
                " -max_target_seqs 100 -perc_identity 95 -outfmt 6 -out {0}/contig_to_ref.blastn".format(
                result_dir, contig_file))

   print('extract contig reference all bests-hits')
   with open(result_dir + '/contig_to_ref.blastn', 'rt') as infile:
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
            contig_name = items[0].strip()
            ref_name = items[1].strip()
            if not genome_contig_list.get(ref_name):
               genome_contig_list[ref_name] = []

            if not contig_name in genome_contig_list[ref_name]:
               genome_contig_list[ref_name].append(contig_name)

         #keep the best score for same contig
         if update:
            oldquery=query
            oldscore=score

   print('create fasta file for each ref')
   for genome in genome_contig_list:
       filename = '{0}/{1}.mc.fasta'.format(result_dir,genome)

       with open(filename, 'wt') as genome_file:
          for contig in genome_contig_list[genome]:
             genome_file.write(">{0}\n{1}\n".format(contig,contig_list[contig]))

      
#=========================================================================================
#put contigs in memory
#=========================================================================================      
def extract_contig_list(contigfilename):
   contig_list = {}
   with open(contigfilename, 'rt') as contigfile:
      header=''
      seq=''
      for line in contigfile:
         if line[0]=='>':
            if header!='':
               contig_list[header] = seq
            header=line.strip().split()[0][1:]
            seq=''
         else:
            seq += line.strip()

      #last sequence
      contig_list[header] = seq

   return contig_list


#=========================================================================================
#MAIN
#=========================================================================================
def main():

   """Main program function
   """
               
   #Extract arguments
   args = config_parameters()

   if not os.path.isdir(args.result_dir):
      os.system("mkdir %s" % args.result_dir)

   list_ref = extract_list_ref(args.ref_folder)

   #extract contig sequence
   contig_list = extract_contig_list(args.contig_file)

   #create metaG contig files
   create_contig_files_shared(args.contig_file, args.ref_folder, list_ref, args.result_dir, contig_list)

   #compare each contigs to ref
   compare_contigs_to_ref(args.contig_file, args.ref_folder, list_ref, args.result_dir)      

   #create report file
   create_report_file(list_ref, args.result_dir, args.output_name)
                       
   print ('end')
                                                                   
if __name__ == '__main__':
    main()

