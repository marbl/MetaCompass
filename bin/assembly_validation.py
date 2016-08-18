#!/usr/bin/env python
# -*- coding: utf-8 -*-

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


def extract_list_ref(ref_folder):

   list_ref = []
   list_ref = [f for f in listdir(ref_folder) if os.path.isfile(join(ref_folder, f))]

   return list_ref

def compare_contigs_to_ref(contig_file, ref_folder, list_ref, result_dir):

   os.system("touch %s/ready.done" % result_dir)
   for ref in list_ref:
      complete_ref = ref_folder + ".".join(ref.split('.')[:-1]) + ".fasta" 
      genome_name = ".".join(ref.split('.')[:-1])
      print (genome_name)
      if not os.path.isfile("%s/%s.done" % (result_dir, genome_name)): 
         #nucmer comparison
         cmd = "tsub 'nucmer -l 35 {0} {1} -p {2}/{3} &>{2}/{3}.log; "\
            "delta-filter -i 95 -l 35 {2}/{3}.delta > {2}/{3}.delta.filt; "\
            "dnadiff -d {2}/{3}.delta.filt -p {2}/{3} &>>{2}/{3}.log 2>&1; "\
            "touch {2}/{3}.done' -q throughput -l nodes=1:ppn=1,mem=10gb"\
            ",walltime=00:30:00".format(
            complete_ref, contig_file, result_dir, genome_name)
         os.system(cmd)


def create_report_file(list_ref, result_dir, output_name):

   #wait till all report files are generated
   nb_genome = len(list_ref)

   process_cmd = subprocess.Popen('ls {0}/*.done | wc -l'.format(result_dir),
                                  shell=True, stdout=subprocess.PIPE)
   
   nb_done = int(process_cmd.communicate()[0])
   
   while nb_done < nb_genome+1:
      print ("waiting for report files...")
      time.sleep(10)
      process_cmd = subprocess.Popen('ls {0}/*.done | wc -l'.format(result_dir),
                                  shell=True, stdout=subprocess.PIPE)
      nb_done = int(process_cmd.communicate()[0])

   with open (output_name, 'wt') as output_file:
      output_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("genome_name", 
                        "nb_contig", "total_ref_bases", "aligned_bases", "%_ref_recovered", 
                        "%_avg_identity", "avg_alignment_len", "total_snps", 
                        "total_gsnps"))
      
      for ref in list_ref:
         genome_name = ".".join(ref.split('.')[:-1])
         report_file_name = result_dir + "/" + genome_name + ".report"
         nb_contig="NA"
         total_ref_bases="NA"
         perc_recovered_ref = 0.0
         aligned_bases = "NA"
         avg_identity = "NA"
         avg_alignment_len = "NA"
         total_snps= "NA"
         total_gsnps = "NA"
 
         with open(report_file_name,'rt') as report_file:
            for line in report_file:
               if "TotalSeqs" in line:
                  nb_contig = line.split()[2].strip()
               elif "TotalBases" in line:
                  total_ref_bases = line.split()[1].strip()
               elif "AlignedBases" in line:
                  aligned_bases = line.split()[1].strip().split('(')[0]
               elif ("AvgLength" in line) and (avg_alignment_len != "NA"):
                  avg_alignment_len = line.split()[1].strip()
               elif ("AvgIdentity" in line) and (avg_identity != "NA"):
                  avg_identity = line.split()[2].strip()
               elif "TotalSNPs" in line:
                  total_snps = line.split()[2].strip()
               elif "TotalGSNPs" in line:
                  total_gsnps = line.split()[2].strip()
            perc_recovered_ref = (float(aligned_bases) / float(total_ref_bases)) * 100
            output_file.write("%s\t%s\t%s\t%s\t%.4f\t%s\t%s\t%s\t%s\n" % (
                              genome_name, nb_contig, total_ref_bases, aligned_bases, 
                              perc_recovered_ref, avg_identity,
                              avg_alignment_len, total_snps, total_gsnps))

def main():

   """Main program function
   """
               
   #Extract arguments
   args = config_parameters()

   if not os.path.isdir(args.result_dir):
      os.system("mkdir %s" % args.result_dir)

   list_ref = extract_list_ref(args.ref_folder)

   #compare each contigs to ref
   compare_contigs_to_ref(args.contig_file, args.ref_folder, list_ref, args.result_dir)      

   #create report file
   create_report_file(list_ref, args.result_dir, args.output_name)
                       
   print ('end')
                                                                   
if __name__ == '__main__':
    main()

