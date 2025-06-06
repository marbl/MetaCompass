import glob, os
import sys
import shutil
import argparse
import pathlib



parser = argparse.ArgumentParser(description="Reference culling")
parser.add_argument('--output_folder', '-o', required=True,
                    metavar='output',
                    help="Metacompass output folder")  # input file
parser.add_argument('--MIN_REFS', '-mr', required=True)
parser.add_argument('--REFGUIDED_CONTIGS', '-rc', required=True)
parser.add_argument('--ALIGNED_READS', '-ar', required=True)
parser.add_argument('--UNALIGNED_READS', '-ur', required=True) # may not be needed

print("Start Creating Final Output Folder")
args = parser.parse_args(sys.argv[1:])
print(args)


if not os.path.exists(args.output_folder + "/MetaCompass_Output"):
    os.makedirs(args.output_folder + "/MetaCompass_Output")
    print("The output directory is created!")
    if not os.path.exists(args.output_folder + "/MetaCompass_Output/contigs"):
        os.makedirs(args.output_folder + "/MetaCompass_Output/contigs")
        print("The contigs directory is created!")
    if not os.path.exists(args.output_folder + "/MetaCompass_Output/read_mapping"):
        os.makedirs(args.output_folder + "/MetaCompass_Output/read_mapping")
        print("The read_mapping directory is created!")

Min_Cand_Ref = open(args.MIN_REFS,"r").read().split("\n")

outf = None
try:
    outf = open(args.output_folder + '/MetaCompass_Output/final_reference_candidates.txt', "w")
except OSError as err:
    print("Cannot open outfile: {e}".format(e=err))
    exit(1)
# Update the minimum candidate list according to which references were skipped in align_reads/pilon
updated_candidate_list = []
count = 0
for i in Min_Cand_Ref:
    if i:
        if os.path.exists(args.REFGUIDED_CONTIGS + "/" + i + "/pilon/"):
            updated_candidate_list.append(i)
            outf.write(i + "\n")
            count = count + 1
# change this code to only run if references were found. Check min_ref_can.txt [wc -l] 
if count:
    print(os.system("ls " + args.REFGUIDED_CONTIGS + "/"))
    for i in updated_candidate_list:
        if i:
            shutil.copyfile(args.REFGUIDED_CONTIGS + "/" + i + "/pilon/contigs.pilon.fasta" , args.output_folder + "/MetaCompass_Output/contigs/"+i+".contigs.pilon.fasta")
            shutil.copyfile(args.REFGUIDED_CONTIGS + "/" + i + "/ref_guided_contig.agp" , args.output_folder + "/MetaCompass_Output/contigs/"+i+".ref_guided_contig.agp")
            shutil.copyfile(args.ALIGNED_READS + "/" + i + "/reads_mapped.fq" , args.output_folder + "/MetaCompass_Output/read_mapping/"+i+".reads_mapped.fq")
            if i==updated_candidate_list[-1]: # for the last reference, move the unmapped read files as well
                shutil.copyfile(args.REFGUIDED_CONTIGS + "/" + i + "/reads_unmapped.1.fq" , args.output_folder + "/MetaCompass_Output/read_mapping/reads_unmapped.1.fq")
                shutil.copyfile(args.REFGUIDED_CONTIGS + "/" + i + "/reads_unmapped.2.fq" , args.output_folder + "/MetaCompass_Output/read_mapping/reads_unmapped.2.fq")
                shutil.copyfile(args.REFGUIDED_CONTIGS + "/" + i + "/reads_unmapped.u.fq" , args.output_folder + "/MetaCompass_Output/read_mapping/reads_unmapped.u.fq")
else:
    print("\n *** No suitable reference candidates remained from reference culling, so no assembly was performed! ***")
