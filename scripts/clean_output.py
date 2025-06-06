import glob, os
import sys
import shutil
import argparse
import pathlib



parser = argparse.ArgumentParser(description="Reference culling")
parser.add_argument('--output_folder', '-o', required=True,
                    metavar='output',
                    help="Metacompass output folder")  # input file
parser.add_argument('--ref_selection', '-rs', action='store_true')
parser.add_argument('--ref_culling', '-rc', action='store_true')
parser.add_argument('--ref_assembly', '-ra', action='store_true')
parser.add_argument('--post_cleaning', '-pc', action='store_true')


args = parser.parse_args(sys.argv[1:])
print(args)

try:
    os.mkdir("tmp/")

    if (args.ref_selection):
        pathlib.Path("tmp/cluster_refs/").mkdir(parents=True, exist_ok=True)
        shutil.copy(args.output_folder + "/reference_selection/cluster_refs/reference_candidates.txt", "tmp/cluster_refs/reference_candidates.txt")
        shutil.copy(args.output_folder + "/reference_selection/filter_refs.log", "tmp/filter_refs.log")
        shutil.copy(args.output_folder + "/reference_selection/reference_selection.log", "tmp/reference_selection.log")
        os.system("rm -r " + args.output_folder + "/reference_selection/*")
        os.system("cp -r tmp/cluster_refs " + args.output_folder + "/reference_selection/cluster_refs")
        shutil.copy("tmp/filter_refs.log", args.output_folder + "/reference_selection/filter_refs.log")
        shutil.copy("tmp/reference_selection.log", args.output_folder + "/reference_selection/reference_selection.log")
        print("Delete ref_selection")

    if (args.ref_culling):
        # os.system("rm -r " + args.output_folder + "/reference_culling/collect_kmers/read_kmers/")
        # os.system("rm -r " + args.output_folder + "/reference_culling/collect_kmers/ref_kmers/")
        print("Delete ref_culling")

    if (args.ref_assembly):
        os.system("rm -r " + args.output_folder + "/reference_assembly/VALET/")
        print("Delete ref_assembly")

    if (args.post_cleaning):
        os.system("rm -r " + args.output_folder + "/reference_assembly/align_reads")
        print("Finish Post Process Cleaning")

except:
    os.system("rm -r tmp/")
