import argparse
import time
import csv
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path
from typing import List, Optional
from .tools import ToolsFactory, Mash, Pilon, BBtools, HeapArranger
from .analyze import DownstreamAnalyzer
from .utils import run_shell_cmd, align_with_minimap2, extract_reads, sam_to_sorted_bam, get_max_len, file_empty, \
    delete_file, get_boc, count_reads, concatenate_cluster_refs, write_list_to_file, move_file


class ReadAligner(object):
    """Class responsible for aligning and mapping read data to reference genomes.

    This class handles the process of reading input data, interleaving reads, and mapping clusters of reads
    to reference genomes. It utilizes various tools and methods to perform the alignment and mapping tasks."""

    def __init__(self):
        # Data storage
        self.inputs = {}  # holds input data
        self.sketch_dir = Path()  # directory for sketches
        self.cluster_mapped_reads: Path = Path()  # Path to last processed cluster's mapped reads
        self.cluster_mapped_reads_lines = 0
        self.mapped_clusters: Path = Path()  # Path to the txt file storing a list of mapped clusters
        self.mapped_genomes = []

        # Current alignment state
        self.curr_ref_align = {}  # Current reference alignment data
        self.curr_nz_align = {}  # Current non-zero alignment data
        self.curr_ref_name: str = "NA"  # Name of the current reference genome
        self.curr_ref_path: Path = Path()  # Path to the current reference genome
        self.curr_out_dir: Path = Path()  # Path to current reference's result output dir
        self.curr_unmapped = {}  # Current unmapped alignment data

        self.interleaved_reads: Path = Path()  # Path to interleaved reads

        # Dependencies
        self.tools_factory: Optional[ToolsFactory] = None  # Factory for different tool objects used in the script
        self.mash_obj: Optional[Mash] = None
        self.pilon_obj: Optional[Pilon] = None
        self.bbtools_obj: Optional[BBtools] = None
        self.heap_arranger: Optional[HeapArranger] = None
        self.downstream_analyzer: Optional[DownstreamAnalyzer] = DownstreamAnalyzer()

        # Misc
        self.total_bases = 0

    def debug_output(self, message: str):
        # Print debug messages if debugging is enabled
        if self.inputs["debug"]:
            timestamp = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            print(f"{timestamp},{message}")

    def set_interleaved_reads(self):
        # interleaved_reads check
        print("DEBUG: Checking interleaved reads")
        forward_reads = self.inputs.get("forward", None)
        reverse_reads = self.inputs.get("reverse", None)
        interleaved_reads = self.inputs.get("interleaved_reads", None)

        if forward_reads and reverse_reads:
            if not interleaved_reads:
                self.interleaved_reads = Path(self.inputs["out"]) / "interleaved.fq"  # assign
                self.interleave_reads()  # Generate
            else:
                self.interleaved_reads = Path(interleaved_reads)
        else:
            if not interleaved_reads:
                raise FileNotFoundError("Both forward and reverse reads are missing. Interleaved reads are required.")
            else:
                self.interleaved_reads = Path(interleaved_reads)

    def read_inputs(self):
        # Parse command line arguments and store them in self.inputs
        parser = argparse.ArgumentParser(
            description="Align reads to a list of references sequentially by aligning to first reference, "
                        "then taking unaligned reads and aligning to next reference until no reads or references left.")

        # Create argument groups for better organization
        input_group = parser.add_argument_group("Input Paths")
        dependency_group = parser.add_argument_group("Path to dependencies")
        output_group = parser.add_argument_group("Output Paths")

        # Input Paths
        input_group.add_argument("-f", "--forward", required=False, type=Path, help="path to forward reads")
        input_group.add_argument("-r", "--reverse", required=False, type=Path, help="path to reverse reads")
        input_group.add_argument("-ir", "--interleaved-reads", required=False, type=Path,
                                 help="Path to interleaved forward & reverse reads")
        input_group.add_argument("-u", "--unpaired", required=False, type=Path, help="path to unpaired reads")
        input_group.add_argument("-refs", "--references", required=False, type=Path,
                                 help="path to reference file directory")
        input_group.add_argument("-rs", "--ref-genome-sketch", required=True, type=Path,
                                 help="path to reference genome sketches")
        input_group.add_argument("-as", "--all-reads-sketch", required=True, type=Path,
                                 help="path to all reads sketches")
        input_group.add_argument("-cl", "--cluster-list", required=True, type=Path, help="path to cluster list")
        input_group.add_argument("-boc", "--breadth-of-coverage", required=False, type=float, default=5.0,
                                 help="required breadth of coverage")
        input_group.add_argument("-mcl", "--max-contig-length", required=False, type=float, default=10000,
                                 help="required max contig length")

        # Output Paths
        output_group.add_argument("-o", "--out", required=True, type=Path, help="path to output directory")

        # Additional arguments
        parser.add_argument("-t", "--threads", type=str, default="2", help="input threads arg to bowtie2")
        parser.add_argument("-debug", "--debug", action='store_true', help="toggle input for debug mode")

        args = parser.parse_args()

        for key in vars(args):
            self.inputs[key] = getattr(args, key)

        if "debug" not in self.inputs:
            self.inputs["debug"] = False

        # initialize dependencies
        self.init_dependencies()

        # set interleaved reads
        self.set_interleaved_reads()

        self.cluster_mapped_reads = Path(self.inputs["out"]) / "cluster_mapped.fq"
        self.mapped_clusters = Path(self.inputs["out"]) / "mapped_clusters.txt"

        if (self.inputs["debug"]):
            print(f"""
            Debug mode is on! 
            
            Cut-off parameters for mapping: 
                * Required breadth of coverage : {self.inputs["breadth_of_coverage"]}
                * Required contig max length: {self.inputs["max_contig_length"]}
            
            """)

            self.print_inputs()

        # sys.exit(0)

    def print_inputs(self):
        print("Inputs are:")
        for k in self.inputs:
            print(k, self.inputs[k])

    def init_dependencies(self):
        # Initializes tool objects like Mash, Pilon, etc.

        self.tools_factory = ToolsFactory(
            output_dir=self.inputs["out"],
            threads=str(self.inputs["threads"]),
        )

        try:
            self.mash_obj = self.tools_factory.get_mash(output_dir=self.inputs["out"])
            self.pilon_obj = self.tools_factory.get_pilon()
            self.bbtools_obj = self.tools_factory.get_bbtools()
            self.heap_arranger = self.tools_factory.get_heap_arranger()
        except FileNotFoundError as fe:
            print(fe)

    def init_curr_alignment_vars(self):
        """Used for initializing paths to current iteration's alignment data"""
        self.curr_ref_align = {
            "mapped_sam": Path(self.curr_out_dir) / "ref_align_mapped.sam",
            "mapped_bam": Path(self.curr_out_dir) / "ref_align_mapped.bam",
            "mapped_sorted": Path(self.curr_out_dir) / "ref_align_mapped_sorted.sam",
            "mapped_ids": Path(self.curr_out_dir) / "ref_align_mapped_ids.txt",
            "mapped_fq": Path(self.curr_out_dir) / "ref_align_mapped.fq",
        }

        self.curr_nz_align = {
            "contig_fasta": Path(self.curr_out_dir) / "nz_contig.fasta",
            "cleaned_fasta": Path(self.curr_out_dir) / "nz_contig_cleaned.fasta",
            "mapped_sam": Path(self.curr_out_dir) / "nz_contig_align_mapped.sam",
            "mapped_bam": Path(self.curr_out_dir) / "nz_contig_align_mapped.bam",
            "mapped_sorted": Path(self.curr_out_dir) / "nz_contig_align_mapped_sorted.sam",
            "mapped_ids": Path(self.curr_out_dir) / "nz_contig_align_mapped_ids.txt",
            "mapped_fq": Path(self.curr_out_dir) / "reads_mapped.fq",
            "unmapped_fq": Path(self.curr_out_dir) / "reads_unmapped.fq",
            "unmapped_bam": Path(self.curr_out_dir) / "nz_contig_align_unmapped.bam",
        }

    def init_unmapped_vars(self):
        """Used for initializing path to current iteration's unmapped alignment data."""
        self.curr_unmapped = {
            "unmapped_bam": Path(self.curr_out_dir) / "unmapped.bam",
            "unmapped_fq": Path(self.curr_out_dir) / "unmapped.fq",
        }

    def align_with_reference(self):
        """Performs DNA sequence read alignment with a reference genome using Minimap2."""
        start_time = time.time()

        align_with_minimap2(
            ref_genome_path=self.curr_ref_path.resolve().as_posix(),
            interleaved_reads_path=self.interleaved_reads.resolve().as_posix(),
            outfile=self.curr_ref_align["mapped_sam"],
            threads=self.inputs.get("threads", 1)
        )

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Alignment w Reference Genome,{time_elapsed}")

    def align_with_non_zero_contig(self):
        """Aligns mapped reads with a non-zero length contig of the reference genome using the minimap2 alignment tool."""
        start_time = time.time()

        align_with_minimap2(
            ref_genome_path=self.curr_nz_align["contig_fasta"].resolve().as_posix(),
            interleaved_reads_path=self.curr_ref_align["mapped_fq"].resolve().as_posix(),
            outfile=self.curr_nz_align["mapped_sam"].resolve().as_posix(),
            threads=self.inputs['threads']
        )

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Non-zero Contig Alignment,{time_elapsed}")

    def extract_non_zero_regions(self):

        depth_file = self.curr_out_dir / "depth.txt"
        nz_bed = self.curr_out_dir / "non_zero_region.bed"
        # Define the output file
        cord_file = self.curr_out_dir / "cord.txt"

        # Generate depth file
        cmd = ["samtools", "depth", self.curr_ref_align["mapped_bam"].resolve().as_posix()]

        run_shell_cmd(cmd, depth_file)

        if depth_file.stat().st_size == 0:
            return

            # Generate non-zero region bed file

        bed_cmd = ["awk", '{print $1"\t"($2)-1"\t"$2"\t"$3}', depth_file]
        run_shell_cmd(bed_cmd, nz_bed)

        # Merge non-zero regions and generate coordinates file
        # Define the bedtools command and its arguments
        merge_cmd = ["bedtools", "merge", "-i", nz_bed.resolve().as_posix()]
        # Define the awk command and its arguments
        awk_cmd = ["awk", '{print $1":"($2)+1"-"$3}']
        run_shell_cmd([merge_cmd, awk_cmd], cord_file)

        # Extract non-zero regions from reference contig
        extract_cmd = ["samtools", "faidx", self.curr_ref_path.resolve().as_posix(), "-r",
                       cord_file.resolve().as_posix()]

        temp_file_1 = self.curr_out_dir / "temp_nz_filt_1.fasta"
        run_shell_cmd(extract_cmd, temp_file_1.resolve().as_posix())

        # Modify the headers in the fasta file
        # cmd = ["sed", "-i", 's/:/_/g ; s/-/_/g ;  s/\./_/g', temp_file_1.resolve().as_posix()]
        cmd = ["sed", "-i", 's/:/_/g ; s/-/_/g ;  s/\\./_/g', temp_file_1.resolve().as_posix()]
        run_shell_cmd(cmd)

        temp_file = self.curr_out_dir / "temp_nz_filt.fasta"

        filter_cmd = ["seqkit", "seq", "-m", "500", temp_file_1.resolve().as_posix()]
        print("DEBUG: running: " + " ".join(filter_cmd))

        run_shell_cmd(filter_cmd, temp_file.resolve().as_posix())

        # overwrite
        self.curr_nz_align["contig_fasta"].touch()
        if temp_file.exists():
            self.curr_nz_align["contig_fasta"].write_text(temp_file.read_text())
        # temp_file.unlink()

        # cord_file.unlink()

    def extract_mapped_ref_alignment_reads(self):
        """Extract mapped reads from reference alignment"""
        start_time = time.time()
        # Extract mapped reads
        extract_reads(
            mapped_sam=self.curr_ref_align["mapped_sam"].resolve().as_posix(),
            mapped_ids=self.curr_ref_align["mapped_ids"].resolve().as_posix(),
            interleaved_reads=self.interleaved_reads.resolve().as_posix(),
            outfile=self.curr_ref_align["mapped_fq"].resolve().as_posix(),
            threads=self.inputs["threads"], mapped=True)

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Mapped reads extraction,{time_elapsed}")

    def extract_mapped_nz_alignment_reads(self):
        """Extract mapped reads from non-zero contig alignment"""
        start_time = time.time()
        # outfile=self.curr_nz_align["mapped_fq"].resolve().as_posix(),
        # Extract mapped reads

        # the mapped reads after non-zero contig alignment are appended to cluster_mapped_reads.fq file
        num_lines = extract_reads(
            mapped_sam=self.curr_nz_align["mapped_sam"].resolve().as_posix(),
            mapped_ids=self.curr_nz_align["mapped_ids"].resolve().as_posix(),
            interleaved_reads=self.curr_ref_align["mapped_fq"].resolve().as_posix(),
            outfile=self.cluster_mapped_reads.resolve().as_posix(),
            threads=self.inputs["threads"], mapped=True)

        self.cluster_mapped_reads_lines += num_lines

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Mapped reads extraction,{time_elapsed}")

    def extract_unmapped_nz_alignment_reads(self):
        """Extract unmapped reads from non-zero contig alignment"""
        start_time = time.time()

        # Extract unmapped mapped reads
        extract_reads(
            mapped_sam=self.curr_nz_align["mapped_sam"].resolve().as_posix(),
            mapped_ids=self.curr_nz_align["mapped_ids"].resolve().as_posix(),
            interleaved_reads=self.interleaved_reads.resolve().as_posix(),
            outfile=self.curr_unmapped["unmapped_fq"].resolve().as_posix(),
            threads=self.inputs["threads"], mapped=False)

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Unmapped reads extraction,{time_elapsed}")

    def interleave_reads(self):
        """Interleave input forward and reverse reads using bbtools"""
        print("DEBUG: Interleaving reads")
        self.bbtools_obj.interleave_reads(
            forward_read=Path(self.inputs['forward']),
            reverse_read=Path(self.inputs['reverse']),
            outfile=self.interleaved_reads
        )

    def pilon_correction(self):
        """Performs error correction of a non-zero length contig using the Pilon tool."""

        start_time = time.time()
        curr_pilon_dir = self.curr_out_dir / "pilon"

        self.pilon_obj.correct(
            genome=self.curr_nz_align["cleaned_fasta"],
            pe_frag_bam=self.curr_nz_align["mapped_bam"],
            output_dir=curr_pilon_dir
        )

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Pilon,{time_elapsed}")

    def downstream_analysis(self, seq_list_file: Path):
        """Perform downstream analysis on the aligned reads

        Args:
            seq_list_file (Path): Path to the list of sequences to be included in the AGP file.
        """
        start_time = time.time()

        self.downstream_analyzer.curr_out_dir = self.curr_out_dir
        self.downstream_analyzer.stats_file = self.curr_out_dir / "stats"
        self.downstream_analyzer.analyze(seq_list_file, self.curr_nz_align, self.curr_ref_path)

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"{self.curr_ref_name},Downstream Analysis,{time_elapsed}")

    def cleanup(self, is_last_iteration: bool = False):
        # TODO: Figure out outputs
        """Clears unnecessary data after each iteration."""
        # Delete fq files only if it is not the last iteration
        out_dir = self.curr_out_dir.resolve().as_posix()
        if not is_last_iteration:
            delete_file(self.curr_ref_align['mapped_fq'].resolve().as_posix())
            delete_file(self.curr_nz_align['mapped_fq'].resolve().as_posix())

        # Delete other files
        files_to_remove = ["*.bai", "*.bed", "*.bt1", "*.bt2", "seq_list.txt", "cord.txt", "*.sam", "*.bam"]
        for ext in files_to_remove:
            delete_file(f"{out_dir}/{ext}")

    def get_next_best_cluster(self, cluster_list: List[List]):
        """Utilizes kmc tools to find the cluster which needs to be mapped next.

        Args:
            cluster_list (List[List]): Contains path to reference sequences in each cluster.
        """

        selected_clusters_store = defaultdict(bool)
        total_clusters = len(cluster_list)

        for _ in range(total_clusters - 1):
            start_time = time.time()
            selected_cluster_kmer, cluster_num = self.heap_arranger.select_cluster(
                all_reads_sketches=self.inputs["all_reads_sketch"],
                curr_subset_reads=self.cluster_mapped_reads,
                ref_genome_sketches=self.inputs["ref_genome_sketch"],
                selected_clusters=self.mapped_clusters
            )

            selected_clusters_store[str(cluster_num)] = True

            selected_cluster_refs = cluster_list[cluster_num]

            end_time = time.time()
            time_elapsed = end_time - start_time
            self.debug_output(f"Cluster {cluster_num + 1}, Time to find {time_elapsed}")

            yield cluster_num + 1, selected_cluster_refs

        start_time = time.time()
        for i in range(total_clusters):
            if not selected_clusters_store[str(i)]:
                end_time = time.time()
                selected_cluster_refs = cluster_list[i]
                time_elapsed = end_time - start_time
                self.debug_output(f"Cluster {i + 1}, Time to find {time_elapsed}")
                yield i + 1, selected_cluster_refs
                break

    def map_cluster_references(self, cluster_refs: List[str] = []):
        """Maps interleaved reads to a list of reference genomes iteratively,
        performing alignment, processing, correction, and downstream analysis.

        Args:
            cluster_refs (List): a list of sequences present in the cluster
        """

        total_iter = len(cluster_refs)
        for current_iter, ref_path in enumerate(cluster_refs):
            start_time = time.time()

            # Set current reference-related attributes and directories
            self.curr_ref_name = Path(ref_path).name
            self.curr_ref_path = Path(ref_path)
            self.curr_out_dir = Path(self.inputs["out"]) / self.curr_ref_name
            os.makedirs(self.curr_out_dir.resolve().as_posix(), exist_ok=True)

            self.init_curr_alignment_vars()

            # if the interleaved reads is empty, exit
            if file_empty(self.interleaved_reads):
                self.exit_map_cluster_reference(start_time, current_iter, total_iter, msg="Skip:Interleaved File empty")
                return

            # Align interleaved reads to the reference genome
            self.align_with_reference()

            # Convert mapped sam output from reference alignment to bam & sort it positionally
            sam_to_sorted_bam(
                input_sam=self.curr_ref_align["mapped_sam"].resolve().as_posix(),
                output_bam=self.curr_ref_align["mapped_bam"].resolve().as_posix(),
                threads=self.inputs['threads']
            )

            # Extract mapped reads from the reference alignment
            self.extract_mapped_ref_alignment_reads()

            self.extract_non_zero_regions()

            if file_empty(self.curr_nz_align["contig_fasta"]):
                self.exit_map_cluster_reference(start_time, current_iter, total_iter, msg="Skip: NZ fasta empty")
                return

            # cut-off conditions
            curr_boc = get_boc(self.curr_nz_align["contig_fasta"], self.curr_ref_path)
            is_boc_less = curr_boc < self.inputs["breadth_of_coverage"]
            curr_max_len = get_max_len(self.curr_nz_align["contig_fasta"])
            is_max_len_less = curr_max_len < self.inputs["max_contig_length"]

            self.debug_output(f"{self.curr_ref_name},BOC,{curr_boc}")
            self.debug_output(f"{self.curr_ref_name},MCL,{curr_max_len}")
            continue_processing = not (is_boc_less and is_max_len_less) or current_iter == 0

            if not continue_processing:
                self.exit_map_cluster_reference(start_time, current_iter, total_iter,
                                                msg="Skip: boc/mcl conditions not satisfied")
                return

            # Align mapped reads to a non-zero length contig

            self.align_with_non_zero_contig()
            # Convert mapped sam output from non-zero reference alignment to bam & sort it positionally
            sam_to_sorted_bam(
                input_sam=self.curr_nz_align["mapped_sam"].resolve().as_posix(),
                output_bam=self.curr_nz_align["mapped_bam"].resolve().as_posix(),
                threads=self.inputs['threads']
            )

            self.init_unmapped_vars()

            # Extract mapped and unmapped non-zero alignment reads
            self.extract_mapped_nz_alignment_reads()
            self.extract_unmapped_nz_alignment_reads()

            if file_empty(self.curr_unmapped["unmapped_fq"]):
                self.exit_map_cluster_reference(start_time, current_iter, total_iter,
                                                msg="Complete: No more unmapped reads")
                return

            # Delete the current set of interleaved reads
            self.interleaved_reads.unlink()

            # Set the extracted unmapped reads as interleaved reads for the next iteration
            self.interleaved_reads = self.curr_unmapped["unmapped_fq"]

            # Calculate sequence list

            cmd1 = ["samtools", "depth", self.curr_nz_align["mapped_bam"].resolve().as_posix(), "-J"]
            cmd2 = ["awk", '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}']
            cmd3 = ["sed", 's/_/\t/g']
            cmd4 = ["awk", '{if($5>0 ){print $1"_"$2"_"$3"_"$4}}']
            print("DEBUG: running: " + " ".join(cmd1))
            print("DEBUG: running: " + " ".join(cmd2))
            print("DEBUG: running: " + " ".join(cmd3))
            print("DEBUG: running: " + " ".join(cmd4))
            seq_list_file = self.curr_out_dir / "seq_list.txt"
            run_shell_cmd([cmd1, cmd2, cmd3, cmd4], seq_list_file)

            # Extract contigs with non-zero coverage from the reference FASTA
            cmd5 = ["seqkit", "grep", "-f", seq_list_file.resolve().as_posix(),
                    self.curr_nz_align["contig_fasta"].resolve().as_posix()]

            print("DEBUG: running: " + " ".join(cmd5))

            run_shell_cmd(cmd5, self.curr_nz_align["cleaned_fasta"].resolve().as_posix())

            self.pilon_correction()
            self.downstream_analysis(seq_list_file)

            self.exit_map_cluster_reference(start_time, current_iter, total_iter, msg="Complete")
            # sys.exit(0)

    def exit_map_cluster_reference(self, iter_start_time, current_iter, total_iter, msg: str = ""):
        """Used to exit an iteration in map_cluster_reference.

        Logs the time taken by the iteration and cleans up excess files

        Args:
            iter_start_time (float) : Starting time of iteration.
            current_iter (int): current iteration number
            total_iter (int): total number of possible iterations i.e. total # of cluster refs
        """
        end_time = time.time()
        time_elapsed = end_time - iter_start_time
        self.debug_output(f"{self.curr_ref_name},Iteration:{msg},{time_elapsed}")
        if "Skip" in msg:
            shutil.rmtree(self.curr_out_dir.resolve().as_posix())
        else:
            self.mapped_genomes.append(self.curr_ref_name)
            self.cleanup(current_iter == total_iter - 1)

    def map_cluster(self, cluster_num: int = 0, cluster_refs: List[str] = []):
        """Map reads for a specific cluster

        Args:
            cluster_num (int) : The cluster which is being mapped
            cluster_refs (List): a list of sequences present in the cluster
        """
        start_time = time.time()

        # Clear the mapped reads from previous cluster
        # if (self.cluster_mapped_reads.exists()):
        #    self.cluster_mapped_reads.unlink()

        # Create dirs for storing mash outputs for the current cluster
        _, ref_msh_path, screen_out_path = self.mash_obj.create_cluster_dirs(cluster_num)

        # Sketch and screen references
        self.mash_obj.sketch_references(cluster_refs, ref_msh_path)

        self.mash_obj.screen_references(self.interleaved_reads, ref_msh_path, screen_out_path)
        # Get the mash order and map references accordingly
        ref_mash_order = self.mash_obj.get_mash_order(screen_out_path)

        self.map_cluster_references(ref_mash_order)

        end_time = time.time()
        time_elapsed = end_time - start_time
        self.debug_output(f"Cluster {cluster_num},Total time,{time_elapsed}")

    def reduce_interleave_reads(self, cluster_list: List[List]):

        start = time.time()
        concat_files = {
            "fna": self.inputs["out"] / "concat_refs.fna",
            "mapped_sam": self.inputs["out"] / "concat_refs_mapped.sam",
            "mapped_ids": self.inputs["out"] / "concat_refs_mapped_ids.txt",
            "mapped_fq": self.inputs["out"] / "concat_refs_mapped.fq"
        }
        # Concatenate all references
        #concatenate_cluster_refs(cluster_list, concat_files["fna"])
        # Align with minimap2

        #align_with_minimap2(
        #    ref_genome_path=concat_files["fna"].resolve().as_posix(),
        #    interleaved_reads_path=self.interleaved_reads.resolve().as_posix(),
        #    outfile=concat_files["mapped_sam"],
        #    threads=self.inputs.get("threads", 1)
        #)

        # Extract mapped reads
        #extract_reads(
        #    mapped_sam=concat_files["mapped_sam"].resolve().as_posix(),
        #    mapped_ids=concat_files["mapped_ids"].resolve().as_posix(),
        #    interleaved_reads=self.interleaved_reads.resolve().as_posix(),
        #    outfile=concat_files["mapped_fq"].resolve().as_posix(),
        #    threads=self.inputs["threads"], mapped=True)

        # Delete interleaved files
        #self.interleaved_reads.unlink()
        # Delete the sam file
        #concat_files["mapped_sam"].unlink()
        # Assign new interleaved reads
        self.interleaved_reads = concat_files["mapped_fq"]
        end = time.time()
        self.debug_output(f"Time to reduce interleaved reads {end - start}")

    def map_clusters(self):
        """Map reads for all clusters"""
        os.makedirs(self.sketch_dir.resolve().as_posix(), exist_ok=True)

        #
        total_clusters = 0
        with open(self.inputs["cluster_list"], "r") as cluster_file:
            cluster_list: List[List] = list(csv.reader(cluster_file))
            self.reduce_interleave_reads(cluster_list)

            for selected_cluster_num, selected_cluster_refs in self.get_next_best_cluster(cluster_list):
                total_clusters += 1
                self.map_cluster(cluster_num=selected_cluster_num, cluster_refs=selected_cluster_refs)

        mapped_genomes_file_path = self.inputs["out"] / "mapped_genomes.txt"
        write_list_to_file(self.mapped_genomes, mapped_genomes_file_path)

        # print the total number of clusters
        self.debug_output(f"Total number of clusters: {total_clusters}")
        num_mapped_reads = self.cluster_mapped_reads_lines // 4
        self.debug_output(f"Number of mapped reads in all clusters: {num_mapped_reads}")

        # copy unmapped reads to the src file
        unmapped_reads_src = self.curr_unmapped["unmapped_fq"]
        unmapped_reads_dst = unmapped_reads_src.parent.parent / "unmapped.fq"
        move_file(unmapped_reads_src, unmapped_reads_dst)
