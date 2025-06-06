import os
from typing import List
from pathlib import Path
from .utils import *
from shutil import which


class Mash(object):
    def __init__(self, bin_path: Path, output_dir: Path, threads: str):
        """Initialize a Mash object.

        Args:
            bin_path (Path): Path to the directory containing Mash binaries.
            output_dir (Path): Path to the output directory for Mash-related data.
            threads (str): Number of threads for parallel processing.

        Raises:
            FileNotFoundError: If the Mash binary is not found in the specified path.
        """
        self.bin_folder = bin_path
        self.bin = "mash" # assume it is in path
        self.output_dir = Path(output_dir)
        self.sketch_dir = self.output_dir / "sketches"
        self.threads = threads

        if not which(self.bin):
            raise FileNotFoundError(f"Binary not found: {self.bin}")
#        if not self.bin.exists():
#            raise FileNotFoundError(f"Binary not found: {self.bin}")

    def sketch_references(self, cluster_refs: List[str], ref_msh_path: Path):
        """Generate Mash sketches for reference sequences.

        Args:
            cluster_refs (List[str]): List of paths to reference sequences in the cluster.
            ref_msh_path (Path): Path to the output reference Mash sketch file.
        """

        cmd = [self.bin, "sketch", "-p", self.threads, "-o", ref_msh_path] + cluster_refs

        run_shell_cmd(cmd)

    def screen_references(self, interleaved_path: Path, ref_msh_path: Path, screen_out_path: Path):
        """Screen reference genomes against a Mash sketch to identify similarities.

        Args:
            interleaved_path (Path): Path to interleaved reads (FQ) file
            ref_msh_path (Path): Path to the reference Mash sketch file.
            screen_out_path (Path): Path to the output screen results file.
        """

        print("Reference genomes path", interleaved_path.resolve().as_posix())
        cmd = [self.bin, "screen", "-w", "-p", self.threads, ref_msh_path.resolve().as_posix(),
               interleaved_path.resolve().as_posix()]

        run_shell_cmd(cmd, screen_out_path.resolve().as_posix())

    @staticmethod
    def get_mash_order(screen_out_path: Path):
        """Extract the order of reference genomes from a Mash screen results file.

        Args:
            screen_out_path (Path): Path to the Mash screen results file.

        Returns:
            List[str]: List of paths to reference genomes in the order of appearance.
        """

        mash_order = []

        with open(screen_out_path, 'r') as file:
            for line in file:
                if not line.startswith('#query-ID'):
                    mash_out = line.split()  # Extract 5th column
                    score_and_path = (float(mash_out[0]), mash_out[4])
                    mash_order.append(score_and_path)

        mash_order.sort(key=lambda x: x[0], reverse=True)  # Sort by the first element of each tuple

        ref_order = [item[1] for item in mash_order]

        return ref_order

    def create_cluster_dirs(self, cluster_num: int):
        """Create directories for storing Mash outputs for the current cluster.

        Args:
            cluster_num (int): The cluster number.

        Returns:
            Tuple[Path, Path, Path]: Tuple containing paths to the cluster sketch directory,
            reference Mash sketch file, and screen results file.
        """
        # Create dirs for storing mash outputs for the current cluster
        cluster_sketch_dir = Path(self.sketch_dir) / f"sketch_{cluster_num}"
        os.makedirs(cluster_sketch_dir.resolve().as_posix(), exist_ok=True)
        ref_msh_path = Path(cluster_sketch_dir) / "ref.msh"
        screen_out_path = Path(cluster_sketch_dir) / "screen.txt"

        return cluster_sketch_dir, ref_msh_path, screen_out_path


class Pilon(object):
    def __init__(self, bin_path: Path = Path(), threads=""):
        """Initialize a Pilon object.

        Args:
            bin_path (Path, optional): Path to the directory containing Pilon binaries. Default is current directory.
            threads (str, optional): Number of threads for parallel processing. Default is an empty string.
        """
        self.bin_folder = bin_path
        self.bin_path = "pilon" # assume it's in path
        self.threads = str(threads)

        if not which(self.bin_path):
            raise FileNotFoundError(f"Binary not found: {self.bin_path}")
        #if not self.bin_path.exists():
        #    raise FileNotFoundError(f"Binary not found: {self.bin_path}")

    def correct(self, genome: Path, pe_frag_bam: Path, output_dir: Path):
        # TODO: discuss with Tu (not all files required)
        """Correct errors in a genome assembly using Pilon.

        Args:
            genome (Path): Path to the genome assembly file.
            pe_frag_bam (Path): Path to the paired-end fragment BAM file.
            output_dir (Path): Path to the output directory for corrected results.
        """
        os.makedirs(output_dir.resolve().as_posix(), exist_ok=True)

        cmd = [
            self.bin_path,
            "--flank", "5",
            "--threads", self.threads,
            "--mindepth", "3",
            "--genome", genome.resolve().as_posix(),
            "--frags", pe_frag_bam.resolve().as_posix(),
            "--output", "contigs.pilon",
            "--outdir", output_dir.resolve().as_posix(),
            "--fix", "bases,amb",
            "--tracks",
            "--changes"
        ]
        run_shell_cmd(cmd, "/dev/null")


class BBtools(object):
    def __init__(self, bin_path: Path, threads: str = ""):
        """Initialize a BBtools object.

        Args:
            bin_path (Path): Path to the directory containing BBtools binaries.
            threads (str, optional): Number of threads for parallel processing. Default is an empty string.

        Raises:
            FileNotFoundError: If the BBtools binary is not found in the specified path.
        """
        self.bin_folder = bin_path
        self.bin_path = "bbmap.sh"  # assume it's in path
        self.threads = threads

#        if not self.bin_path.exists():
        if not which(self.bin_path):
            raise FileNotFoundError(f"Binary not found: {self.bin_path}")

    def interleave_reads(self, forward_read: Path, reverse_read: Path, outfile: Path):
        """Interleave paired-end reads from forward and reverse read files.

        Args:
            forward_read (Path): Path to the forward read file.
            reverse_read (Path): Path to the reverse read file.
            outfile (Path): Path to the output interleaved read file.
        """
        print("DEBUG: writing interleave reads to: " + outfile.resolve().as_posix())
        cmd = [
            f"reformat.sh",
            f"t={self.threads}",
            "overwrite=f",
            f"in1={forward_read.resolve().as_posix()}",
            f"in2={reverse_read.resolve().as_posix()}",
            f"out={outfile.resolve().as_posix()}",

        ]
        run_shell_cmd(cmd)


class KMC(object):
    def __init__(self, bin_path=Path()):
        self.bin_path = bin_path
        self.heap_script = self.bin_path / "heap.sh"

    def find_best_kmer_reference(self, all_reads, subset_reads, fasta_dir, selected_refs):
        kmc_cmd = ["bash", self.heap_script.resolve().as_posix(),
                   all_reads.resolve().as_posix(),
                   subset_reads.resolve().as_posix(),
                   fasta_dir.resolve().as_posix(),
                   selected_refs.resolve().as_posix(),
                   self.bin_path.resolve().as_posix()
                   ]

        run_shell_cmd(kmc_cmd)


class HeapArranger(object):
    def __init__(self, bin_path: Path = Path()):
        self.kmc_path = "kmc"  # Path to kmc binary
        self.heap_script = bin_path / "heap.sh"  # Path to the heap script

    def select_cluster(self, all_reads_sketches: Path, curr_subset_reads: Path, ref_genome_sketches: Path,
                       selected_clusters: Path):
        """
        Selects a cluster using the heap script with specified input files.

        Args:
            all_reads_sketches (Path): Path to sketches of all reads.
            curr_subset_reads (Path): Path to the current subset of reads.
            ref_genome_sketches (Path): Path to sketches of reference genomes.
            selected_clusters (Path): Path to store selected clusters information.

        Returns:
            tuple: A tuple containing the selected cluster's kmer and cluster number.
        """
        cmd = ["bash", self.heap_script.resolve().as_posix(),
               all_reads_sketches.resolve().as_posix(),
               curr_subset_reads.resolve().as_posix(),
               ref_genome_sketches.resolve().as_posix(),
               selected_clusters.resolve().as_posix(),
               self.kmc_path]

        print("DEBUG running command: " + " ".join(cmd))

        captured_output = run_shell_cmd(cmd).splitlines()[-1]

        print("DEBUG output: " + captured_output)
        selected_cluster_kmer, _ = map(lambda s: s.split(":")[1].strip(), captured_output.split(","))
        print("DEBUG selected kmer: "  + selected_cluster_kmer)
        cluster_num = int(selected_cluster_kmer.split(".")[0].split("_")[1]) - 1

        return selected_cluster_kmer, cluster_num


class ToolsFactory(object):
    def __init__(self, output_dir: Path, threads: str = ""):
        """Initialize a ToolsFactory object.

        Args:
            bin_path (Path): Path to the directory containing tool binaries.
            output_dir (Path): Default output directory for tool-related data.
            threads (str, optional): Default number of threads for parallel processing. Default is an empty string.
        """
        self.default_threads = threads
        self.bin_path = Path(__file__).parent.parent.parent / "bin"
        self.default_output_dir = output_dir

    def get_mash(self, output_dir: Path = Path(), threads=""):
        """Get a Mash object with specified or default settings.

        Args:
            output_dir (Path, optional): Output directory for Mash-related data. Default is default_output_dir.
            threads (str, optional): Number of threads for parallel processing. Default is default_threads.

        Returns:
            Mash: A Mash object configured with the provided or default settings.
        """
        out_dir = output_dir if output_dir else self.default_output_dir
        t = threads if threads else self.default_threads

        return Mash(
            bin_path=self.bin_path,
            output_dir=out_dir,
            threads=t
        )

    def get_bbtools(self, threads=""):
        """Get a BBtools object with specified or default settings.

        Args:
            threads (str, optional): Number of threads for parallel processing. Default is default_threads.

        Returns:
            BBtools: A BBtools object configured with the provided or default settings.
        """
        t = threads if threads else self.default_threads

        return BBtools(
            bin_path=self.bin_path,
            threads=t
        )

    def get_pilon(self, threads=""):
        """
        Get a Pilon object with specified or default settings.

        Args:
            threads (str, optional): Number of threads for parallel processing. Default is default_threads.

        Returns:
            Pilon: A Pilon object configured with the provided or default settings.
        """
        t = threads if threads else self.default_threads

        return Pilon(
            bin_path=self.bin_path,
            threads=t
        )

    def get_heap_arranger(self):
        return HeapArranger(bin_path=self.bin_path)
