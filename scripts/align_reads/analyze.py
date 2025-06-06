from pathlib import Path
from .utils import run_shell_cmd


class DownstreamAnalyzer(object):
    """Class for performing downstream analysis on aligned reads.

    This class provides methods for running various analysis scripts on aligned reads,
    such as generating AGP files and calculating N50 statistics."""

    def __init__(self, scripts=Path()):
        """Initialize a DownstreamAnalyzer object.

        Args:
            scripts (Path, optional): Path to the directory containing analysis scripts. Default is an empty path.
        """
        self.scripts = Path(__file__).parent.parent
        self.agp_script = scripts / "AGP_generator.py"
        self.n50_script = scripts / "n50.py"
        self.curr_out_dir = Path()
        self.stats_file = Path()

    def generate_agp_file(self, seq_list_file: Path):
        """Generate AGP file using the AGP generator script.

        Args:
            seq_list_file (Path): Path to the list of sequences to be included in the AGP file.

        Note:
            This method sorts the sequence list, runs the AGP generator script, and removes the temporary sorted file.
        """
        sort_cmd = ["sort", "-t_", "-k1,1", "-k2,2", "-k3,3n", seq_list_file.resolve().as_posix()]
        seq_list_sorted = self.curr_out_dir / "seq_list_sorted.txt"
        run_shell_cmd(sort_cmd, seq_list_sorted.resolve().as_posix())

        agp_cmd = ["python", self.agp_script.resolve().as_posix(), "-o", self.curr_out_dir.resolve().as_posix(), "-c",
                   seq_list_sorted.resolve().as_posix()]
        run_shell_cmd(agp_cmd)
        seq_list_sorted.unlink()

    def calculate_average_depth(self, curr_nz_align: dict):
        """Calculate the average depth of coverage for the mapped reads.

        Args:
            curr_nz_align (dict): Dictionary containing non-zero alignment data.

        Returns:
            float: The calculated average depth of coverage.
        """
        # First command

        cmd1 = ["samtools", "depth", curr_nz_align["mapped_bam"].resolve().as_posix(), "-J"]

        # Second command
        cmd2 = ["awk", '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}']

        # Third command
        cmd3 = ["sed", "s/_/\t/g"]

        # Fourth command
        cmd4 = ["awk", '{s=s+$4-$3+1;g=g+$5} END {print g/s}']

        # Combining all commands for the pipeline
        cmds = [cmd1, cmd2, cmd3, cmd4]

        run_shell_cmd(cmds, self.stats_file.resolve().as_posix())

    def calculate_non_zero_coverage(self, curr_nz_align: dict):
        """Calculate the coverage for non-zero regions in the mapped reads.

        Args:
            curr_nz_align (dict): Dictionary containing non-zero alignment data.

        Note:
            This method executes a pipeline of commands using bedtools, awk, sed, and calculates the average coverage
            for non-zero regions in the mapped reads.
        """

        # First command: bedtools genomecov
        cmd1 = ["bedtools", "genomecov", "-d", "-ibam", curr_nz_align["mapped_bam"].resolve().as_posix()]

        # Second command: awk to sum coverage
        cmd2 = ["awk", '{arr[$1]+=$3;} END {for (i in arr) print i, arr[i]}']

        # Third command: sed to replace underscores with tabs
        cmd3 = ["sed", "s/_/\t/g"]

        # Fourth command: awk to calculate average for non-zero regions
        cmd4 = ["awk", '{if($5!=0){print $1"_"$2"_"$3"_"$4"\t" $5/($4-$3+1)}}']

        # Combining all commands for the pipeline
        cmds = [cmd1, cmd2, cmd3, cmd4]

        # Call run_shell_cmd
        run_shell_cmd(cmds, self.stats_file.resolve().as_posix())

    def extract_n50_stats(self, curr_nz_align: dict):
        """Extract N50 statistics from the cleaned non-zero length contig FASTA file.

        Args:
            curr_nz_align (dict): Dictionary containing non-zero alignment data.
        """
        # First command: n50 tool
        cmd1 = ["n50", "--format", "tsv", curr_nz_align["cleaned_fasta"].resolve().as_posix()]

        # Second command: sed to pick the second line
        cmd2 = ["sed", "-n", "2 p"]

        # Third command: awk to format output
        cmd3 = ["awk", "{print \"#seqs:\\n\" $2 \"\\nTotal length:\\n\" $3 \"\\nMax:\\n\" $6 }"]

        # Combining all commands for the pipeline
        cmds = [cmd1, cmd2, cmd3]

        # Call run_shell_cmd
        run_shell_cmd(cmds, self.stats_file.resolve().as_posix())

    def calculate_n50_from_fasta(self, curr_nz_align: dict, curr_ref_path: Path):
        """Convert contig and reference FASTA files to tabular format and calculate N50 statistics.

        Args:
            curr_nz_align (dict): Dictionary containing non-zero alignment data.
            curr_ref_path (Path): Path to the current reference genome FASTA file.
        """

        # Convert assembly FASTA to tabular format
        cmd1 = ["seqkit", "fx2tab", "--length", "--name", curr_nz_align["cleaned_fasta"].resolve().as_posix()]
        length_assembly = self.curr_out_dir / "length_assembly.txt"
        run_shell_cmd(cmd1, length_assembly.resolve().as_posix())

        # Convert reference FASTA to tabular format
        cmd2 = ["seqkit", "fx2tab", "--length", "--name", curr_ref_path.resolve().as_posix()]
        length_ref = self.curr_out_dir / "length_ref.txt"
        run_shell_cmd(cmd2, length_ref.resolve().as_posix())

        # Calculate N50 using the Python script

        cmd3 = ["python", self.n50_script.resolve().as_posix(), "-i", length_assembly.resolve().as_posix(), "-r",
                length_ref.resolve().as_posix()]
        run_shell_cmd(cmd3, self.stats_file.resolve().as_posix())

    def analyze(self, seq_list_file: Path, cur_nz_align: dict, curr_ref_path: Path):
        """Perform downstream analysis on aligned reads.

        Args:
            seq_list_file (Path): Path to the list of sequences.
            cur_nz_align (dict): Dictionary containing non-zero alignment data.
            curr_ref_path (Path): Path to the current reference genome FASTA file.
        """
        self.generate_agp_file(seq_list_file)
        self.calculate_average_depth(cur_nz_align)
        self.calculate_non_zero_coverage(cur_nz_align)
        self.extract_n50_stats(cur_nz_align)
        self.calculate_n50_from_fasta(cur_nz_align, curr_ref_path)
