import subprocess
from pathlib import Path
from shutil import copy
from shutil import copyfileobj
from typing import List
from Bio import SeqIO
import sys
import traceback


def delete_file(f: [str | Path]):
    if isinstance(f, Path):
        f.unlink()
        return

    remove_string = f"rm -rf {f}"
    try:
        subprocess.Popen(remove_string, shell=True).wait()
    except Exception as e:
        print("Error while removing", str(e))

    if Path(f).exists():
        raise FileExistsError(f"Unable to delete file with path {Path(f).resolve().as_posix()}")


def run_shell_cmd(cmds: List, outfile=None, count_lines=False):
    """
    Run one or more shell commands, optionally piping them together.

    Args
        cmds (List): A single command as a list or multiple commands as a list of lists.
        outfile (Path): Optional output file to redirect the final command's output.
    """

    try:
        if not isinstance(cmds[0], list):  # If cmds is a single command
            cmds = [cmds]

        # Start the first process
        processes = [subprocess.Popen(cmds[0], stdout=subprocess.PIPE)]

        # Start the rest of the processes, connecting each one's input to the previous one's output
        for cmd in cmds[1:]:
            processes.append(subprocess.Popen(cmd, stdin=processes[-1].stdout, stdout=subprocess.PIPE))
            processes[-2].stdout.close()  # Allow previous process to receive a SIGPIPE if the next process exits

        # Wait for the final process to complete and capture its output
        final_process = processes[-1]
        stdout, stderr = final_process.communicate()

        # Check for errors
        if final_process.returncode != 0:
            raise subprocess.CalledProcessError(final_process.returncode, final_process.args)

        # Write the output to a file if given
        if outfile:
            with open(outfile, 'a') as f:
                f.write(stdout.decode())

        if count_lines:
            return stdout.decode(), len(stdout.decode().splitlines())
        else:
            return stdout.decode()

    except subprocess.CalledProcessError as e:
        error_message = e.output.decode() if e.output else "No output"
        print(f"Command {e.cmd} failed with return code {e.returncode} - {error_message}. {e}")
        traceback.print_exc()
    except Exception as e:
        print(f"An error occurred - {str(e)}.")
        traceback.print_exc()


def align_with_minimap2(ref_genome_path: str, interleaved_reads_path: str,
        outfile: Path, threads: str = "1",
                        mapped_only: bool = True):
    """Used to align data using minimap2

    Args:
        ref_genome_path (Path): Path to the reference genome file.
        interleaved_reads_path (Path): Path to the interleaved reads file.
        outfile (Path): Path to the output alignment file.
        threads (str, optional): Number of threads for parallel processing. Default is "1".
        mapped_only (bool, optional): Flag to indicate whether to output only mapped reads. Default is True.

    """
    ref_genome = Path(ref_genome_path)
    ref_genome_size = get_file_size(ref_genome)
    size_limit_flag = ""

    if ref_genome_size > 4:
        size_limit_flag = f"-I{int(ref_genome_size) + 4}g"  # set the size limit 4 more than the size of the file in GB

        cmd = ["minimap2", "-t", str(threads), "--heap-sort=yes",
               size_limit_flag, "-ax", "sr",
               ref_genome_path, interleaved_reads_path]
    else:
        cmd = ["minimap2", "-t", str(threads), "--heap-sort=yes", "-ax", "sr",
               ref_genome_path, interleaved_reads_path]

    if file_empty(Path(ref_genome_path)):
        print(f"File {ref_genome_path} is empty")
        sys.exit(1)

    if file_empty(Path(interleaved_reads_path)):
        print(f"File {interleaved_reads_path} is empty")
        sys.exit(1)

    if mapped_only:
        cmd.insert(3, "--sam-hit-only")

    print("DEBUG: running: " + " ".join(cmd))

    run_shell_cmd(cmd, outfile)


def sam_to_sorted_bam(input_sam: Path, output_bam: Path, threads: str = "1"):
    """Used to generate a sorted bam file from a sam file using samtools

    Args:
        input_sam (Path): Path to the input SAM file.
        output_bam (Path): Path to the output sorted BAM file.
        threads (str, optional): Number of threads for parallel processing. Default is "1".
    """
    view_cmd = ["samtools", "view", "-Sb", "-@", str(threads), input_sam]
    sort_cmd = ["samtools", "sort", "--threads", str(threads), "-o", output_bam, "-"]
    # view_cmd | sort_cmd
    run_shell_cmd([view_cmd, sort_cmd])
    bam_index_cmd = ["samtools", "index", output_bam]
    run_shell_cmd(bam_index_cmd)


def extract_read_ids(sam_file: str, outfile: str, threads: str = "1"):
    """Used to extract read ids from a samfile

    Args:
        sam_file (str): Path to the input SAM file.
        outfile (str): Path to the output file to store extracted read IDs.
        threads (str, optional): Number of threads for parallel processing. Default is "1".
    """
    # Convert SAM to FASTQ
    cmd1 = ["samtools", "fastq", "-@", str(threads), sam_file]

    # Process with awk to extract read IDs
    cmd2 = ["awk", 'NR%4==1 {$0=gensub(/^@(.*)\/.*/, "\\1", 1); print}']

    run_shell_cmd([cmd1, cmd2], outfile)


def extract_reads(mapped_sam: str, mapped_ids: str, interleaved_reads: str, outfile: str, threads="1",
                  mapped: bool = True):
    """Used to extract reads present in a sam file from an interleaved reads file using seqkit & samtools

    Args:
        mapped_sam (str): Path to the mapped SAM file.
        mapped_ids (str): Path to the file containing mapped read IDs.
        interleaved_reads (str): Path to the interleaved reads file.
        outfile (str): Path to the output file to store the extracted reads.
        threads (str, optional): Number of threads for parallel processing. Default is "1".
        mapped (bool, optional): Flag to indicate whether to extract mapped reads. Default is True.

    Returns:
        num_lines (int): total number of lines of extracted reads written to the output file.
    """

    if not Path(mapped_ids).exists():
        extract_read_ids(sam_file=mapped_sam,
                         outfile=mapped_ids,
                         threads=threads)

    cmd = ["seqkit", "grep", "-I", "-j", str(threads), "-f", mapped_ids, interleaved_reads]
    print("DEBUG: running:" + " ".join(cmd))

    if not mapped:
        cmd.insert(5, "-v")
        print("unmapped reads extracted to", outfile)

    _, num_lines = run_shell_cmd(cmd, outfile, count_lines=True)
    return num_lines

    # Path(interleaved_reads).unlink() # delete the interleaved reads


def count_bases(file_path: Path):
    """
    Count the total number of bases in a FASTA or FASTQ file, excluding header lines.

    Args:
        file_path (Path): The path to the input FASTA or FASTQ file.

    Returns:
        int: The total number of bases in the file.
    """
    total_bases = 0
    first_line = None

    with open(file_path.resolve().as_posix(), "r") as file:
        first_line = file.readline().strip()

    if first_line.startswith(">"):
        # Assume it's a FASTA file
        for record in SeqIO.parse(file_path.resolve().as_posix(), "fasta"):
            total_bases += len(record.seq)
    elif first_line.startswith("@"):
        # Assume it's a FASTQ file
        for record in SeqIO.parse(file_path.resolve().as_posix(), "fastq"):
            total_bases += len(record.seq)
    elif first_line.startswith(";"):  # You can add more conditions for other formats
        # Assume it's a FASTA (fna) file with a different header format
        for record in SeqIO.parse(file_path.resolve().as_posix(), "fasta"):
            total_bases += len(record.seq)
    else:
        raise ValueError("Unrecognized file format: not FASTA or FASTQ.")

    return total_bases


def count_reads(file_path: Path):
    """
    Count the total number of reads in a FASTQ file

    Args:
        file_path (Path): The path to the input FASTQ file.

    Returns:
        int: The total number of readss in the file.
    """
    count = 0
    with open(file_path.resolve().as_posix(), "r") as f:
        for _ in SeqIO.parse(f, "fastq"):
            count += 1

    return count


def get_boc(nz_contig_file: Path, ref_genome_file: Path):
    """
   Finds the bread of coverage using the depth file generated by samtools

   Args:
       nz_contig_file (Path): The path to a non-zero contigs FASTA file
       ref_genome_file(Path): The path to the corresponding reference genome file

   Returns:
       float: The breadth of coverage
   """

    nzero_bases = count_bases(nz_contig_file)
    total_bases = count_bases(ref_genome_file)

    print("N-zero bases", nzero_bases)
    print("Total bases", total_bases)

    boc = float(nzero_bases * 100) / total_bases
    print("Calculated breadth of coverage", boc)

    return boc


## Takes a set of file names and concatenates them into a single one
def concatenate_cluster_refs(cluster_list: List, concat_out_file: Path):
    Path(concat_out_file).parent.mkdir(parents=True, exist_ok=True)

    all_refs = []
    concat_out_file_path = concat_out_file.resolve().as_posix()
    for cluster_refs in cluster_list:
        all_refs.extend(cluster_refs)

    total_refs = len(all_refs)

    if total_refs == 0:
        return

    with open(concat_out_file_path, "wb") as outf:
        for i in range(0, total_refs):
            with open(Path(all_refs[i]).resolve().as_posix(), 'rb') as inf:
                copyfileobj(inf, outf)
        outf.close()

# Copies a file
def concat_single_file(src, dst):
    try:
        with open(src, "r") as source_file, open(dst, "a") as destination_file:
            source_content = source_file.read()
            destination_file.write(source_content)
    except FileNotFoundError:
        print(f"Error: Source file '{src}' not found.")
    except PermissionError:
        print(f"Error: Permission denied when accessing files.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


def get_max_len(fasta_file: Path):
    cmd = ["seqkit", "stats", fasta_file.resolve().as_posix()]
    output = run_shell_cmd(cmd)
    max_len = int(output.split('\n')[1].split()[-1].replace(',', ''))
    return max_len


def file_empty(file_path: Path):
    return not file_path.exists() or file_path.stat().st_size == 0


def write_list_to_file(arr, fpath):
    with open(fpath, "w") as f:
        for l in arr:
            f.write(str(l) + "\n")


def get_file_size(file: Path):
    """
    Returns the size of the file in GB
    Args:
       file (Path): Contains the path to the file

    Returns:
       float: size of the file in GB
    """
    if file.exists():
        file_bytes = file.stat().st_size
        file_size_gb = float(file_bytes) / (1024 ** 3)
        return file_size_gb
    else:

        return 0

def move_file(src:Path, dst:Path):
    if src:
        src.rename(dst)
