import os
import sys
import subprocess
import pysam
from filecmp import cmp
from unittest.mock import patch, MagicMock, call
from pathlib import Path

# from ....scripts.align_reads.utils import *
project_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, project_root.as_posix())
from scripts.align_reads.utils import *

project_name = Path(__file__).resolve().parent.parent.parent.parent.name

subprocess_patch = f"{project_name}.scripts.align_reads.utils.subprocess.Popen"


def count_sam_reads(sam_file: Path):
    read_count = 0
    with pysam.AlignmentFile(sam_file, "r") as samfile:
        for _ in samfile:
            read_count += 1

    return read_count


def is_bam_sorted(bam_file: Path):
    try:
        with pysam.AlignmentFile(bam_file.resolve().as_posix(), 'rb') as bf:
            return bf.header.get('HD', {}).get('SO') == 'coordinate'
    except Exception as e:
        print(f"Error: {e}")
        return False


def test_run_shell_cmd_single_command():
    with patch(subprocess_patch) as mock_popen:
        # Mocking the subprocess to return a simple output
        process_mock = mock_popen.return_value
        process_mock.communicate.return_value = (b"test output", b"")
        process_mock.returncode = 0  # Set the return code to 0

        run_shell_cmd(["echo", "hello"])

        # Check if the command was called correctly
        mock_popen.assert_called_once_with(["echo", "hello"], stdout=subprocess.PIPE)


def test_run_shell_cmd_with_mocked_piping():
    with patch(subprocess_patch) as mock_popen:
        # Mock the communicate method to return output for echo and grep respectively
        mock_process1 = MagicMock()
        mock_process1.communicate.return_value = (b"hello world", b"")
        mock_process1.returncode = 0

        mock_process2 = MagicMock()
        mock_process2.communicate.return_value = (b"hello", b"")
        mock_process2.returncode = 0

        # Make the mock_popen return our two mock processes in order
        mock_popen.side_effect = [mock_process1, mock_process2]

        # Call the function with piping scenario
        result = run_shell_cmd([["echo", "hello world"], ["grep", "hello"]])

        # Check the commands were called correctly
        mock_popen.assert_has_calls([
            call(["echo", "hello world"], stdout=subprocess.PIPE),
            call(["grep", "hello"], stdin=mock_process1.stdout, stdout=subprocess.PIPE)
        ])

        # Check the final output
        # assert result == "hello"


def test_run_shell_cmd_with_real_piping():
    # Use the function to run a real echo and grep command and write to an outfile
    outfile_name = "test_output.txt"
    run_shell_cmd([["echo", "hello world"], ["grep", "hello"]], outfile=outfile_name)

    # Read the contents of the output file
    with open(outfile_name, 'r') as f:
        result = f.read()

    # Cleanup: Optionally delete the outfile after reading its contents
    os.remove(outfile_name)

    # Check that the output is as expected
    assert result == "hello world\n"


def test_align_with_minimap2_all(test_reads_unit, tmp_dir):
    refs_path = test_reads_unit.ref_fna.resolve().as_posix()
    interleaved_path = test_reads_unit.interleaved_reads.resolve().as_posix()

    out_sam = tmp_dir / "all_test.sam"

    # imported from utils
    align_with_minimap2(
        ref_genome_path=refs_path,
        interleaved_reads_path=interleaved_path,
        outfile=out_sam.resolve().as_posix(),
        threads=1,
        mapped_only=False
    )

    assert out_sam.exists()
    out_sam_reads = count_sam_reads(out_sam)
    assert out_sam_reads == 65925
    out_sam.unlink()


def test_align_with_minimap2_mapped(test_reads_unit, tmp_dir):
    refs_path = test_reads_unit.ref_fna.resolve().as_posix()
    interleaved_path = test_reads_unit.interleaved_reads.resolve().as_posix()

    out_sam = tmp_dir / "mapped_test.sam"

    # imported from utils
    align_with_minimap2(
        ref_genome_path=refs_path,
        interleaved_reads_path=interleaved_path,
        outfile=out_sam.resolve().as_posix(),
        threads=1,
        mapped_only=True
    )

    assert out_sam.exists()
    out_sam_reads = count_sam_reads(out_sam)
    assert out_sam_reads == 64845
    out_sam.unlink()


#
#
def test_sam_to_sorted_bam(test_reads_unit, tmp_dir):
    input_sam = test_reads_unit.mapped_sam
    output_bam = tmp_dir / "sorted.bam"
    threads = "1"

    # imported from utils
    sam_to_sorted_bam(input_sam, output_bam.resolve().as_posix(), threads)
    assert is_bam_sorted(output_bam)


#
#
def test_extract_reads_id(test_reads_unit, tmp_dir):
    out_reads_id = tmp_dir / "read_ids.txt"
    input_sam = test_reads_unit.mapped_sam
    mapped_reads_ids = test_reads_unit.mapped_reads_ids

    # imported from utils
    extract_read_ids(input_sam.resolve().as_posix(),
                     out_reads_id.resolve().as_posix(), threads=1)

    assert cmp(mapped_reads_ids, out_reads_id, shallow=False)
    out_reads_id.unlink()


#
#
def test_extract_reads_mapped(test_reads_unit, tmp_dir):
    out_reads_fq = tmp_dir / "mapped_test.fq"
    input_sam = test_reads_unit.mapped_sam
    read_ids = test_reads_unit.mapped_reads_ids
    interleaved_reads = test_reads_unit.interleaved_reads
    expected_mapped_fq = test_reads_unit.mapped_fq

    # imported from utils
    num_lines = extract_reads(
        input_sam.resolve().as_posix(),
        read_ids.resolve().as_posix(),
        interleaved_reads.resolve().as_posix(),
        out_reads_fq.resolve().as_posix(),
        threads=1, mapped=True)

    assert cmp(expected_mapped_fq, out_reads_fq, shallow=False)
    print("Number of extracted reads are", num_lines)

    out_reads_fq.unlink()


def test_extract_reads_unmapped(test_reads_unit, tmp_dir):
    out_reads_fq = tmp_dir / "unmapped.fq"
    input_sam = test_reads_unit.all_sam
    read_ids = test_reads_unit.all_reads_ids
    interleaved_reads = test_reads_unit.interleaved_reads
    expected_fq = test_reads_unit.unmapped_fq

    # imported from utils
    extract_reads(
        input_sam.resolve().as_posix(),
        read_ids.resolve().as_posix(),
        interleaved_reads.resolve().as_posix(),
        out_reads_fq.resolve().as_posix(),
        threads=1, mapped=False)

    assert cmp(expected_fq, out_reads_fq, shallow=False)

    out_reads_fq.unlink()
