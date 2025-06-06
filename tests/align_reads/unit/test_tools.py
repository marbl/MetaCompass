import os
import sys
import inspect
from pathlib import Path
from Bio import SeqIO

project_root = Path(__file__).resolve().parent.parent.parent.parent
sys.path.insert(0, project_root.as_posix())
from scripts.align_reads.tools import *


def is_interleaved(fastq_file):
    fastq_iter = SeqIO.parse(fastq_file, "fastq")
    try:
        while True:
            first_read = next(fastq_iter)
            second_read = next(fastq_iter)
            # Check if the read IDs indicate that they are paired
            # This condition may need to be adjusted depending on the naming convention of your reads
            if not (first_read.id[:-2] == second_read.id[:-2] and
                    first_read.id[-1] == "1" and second_read.id[-1] == "2"):
                return False
    except StopIteration:
        # Reached the end of the file without finding unpaired reads
        return True


def test_tools_factory(tmp_dir, project_root):
    tf = ToolsFactory(
        output_dir=tmp_dir,
        threads="1"
    )
    expected_bin_path = project_root / "bin"
    assert tf.bin_path == expected_bin_path
    assert tf.default_output_dir == tmp_dir
    assert tf.default_threads == "1"

    actual_methods = [name for name, _ in inspect.getmembers(tf, callable) if not name.startswith('__')]

    callable_methods = ["get_mash", "get_bbtools", "get_pilon", "get_heap_arranger"]

    assert sorted(callable_methods) == sorted(actual_methods)


def test_mash_tool(tools_factory):
    tf = tools_factory
    mash_obj = tf.get_mash()
    assert isinstance(mash_obj, Mash)


def test_pilon_tool(tools_factory):
    tf = tools_factory
    pilon_obj = tf.get_pilon()
    assert isinstance(pilon_obj, Pilon)


def test_bbtool(tools_factory):
    tf = tools_factory
    bbtools_obj = tf.get_bbtools()
    assert isinstance(bbtools_obj, BBtools)


def test_bbmap_interleave(test_reads_unit, tmp_dir, tools_factory):
    small_fwd = test_reads_unit.small_forward_read
    small_rev = test_reads_unit.small_reverse_read
    small_inter = tmp_dir / "small_interleaved.fq"
    tf = tools_factory
    bbtools_obj = tf.get_bbtools()
    bbtools_obj.interleave_reads(small_fwd, small_rev, small_inter)

    assert is_interleaved(small_inter.resolve().as_posix())


def test_heap_arranger(tools_factory):
    tf = tools_factory
    heap_arranger = tf.get_heap_arranger()
    assert isinstance(heap_arranger, HeapArranger)
