import pytest
import sys
import os
import shutil
import pysam

from pathlib import Path

# Get the root of the project (where both 'scripts' and 'tests' directories are)
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, project_root.as_posix())

# Now you should be able to import your modules
from scripts.align_reads.read_aligner import ReadAligner
from scripts.align_reads.tools import ToolsFactory


# from ...scripts.align_reads.read_aligner import ReadAligner
# from ...scripts.align_reads.tools import ToolsFactory

# Common fixtures
@pytest.fixture
def tmp_dir():
    t = Path.cwd() / "tmp"
    t.mkdir(parents=True, exist_ok=True)  # setup
    yield t
    shutil.rmtree(t.resolve().as_posix())  # teardown


@pytest.fixture
def project_root():
    yield Path(__file__).resolve().parent.parent.parent


# Fixtures for unit tests
class TestReadsUnit(object):
    def __init__(self):
        self.dir = Path(__file__).resolve().parent.parent / 'test_data' / 'align_reads' / 'unit'
        self.forward_read = self.dir / "forward.fq"
        self.reverse_read = self.dir / "reverse.fq"
        self.small_forward_read = self.dir / "small_forward.fq"
        self.small_reverse_read = self.dir / "small_reverse.fq"
        self.interleaved_reads = self.dir / "interleaved.fq"
        self.mapped_sam = self.dir / "mapped.sam"
        self.all_sam = self.dir / "all.sam"
        self.ref_fna = self.dir / "refs.fna"
        self.mapped_reads_ids = self.dir / "mapped_reads_ids.txt"
        self.all_reads_ids = self.dir / "all_reads_ids.txt"
        self.mapped_fq = self.dir / "mapped.fq"
        self.all_fq = self.dir / "all.fq"
        self.unmapped_fq = self.dir / "unmapped.fq"


@pytest.fixture
def test_reads_unit():
    yield TestReadsUnit()


@pytest.fixture
def tools_factory(tmp_dir):
    tf = ToolsFactory(
        output_dir=tmp_dir,
        threads="1"
    )
    yield tf


#
# @pytest.fixture
# def tool_paths():
#     curr_file = Path(__file__).resolve()
#     project_src = curr_file.parent.parent.parent
#
#     paths = {
#         "bin_path": project_src / "bin",
#         "scripts_path": project_src / "scripts"
#     }
#
#     yield paths


# Fixtures for integration tests
class TestReadsIntegration(object):
    def __init__(self):
        self.dir = Path(__file__).resolve().parent.parent / 'test_data' / 'align_reads' / 'integration'
        self.forward_read = self.dir / "forward.fq"
        self.reverse_read = self.dir / "reverse.fq"
        self.interleaved_read = self.dir / "interleaved.fq"
        self.actual_culling = self.dir / "actual_reference_culling"
        self.actual_selection = self.dir / "actual_reference_selection"
        self.expected_assembly = self.dir / "expected_reference_assembly"
        self.ref_db = Path(__file__).resolve().parent.parent / 'test_data' / "ref_db"

        ref_fna_path = "GCA_009556455.1/GCA_009556455.1_ASM955645v1_genomic.fna"

        # create clusters.txt in actual_culling
        with open(self.actual_culling / "clusters.txt", "w") as f:
            txt_to_write = f"{self.ref_db.resolve().as_posix()}/{ref_fna_path}"
            f.write(txt_to_write)

        # create cluster_content.csv
        with open(self.actual_culling / "cluster_content.csv", "w") as f:
            txt_to_write = f"cluster_1,{self.ref_db.resolve().as_posix()}/{ref_fna_path}"
            f.write(txt_to_write)


@pytest.fixture
def test_reads_integration():
    yield TestReadsIntegration()


@pytest.fixture
def inputs_integration(test_reads_integration, tmp_dir):
    # copy expected culling and selection into temporary directory
    tmp_culling = tmp_dir / "actual_reference_culling"
    tmp_selection = tmp_dir / "actual_reference_selection"
    tmp_assembly = tmp_dir / "actual_reference_assembly"

    shutil.copytree(test_reads_integration.actual_culling, tmp_culling)
    shutil.copytree(test_reads_integration.actual_selection, tmp_selection)

    t = {
        "interleaved_reads": tmp_culling / "interleaved_reads.fq",
        "cluster_list": tmp_culling / "clusters.txt",
        "ref_genome_sketch": tmp_culling,
        "all_reads_sketch": tmp_culling / "reads.base_2.kmers",
        "out": tmp_assembly,
        "max_contig_length": 2000,
        "breadth_of_coverage": 5.0,
        "threads": 1,
        "debug": True,
    }

    yield t


def initialize_aligner(inputs):
    aligner = ReadAligner()
    aligner.inputs = inputs
    aligner.init_dependencies()
    aligner.set_interleaved_reads()

    aligner.cluster_mapped_reads = Path(aligner.inputs["out"]) / "cluster_mapped.fq"
    aligner.mapped_clusters = Path(aligner.inputs["out"]) / "mapped_clusters.txt"

    return aligner


@pytest.fixture
def read_aligner_integration(inputs_integration):
    aligner = initialize_aligner(inputs_integration)
    yield aligner
