import pytest
import shutil
from pathlib import Path

# Get the root of the project (where both 'scripts' and 'tests' directories are)
project_root = Path(__file__).resolve().parent.parent


@pytest.fixture
def tmp_dir():
    t = Path.cwd() / "tmp"
    t.mkdir(parents=True, exist_ok=True)  # setup
    yield t
    shutil.rmtree(t.resolve().as_posix())  # teardown


@pytest.fixture
def project_root():
    yield Path(__file__).resolve().parent.parent.parent


class SelectionCullingRefDb(object):
    def __init__(self):
        self.dir = Path(__file__).resolve().parent.parent / 'test_data' / 'sel_cull_db'
        self.ref_db = self.dir / "RefSeq_V2_db"
        self.forward_read = self.dir / "reads" / "forward.fq"
        self.reverse_read = self.dir / "reads" / "reverse.fq"


@pytest.fixture
def sel_cull_ref_db():
    yield SelectionCullingRefDb()


@pytest.fixture
def metacompass_runner(project_root):
    r = project_root / "tests" / "utils" / "metacompass_run.sh"
    yield r


@pytest.fixture
def metacompass_path(project_root):
    r = project_root / "metacompass.nf"
    yield r
