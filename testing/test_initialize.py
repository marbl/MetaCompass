"""
In root directory of project:

pytest -v -s testing/initialize_test.py
"""

import subprocess
import os
import pytest

# Script: check_exist.sh
validPrograms = "nextflow,samtools,bedtools,python3,bowtie2,kmer-mask,datasets"
invalidPrograms = "program1,programZ"

validFiles = "testing/initialize/valid_file.fasta"
invalidFiles = "testing/initialize/invalid_file.fasta"

# Script: initialize.nf
validReference = "testing/initialize/valid_reference"
invalidReference = "testing/initialize/invalid_reference"

validOut = "testing/initialize/valid_output"
invalidOut = "testing/initialize/invalid_output"


########################################################################################################################

############# Script: check_exist.sh #############

testCheckExist = [
    (f"--help", 0),
    (f"-b {validPrograms}", 0),
    (f"-f {validFiles}", 0),
    (f"-b {invalidPrograms}", 1),
    (f"-f {invalidFiles}", 1),
]

@pytest.mark.parametrize("params, expected_status", testCheckExist)
def test_return_status(params, expected_status):
    status = subprocess.call(
        f"bash scripts/check_exist.sh {params}", shell=True)
    assert expected_status == status

############# Script: initialize.nf #############

testInitialize = [
    (f"--reference_db {validReference} --output {validOut}", 0),
    (f"--help", 0),
    (f"--output {validOut}", 0), # case when user only wants to run denovo assembly
    (f"--output {invalidOut}", 1),
    (f"--reference_db {validReference}", 1),
    (f"--reference_db {invalidReference}", 1),
    (f"--reference_db {invalidReference} --output {validOut}", 1),
    (f"--reference_db {validReference} --output {invalidOut}", 1),
]

print("\nDebug: validOut = " + validOut + "\n")

@pytest.mark.parametrize("params, expected_status", testInitialize)
def test_initialize_return_status(params, expected_status):
    status = subprocess.call(
        f"nextflow pipeline/initialize.nf {params}", shell=True)
    assert expected_status == status
