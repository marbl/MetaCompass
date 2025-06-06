"""
pytest -v -s test_refcull.py

bash /mnt/f/ref-culling/testing/align_reads.sh --forward /mnt/f/ref-culling/testing/reads/shakya4k.1.fq --reverse /mnt/f/ref-culling/testing/reads/shakya4k.2.fq -refs /mnt/f/ref-culling/testing/reads/data/                                            -min-list /mnt/f/ref-culling/testing/reads/min_reference_candidates.txt                                            -o /mnt/f/ref-culling/testing/output   

bash /mnt/f/ref-culling/testing/align_reads.sh --forward /mnt/f/ref-culling/testing/reads/test.1.fq --reverse /mnt/f/ref-culling/testing/reads/test.2.fq -refs /mnt/f/ref-culling/testing/reads/testData/                                            -min-list /mnt/f/ref-culling/testing/reads/test.txt                                            -o /mnt/f/ref-culling/testing/output 

"""

import subprocess
import os
import pytest
import sys

cwd = os.getcwd() + "/"
nexistForward = cwd  + "reads/shakya4k.1.f"
nexistReverse = cwd + "reads/shakya4k.2.f"
nexistRefs = cwd + "reads/dat"
nexistMinList = cwd + "reads/min_reference_candi"
nexistOutput = cwd + "outputt"

emptyForward = cwd  + "empty.fq"
emptyReverse = cwd + "empty.fq"
emptyRefs = cwd + "emptyfolder"
emptyMinList = cwd + "empty.txt"

validForward = cwd  + "reads/shakya4k.1.fq"
validReverse = cwd + "reads/shakya4k.2.fq"
validRefs = cwd + "reads/data"
validMinList = cwd + "reads/min_reference_candidates.txt"
validOutput = cwd + "output"
########################################################################################################################


testAlignReadsExitStatus = [
    (f'--help', 0), # test --help
    (f'--forward "{nexistForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{nexistReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{nexistRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{nexistMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{nexistOutput}"', 1),
    (f'--forward "{emptyForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{emptyReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{emptyRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{emptyMinList}" -o "{validOutput}"', 1),
    (f'--reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -min-list "{validMinList}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -o "{validOutput}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}"', 1),
    (f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"', 0) # test valid inputs
]

@pytest.mark.parametrize("params, expected_status", testAlignReadsExitStatus)
def test_return_status(params, expected_status):
    status = subprocess.call(f'./align_reads.sh {params}', shell=True)
    assert expected_status == status

outputMapped1 = validOutput + "/GCF_000016785.1/reads_mapped.1.fq"
outputUnmapped1 = validOutput + "/GCF_000016785.1/reads_unmapped.1.fq"
outputMapped2 = validOutput + "/GCF_000016785.1/reads_mapped.2.fq"
outputUnmapped2 = validOutput + "/GCF_000016785.1/reads_unmapped.2.fq"

testAlignReadsOutputSize = [
    (validForward, outputMapped1, outputUnmapped1),
    (validReverse, outputMapped2, outputUnmapped2)]


@pytest.mark.parametrize("inputF, mappedO, unmappedO", testAlignReadsOutputSize)
def test_output_size(inputF, mappedO, unmappedO):
    if os.path.exists(inputF) and os.path.exists(mappedO) and os.path.exists(unmappedO):
        assert os.stat(inputF).st_size == os.stat(mappedO).st_size + os.stat(unmappedO).st_size
    else:
        assert False

validForward = cwd  + "reads/test.1.fq"
validReverse = cwd + "reads/test.2.fq"
validRefs = cwd + "reads/testData"
validMinList = cwd + "reads/test.txt"
validOutput = cwd + "testOutput"
outputMapped1 = validOutput + "/GCF_000016785.1/reads_mapped.1.fq"
outputUnmapped1 = validOutput + "/GCF_000016785.1/reads_unmapped.1.fq"
outputMapped2 = validOutput + "/GCF_000016785.1/reads_mapped.2.fq"
outputUnmapped2 = validOutput + "/GCF_000016785.1/reads_unmapped.2.fq"

def test_functionality():
    params = f'--forward "{validForward}" --reverse "{validReverse}" -refs "{validRefs}" -min-list "{validMinList}" -o "{validOutput}"'
    status = subprocess.call(f'./align_reads.sh {params}', shell=True)
    assert 0 == status
    assert os.path.exists(outputMapped1) and os.path.exists(outputUnmapped1) and os.path.exists(outputMapped2) and os.path.exists(outputUnmapped2)
    assert os.stat(outputUnmapped1).st_size != 0 and os.stat(outputUnmapped2).st_size != 0
    assert os.stat(outputMapped1).st_size == 0 and os.stat(outputMapped2).st_size == 0



        