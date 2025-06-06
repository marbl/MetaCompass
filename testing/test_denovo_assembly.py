"""
pytest -v -s test_name.py
"""

import subprocess
import os
import pytest
import sys

cwd = os.getcwd() + "/testing/" 
validWorkDir = cwd + "denovo_assembly"

validLog = cwd + 'denovo_assembly/output'
validUnpairedFasta = f'--unpaired {cwd}denovo_assembly/input/unpaired/mc.sam.unmapped.u.fq'
validUnpairedOutput = cwd + 'denovo_assembly/output/unpaired'

validPairedFasta = f'--forward {cwd}denovo_assembly/input/paired/mc.sam.unmapped.1.fq --reverse {cwd}denovo_assembly/input/paired/mc.sam.unmapped.2.fq'
validPairedOutput = cwd + 'denovo_assembly/output/paired'

validPairedUnpairedFasta = f'--forward {cwd}denovo_assembly/input/paired_unpaired/mc.sam.unmapped.1.fq --reverse {cwd}denovo_assembly/input/paired_unpaired/mc.sam.unmapped.2.fq --unpaired {cwd}denovo_assembly/input/paired_unpaired/mc.sam.unmapped.u.fq'
validPairedUnpairedOutput = cwd + 'denovo_assembly/output/paired_unpaired'

invalidOnlyForwardFasta = f'--forward {cwd}denovo_assembly/input/paired/mc.sam.unmapped.1.fq'
invalidOnlyReverseFasta = f'--reverse {cwd}denovo_assembly/input/paired/mc.sam.unmapped.2.fq'
invalidEmptyFasta = '""'
invalidUnpairedOutput = cwd + "blabla"
invalidLog = "blabla"
invalidWorkDir = "blabla"

pairedFinalContigs = cwd + "denovo_assembly/output/paired/denovo_assembly/megahit/final.contigs.fa"
expectedPairedFinalContigs = cwd + "denovo_assembly/expected/final.contigs_paired.fa"

pairedUnpairedFinalContigs = cwd + "denovo_assembly/output/paired_unpaired/denovo_assembly/megahit/final.contigs.fa"
expectedPairedUnpairedFinalContigs = cwd + "denovo_assembly/expected/final.contigs_paired_unpaired.fa"

"""
unpairedFinalContigs = cwd + "valid/unpaired/denovo_assembly/megahit/final.contigs.fa"
expectedUnpairedFinalContigs = cwd + "expected/final.contigs_pairedunpaired.fa"
"""

########################################################################################################################


testDeNovoExitStatus = [
    (f'--reads "{validUnpairedFasta}" --output {invalidUnpairedOutput} --log {validLog}', 1), # not exist output
    (f'--reads "{validUnpairedFasta}" --output {validUnpairedOutput} --log {invalidLog}', 1), # not exist log
    (f'--reads "{validUnpairedFasta}" --output {validUnpairedOutput} --log {validLog} --workdir {invalidWorkDir}', 1), # not exist work dir
    

    (f'--output "{validUnpairedOutput}" --log {validLog}', 1), # no reads
    (f'--reads "{validUnpairedFasta}" --log {validLog}', 1), # no output

    (f'--reads "{invalidOnlyForwardFasta}" --output {validPairedOutput} --log {validLog}', 1), # only forward
    (f'--reads "{invalidOnlyReverseFasta}" --output {validPairedOutput} --log {validLog}', 1), # only reverse
    (f'--reads "{invalidOnlyForwardFasta} {validUnpairedFasta}" --output {validPairedOutput} --log {validLog}', 1), # only forward and unpaired
    (f'--reads "{invalidOnlyReverseFasta} {validUnpairedFasta}" --output {validPairedOutput} --log {validLog}', 1), # only reverse and unpaired
    (f'--reads "{invalidEmptyFasta}" --output {validPairedOutput} --log {validLog}', 1), # empty for reads
    
    (f'--help', 0),
    (f'--reads "{validUnpairedFasta}" --output {validUnpairedOutput} --log {validLog} --workdir {validWorkDir} --threads 1', 0), # valid unpaired reads
    (f'--reads "{validPairedFasta}" --output {validPairedOutput} --log {validLog} --workdir {validWorkDir} --threads 1', 0), # valid paired reads
    
    (f'--reads "{validPairedUnpairedFasta}" --output {validPairedUnpairedOutput} --log {validLog} --workdir {validWorkDir} --threads 1', 0) # valid unpaired and paired reads
]



@pytest.mark.parametrize("params, expected_status", testDeNovoExitStatus)
def test_return_status(params, expected_status):
    status = subprocess.call(f'nextflow run pipeline/denovo_assembly.nf {params}', shell=True)
    assert expected_status == status


########################################################################################################################


def identical(f1, f2):
  file1 = open(f1, "r")
  lines1 = file1.readlines()

  file2 = open(f2, "r")
  lines2 = file2.readlines()

  seqs = []
  seqs2 = []
  for line in lines1:
    if not line.startswith(">"):
      seqs.append(line[:-1])

  for line in lines2:
    if not line.startswith(">"):
      tmp = line[:-1]
      if tmp in seqs:
        seqs.remove(tmp) 
      else:
        seqs2.append(tmp)

  if len(seqs) == 0 and len(seqs2) == 0:
    return True
  else:
    return False

testFinalContigs = [
	(pairedFinalContigs, expectedPairedFinalContigs),
	(pairedUnpairedFinalContigs, expectedPairedUnpairedFinalContigs)
]

@pytest.mark.parametrize("contigs, expectedContigs", testFinalContigs)
def test_final_contigs(contigs, expectedContigs):
	assert identical(contigs, expectedContigs) is True
