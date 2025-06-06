# """
# pytest -v -s test_name.py
# """

#Scripts: collect_kmers.sh

import subprocess
import os
import shutil
import pytest
import sys
from pathlib import Path

cwd = os.getcwd() # testing dir
workdir = Path(cwd) # metacompass workdir is parent dir of testing dir

# READS: tutorial sample
valid_forward = f'{workdir}/tutorial/thao2000.1.fq'
valid_reverse = f'{workdir}/tutorial/thao2000.2.fq'
valid_unpaired = f'{workdir}/tutorial/c_rudii.fastq'
invalid_read = f'{workdir}/tutorial/file_not_exist'

# REFERENCES:
valid_1ref_folder = f'{workdir}/testing/collect_kmers_data/1ref'
valid_5ref_folder = f'{workdir}/testing/collect_kmers_data/5ref'
invalid_0ref_folder = f'{workdir}/testing/collect_kmers_data/0ref'

# Create output dir for pytest
test_output = f'{workdir}/testing/collect_kmers_output'

if os.path.exists(test_output):
    print(f"Output folder {test_output} exists and would be overwritten.") 
    shutil.rmtree(test_output)
os.makedirs(test_output)

# Create 1 output folder for all invalid cases. This folder includes log file of all invalid cases
invalid_cases = f'{test_output}/invalid'
os.mkdir(invalid_cases)

########################## PARAMS ########################
# ----- 0. PARAM: HELP -----
params_help = f'--help'

# ----- 1. PROVIDE VALID READ AND VALID REF FOLDER -----
# Reads: valid paired + unpaired, ref folder : has 1 ref
out_paired_unpaired_1ref = f'{test_output}/paired_unpaired_1ref' #output folder
os.mkdir(out_paired_unpaired_1ref) #create output folder
params_paired_unpaired_1ref = f'-f {valid_forward} -r {valid_reverse} -u {valid_unpaired}\
                                -refs {valid_1ref_folder} -o {out_paired_unpaired_1ref}\
                                -l {out_paired_unpaired_1ref}/log'

# Reads: valid paired + unpaired, ref folder : has 5 refs
out_paired_unpaired_5ref = f'{test_output}/paired_unpaired_5ref'
os.mkdir(out_paired_unpaired_5ref)
params_paired_unpaired_5ref = f'-f {valid_forward} -r {valid_reverse} -u {valid_unpaired}\
                                -refs {valid_5ref_folder} -o {out_paired_unpaired_5ref}\
                                -l {out_paired_unpaired_5ref}/log'



# Reads: valid paired, ref folder: has 1 ref
out_paired_1ref = f'{test_output}/paired_1ref'
os.mkdir(out_paired_1ref)
params_paired_1ref = f'-f {valid_forward} -r {valid_reverse}\
                                -refs {valid_1ref_folder} -o {out_paired_1ref}\
                                -l {out_paired_1ref}/log'

# Reads: valid paired, ref folder: has 5 refs
out_paired_5ref = f'{test_output}/paired_5ref'
os.mkdir(out_paired_5ref)
params_paired_5ref = f'-f {valid_forward} -r {valid_reverse}\
                        -refs {valid_5ref_folder} -o {out_paired_5ref}\
                        -l {out_paired_5ref}/log'

# Reads: valid unpaired, ref folder: has 1 ref
out_unpaired_1ref = f'{test_output}/unpaired_1ref'
os.mkdir(out_unpaired_1ref)
params_unpaired_1ref = f'-u {valid_unpaired}\
                        -refs {valid_1ref_folder} -o {out_unpaired_1ref}\
                        -l {out_unpaired_1ref}/log'

# Reads: valid unpaired, ref folder: has 5 refs
out_unpaired_5ref = f'{test_output}/unpaired_5ref'
os.mkdir(out_unpaired_5ref)
params_unpaired_5ref = f'-u {valid_unpaired}\
                        -refs {valid_5ref_folder} -o {out_unpaired_5ref}\
                        -l {out_unpaired_5ref}/log'

## ----- 2. PROVIDE READ AND REF BUT 1 OF THEM NOT VALID -----
# Reads: valid paired + unpaired, invalid ref folder : ref folder has 0 ref
params_paired_unpaired_0ref = f'-f {valid_forward} -r {valid_reverse} -u {valid_unpaired}\
                                -refs {invalid_0ref_folder} -o {invalid_cases}\
                                -l {invalid_cases}/paired_unpaired_0ref.log'

# Reads: valid paired, invalid ref folder : ref folder has 0 ref
params_paired_0ref = f'-f {valid_forward} -r {valid_reverse} \
                        -refs {invalid_0ref_folder} -o {invalid_cases}\
                        -l {invalid_cases}/paired_0ref.log'

# Reads: valid unpaired, invalid ref folder : ref folder has 0 ref
params_unpaired_0ref = f'-u {valid_unpaired}\
                        -refs {invalid_0ref_folder} -o {invalid_cases}\
                        -l {invalid_cases}/unpaired_0ref.log'


# Reads: valid paired, invalid ref folder : folder provided but not exist
params_invalid_refs = f'-f {valid_forward} -r {valid_reverse}\
                        -refs /folder/not/exist -o {invalid_cases}\
                        -l {invalid_cases}/invalid_refs.log'

# Reads: provide valid paired and invalid unpaired, ref folder: has 5 refs
params_valid_paired_invalid_unpaired = f'-f {valid_forward} -r {valid_reverse} -u {invalid_read}\
                                    -refs {valid_5ref_folder} -o {invalid_cases}\
                                    -l {invalid_cases}/valid_paired_invalid_unpaired.log'

# Reads: provide invalid paired and invalid unpaired, ref folder: has 5 refs
params_invalid_paired_invalid_unpaired = f'-f {invalid_read} -r {invalid_read} -u {invalid_read}\
                                    -refs {valid_5ref_folder} -o {invalid_cases}\
                                    -l {invalid_cases}/invalid_paired_invalid_unpaired.log'


## ----- 3. NOT PROVIDE REF OR READ -----
# Reads: not provided, ref folder: has 5 refs
params_read_not_provided = f'-refs {valid_5ref_folder} -o {invalid_cases}\
                            -l {invalid_cases}/not_provided_read.log'

# Reads: valid paired, ref folder: not provided
params_refs_not_provided = f'-f {valid_forward} -r {valid_reverse} \
                            -o {invalid_cases} -l {invalid_cases}/not_provided_refs.log'

### ----- 4. PROBLEMS WITH OUTPUT FOLDER -----
# Output folder is provided but not exist
params_invalid_outfolder = f'-u {valid_unpaired}\
                            -refs {valid_5ref_folder} -o folder/not/exist\
                            -l {invalid_cases}/invalid_outfolder.log'


# Output folder is not provided
params_out_folder_not_provided = f'-u {valid_unpaired}\
                                -refs {valid_5ref_folder} \
                                -l {invalid_cases}/not_provided_outfolder.log'

########################## PYTEST ########################
testCollectKmers = [
    (params_help, 0),
    (params_paired_unpaired_1ref, 0),
    (params_paired_unpaired_5ref, 0),
    (params_paired_1ref, 0),
    (params_paired_5ref, 0),
    (params_unpaired_1ref, 0),
    (params_unpaired_5ref, 0),
    (params_paired_unpaired_0ref, 1),
    (params_paired_0ref, 1),
    (params_unpaired_0ref, 1),
    (params_invalid_refs, 1),
    (params_valid_paired_invalid_unpaired, 1),
    (params_invalid_paired_invalid_unpaired, 1),
    (params_read_not_provided, 1),
    (params_refs_not_provided, 1),
    (params_invalid_outfolder, 1),
    (params_out_folder_not_provided, 1)
]

@pytest.mark.parametrize("params, expected_status", testCollectKmers)
def test_return_status_refsel(params, expected_status):
    status = subprocess.call(f'{workdir}/scripts/collect_kmers.sh -ms 28 {params}', shell=True)
    assert expected_status == status
    