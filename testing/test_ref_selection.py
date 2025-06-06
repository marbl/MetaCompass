# """
# pytest -v -s test_name.py
# """

import subprocess
import os
import shutil
import pytest
import sys
from pathlib import Path

#Valid vars
cwd = os.getcwd() # testing dir
#workdir = Path(cwd).parent # metacompass workdir is parent dir of testing dir
workdir = cwd
print(f"workdir = {workdir}")
valid_DB = '/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/RefSeq_V1_db' 
valid_forward = f'{workdir}/tutorial/thao2000.1.fq' #tutorial sample
valid_reverse = f'{workdir}/tutorial/thao2000.2.fq'
valid_unpaired = f'{workdir}/tutorial/c_rudii.fastq'
# valid_forward = '/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/hmp_stool/SRS011302.denovo_duplicates_marked.trimmed.1.fastq'
# valid_reverse = '/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/hmp_stool/SRS011302.denovo_duplicates_marked.trimmed.2.fastq'
# valid_unpaired = '/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/hmp_stool/SRS011302.denovo_duplicates_marked.trimmed.singleton.fastq'
test_output = f'{workdir}/testing/ref_selection_output'

if os.path.exists(test_output):
    print(f"Output folder {test_output} exists and would be overwritten.") 
    shutil.rmtree(test_output)
os.makedirs(test_output)

#Invalid vars
invalid_read = f"{workdir}/testing/ref_selection_data/file_not_exist.txt"
invalid_DB = "/fs/cbcb-lab/mpop/Fraunhofer_metagenomics/folder_not_exist"

#Valid Parameters
params_valid_paired = f'--forward {valid_forward} --reverse {valid_reverse}'
params_valid_unpaired = f'--unpaired {valid_unpaired}'
params_valid_paired_unpaired = f'--forward {valid_forward} --reverse {valid_reverse} --unpaired {valid_unpaired}'
params_valid_DB = f'--reference {valid_DB}'
params_filter =  '--ms 28 \
                --clean 0.0 \
                --match 0.01 \
                --readlen 200 '
params_cluster = '--deapth_of_coverage 5 \
                --breadth_of_coverage 0.9 \
                --percent_markers_covered 75 \
                --threads 12 '
params_help = '--help'


################################################################################
#Pipeline: ref_selection.nf

#refsel_paired: Valid, both filtering and clustering, paired reads
valid_out_1 = f'{test_output}/refsel_paired'
valid_log_1 = f'{test_output}/refsel_paired/refsel_paired.log'
os.mkdir(valid_out_1)
params_refsel_paired = f'{params_valid_paired} {params_valid_DB} \
                        --filter_refs true {params_filter} {params_cluster} \
                        --output {valid_out_1} --workdir {workdir} \
                        &>> {valid_log_1}'

#refsel_unpaired: Valid, both filtering and clustering, unpaired reads
valid_out_2 = f'{test_output}/refsel_unpaired'
valid_log_2 = f'{test_output}/refsel_unpaired/refsel_unpaired.log'
os.mkdir(valid_out_2)
params_refsel_unpaired = f'{params_valid_unpaired} {params_valid_DB} \
                        --filter_refs true {params_filter} {params_cluster} \
                        --output {valid_out_2} --workdir {workdir} \
                        &>> {valid_log_2}'

#refsel_paired_unpaired: Valid, both filtering and clustering, both paired and unpaired reads
valid_out_3 = f'{test_output}/refsel_paired_unpaired'
valid_log_3 = f'{test_output}/refsel_paired_unpaired/refsel_paired_unpaired.log'
os.mkdir(valid_out_3)
params_refsel_paired_unpaired = f'{params_valid_paired_unpaired} {params_valid_DB} \
                                --filter_refs true {params_filter} {params_cluster} \
                                --output {valid_out_3} --workdir {workdir} \
                                &>> {valid_log_3}'

#refsel_paired_nofilter: Valid, only clustering, paired reads
valid_out_4 = f'{test_output}/refsel_paired_nofilter'
valid_log_4 = f'{test_output}/refsel_paired_nofilter/refsel_paired_nofilter.log'
os.mkdir(valid_out_4)
params_refsel_paired_nofilter = f'{params_valid_paired} {params_valid_DB} \
                                --filter_refs false {params_filter} {params_cluster} \
                                --output {valid_out_4} --workdir {workdir} \
                                &>> {valid_log_4}'

#refsel_unpaired_nofilter: Valid, only clustering, unpaired reads
valid_out_5 = f'{test_output}/refsel_unpaired_nofilter'
valid_log_5 = f'{test_output}/refsel_unpaired_nofilter/refsel_unpaired_nofilter.log'
os.mkdir(valid_out_5)
params_refsel_unpaired_nofilter = f'{params_valid_unpaired} {params_valid_DB} \
                                --filter_refs false {params_filter} {params_cluster} \
                                --output {valid_out_5} --workdir {workdir} \
                                &>> {valid_log_5}'

#refsel_paired_unpaired_nofilter: Valid, only clustering, both paired and unpaired reads
valid_out_6 = f'{test_output}/refsel_paired_unpaired_nofilter'
valid_log_6 = f'{test_output}/refsel_paired_unpaired_nofilter/refsel_paired_unpaired_nofilter.log'
os.mkdir(valid_out_6)
params_refsel_paired_unpaired_nofilter = f'{params_valid_paired_unpaired} \
                                {params_valid_DB} \
                                --filter_refs false {params_filter} {params_cluster} \
                                --output {valid_out_6} --workdir {workdir} \
                                &>> {valid_log_6}'

#refsel_invalid_read: Invalid path to read file, valid output and log file
valid_out_7 = f'{test_output}/refsel_invalid_read'
valid_log_7 = f'{test_output}/refsel_invalid_read/refsel_invalid_read.log'
os.mkdir(valid_out_7)
params_refsel_invalid_read = f'--forward {invalid_read} --reverse {invalid_read} \
                            {params_valid_DB} \
                            --filter_refs true {params_filter} {params_cluster} \
                            --output {valid_out_7} --workdir {workdir} \
                            &>> {valid_log_7}'

#refsel_invalid_DB: Invalid Database, valid output and log file
valid_out_8 = f'{test_output}/refsel_invalid_DB'
valid_log_8 = f'{test_output}/refsel_invalid_DB/refsel_invalid_DB.log'
os.mkdir(valid_out_8)
params_refsel_invalid_DB = f'{params_valid_paired_unpaired} --reference {invalid_DB} \
                            --filter_refs true {params_filter} {params_cluster} \
                            --output {valid_out_8} --workdir {workdir} \
                            &>> {valid_log_8}'

#refsel_invalid_workdir: Workdir is incorrect
valid_out_9 = f'{test_output}/refsel_invalid_workdir'
valid_log_9 = f'{test_output}/refsel_invalid_workdir/refsel_invalid_workdir.log'
os.mkdir(valid_out_9)
params_refsel_invalid_workdir = f'{params_valid_paired_unpaired} {params_valid_DB} \
                            --filter_refs true {params_filter} {params_cluster} \
                            --output {valid_out_9} --workdir /folder/not/exist \
                            &>> {valid_log_9}'

#refsel_invalid_outfolder: Valid read & Database, invalid output folder
valid_log_10 = f'{test_output}/refsel_invalid_outfolder.log'
params_refsel_invalid_outfolder = f'{params_valid_paired_unpaired} {params_valid_DB} \
                --filter_refs true {params_filter} {params_cluster} \
                --output /folder/not/exist --workdir {workdir} \
                &>> {valid_log_10}'

testRefSelection = [
    (params_help, 0),
    (params_refsel_paired_unpaired, 0),
    (params_refsel_paired, 0),
    (params_refsel_unpaired, 0),
    (params_refsel_paired_unpaired_nofilter, 0),
    (params_refsel_paired_nofilter, 0),
    (params_refsel_unpaired_nofilter, 0),
    (params_refsel_invalid_read, 1),
    (params_refsel_invalid_DB, 1),
    (params_refsel_invalid_workdir, 1),
    (params_refsel_invalid_outfolder, 1)
]

@pytest.mark.parametrize("params, expected_status", testRefSelection)
def test_return_status_refsel(params, expected_status):
    print(f'nextflow run {workdir}/pipeline/ref_selection.nf {params}')
    status = subprocess.call(f'nextflow run {workdir}/pipeline/ref_selection.nf {params}', shell=True)
    assert expected_status == status

################################################################################
#Script: filter_refs.sh

#Both paired and unpaired reads
valid_out_11 = f'{test_output}/filter_paired_unpaired'
valid_log_11 = f'{test_output}/filter_paired_unpaired/filter_paired_unpaired.log'
os.mkdir(valid_out_11)
params_filter_paired_unpaired = f'{params_valid_paired_unpaired} {params_valid_DB} \
                                {params_filter} --output {valid_out_11} \
                                --log {valid_log_11} \
                                &>> {valid_log_11}'

#Only paired reads
valid_out_12 = f'{test_output}/filter_paired'
valid_log_12 = f'{test_output}/filter_paired/filter_paired.log'
os.mkdir(valid_out_12)
params_filter_paired = f'{params_valid_paired} {params_valid_DB} \
                    {params_filter} --output {valid_out_12} \
                    --log {valid_log_12} \
                    &>> {valid_log_12}'

#Only unpaired read
valid_out_13 = f'{test_output}/filter_unpaired'
valid_log_13 = f'{test_output}/filter_unpaired/filter_unpaired.log'
os.mkdir(valid_out_13)
params_filter_unpaired = f'{params_valid_unpaired} {params_valid_DB} \
                        {params_filter} --output {valid_out_13} \
                        --log {valid_log_13} \
                        &>> {valid_log_13}'

#filter_invalid_read
valid_out_14 = f'{test_output}/filter_invalid_read'
valid_log_14 = f'{test_output}/filter_invalid_read/filter_invalid_read.log'
os.mkdir(valid_out_14)
params_filter_invalid_read = f'--forward {invalid_read} --reverse {invalid_read} {params_valid_DB} \
                            {params_filter} --output {valid_out_14} \
                            --log {valid_log_14} \
                            &>> {valid_log_14}'

#filter_invalid_DB
valid_out_15 = f'{test_output}/filter_invalid_DB'
valid_log_15 = f'{test_output}/filter_invalid_DB/filter_invalid_DB.log'
os.mkdir(valid_out_15)
params_filter_invalid_DB = f'{params_valid_unpaired} --reference {invalid_DB} \
                            {params_filter} --output {valid_out_15} \
                            --log {valid_log_15} \
                            &>> {valid_log_15}'

#filter_invalid_outfolder
valid_log_16 = f'{test_output}/filter_invalid_outfolder.log'
params_filter_invalid_outfolder = f'{params_valid_unpaired} {params_valid_DB} \
                            {params_filter} --output /folder/not/exist \
                            --log {valid_log_16} \
                            &>> {valid_log_16}'


testFilterRef = [
    (params_help, 0),
    (params_filter_paired_unpaired, 0),
    (params_filter_paired, 0),
    (params_filter_unpaired, 0),
    (params_filter_invalid_read, 1),
    (params_filter_invalid_DB, 1),
    (params_filter_invalid_outfolder, 1)
]

@pytest.mark.parametrize("params, expected_status", testFilterRef)
def test_return_status_filtering(params, expected_status):
    status = subprocess.call(f'{workdir}/scripts/filter_refs.sh {params}', shell=True)
    assert expected_status == status

################################################################################
#Script: cluster_refs.sh

valid_inputs_refs = f"{test_output}/filter_paired_unpaired"

#Both paired and unpaired reads, filtering step was done before
valid_out_17 = f'{test_output}/cluster_paired_unpaired'
valid_log_17 = f'{test_output}/cluster_paired_unpaired/cluster_paired_unpaired.log'
os.mkdir(valid_out_17)
params_cluster_paired_unpaired = f'{params_valid_paired_unpaired} {params_valid_DB} \
                            --inputs {valid_inputs_refs} \
                            --filter_refs true \
                            {params_cluster} \
                            --scripts {workdir}/scripts --output {valid_out_17} \
                            --log {valid_log_17} \
                            &>> {valid_log_17}'

#Only paired reads, filtering step was not done before
valid_out_18 = f'{test_output}/cluster_paired'
valid_log_18 = f'{test_output}/cluster_paired/cluster_paired.log'
os.mkdir(valid_out_18)
params_cluster_paired = f'{params_valid_paired} {params_valid_DB} \
                        --inputs {valid_inputs_refs} \
                        --filter_refs false \
                        {params_cluster} \
                        --scripts {workdir}/scripts --output {valid_out_18} \
                        --log {valid_log_18} \
                        &>> {valid_log_18}'

#Only unpaired read, filtering step was done before
valid_out_19 = f'{test_output}/cluster_unpaired'
valid_log_19 = f'{test_output}/cluster_unpaired/cluster_unpaired.log'
os.mkdir(valid_out_19)
params_cluster_unpaired = f'{params_valid_unpaired} {params_valid_DB} \
                        --inputs {valid_inputs_refs} \
                        --filter_refs true \
                        {params_cluster} \
                        --scripts {workdir}/scripts --output {valid_out_19} \
                        --log {valid_log_19} \
                        &>> {valid_log_19}'

# cluster_invalid_read
valid_out_20 = f'{test_output}/cluster_invalid_read'
valid_log_20 = f'{test_output}/cluster_invalid_read/cluster_invalid_read.log'
os.mkdir(valid_out_20)
params_cluster_invalid_read = f'--forward {invalid_read} --reverse {invalid_read} {params_valid_DB} \
                        --inputs {valid_inputs_refs} \
                        --filter_refs true \
                        {params_cluster} \
                        --scripts {workdir}/scripts --output {valid_out_20} \
                        --log {valid_log_20} \
                        &>> {valid_log_20}' 

# cluster_invalid_DB
valid_out_21 = f'{test_output}/cluster_invalid_DB'
valid_log_21 = f'{test_output}/cluster_invalid_DB/cluster_invalid_DB.log'
os.mkdir(valid_out_21)
params_cluster_invalid_DB = f'{params_valid_paired_unpaired} --reference {invalid_DB} \
                        --inputs {valid_inputs_refs} \
                        --filter_refs true \
                        {params_cluster} \
                        --scripts {workdir}/scripts --output {valid_out_21} \
                        --log {valid_log_21} \
                        &>> {valid_log_21}'

# cluster_invalid_scripts
valid_out_22 = f'{test_output}/cluster_invalid_scripts'
valid_log_22 = f'{test_output}/cluster_invalid_scripts/cluster_invalid_scripts.log'
os.mkdir(valid_out_22)
params_cluster_invalid_scripts = f'{params_valid_paired_unpaired} {params_valid_DB} \
                        --inputs {valid_inputs_refs} \
                        --filter_refs true \
                        {params_cluster} \
                        --scripts /incorrect/path --output {valid_out_22} \
                        --log {valid_log_22} \
                        &>> {valid_log_22}'

# cluster_invalid_inputs_refs
valid_out_23 = f'{test_output}/cluster_invalid_inputs_refs'
valid_log_23 = f'{test_output}/cluster_invalid_inputs_refs/cluster_invalid_inputs_refs.log'
os.mkdir(valid_out_23)
params_cluster_invalid_inputs_refs = f'{params_valid_paired_unpaired} {params_valid_DB} \
                        --inputs /folder/not/exist \
                        --filter_refs true \
                        {params_cluster} \
                        --scripts {workdir}/scripts --output {valid_out_23} \
                        --log {valid_log_23} \
                        &>> {valid_log_23}'

# cluster_invalid_outfolder
valid_log_24 = f'{test_output}/cluster_invalid_outfolder.log'
params_cluster_invalid_outfolder = f'{params_valid_paired_unpaired} {params_valid_DB} \
                                --inputs {valid_inputs_refs} \
                                --filter_refs true \
                                {params_cluster} \
                                --scripts {workdir}/scripts --output /folder/not/exist \
                                --log {valid_log_24} \
                                &>> {valid_log_24}'

testClusterRef = [
    (params_help, 0),
    (params_cluster_paired_unpaired, 0),
    (params_cluster_paired, 0),
    (params_cluster_unpaired, 0),
    (params_cluster_invalid_read, 1),
    (params_cluster_invalid_DB, 1),
    (params_cluster_invalid_scripts, 1),
    (params_cluster_invalid_inputs_refs, 1),
    (params_cluster_invalid_outfolder, 1)
]

@pytest.mark.parametrize("params, expected_status", testClusterRef)
def test_return_status_clustering(params, expected_status):
    status = subprocess.call(f'{workdir}/scripts/cluster_refs.sh {params}', shell=True)
    assert expected_status == status

################################################################################
