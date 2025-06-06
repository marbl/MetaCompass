#!/usr/bin/env bash

echo "##################################"
echo "##### Installing MetaCompass #####"
echo "##################################"


echo
echo
echo "##### Creating Executables #####"
echo

command -v scripts/cull_reference_candidates/refcull &> /dev/null
mod=$?
if [[ $mod != 0 ]]; then
    cmd="cd scripts/cull_reference_candidates"
    $cmd
    cmd="g++ -std=c++11 cmdopt.cpp cull_reference_candidate.cpp process_map.cpp -o refcull"
    $cmd
    cmd=$?
    if [[ $cmd != 0 ]]; then
        echo "refcull could not be created."
    else
        echo "refcull created."
    fi
    cd -
else
    echo "refcull already created."
fi

command -v ./scripts/FRC/get_FRC_data &> /dev/null
mod=$?
if [[ $mod != 0 ]]; then
    cmd="cd scripts/FRC"
    $cmd
    cmd="make"
    $cmd
    cmd=$?
    if [[ $cmd != 0 ]]; then
        echo "FRC could not be created."
    else
        echo "FRC created."
    fi
    cd -
else
    echo "get_FRC_data already created."
fi

# check that bins are in path blastn, bowtie2, samtools, bedtools, java, jellyfish, 
check_bins(){

	for b in $to_check
	do
		t1=$(which $b)
		if [[ $? == 0 ]]; then
			echo "$b found."
		else
			echo "$b not found."
			exit 1
		fi
	done
}

echo
echo
echo "##### Checking MetaCompass Dependencies #####"
echo

# Need to check separately: MetaCarvel, Python3.9
# NOTE: No blastn?
to_check="nextflow kmer-mask datasets bowtie2 samtools java jellyfish bedtools samtools megahit"
check_bins

echo
echo
echo "##### Checking Python Dependencies #####"
echo

echo "Python version is $(python -c 'import sys; print(".".join(map(str, sys.version_info[:3])))'), make sure it is at least 3.9."
    
python -c "import numpy" &> /dev/null
mod=$?
if [[ $mod != 0 ]]; then
    echo "numpy not installed."
else
    echo "numpy installed."
fi

python -c "import pandas" &> /dev/null
mod=$?
if [[ $mod != 0 ]]; then
    echo "pandas not installed."
else
    echo "pandas installed."
fi

python -c "import matplotlib" &> /dev/null
mod=$?
if [[ $mod != 0 ]]; then
    echo "matplotlib not installed."
else
    echo "matplotlib installed."
fi
echo
echo
echo "##### Reminder to set paths correctly #####"
echo

echo "Please make sure frc_path is set to $(pwd)/get_FRC_data in the nextflow.config file."
echo "Please make sure metacarvel_path is set correctly in the nextflow.config file."
echo "Please make sure reference_db is set correctly in the nextflow.config file."
echo
