#!/bin/bash

# Provide list of comma seperated bins or files to look for

check_bins(){

	for b in $to_check
	do
		t1=$(which $b)
		if [ $? == 0 ]; then
			echo "$b found."
		else
			echo "$b not found."
			exit 1
		fi
	done
}


check_files(){

	for f in $to_check
	do
		if [ -r $f ]; then
			echo "$f found and readable."
		else
			echo "$f doesn't exist or is not readable."
			exit 1
		fi
	done

}

usage(){
	echo "Usage: ./check_exist.sh [ -b <bins in PATH> | -f <filepaths> ]"
	echo "Ex: ./check_exist.sh -b nextflow,bowtie2"
	echo "Ex: ./check_exist.sh -f /home/example.txt,./example2.txt"
}

to_check=$(echo $2 | tr "," " ")
echo "Checking for: $to_check"

if [ $1 == "-b" ]; then

	check_bins

elif [ $1 == "-f" ]; then

	check_files

else

	usage

fi