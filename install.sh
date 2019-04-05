#!/usr/bin/env bash
 
echo "#Installing MetaCompass"
cmd="g++ -Wall -W -O2 -o ./bin/extractSeq ./src/utils/extractSeq.cpp"
echo $cmd
$cmd

cmd="g++ -Wall -W -O2 -o ./bin/formatFASTA ./src/utils/formatFASTA.cpp"
echo $cmd
$cmd

cmd="g++ -Wall -W -O2 -o ./bin/buildcontig ./src/buildcontig/buildcontig.cpp ./src/buildcontig/cmdoptions.cpp ./src/buildcontig/memory.cpp ./src/buildcontig/procmaps.cpp ./src/buildcontig/outputfiles.cpp"
echo $cmd
$cmd

#downloaded on Jan 11th 2019
cmd="wget --no-check-certificate https://obj.umiacs.umd.edu/metacompassdb/2019/01/refseq.tar.gz"
echo $cmd
$cmd

cmd="tar -xzvf refseq.tar.gz"
echo $cmd
$cmd

#optional
cmd="rm refseq.tar.gz"
echo $cmd
$cmd

