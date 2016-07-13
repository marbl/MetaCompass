#!/usr/bin/bash
 
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

cmd="wget https://gembox.cbcb.umd.edu/metacompass/refseq.tar.gz"
echo $cmd
$cmd

cmd="tar -xzvf refseq.tar.gz"
echo $cmd
$cmd
