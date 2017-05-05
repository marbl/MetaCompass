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

#wget https://gembox.cbcb.umd.edu/metacompass/refseq.tar.gz
#cd wgs-8.3rc2
#cd kmer && make install && cd ..
#cd src && make && cd ..
#cd ..
#cmd="mkdir -p ./src/metaphyler/markers"
#echo $cmd
#$cmd

cmd="wget --no-check-certificate https://gembox.cbcb.umd.edu/metacompass/markers.tar.gz -P ./src/metaphyler"
echo $cmd
$cmd

cmd="tar -xzvf ./src/metaphyler/markers.tar.gz"
echo $cmd
$cmd


#cmd="tar -xzvf mc.tar.gz"
#echo $cmd
#$cmd

cmd="wget https://gembox.cbcb.umd.edu/metacompass/refseq.tar.gz"
echo $cmd
$cmd

cmd="tar -xzvf refseq.tar.gz"
echo $cmd
$cmd

cmd="wget https://gembox.cbcb.umd.edu/metacompass/test.tar.gz"
echo $cmd
$cmd

cmd="tar -xzvf test.tar.gz"
echo $cmd
$cmd

