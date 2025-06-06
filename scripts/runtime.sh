#!/bin/bash
N=`date +%s%N`
export PS4='+[$(((`date +%s%N`-$N)/1000000))ms][${BASH_SOURCE}:${LINENO}]: ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'; 
exec 2>runtime.log; set -x;
for ((i=0;i<10;i++)); do
    o=$(($RANDOM*$RANDOM/$RANDOM))
    echo $o
    sleep 0.$o
done

