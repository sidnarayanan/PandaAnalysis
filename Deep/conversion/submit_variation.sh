#!/bin/bash

executable=$1

if [[ "$executable" == "" ]]; then 
    echo "Usage: $0 <executable>" >&2
    exit 1
fi

rm -rf cache
mkdir cache

#for v in Nominal ISRRenormUp ISRRenormDown; do
for v in FSRRenormUp FSRRenormDown; do
    submit --exec $executable --arglist partitions/QCD_${v}.txt --cache $(readlink -f cache/QCD_$v) 
done
