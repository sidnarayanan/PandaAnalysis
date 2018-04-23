#!/bin/bash

executable=$1

if [[ "$executable" == "" ]]; then 
    echo "Usage: $0 <executable>" >&2
    exit 1
fi

rm -rf ${SUBMIT_NPY}/cache
mkdir ${SUBMIT_NPY}/cache

#for f in ZpWW ZpTT ZpA0h QCD; do
for f in Top_lo QCD; do
#for f in VBFHbb ZvvHbb ggZvvHbb; do
    submit --exec $executable --arglist ${SUBMIT_NPY}/partitions/${f}.txt --cache $(readlink -f ${SUBMIT_NPY}/cache/$f) 
done
