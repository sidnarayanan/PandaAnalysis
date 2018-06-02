#!/bin/bash 

for e in convert.py convert_higgs.py; do 
    echo $e
    check --exec $e $@ 
done
