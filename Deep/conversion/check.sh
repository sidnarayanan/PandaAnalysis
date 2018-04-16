#!/bin/bash

for f in ${SUBMIT_NPY}/cache/*; do
    PInfo -n check.sh $f
    check --cache $f $@
done
