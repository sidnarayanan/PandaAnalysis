#!/bin/bash

for f in ${SUBMIT_NPY}/cache/*; do
    logger.info -n check.sh $f
    check --cache $f $@
done
