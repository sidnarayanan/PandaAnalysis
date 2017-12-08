#!/bin/bash

LAST=1

while true
do
    if [[ $(ls -A $SUBMIT_OUTDIR/locks/ | wc -l) > 0 ]]; then
        RECENT=$(stat -c%Z $(ls -ltr $SUBMIT_OUTDIR/locks/* | tail -n1 | awk '{ print $9; }' 2>/dev/null) 2>/dev/null)
    else
        RECENT=1
    fi
    if [[ "$LAST" -le "$RECENT" || $(expr $(date +%s) - $LAST) > 60 ]]; 
    then
        clear
        ${CMSSW_BASE}/src/PandaAnalysis/T3/bin/checkMissingFiles.py
        LAST=$(date +%s) # this takes care of the 60s refresh
        date
    fi
    sleep 5
done
