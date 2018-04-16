#!/bin/bash

CATALOG=${HOME}/catalog/t3mit/pandaf/008/

echo "cataloging snarayan"

cd /mnt/hadoop/scratch/snarayan/gridpanda/

for f in * ; 
do 
    pd=$(echo $f | sed -e "s?ZPrimeToTTJets?ZprimeToTT?" -e "s?GeV??g" -e "s?_M?_M-?" -e "s?_W?_W-?" -e "s?\$?_TuneCUETP8M1_13TeV-madgraphMLM-pythia8?")
    echo -n "$f -> $pd "
    mkdir -p ${CATALOG}/${pd}
    ls $PWD/$f/*root 2>/dev/null | sort -u |  sed "s?/mnt/hadoop?root://t3serv006.mit.edu/?" > ${CATALOG}/${pd}/RawFiles.00  
    echo "-> $(cat ${CATALOG}/${pd}/RawFiles.00 | wc -l) files"
done

echo "cataloging bmaier"

cd /mnt/hadoop/scratch/bmaier/gridpanda/

for f in * ; 
do 
    pd=$(echo $f | sed -e "s?ZPrimeToTTJets?ZprimeToTT?" -e "s?GeV??g" -e "s?_M?_M-?" -e "s?_W?_W-?" -e "s?\$?_TuneCUETP8M1_13TeV-madgraphMLM-pythia8?")
    echo -n "$f -> $pd "
    mkdir -p ${CATALOG}/${pd}
    ls $PWD/$f/*root 2>/dev/null | sort -u |  sed "s?/mnt/hadoop?root://t3serv006.mit.edu/?" >> ${CATALOG}/${pd}/RawFiles.00  
    echo "-> $(cat ${CATALOG}/${pd}/RawFiles.00 | wc -l) files"
done

echo "total of $(cat ${CATALOG}/ZprimeToTT*/RawFiles.00 | wc -l) files"
