#!/bin/bash

now=$( date '+%Y%m%d%H%M%S' )
echo Starting on $( date )

for j in $( seq -w 1 100 ); do
  for i in $( seq 0 8 ); do
    beginCent=$[$i]
    endCent=$[$beginCent]
    echo '******************************************************'
    echo Submitting for centralities $beginCent and $endCent
    echo '******************************************************'
    star-submit-template -u ie -template runjob.xml -entities nEvts=10853300000,startCent=$beginCent,endCent=$endCent,outName=LC.toyMC-${beginCent}-${endCent}.$j.root,home=/gpfs/mnt/gpfs01/star/pwg/msimko/LambdaCdataDrivenSImulator,JobId=$now  
  done
done
