#!/bin/bash

now=$( date '+%Y%m%d%H%M%S' )
echo Starting on $( date )

for j in $( seq -w 1 100 ); do
  for i in $( seq 0 7 ); do
    beginCent=$[$i]
    endCent=$[$beginCent]
    echo '******************************************************'
    echo Submitting for centralities $beginCent and $endCent
    echo '******************************************************'
    star-submit-template -template runjob.xml -entities nEvts=10853300000,startCent=$beginCent,endCent=$endCent,outName=LC.toyMC-${beginCent}-${endCent}.$j.root,home=/global/project/projectdirs/star/pwg/starhf/simkomir/LambdaC/DataDrivenFastSim,JobId=$now.$j  
  done
done
