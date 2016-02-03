#!/bin/bash

now=$( date '+%Y%m%d%H%M%S' )
echo Starting on $( date )

for i in $( seq 0 3 ); do
  beginCent=$[$i*2]
  endCent=$[$beginCent+1]
  echo '******************************************************'
  echo Submitting for centralities $beginCent and $endCent
  echo '******************************************************'
  star-submit-template -u ie -template runjob.xml -entities nEvts=10853300000,startCent=$beginCent,endCent=$endCent,outName=LC.toyMC-${beginCent}-${endCent}.root,home=/gpfs/mnt/gpfs01/star/pwg/msimko/LambdaCdataDrivenSImulator,JobId=$now  
done
