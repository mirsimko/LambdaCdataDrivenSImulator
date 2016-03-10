#!/bin/bash

count=0
for iCent in $( seq 0 7 ); do
  for iFile in $( seq -w 1 100 ); do
    if [ ! -f LC.toyMC-$iCent-$iCent.$iFile.root.tar ]; then
      echo LC.toyMC-$iCent-$iCent.$iFile.root.tar is missing
      let "count++"
    fi
  done
done

if [ $count == 0 ]; then
  echo No file is missing ... life is good 
fi
