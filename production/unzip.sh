#!/bin/bash

count=1

mkdir temp

for tarFile in LC.toyMC-*.*.root.tar; do
  # echo unzipping $tarFile
  tar -xzf $tarFile
  if [ "$(( $count % 100 ))" -eq 0 ]; then
    echo Hadd to temp/LC.toyMC.$count.temp.root
    hadd temp/LC.toyMC.$count.temp.root *.root >> /dev/null
    rm -f *.root
  fi
  let "count++" 
done

if [ -f *.root ]; then
  echo Hadd to temp/LC.toyMC.$count.temp.root
  hadd temp/LC.toyMC.$count.temp.root *.root >> /dev/null
  rm -f *.root
fi

hadd LC.toyMC.root temp/*.temp.root

rm -rf temp
