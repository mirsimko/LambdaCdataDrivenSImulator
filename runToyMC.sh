#!/bin/bash

nPart=${1:-1085330000} # number of produced events after vz cut in p+K+pi channel
startCent=${2:-0}
endCent=${3:-2}
result=${4:-test.root}
root -l -b -q 'toyMcEffLc.C++g('$nPart','$startCent','$endCent',"'$result'")'
