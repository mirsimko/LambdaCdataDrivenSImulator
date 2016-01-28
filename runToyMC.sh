#!/bin/bash

nPart=${1:-2000000}
root -l -b -q 'toyMcEffLc.C++g('$nPart',0,1,"test.root")'
