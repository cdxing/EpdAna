#!/bin/bash
#  This script is for test run the PicoAnalyzer.cxx

#This is argument for iteration, if supplied
iter=1
if [ $# -ge 1 ]
then
        iter=$1
fi
stardev
wait
root4star -b -q -l RunAnalyzer.C+
wait
root4star -b -q -l PicoAnalyzer.cxx+\(\"./st_physics_20179040_raw_2500006.picoDst.root\",\"test\",1,0,0,$iter\)
wait
# chmod u+x
