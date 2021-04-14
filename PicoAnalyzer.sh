#!/bin/bash
#  This script is for test run the PicoAnalyzer.cxx

#This is argument for iteration, if supplied
iter=1
if [ $# -ge 1 ]
then
        iter=$1
fi
stardev
root4star -b -q -l RunAnalyzer.C+
root4star -b -q -l PicoAnalyzer.cxx+'("st_physics_adc_19158057_raw_3500002.picoDst.root","test",1,0,0,$iter)'
# chmod u+x
