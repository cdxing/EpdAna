#!/bin/sh
#  This script is for test run the PicoAnalyzer.cxx

#This is argument for iteration, if supplied
cut=0
var=0
iter=1
if [ $# -ge 1 ]
then
        cut=$1
	var=$2
	iter=$3
fi
source /afs/rhic/rhstar/group/star_cshrc.csh
wait
#echo "Hello, world!"
stardev
wait
root4star -b -q -l RunAnalyzer.C+
wait
root4star -b -q -l PicoAnalyzer.cxx+\(\"st_physics_adc_19158057_raw_3500002.picoDst.root\",\"test\",1,$cut,$var,$iter\)
wait
# chmod u+x
