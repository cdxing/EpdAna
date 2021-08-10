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
wait
root4star -b -q -l RunAnalyzer.C+
wait
root4star -b -q -l PicoAnalyzer.cxx+\(\"hlt_22031042_10_01_000.picoDst.root\",\"test\",1,$cut,$var,$iter\)
wait
# chmod u+x
