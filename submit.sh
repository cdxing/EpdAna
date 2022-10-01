#!/bin/bash

nJob1=$1
# sys_cutN

nJob2=$2
# sys_varN

nJob3=$3
# inter

#Filelist=$4
run0=$4

run1=$5

star-submit-template -template submit_template.xml -entities cutID=$nJob1,varID=$nJob2,iterID=$nJob3,run_start=$run0,run_end=$run1
#,listOfFiles=$Filelist
