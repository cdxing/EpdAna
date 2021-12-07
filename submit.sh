#!/bin/bash

nJob1=$1
# sys_cutN

nJob2=$2
# sys_varN

nJob3=$3
# inter

Filelist=$4

star-submit-template -template submit_template.xml -entities cutID=$nJob1,varID=$nJob2,iterID=$nJob3,listOfFiles=$Filelist
