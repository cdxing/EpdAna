#!/bin/bash
for name in merged_EpCorrection_OUTPUT*
do
    newname="$(echo "$name" | cut -9- )"
    mv "$name" "$newname"
done
for name in orrection_OUTPUT*
do
    newname="$(echo "$name" | cut -9- )"
    mv "$name" "$newname"
done
for name in _OUTPUT*
do
    newname=EpCorrection_INPUT"$(echo "$name" | cut -7- )"
    mv "$name" "$newname"
done
