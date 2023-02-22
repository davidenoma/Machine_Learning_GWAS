#!/bin/bash
#remove first column
cut -f2- -d"," $1  > hold
mv hold $1

#remove first row
awk NR\>1 $1 > hold
mv hold $1
#replace comma with space
sed -i 's/,/ /g' $1