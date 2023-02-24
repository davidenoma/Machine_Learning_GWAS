#!/bin/bash
#remove first column
cut -f2- -d"," $1  > hold
mv hold $1

#remove first row
awk NR\>1 $1 > hold
mv hold $1
#replace comma with space
sed -i 's/,/ /g' $1

#Generate the recode
java -Xmx2g -jar Jawamix5/jawamix5.jar char2num -ig call_method_54.tair9.FT10.csv -o char2num

#convert genotye tp hdf5
java -Xmx2g -jar Jawamix5/jawamix5.jar import -ig char2num -o genohdf5
#reate kinship matrix
java -Xmx2g -jar Jawamix5/jawamix5.jar kinship -ig genohdf5 -o kinship/

#Run Genome wide linear regression for phenotype(s)
java -Xmx2g -jar Jawamix5/jawamix5.jar lm -ig genohdf5 -ip FT10.txt -o emmax_/ -ik kinship/.rescaled.IBS

