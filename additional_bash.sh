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
java -Xmx2g -jar Jawamix5/jawamix5.jar char2num -ig genotype.csv -o genotype2num

#convert genotye tp hdf5
java -Xmx2g -jar Jawamix5/jawamix5.jar import -ig genotype2num -o genohdf5
#reate kinship matrix
java -Xmx2g -jar Jawamix5/jawamix5.jar kinship -ig genohdf5 -o kinship_clean/

#Run Genome wide linear regression for phenotype(s)
java -Xmx2g -jar Jawamix5/jawamix5.jar lm -ig genohdf5 -ip phenotype.csv -o lm_clean/ -ik kinship_clean/.rescaled.IBS
java -Xmx2g -jar Jawamix5/jawamix5.jar emmax_stepwise -ig genohdf5 -ip phenotype.csv -o emmax_stepwise_clean/ -ik kinship_clean/.rescaled.IBS

