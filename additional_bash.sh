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
java -Xmx2g -jar Jawamix5/jawamix5.jar kinship -ig genohdf5_final -o kinship_final/

#Run Genome wide linear regression for phenotype(s)
java -Xmx2g -jar Jawamix5/jawamix5.jar lm -ig genohdf5_final -ip phenotype_final.txt -o lm_final/ -ik kinship_final/.rescaled.IBS
java -Xmx2g -jar Jawamix5/jawamix5.jar emmax_stepwise -ig genohdf5_final -ip phenotype_final.txt -o emmax_final/ -ik kinship_final/.rescaled.IBS

#Annotation of output
awk 'BEGIN{FS=OFS="\t"}{if(NR>1){ print $1,$2-1,$2 }}'  emmax_snps_for_bedtools_final.csv | bedtools intersect -wb -a gene_model_update.gff -b stdin  > annotation_of_top_snps