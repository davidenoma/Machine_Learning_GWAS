#remove first column
cut -f2- -d"," genotype.tped  > trans.tped

#remove first row
 awk NR\>1 final_genotype.csv > genotype.tped
#replace comma with space
sed -i 's/,/ /g' trans.tped