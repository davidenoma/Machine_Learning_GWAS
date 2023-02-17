import pandas as pd

#Loading the data files

test_genotype_file = "/home/davidenoma/Documents/MDGE612/individual_ids"
genotype_path = "/home/davidenoma/Documents/MDGE612/call_method_54.tair9.FT10.csv"
phenotype_path = "/home/davidenoma/Documents/MDGE612/FT10.txt"

#creation of ped and fam files

genotype = pd.read_csv(test_genotype_file)
snpids = []
for i in range(genotype.shape[0]):
    snpids.append("SNP"+str(i))
genotype.insert(1,'SNPIDS',snpids)
# genotype.to_csv('/home/davidenoma/Documents/MDGE612/update_to_ind_id_file.csv')
chrom_snp_id_gen_distance_bp = genotype.iloc[:,:3]
gene_T = genotype.iloc[:,3:].transpose()
print(chrom_snp_id_gen_distance_bp,gene_T)



