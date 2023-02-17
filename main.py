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
gen_distance = [0]* genotype.shape[0]
chrom_snp_id_gen_distance_bp.insert(2,"gen distance",gen_distance)
gene_T = genotype.iloc[:,3:].transpose()

print(chrom_snp_id_gen_distance_bp,"\n",gene_T)
phenotype = pd.read_csv(phenotype_path,index_col=None)
#preparing the phenotype file
family_id = [1] * phenotype.shape[0]

individual_id = [0] * phenotype.shape[0]
paternal_id = [0] * phenotype.shape[0]
maternal_id = [1] * phenotype.shape[0]

phenotype.insert(0,'maternal_id',maternal_id)
phenotype.insert(0,'paternal_id',paternal_id)
phenotype.insert(0,'invidual_id',individual_id)
phenotype.insert(0,'family_id',family_id)

phenotype.to_csv('/home/davidenoma/Documents/MDGE612/phenotype_recode.csv',index=False)


