import pandas as pd

#Loading the data files
ospath = "C:/Users/HP/OneDrive/Desktop/PhD. BMB/COURSE_WORK/MDGE612"

test_genotype_file = ospath+"/individual_ids"
genotype_path = ospath+"/call_method_54.tair9.FT10.csv"
phenotype_path = ospath+"/FT10.txt"

#creation of ped and fam files
#File specification at: https://zzz.bwh.harvard.edu/plink/data.shtml
def genotype_cleaning():
    genotype = pd.read_csv(genotype_path)
    snpids = []
    for i in range(genotype.shape[0]):
        snpids.append("SNP" + str(i))
    genetic_distance = [0] * genotype.shape[0]
    genotype.insert(1, 'SNPIDS', snpids)

    # genotype.to_csv(ospath+'/final_genotype.csv')

    # chrom_snp_id_gen_distance_bp = genotype.iloc[:,:3]
    # gen_distance = [0]* genotype.shape[0]
    # chrom_snp_id_gen_distance_bp.insert(2,"gen distance",gen_distance)
    # gene_T = genotype.iloc[:,3:].transpose()
    #
    # print(chrom_snp_id_gen_distance_bp,"\n",gene_T)


def phenotype_cleaning():
    phenotype = pd.read_csv(phenotype_path, index_col=None)
    # #preparing the phenotype file
    family_id = [1] * phenotype.shape[0]

    individual_id = [0] * phenotype.shape[0]
    paternal_id = [0] * phenotype.shape[0]
    maternal_id = [1] * phenotype.shape[0]
    #
    phenotype.insert(0, 'maternal_id', maternal_id)
    phenotype.insert(0, 'paternal_id', paternal_id)
    phenotype.insert(0, 'invidual_id', individual_id)
    phenotype.insert(0, 'family_id', family_id)

    # phenotype.to_csv(ospath+'/phenotype_recode.csv',index=False)

#try to replace info on the tped file tab
example_ped = pd.read_csv(ospath+'/example.tped',header=None,index_col=None,sep=' ')
# print(example_ped)

for i in range(example_ped.shape[0]*2):
    (example_ped.iloc[i,:])





#Multiply all alleles at each postion