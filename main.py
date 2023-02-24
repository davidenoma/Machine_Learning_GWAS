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
example_ped = pd.read_csv(ospath+'/trans.tped',header=None,sep=' ')

#Create a new tped file with duplicated alleles at each position
#Multiply all alleles at each postion
new_example_ped = pd.DataFrame()
def mulitply_alleles():
    for i in range(example_ped.shape[0]):
        print(i)
        x = example_ped.iloc[i, 4:]
        new_list = list()
        [new_list.extend([z] * 2) for z in list(x)]
        if i == 0:
            new_example_ped = pd.DataFrame(new_list)
            new_example_ped = new_example_ped.transpose()
        else:
            update_ped = pd.DataFrame(new_list)
            update_ped = update_ped.transpose()
            new_example_ped = new_example_ped.append(update_ped, ignore_index=True)
    print(example_ped.iloc[:, :4].shape, new_example_ped.shape)
    combined_final_tped = pd.merge(example_ped.iloc[:, :4], new_example_ped, left_index=True, right_index=True)
    combined_final_tped.columns = range(combined_final_tped.shape[1])
    print(combined_final_tped.shape)
    combined_final_tped.to_csv(ospath + '/new_trans.tped')


#Trying to insert genetic distance to tped file
to_update_trans = pd.read_csv(ospath+'/trans.tped', sep=" ",header=None)

to_update_trans.insert(2,'2',[0]*to_update_trans.shape[0])
to_update_trans.columns = range(to_update_trans.shape[1])
to_update_trans.to_csv(ospath+"/update_trans.tped",sep =" ",header=None)
