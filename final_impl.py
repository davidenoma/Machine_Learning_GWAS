from multiprocessing import freeze_support
from statistics import mean

import numpy as np
import pandas as pd
from sklearn import linear_model, svm
from sklearn.decomposition import PCA
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor
import statsmodels.api as sm
from scipy import stats
import multiprocessing


ospath = "C:/Users/HP/OneDrive/Desktop/PhD._BMB/COURSE_WORK/MDGE612/"
# /home/davidenoma/Documents/MDGE612

genotype_path = ospath+"/call_method_54.tair9.FT10.csv"
phenotype_path = ospath+"/FT10.txt"

def QConPreditionFiles(genotype,phenotype):
    genotype = pd.read_csv(genotype,index_col=None)
    phenotype = pd.read_csv(phenotype, sep="\t")
    print(genotype.shape, phenotype.shape)

    #Remove NA from the phenotype file
    phenotype = phenotype.dropna()
    #get the phenotype ids for each of the plants
    true_pheno = phenotype.loc[:, 'ecotype_id'].values
    index = [str(x) for x in true_pheno]

    #Starting from index 2 to skip the chromosome and postions columns
    i = 2
    #select columns in the genotype file if they are in the index create of the phenotype ids
    genotype = genotype.iloc[:, i:][[col for col in genotype.columns[i:] if col in index]]

    i = 0
    #Removing the  phenotype ids that are not in the genotype file.
    while i < len(index):
        if index[i] not in genotype.columns:
            phenotype=phenotype[phenotype.ecotype_id != true_pheno[i]]
        i = i + 1



    print(genotype.shape, phenotype.shape)

    return genotype,phenotype

def removeMAF(genotype_path):
    genotype = pd.read_csv(genotype_path, index_col=None)
    print('Before removing single snps with only one allele across all samples:',genotype.shape)
    #remove snps with MAF less than 0.05
    i = 0
    while i < genotype.shape[0]:
        snps = genotype.iloc[i,2:].values
        unique_items = set(snps)
        if (len(unique_items) != 2):
            genotype = genotype.drop(i)
        i = i + 1
    print('After removing single snps with only one allele accross all samples:',genotype.shape)
    print('Before removing snps with MAF less than 0.05:',genotype.shape)
    genotype.set_axis([x for x in range(genotype.shape[0])], axis=0, inplace=True)

    i=0
    while i < genotype.shape[0]:
        snps = genotype.iloc[i, 2:].values
        unique_items = set(snps)
        individual_snps = len(snps)
        if (len(unique_items) == 2):
            snp1 = list(snps).count(list(unique_items)[0])
            snp2 = list(snps).count(list(unique_items)[1])
            if (snp1 / individual_snps < 0.05) or (snp2 / individual_snps < 0.05):
                genotype = genotype.drop(i)
                print(i, snp1 / individual_snps, snp2 / individual_snps)
        i = i + 1
    print('After removing snps with MAF less than 0.05:',genotype.shape)

    return genotype
def prepare_files_for_modelling_in_R():
    genotype = pd.read_csv(ospath + '/genotype_final_after_maf', index_col=None)
    print(genotype.shape)
    genotype = genotype.transpose()
    #remove columns chromosome and positions
    genotype.set_axis([str(x) for x in genotype.iloc[1, :].values], axis="columns", inplace=True)
    genotype.drop(['Chromosome', 'Positions'], inplace=True)
    # read in coordinates after genotype mapping
    emmax = pd.read_csv(ospath + '/emmax_maf/EMMAX.0_5_FT10.top.csv')
    ##The sorting was done on excel in ascending order from smallest p-value
    emmax_coord = [str(x) for x in emmax.iloc[:, 1].values]
    emmax_geno_for_R = genotype.loc[:,emmax_coord]
    print(emmax_geno_for_R.shape)

def write_to_file(genotype,phenotype):
    #This is neccessary to perform genome-wide association mapping  with Jawamax
    genotype.to_csv(ospath+'/genotype_final.csv',index_label=None,index=False)
    phenotype.to_csv(ospath+'/phenotype_final.txt',sep='\t',index_label=None,index=False)
#
# genotype,phenotype = QConPreditionFiles(genotype_path,phenotype_path)
# write_to_file(genotype,phenotype)
# removeMAF(ospath+'/genotype_final.csv')
prepare_files_for_modelling_in_R()