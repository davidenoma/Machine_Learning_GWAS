import pandas as pd
from sklearn import linear_model, svm
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
#Loading the data files
from xgboost import XGBRegressor

ospath = "C:/Users/HP/OneDrive/Desktop/PhD._BMB/COURSE_WORK/MDGE612"

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
    gen_distance = [0]* genotype.shape[0]
    # chrom_snp_id_gen_distance_bp.insert(2,"gen distance",gen_distance)
    gene_T = genotype.iloc[:,3:].transpose()
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

# #try to replace info on the tped file tab
# example_ped = pd.read_csv(ospath+'/trans.tped',header=None,sep=' ')
#
# #Create a new tped file with duplicated alleles at each position
# #Multiply all alleles at each postion
# new_example_ped = pd.DataFrame()
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
# to_update_trans = pd.read_csv(ospath+'/trans.tped', sep=" ",header=None)
# to_update_trans.insert(2,'2',[0]*to_update_trans.shape[0])
# to_update_trans.columns = range(to_update_trans.shape[1])
# to_update_trans.to_csv(ospath+"/update_trans.tped",sep =" ",header=None)
def saveToDisk(genotype,phenotype):
    genotype.to_csv(ospath + '/genotype.csv', index=False)
    phenotype.to_csv(ospath + '/phenotype.csv', index=False, sep=" ")

def QConPreditionFiles(genotype,phenotype):
    genotype = pd.read_csv(genotype_path)
    phenotype = pd.read_csv(phenotype_path, sep="\t")
    phenotype = phenotype.dropna()

    true_pheno = phenotype.iloc[:, 0].values
    index = [str(x) for x in true_pheno]
    i = 2

    while i < len(genotype.columns):
        if genotype.columns[i] not in index:
            genotype.drop(columns=genotype.columns[i], inplace=True)
        i = i + 1
    while i < len(index):
        if index[i] not in genotype.columns:
            # final_index.append(index[i])
            print(index)
            index.remove(index[i])
        i = i + 1
    print(genotype.shape, phenotype.shape)
    return genotype,phenotype


# Split the data into training and testing sets

genotype = pd.read_csv(ospath + '/genotype2num', index_col=None)
genotype = genotype.transpose()
phenotype = pd.read_csv(ospath + '/phenotype.csv', index_col=None,sep="\t")
genotype.set_axis([str(x) for x in genotype.iloc[1,:].values],axis="columns",inplace=True)
# genotype.columms = genotype.iloc[1,:].values
genotype.drop(['Chromosome','Positions'],inplace=True)

#read in coordinates after genotype mapping
emmax = pd.read_csv(ospath+'/emmax_stepwise_clean/EMMAX.0_5_FT10.top')
emmax_coord = [str(x) for x in emmax.iloc[:,1].values]

lm = pd.read_csv(ospath+'/lm_clean/LM.0_5_FT10.top')
lm_coord = [str(x) for x in lm.iloc[:,1].values]

#Linear Regression
# X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:,emmax_coord], phenotype.iloc[:,1], test_size=0.2, random_state=0)

X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:,lm_coord], phenotype.iloc[:,1].values, test_size=0.2, random_state=0)
# Fit the linear regression model


reg = LinearRegression().fit(X_train, y_train)
# Predict the phenotype using the genotypes
y_pred = reg.predict(X_test)
# Calculate the mean squared error
mse = mean_squared_error(y_test, y_pred)
# Print the mean squared error
print("Mean Squared Error:", mse)

#Ridge regression

ridge_model = Ridge(alpha=1.0).fit(X_train, y_train)
# Predict the phenotype using the genotypes
y_pred = ridge_model.predict(X_test)
# Calculate the mean squared error
mse = mean_squared_error(y_test, y_pred)
# Print the mean squared error
print("Mean Squared Error:", mse)
print(genotype.loc[:,emmax_coord].shape, len(phenotype.iloc[:,1].values))


lasso_model = linear_model.Lasso(alpha=0.1).fit(X_train, y_train)
# Predict the phenotype using the genotypes
y_pred = lasso_model.predict(X_test)
# Calculate the mean squared error
mse = mean_squared_error(y_test, y_pred)
# Print the mean squared error
print("Mean Squared Error:", mse)


#SVR
svm_model = svm.SVR(epsilon=0.4).fit(X_train,y_train)
# Predict the phenotype using the genotypes
y_pred = svm_model.predict(X_test)
# Calculate the mean squared error
mse = mean_squared_error(y_test, y_pred)
print("Mean Squared Error:", mse)

#
model = XGBRegressor().fit(X_train, y_train)
y_pred = model.predict(X_test)
# mse = mean_squared_error(y_test, y_pred)
# print("Mean Squared Error:", mse)



