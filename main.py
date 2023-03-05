import numpy as np
import pandas as pd
from sklearn import linear_model, svm
from sklearn.linear_model import LinearRegression, Ridge
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor

ospath = "C:/Users/HP/OneDrive/Desktop/PhD._BMB/COURSE_WORK/MDGE612"
#/home/davidenoma/Documents/MDGE612

genotype_path = ospath+"/call_method_54.tair9.FT10.csv"
phenotype_path = ospath+"/FT10.txt"

#creation of ped and fam files
#File specification at: https://zzz.bwh.harvard.edu/plink/data.shtml
def genotype_cleaning():
    genotype = pd.read_csv(genotype_path)
    snpids = []
    for i in range(genotype.shape[0]):
        snpids.append("SNP" + str(i))
    genotype.insert(1,"snpid",snpids)

    gen_distance = [0] * genotype.shape[0]
    genotype.insert(2,"gen_distance",gen_distance)

    genotype.to_csv(ospath+'/new_genotype.tped',sep=" ",index=False)



def phenotype_cleaning():
    phenotype = pd.read_csv(phenotype_path, index_col=None, sep="\t")
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

    phenotype.to_csv(ospath+'/new_phenotype.tfam',index=False,sep=" ")

def saveToDisk(genotype,phenotype):
    genotype.to_csv(ospath+'/new_genotype.tped', sep=" ",index=False)
    phenotype.to_csv(ospath+'/new_phenotype.tfam', index=False, sep=" ")

def QConPreditionFiles(genotype,phenotype):
    genotype = pd.read_csv(genotype,sep=" ",index_col=None)
    phenotype = pd.read_csv(phenotype, sep=" ")
    phenotype = phenotype.dropna()

    true_pheno = phenotype.loc[:, 'ecotype_id'].values
    index = [str(x) for x in true_pheno]
    i = 4
    print(genotype.shape,phenotype.shape)

    while i < len(genotype.columns):
        if genotype.columns[i] not in index:
            print('found genotype missing')
            genotype.drop(columns=genotype.columns[i], inplace=True)
        i = i + 1

    i = 0
    while i < len(index):
        if index[i] not in genotype.columns:
            phenotype=phenotype[phenotype.ecotype_id != true_pheno[i]]
            print('found pheno missing ')
        i = i + 1
    print(genotype.shape, phenotype.shape)
    return genotype,phenotype


def modelling():
    # Split the data into training and testing sets

    genotype = pd.read_csv(ospath + '/updated_recode_file', index_col=None)
    genotype = genotype.transpose()
    phenotype = pd.read_csv(ospath + '/phenotype_for_jawa', index_col=None, sep="\t")
    genotype.set_axis([str(x) for x in genotype.iloc[1, :].values], axis="columns", inplace=True)
    # genotype.columms = genotype.iloc[1,:].values
    genotype.drop(['Chromosome', 'Positions'], inplace=True)
    # read in coordinates after genotype mapping
    emmax = pd.read_csv(ospath + '/emmax_step_wise_manual_recode/EMMAX.0_5_FT10.top.csv')

    # Selecting the top 100 Features affter sorting in ascending order based on p-value
    ##The sorting was done on excel.
    emmax_coord = [str(x) for x in emmax.iloc[:19, 1].values]

    lm = pd.read_csv(ospath+'/lm_manual_recode/LM.0_5_FT10.top.csv')
    lm_coord = [str(x) for x in lm.iloc[:49,1].values]

    # Split into train and test set
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, emmax_coord], phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)
    # Fit the linear regression model

    reg = LinearRegression().fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = reg.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    print("Mean Squared Error: Linear Regression", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for : Linear Regression", r_squared)
    print("\n")
    # Ridge regression

    ridge_model = Ridge(alpha=1.0).fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = ridge_model.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    print("Mean Squared Error Ridge Regression :L2 ", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for Ridge Regression L2: ", r_squared)
    print("\n")

    lasso_model = linear_model.Lasso(alpha=0.1).fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = lasso_model.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    print("Mean Squared Error Lasso L1 Regression:", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for Lasso L1 Regression:", r_squared)
    print("\n")

    # SVR
    svm_model = svm.SVR(epsilon=0.4).fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = svm_model.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    print("Mean Squared Error(SVR):", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for SVR:", r_squared)
    print("\n")




    #For complex models, it might be useful to have more features
    # Split into train and test set
    emmax_coord = [str(x) for x in emmax.iloc[:, 1].values]
    lm_coord = [str(x) for x in lm.iloc[:1999, 1].values]
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, emmax_coord], phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)

    model = XGBRegressor().fit(X_train.to_numpy(), y_train)
    y_pred = model.predict(X_test.to_numpy())
    mse = mean_squared_error(y_test,y_pred)
    print("Mean Squared Error (XGBOOST):", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for XGBOOST:", r_squared)
    print("\n")

    # For complex models, it might be useful to have more features
    # Split into train and test set
    emmax_coord = [str(x) for x in emmax.iloc[:, 1].values]
    lm_coord = [str(x) for x in lm.iloc[:1999, 1].values]
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, lm_coord], phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)

    # hidden_layer_sizes = (50, 10, 2)
    model = MLPRegressor(hidden_layer_sizes = (100,50, 30, 2),random_state=1, activation="relu", solver="adam", max_iter=1000)
    # Train model
    model.fit(X_train, y_train)

    # Make predictions on test data
    y_pred = model.predict(X_test)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for Neural network:", r_squared)
    mse = mean_squared_error(y_test, y_pred)
    print("Mean Squared Error ( Neural network:):", mse)
    print("\n")

#Implementation line
modelling()
