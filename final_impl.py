from multiprocessing import freeze_support

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
    # phenotype = phenotype.dropna()
    phenotype = phenotype.dropna()
    true_pheno = phenotype.loc[:, 'ecotype_id'].values
    index = [str(x) for x in true_pheno]
    chrom = genotype.iloc[:,0].values
    pos = genotype.iloc[:,1].values

    i = 2
    genotype = genotype.iloc[:, i:][[col for col in genotype.columns[i:] if col in index]]
    i = 0
    while i < len(index):
        if index[i] not in genotype.columns:
            phenotype=phenotype[phenotype.ecotype_id != true_pheno[i]]
        i = i + 1

    genotype.insert(0, 'Chromosome',chrom,allow_duplicates=True)
    genotype.insert(1, 'Positions',pos)


    print(genotype.shape, phenotype.shape)

    return genotype,phenotype

def removeMAF(genotype_path):

    genotype = pd.read_csv(genotype_path, index_col=None)
    print('Before removing single snps with only one allele accross all samples:',genotype.shape)
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
    genotype.to_csv(ospath + '/genotype_final_after_maf.csv', index_label=None, index=False)
    return genotype


def write_to_file(genotype,phenotype):
    #This is neccessary to perform genome-wide association mapping  with Jawamax
    genotype.to_csv(ospath+'/genotype_final.csv',index_label=None,index=False)
    phenotype.to_csv(ospath+'/phenotype_final.txt',sep='\t',index_label=None,index=False)


def OLS_GWAS():
    genotype = pd.read_csv(ospath + '/genotype2num_final', index_col=None)
    # genotype = pd.read_csv(ospath + '/genotype2num_final', index_col=None)
    genotype = genotype.transpose()
    phenotype = pd.read_csv(ospath + '/phenotype_final.txt', index_col=None, sep="\t")
    phenotype = phenotype.iloc[:, 1].values
    genotype.set_axis([str(x) for x in genotype.iloc[1, :].values], axis="columns", inplace=True)
    genotype.drop(['Chromosome', 'Positions'], inplace=True)
    # read in coordinates after genotype mapping
    # add constant to X for intercept term
    print(phenotype.shape,genotype.shape)
    genotype = sm.add_constant(genotype)
    # fit OLS model
    model = sm.OLS(phenotype, genotype).fit()

    print(model.summary())
    output = pd.DataFrame(
        {'coefficients': model.params, 'std_error': model.bse, 't_values': model.tvalues, 'p_values': model.pvalues})
    output.index = genotype.columns
    # print output
    print(output)


def modelling():
    # Split the data into training and testing sets
    genotype = pd.read_csv(ospath + '/genotype2num_final', index_col=None)
    genotype = genotype.transpose()
    phenotype = pd.read_csv(ospath + '/phenotype_final.txt', index_col=None, sep="\t")
    genotype.set_axis([str(x) for x in genotype.iloc[1, :].values], axis="columns", inplace=True)
    genotype.drop(['Chromosome', 'Positions'], inplace=True)
    # read in coordinates after genotype mapping

    print(genotype.head(), genotype.shape)
    emmax = pd.read_csv(ospath + '/emmax_final/EMMAX.0_5_FT10.top.csv')
    # Selecting the top 100 Features affter sorting in ascending order based on p-value
    ##The sorting was done on excel.
    emmax_coord = [str(x) for x in emmax.iloc[:49, 1].values]



    lm = pd.read_csv(ospath+'/lm_final/LM.0_5_FT10.top.csv')
    lm_coord = [str(x) for x in lm.iloc[:49,1].values]
    # Split into train and test set
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, emmax_coord], phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)
    # # Fit the linear regression model
    pca = PCA(n_components=10)
    X_train_pca = pca.fit_transform(X_train)
    X_test_pca = pca.transform(X_test)
    print(pca.explained_variance_ratio_)
    # Create a LinearRegression object and fit it to the transformed training data
    lr = LinearRegression()
    lr.fit(X_train_pca, y_train)

    # Make predictions on the transformed testing data
    y_pred = lr.predict(X_test_pca)

    # Calculate the mean squared error of the predictions
    mse = mean_squared_error(y_test, y_pred)
    print(f"Mean Squared Error: {mse:.2f}")



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
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, emmax_coord], phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)

    # hidden_layer_sizes = (50, 10, 2)
    #100,50,30,2
    model = MLPRegressor(hidden_layer_sizes = (100,50, 30, ),random_state=1, activation="relu", solver="adam", max_iter=1000)
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
# modelling()
# genotype,phenotype = QConPreditionFiles(genotype_path,phenotype_path)
# write_to_file(genotype,phenotype)


OLS_GWAS()

