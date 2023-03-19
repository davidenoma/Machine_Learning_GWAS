from sklearn.decomposition import PCA
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

def recode_genotype(genotype):

    return

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


def OLS_GWAS():
    genotype = pd.read_csv(ospath + '/genotype_final_after_maf_recode', index_col=None)
    # genotype = pd.read_csv(ospath + '/genotype2num_final', index_col=None)
    genotype = genotype.transpose()
    phenotype = pd.read_csv(ospath + '/phenotype_final.txt', index_col=None, sep="\t")
    phenotype = phenotype.iloc[:, 1].values
    genotype.set_axis([str(x) for x in genotype.iloc[1, :].values], axis="columns", inplace=True)
    genotype.drop(['Chromosome', 'Positions'], inplace=True)
    # read in coordinates after genotype mapping
    # add constant to X for intercept term

    emmax = pd.read_csv(ospath + '/emmax_maf/EMMAX.0_5_FT10.top.csv')
    # Selecting the top 100 Features affter sorting in ascending order based on p-value
    ##The sorting was done on excel.
    emmax_coord = [str(x) for x in emmax.iloc[:300, 1].values]

    genotype = genotype.loc[:,emmax_coord].to_numpy()
    print(genotype.shape,phenotype.shape)
    genotype= sm.add_constant(genotype)
    # fit OLS model
    model = sm.OLS(phenotype, genotype).fit()
    print(model.summary())
    output = pd.DataFrame(
        {'coefficients': model.params, 'std_error': model.bse, 't_values': model.tvalues, 'p_values': model.pvalues})
    # output.index = genotype.columns
    # print output
    print(output.head())


def modelling():
    # Split the data into training and testing sets
    genotype = pd.read_csv(ospath + '/genotype_final_after_maf', index_col=None)
    genotype = genotype.transpose()
    phenotype = pd.read_csv(ospath + '/phenotype_final.txt', index_col=None, sep="\t")
    genotype.set_axis([str(x) for x in genotype.iloc[1, :].values], axis="columns", inplace=True)
    genotype.drop(['Chromosome', 'Positions'], inplace=True)
    # read in coordinates after genotype mapping


    emmax = pd.read_csv(ospath + '/emmax_maf/EMMAX.0_5_FT10.top.csv')
    # Selecting the top 100 Features affter sorting in ascending order based on p-value
    ##The sorting was done on excel.
    emmax_coord = [str(x) for x in emmax.iloc[:, 1].values]
    lm = pd.read_csv(ospath+'/lm_maf/LM.0_5_FT10.top.csv')
    lm_coord = [str(x) for x in lm.iloc[:4200,1].values]


    # emmax_geno_for_R = genotype.loc[:,emmax_coord]
    # print(emmax_geno_for_R.head(),emmax_geno_for_R.to_csv(ospath+'emmax_index_genotype.csv',index_label=None,index=False))

    # Split into train and test set
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:,emmax_coord], phenotype.iloc[:, 1].values,
                                                   test_size=0.2, random_state=0)

    max_index, pcs = PCA_impl(X_test, X_train, y_test, y_train)
    pca = PCA(n_components=pcs[max_index], svd_solver='full')
    X_train_pca = pca.fit_transform(X_train)
    X_test_pca = pca.transform(X_test)
    # print(pca.explained_variance_ratio_)

    reg = LinearRegression().fit(X_train_pca, y_train)
    # Predict the phenotype using the genotypes
    y_pred = reg.predict(X_test_pca)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    # print("Mean Squared Error: Linear Regression", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    # print("R-squared for : Linear Regression", r_squared)
    print("\n")

    alphas =  np.arange(0, 10.001, 0.001).tolist()
    r2s_l1 =[]
    r2s_l2 = []
    r2s_l1,r2s_l2 = OptimizeAlpha(X_test_pca, X_train_pca, alphas,r2s_l1,r2s_l2,y_test, y_train)
    print(max(r2s_l1),alphas[r2s_l1.index(max(r2s_l1))])
    print(max(r2s_l2),alphas[r2s_l2.index(max(r2s_l2))])

    #Lasso
    lasso_r2s_scores = []
    for i in range(0,1000):
        lasso_model = linear_model.Lasso(alpha=alphas[r2s_l1.index(max(r2s_l1))]).fit(X_train_pca, y_train)
        # Predict the phenotype using the genotypes
        y_pred = lasso_model.predict(X_test_pca)
        # Calculate the mean squared error
        mse = mean_squared_error(y_test, y_pred)
        # Print the mean squared error
        # print("Mean Squared Error Lasso L1 Regression:", mse)
        # Calculate R-squared
        r_squared = r2_score(y_test, y_pred)
        lasso_r2s_scores.append(r_squared)
        # print("R-squared for Lasso L1 Regression:", r_squared)

    # Ridge regression
    ridge_r2s_scores = []

    for i in range(0,1000):
        ridge_model = Ridge(alpha=alphas[r2s_l2.index(max(r2s_l2))]).fit(X_train_pca, y_train)
        # Predict the phenotype using the genotypes
        y_pred = ridge_model.predict(X_test_pca)
        # Calculate the mean squared error
        mse = mean_squared_error(y_test, y_pred)
        # Print the mean squared error
        # print("Mean Squared Error Ridge Regression :L2 ", mse)
        # Calculate R-squared
        r_squared = r2_score(y_test, y_pred)
        ridge_r2s_scores.append(r_squared)
        # print("R-squared for Ridge Regression L2: ", r_squared)
        # print("\n")
    print( mean(lasso_r2s_scores),mean(ridge_r2s_scores),max(lasso_r2s_scores),max(ridge_r2s_scores))
    return


    # SVR
    svm_model = svm.SVR(epsilon=0.2).fit(X_train_pca, y_train)
    # Predict the phenotype using the genotypes
    y_pred = svm_model.predict(X_test_pca)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    print("Mean Squared Error(SVR):", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for SVR:", r_squared)
    print("\n")

    #For complex models, it might be useful to have more features
    # Split into train and test set
    # emmax_coord = [str(x) for x in emmax.iloc[:, 1].values]
    # lm_coord = [str(x) for x in lm.iloc[:1999, 1].values]
    # X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, emmax_coord], phenotype.iloc[:, 1].values,
    #                                                     test_size=0.2, random_state=0)

    model = XGBRegressor().fit(X_train_pca, y_train)
    y_pred = model.predict(X_test_pca)
    mse = mean_squared_error(y_test,y_pred)
    print("Mean Squared Error (XGBOOST):", mse)
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print("R-squared for XGBOOST:", r_squared)
    print("\n")
    return
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


def OptimizeAlpha(X_test_pca, X_train_pca, alphas, r2s_l1,r2s_l2, y_test, y_train):

    for i in range(len(alphas)):
        lasso_model = linear_model.Lasso(alpha=alphas[i]).fit(X_train_pca, y_train)
        # Predict the phenotype using the genotypes
        y_pred = lasso_model.predict(X_test_pca)
        # Calculate the mean squared error
        mse = mean_squared_error(y_test, y_pred)
        # Print the mean squared error
        # print("Mean Squared Error Lasso L1 Regression:", mse)
        # Calculate R-squared
        r_squared = r2_score(y_test, y_pred)
        r2s_l1.append(r_squared)
        # print("R-squared for Lasso L1 Regression:", r_squared)
    # Ridge regression

    for i in range(len(alphas)):
        ridge_model = Ridge(alpha=alphas[i]).fit(X_train_pca, y_train)
        # Predict the phenotype using the genotypes
        y_pred = ridge_model.predict(X_test_pca)
        # Calculate the mean squared error
        mse = mean_squared_error(y_test, y_pred)
        # Print the mean squared error
        # print("Mean Squared Error Ridge Regression :L2 ", mse)
        # Calculate R-squared
        r_squared = r2_score(y_test, y_pred)
        r2s_l2.append(r_squared)
        # print("R-squared for Ridge Regression L2: ", r_squared)
        # print("\n")
    return  r2s_l1,r2s_l2


def PCA_impl(X_test, X_train, y_test, y_train):
    r2s = []
    pcs = []
    for i in range(10, 120):
        # # Fit the linear regression model
        pca = PCA(n_components=i, svd_solver='full')
        X_train_pca = pca.fit_transform(X_train)
        X_test_pca = pca.transform(X_test)
        lr = LinearRegression()
        lr.fit(X_train_pca, y_train)
        y_pred = lr.predict(X_test_pca)
        mse = mean_squared_error(y_test, y_pred)
        r_squared = r2_score(y_test, y_pred)
        r2s.append(r_squared)
        pcs.append(i)
    max_index = r2s.index(max(r2s))
    print(len(pcs), len(r2s), max(r2s), pcs[max_index])
    return max_index, pcs


#Implementation line
modelling()
