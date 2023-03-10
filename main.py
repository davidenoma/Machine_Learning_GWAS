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

def modelling():
    # Split the data into training and testing sets
    genotype = pd.read_csv(ospath + '/updated_recode_file', index_col=None)
    genotype = genotype.transpose()


    # Print the explained variance ratios of the principal components

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
    X_train, X_test, y_train, y_test = train_test_split(genotype, phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)
    pca = PCA(n_components=100)
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
    r_squared = r2_score(y_test, y_pred)
    print(f"Mean Squared Error PCA: {mse:.2f}","Rsquared: ",r_squared)

    # Split into train and test set
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:,emmax_coord], phenotype.iloc[:, 1].values,
                                                        test_size=0.2, random_state=0)

    # Fit the linear regression model
    reg = LinearRegression().fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = reg.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    print(f"Mean Squared Error: Linear Regression: {mse:.2f}")
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print(f"R-squared for : Linear Regression {r_squared:.2f}")
    print("\n")
    # Ridge regression

    ridge_model = Ridge(alpha=1.0).fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = ridge_model.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    print(f"Mean Squared Error: Ridge Regression L2: {mse:.2f}")

    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print(f"R-squared for Ridge Regression L2:  {r_squared:.2f}")
    print("\n")

    lasso_model = linear_model.Lasso(alpha=0.1).fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = lasso_model.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    # Print the mean squared error
    print(f"Mean Squared Error Lasso L1 Regression: { mse:.2f}")
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print(f"R-squared for Lasso L1 Regression: {r_squared:.2f}")
    print("\n")

    # SVR
    svm_model = svm.SVR(epsilon=0.4).fit(X_train, y_train)
    # Predict the phenotype using the genotypes
    y_pred = svm_model.predict(X_test)
    # Calculate the mean squared error
    mse = mean_squared_error(y_test, y_pred)
    print(f"Mean Squared Error(SVR): {mse:.2f}")
    # Calculate R-squared
    r_squared = r2_score(y_test, y_pred)
    print(f"R-squared for SVR: {r_squared:.2f}")
    print("\n")




    #For complex models, it might be useful to have more features
    # Split into train and test set
    emmax_coord = [str(x) for x in emmax.iloc[:, 1].values]
    lm_coord = [str(x) for x in lm.iloc[:1999, 1].values]
    X_train, X_test, y_train, y_test = train_test_split(genotype.loc[:, emmax_coord], phenotype.iloc[:, 1].values,                                                        test_size=0.2, random_state=0)

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
modelling()
