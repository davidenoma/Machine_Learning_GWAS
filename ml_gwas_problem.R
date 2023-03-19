#REQUIRED LIBRARIES
library(readr)
library(qqman)
# library(colortools)
library(tidyverse)
library(caret)
library(glmnet)
# library(OIdata) 
#Setting up the predicition model
ospath = ""
setwd(ospath)
set.seed(0)
#LOAD filtered gene file
gene <- read_csv("genotype_final_after_maf")
#LOAD filtered phenotype file
FT10 <- read_delim("phenotype_final.txt", "\t", escape_double = FALSE, trim_ws = TRUE)
EMMAX_0_5_FT10  = read_csv("emmax_maf/EMMAX.0_5_FT10.top.csv")

all = EMMAX_0_5_FT10[1:100,]
#Get p value of SNPs
pval = all$pvalue

filter=all
#Get the positions of SNPs
positions = filter$location

#Get the genomic data of significant SNPs from GWAS analysis
filtered_genes = gene[is.element(gene$Positions, positions),]

#Transpose data
trans = t(filtered_genes)

#Formulate the SNP ID as a Combination of the Chromosome and SNP position
locations = paste0(filtered_genes$Chromosome,".",filtered_genes$Positions)
#Rename the columns with the new SNP ID
colnames(trans)=locations
#Remove row one and two to elminate unwanted data
trans = trans[c(-1,-2),]
ecotype_id = rownames(trans)

#Combine the plant IDs to the main table to be used in the regression analysis
Full = cbind(ecotype_id,trans)
names(FT10)[2] <- "Phenotype"
Complete = merge(Full, FT10,by="ecotype_id",sort=FALSE)

sample_Size = nrow(Complete)


for(i in 1:1000){
train_Size = as.integer(sample_Size*0.8)
#Randomly select the plant IDs
train_Sample = sample(nrow(Complete),train_Size)
train_Data = Complete[train_Sample,-1]
test_Data = Complete[-train_Sample,-1]


x <- model.matrix(Phenotype~., train_Data)[,-1]
y <- train_Data$Phenotype

x.test <- model.matrix(Phenotype~., test_Data)[,-1]

#Ridge regression model alpha = 0
ridge <- cv.glmnet(x, y, alpha = 0)
ridge$lambda.min

model <- glmnet(x, y, alpha = 0, lambda = ridge$lambda.min)
predictions <- model %>% predict(x.test) %>% as.vector()

ridge_results = data.frame(
  MSE = RMSE(predictions, test_Data$Phenotype),
  Rsquare = R2(predictions, test_Data$Phenotype)
)

#Lasso regression model alpha = 1
lasso <- cv.glmnet(x, y, alpha = 1)
lasso$lambda.min
model <- glmnet(x, y, alpha = 1, lambda = lasso$lambda.min)
predictions <- model %>% predict(x.test) %>% as.vector()
lasso_results = data.frame(
  MSE = RMSE(predictions, test_Data$Phenotype),
  Rsquare = R2(predictions, test_Data$Phenotype)
)
}
lasso_results
ridge_results

