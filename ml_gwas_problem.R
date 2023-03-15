#REQUIRED LIBRARIES
library(readr)
library(qqman)
# library(colortools)
library(tidyverse)
library(caret)
library(glmnet)
# library(OIdata) 
#Setting up the predicition model
ospath = "C:/Users/HP/OneDrive/Desktop/PhD._BMB/COURSE_WORK/MDGE612/"
setwd(ospath)
set.seed(0)
##MANHATTAN PLOT GENERATION

#LOAD unifiltered EMMAX file to get all the SNPS
# ALL_EMMAX <- read_csv("Data/ALL_EMMAX.top")
# #LOAD filtered EMMAX file to get the significant SNPS
# EMMAX_0_5_FT10 <- read_csv("Data/EMMAX.0_5_FT10.top")
# 
# 
# snp_Name = paste0(ALL_EMMAX$`#chr`,".",ALL_EMMAX$location)
# FULL = cbind(ALL_EMMAX,snp_Name)
# interest = paste0(EMMAX_0_5_FT10$`#chr`,".",EMMAX_0_5_FT10$location)
# 
# colour = c("#D42930","#D4CD29","#29D4CD","#2930D4","#CD29D4")
# colour=colour[-3]
# 
# #THE SIGNIFICANT SNPS USED FOR THE ANALYSIS WILL BE HIGHLIGHTED IN GREEN
# manhattan(FULL, chr="#chr", bp="location", snp="snp_Name", p="pvalue",highlight = interest,col= colour,
#           suggestiveline = F, genomewideline = F)
#=


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

R_ridge = vector()
R_lasso = vector()
R2_linear_regression = vector()

# for(i in 1:100){
  # print(i)
  #Number of plants to include in the training set which is 80% if the full data size
  train_Size = as.integer(sample_Size*0.8)
  #Randomly select the plant IDs
  train_Sample = sample(nrow(Complete),train_Size)
  train_Data = Complete[train_Sample,-1]
  test_Data = Complete[-train_Sample,-1]
  
  
  x <- model.matrix(Phenotype~., train_Data)[,-1]
  y <- train_Data$Phenotype
  
  x.test <- model.matrix(Phenotype~., test_Data)[,-1]
  
  # #trying the linear model
  # model <- glm(y ~ x, family = gaussian)
  # 
  # predictions <- model %>% predict(x.test) %>% as.vector()
  # 
  # data.frame(
  #   RMSE = RMSE(predictions, test_Data$Phenotype),
  #   Rsquare = R2(predictions, test_Data$Phenotype)
  # )
  # 
  # Rsquare_LM = R2(predictions, test_Data$Phenotype)
  # R2_linear_regression=c(R2_linear_regression,Rsquare_LM)
  
  
  
  #Ridge regression model depicted by alpha = 0
  ridge <- cv.glmnet(x, y, alpha = 0)
  ridge$lambda.min
  
  model <- glmnet(x, y, alpha = 0, lambda = ridge$lambda.min)
  #coef(model)
  
   predictions <- model %>% predict(x.test) %>% as.vector()
  data.frame(
    RMSE = RMSE(predictions, test_Data$Phenotype),
    Rsquare = R2(predictions, test_Data$Phenotype)
  )
  
  Rsquare_ridge = R2(predictions, test_Data$Phenotype)
  R_ridge=c(R_ridge,Rsquare_ridge)
  
  
  
  #Lasso regression model depicted by alpha = 1
  lasso <- cv.glmnet(x, y, alpha = 1)
  lasso$lambda.min
  model <- glmnet(x, y, alpha = 1, lambda = lasso$lambda.min)
  predictions <- model %>% predict(x.test) %>% as.vector()
  
  data.frame(
    RMSE = RMSE(predictions, test_Data$Phenotype),
    Rsquare = R2(predictions, test_Data$Phenotype)
  )
  
  Rsquare_lasso = R2(predictions, test_Data$Phenotype)
  R_lasso=c(R_lasso,Rsquare_lasso)
  
# }

print("Mean R Square of Ridge Regression:")
mean(R_ridge)
print("Mean R Square of Lasso Regression:", mean(R_lasso))
mean(R_lasso)

