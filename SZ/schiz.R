setwd("C:/Academic/Research/01 SZ")

library(tidyr)
library(caret)
require(ABAData)
library(ABAEnrichment)
library(mygene)
library(splitstackshape)
library(glmnet)


data("dataset_adult")


df = spread(dataset_adult, structure, signal)
df = df[,-2]
df = df[complete.cases(df),]
duplicate_index = duplicated(df$hgnc_symbol)
df = df[!duplicate_index,]


#using genes from PMC4754044
ngenes = read.csv("PMC4754044_table3.csv", sep=",")
names(ngenes) = c("hgnc_symbol", "pval", "variants")
#index = df$hgnc_symbol %in% ngenes$hgnc_symbol
index = df$hgnc_symbol %in% sz$REPORTED.GENE.S.
df1 = df[index,]
df1 = merge(df1, sz, by.x ="hgnc_symbol", by.y="REPORTED.GENE.S.")
#df1 = df1[,c(1:3, 418, 5:417)]
#df2 = df1[,-c(1:3)]

x = model.matrix(pval~., df2)[,-1]
y=df2$pval
grid = 10^seq(10, -2, length=100)
set.seed(1)
train=sample(1:nrow(x), nrow(x)*0.7)
test=(-train)
lasso.mod = glmnet(x[train,], y[train], alpha=1, lambda=grid)
