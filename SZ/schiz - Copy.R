setwd("C:/Academic/Research/01 SZ")

library(tidyr)
library(caret)
require(ABAData)
library(ABAEnrichment)
library(mygene)
library(splitstackshape)


data("dataset_adult")

data = read.csv2("nature13595-s3.csv", sep=",")
genes = data$eQTL.gene
data = read.csv2("nature13595-s3-blood.csv", sep=",")
blood= data$eQTL.gene

ng = as.data.frame(as.character(genes))
ng1 = as.data.frame(as.character(blood))
names(ng1) = names(ng)
ng = rbind(ng, ng1)
ng = unlist(ng)

df = spread(dataset_adult, structure, signal)
df = df[complete.cases(df),]
duplicate_index = duplicated(df$entrezgene)
df = df[!duplicate_index,]
df$schiz = df$hgnc_symbol %in% ng
df = df[,c(1:4, 419, 5:418)]
genes1 = df$hgnc_symbol

d1 = df[, c(5:419)]
set.seed(200)
q = createDataPartition(d1$schiz, p=0.7, list = F)
training = d1[q,]
testing = d1[-q,]
mod.glm = glm(schiz~., data=training, family=binomial())

glm.probs = predict(mod.glm, type = "response")
glm.pred = rep(F, dim(training)[1])
limit = median(glm.probs)
glm.pred[glm.probs>limit] = T
z=table(glm.pred, training$schiz)
z
sensitivity(as.factor(glm.pred), as.factor(training$schiz))
specificity(as.factor(glm.pred), as.factor(training$schiz))

test.probs = predict(mod.glm, testing, type="response")
test.pred = rep(F, dim(testing)[1])
test.pred[test.probs>limit] = T
z=table(test.pred, testing$schiz)
z
sensitivity(as.factor(test.pred), as.factor(testing$schiz))
specificity(as.factor(test.pred), as.factor(testing$schiz))


res = queryMany(genes1, scopes='symbol', fields=c('entrezgene', 'hgnc', 'go'), species='human')
# 
# z.BP = rep(NA, dim(res)[1])
# for(index in 1:dim(res)[1]){
#   z = res$go.BP[[index]]$term[1]
#   if(class(z) == "character"){
#     z.BP[index] = z
#   }
# }
# z.BP = as.data.frame(z.BP)
# 
# z.CC = rep(NA, dim(res)[1])
# for(index in 1:dim(res)[1]){
#   z = res$go.CC[[index]]$term[1]
#   if(class(z) == "character"){
#     z.CC[index] = z
#   }
# }
# z.CC = as.data.frame(z.CC)
# 
# 
# z.MF = rep(NA, dim(res)[1])
# for(index in 1:dim(res)[1]){
#   z = res$go.MF[[index]]$term[1]
#   if(class(z) == "character"){
#     z.MF[index] = z
#   }
# }
# z.MF = as.data.frame(z.MF)
# 
# 
# GO.frame = data.frame(entrezgene = res$entrezgene, BP = z.BP$z.BP, CC = z.CC$z.CC, MF = z.MF$z.MF)
# GO.frame = GO.frame[complete.cases(GO.frame),]
# GO_duplicate_index=duplicated(GO.frame$entrezgene)
# GO.frame = GO.frame[!GO_duplicate_index,]
# 
# df = merge(GO.frame, df)
# df = df[,c(1, 5:8, 2:4, 9:422)]
# 
z = res
z1=sapply(z$go.BP, "[[", "term")
lenz1 = lapply(z1, length)
z3 = data.frame(genes = z$query, freq = unlist(lenz1))
z4 = expandRows(z3, "freq")
z5 = cbind(genes=z4, MF=unlist(z1))
z6 = z5[z5$genes %in% ng,]
length(unique(z$genes))
# 
