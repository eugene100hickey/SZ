setwd("~/Desktop/Academic/Research/Schiz/SZ")

library(tidyverse)
library(limma)
require(ABAData)
library(ABAEnrichment)
library(colorspace)

#reads in list of SZ genes from Nature publication
#genes stored in object ng
gene1 <- suppressWarnings(read_csv("nature13595-s3.csv")$eQTL)
gene2 <- suppressWarnings(read_csv("nature13595-s3-blood.csv")$eQTL)
ng <- c(gene1, gene2)



#reads in ABA Data
data("dataset_adult")

df  <-  spread(dataset_adult, structure, signal)
df  <-  df[complete.cases(df),]
duplicate_index  <-  duplicated(df$entrezgene)
df  <-  df[!duplicate_index,]

d1_all_genes = df[, c(5:418)]

df$schiz = df$hgnc_symbol %in% ng

#function that gets all brain superstructure areas from a given brain code
brainX <-  function(x){
  get_name(get_superstructures(x))
}

#gets the brain codes from the dataframe column names
brain_labels <- names(d1_all_genes)

#runs the function brainX above to get the brain areas,
#builds a dataframe called brain with the codes and the areas
brain <- sapply(brain_labels, brainX) %>% 
  unlist()
brain = data.frame(code = names(brain), area = brain)

#groups brain areas at a level two below overall brain
#we can change this if we think we have too many groups
#the telencephalon accounts for half the brain areas, so we subdivide further
#likewise the metencephalon
#finish with 17 groups
index = which(brain$area == "Br_Brain")
#index = index + ifelse(brain$area[index+2] == "Tel_Telencephalon", 2,0)
#index = index + ifelse(brain$area[index+2] == "MET_Metencephalon", 1,0)
brain = brain[index+2,]
brain$area <- as.character(brain$area)

brain$telenceph <- brain$area == "Tel_Telencephalon"
design <- model.matrix(~0+schiz, data = df)
cm <- makeContrasts(status = schizTRUE - schizFALSE, levels = design)
fit <- lmFit(dft, design)
fit1 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit1)
results <- decideTests(fit2)
summary(results)

plotDensities(d1_all_genes, legend = F)

design <- model.matrix(~0+schiz, data = df)
cm <- makeContrasts(status = telencephTRUE - telencephFALSE, levels = design)
fit <- lmFit(d1_all_genes, design)
fit1 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit1)
results <- decideTests(fit2)
summary(results)

plotDensities(d1_all_genes, legend = F)
plotMDS(d1_all_genes, labels = brain$area, 
        col = rainbow_hcl(7)[c(as.factor(brain$area))], 
        cex = 0.5, gene.selection = "common")
top_areas <- topTable(fit2, number = 6)
dft[rownames(top_areas),1:5]

stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
hist(stats[, "P.Value"])


z <- normalizeBetweenArrays(d1_all_genes)
plotDensities(normalizeBetweenArrays(d1_all_genes), legend = F)
keep <- rowMeans(z) > 5
z <- z[keep,]
plotDensities(z, legend = F)

plotMDS(d1_all_genes, labels = brain$area, 
        col = rainbow_hcl(7)[c(as.factor(brain$area))], 
        cex = 0.5, gene.selection = "common")
# legend("topleft", pch=16, cex = 0.4, col=rainbow_hcl(7),
#        legend=unique(brain$area))

top_genes <- topTable(fit2, number = 3)
df[rownames(top_genes),1:5]

stats <- topTable(fit2, number = nrow(fit2), sort.by = "none")
hist(stats[, "P.Value"])

volcanoplot(fit2, highlight = 10, names = df$hgnc_symbol)

enrich_kegg <- kegga(fit2, geneid = df$entrezgene, species = "Hs")
topKEGG(enrich_kegg, number = 10)

enrich_go <- goana(fit2, geneid = df$entrezgene, species = "Hs")
topGO(enrich_go, number = 10, ontology = "BP")
