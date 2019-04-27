setwd("~/Desktop/Academic/Research/Schiz/SZ")

library(tidyverse)
require(ABAData)
library(ABAEnrichment)
library(WGCNA)
library(ggdendro)

#read in ABAData
data("dataset_adult")

#read in SZ genes from Ripke
gene1 <- suppressWarnings(read_csv("nature13595-s3.csv")$eQTL)
gene2 <- suppressWarnings(read_csv("nature13595-s3-blood.csv")$eQTL)
sz_genes <- c(gene1, gene2)

#place ABAData in wide format
allen <- dataset_adult %>% 
  filter(hgnc_symbol %in% sz_genes) %>% 
  select(hgnc_symbol, structure, signal) %>% 
  distinct() %>% 
  spread(key = structure, value = signal) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "hgnc_symbol")

brainX <-  function(x){
  get_name(get_superstructures(x))
}

# make "brain" with area descriptors for brain ABA areas
brain_labels <- names(allen)
brain <- sapply(brain_labels, brainX) %>% 
  unlist()
brain = data.frame(code = names(brain), area = brain)

# can modify this paragraph to get different refinements of brain area
index <- which(brain$area == "Br_Brain")
index <- index + ifelse(brain$area[index + 2] == "Tel_Telencephalon", 2, 0)
index <- index + ifelse(brain$area[index + 2] == "MET_Metencephalon", 1, 0)
brain <- brain[index + 2, ]
brain$area <- as.character(brain$area)
brain$code <- gsub("\\.....|\\......", "", brain$code)

# back to WCGNA
datExpr0 <- t(allen)

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree
ggdendrogram(sampleTree, rotate = T) +
  geom_hline(yintercept = 15, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust == 1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Form a data frame analogous to expression data that will hold 
# the clinical traits.
datTraits <- brain %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "code")
# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
cols <- datTraits %>% 
  select(area) %>% 
  as.factor() %>% 
  levels()
cols <- str_split(string = cols, pattern = "\",")[[1]] %>% 
  length() %>% 
  rainbow()
traitColors <- cols[as.factor(datTraits$area)]
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = rownames(datTraits),
                    main = "Sample dendrogram and trait heatmap")





net = blockwiseModules(datExpr, power = 6,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3)
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
