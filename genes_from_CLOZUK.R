require(tidyverse)
require(httr)
require(readxl)
require(ABAData)
require(factoextra)
require(fpc)
require(WGCNA)


url <- "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5918692/bin/NIHMS958804-supplement-Supplementary_Table.xlsx"
GET(url, write_disk(temp_file <- tempfile(fileext = ".xlsx"))) # downloads the .xlsx file
df <- read_excel(temp_file, sheet = 5, skip = 6) # reads into a dataframe. First six rows of the excel file are just header
unlink(temp_file)     # deletes the temporary file


##################################################################
# makes sz_genes, a dataframe with a single column of the CLOZUK genes
#######################################################################

sz_genes <- df$`Gene(s) tagged` %>% 
  str_split(",") %>% 
  unlist() %>% 
  as.data.frame() %>% 
  distinct()
names(sz_genes) <- "genes"

data("dataset_adult")


########################################################################################
# makes wide_allen_sz, a dataframe which has sz_genes as rownames, brain areas as columns
# it's wide format and scaled
########################################################################################

wide_allen_sz <-
  dataset_adult %>% 
  select(hgnc_symbol, structure, signal) %>% 
  filter(hgnc_symbol %in% sz_genes$genes) %>% 
  distinct() %>% 
  spread(key = structure, value = signal) %>% 
  remove_rownames() %>% 
  column_to_rownames(var = "hgnc_symbol") %>% 
  t() %>% 
  scale() %>% 
  t() %>% 
  as.data.frame()


#########################################################
# a couple of toy examples. The first one plots a heatmap.
# striping here would be indicative of clusters
#########################################################
distance <- get_dist(wide_allen_sz, method = "euclidean", stand = F)
fviz_dist(distance, 
          gradient = list(low="#00AFBB", mid = "white", high = "#FC4E07"), 
          lab_size = 6)


##########################################################
# cluster plot, choose 7 clusters, used kmeans
###########################################################
km <- kmeans(wide_allen_sz, centers = 7, nstart = 25)
chooseTopHubInEachModule(t(wide_allen_sz), km$cluster, power = 9)
fviz_cluster(km, wide_allen_sz, repel = T)

########################################################
# clusterboot measures how stable each cluster is
# values above 0.75 usually mean a real cluster
# values below 0.5 usually mean a rag-tag grouping
#######################################################
cboot <- clusterboot(wide_allen_sz, 
                     clustermethod = kmeansCBI, 
                     bootmethod="boot", 
                     k=7, B=100, 
                     count=F)

cboot$bootmean
