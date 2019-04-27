# linking our SZ genes to pathways

library(Homo.sapiens)
library(AnnotationDbi)
library(tidyverse)
library(magrittr)
library(KEGGprofile)
library(rvest)
library(topGO)
library(limma)
library(ABAData)
library(beepr)
# library(hgu95av2.db)
# library(org.Hs.eg.db)


setwd("~/Desktop/Academic/Research/Schiz/SZ")

sz_genes <- c(suppressWarnings(read_csv("nature13595-s3.csv")$eQTL), 
              suppressWarnings(read_csv("nature13595-s3-blood.csv")$eQTL)
)

data("dataset_adult")

sz_entrezgenes <- dataset_adult %>% 
  filter(hgnc_symbol %in% sz_genes) %>% 
  select(entrezgene) %>% 
  distinct() %>% 
  pull(entrezgene) %>% 
  as.character()

sz_KEGGresult <- find_enriched_pathway(sz_entrezgenes,species = 'hsa', 
                                      download_latest = T, 
                                      returned_pvalue = 0.05, 
                                      returned_adjpvalue = 0.2, 
                                      returned_genenumber = 2)
sz_KEGGresult <- sz_KEGGresult[[1]] 
sz_KEGGresult %<>% rownames_to_column(var = "PATHWAY")
View(sz_KEGGresult)
beep(9)

# setwd("~/Desktop/Academic/Science4/Projects/2019")
# autism_genes <- read_csv("autism_genes_sfari.csv")
# sz_genes <- autism_genes$Autism_Genes



z <- AnnotationDbi::select(org.Hs.eg.db,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "PATH", "MAP"), keys = sz_genes) %>%
  distinct()

chroms <- AnnotationDbi::select(org.Hs.eg.db,
                           keytype = "SYMBOL",
                           columns = c("SYMBOL", "MAP"), keys = sz_genes) %>%
  distinct()
# 
# chroms$CHROM <- str_split(chroms$MAP, 
#                           pattern = "[p,q]", 
#                           n = 2, 
#                           simplify = T)[,1]
chroms$CHROM <- ifelse(str_detect(chroms$SYMBOL, pattern = "C.{1,}orf"),
                       parse_number(chroms$SYMBOL),
                       parse_number(chroms$MAP))

chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "TMEM22", 3))
chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "SEPW1", 19))
chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "RP5-874C20.3", 20))
chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "USMG5", 10))
chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "FAM63A", 1))
chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "DUS2L", 16))
chroms <- chroms %>% mutate(CHROM = replace(CHROM, SYMBOL == "ZNFX1-AS1", 20))



url <- "https://www.genome.jp/kegg/pathway.html"
w <- read_html(url)
pathway_data_html <- html_nodes(w, "dd a , b+ .list dt")
pathway_data <- html_text(pathway_data_html)
pathway_numbers <- pathway_data[c(T, F)]
pathway_names <- pathway_data[c(F, T)]
pathways <- data_frame(PATH = pathway_numbers, names = pathway_names)

pathway_all_html <- html_nodes(w, ".clear+ b , h4+ b , dd a , b+ .list dt , a+ b")
pathway_all <- html_text(pathway_all_html)
pathway_all <- data_frame(entry = pathway_all) %>%
  filter(!entry %in% pathways$PATH) %>%
  rownames_to_column() %>%
  filter(!entry %in% pathways$names) %>%
  mutate(rowname = as.numeric(rowname))

group_lengths <- c(
  as.numeric(pathway_all$rowname),
  dim(pathways)[1] + dim(pathway_all)[1] + 1
)
group_lengths <- group_lengths[-1]
group_lengths <- group_lengths - pathway_all$rowname - 1
pathways <- cbind(pathways, category = rep(pathway_all$entry, group_lengths))

pathway_groups_html <- html_nodes(w, "a+ h4")
pathway_groups <- html_text(pathway_groups_html)
pw <- strsplit(x = as.character(pathways$category), split = ".", fixed = T) %>%
  unlist()
pw <- pw[c(T, F)]
pathways <- cbind(pathways, cat_number = pw)

pwg <- strsplit(pathway_groups, ".", fixed = T)
pwg <- unlist(pwg)
pwg <- pwg[-length(pwg)]
pwg <- data_frame(cat_number = pwg[c(T, F)], cat_name = pwg[c(F, T)])
pathways <- left_join(pathways, pwg)

z1 <- getKEGGPathwayNames("hsa", remove.qualifier = T)

sz_paths <- left_join(z, pathways)