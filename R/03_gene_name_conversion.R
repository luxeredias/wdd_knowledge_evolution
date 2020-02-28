#script to convert all gene names in WDD files to official gene symbol
#and all gene names in Reactome GMT file to official gene symbol

#packages
library(dplyr)
library(clusterProfiler)
source("~/Área de Trabalho/GitHub_projects/alias2official/scripts/alias2official_gene_symbol.R")

#set working directory
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#load edges file
data("all_edges")

all_edges <- bind_rows(edges_list)
genes_diseases <- unique(all_edges$Source)

genes_diseases_converted <- alias2official_gene_symbol(genes = genes_diseases,
                                                       manualSearch = T)

all_edges_converted <- merge(all_edges,genes_diseases_converted,
                             by.x="Source",by.y="original") %>%
  dplyr::mutate(Source=converted) %>%
  dplyr::select(-converted)

save(all_edges_converted,file="data/all_edges_converted.RData")

reactome_GMT <- read.gmt("data/Filtered_Reactome_2016.txt")
genes_GMT <- unique(reactome_GMT$gene)

genes_GMT_converted <- alias2official_gene_symbol(genes = genes_GMT,
                                                  manualSearch = T)


reactome_GMT_converted <- merge(reactome_GMT,genes_GMT_converted,
                                by.x="gene",by.y="original") %>%
  dplyr::mutate(gene=converted) %>%
  dplyr::select(-converted) %>%
  dplyr::select(ont,gene)

save(reactome_GMT_converted,file="data/reactome_GMT_converted.RData")


cannonical_GMT <- read.gmt("data/cannonical_all_msigdb.gmt")
genes_cannonical <- unique(cannonical_GMT$gene)

genes_cannonical_converted <- alias2official_gene_symbol(genes = genes_cannonical,
                                                       manualSearch = T)

cannonical_GMT_converted <- merge(cannonical_GMT,genes_cannonical_converted,
                                  by.x="gene",by.y="original") %>%
  dplyr::mutate(gene=converted) %>%
  dplyr::select(-converted) %>%
  dplyr::select(ont,gene)

save(cannonical_GMT_converted,file="data/cannonical_GMT_converted.RData")

reactome_full_GMT <- read.gmt("data/ReactomePathwaysAll.gmt")
genes_reactome_all <- unique(reactome_full_GMT$gene)

gene_reactome_all_converted <- alias2official_gene_symbol(genes = genes_reactome_all,
                                                          manualSearch = T)
gene_reactome_all_converted$converted <- ifelse(gene_reactome_all_converted$converted=="",
                                                gene_reactome_all_converted$original,
                                                gene_reactome_all_converted$converted)
reactome_full_GMT_converted <- merge(reactome_full_GMT,
                                     gene_reactome_all_converted,
                                     by.x="gene",by.y="original") %>%
  dplyr::mutate(gene=converted) %>%
  dplyr::select(-converted) %>%
  dplyr::select(ont,gene)

save(reactome_full_GMT_converted,file="data/reactome_full_GMT_converted.RData")
