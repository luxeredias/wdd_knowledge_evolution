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

reactome_GMT <- read.gmt("data/Filtered_Reactome_2016.txt") 
genes_GMT <- unique(reactome_GMT$gene)

genes_diseases_converted <- alias2official_gene_symbol(genes = genes_diseases,
                                                       manualSearch = T)

genes_GMT_converted <- alias2official_gene_symbol(genes = genes_GMT,
                                                  manualSearch = T)

all_edges_converted <- merge(all_edges,genes_diseases_converted,
                             by.x="Source",by.y="original") %>%
  dplyr::mutate(Source=converted) %>%
  dplyr::select(-converted)

reactome_GMT_converted <- merge(reactome_GMT,genes_GMT_converted,
                                by.x="gene",by.y="original") %>%
  dplyr::mutate(gene=converted) %>%
  dplyr::select(-converted) %>%
  dplyr::select(ont,gene)

save(all_edges_converted,file="data/all_edges_converted.RData")
save(reactome_GMT_converted,file="data/reactome_GMT_converted.RData")
