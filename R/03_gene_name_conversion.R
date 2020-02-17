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
data("human_gene_info")

all_edges <- bind_rows(edges_list)
genes_diseases <- unique(all_edges$Source)

reactome_GMT <- read.gmt("data/Filtered_Reactome_2016.txt") 
genes_GMT <- unique(reactome_GMT$gene)

genes_diseases_converted <- alias2official_gene_symbol(genes = genes_diseases,
                                                       manualSearch = T)











genes_human <- unique(human_gene_info$Symbol)

dis_GMT <- intersect(genes_diseases,genes_GMT)
dis_notGMT <- genes_diseases[!genes_diseases %in% dis_GMT]
GMT_notdis <- genes_GMT[!genes_GMT %in% dis_GMT]

dis_genes <- data.frame(to="dis",from=genes_diseases)
GMT_genes <- data.frame(to="GMT",from=genes_GMT)
human_genes <- data.frame(to="human",from=genes_human)

all_edges <- as.data.frame(rbind(dis_genes,GMT_genes,human_genes))

graph <- graph_from_data_frame(all_edges,directed = F)
degree <- as.data.frame(igraph::degree(graph))

degree <- degree %>%
  rownames_to_column(var="node")

exclusive_genes_dis <- degree %>%
  filter(`igraph::degree(graph)`==1 & 
           !(node %in% human_genes$from) &
           !(node %in% GMT_genes$from)) %>%
  pull(node) %>%
  unique()

exclusive_genes_GMT <- degree %>%
  filter(`igraph::degree(graph)`==1 & 
           !(node %in% human_genes$from) &
           !(node %in% dis_genes$from)) %>%
  pull(node) %>%
  unique()

exclusive_genes_dis_converted <- alias2official_gene_symbol(exclusive_genes_dis,
                                                            manualSearch = T)

exclusive_genes_GMT_converted <- alias2official_gene_symbol(exclusive_genes_GMT,
                                                            manualSearch = T)

table(exclusive_genes_dis_converted$converted %in% solved_genes)

solved_genes_dis_GMT <- degree %>%
  filter(`igraph::degree(graph)`==2 & !(node %in% human_genes$from)) %>%
  pull(node) %>%
  unique()

solved_genes_dis_GMT_human <- degree %>%
  filter(`igraph::degree(graph)`==3) %>%
  pull(node) %>%
  unique()

genes_dis_human <- degree %>%
  filter(`igraph::degree(graph)`==2 & !(node %in% GMT_genes$from)) %>%
  pull(node) %>%
  unique()
  
genes_GMT_human <- degree %>%
  filter(`igraph::degree(graph)`==2 & !(node %in% dis_genes$from)) %>%
  pull(node) %>%
  unique()

solved_genes <- c(solved_genes_dis_GMT,solved_genes_dis_GMT_human)



write.csv(all_edges,"all_edges.csv",quote = F,row.names = F)

converted <- alias2official_gene_symbol(genes = genes_to_convert,
                                        manualSearch = T)

all_genes_2 <- unique(c(all_genes_in_GMT,converted$converted))
write.table(all_genes_2,file="all_genes_2.txt",quote = F,row.names = F)

table(all_genes_2 %in% unique(reactome_GMT$gene))

converted_in_GMT <- converted$converted[converted$converted %in% unique(reactome_GMT$gene)]

converted_not_in_GMT <- converted$converted[!converted$converted %in% converted_in_GMT]


mart <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
listFilters(mart)
BM <- getBM(values=list(external_gene_name=all_genes,
                        external_synonym=all_genes),
            attributes = c("external_gene_name","ensembl_gene_id","entrezgene_id"),
            filters = c("external_gene_name","external_synonym"),
            mart = mart)

listAttributes(mart)

conveted_GMT <- alias2official_gene_symbol(gene = unique(reactome_GMT$gene),
                                           manualSearch = T)
