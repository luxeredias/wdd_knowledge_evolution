#script to convert all gene names in WDD files to official gene symbol
#and all gene names in Reactome GMT file to official gene symbol

#packages
library(dplyr)
source("~/Área de Trabalho/GitHub_projects/alias2official/scripts/alias2official_gene_symbol.R")

#load edges file
load("edges.RData")

all_genes <- c()
for (i in 1:3){
  genes <- unique(edges_list[[i]]$Source)
  all_genes <- unique(c(all_genes,genes))
}

reactome_GMT <- read.gmt("GMT_files/all_filtered/Filtered_Reactome_2016.txt") 

all_genes_in_GMT <- all_genes[all_genes %in% unique(reactome_GMT$gene)]

genes_to_convert <- all_genes[!all_genes %in% unique(reactome_GMT$gene)]

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
