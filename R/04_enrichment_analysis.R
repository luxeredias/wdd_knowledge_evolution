#script to perform gene set enrichment analysis of the genes associated to
#infectious, inflammatory and psychiatric disorders from 1990 to 2018

#libraries used
library(tidyr)
library(purrr)
library(patchwork)
library(msigdbr)
library(clusterProfiler)
library(tibble)
library(dplyr)
library(readr)
library(igraph)
library(ggplot2)
library(stringr)
library(gplots)
library(reshape2)
library(lattice)
library(corrplot)
library(ggridges)
library(enrichR)
library(rvest)
library(httr)
library(heatmaply)
library(VennDiagram)
options(stringsAsFactors = F)

#set working directory
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#load necessary data
data("all_edges_converted")
data("reactome_GMT_converted")
data("diseases")

#useful objects
classes <- c("inflammatory","psychiatric","infectious")
color <- c("#ff743e","#00b7da","#68ad36")

#edit Reactome GMT file
reactome_GMT_converted$reactome_Id <- as.character(str_extract_all(string = reactome_GMT_converted$ont,
                                                      pattern = "R-HSA-\\d+"))
reactome_GMT_converted$ont <- str_remove_all(string = reactome_GMT_converted$ont,
                                             pattern = " Homo sapiens R-HSA-\\d+")
selected_terms <- read.csv("data/picked_terms_reactome.txt",col.names = "term")
reactome_GMT_selected <- reactome_GMT_converted %>%
  filter(ont %in% selected_terms$term) %>%
  select(1,2)

all_genes <- all_edges_converted %>%
  select(1,2,31) %>%
  rename(docs=3) %>%
  filter(docs > 0) %>%
  pull(Source) %>%
  unique()

#perform enrichment for all diseases in each year
for (i in 1:27){
  dis <- top_9_all[i]
  df <- all_edges_converted %>%
    filter(Target==dis)
  cols <- colnames(df)[-c(1,2)]
  for (j in 1:29){
    col <- which(colnames(df)==cols[j])
    
    genes <- df %>%
      select(1,2,31) %>%
      rename(docs=3) %>%
      filter(docs > 0) %>%
      pull(Source) %>%
      unique()
    #reactome enrichment
    if (length(genes)>0){
      tryCatch({
        reactome_results <- enricher(gene = genes,pvalueCutoff = 0.05,
                                     pAdjustMethod = "BH",
                                     universe = all_genes,
                                     TERM2GENE = reactome_GMT_converted[,-3])
        reactome_results <- reactome_results@result
        reactome_results <- reactome_results %>%
          filter(pvalue < 0.05,Count>0) %>%
          dplyr::select(ID,p.adjust)
        colnames(reactome_results) <- c("Term",dis)
        if (j == 1){
          reactome_results_all <- reactome_results
        }else{
          reactome_results_all <- reactome_results_all %>% 
            full_join(reactome_results,by="Term")
        }},error=function(err){
          reactome_results <- data.frame(Term="",x=as.numeric(""))
          colnames(reactome_results) <- c("Term",dis)
          if (j == 1){
            reactome_results_all <- reactome_results
          }else{
            reactome_results_all <- reactome_results_all %>% 
              full_join(reactome_results,by="Term")
          }
        })
    }else{
      reactome_results <- data.frame(Term="",x=as.numeric(""))
      colnames(reactome_results) <- c("Term",dis)
      if (j == 1){
        reactome_results_all <- reactome_results
      }else{
        reactome_results_all <- reactome_results_all %>% 
          full_join(reactome_results,by="Term")
      }
    }
  }
  }#end of enrichment loop for each disease in the given year
  rownames(reactome_results_all) <- reactome_results_all$Term
  reactome_results_all$Term <- NULL
  reactome_results_all[is.na(reactome_results_all)] <- 1
  reactome_results_all <- -log(reactome_results_all)
  
  evolution_reactome_list[[k]] <- reactome_results_all
}

#select the top 9 most connected diseases in 2018
top_9_all <- c()
for (i in 1:3){
  class <- classes[i]
  dis <- diseases$disease[diseases$class==class]
  df <- all_edges_converted %>%
    filter(Target %in% dis)
  degrees <- table(df$Target)
  top_9 <- as.character(as.data.frame(degrees[order(degrees,decreasing = T)])[1:9,1])
  top_9_all <- c(top_9_all,top_9)
}



