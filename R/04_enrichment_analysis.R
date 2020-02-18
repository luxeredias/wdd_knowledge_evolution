#script to perform gene set enrichment analysis of the genes associated to
#infectious, inflammatory and psychiatric disorders from 1990 to 2018

#libraries used
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(reshape2)
library(tibble)
options(stringsAsFactors = F)

library(tidyr)
library(purrr)
library(patchwork)
library(msigdbr)

library(readr)
library(igraph)


library(gplots)

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
reactome_GMT_converted <- read.gmt("data/ReactomePathwaysLevel3.gmt")

picked_terms <- read.csv("data/picked_terms_reactome.txt",col.names = "term")
reactome_GMT_selected <- reactome_GMT_converted %>%
  dplyr::filter(ont %in% picked_terms$term) %>%
  dplyr::select(1,2)

all_edges_converted <- all_edges_converted %>%
  dplyr::filter(!Source=="")

all_genes <- all_edges_converted %>%
  dplyr::select(1,2,31) %>%
  dplyr::rename(docs=3) %>%
  dplyr::filter(docs > 0) %>%
  dplyr::pull(Source) %>%
  unique()

#select the top 9 most connected diseases in 2018
top_9_all <- c()
Reactome_results_all_list <- list()
selected_terms_list <- list()
selected_terms_years_list <- list()
for (i in 1:3){
  #selecte edges of each disease type
  class <- classes[i]
  diss <- diseases$disease[diseases$class==class]
  edge <- all_edges_converted %>%
    dplyr::filter(Target %in% diss)
  degrees <- table(edge$Target)
  top_9 <- as.character(as.data.frame(degrees[order(degrees,decreasing = T)])[1:9,1])
  top_9_all <- c(top_9_all,top_9)
  
  #perform functional enrichment against Reactome for each of the top9 disease's genes
  for (j in 1:length(top_9)){
    dis <- top_9[j]
    genes <- edge %>%
      dplyr::filter(Target == dis) %>%
      dplyr::pull(Source) %>%
      unique()
    #Reactome enrichment
    tryCatch({reactome_results <- enricher(gene = genes,pvalueCutoff = 0.05,
                                           pAdjustMethod = "none",
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
      reactome_results <- data.frame(Term="",dis=as.character(""))
      if (j == 1){
        reactome_results_all <- reactome_results
      }else{
        reactome_results_all <- reactome_results_all %>% 
          full_join(reactome_results,by="Term")
      }})
  }
  
  #save aggregated results to Reactome_results_all list
  Reactome_results_all_list[[i]] <- reactome_results_all
  
  #plot results
  data_plot <- reactome_results_all
  data_plot[is.na(data_plot)] <- 1
  data_plot <- data_plot %>%
    column_to_rownames("Term")
  data_plot <- -log(data_plot)
  pal <- colorRampPalette(c("white",color[i]))(1000)
  filename <- paste0("figures/",classes[i],"_Reactome.svg")
  svg(filename,width = 10,height = 10)
  heatmap(as.matrix(data_plot),
          col=pal,scale = "none",margins = c(10,30), #margin = c(col,row)
          cexRow = .6,
          cexCol = .5)
  dev.off()
  
  #select terms to be used in evolution analysis: terms that were enriched
  #in the majority of the top 9 diseases in each category
  
  #get the 20 terms that are enriched in the highest number of diseases
  terms <- data.frame(term=unique(rownames(data_plot)),
                      enriched_in=as.numeric(""))
  
  for (t in 1:length(terms$term)){
    term <- terms$term[t]
    df <- data_plot %>%
      filter(rownames(data_plot)==term) %>%
      t()
    diseases_with_enrichment <- table(df>0)["TRUE"]
    terms$enriched_in[t] <- diseases_with_enrichment
  }
  
  selected_terms <- terms %>%
    arrange(desc(enriched_in)) %>%
    pull(term) %>%
    .[1:20]
  
  selected_terms_list[[i]] <- selected_terms
  
  data_selected <- data_plot[rownames(data_plot) %in% selected_terms,]
  filename <- paste0("figures/",classes[i],"_selected_terms.svg")
  svg(filename,width = 5,height = 5)
  heatmap(as.matrix(data_selected[,colSums(data_selected)>0]),
          col=pal,scale = "none",margins = c(8,15), #margin = c(col,row)
          cexRow = .6,
          cexCol = .5)
  dev.off()
  print(paste0(classes[i]," done"))
}

all_selected_terms <- unique(unlist(selected_terms_list))

common_terms <- Reduce(intersect,selected_terms_list)

inflammatory_unique_terms <- setdiff(selected_terms_list[[1]],common_terms)
psichiatric_unique_terms <- setdiff(selected_terms_list[[2]],common_terms)
infectious_unique_terms <- setdiff(selected_terms_list[[3]],common_terms)


#perform enrichment for all diseases in each year for the selected terms
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





