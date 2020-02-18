#script to perform gene set enrichment analysis of the genes associated to
#infectious, inflammatory and psychiatric disorders from 1990 to 2018

#libraries used
library(stringr)
library(dplyr)
library(clusterProfiler)
library(ggplot2)

library(tidyr)
library(purrr)
library(patchwork)
library(msigdbr)

library(tibble)

library(readr)
library(igraph)


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
  dplyr::filter(ont %in% selected_terms$term) %>%
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
  #save aggregated results to results_all lists
  Reactome_results_all_list[[i]] <- reactome_results_all
  print(paste(dis_type[i],"done"))
  #plot reactome enrichment for each disease class
  pal <- colorRampPalette(c("white",color[i]))(1000)
  data <- reactome_results_all
  data[is.na(data)] <- 1
  rownames(data) <- data$Term
  data$Term <- NULL
  data <- -log(data)
  filename <- paste0("figures/",dis_type[i],"_reactome.svg")
  svg(filename,width = 10,height = 10)
  heatmap(as.matrix(data[-1,]),
          col=pal,scale = "none",margins = c(15,40), #margin = c(col,row)
          cexRow = .6,
          cexCol = .5)
  dev.off()
}

selected_terms_list <- list()
selected_terms_years_list <- list()
for (i in 1:3){
  data <- Reactome_results_all_list[[i]]
  data[is.na(data)] <- 1
  #data$Term <- gsub(x = data$Term,pattern = "KEGG_",replacement = "")
  #data$Term <- gsub(x = data$Term,pattern = "_",replacement = " ")
  rownames(data) <- data$Term
  data$Term <- NULL
  data <- -log(data)
  pal <- colorRampPalette(c("white",color[i]))(1000)
  filename <- paste0("~/Área de Trabalho/Papers_Helder/Evolution_all/Figures/Figure_02/Version_02/",dis_type[i],"_Reactome_v2.0.svg")
  svg(filename,width = 10,height = 10)
  heatmap(as.matrix(data),
          col=pal,scale = "none",margins = c(10,30), #margin = c(col,row)
          cexRow = .6,
          cexCol = .5)
  dev.off()
  PA <- data
  PA[PA>0] <- 1
  nterms <- c()
  for (g in 1:ncol(PA)){
    selected_terms <- rownames(PA[rowSums(PA)>g,])
    nterms <- c(nterms,length(selected_terms))
  }
  nterms <- abs(10-nterms)
  g<-which(nterms==min(nterms))
  selected_terms <- c(rownames(PA[rowSums(PA)>g,]),picked_terms)  
  selected_terms_list[[i]] <- selected_terms
  data_selected <- data[rownames(data) %in% selected_terms,]
  filename <- paste0("~/Área de Trabalho/Papers_Helder/Evolution_all/Figures/Figure_02/Version_02/",dis_type[i],"_selected_terms.svg")
  svg(filename,width = 5,height = 5)
  heatmap(as.matrix(data_selected[,colSums(data_selected)>0]),
          col=pal,scale = "none",margins = c(8,15), #margin = c(col,row)
          cexRow = .6,
          cexCol = .5)
  dev.off()
}


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





