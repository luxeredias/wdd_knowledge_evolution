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

#set working directory
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#load necessary data
data("all_edges_converted")
all_edges_converted <- all_edges_converted %>%
  dplyr::filter(!Source=="")
data("reactome_GMT_converted")
#data("cannonical_GMT_converted")
#data("reactome_full_GMT_converted")
data("all_edges")
data("diseases")

#useful objects
classes <- c("inflammatory","psychiatric","infectious")
color <- c("#ff743e","#00b7da","#68ad36")
inflammatory_diseases <- unique(edges_list[[1]]$Target)
psychiatric_diseases <- unique(edges_list[[2]]$Target)
infectious_diseases <- unique(edges_list[[3]]$Target)

#edit Reactome GMT file
reactome_GMT_converted$reactome_Id <- as.character(str_extract_all(string = reactome_GMT_converted$ont,
                                                      pattern = "R-HSA-\\d+"))
reactome_GMT_converted$ont <- str_remove_all(string = reactome_GMT_converted$ont,
                                             pattern = " Homo sapiens R-HSA-\\d+")

all_genes <- all_edges_converted %>%
  dplyr::select(1,2,31) %>%
  dplyr::rename(docs=3) %>%
  dplyr::filter(docs > 0) %>%
  dplyr::pull(Source) %>%
  unique()

#select the top 9 most connected diseases in 2018
top_9_list <- list()
Reactome_results_all_list <- list()
selected_terms_list <- list()
all_terms_list <- list()
GMT <- reactome_GMT_converted[,-3]
save(GMT,file="data/GMT.RData")
for (i in 1:3){
  #selecte edges of each disease type
  class <- classes[i]
  diss <- diseases$disease[diseases$class==class]
  edge <- all_edges_converted %>%
    dplyr::filter(Target %in% diss)
  degrees <- table(edge$Target)
  top_9 <- as.character(as.data.frame(degrees[order(degrees,decreasing = T)])[1:9,1])
  top_9_list[[i]] <- top_9
  
  #perform functional enrichment against Reactome for each of the top9 disease's genes
  for (j in 1:length(diss)){
    dis <- diss[j]
    #dis <- top_9[j]
    genes <- edge %>%
      dplyr::filter(Target == dis) %>%
      dplyr::pull(Source) %>%
      unique()
    #Reactome enrichment
    tryCatch({reactome_results <- enricher(gene = genes,pvalueCutoff = 0.05,
                                           pAdjustMethod = "none",
                                           TERM2GENE = GMT)
    reactome_results <- reactome_results@result
    reactome_results <- reactome_results %>%
      filter(pvalue < 0.01,Count>0) %>%
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
  all_terms <- unique(reactome_results_all$Term)
  all_terms_list[[i]] <- all_terms
  enrichment_filename <- paste0("data/all_",class,"_reactome_enrichment.csv")
  write.csv(reactome_results_all,file = enrichment_filename,quote = F,row.names = F)
  
  # #plot results
  # data_plot <- reactome_results_all
  # data_plot[is.na(data_plot)] <- 1
  # data_plot <- data_plot %>%
  #   column_to_rownames("Term")
  # data_plot <- -log(data_plot)
  # filename <- paste0("figures/",classes[i],"_Reactome.svg")
  # svg(filename,width = length(diss),height = length(rownames(data_plot)))
  # heatmap(as.matrix(data_plot),
  #         col=pal,scale = "none",margins = c(10,30), #margin = c(col,row)
  #         cexRow = .6,
  #         cexCol = .5)
  # dev.off()
  
  #get the 20 terms that are enriched in the highest number of diseases
  #transform enrichment results in -log10pvalue
  data_plot <- reactome_results_all
  data_plot[is.na(data_plot)] <- 1
  data_plot <- data_plot %>%
    column_to_rownames("Term")
  data_plot <- -log(data_plot)
  
  #enriched terms dataframe
  terms <- data.frame(term=unique(rownames(data_plot)),
                      enriched_in=as.numeric(""))
  
  #calculate the number of diseases with enrichment for each term
  for (t in 1:length(terms$term)){
    term <- terms$term[t]
    df <- data_plot %>%
      filter(rownames(data_plot)==term) %>%
      t()
    diseases_with_enrichment <- table(df>0)["TRUE"]
    terms$enriched_in[t] <- diseases_with_enrichment
  }
  
  #select the top 20 terms enriched in more diseases
  selected_terms <- terms %>%
    arrange(desc(enriched_in)) %>%
    pull(term) %>%
    .[1:20]
  
  selected_terms_list[[i]] <- selected_terms
  
  #subset enrichment table with top 20 terms and top 9 diseases
  data_selected <- data_plot[rownames(data_plot) %in% selected_terms,top_9]
  
  #plot
  filename <- paste0("figures/",classes[i],"_selected_terms.svg")
  pal <- colorRampPalette(c("white",color[i]))(1000)
  svg(filename,width = 5,height = 5)
  heatmap(as.matrix(data_selected[,colSums(data_selected)>0]),
          col=pal,scale = "none",margins = c(8,15), #margin = c(col,row)
          cexRow = .6,
          cexCol = .5)
  dev.off()
  print(paste0(classes[i]," done"))
}

all_selected_terms <- unique(unlist(selected_terms_list))
all_terms <- unique(unlist(all_terms_list))

common_selected_terms <- Reduce(intersect,selected_terms_list)
common_all_terms <- Reduce(intersect,all_terms_list)

# inflammatory_unique_terms <- setdiff(selected_terms_list[[1]],common_terms)
# psichiatric_unique_terms <- setdiff(selected_terms_list[[2]],common_terms)
# infectious_unique_terms <- setdiff(selected_terms_list[[3]],common_terms)

#perform enrichment for all diseases in each year for the selected terms
evolution_reactome_list <- list()
evolution_reactome_melt_list <- list()
for (i in 1:99){
  # dis <- unlist(top_9_list)[i]
  dis <- unique(all_edges_converted$Target)[i]
  df <- all_edges_converted %>%
    filter(Target==dis)
  cols <- colnames(df)[-c(1,2)]
  #perform enrichment for disease in each year
  for (j in 1:29){
    col <- which(colnames(df)==cols[j])
    genes_year <- all_edges_converted %>%
      dplyr::select(1,2,col) %>%
      dplyr::rename(docs=3) %>%
      dplyr::filter(docs > 0) %>%
      dplyr::pull(Source) %>%
      unique()
    
    genes <- df %>%
      dplyr::select(1,2,col) %>%
      dplyr::rename(docs=3) %>%
      dplyr::filter(docs > 0) %>%
      dplyr::pull(Source) %>%
      unique()
    
    #reactome enrichment
    if (length(genes)>0){
        reactome_results <- enricher(gene = genes,pvalueCutoff = 0.05,
                                     pAdjustMethod = "none",
                                     TERM2GENE = GMT)
        if (length(reactome_results)==0){
          reactome_results <- data.frame(Term="",x=as.numeric(""))
          colnames(reactome_results) <- c("Term",dis)
          if (j == 1){
            reactome_results_all <- reactome_results
          }else{
            reactome_results_all <- reactome_results_all %>% 
              full_join(reactome_results,by="Term")
          }
        }else{
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
          }
        }
      
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
  }#end of enrichment loop for each disease in the given year
  
  reactome_results_all <- reactome_results_all[-1,]
  reactome_results_all[,-1][is.na(reactome_results_all[,-1])] <- 1
  reactome_results_all[,-1] <- -log(reactome_results_all[,-1])
  colnames(reactome_results_all) <- c("Term",seq(1990,2018,1))
  evolution_reactome_list[[i]] <- reactome_results_all
  
  reactome_results_melt <- melt(reactome_results_all)
  if (length(reactome_results_all$Term)==0){
    evolution_reactome_melt_list[[i]] <- reactome_results_melt
  }else{
    reactome_results_melt$dis <- dis
    evolution_reactome_melt_list[[i]] <- reactome_results_melt
  }
}

#names(evolution_reactome_list) <- unlist(top_9_list)
names(evolution_reactome_list) <- unique(all_edges_converted$Target)
evolution_reactome_melt <- bind_rows(evolution_reactome_melt_list)

evolution_reactome_melt$class <- ifelse(evolution_reactome_melt$dis %in% inflammatory_diseases,yes = "inflammatory",no="")
evolution_reactome_melt$class <- ifelse(evolution_reactome_melt$dis %in% psychiatric_diseases,yes = "psychiatric",no=evolution_reactome_melt$class)
evolution_reactome_melt$class <- ifelse(evolution_reactome_melt$dis %in% infectious_diseases,yes = "infectious",no=evolution_reactome_melt$class)

evolution_reactome_cast <- dcast(evolution_reactome_melt,formula = Term~variable+dis)
evolution_reactome_cast[is.na(evolution_reactome_cast)] <- 0

evolution_reactome_melt_2 <- melt(evolution_reactome_cast)
evolution_reactome_melt_2$year <- as.numeric(substr(evolution_reactome_melt_2$variable,start = 1,stop = 4))
evolution_reactome_melt_2$variable <- str_sub(evolution_reactome_melt_2$variable,start = 6,end = str_length(evolution_reactome_melt_2$variable))
evolution_reactome_melt_2 <- evolution_reactome_melt_2[,c(1,4,3,2)]
colnames(evolution_reactome_melt_2) <- c("Term","year","value","dis")

evolution_reactome_melt_2$class <- ifelse(evolution_reactome_melt_2$dis %in% inflammatory_diseases,yes = "inflammatory",no="")
evolution_reactome_melt_2$class <- ifelse(evolution_reactome_melt_2$dis %in% psychiatric_diseases,yes = "psychiatric",no=evolution_reactome_melt_2$class)
evolution_reactome_melt_2$class <- ifelse(evolution_reactome_melt_2$dis %in% infectious_diseases,yes = "infectious",no=evolution_reactome_melt_2$class)

top_9_all <- unlist(top_9_list)
evolution_reactome_melt <- evolution_reactome_melt_2

#plot evolution of all terms for all diseases
pdf("figures/all_terms_evolution.pdf")
for (i in 1:length(unique(evolution_reactome_melt$Term))){
  term <- unique(evolution_reactome_melt$Term)[i]
  p<-evolution_reactome_melt %>%
    filter(Term==term) %>%
    ggplot(aes(x=year,y=dis,fill=class,alpha=log(value+1)))+
      geom_tile(show.legend = F)+
      theme_minimal(base_line_size = .1)+
      theme(axis.title.y = element_blank(),
            #axis.text.y = element_blank(),
            axis.text.x = element_blank())+
      scale_x_discrete(limits = c(1990,2000,2010,2018),name="Year")+
      facet_grid(rows = vars(class),cols = vars(Term), scales = "free", space = "free")+
      #facet_wrap(facets = ~term,scales = "free")+
      scale_fill_manual(values = c(color[3],color[1],color[2]))+
      theme(strip.text = element_text(size = 8))
  print(p)
}
dev.off()

#plot evolution of common terms for all diseases
pdf("figures/common_terms_all_dis_evolution.pdf")
for (i in 1:length(common_all_terms)){
  term <- common_all_terms[i]
  p<-evolution_reactome_melt %>%
    filter(Term==term) %>%
    ggplot(aes(x=year,y=dis,fill=class,alpha=log(value+1)))+
    geom_tile(show.legend = F)+
    theme_minimal(base_line_size = .1)+
    theme(axis.title.y = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_blank())+
    scale_x_discrete(limits = c(1990,2000,2010,2018),name="Year")+
    facet_grid(rows = vars(class),cols = vars(Term), scales = "free", space = "free")+
    #facet_wrap(facets = ~term,scales = "free")+
    scale_fill_manual(values = c(color[3],color[1],color[2]))+
    theme(strip.text = element_text(size = 8))
  print(p)
}  
dev.off()

#plot evolution of common terms for top 9 diseases in each category
pdf("figures/common_terms_top_9_dis_evolution.pdf")
for (i in 1:length(common_all_terms)){
  term <- common_all_terms[i]
  p<-evolution_reactome_melt %>%
    filter(Term==term & dis %in% top_9_all) %>%
    ggplot(aes(x=year,y=dis,fill=class,alpha=log(value+1)))+
    geom_tile(show.legend = F)+
    theme_minimal(base_line_size = .1)+
    theme(axis.title.y = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_blank())+
    scale_x_discrete(limits = c(1990,2000,2010,2018),name="Year")+
    facet_grid(rows = vars(class),cols = vars(Term), scales = "free", space = "free")+
    #facet_wrap(facets = ~term,scales = "free")+
    scale_fill_manual(values = c(color[3],color[1],color[2]))+
    theme(strip.text = element_text(size = 8))
  print(p)
}  
dev.off()





