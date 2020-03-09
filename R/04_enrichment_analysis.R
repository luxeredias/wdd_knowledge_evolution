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
    .[1:10]
  
  selected_terms_list[[i]] <- selected_terms
  
  #subset enrichment table with top 20 terms and top 9 diseases
  data_selected <- data_plot %>%
    rownames_to_column("Term") %>%
    filter(Term %in% selected_terms) %>%
    dplyr::select(c(1,which(colnames(.) %in% top_9)))
  
  #calculate distance between terms and diseases
  dist_r <- dist(data_selected %>% column_to_rownames("Term"))
  clust_r <- hclust(dist_r)
  order_r <- clust_r$order
  
  dist_c <- dist(t(data_selected %>% column_to_rownames("Term")))
  clust_c <- hclust(dist_c)
  order_c <- clust_r$order
  
  #plot heatmap of selected terms in top9 diseases
  data_selected <- melt(data_selected)
  data_selected$Term <- factor(data_selected$Term,levels = unique(data_selected$Term)[rev(order_r)])
  data_selected$variable <- factor(data_selected$variable,levels = unique(data_selected$variable)[rev(order_c)])
  
  filename <- paste0("figures/",classes[i],"_selected_terms.svg")
  svg(filename,width = 8,height = 5)
  p <- ggplot(data_selected,aes(x=variable,y=Term,fill=color[i],alpha=value))+
    geom_tile(show.legend = F)+
    theme_minimal(base_line_size = .1)+
    theme(axis.text.x = element_text(angle = 45, hjust = 0,size = 7.5),
          axis.text.y = element_text(size = 7.5),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.margin=unit(c(t=1,r=4,b=1,l=1),"cm"))+
    scale_fill_manual(values = color[i])+
    scale_x_discrete(position = "top")
  print(p)
  dev.off()
  
  # #plot
  # filename <- paste0("figures/",classes[i],"_selected_terms.svg")
  # pal <- colorRampPalette(c("white",color[i]))(1000)
  # svg(filename,width = 5,height = 5)
  # heatmap(as.matrix(data_selected[,colSums(data_selected)>0]),
  #         col=pal,scale = "none",margins = c(8,15), #margin = c(col,row)
  #         cexRow = .6,
  #         cexCol = .5)
  # dev.off()
  print(paste0(classes[i]," done"))
}

all_selected_terms <- unique(unlist(selected_terms_list))
all_terms <- unique(unlist(all_terms_list))

common_selected_terms <- Reduce(intersect,selected_terms_list)
common_all_terms <- Reduce(intersect,all_terms_list)

#plot all common terms in top9 diseases in one heatmap
for (i in 1:3){
  reactome_results_common_terms <- Reactome_results_all_list[[i]] %>%
    dplyr::filter(Term %in% common_all_terms) %>%
    dplyr::select(c(1,which(colnames(.) %in% top_9_list[[i]]))) %>%
    column_to_rownames("Term")
  if (i==1){
    reactome_results_common_terms_all <- reactome_results_common_terms
  }else{
    reactome_results_common_terms_all <- cbind(reactome_results_common_terms_all,
                                               reactome_results_common_terms)
  }
}

data_plot <- reactome_results_common_terms_all
data_plot[is.na(data_plot)] <- 1

#determine distance between terms to plot in dendrogram order
dist <- dist(data_plot)
clust <- hclust(dist)
order <- clust$order

#create melt table to plot in ggplot2
data_plot <- -log(data_plot)
data_plot <- data_plot %>%
  rownames_to_column("Term") %>%
  melt()
data_plot$class <- c(rep("inflammatory",288),
                     rep("psychiatric",288),
                     rep("infectious",288))

data_plot$Term <- factor(data_plot$Term,levels = data_plot$Term[order])

#plot heatmap with all common terms in top9 diseases of each category
svg("figures/reactome_common_terms_top9_diseaes.svg",width = 14,height = 7)
p <- ggplot(data_plot,aes(x=variable,y=Term,fill=class,alpha=value))+
  geom_tile(show.legend = F)+
  theme_minimal(base_line_size = .1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0,size = 7.5),
        axis.text.y = element_text(size = 7.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin=unit(c(t=1,r=4,b=1,l=1),"cm"))+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = c(color[3],color[1],color[2]))
print(p)
dev.off()

#plot mostly enriched common terms in top9 diseases in each category
svg("figures/reactome_selected_common_terms_top9_diseaes.svg",width = 12,height = 5)
selected_terms <- clust$labels[rev(order)][1:10]
p <- data_plot %>%
  filter(Term %in% selected_terms) %>%
  ggplot(aes(x=variable,y=Term,fill=class,alpha=value))+
  geom_tile(show.legend = F)+
  theme_minimal(base_line_size = .1)+
  theme(axis.text.x = element_text(angle = 45, hjust = 0,size = 7.5),
        axis.text.y = element_text(size = 7.5),
        axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        plot.margin=unit(c(t=1,r=4,b=1,l=1),"cm"))+
  scale_x_discrete(position = "top")+
  scale_fill_manual(values = c(color[3],color[1],color[2]))
print(p)
dev.off()

evolution_terms <- selected_terms[c(1,3,4,5,8)]

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
    geom_tile(show.legend = T)+
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

#plot evolution of selected common terms for top9 diseases in one image
evolution_terms <- all_terms[c(which(all_terms %in% evolution_terms),
                               14,68,101,24,67,34,28,18,6)]

selected_evolution_terms <- evolution_terms[c(1,3,22,5,11,10)]

svg("figures/selected_terms_top9_dis_evolution.svg",width = 20,height = 5)
p <- evolution_reactome_melt %>%
  filter(dis %in% top_9_all & Term %in% selected_evolution_terms) %>%
  ggplot(aes(x=year,y=dis,fill=class,alpha=value))+
  geom_tile(show.legend = F)+
  theme_minimal(base_line_size = .2)+
  theme(axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank()
        )+
  scale_x_discrete(limits = c(1990,2000,2010,2018),name="Year (1990-2018)")+
  facet_grid(rows = vars(class),cols = vars(Term), scales = "free", space = "free")+
  scale_fill_manual(values = c(color[3],color[1],color[2]))+
  theme(strip.text = element_text(size = 8))
print(p)
dev.off()

#plot only legend of each disease class heatmap
for (i in 1:3){
  file <- paste0("figures/legend_",i,".svg")
  svg(filename = file,
      width = 1,height = 1.7)
  
  p<- evolution_reactome_melt %>%
    filter(dis %in% top_9_all & 
             Term %in% selected_evolution_terms) %>%
    ggplot(aes(x=1,y=value,color=value))+
    geom_point()+
    scale_y_discrete(limits=c(0,max(x$value)))+
    scale_color_gradient(low="white",
                         high=color[i])
  legend <- get_legend(p)
  plot(legend)
  dev.off()
}  

#calculate the percentage of diseases enriched for each selected term
#in each year (separated by category of disease)
for (i in 1:3){
  class <- classes[i]
  nsig_all_years_percent <- as.data.frame(matrix(nrow = 6,ncol = 29))
  colnames(nsig_all_years_percent) <- seq(1990,2018,1)
  rownames(nsig_all_years_percent) <- selected_evolution_terms
  nsig_all_years_absolute <- nsig_all_years_percent
  for (k in 1:6){
    nsigs <- c()
    termo <- selected_evolution_terms[k]
    for (j in 1:29){
      year <- seq(1990,2018,1)[j]
      nsigs_year <- evolution_reactome_melt %>%
        filter(year==!!year &
                 Term==termo &
                 class==!!class &
                 value > 4.6)
      nsigs_year <- length(unique(nsigs_year$dis))
      nsig_all_years_absolute[k,j] <- nsigs_year
      nsigs_year <- nsigs_year/length(unique(edges_list[[i]]$Target))*100
      nsig_all_years_percent[k,j] <- nsigs_year
    }
  }
  nsig_all_years_percent <- nsig_all_years_percent %>%
    rownames_to_column("Term")
  nsig_all_years_percent_melt <- melt(nsig_all_years_percent)
  nsig_all_years_percent_melt$class <- class
  
  nsig_all_years_absolute <- nsig_all_years_absolute %>%
    rownames_to_column("Term")
  nsig_all_years_absolute_melt <- melt(nsig_all_years_absolute)
  nsig_all_years_absolute_melt$class <- class
  
  if (i==1){
    nsig_all_years_all_percent <- nsig_all_years_percent_melt
    nsig_all_years_all_absolute <- nsig_all_years_absolute_melt
  }else{
    nsig_all_years_all_percent <- as.data.frame(rbind(nsig_all_years_all_percent,nsig_all_years_percent_melt))
    nsig_all_years_all_absolute <- as.data.frame(rbind(nsig_all_years_all_absolute,nsig_all_years_absolute_melt))
  }
}

nsig_all_years_all_percent$variable <- as.numeric(as.character(nsig_all_years_all_percent$variable))
nsig_all_years_all_absolute$variable <- as.numeric(as.character(nsig_all_years_all_absolute$variable))

#plot line graphs of the percent of diseases that present enrichment
#for each of the selected terms
svg("figures/nsig_years_percent_selected_terms_top9.svg",width = 20,height = 5)
p<-ggplot(nsig_all_years_all_percent,aes(x=variable,y=value,color=class))+
  geom_line(show.legend = T,size=1)+
  scale_color_manual(values = c(color[3],color[1],color[2]))+
  scale_y_continuous(limits=c(0,100))+
  facet_wrap(facets = ~Term,nrow = 1,ncol=6)+
  theme_minimal()
print(p)
dev.off()

#plot line graphs of the absolute number of diseases that present enrichment
#for each of the selected terms
svg("figures/nsig_years_absolute_selected_terms_top9.svg",width = 20,height = 5)
p<-ggplot(nsig_all_years_all_absolute,aes(x=variable,y=value,color=class))+
  geom_line(show.legend = T,size=1)+
  scale_color_manual(values = c(color[3],color[1],color[2]))+
  scale_y_continuous(limits=c(0,20))+
  facet_wrap(facets = ~Term,nrow = 1,ncol=6)+
  theme_minimal()
print(p)
dev.off()

#plot all terms in top9 dis evolution (similar to Fig. 2)
svg("figures/all_common_terms_top9_dis_evolution.svg",width = 100,height = 5)
p <- evolution_reactome_melt %>%
  filter(dis %in% top_9_all & Term %in% common_all_terms) %>%
  ggplot(aes(x=year,y=dis,fill=class,alpha=log2(value)))+
  geom_tile(show.legend = F)+
  theme_minimal(base_line_size = .2)+
  theme(axis.title.y = element_blank(),
        #axis.text.y = element_blank(),
        axis.text.x = element_blank()
  )+
  scale_x_discrete(limits = c(1990,2000,2010,2018),name="Year (1990-2018)")+
  facet_grid(rows = vars(class),cols = vars(Term), scales = "free", space = "free")+
  scale_fill_manual(values = c(color[3],color[1],color[2]))+
  theme(strip.text = element_text(size = 8))
print(p)
dev.off()

#plot all terms in top9 dis evolution in groups of 4
intervals <- list(c(1:4),c(5:8),c(9:12),c(13:16),
               c(17:20),c(21:24),c(25:28),c(29:32))

for (i in 1:8){
  file <- paste0("figures/all_common_terms_top9_dis_evolution_",i,".svg")
  int <- intervals[[i]]
  svg(file,width = 12,height = 5)
  terms <- common_all_terms[int]
  p <- evolution_reactome_melt %>%
    filter(dis %in% top_9_all & Term %in% terms) %>%
    ggplot(aes(x=year,y=dis,fill=class,alpha=value))+
    geom_tile(show.legend = F)+
    theme_minimal(base_line_size = .2)+
    theme(axis.title.y = element_blank(),
          #axis.text.y = element_blank(),
          axis.text.x = element_blank()
    )+
    scale_x_discrete(limits = c(1990,2000,2010,2018),name="Year (1990-2018)")+
    facet_grid(rows = vars(class),cols = vars(Term), scales = "free", space = "free")+
    scale_fill_manual(values = c(color[3],color[1],color[2]))+
    theme(strip.text = element_text(size = 8))
  print(p)
  dev.off()
}  
