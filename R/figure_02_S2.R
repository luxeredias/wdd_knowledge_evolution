#With this scritp, you can reproduce the analysis and plot the images that
#compose panel A-C in figure 02 of the paper

#load necessary packages
#file management
library(data.table)
#data manipulation
library(dplyr)
library(tidyverse)
library(reshape2)
library(stringr)
#network/gene analysis
library(igraph)
library(clusterProfiler)
#literature scraping
library(easyPubMed)
library(MeSH.db)
#plotting
library(ggplot2)
library(ggrepel)
library(ggraph)
library(ggridges)

options(stringsAsFactors = F)

#load data and create necessary objects
data("all_edges")
all_edges <- do.call(rbind,edges_list)
dis_dis_nodes <- as.data.frame(fread("data/dis_dis_nodes.csv"))
dis_class_df <- dis_dis_nodes %>%
  dplyr::select(1,3) %>%
  dplyr::rename(disease=1,class=2)
diseases <- unique(all_edges$Target)
genes <- unique(all_edges$Source)
years <- seq(1990,2018)

#calculate dis-dis similarity between diseases of combinations of classes
#in each year
#Method: run enricher() function from clusterProfiler package to obtain
#the pvalue associated with the gene sharing between disease A and B
#in each year. Disease A: source of "genes". Diases B: in gmt file
#obs: gmt file for genes associated with each disease in each year created
#on the fly for each year
diseases <- unique(all_edges$Target)
all_genes <- unique(all_edges$Source)
years <- seq(1990,2018,1)
for(i in 1:length(years)){
  yr <- years[i]
  col <- which(grepl(yr,colnames(all_edges)))
  df <- all_edges %>%
    dplyr::select(1,2,col) %>%
    rename(docs=3) %>%
    filter(docs!=0)
  gmt <- df[,c(2,1)] %>% rename(dis=1,gene=2) %>%
    arrange(dis)
  gene_list <- split(df$Source,df$Target)
  enr <- lapply(gene_list,enricher,
                pAdjustMethod = "BH",
                universe = all_genes,
                TERM2GENE = gmt)
  enr <- enr[enr!="NULL"]
  enr_out <- do.call(rbind,lapply(enr,function(x){x<-x@result})) %>%
    rownames_to_column("edge") %>%
    separate(edge,sep = "\\.",into = c("Source","Target")) %>%
    mutate(edge=paste0(Source,"_",ID)) %>%
    dplyr::select(1,3,12,7,8,10,11) %>%
    dplyr::mutate(year=yr)
  if(i==1){
    enr_all_out <- enr_out
  }else{
    enr_all_out <- rbind(enr_all_out,enr_out)
  }
}

saveRDS(enr_all_out,file = "intermediate/dis_dis_similarity_clusterProfiler.RDS")
enr_all_out <- readRDS("intermediate/dis_dis_similarity_clusterProfiler.RDS")

#annotate comparisons acording to disease classes
enr_all_out_2 <- enr_all_out %>%
  filter(Source!=ID) %>%
  mutate(logpadj=-log(p.adjust),
         logpval=-log(pvalue)) %>%
  rename(Target=ID) %>%
  left_join(dis_class_df,by=c("Source"="disease")) %>%
  rename(class_source=class) %>%
  left_join(dis_class_df,by=c("Target"="disease")) %>%
  rename(class_target=class) %>%
  mutate(comparison=ifelse(class_target!=class_source,
                           "between_classes","within_classes")) %>%
  mutate(between=paste0(class_source,"_vs_",class_target))

####figure S2 - line graphs of dis-dis similarity through years####
#select the top5 most similar diseases to each of the 9 most connected
#diseases in each cateogry (in the 2018 network)
top_9_2018 <- enr_all_out_2 %>%
  filter(year==2018) %>%
  filter(Source %in% top_9) %>%
  group_by(Source,between) %>%
  top_n(n = 5,wt = logpval)

#plot line graphs for dis-dis similarity through year for all category pairs
selected_between <- unique(enr_all_out_2$between)[c(4,8,2,7,5,6)]
for(i in 1:length(selected_between)){
  sel <- selected_between[i]
  p <- enr_all_out_2 %>%
    filter(edge %in% top_9_2018$edge,
           class_source!=class_target,
           between==sel) %>%
    ggplot(aes(x=year,y=logpval,color=edge))+
    geom_line(size=1)+
    scale_x_continuous(breaks = c(1990,2000,2010,2018))+
    geom_text_repel(data = enr_all_out_2 %>%
                      filter(edge %in% top_9_2018$edge,
                             class_source!=class_target,
                             between==sel,
                             year==2018),
                    aes(label=Target),
                    nudge_x = -10,nudge_y = 5,
                    color="black",size=2.5,segment.size = .1)+
    facet_wrap(facets = ~Source)+
    theme_minimal()+
    ylab("similarity")+xlab("year")+
    theme(legend.position = "none")
  fil <- paste0("figures/figure_S2/",sel,"_similarity_evolution.pdf")
  pdf(file = fil,width = 13.70,height = 10.15)
  print(p)
  dev.off()
}

####figure 2A - heatmaps of dis-dis similarity between categories in 2018####
#plot heatmaps of dis-dis similarity for:
#psychiatric vs inflammatory
#psychiatric vs infectious
#inflammatory vs infectious
myColor_list <- list(psy_inflam=colorRampPalette(colors = c("white","pink","red3"))(n=1000),
                     psy_infect=colorRampPalette(colors = c("white","lightgreen","darkgreen"))(n=1000),
                     inflam_infect=colorRampPalette(colors = c("white","plum1","purple4"))(n=1000))

folders <- list.dirs(path = "figures/figure_02/panel_A/",recursive = F,full.names = F)
folders <- folders[c(3,2,1)]
for(i in 1:length(selected_between[c(1,2,4)])){
  sel <- selected_between[c(1,2,4)][i]
  df_sel <- enr_all_out_2 %>%
    filter(edge %in% top_9_2018$edge,
           class_source!=class_target,
           between==sel,
           year==2018)
  
  df_2018 <- enr_all_out_2 %>%
    filter(Source %in% unique(df_sel$Source),
           Target %in% unique(df_sel$Target),
           class_source!=class_target,
           between==sel,
           year==2018)
  
  mx_2018 <- reshape2::dcast(df_2018,
                             formula = Source~Target,
                             fun.aggregate = sum,
                             value.var = "logpval") %>%
    column_to_rownames("Source") %>% as.matrix()
  
  max <- max(mx_2018)
  
  rownames(mx_2018) <- str_to_title(str_to_lower(rownames(mx_2018)))
  colnames(mx_2018) <- str_to_title(str_to_lower(colnames(mx_2018)))
  
  myColor <- myColor_list[[i]]
  fil_2018 <- paste0("figures/figure_02/panel_A/dis_dis_similarity_2018_",sel,".pdf")
  pdf(fil_2018,width = 7.8,height = 4.5)
  p <- pheatmap::pheatmap(mx_2018,
                        cluster_rows = T,
                        cluster_cols = T,
                        color = myColor,breaks = seq(0,max,max/1000))
  dev.off()
  
  col_order <- p$tree_col$order
  row_order <- p$tree_row$order
  
  years_plot <- seq(1990,2018,1)
  for(j in 1:length(years_plot)){
    yr <- years_plot[j]
    good_edges <- top_9_2018 %>%
      filter(between==sel) %>%
      pull(edge)
    
    plot_mx <- enr_all_out_2 %>%
      filter(year==yr,
             edge %in% top_9_2018$edge,
             between==sel,
             class_source!=class_target)
    
    plot_mx <- enr_all_out_2 %>%
      filter(year==yr,
             Source %in% plot_mx$Source,
             Target %in% plot_mx$Target,
             between==sel,
             class_source!=class_target)
    
    edges_out <- good_edges[!good_edges %in% plot_mx$edge]
    
    if(length(edges_out)==0){
      plot_mx_2 <- plot_mx %>%
        reshape2::dcast(formula=Source~Target,
                        fun.aggregate = sum,
                        value.var = "logpval") %>%
        column_to_rownames("Source") %>%
        as.matrix()
    }else{
      plot_mx_sup <- data.frame(edge=edges_out,logpval=0) %>%
        separate(edge,sep = "_",into = c("Source","Target"))
      
      plot_mx_2 <- plot_mx %>%
        dplyr::select(Source,Target,logpval) %>%
        rbind(plot_mx_sup) %>%
        reshape2::dcast(formula = Source~Target,
                        fun.aggregate = sum,
                        value.var = "logpval") %>%
        column_to_rownames("Source") %>%
        as.matrix()
    }
    
    rownames(plot_mx_2) <- str_to_title(str_to_lower(rownames(plot_mx_2)))
    colnames(plot_mx_2) <- str_to_title(str_to_lower(colnames(plot_mx_2)))
    
    plot_mx_2 <- plot_mx_2[row_order,col_order]
    
    fil <- paste0("figures/figure_02/panel_A/",folders[i],"/",sel,"_top_pairs_",yr,".png")
    png(filename = fil,res = 300,width = 750*3,height = 450*3)
    pheatmap::pheatmap(plot_mx_2,
                       cluster_rows = F,
                       cluster_cols = F,
                       color = myColor,
                       breaks = seq(0,max,max/1000),main = yr)
    dev.off()
  }
}

#make gifs from heatmaps from 1990 to 2018 (yearly)
for(i in 1:length(folders)){
  path <- paste0("figures/similarity_evolution/",folders[i],"/")
  imgs <- list.files(path,full.names = T)
  img_list <- lapply(imgs, image_read)
  img_joined <- image_join(img_list)
  img_animated <- image_animate(img_joined, delay = 50,loop = 1)
  path_out <- paste0(path,folders[i],"_evolution.gif")
  image_write(image = img_animated,
              path = path_out)
  print(i)
}

####figure 2B and 2C - similarity-to-paper ratio####
#Search pubmed for the number of papers with each disease pair from the
#figures above
#load MeSH terms for each disease (MeSH terms provided in data/ folder)
dis_MeSH_all <- readRDS("data/dis_MeSH_all.rds")

dis_MeSH_all <- dis_MeSH_all %>%
  dplyr::filter(!is.na(MESHTERM)) %>%
  dplyr::mutate(MESHTERM=paste0(MESHTERM,"[MeSH Terms]")) %>%
  dplyr::group_by(disease,class) %>%
  dplyr::summarise(MESHTERM=paste0(MESHTERM,collapse = " OR ")) %>%
  dplyr::mutate(MESHTERM=paste0("(",MESHTERM,")"))

top_pairs <- top_9_2018 %>%
  dplyr::filter(between %in% selected_between[c(1,2,4)]) %>%
  dplyr::pull(edge) %>%
  unique()

#create terms for disease pairs that will be searched in PubMed
search_grid <- data.frame(edge=top_pairs) %>%
  separate(col = "edge",sep = "_",into = c("a","b"),remove = F) %>%
  dplyr::left_join(dis_MeSH_all,by=c("a"="disease")) %>%
  dplyr::left_join(dis_MeSH_all,by=c("b"="disease")) %>%
  dplyr::mutate(search=paste0(MESHTERM.x," AND ",MESHTERM.y)) %>%
  dplyr::select(edge,search)

#search queries in PubMed using easyPubMed packege
#api_key <- "" INCLUDE YOUR NCBI API KEY HERE TO GO FASTER
for(i in 1:nrow(search_grid)){
  my_query <- search_grid$search[i]
  my_query <- get_pubmed_ids(as.character(my_query),
                             api_key = api_key)
  if(my_query$Count==0){
    next()
  }
  
  source("R/fetch_pubmed_data_mod.R")
  retmx <- my_query$Count
  my_abstracts_xml <- fetch_pubmed_data_2(my_query,retmax = retmx)
  all_xml <- articles_to_list(my_abstracts_xml)
  final_df <- do.call(rbind, lapply(all_xml, article_to_df,
                                    max_chars = -1, getAuthors = FALSE))
  papers_per_year <- final_df %>%
    dplyr::group_by(year) %>%
    dplyr::summarise(count=length(doi)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(edge=search_grid$edge[i]) %>%
    dplyr::select(edge,year,count)
  
  if(i==1){
    papers_per_year_all <- papers_per_year
  }else{
    papers_per_year_all <- rbind(papers_per_year_all,papers_per_year)
  }
  print(i)
}

papers_per_year_all$year <- as.numeric(papers_per_year_all$year)
saveRDS(papers_per_year_all,file = "intermediate/easyPubMed_top_dis_dis_clusterProfiler.RDS")
papers_per_year_all <- readRDS("intermediate/easyPubMed_top_dis_dis_clusterProfiler.RDS")

#plot number of papers each year between top pairs of diseases
pubmed_search <- reshape2::dcast(papers_per_year_all,
                                   formula = edge~year,
                                   value.var = "count",
                                   fun.aggregate = sum) %>%
  reshape2::melt() %>%
  dplyr::mutate(variable=as.numeric(as.character(variable))) %>%
  dplyr::rename(year=variable,count=value) %>%
  dplyr::group_by(edge) %>%
  dplyr::mutate(sum=sum(count)) %>%
  dplyr::filter(sum>0) %>%
  dplyr::ungroup() %>%
  separate(col = edge,sep = "_",into = c("Source","Target"),remove = F) %>%
  left_join(dis_class_df,by=c("Source"="disease")) %>%
  rename(class_source=class) %>%
  left_join(dis_class_df,by=c("Target"="disease")) %>%
  rename(class_target=class) %>%
  mutate(between=paste0(class_source,"_vs_",class_target))

#plot number of papers per year for top pairs in each comparison
for(i in 1:length(selected_between[c(1,2,4)])){
  sel <- selected_between[c(1,2,4)][i]
  p <- pubmed_search %>%
    filter(year<2019,
           between==sel) %>%
    ggplot(aes(x=year,y=count,color=edge))+
    geom_line()+
    geom_text_repel(data = pubmed_search %>%
                      separate(col = edge,sep = "_",into = c("Source","Target"),remove = F) %>%
                      filter(between==sel,
                             year==2018),
                    aes(label=Target,x=year,y=count,color=edge),
                    segment.size = .4,
                    segment.alpha = 0.5,
                    segment.colour = "black",
                    size=2,hjust=1,force = 1,
                    nudge_x = -10,nudge_y = 30)+
    scale_x_continuous(limits = c(1990,2018),breaks = c(seq(1990,2018,10),2018))+
    facet_wrap(facets = ~Source,scales = "fixed")+
    theme_minimal()+
    ylab("papers")+
    theme(legend.position = "none")
  fil <- paste0("figures/figure_S2/",sel,"_papers_per_year.pdf")
  pdf(file = fil,
      width = 9,height = 7)
  print(p)
  dev.off()
}

#calculate cumulative paper number per disease pair
all_years <- unique(pubmed_search$year)
for (i in 1:length(all_years)){
  yr <- all_years[i]
  pm <- pubmed_search %>%
    dplyr::filter(year <= yr) %>%
    dplyr::group_by(edge) %>%
    dplyr::summarise(count=sum(count)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(year=yr) %>%
    dplyr::select(edge,year,count)
  if(i==1){
    pubmed_search_cummulative <- pm
  }else{
    pubmed_search_cummulative <- rbind(pubmed_search_cummulative,pm)
  }
}

#calculate similarity score
pubmed_network <- pubmed_search_cummulative %>%
  separate(col = edge,sep = "_",into = c("Source","Target"),remove = F) %>%
  full_join(enr_all_out_2,by=c("Source","Target","edge","year")) %>%
  filter(edge %in% top_9_2018$edge)

pubmed_network$count[is.na(pubmed_network$count)] <- 0

pubmed_network$sim_paper <- pubmed_network$logpval/(pubmed_network$count+1)

pubmed_network_scored <- pubmed_network %>%
  filter(year==2018) %>%
  filter(!is.na(between))

write.csv(pubmed_network_scored,file="tables/table_S4.csv",row.names = F)

#Figure 2B - plot similarity-to-paper ratio histogram
p <- pubmed_network_scored %>%
  filter(between %in% selected_between[c(1,2,4)]) %>%
  mutate(between=factor(between,levels = selected_between[c(1,2,4)])) %>%
  ggplot(aes(sim_paper,fill=between))+
  geom_histogram(color="black")+
  scale_fill_manual(values = rev(c("#cb8adb","#72d172","#ef848b")))+
  xlab("similarity-to-paper ratio")+
  ylab("number of disease pairs")+
  facet_wrap(facets = ~between,ncol = 1)+
  theme_minimal()+
  theme(legend.position = "none")
pdf(file = "figures/figure_02/panel_B/sim_paper_score_histogram.pdf",
    width = 4,height = 10)
print(p)
dev.off()

#Figure 02C - papers x similarity plots of all top dis-dis pairs
top_9_2018_between_classes <- top_9_2018 %>%
  filter(between %in% selected_between[c(1,2,4)])

for(i in 1:length(top_9_2018_between_classes$edge)){
  ed <- top_9_2018_between_classes$edge[i]
  title <- unlist(str_split(ed,pattern = "_"))
  title <- paste0(str_to_title(tolower(title)),collapse = " vs ")
  p <- pubmed_network_2 %>%
    filter(edge==ed) %>%
    filter(year < 2019 & year >= 1990) %>%
    mutate(Source=str_to_title(tolower(Source)),
           Target=str_to_title(tolower(Target))) %>%
    ggplot(aes(x=year,y=count))+
    geom_line(color="red",size=1)+
    geom_line(data = pubmed_network_2 %>%
                filter(edge==ed) %>%
                filter(!is.na(logpval)) %>%
                mutate(Source=str_to_title(tolower(Source)),
                       Target=str_to_title(tolower(Target))),
              aes(x=year,y=logpval*10),color="blue",size=1)+
    scale_y_continuous(name = "papers",limits = c(0,1000),breaks = seq(0,1000,250),
                       sec.axis = sec_axis(~./10,name="genetic similarity",
                                           breaks = seq(0,90,20)))+
    scale_x_continuous(limits = c(1990,2018),
                       breaks = c(1990,2000,2010,2018))+
    theme_minimal()+
    theme(
      axis.title.y = element_text(color = "red", size=8),
      axis.title.y.right = element_text(color = "blue", size=8)
    )+
    ggtitle(label = title)+
    theme(plot.title = element_text(size = 8))
  
  title_fil <- gsub(title,pattern = " ",replacement = "_")
  fil <- paste0("figures/figure_02/panel_C/",title_fil,".pdf")
  pdf(fil,width = 2.33,height = 2.04)
  print(p)
  dev.off()
}
