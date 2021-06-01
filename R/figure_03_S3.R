#With this scritp, you can reproduce the analysis and plot the images that
#compose panel A-C in figure 03 of the paper and figure S3

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
library(limma)
#plotting
library(ggplot2)
library(ggrepel)
library(ggraph)

options(stringsAsFactors = F)

#load and create necessary objects
#pathways gmt
pathways_gmt <- fread("data/all_pathways_and_genes.csv")
reactome_gmt <- pathways_gmt %>%
  filter(db=="REACTOME")

#gene disease data
data("all_edges")
all_edges <- do.call(rbind,edges_list)
dis_dis_nodes <- as.data.frame(fread("data/dis_dis_nodes.csv"))
dis_class_df <- dis_dis_nodes %>%
  dplyr::select(1,3) %>%
  dplyr::rename(disease=1,class=2)
diseases <- unique(all_edges$Target)

#conversion table for genes in all_edges
load("data/gene_convertion_table.RData")
all_edges_converted <- all_edges %>%
  dplyr::left_join(gene_convertion_table,by=c("Source"="alias")) %>%
  dplyr::select(-1) %>%
  dplyr::select(converted,Target,everything())

#top 9 diseases
top_9_df <- readRDS("intermediate/top9_diseases_df.RDS")
top_9 <- top_9_df$disease

#run enrichment for genes related to top9 diseases in 2018
all_edges_top_9 <- all_edges_converted %>%
  filter(Doc.2018 > 0) %>%
  filter(Target %in% top_9)

for(i in 1:length(top_9)){
  genes <- all_edges_top_9 %>%
    filter(Target==top_9[i]) %>%
    pull(converted) %>% unique()
  enr <- enricher(genes,pAdjustMethod = "BH",TERM2GENE = reactome_gmt)
  enr_res <- enr@result
  enr_res$dis <- top_9[i]
  if(i==1){
    enr_res_all <- enr_res
  }else{
    enr_res_all <- rbind(enr_res_all,enr_res)
  }
}

#filter significant terms
enr_res_all_sig <- enr_res_all %>%
  filter(p.adjust < 0.01)

saveRDS(enr_res_all_sig,file = "intermediate/ORA_results_top9_diseases_2018.RDS")
enr_res_all_sig <- readRDS(file = "intermediate/ORA_results_top9_diseases_2018.RDS")

#run term-term network analysis to find term domains
enrichment_term_gene <- enr_res_all_sig %>%
  dplyr::select(ID,geneID) %>%
  separate_rows(geneID,sep = "/") %>%
  dplyr::rename(term=1,gene=2)

#create a list containing all term genes (also in edges_all)
term_gene_list <- split(enrichment_term_gene$gene,
                        enrichment_term_gene$term)

#run enrichment using lapply and create enrichment dataframe for term-term
term_term_enrichment <- lapply(term_gene_list,enricher,
                               pAdjustMethod = "BH",TERM2GENE=enrichment_term_gene)

for(i in 1:length(term_term_enrichment)){
  obj <- term_term_enrichment[[i]]
  if(class(obj)=="NULL"){
    next()
  }
  res <- obj@result
  res$term <- names(term_term_enrichment)[i]
  if(i==1){
    term_term_all <- res
  }else{
    term_term_all <- rbind(term_term_all,res)
  }
}

term_term_network <- term_term_all %>%
  dplyr::select(ID,term,p.adjust) %>%
  dplyr::mutate(weight=-log(p.adjust)) %>%
  dplyr::select(-3) %>%
  dplyr::rename(from=1,to=2) %>%
  filter(from != to,
         weight > -log(0.01))
rownames(term_term_network) <- NULL

#create igraph object and run louvain algorithm to detect clusters
g <- graph_from_data_frame(term_term_network,directed = F,
                           vertices = unique(c(term_term_network$from,
                                               term_term_network$to)))
louv <- cluster_louvain(graph = g,weights = term_term_network$weight)
mods <- data.frame(node=louv$names,mod=louv$membership)
degree <- strength(graph = g,vids = V(g),mode = "all",weights = E(g)$weight)
nodes <- mods %>%
  mutate(degree = degree) %>%
  mutate(node=ifelse(str_sub(node,
                             start = nchar(node),
                             end = nchar(node))==" ",
                     no=node,
                     yes=str_sub(node,
                                 start = 1,
                                 end = nchar(node)-1)
  ))

#load manual annotation of clusters (made outside R by authors. In data/)
nodes_annot <- fread("data/term_term_nodes_annot.csv")
nodes_2 <- nodes %>%
  dplyr::select(node) %>%
  left_join(nodes_annot,by="node")

label_terms <- nodes_2 %>%
  group_by(mod) %>%
  top_n(degree,n = 2) %>%
  ungroup()

nodes_2 <- nodes_2 %>%
  mutate(label=ifelse(node %in% label_terms$node,
                      node,""))

table_S5 <- enr_res_all_sig %>%
  dplyr::left_join(nodes_annot,by = c("ID"="node")) %>%
  dplyr::filter(!is.na(annot)) %>%
  dplyr::select(-degree)

#Save table S5 in .csv format (enriched terms for top 9 diseases in 2018)
write.csv(table_S5,file = "tables/table_S5.csv",row.names = F)

#update igraph object with degree, annot and labels
V(g)$degree <- nodes_2$degree
V(g)$mod <- nodes_2$annot
V(g)$label <- nodes_2$label

#remove smaller component (contains only 3 terms)
components <- igraph::clusters(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(g)[components$membership == biggest_cluster_id]
g <- igraph::induced_subgraph(g, vert_ids)

####Figure 3A network####
#plot term-term network with mod annotations
load("data/term_term_network_palette.RData")
p <- ggraph(graph = g,layout = "auto")+
  geom_edge_link(aes(color=log(weight)))+
  geom_node_point(aes(size=degree,fill=mod),color="black",pch=21)+
  #geom_node_text(aes(label = label),repel=TRUE)+
  scale_fill_manual(values = pal)+
  scale_edge_color_continuous(low = "white",high = "grey40")+
  theme_void()
pdf(file = "figures/figure_03/panel_A/term_term_network_mods.pdf",
    width = 10.3,height = 6.6)
print(p)
dev.off()

pdf(file = "figures/figure_03/panel_A/term_term_network_mod_facets.pdf",
    width = 12.38,height = 8.06)
p <- ggraph(graph = g,layout = "auto")+
  geom_edge_link(aes(color=log(weight)))+
  geom_node_point(aes(size=degree,fill=mod),color="black",pch=21)+
  scale_fill_manual(values = pal)+
  scale_edge_color_continuous(low = "white",high = "grey40")+
  theme_void()+
  theme(legend.position = "none")+
  facet_nodes(facets = ~mod)+
  theme(strip.text.x = element_text(size = 7))
print(p)
dev.off()

#run enrichment for genes in each top 9 disease per year and plot results in network
years <- 1990:2018
for(i in 1:length(years)){
  yr <- years[i]
  col <- which(grepl(yr,colnames(all_edges_converted)))
  df <- all_edges_converted %>%
    dplyr::filter(Target %in% top_9) %>%
    dplyr::select(1,2,col) %>%
    dplyr::rename(gene=1,dis=2,docs=3) %>%
    dplyr::filter(docs > 0)
  gene_list <- split(df$gene,df$dis)
  enr <- lapply(gene_list,enricher,pAdjustMethod = "BH",TERM2GENE = reactome_gmt)
  for(j in 1:length(enr)){
    enr_res <- enr[[j]]
    if(class(enr_res)=="NULL"){next}
    enr_res <- enr_res@result
    enr_res$year <- yr
    enr_res$dis <- names(enr)[j]
    if(i==1 & j==1){
      enr_res_all_dis <- enr_res
    }else{
      enr_res_all_dis <- rbind(enr_res_all_dis,enr_res)
    }
  }
}  

saveRDS(enr_res_all_dis,file = "intermediate/ORA_results_top9_diseases_all_years.RDS")

enr_res_all_dis_selected <- enr_res_all_dis %>%
  filter(p.adjust < 0.01) %>%
  mutate(ID=ifelse(str_sub(ID,
                           start = nchar(ID),
                           end = nchar(ID))==" ",
                   no=ID,
                   yes=str_sub(ID,
                               start = 1,
                               end = nchar(ID)-1)
  )) %>%
  left_join(dis_class_df,by=c("dis"="disease")) %>%
  filter(ID %in% terms_in_network) %>%
  left_join(nodes_2 %>% dplyr::select(node,annot),by=c("ID"="node"))

enrichment_df_plot <- enr_res_all_dis_selected %>%
  dplyr::select(ID,p.adjust,year,dis,class,annot) %>%
  dplyr::mutate(score=-log(p.adjust))

#plot enrichment results for each module in boxplots
enrichment_boxplot_plot <- enrichment_df_plot %>%
  filter(year==2018)

cols <- c("#68ad36","#ff743e", "#00b7da")
classes <- unique(enrichment_df_plot$class)

####Figure 3A boxplots####
p <- enrichment_boxplot_plot %>%
  #filter(dis %in% top_9[c(1,2,6,12,11,10,25,26,20,19)]) %>%
  # group_by(dis,class,annot) %>%
  # summarise(score=mean(score)) %>%
  # arrange(class,dis,desc(score)) %>%
  # mutate(annot=factor(annot,levels = unique(annot))) %>%
  mutate(class=factor(class,levels = c("INFECTIOUS","INFLAMMATORY","PSYCHIATRIC"))) %>%
  ggplot(aes(x=class,y=log(score),fill=class))+
  geom_jitter(aes(color=class),shape=1)+
  geom_boxplot(outlier.color = NA,alpha=0.3)+
  #scale_fill_gradient(low = "white",high = "red3")+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  facet_wrap(facets = ~annot,nrow = 1)+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        strip.text.x = element_text(size = 4))
pdf(file = "figures/figure_03/panel_A/annot_enrichment_boxplots_2018.pdf",
    width = 13.92,height = 3.22)
print(p)
dev.off()

#####Figure 3B - networks####
#plot term-term network with nodes colored according to enrichment values
#in each disease in 2018
terms_in_network <- V(g)$name
g2 <- g
dis_enrichment <- enr_res_all_sig %>%
  mutate(score=-log(p.adjust)) %>%
  dcast(formula = ID~dis,value.var = "score") %>%
  filter(ID %in% terms_in_network)

dis_enrichment[is.na(dis_enrichment)] <- 0
dis_enrichment <- dis_enrichment[match(V(g2)$name,dis_enrichment$ID),]
dis_enrichment <- dis_enrichment %>%
  mutate(ID=ifelse(str_sub(ID,
                           start = nchar(ID),
                           end = nchar(ID))==" ",
                   no=ID,
                   yes=str_sub(ID,
                               start = 1,
                               end = nchar(ID)-1)
  ))

nodes_3 <- nodes_2 %>%
  dplyr::select(-label) %>%
  left_join(dis_enrichment,by=c("node"="ID"))

for(i in 1:length(top_9)){
  dis <- top_9[i]
  dis_fil <- gsub(dis,pattern = " ",replacement = "_")
  col <- which(colnames(dis_enrichment)==dis)
  V(g2)$enrichment <- log(dis_enrichment[,col]+1)
  fil1 <- paste0("figures/figure_03/panel_B/term_term_network_full_enrichment_",dis_fil,".pdf")
  pdf(file = fil1,
      width = 7.3,height = 6.6)
  p <- ggraph(graph = g2,layout = "auto")+
    geom_edge_link(aes(color=log(weight)))+
    geom_node_point(aes(size=degree,fill=enrichment),color="black",pch=21)+
    scale_fill_gradient(low = "white",high = "red3")+
    scale_edge_color_continuous(low = "white",high = "grey40")+
    theme_void()+
    ggtitle(dis)
  print(p)
  dev.off()
  
  fil2 <- paste0("figures/figure_03/panel_B/term_term_network_mods_enrichment_",dis_fil,".pdf")
  pdf(file = fil2,
      width = 12.38,height = 8.06)
  p <- ggraph(graph = g2,layout = "auto")+
    geom_edge_link(aes(color=log(weight)))+
    geom_node_point(aes(size=degree,fill=enrichment),color="black",pch=21)+
    scale_fill_gradient(low = "white",high = "red3")+
    scale_edge_color_continuous(low = "white",high = "grey40")+
    theme_void()+
    theme(legend.position = "none")+
    facet_nodes(facets = ~mod)+
    theme(strip.text.x = element_text(size = 7))+
    ggtitle(dis)
  print(p)
  dev.off()
}

####Figure 3B boxplots####
p <- enrichment_heatmap_plot %>%
  #filter(dis %in% top_9[c(25,26,20,19,1,2,6,12,11,10)]) %>%
  filter(dis %in% top_9[c(18,26,4,17,10,2)]) %>%
  # group_by(dis,class,annot) %>%
  # summarise(score=mean(score)) %>%
  arrange(class,dis,desc(score)) %>%
  #mutate(dis=factor(dis,levels = top_9[c(25,26,20,19,1,2,6,12,11,10)])) %>%
  mutate(dis=factor(dis,levels = top_9[c(18,26,4,17,10,2)])) %>%
  mutate(annot=factor(annot,levels = unique(annot))) %>%
  ggplot(aes(x=annot,y=log(score),fill=class))+
  geom_jitter(aes(color=class),shape=1)+
  geom_boxplot(outlier.color = NA,alpha=0.3)+
  #scale_fill_gradient(low = "white",high = "red3")+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  coord_flip()+
  facet_wrap(facets = ~dis,nrow = 1)+
  theme_minimal()+
  theme(axis.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 5),
        legend.position = "none")
pdf(file = "figures/figure_03/panel_B/annot_enrichment_boxplots_per_dis.pdf",
    width = 9,height = 3)
print(p)
dev.off()

####Figure S3 - boxplots####
p <- enrichment_boxplot_plot %>%
  #filter(dis %in% top_9[c(25,26,20,19,1,2,6,12,11,10)]) %>%
  filter(!dis %in% top_9[c(18,26,4,17,10,2)]) %>%
  # group_by(dis,class,annot) %>%
  # summarise(score=mean(score)) %>%
  arrange(class,dis,desc(score)) %>%
  #mutate(dis=factor(dis,levels = top_9[c(25,26,20,19,1,2,6,12,11,10)])) %>%
  mutate(dis=factor(dis,levels = unique(dis))) %>%
  mutate(annot=factor(annot,levels = unique(annot))) %>%
  ggplot(aes(x=annot,y=log(score),fill=class))+
  geom_jitter(aes(color=class),shape=1)+
  geom_boxplot(outlier.color = NA,alpha=0.3)+
  #scale_fill_gradient(low = "white",high = "red3")+
  scale_fill_manual(values = cols)+
  scale_color_manual(values = cols)+
  coord_flip()+
  facet_wrap(facets = ~dis,nrow = 3)+
  theme_minimal()+
  scale_x_discrete(name="cluster")+
  scale_y_continuous(name="log(enrichment score)")+
  theme(axis.text.y = element_text(size = 7),
        strip.text.x = element_text(size = 5),
        legend.position = "none")
pdf(file = "figures/figure_S3/annot_enrichment_boxplots_all_dis.pdf",
    width = 10,height = 6.6)
print(p)
dev.off()

####Figure 3C####
#plot yearly enrichment by mod and dis category
annots <- unique(enrichment_df_plot$annot)
p <- enrichment_df_plot %>%
  group_by(dis,class,year,annot) %>%
  summarise(score=mean(score)) %>%
  ungroup() %>%
  mutate(class=factor(class,levels = c("INFECTIOUS","INFLAMMATORY","PSYCHIATRIC"))) %>%
  arrange(class) %>%
  mutate(dis=factor(dis,unique(dis))) %>%
  mutate(annot=factor(annot,unique(annot))) %>%
  ggplot(aes(x=year,y=dis,fill=class,color=class))+
  geom_ridgeline(aes(height=score),scale=0.05,alpha=0.274,size=.2)+
  scale_fill_manual(values=c("#68ad36", "#ff743e", "#00b7da"))+
  scale_color_manual(values=c("#68ad36", "#ff743e", "#00b7da"))+
  scale_x_continuous(breaks=c(1990,2005,2018))+
  theme_ridges(font_size = 5)+
  theme_minimal()+
  theme(legend.position = 'none',
        strip.text.x = element_text(size = 6,angle = 30),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6))+
  facet_wrap(facets = ~annot,nrow=1)
pdf(file = "figures/figure_03/panel_C/annot_enrichment_evolution_top9.pdf",
    width = 9.4,height = 5.4)
print(p)
dev.off()
