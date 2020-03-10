#script to select common hubs between disease categories and plot line
#graphs of the number of diseases associated to each of them

#necessary packages
library(dplyr)
library(igraph)
library(tibble)
library(clusterProfiler)
library(ggplot2)
options(stringsAsFactors = F)

#set working directory
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#load necessary data
data("all_edges")
data("GMT")

# data("top_9_list")
# top_9_all <- unlist(top_9_list)

#importante objects
classes <- c("inflammatory","psychiatric","infectious")
color <- c("#ff743e","#00b7da","#68ad36")

#calculate degree of genes associated with diseaes in each cateogry in 2018
hubs <- list()
for (i in 1:3){
  all_edges_2018 <- edges_list[[i]] %>%
    dplyr::select(c(1,2,31)) %>%
    dplyr::filter(Doc.2018 > 0) %>%
    dplyr::select(-3)
  
  all_nodes_2018 <- data.frame(Id=unique(c(all_edges_2018$Source,all_edges_2018$Target)),
                               Class=c(rep("GENE",length(unique(all_edges_2018$Source))),
                                       rep("DISEASE",length(unique(all_edges_2018$Target))))
                               )
  
  graph <- graph_from_data_frame(d = all_edges_2018,
                                 directed = F,
                                 vertices = all_nodes_2018)
  
  degree <- as.data.frame(degree(graph)) %>%
    rownames_to_column("node")
  
  colnames(degree) <- c("Id","degree")
  degree <- merge(degree,all_nodes_2018,
                  by = "Id")
  genes <- degree %>%
    filter(Class=="GENE") %>%
    arrange(desc(degree)) %>%
    select(-3) %>%
    mutate(dis_class=classes[i])
  if (i==1){
    all_genes <- genes
  }else{
    all_genes <- rbind(all_genes,genes)
  }
  
  hub_genes <- genes$Id[1:50]
  hubs[[i]] <- hub_genes
  names(hubs)[i] <- classes[i]
  
  #perform enrichment of hub genes and plot results
  gsea_results <- enricher(gene = hub_genes,pvalueCutoff = 0.05,
                           pAdjustMethod = "none",
                           TERM2GENE = GMT)
  
  gsea_plot <- gsea_results@result %>%
    dplyr::select(ID,p.adjust) %>%
    arrange(desc(-log(p.adjust))) %>%
    mutate(logp = -log(p.adjust)) %>%
    top_n(wt = logp,n = 5)
  
  gsea_plot$ID <- factor(gsea_plot$ID,levels = rev(gsea_plot$ID))
  pallete <- colortools::sequential(color[i],plot = F)
  filename <- paste0("figures/",classes[i],"_hub_genes_enrichment.svg")
  svg(filename,width = 6.5,height = 1.7)
  p<-ggplot(gsea_plot,aes(x=ID,y=logp,fill=logp))+
    coord_flip()+
    geom_bar(stat = "identity")+
    scale_fill_gradient(low  = pallete[5],
                        high = pallete[21])+
    geom_hline(yintercept = -log(0.01),
               linetype="dashed",
               size=.3)+
    theme_minimal()
  print(p)
  dev.off()
}

#plot venn diagram of the hub genes shared between disease categories
common_hubs <- intersect(intersect(hubs[[1]],hubs[[2]]),hubs[[3]])
venn.diagram(x = hubs,
             category.names = classes,
             filename = "figures/common_hubs.svg",imagetype = "svg",
             output=T,resolution = 300,force.unique = T,height = 10,width = 10,
             
             # Circles
             lwd = 2,
             lty = 'blank',
             fill = color,
             
             # Numbers
             cex = 2,
             fontface = "bold",
             fontfamily = "sans",
             
             # Set names
             cat.cex = 2,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-27, 27, 135),
             cat.dist = c(0.055, 0.055, 0.085),
             cat.fontfamily = "sans",
             rotation = 1)

all_genes_edges <- all_genes %>%
  rename(Source=Id,Target=dis_class) %>%
  select(c(1,3,2))

write.csv(all_genes_edges,"data/all_genes_degrees.csv",row.names = F)
