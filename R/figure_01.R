#With this scritp, you can reproduce the analysis and plot the images that
#compose panel A-E in figure 01 of the paper

#load necessary packages
#file management
library(data.table)
#data manipulation
library(dplyr)
library(tidyverse)
library(reshape2)
#network/gene analysis
library(igraph)
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

#write Table S1 in .csv format (studied diseases in each class)
write.csv(dis_class_df,file = "tables/table_S1.csv",row.names = F)

####Figure 01A - disease-disease networks from 1990 to 2018####
#calculate disease-disease similarity in each year according to gene
#sharing between diseases
#calculate number of genes per disease per year
dis_dis_edges <- as.data.frame(t(combn(diseases,m = 2))) %>%
  dplyr::rename(dis_1=1,dis_2=2)
years_df <- as.data.frame(matrix(nrow = nrow(dis_dis_edges),ncol = 29))
colnames(years_df) <- seq(1990,2018)
dis_dis_edges <- cbind(dis_dis_edges,years_df)

for(i in 1:length(years)){
  yr <- years[i]
  col <- which(grepl(yr,colnames(all_edges)))
  df <- all_edges %>%
    dplyr::select(c(1,2,col)) %>%
    dplyr::rename(docs=3) %>%
    dplyr::filter(docs > 1)
  genes_per_dis_yr <- as.data.frame(table(df$Target)) %>%
    rename(disease=1)
  colnames(genes_per_dis_yr)[2] <- yr
  if(i==1){
    genes_per_year <- genes_per_dis_yr
  }else{
    genes_per_year <- genes_per_year %>%
      full_join(genes_per_dis_yr,by="disease")
  }
  
  for(j in 1:nrow(dis_dis_edges)){
    dis_1 <- dis_dis_edges[j,1]
    dis_2 <- dis_dis_edges[j,2]
    if((!dis_1 %in% df$Target) | (!dis_2 %in% df$Target)){next()}
    
    dis_1_genes <- df %>%
      dplyr::filter(Target==dis_1) %>%
      dplyr::pull(Source) %>% unique()
    dis_2_genes <- df %>%
      dplyr::filter(Target==dis_2) %>%
      dplyr::pull(Source) %>% unique()
    comm <- dplyr::intersect(dis_1_genes,dis_2_genes)
    pval <- round(phyper(q = length(comm)-1,
                         m = length(dis_1_genes),
                         n = length(unique(df$Source))-length(dis_1_genes),
                         k = length(dis_2_genes),
                         lower.tail = F),
                  digits = 10000)
    log_pval <- -log(pval)
    dis_dis_edges[j,col] <- log_pval
  }
}

genes_per_year[is.na(genes_per_year)] <- 0
saveRDS(object = genes_per_year,file = "intermediate/genes_per_year.RDS")
genes_per_year <- readRDS("intermediate/genes_per_year.RDS")

dis_dis_edges[is.na(dis_dis_edges)] <- 0
saveRDS(object = dis_dis_edges,file = "intermediate/dis_dis_edges.RDS")
dis_dis_edges <- readRDS("intermediate/dis_dis_edges.RDS")

#write Table S2 in .csv format (cummulative genes per disease per year)
write.csv(genes_per_year,file = "tables/table_S2.csv",row.names = F)

#get top9 diseases connected to more genes in each disease class in 2018
top9_diseases_df <- genes_per_year %>%
  dplyr::select(1,30) %>%
  dplyr::rename(genes=2) %>%
  dplyr::left_join(dis_class_df,by="disease") %>%
  group_by(class) %>%
  top_n(n=9,wt = genes) %>%
  ungroup() %>%
  filter(disease!="MEASLES")

saveRDS(top9_diseases_df,file = "intermediate/top9_diseases_df.RDS")

#Plot dis-dis network in each year
#transform dis_dis_edges into weights matrix and remove not significant edges
dis_dis_weights <- dis_dis_edges %>%
  dplyr::mutate(edge=paste0(dis_1,"_",dis_2)) %>%
  column_to_rownames("edge") %>%
  dplyr::select(-c(1,2)) %>%
  as.matrix()

dis_dis_weights[dis_dis_weights < -log(0.05)] <- NA

#determine years to plot/analyze
years_plot <- c(1990,2000,2010,2018)

#determine category colors and build pallete for weights
colors <- c("#ff743e","#68ad36","#00b7da")
names(colors) <- c("INFLAMMATORY","INFECTIOUS","PSYCHIATRIC")

saveRDS(colors,file = "intermediate/colors.RDS")

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

weights_log_vec <- c()
for (i in 1:length(years_plot)){
  yr <- as.character(years_plot[i])
  ed <- dis_dis_edges %>%
    dplyr::select(1,2,yr) %>%
    dplyr::rename(from=1,to=2,weight=3) %>%
    dplyr::filter(weight >= -log(0.05))
  ed_out <- ed
  colnames(ed_out) <- c("Source","Target",yr)
  ed_out$edge <- paste0(ed_out$Source,"_",ed_out$Target)
  ed_out <- ed_out %>%
    dplyr::select(edge,3)
  weights_log_vec <- c(weights_log_vec,log(ed_out[,2]+1))
}  

#create palette for edge weights to keep constant from year to year
#(the edge color will be proportional to the max overall value)
mypal <- colorRampPalette(c("grey88", "black"))(10000)
weights_log_vec_sorted <- sort(unique(weights_log_vec))
pal <- map2color(weights_log_vec_sorted,mypal)
pal_final <- weights_log_vec_sorted
names(pal_final) <- pal

#create a network for each year with dis-dis relationships
#edge weights = -log(pval) of the dis-dis gene sharing
#node size = number of genes associated to disease in each year
p_list <- list()
top_9 <- top9_diseases_df$disease
for (i in 1:length(years_plot)){
  yr <- as.character(years_plot[i])
  #subset dis_dis_edges for year
  ed <- dis_dis_edges %>%
    dplyr::select(1,2,yr) %>%
    dplyr::rename(from=1,to=2,weight=3) %>%
    dplyr::filter(weight >= -log(0.05))
  
  #get nodes of the year
  md <- dis_dis_nodes %>%
    dplyr::filter(Id %in% unique(c(ed$from,ed$to)))
  
  #build igraph object
  g <- graph_from_data_frame(ed,directed = F,vertices = md)
  
  #set node size by number of genes each year
  gpy_cols <- which(colnames(genes_per_year)==yr)
  order <- quos(V(g)$name)
  
  gpy <- genes_per_year %>%
    dplyr::select(1,gpy_cols) %>%
    filter(disease %in% V(g)$name)
  gpy <- gpy[match(V(g)$name,gpy$dis),]
  
  #set nodes class, lables (top 9 only), color (class) and size (genes per year)
  V(g)$class <- dis_dis_nodes$Class[dis_dis_nodes$Id %in% V(g)$name]
  V(g)$label <- ifelse(V(g)$name %in% top_9,V(g)$name,"")
  V(g)$color <- colors[match(V(g)$class,names(colors))]
  V(g)$size <- gpy[,2]
  
  #set edge weight (-log(pval) dis-dis gene sharing)
  E(g)$weight <- log(ed$weight+1)
  E(g)$width <- E(g)$weight
  E(g)$color <- "grey80"
  
  pal_in <- names(pal_final[pal_final %in% E(g)$weight])
  
  #detect number of components in year network
  comps <- components(g,mode = "strong")
  #create one network per component (only 1 if ncomps = 1)
  #add plot to list acording to year and component number (e.g.:1990_comp1)
  if(comps$no > 1 & yr != "1990"){
    for(k in 1:comps$no){
      sub_v <- names(comps$membership)[comps$membership==k]
      g2 <- subgraph(g,v=sub_v)
      plot_name <- paste0(yr,"_comp_",k)
      lay <- create_layout(g2, layout = 'igraph', algorithm = 'nicely')
      p <- ggraph(lay)+
        geom_edge_link(aes(color=weight,alpha=weight),width=.5)+
        geom_node_point(aes(size=size),color="black")+  
        geom_node_point(aes(color=class,size=size))+
        scale_edge_color_gradient(low = "grey85",high = "grey70")+
        scale_edge_color_gradient(low = pal_in[1],high = pal_in[length(pal_in)])+
        scale_size_continuous(limits = c(0,900))+
        scale_edge_alpha_continuous(range = c(0,1))+
        scale_color_manual(values = colors)+
        theme_graph()
      #print(p)
      p_list[[plot_name]] <- p
    }
  }else{
    lay <- create_layout(g, layout = 'igraph', algorithm = 'nicely')
    plot_name <- paste0(yr,"_comp_1")
    p <- ggraph(lay)+
      geom_edge_link(aes(color=weight,alpha=weight),width=.5)+
      geom_node_point(aes(size=size),color="black")+
      geom_node_point(aes(color=class,size=size))+
      scale_edge_color_gradient(low = pal_in[1],high = pal_in[length(pal_in)])+
      scale_size_continuous(limits = c(0,900),range = c(0.01,9))+
      scale_edge_alpha_continuous(range = c(0,1))+
      scale_color_manual(values = colors)+
      theme_graph()
    #print(p)
    p_list[[plot_name]] <- p
  }
}

#plot each network in svg (with legend) and png (without legend)
for(i in 1:length(p_list)){
  p <- p_list[[i]]
  plot_name <- names(p_list)[i]
  
  name_svg <- paste0("figures/figure_01/panel_A/",plot_name,"_ggraph.svg")
  svg(name_svg,width = 11.21,height = 9.53)
  print(p)
  dev.off()
  
  name_png <- paste0("figures/figure_01/panel_A/",plot_name,"_ggraph.png")
  png(name_png,width = 3000,height = 3000,units = "px",res = 500)
  print(p+theme(legend.position = "none"))
  dev.off()
}

####Figure 01B - number of genes in each category per year####
genes_per_year_melt <- melt(genes_per_year)
genes_per_year_melt$variable <- as.numeric(as.character(genes_per_year_melt$variable))

genes_per_year_all <- genes_per_year_melt %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(disease="all",
            value=sum(value)) %>%
  dplyr::ungroup() %>%
  dplyr::rename(class=disease) %>%
  dplyr::select(class,variable,value)

colors_plot <- c("black",colors)
names(colors_plot)[1] <- "all"
p <- genes_per_year_melt %>%
  dplyr::left_join(dis_class_df,by="disease") %>%
  dplyr::group_by(class,variable) %>%
  dplyr::summarise(value=sum(value)) %>%
  dplyr::ungroup() %>%
  dplyr::select(class,variable,value) %>%
  rbind(genes_per_year_all) %>%
  ggplot(aes(x=variable,y=value,color=class))+
  geom_line()+
  scale_color_manual(values = colors_plot)+
  scale_x_continuous(name = "year",breaks = years_plot)+
  scale_y_continuous(name = "genes")+
  theme_classic()
pdf(file = "figures/figure_01/panel_B/genes_per_year.pdf",
    width = 5,height = 3.6)
print(p)
dev.off()

####Figure 01C - distribution of genes in 2018 per disease category####
genes_2018 <- genes_per_year_melt %>%
  filter(variable==2018)

#violin plot of genes per disease in each disease class
p <- genes_2018 %>%
  left_join(dis_class_df,by="disease") %>%
  ggplot(aes(x=class,y=value))+
  geom_violin(aes(fill=class),show.legend = F)+
  geom_jitter(aes(fill=class),width = .3,shape=21,colour="black",show.legend = F)+
  scale_fill_manual(values=colors)+
  scale_y_continuous(name = "genes")+
  scale_y_continuous(name = "category")+
  theme_minimal()
pdf("figures/figure_01/panel_C/genes_in_2018_violin_plot.pdf",
    width = 3.63,height = 3.57)
print(p)
dev.off()

####Figure 01D - new genes per year in each disease class####
#calculate the yearly variation of genes
gene_variation <- as.data.frame(t(apply(genes_per_year[,-1], 1, diff))) %>%
  mutate(disease=genes_per_year$disease) %>%
  select(disease,everything())

saveRDS(gene_variation,file = "intermediate/gene_variation_per_year.RDS")

#write Table S3 in .csv format (new genes per disease per year)
write.csv(gene_variation,file = "tables/table_S3.csv",row.names = F)

gene_variation_melt <- melt(gene_variation)
gene_variation_melt$variable <- as.numeric(as.character(gene_variation_melt$variable))

#plot the yearly variation of genes in each class
#line plot
p <- gene_variation_melt %>%
  left_join(dis_class_df,by="disease") %>%
  dplyr::group_by(class,variable) %>%
  dplyr::summarise(value=sum(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(variable != 2018) %>%
  ggplot(aes(x=variable,y=value,color=class))+
  geom_line()+
  scale_color_manual(values = colors)+
  facet_wrap(facets = ~class)+
  theme_minimal()
pdf("figures/figure_01/panel_D/new_genes_per_year_class_line.pdf",
    width = 7.5,height = 3.6)
print(p)
dev.off()

#ridge plot
p <- gene_variation_melt %>%
  left_join(dis_class_df,by="disease") %>%
  dplyr::group_by(class,variable) %>%
  dplyr::summarise(value=sum(value)) %>%
  dplyr::ungroup() %>%
  dplyr::filter(variable != 2018) %>%
  dplyr::mutate(class=factor(class,
                             levels = rev(c("INFECTIOUS","INFLAMMATORY","PSYCHIATRIC")))) %>%
  ggplot(aes(x=variable,y=class,fill=class,color=class))+
  geom_ridgeline(aes(height=value),scale=.005,alpha=0.274,size=.2)+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  scale_x_discrete(limits=c(1991,2000,2010,2017))+
  theme_ridges(font_size = 10)+
  theme(legend.position = 'none')
pdf(file = "figures/figure_01/panel_D/new_genes_per_year_class_ridge.pdf",
    width = 3.43,height = 3.57)
print(p)
dev.off()

####Figure 01E - new genes per year in each disease####
#plot the yearly variation of genes in all diseases
p <- gene_variation_melt %>%
  dplyr::left_join(dis_class_df,by="disease") %>%
  ggplot(aes(x=variable,y=disease,fill=class,color=class))+
  geom_ridgeline(aes(height=value,color=colors_plot_2),
                 scale=0.05,alpha=0.274,color="transparent")+
  scale_fill_manual(values=colors)+
  scale_y_discrete(limits=diseases)+
  scale_x_discrete(limits=c(1991,2000,2010,2017))+
  theme_ridges(font_size = 7)+
  theme(legend.position = 'none')
pdf(file = "figures/figure_01/panel_E/new_genes_per_year_per_disease_all.pdf",
    width = 4.35,height = 9.5)
print(p)
dev.off()
