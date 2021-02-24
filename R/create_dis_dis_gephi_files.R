#packages
library(tidyverse)
library(data.table)
library(dplyr)
library(readr)
library(igraph)
library(ggraph)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = F)

#set working directory
dir <- "~/Área de Trabalho/Papers_Helder/Evolution_all/"
setwd(dir)

#load dis-dis edges
dis_dis_edges <- fread("dis_dis_edges.csv",header = T)
dis_dis_nodes <- fread("dis_dis_nodes.csv")

dis_dis_weights <- dis_dis_edges %>%
  mutate(edge=paste0(Source,"_",Target)) %>%
  column_to_rownames("edge") %>%
  dplyr::select(-c(1,2)) %>%
  as.matrix()

dis_dis_weights[dis_dis_weights < -log(0.05)] <- NA

dis_dis_weights_scaled <- scale(dis_dis_weights,center = F)
years <- seq(1990,2018,4)
years <- c(1990,2000,2010,2018)

wdegs_df <- as.data.frame(matrix(nrow=99,ncol=1))
colnames(wdegs_df) <- "dis"
wdegs_df$dis <- unique(c(dis_dis_edges$Source,dis_dis_edges$Target))
wdegs_log_df <- wdegs_df
glist <- list()
edlist <- list()
mdlist <- list()
colors <- c("#ff743e","#68ad36","#00b7da")
names(colors) <- c("INFLAMMATORY","INFECTIOUS","PSYCHIATRIC")

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

weights_log_vec <- c()
for (i in 1:length(years)){
  yr <- as.character(years[i])
  ed <- dis_dis_edges %>%
    select(1,2,yr) %>%
    rename(from=1,to=2,weight=3) %>%
    filter(weight >= -log(0.05))
  ed_out <- ed
  colnames(ed_out) <- c("Source","Target",yr)
  ed_out$edge <- paste0(ed_out$Source,"_",ed_out$Target)
  ed_out <- ed_out %>%
    dplyr::select(edge,3)
  weights_log_vec <- c(weights_log_vec,log(ed_out[,2]+1))
}  

mypal <- colorRampPalette(c("grey88", "darkgreen"))(10000)
x <- sort(unique(weights_log_vec))
pal <- map2color(x,mypal)
pal_final <- x
names(pal_final) <- pal

load("~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/data/genes_per_year.RData")

p_list <- list()
#weights_log_vec <- c()
for (i in 1:length(years)){
  yr <- as.character(years[i])
  ed <- dis_dis_edges %>%
    select(1,2,yr) %>%
    rename(from=1,to=2,weight=3) %>%
    filter(weight >= -log(0.05))
  
  md <- dis_dis_nodes %>%
    filter(Id %in% unique(c(ed$from,ed$to)))
  
  edlist[[i]] <- ed
  
  ed_out <- ed
  colnames(ed_out) <- c("Source","Target",yr)
  ed_out$edge <- paste0(ed_out$Source,"_",ed_out$Target)
  ed_out <- ed_out %>%
    dplyr::select(edge,3)
  
  #weights_log_vec <- c(weights_log_vec,log(ed_out[,2]+1))
  
  if(i==1){
    ed_out_all <- ed_out
  }else{
    ed_out_all <- full_join(ed_out_all,ed_out,by="edge")
  }
  
  mdlist[[i]] <- md
  
  g <- graph_from_data_frame(ed,directed = F,vertices = md)
  glist[[i]] <- g
  
  # wdegs <- as.data.frame(strength(graph = g,vids = V(g),
  #                                 weights = E(g)$weight,mode = "total")) %>%
  #   rownames_to_column("dis")
  # 
  # colnames(wdegs) <- c("dis",yr)
  # wdegs_log <- wdegs
  # wdegs_log[,2] <- log(wdegs_log[,2]+1)
  # wdegs_log_df <- full_join(wdegs_log_df,wdegs_log,by="dis")
  # wdegs_df <- full_join(wdegs_df,wdegs,by="dis")
  
  gpy_cols <- which(colnames(genes_per_year_all)==yr)
  order <- quos(V(g)$name)
  gpy <- genes_per_year_all %>%
    dplyr::select(1,all_of(gpy_cols)) %>%
    filter(dis %in% V(g)$name)
  gpy <- gpy[match(V(g)$name,gpy$dis),]
  
  V(g)$class <- dis_dis_nodes$Class[dis_dis_nodes$Id %in% V(g)$name]
  V(g)$color <- colors[match(V(g)$class,names(colors))]
  V(g)$size <- gpy[,2]
  #V(g)$size <- wdegs_log[,2]
  E(g)$weight <- log(ed_out[,2]+1)
  E(g)$width <- E(g)$weight
  E(g)$color <- "grey80"
  name <- paste0("Figures/Figure_01/dis_dis_networks_v2/",yr,"_ggraph.pdf")
  pal_in <- names(pal_final[pal_final %in% E(g)$weight])
  lay <- create_layout(g, layout = "fr")
  p <- ggraph(lay)+
      geom_edge_link(aes(color=weight,alpha=weight))+
      geom_node_point(aes(color=class,size=size))+
      #scale_edge_color_gradient(low = "grey80",high = "darkgreen")+
      scale_edge_color_gradient(low = pal_in[1],high = pal_in[length(pal_in)])+
      scale_size_continuous(limits = c(0,900))+
      scale_color_manual(values = colors)+
      theme_graph()
  #print(p)
  #pdf(name,width = 10,height = 10)
  p_list[[i]] <- p
  #dev.off()
}

for(i in 1:4){
  p <- p_list[[i]]
  name_svg <- paste0("Figures/Figure_01/dis_dis_networks_v2/",years[i],"_ggraph.svg")
  svg(name_svg,width = 6.12,height = 5.95)
  print(p)
  dev.off()
  
  name_png <- paste0("Figures/Figure_01/dis_dis_networks_v2/",years[i],"_ggraph.png")
  png(name_png,width = 3000,height = 3000,units = "px",res = 500)
  print(p+theme(legend.position = "none"))
  dev.off()
}


wdegs_df[,-1][is.na(wdegs_df[,-1])] <- 0
wdegs_log_df[,-1][is.na(wdegs_log_df[,-1])] <- 0

dis_dis_nodes_2 <- dis_dis_nodes %>%
  left_join(wdegs_log_df,by=c("Id"="dis"))

dis_dis_edges_2 <- ed_out_all %>%
  separate(edge,into = c("Source","Target"),sep = "_")

dis_dis_edges_2[,-c(1,2)][is.na(dis_dis_edges_2[,-c(1,2)])] <- 0

write.csv(dis_dis_nodes_2,file = "dis_dis_nodes_2.csv",row.names = F)
write.csv(dis_dis_edges_2,file = "dis_dis_edges_2.csv",row.names = F)

wdegs_mx <- as.matrix(wdegs_df[,-1])
min_wd <- min(wdegs_mx[wdegs_mx>0])
max_wd <- max(wdegs_mx)

weights_mx <- as.matrix(dis_dis_edges[,-c(1,2)])
min_weights <- 0
max_weights <- max(weights_mx)
maxs <- rep(max_weights,)
mins <- apply(weights_mx, 2, min)
weights_scale <- scale(weights_mx,center = rep(min_weights,15),scale = rep(max_weights-min_weights,15))
weights_df <- cbind(dis_dis_edges[,c(1,2)],weights_scale)

i=1
for(i in 1:8){
  yr <- as.character(years[i])
  ed <- edlist[[i]]
  g <- graph_from_data_frame(ed,directed = F)
  md <- md[match(V(g)$name,md$Id),]
  wds <- wdegs_df %>%
    select(1,yr)
  wds <- wds[match(V(g)$name,wds$dis),]
  wd_n <- wds[,2]
  g <- set_vertex_attr(graph = g,name = "class",index = V(g),md$Class)
  g <- set_vertex_attr(graph = g,name = "wdeg",index = V(g),wd_n)
  
  p <- ggraph(g, layout = 'fr') + 
        geom_edge_link(aes(alpha=weight),color="green")+
        geom_node_point(aes(size = wdeg,color=class)) + 
        theme(legend.position = 'bottom')+
        scale_size_continuous(limits = c(min_wd,max_wd))
  print(p)
}


