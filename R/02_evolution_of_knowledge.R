#script to analyze the evolution of knowledge and produce figures for the
#first section of results (Figure 01A-D)

#libraries
library(dplyr)
library(reshape2)
library(ggplot2)

#set working directory
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#load relevant files
data("all_edges")

#define colors that will be used for each disease class
colors <- c("#68ad36", "#ff743e", "#00b7da")

#get edges files
all_edges <- bind_rows(edges_list)
inflammatory_edges <- edges_list[[1]]
infectious_edges <- edges_list[[2]]
psychiatric_edges <- edges_list[[3]]

#get disease names
diseases <- data.frame(disease=unique(all_edges$Target))
diseases$class <- c(rep("inflammatory",27),rep("infectious",63),rep("psychiatric",9))

inflammatory_diseases <- diseases$disease[diseases$class=="inflammatory"]
infectious_diseases <- diseases$disease[diseases$class=="infectious"]
psychiatric_diseases <- diseases$disease[diseases$class=="psychiatric"]

#####calculate the number of genes per year per disease and class#####
for (i in 1:length(diseases$disease)){
  dis <- diseases$disease[i]
  class <- diseases$class[i]
  cols <- colnames(all_edges)[-c(1,2)]
  df <- all_edges %>%
    filter(Target==dis)
  for (j in 1:29){
    year <- as.numeric(gsub(x = cols[j],pattern = "Doc.",replacement = ""))
    col <- which(colnames(all_edges)==cols[j])
    genes <- df %>%
      select(1,2,col) %>%
      rename(docs=3) %>%
      filter(docs > 0) %>%
      pull(Source) %>%
      unique() %>%
      length()
    gen_y <- data.frame(dis=dis,genes=genes)
    colnames(gen_y) <- c("dis",year)
    if (j==1){
      genes_year <- gen_y
    }else{
      genes_year <- cbind(genes_year,gen_y[,-1])
      colnames(genes_year) <- c(colnames(genes_year)[1:j],year)
    }
  }
  if(i==1){
    genes_per_year_all <- genes_year
  }else{
    genes_per_year_all <- as.data.frame(rbind(genes_per_year_all,genes_year))
  }
}

save(genes_per_year_all,file="data/genes_per_year.RData")
write.csv(genes_per_year_all,file = "data/genes_per_year.csv",quote = F,row.names = F)

genes_per_year_all_melt <- melt(genes_per_year_all)
genes_per_year_all_melt$variable <- as.numeric(as.character(genes_per_year_all_melt$variable))
genes_per_year_all_melt$class <- c(rep("inflammatory",27),
                                   rep("psychiatric",9),
                                   rep("infectious",63))

colnames(genes_per_year_all_melt) <- c("dis","year","genes","class")

#####plot results of genes per year#####
#line graph of genes per year in each disease class
svg(filename = "figures/genes_per_year_lines.svg")
ggplot(genes_per_year_all_melt,aes(x=year,y=genes))+
  stat_summary(aes(color=class),fun.y = sum,
               geom = "line",size=1,show.legend = T)+
  scale_x_discrete(limits=c(1990,2000,2010,2018))+
  scale_color_manual(values = colors)+
  theme_minimal()
dev.off()

#boxplot of genes per disease in each disease class in 2018
genes_2018 <- genes_per_year_all_melt %>%
  filter(year==2018)

svg(filename = "figures/genes_per_class_boxplot.svg")
ggplot(genes_2018,aes(x=class,y=genes))+
  geom_boxplot(aes(fill=class),outlier.shape = NA,show.legend = F)+
  geom_point(aes(fill=class),shape = 21, colour = "black",show.legend = F)+
  scale_fill_manual(values=colors)+
  theme_minimal()
dev.off()

#violin plot of genes per disease in each disease class
svg(filename = "figures/genes_per_class_violin.svg")
ggplot(genes_2018,aes(x=class,y=genes))+
  geom_violin(aes(fill=class),show.legend = F)+
  geom_jitter(aes(fill=class),width = .3,shape=21,colour="black",show.legend = F)+
  scale_fill_manual(values=colors)+
  theme_minimal()
dev.off()

#####calculate the yearly variation of genes#####
gene_variation <- as.data.frame(t(apply(genes_per_year_all[,-1], 1, diff)))
rownames(gene_variation) <- genes_per_year_all$dis
gene_variation$disease <- rownames(gene_variation)
gene_variation <- gene_variation[,c(29,1:28)]

save(gene_variation,file = "data/gene_variation.RData")
write.csv(gene_variation,file = "data/gene_variation.csv",quote = F,row.names = F)

gene_variation_matrix <- as.matrix(gene_variation[,-1])
rownames(gene_variation_matrix) <- diseases$disease

gene_variation_melt <- melt(gene_variation)
gene_variation_melt$variable <- as.numeric(as.character(gene_variation_melt$variable))
gene_variation_melt$dis_class <- c(rep("inflammatory",27),
                               rep("psychiatric",9),
                               rep("infectious",63))

#####plot the yearly variation of genes in each class#####
genes_each_year_dis_class <- gene_variation_melt %>%
  group_by(dis_class,variable) %>%
  summarise(sum=sum(value))
genes_each_year_dis_class$variable <- as.numeric(as.character(genes_each_year_dis_class$variable))
#remove the year 2018. Searches were made June-August of 2018 so data for
#the year is incomplete
genes_each_year_dis_class <- genes_each_year_dis_class[!genes_each_year_dis_class$variable==2018,]

#genes each year line plot and ridge plot
#line plot
svg(filename = "figures/genes_each_year_line.svg",width = 12,height = 4)
ggplot(genes_each_year_dis_class,aes(x=variable,y = sum,color=dis_class))+
  geom_line()+
  facet_wrap(facets = as.factor(genes_each_year_dis_class$dis_class))+
  scale_color_manual(values = colors)+
  theme_minimal()
dev.off()

#ridge plot
svg(filename = "figures/genes_each_year_ridgeline.svg",width = 12,height = 4)
ggplot(genes_each_year_dis_class,aes(x=variable,y=dis_class,fill=dis_class,color=dis_class))+
  geom_ridgeline(aes(height=sum),scale=.005,alpha=0.274,size=.2)+
  scale_fill_manual(values=c("#68ad36", "#ff743e", "#00b7da"))+
  scale_color_manual(values=c("#68ad36", "#ff743e", "#00b7da"))+
  scale_x_discrete(limits=c(1991,2000,2010,2017))+
  theme_ridges(font_size = 10)+
  theme(legend.position = 'none')
dev.off()

#####plot the yearly variation of genes in all and selected diseases#####
#all diseases by class (separate files)
dis_classes <- c("infectious","inflammatory","psychiatric")
for (i in 1:3){
  class <- dis_classes[i]
  
  #all one file per class
  scales <- c(0.1,0.05,0.05)
  heights <- c(10,6.5,5)
  filename <- paste0("figures/genes_each_year_",class,"_diseases.svg")
  svg(filename,width = heights[i]+(heights[i]/2),height = heights[i])
  p <- gene_variation_melt %>%
    filter(dis_class==class) %>%
    ggplot(aes(x=variable,y=disease,fill=dis_class,color=dis_class))+
    geom_ridgeline(aes(height=value),scale=scales[i],alpha=0.274,color="transparent")+
    scale_fill_manual(values=colors[i])+
    scale_x_discrete(limits=c(1991,2000,2010,2017))+
    theme_ridges(font_size = 10)+
    theme(legend.position = 'none')
  print(p)
  dev.off()

  #top 9
  top_9 <- gene_variation_melt %>%
    filter(dis_class==class) %>%
    arrange(desc(value)) %>%
    pull(disease) %>%
    unique() %>%
    .[1:9]
  
  order <- gene_variation_melt %>%
    filter(dis_class==class & disease %in% top_9) %>%
    dcast(formula = disease~variable,value.var = "value") %>%
    .[,-1] %>% as.matrix() %>% dist() %>% hclust() %>%
    .$order
  
  filename <- paste0("figures/genes_each_year_",class,"_top_9.svg")
  svg(filename,width = heights[i]+(heights[i]/2),height = heights[i])
  scales <- c(0.03,0.03,0.03)
  p <- gene_variation_melt %>%
    filter(dis_class==class & disease %in% top_9) %>%
    ggplot(aes(x=variable,y=disease,fill=dis_class,color=dis_class))+
    geom_ridgeline(aes(height=value),scale=scales[i],alpha=0.274)+
    scale_y_discrete(limits=top_9[order])+
    scale_x_discrete(limits=c(1991,2000,2010,2017))+
    scale_color_manual(values = colors[i])+
    scale_fill_manual(values = colors[i])+
    theme_ridges(font_size = 10)+
    theme(legend.position = 'none')
  print(p)
  dev.off()
}

#all in one file
filename <- paste0("figures/genes_each_year_all_diseases.svg")
order <- diseases$disease
svg(filename,width = 20,height = 20)
gene_variation_melt %>%
  ggplot(aes(x=variable,y=disease,fill=dis_class,color=dis_class))+
  geom_ridgeline(aes(height=value),scale=0.05,alpha=0.274,color="transparent")+
  scale_fill_manual(values=colors)+
  scale_y_discrete(limits=order)+
  scale_x_discrete(limits=c(1991,2000,2010,2017))+
  theme_ridges(font_size = 12)+
  theme(legend.position = 'none')
dev.off()

