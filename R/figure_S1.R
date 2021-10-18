#With this scritp, you can reproduce the analysis and plot the images that
#compose panel A-B in figure S1 of the paper

#load necessary packages
#file management
library(data.table)
#data manipulation
library(dplyr)
library(tidyverse)
library(reshape2)
#plotting
library(ggplot2)

options(stringsAsFactors = F)

#load data and create necessary objects
data("all_edges")
all_edges <- do.call(rbind,edges_list)
dis_dis_edges <- readRDS("intermediate/dis_dis_edges.RDS")
dis_dis_nodes <- as.data.frame(fread("data/dis_dis_nodes.csv"))
dis_class_df <- dis_dis_nodes %>%
  dplyr::select(1,3) %>%
  dplyr::rename(disease=1,class=2)
diseases <- unique(dis_class_df$disease)
years <- seq(1990,2018)
top_9_df <- readRDS("intermediate/top9_diseases_df.RDS")
top_9 <- top_9_df$disease
colors <- readRDS("intermediate/colors.RDS")

####Figure S1 A and B####
#Find closest diseases from each pair of classes
df <- dis_dis_edges %>%
  filter(dis_1 != dis_2) %>%
  left_join(dis_dis_nodes[,-1],by=c("dis_1"="Label")) %>%
  left_join(dis_dis_nodes[,-1],by=c("dis_2"="Label")) %>% 
  mutate(comparison=paste0(Class.x,"_x_",Class.y),
         same_class=ifelse(Class.x==Class.y,
                           "within_classes",
                           "between_classes"))

#create melt dataframe to plot in ggplot2
df_melt <- reshape2::melt(df)
df_melt$variable <- as.numeric(as.character(df_melt$variable))
df_melt$edge <- paste0(df_melt$Source,"_",df_melt$Target)

#plot distance between classes through the years
#distance between and within classes
colors_plot <- c("#cd0000","#006400","#551a8b",colors)
names(colors_plot) <- unique(df_melt$comparison)[c(2,5,3,1,6,4)]

p <- df_melt %>%
  filter(dis_1 %in% top_9,dis_2 %in% top_9) %>%
  #filter(Class.x != Class.y) %>%
  ggplot(aes(x=variable,y=value,color=comparison))+
  geom_point(size=.1,alpha=.2)+
  geom_smooth()+
  theme_minimal()+
  scale_color_manual(values = colors_plot)+
  scale_x_continuous(name = "year",
                     limits = c(1990,2018),
                     breaks = c(1990,2000,2010,2018))+
  scale_y_continuous(name = "similarity")+
  facet_wrap(facets = ~same_class,scales = "fixed")

pdf(file = "figures/figure_S1/between_and_within_class_similarity_evolution.pdf",
    width = 8,height = 2.7)
print(p)
dev.off()

####Figure S1C - number of paper of diseases with <100 or >100 genes####
dis_MeSH_all <- readRDS("data/dis_MeSH_all.rds")

#calculate number of genes per disease in 2018
edges_2018 <- all_edges %>%
  dplyr::filter(Doc.2018>0) %>%
  dplyr::select(1,2)

genes_per_disease <- as.data.frame(table(edges_2018$Target)) %>%
  rename(disease=1,genes=2) %>%
  left_join(dis_MeSH_all,by="disease")

mesh_terms <- unique(genes_per_disease$MESHTERM)

#search queries in PubMed using easyPubMed packege
count_vec <- c()
for(i in 1:length(mesh_terms)){
  my_query <- mesh_terms[i]
  my_query <- get_pubmed_ids(as.character(my_query),
                             api_key = "42f9ff230d26b4a045680973c38549d5e708")
  count_vec <- c(count_vec,my_query$Count)
}

papers_per_disease <- data.frame(MESHTERM=mesh_terms,
                                 papers=as.numeric(count_vec))

saveRDS(papers_per_disease,file = "intermediate/papers_per_disease_2018.RDS")
papers_per_disease <- readRDS("intermediate/papers_per_disease_2018.RDS")

genes_papers_per_disease <- genes_per_disease %>%
  left_join(papers_per_disease,by="MESHTERM") %>%
  mutate(number_genes=ifelse(genes>=100,">100","<100"))

p <- genes_papers_per_disease %>%
  ggplot(aes(x=number_genes,y=papers))+
  #geom_boxplot(outlier.color = NA)+
  geom_violin()+
  geom_jitter(aes(color=class),width = .2)+
  theme_minimal()+
  scale_color_manual(values = c("#68ad36","#ff743e","#00b7da"))+
  scale_y_continuous(labels = scales::comma)
pdf(file = "figures/figure_S1/genes_vs_papers_2018.pdf",
    width = 3.23,height = 4)
print(p)
dev.off()

papers_100_more <- genes_papers_per_disease %>%
  filter(number_genes == ">100") %>%
  pull(papers) %>% as.numeric()

paper_100_less <- genes_papers_per_disease %>%
  filter(number_genes == "<100") %>%
  pull(papers) %>% as.numeric()

t.test(x = paper_100_less,y = papers_100_more)

#Make correlation between number of papers and number of genes for each disease
#(new figure S1C after iScience revision)

p <- genes_papers_per_disease %>%
  group_by(disease,class) %>%
  summarise(genes=mean(genes),papers=sum(papers)) %>%
  ungroup() %>%
  ggplot(aes(x=papers,y = genes))+
  geom_point(aes(color=class))+
  geom_smooth(method = "glm")+
  theme_minimal()+
  scale_color_manual(values = c("#68ad36","#ff743e","#00b7da"))+
  scale_x_continuous(labels = scales::comma)
  facet_wrap(facets = ~class)
pdf(file = "figures/figure_S1/genes_vs_papers_2018_correlation.pdf",
    width = 4.7,height = 3.7)
print(p)
dev.off()

cor.df <- genes_papers_per_disease %>%
  group_by(disease,class) %>%
  summarise(genes=mean(genes),papers=sum(papers)) %>%
  ungroup()

cor.test(x = cor.df$papers,
         y = cor.df$genes)
