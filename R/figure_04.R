#With this scritp, you can reproduce the analysis and plot the images that
#compose panel A-E in figure 04 of the paper

#necessary packages
library(dplyr)
library(igraph)
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
library(data.table)
library(reshape2)
library(ggraph)
library(ggrepel)
options(stringsAsFactors = F)

#load necessary data
data("all_edges")
all_edges <- do.call(rbind,edges_list)
load("data/gene_convertion_table.RData")
all_edges_converted <- all_edges %>%
  dplyr::left_join(gene_convertion_table,by=c("Source"="alias")) %>%
  dplyr::select(-1) %>%
  dplyr::select(converted,Target,everything())
dis_dis_nodes <- as.data.frame(fread("data/dis_dis_nodes.csv"))
dis_class_df <- dis_dis_nodes %>%
  dplyr::select(1,3) %>%
  dplyr::rename(disease=1,class=2)
diseases <- unique(all_edges$Target)
top_9_df <- readRDS("intermediate/top9_diseases_df.RDS")
top_9 <- top_9_df$disease

#get the genes that are shared between all disease categories
edges_2018 <- all_edges_converted %>%
  dplyr::filter(Doc.2018 > 0) %>%
  dplyr::select(converted,Target) %>%
  dplyr::rename(from=1,to=2) %>%
  filter(!is.na(from))

nodes <- data.frame(Id=unique(c(edges_2018$from,
                                edges_2018$to)),
                    Type=c(rep("gene",length(unique(edges_2018$from))),
                           rep("diseases",length(unique(edges_2018$to)))))

classes <- unique(dis_class_df$class)
for(i in 1:3){
  cls <- classes[i]
  
  edges_class <- edges_2018 %>%
    left_join(dis_class_df,by=c("to"="disease")) %>%
    filter(class==cls) %>%
    dplyr::select(-3)
  
  nodes_class <- nodes %>%
    filter(Id %in% unique(c(edges_class$from,edges_class$to)))
  
  g <- igraph::graph_from_data_frame(edges_class,directed = F)
  degree <- as.data.frame(igraph::degree(g)) %>%
    rownames_to_column("Id") %>%
    rename(degree=2)
  
  nodes_class_2 <- nodes_class %>%
    left_join(degree,by = "Id")
  
  top_genes <- nodes_class_2 %>%
    filter(Type=="gene") %>%
    top_n(wt = degree,n = 100) %>%
    mutate(dis_class=cls)
  if(i==1){
    top_genes_all <- top_genes
  }else{
    top_genes_all <- rbind(top_genes_all,top_genes)
  }
}

top_genes_list <- split(top_genes_all$Id,top_genes_all$dis_class)
class_genes <- edges_2018 %>%
  left_join(dis_class_df,by=c("to"="disease")) %>%
  dplyr::select(-2) %>%
  rename(to=2) %>%
  filter(!duplicated(.))
class_genes_list <- split(class_genes$from,class_genes$to)
names(class_genes_list) <- c("Infectious","Inflammatory","Psychiatric")

####Figure 3A####
pdf("figures/figure_04/panel_A/upset_shared_genes_category.pdf",width = 4.55,height = 3.5)
UpSetR::upset(UpSetR::fromList(class_genes_list),
              sets.bar.color = c("#68ad36","#00b7da","#ff743e"),
              mainbar.y.label = "genes",sets.x.label = "class genes")

dev.off()

g2 <- graph_from_data_frame(d = class_genes,directed = F)
degree_classes <- as.data.frame(igraph::degree(g2)) %>%
  rownames_to_column("Id") %>%
  rename(degree=2)

genes_in_all_classes <- degree_classes %>%
  filter(degree==3)

#load drug-diseases and drug-gene information from CTD (download files)
download.file(url = "http://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz",
              destfile = "intermediate/CTD_chemicals_diseases.csv.gz")
CTD_drug_disease <- fread("intermediate/CTD_chemicals_diseases.csv.gz")
colnames(CTD_drug_disease) <- c("ChemicalName","ChemicalID","CasRN","DiseaseName",
                                "DiseaseID","DirectEvidence","InferenceGeneSymbol",
                                "InferenceScore","OmimIDs","PubMedIDs")

download.file(url = "http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz",
              destfile = "intermediate/CTD_chem_gene_ixns.csv.gz")
CTD_drug_gene <- fread("intermediate/CTD_chem_gene_ixns.csv.gz")
colnames(CTD_drug_gene) <- c("ChemicalName","ChemicalID","CasRN","GeneSymbol",
                             "GeneID","GeneForms","Organism","OrganismID","Interaction",
                             "InteractionActions","PubMedIDs")

#load diseases and MeSH IDs
diseases_MeSH <- readRDS("data/dis_MeSH_all.rds")
diseases_MeSH <- diseases_MeSH %>%
  separate_rows(MeSH,sep = "\\|") %>%
  filter(!is.na(MeSH))

all_MeSH <- paste0(diseases_MeSH$MeSH,collapse = "|")

#get drug-disease interactions for the top9 diseases in each category
top_9_MeSH <- diseases_MeSH %>%
  filter(disease %in% top_9)

top_9_MeSH_search <- paste0(top_9_MeSH$MeSH,collapse = "|")

drug_disease_top_9 <- CTD_drug_disease %>%
  filter(grepl(top_9_MeSH_search,DiseaseID),
         DirectEvidence=="therapeutic")

drug_top_genes <- CTD_drug_gene %>%
  filter(GeneSymbol %in% genes_in_all_classes$Id,
         OrganismID==9606) %>%
  filter(!duplicated(paste0(ChemicalID,GeneSymbol))) %>%
  filter(ChemicalID %in% unique(drug_disease_top_9$ChemicalID))

drug_top_genes_edges <- drug_top_genes %>%
  dplyr::select(ChemicalName,GeneSymbol) %>%
  dplyr::rename(from=1,to=2)

drug_top_genes_nodes <- data.frame(node=unique(c(drug_top_genes_edges$from,
                                                 drug_top_genes_edges$to)),
                                   type=c(rep("drug",length(unique(drug_top_genes_edges$from))),
                                          rep("gene",length(unique(drug_top_genes_edges$to)))))

ngenes_per_drug <- as.data.frame(table(drug_top_genes$ChemicalName)) %>%
  rename(ChemicalName=1,ngenes=2)

####Figure 3D####
p <- ngenes_per_drug %>%
  top_n(wt = ngenes,n = 20) %>%
  arrange(desc(ngenes)) %>%
  mutate(ChemicalName=factor(ChemicalName,levels = rev(ChemicalName))) %>%
  ggplot(aes(y=ChemicalName,x=ngenes))+
  geom_segment(aes(x=0, xend=ngenes, y=ChemicalName, yend=ChemicalName), color="skyblue")+
  geom_point(aes(size=ngenes))+
  theme_minimal()
pdf(file = "figures/figure_04/panel_D/hub_genes_per_drug.pdf",width = 4.5,height = 4.2)
print(p)
dev.off()

ndrugs_per_gene <- as.data.frame(table(drug_top_genes$GeneSymbol)) %>%
  rename(GeneSymbol=1,ndrugs=2)

####Figure 3B####
p <- ndrugs_per_gene %>%
  top_n(wt = ndrugs,n = 20) %>%
  arrange(desc(ndrugs)) %>%
  mutate(GeneSymbol=factor(GeneSymbol,levels = rev(GeneSymbol))) %>%
  ggplot(aes(y=GeneSymbol,x=ndrugs))+
  geom_segment(aes(x=0, xend=ndrugs, y=GeneSymbol, yend=GeneSymbol), color="pink")+
  geom_point(aes(size=ndrugs))+
  theme_minimal()
pdf(file = "figures/figure_04/panel_B/drugs_per_hub_gene.pdf",width = 4.5,height = 4.2)
print(p)
dev.off()

####Figure 3E####
top_20_drugs <- ngenes_per_drug %>%
  top_n(wt = ngenes,n = 20)
top_20_genes <- ndrugs_per_gene %>%
  top_n(wt = ndrugs,n = 20)

top_20_drugs_genes <- drug_top_genes_edges %>%
  filter(from %in% top_20_drugs$ChemicalName,
         to %in% top_20_genes$GeneSymbol)

g_top_drugs <- graph_from_data_frame(top_20_drugs_genes,directed = F)
V(g_top_drugs)$degree <- degree(g_top_drugs)
V(g_top_drugs)$type <- c(rep("drug",20),rep("gene",20))
p <- ggraph(graph = g_top_drugs,layout = "kk")+
  geom_edge_link(color="grey90")+
  geom_node_point(aes(size=degree,fill=type),color="black",pch=21)+
  geom_node_text(aes(label = name),repel=TRUE)+
  theme_void()
# theme(legend.position = "none")
pdf(file = "figures/figure_04/panel_E/drug_hub_network.pdf",width = 7.62,height = 5.95)
print(p)
dev.off()

#make timeline of top20 genes per category
top20_timeline <- all_edges_converted %>%
  dplyr::filter(converted %in% as.character(top_20_genes$GeneSymbol)) %>%
  dplyr::left_join(dis_class_df,by=c("Target"="disease")) %>%
  dplyr::select(-Target) %>%
  dplyr::select(converted,class,everything()) %>%
  dplyr::group_by(converted,class) %>%
  dplyr::summarise_all(sum) %>%
  dplyr::ungroup()

gn_vec <- c()
cls_vec <- c()
yr_vec <- c()
for(i in 1:20){
  gn <- unique(top20_timeline$converted)[i]
  for(j in 1:3){
    cls <- classes[j]
    df <- top20_timeline %>%
      filter(converted==gn,
             class==cls) %>%
      melt()
    yr <- as.character(df$variable)[which.max(df$value > 0)]
    yr <- as.numeric(gsub(yr,pattern = "Doc.",replacement = ""))
    
    gn_vec <- c(gn_vec,gn)
    cls_vec <- c(cls_vec,cls)
    yr_vec <- c(yr_vec,yr)
  }
}

timeline_df <- data.frame(gene=gn_vec,class=cls_vec,year=yr_vec)
class_colors <- c("#68ad36","#00b7da","#ff743e")

####Figure 3C####
p <- timeline_df %>%
  arrange(desc(year)) %>%
  mutate(gene=factor(gene,levels = unique(gene))) %>%
  mutate(class=factor(class,levels = classes)) %>%
  ggplot(aes(x=year,y=gene,color=class))+
  geom_point(size=4)+
  scale_color_manual(values = rev(class_colors))+
  scale_x_continuous(breaks = 1990:2018)+
  #geom_text_repel(aes(label=gene),force_pull = -.1)+
  theme_minimal()
pdf(file = "figures/figure_04/panel_C/top20_genes_timeline_v2.0.pdf",width = 13.01,height = 4.39)
print(p)
dev.off()
