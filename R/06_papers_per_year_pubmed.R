#script to query pubmed and get number of papers per year per disease/gene
library(dplyr)
library(tidyverse)
library(ggplot2)
library(MeSH.db)
library(easyPubMed)
source("R/PubMedTrend.R")

#set working dir
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#get MeSH terms of diseases
load("data/diseases.RData")
dis_MeSH <- read.csv("data/diseases.csv")

MeSH_ids <- dis_MeSH %>%
  separate_rows(MeSH,sep = "\\|") %>%
  pull(MeSH) %>%
  unique()

MeSH_terms <- select(MeSH.db,keys = MeSH_ids,columns = c("MESHID","MESHTERM"),keytype = "MESHID")

dis_MeSH_all <- dis_MeSH %>%
  separate_rows(MeSH,sep = "\\|") %>%
  #filter(!is.na(MeSH)) %>%
  full_join(MeSH_terms,by=c("MeSH"="MESHID"))

dis_MeSH_all$MESHTERM[67] <- "Zika Virus Infection"
dis_MeSH_all$MESHTERM[95] <- "Autism Spectrum Disorder"

#search number of papers per year per disease using PubMedTrend.R
queries <- dis_MeSH_all$MESHTERM
names(queries) <- dis_MeSH_all$MESHTERM

# pubmed_trend <- PubMedTrend(query = queries,yrStart = 1990,yrMax = 2018)
# colnames(pubmed_trend) <- c("MESHTERM","year","count")

# save(pubmed_trend,file = "data/pubmed_trend_all_dis.RData")
load("data/pubmed_trend_all_dis.RData")
colnames(pubmed_trend) <- c("MESHTERM","year","count")

pubmed_trend <- pubmed_trend %>%
  filter(MESHTERM!="BONE INFECTION")

top_9_dis_class <- pubmed_trend %>%
  left_join(dis_MeSH_all,by="MESHTERM") %>%
  group_by(MESHTERM,class) %>%
  dplyr::summarise(count=sum(count))

p <- pubmed_trend %>%
  left_join(dis_MeSH_all,by="MESHTERM") %>%
  filter(class=="infectious") %>%
  ggplot(aes(x=year,y=count,color=MESHTERM))+
  geom_line()
p

dis_vec <- c()
max_vec <- c()
yr_max_vec <- c()
years <- c(1990:2018)
for(i in 1:length(dis_MeSH_all$MESHTERM[-25])){
  dis <- dis_MeSH_all$MESHTERM[-25][i]
  cnts <- pubmed_trend %>%
    filter(MESHTERM==dis) %>%
    as.data.frame()
  if(any(!years %in% cnts$year)){
    yr <- years[which(!years %in% cnts$year)]
    df <- data.frame(MESHTERM=dis,year=as.character(yr),count=0)
    cnts <- rbind(cnts,df)
  }
  dif_vec <- c()
  for(k in 1:length(c(1991:2018))){
    yr <- c(1991:2018)[k]
    cnt_yr_1 <- cnts$count[cnts$year==yr]
    cnt_yr_0 <- cnts$count[cnts$year==yr-1]
    dif <- cnt_yr_1 - cnt_yr_0
    names(dif) <- yr
    dif_vec <- c(dif_vec,dif)
  }
  
  max_dif <- max(dif_vec)
  yr_max_dif <- names(dif_vec[dif_vec==max_dif])
  if(length(yr_max_dif)>1){
    dis_vec <- c(dis_vec,rep(dis,length(yr_max_dif)))
    max_vec <- c(max_vec,rep(max_dif,length(yr_max_dif)))
    yr_max_vec <- c(yr_max_vec,yr_max_dif) 
  }else{
    dis_vec <- c(dis_vec,dis)
    max_vec <- c(max_vec,max_dif)
    yr_max_vec <- c(yr_max_vec,yr_max_dif)
  }
}

peak_df <- data.frame(disease=dis_vec,
                      peak=max_vec,
                      year_peak=as.numeric(yr_max_vec))

pubmed_trend_names <- pubmed_trend %>%
  left_join(dis_MeSH_all,by="MESHTERM") %>%
  dplyr::group_by(disease,year,class) %>%
  dplyr::summarise(count=sum(count)) %>%
  dplyr::select(disease,year,count,class) %>%
  dplyr::rename(papers=3)

gene_variation_melt_2 <- gene_variation_melt %>%
  dplyr::select(disease,variable,value,dis_class) %>%
  dplyr::rename(year=2,genes=3,class=4)

papers_vs_genes <- gene_variation_melt_2 %>%
  dplyr::full_join(pubmed_trend_names,by=c("disease","year","class"))

papers_vs_genes$genes[is.na(papers_vs_genes$genes)] <- 0
papers_vs_genes$papers[is.na(papers_vs_genes$papers)] <- 0

papers_vs_genes$gene_per_paper <- papers_vs_genes$genes/papers_vs_genes$papers
papers_vs_genes$gene_per_paper[is.infinite(papers_vs_genes$gene_per_paper)] <- 0

p <- papers_vs_genes %>%
  filter(disease=="PSORIASIS") %>%
  #filter(class=="psychiatric") %>%
  ggplot(aes(x=papers,y=genes,color=disease))+
  geom_point()+
  geom_smooth(method = lm, se = FALSE)
p

p <- papers_vs_genes %>%
  filter(disease=="PSORIASIS") %>%
  #filter(class=="psychiatric") %>%
  ggplot(aes(x=year,fill=disease))+
  geom_col(aes(y=papers),alpha=0.2)+
  geom_line(aes(y=gene_per_paper*100*100*10,color=disease))+
  scale_y_continuous(name = "Papers",
                     sec.axis = sec_axis(trans = ~./100000,name = "Genes per paper"))
p

p <- papers_vs_genes %>%
  filter(disease=="PSORIASIS") %>%
  #filter(class=="psychiatric") %>%
  ggplot(aes(x=year,fill=disease))+
  geom_col(aes(y=genes),alpha=0.2)+
  geom_line(aes(y=gene_per_paper*1000,color=disease))+
  scale_y_continuous(name = "Genes",
                     sec.axis = sec_axis(trans = ~./100000,name = "Genes per paper"))
p

p <- papers_vs_genes %>%
  filter(disease=="PSORIASIS") %>%
  #filter(class=="psychiatric") %>%
  ggplot(aes(x=year,fill=disease))+
  geom_col(aes(y=genes*100),alpha=0.2)+
  geom_col(aes(y=papers),alpha=0.2)+
  scale_y_continuous(name = "Papers",
                     sec.axis = sec_axis(trans = ~./100,name = "Genes"))
p
