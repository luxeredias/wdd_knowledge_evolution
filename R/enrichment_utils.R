#function to run functional enrichment against Reactome for genes in Evolution
#network

#make a list of genes per disease per year (remove relationships with zero docs)
make_gene_lists <- function(edges){
  diseases <- unique(edges$Target)
  gene_list <- sapply(diseases,simplify = F,function(dis){
    ed <- edges %>%
      dplyr::filter(Target==dis) %>%
      dplyr::select(-2) %>%
      dplyr::group_by(Source) %>%
      dplyr::summarise_all(sum) %>%
      dplyr::ungroup() %>%
      column_to_rownames("Source") %>%
      as.matrix()
    
    dis_genes_year <- apply(ed,MARGIN = 2,function(yr){
      yr <- yr[yr>0]
      return(names(yr))
    })
  })
}

# gene_list <- make_gene_lists(edges = all_edges_converted)
# gmt <- "~/Área de Trabalho/Papers_Helder/COVID-19_genome_expression/data/c2.cp.reactome.v7.0.symbols.gmt"

run_enrichment <- function(gene_list,gmt){
  require(clusterProfiler)
  require(stringr)
  
  #load gmt
  gmt <- read.gmt(gmt)
  
  #edit gmt to make it nice
  gmt$ont <- as.character(gmt$ont)
  gmt$ont <- stringr::str_replace_all(stringr::str_remove(gmt$ont,"REACTOME_"),"_"," ")
  gmt$ont <- tolower(gmt$ont)
  substr(gmt$ont,1,1) <- toupper(substr(gmt$ont,1,1))
  
  enrichment <- lapply(gene_list,function(dis){
    lapply(dis,function(year_genes){
      enr <- clusterProfiler::enricher(year_genes,
                                       pAdjustMethod = "BH",
                                       TERM2GENE = gmt)
      enr_result <- ifelse(test = class(enr)=="NULL",
                           yes="NULL",
                           no=enr@result)
    })
  })
}

# enr_out <- run_enrichment(gene_list,gmt = gmt)
# save(enr_out,file="~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/data/enrichment_years_new.RData")
load("data/enrichment_years_new.RData")

z<-Filter(function(x) x$p.adjust < 0.05,enr_out)
enrichment_list <- enr_out
padj=0.05
make_network_from_enrichment <- function(enrichment_list,
                                         padj=0.05){
  dis <- enrichment_list[[1]]
  year_enr <- dis[[28]]
  enr_df <- lapply(enrichment_list,function(dis){
    lapply(dis,function(year_enr) {
      if(year_enr=="NULL"){
        enr <- "NULL"
      }else{
        enr <- year_enr %>% filter(p.adjust < 0.05)
        
      }
    })
  })
}



