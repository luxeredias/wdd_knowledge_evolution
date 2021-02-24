#function to run functional enrichment against Reactome for genes in Evolution
#network

#make a list of genes per disease per year (remove relationships with zero docs)
edges <- all_edges_converted
years <- seq(1990,2018,1)

make_term_dis_network <- function(edges,years,gmt,enr_pval){
  diseases <- unique(edges$Target)
  univ <- unique(edges$Source)
  for (i in 1:length(diseases)){
    for(j in 1:length(years)){
      
      dis <- diseases[i]
      yr <- years[j]
      
      col <- which(grepl(yr,colnames(edges)))
      
      ed <- edges %>%
        dplyr::filter(Target==dis) %>%
        dplyr::select(1,2,col) %>%
        dplyr::rename(docs=3) %>%
        dplyr::filter(docs > 0)
      genes <- unique(ed$Source)
      enr <- clusterProfiler::enricher(genes,pAdjustMethod = "BH",
                                       TERM2GENE = gmt,universe = univ)
      if(class(enr)=="NULL"){
        term_dis_net <- data.frame(ID="",disease=dis,yr="")
        colnames(term_dis_net)[3] <- yr
        if(j==1){
          term_dis_net_a <- term_dis_net
        }else{
          term_dis_net_a <- term_dis_net_a %>%
            full_join(term_dis_net,by=c("ID","disease"))
        }
      }else{
        enr_result <- enr@result
        term_dis_net <- enr_result %>%
          dplyr::mutate(disease=dis) %>%
          dplyr::filter(pvalue < enr_pval) %>%
          dplyr::select(ID,disease,pvalue)
        colnames(term_dis_net)[3] <- yr
        if(j==1){
          term_dis_net_a <- term_dis_net
        }else{
          term_dis_net_a <- term_dis_net_a %>%
            full_join(term_dis_net,by=c("ID","disease"))
        }
      }
    }
    if(i==1){
      term_dis_net_all <- term_dis_net_a
    }else{
      term_dis_net_all <- rbind(term_dis_net_all,term_dis_net_a)
    }
  }
  term_dis_net_all <- term_dis_net_all %>%
    filter(ID!="")
  return(term_dis_net_all)
}  

#save(term_dis_net_all_2,file = "data/terms_dis_net_all.RData")

years <- seq(1990,2018,1)
plot_term_dis_network <- function(term_dis_net,years){
  require(igraph)
  for(i in 1:length(years)){
    yr <- years[i]
    col <- which(colnames(term_dis_net)==as.character(yr))
    df <- term_dis_net %>%
      dplyr::select(1,2,col) %>%
      dplyr::rename(from=1,to=2,weight=3) %>%
      dplyr::filter(!is.na(weight)) %>%
      dplyr::mutate(weight=as.numeric(weight))
    g <- igraph::graph_from_data_frame(df,directed = F)
    
  }
}

