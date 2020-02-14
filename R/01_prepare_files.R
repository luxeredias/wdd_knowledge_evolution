#script to manipulate evolution files and create Gephi files
#for inflammatory, infectious and psychiatric diseases
library(dplyr)
options(stringsAsFactors = F)

#set working directory
dir <- "~/Área de Trabalho/GitHub_projects/wdd_knowledge_evolution/"
setwd(dir)

#load edges files
#inflammatory diseases
inflammatory_edges <- read.csv("data/inflammatory_edges.csv")
inflammatory_edges <- inflammatory_edges[,c(1,2,25:53)]

inflammatory_edges_gephi <- inflammatory_edges[,c(1,2,seq(3,31,2))]
inflammatory_nodes_gephi <- data.frame(Id=c(unique(inflammatory_edges_gephi$Source),unique(inflammatory_edges_gephi$Target)),
                                       Label=c(unique(inflammatory_edges_gephi$Source),unique(inflammatory_edges_gephi$Target)),
                                       Class=c(rep("GENE",length(unique(inflammatory_edges_gephi$Source))),
                                               rep("CONDITION",length(unique(inflammatory_edges_gephi$Target)))),
                                       Size=c(rep(1,length(unique(inflammatory_edges_gephi$Source))),
                                              rep(10,length(unique(inflammatory_edges_gephi$Target)))),
                                       dis_class="INFLAMMATORY")

write.csv(inflammatory_edges_gephi,"data/inflammatory_edges_gephi.csv",quote = F,row.names = F)
write.csv(inflammatory_nodes_gephi,"data/inflammatory_nodes_gephi.csv",quote = F,row.names = F)

#infectious diseases
infectious_edges <- read.csv("data/infectious_edges.csv")
infectious_edges_gephi <- infectious_edges

#change names with "," to avoid errors in Gephi
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="INFECTION, MYCOBACTERIUM"] <- "MYCOBACTERIUM INFECTION"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="INFECTIONS, TOXOPLASMA GONDII"] <- "TOXOPLASMA GONDII INFECTION"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="FEVER, YELLOW"] <- "YELLOW FEVER"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="INFECTION, YERSINIA"] <- "YERSINIA INFECTION"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="INFECTION, LISTERIA"] <- "LISTERIA INFECTION"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="FASCIITIDES, NECROTIZING"] <- "NECROTIZING FASCIITIDES"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="INFECTION, STREPTOCOCCAL"] <- "STREPTOCOCCAL INFECTION"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="DISEASE, FUNGUS"] <- "FUNGUS DISEASE"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="FLU, HUMAN"] <- "HUMAN FLU"
infectious_edges_gephi$Target[infectious_edges_gephi$Target=="INFECTION, SALMONELLA"] <- "SALMONELLA INFECTION"

infectious_edges <- infectious_edges_gephi
infectious_edges_gephi <- infectious_edges_gephi[,c(1,2,seq(3,31,2))]

infectious_nodes_gephi <- data.frame(Id=c(unique(infectious_edges_gephi$Source),unique(infectious_edges_gephi$Target)),
                                     Label=c(unique(infectious_edges_gephi$Source),unique(infectious_edges_gephi$Target)),
                                     Class=c(rep("GENE",length(unique(infectious_edges_gephi$Source))),
                                             rep("CONDITION",length(unique(infectious_edges_gephi$Target)))),
                                     Size=c(rep(1,length(unique(infectious_edges_gephi$Source))),
                                            rep(10,length(unique(infectious_edges_gephi$Target)))),
                                     dis_class="INFECTIOUS")

write.csv(infectious_edges_gephi,"data/infectious_edges_gephi.csv",quote = F,row.names = F)
write.csv(infectious_nodes_gephi,"data/infectious_nodes_gephi.csv",quote = F,row.names = F)

#psychiatric disorders
psychiatric_edges <- read.csv("data/psychiatric_edges.csv")
psychiatric_edges <- psychiatric_edges[,c(1,2,25:53)]

#change names with "," to avoid errors in Gephi
psychiatric_edges$Target[psychiatric_edges$Target=="DISORDER, SCHIZOPHRENIC"] <- "SCHIZOPHRENIA"

psychiatric_edges_gephi <- psychiatric_edges[,c(1,2,seq(3,31,2))]
psychiatric_nodes_gephi <- data.frame(Id=c(unique(psychiatric_edges_gephi$Source),unique(psychiatric_edges_gephi$Target)),
                                      Label=c(unique(psychiatric_edges_gephi$Source),unique(psychiatric_edges_gephi$Target)),
                                      Class=c(rep("GENE",length(unique(psychiatric_edges_gephi$Source))),
                                              rep("CONDITION",length(unique(psychiatric_edges_gephi$Target)))),
                                      Size=c(rep(1,length(unique(psychiatric_edges_gephi$Source))),
                                             rep(10,length(unique(psychiatric_edges_gephi$Target)))),
                                      dis_class="PSYCHIATRIC")

write.csv(psychiatric_edges_gephi,"data/psychiatric_edges_gephi.csv",quote = F,row.names = F)
write.csv(psychiatric_nodes_gephi,"data/psychiatric_nodes_gephi.csv",quote = F,row.names = F)

#make a list of processed edges files and save in R object
edges_list <- list(inflammatory_edges,psychiatric_edges,infectious_edges)
save(edges_list,file = "data/all_edges.RData")

all_edges <- bind_rows(inflammatory_edges,infectious_edges,psychiatric_edges)
diseases <- data.frame(disease=unique(all_edges$Target))
diseases$class <- c(rep("inflammatory",27),rep("infectious",63),rep("psychiatric",9))

all_nodes <- data.frame(Id=c(unique(all_edges$Source),unique(all_edges$Target)))
all_nodes$Label <- c(rep(" ",3723),all_nodes$Id[3724:3822])
all_nodes$Class <- c(rep("GENE",3723),rep("CONDITION",99))
all_nodes$Size <- c(rep(1,3723),rep(10,99))
all_nodes$Dis_class <- c(rep(" ",3723),rep("INFLAMMATORY",27),rep("INFECTIOUS",63),rep("PSYCHIATRIC",9))

save(diseases,file = "data/diseases.RData")
write.csv(all_edges,"data/all_edges.csv",quote = F,row.names = F)
write.csv(all_nodes,"data/all_nodes.csv",quote = F,row.names = F)
