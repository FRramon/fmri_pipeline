# data_path = "/Volumes/LaCie/"
#labels <- colnames(read.csv("/Users/francoisramon/Desktop/These/fMRI_pipeline/labels_schaefer.csv"))
# 
# 
# M <- read.csv(data_path,header=F)
# 
# 
# grep("Default",labels)

library(igraph)
library(NetworkConhect)
library(gplots)

# path <- "/Volumes/LaCie/derivatives/grouped_results"
# data_path = "/Volumes/LaCie/derivatives/grouped_results/FC_schaefer.csv"
# DF <- read.csv(data_path)
# colnames(DF) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
# DF <- DF[c("subject_id","visit_id","from","to","weight")]
# 
# 
# DF$subject_id <- as.character(DF$subject_id)
# DF$visit_id <- as.character(DF$visit_id)
# 
# DF$visit_id[DF$visit_id == "1"] <- "V1"
# DF$visit_id[DF$visit_id == "2"] <- "V2"
# DF$visit_id[DF$visit_id == "3"] <- "V3"
# 
# weight <- "PearsonCorrel"
# real_weight <- "PearsonCorrel"

makeGraph <- function(data,s_id,v_id,WM_metric=c('FBC','ODI','GFA','FA','Fintra','PearsonCorrel'),threshold = 0){
  WM_metric<-match.arg(WM_metric)
  df <- subset(data, subject_id == s_id & visit_id == v_id)
  edgelist <- df[, c('from','to','weight')]
  # edgelist<-data.frame(df[,3],df[,4],df[,5])
  # colnames(edgelist) <- c('from','to','weight')
  elist<-subset(edgelist,weight>threshold)
  #finv <- function(x) 1/x
  #elist$weight <- as.numeric(lapply(elist$weight,finv))
  g<-graph_from_data_frame(elist,directed=FALSE)
  #print(is.weighted(g))
  g
}

computeMetric <- function(data,
                      v_id,
                      WM_metric = c('FBC','ODI','GFA','Fintra','FA','PearsonCorrel'),
                      eval = c('clust_coeff','characteristic_path','global_eff','local_eff','smallworldeness','richcore','strength','betweenness'),
                      thresh_method,
                      tvalue,
                      rsnet = "All"
){
  #LrsNet <- c("Default","SomMot","SalVentAttn","Vis","Limbic","Cont","DorsAttn")
  #rsNet <- LrsNet[5]
  
  
  WM_metric<-match.arg(WM_metric)
  eval<-match.arg(eval)
  subject_ids <- get_subject_ids(data,v_id)
  RES <- vector("numeric", length(subject_ids)) # créer une nouvelle liste vide
  
  for(j in 1:length(subject_ids)){
    s_id <- subject_ids[j]
    if(thresh_method == "threshold"){
      g <- makeGraph(data,s_id,v_id,WM_metric,tvalue)
      if (rsnet != "All"){ 
        i_DMN <- grep(rsnet,V(g)$name)
        nodes_DMN <- V(g)$name[i_DMN]
        g <- induced_subgraph(
          g,
          i_DMN 
        )
      } 
    } else if(thresh_method == "density"){
      gT <- makeGraph(data,s_id,v_id,WM_metric,0)
      g <- sparseThresh(gT,tvalue)
      i_DMN <- grep(rsnet,V(g)$name)
      if (rsnet != "All"){ 
        nodes_DMN <- V(g)$name[i_DMN]
        g <- induced_subgraph(
          g,
          i_DMN)
      }
    }
    if(eval == 'clust_coeff'){
      value <- transitivity(g,"global")
      #value <- normalized_cluster_coeff(g)
    } else if (eval == 'characteristic_path'){
      #value <- mean_distance(g)
      value <- mean_distance(g,weights = 1/(E(subg)$weight))
      #value <- normalized_shortest_path(g)
    } else if (eval == 'global_eff'){
      # value <- global_efficiency(g)
      value <- global_efficiency(g,weights = 1/(E(subg)$weight))
    } else if (eval == 'local_eff'){
      #value <- average_local_efficiency(g)
      value <- average_local_efficiency(g,weights = 1/(E(subg)$weight))
    } else if (eval == 'smallworldeness'){
      value <- smallworldeness(g)
    } else if (eval == 'richcore'){
      value <- rich_core(g)
    } else if (eval =='strength'){
      value <- mean(strength(g))
    } else if (eval =='betweenness'){
      value <- mean(betweenness(g))
    }
    RES[j] <- value
  }
  RES
}


getResultsInSubNetworks <- function(DF,WM_metric,eval,thresh_method,threshold,rsnet){
  M1 <- computeMetric(DF,'V1',WM_metric,eval,thresh_method,threshold,rsnet)
  M2 <- computeMetric(DF,'V2',WM_metric,eval,thresh_method,threshold,rsnet)
  M3 <- computeMetric(DF,'V3',WM_metric,eval,thresh_method,threshold,rsnet)
  
  idsV1 <- get_subject_ids(DF,"V1")
  idsV2 <- get_subject_ids(DF,"V2")
  idsV3 <- get_subject_ids(DF,"V3")
  
  Gdata <- data.frame(
    ids = c(idsV1,idsV2,idsV3),
    Y=c(M1, M2, M3),
    Group =factor(rep(c("V1", "V2", "V3"), times=c(length(M1), length(M2), length(M3))))
  )
  
  Gdata
}
  
# LIMBIC CLUSTERING COEF ????
# R <- getResultsInSubNetworks(DF,"PearsonCorrel","clust_coeff","threshold",0,"Limbic")
# 
# g <- ggviolin(Gdata, x = "Group", y = "Y", fill = "Group",
#               palette = "npg",#color_palette,
#               alpha = .6,
#               width = .5,
#               linetype = 'blank',
#               add= "jitter",color = "Group")
# print(g)
# 
# llm  <- lme4::lmer(Y ~ Group + (1|ids),data = R)#,family = Gamma(link = "log"))
# car::Anova(llm,type=2)
# #   
# library(sjPlot)
# tab_model(llm,show.p=T)
# broom.mixed::tidy(llm)
# 
# shapiro.test(R$Y)
# shapiro.test(R[R$Group=="V1",]$Y)
# shapiro.test(R[R$Group=="V2",]$Y)
# shapiro.test(R[R$Group=="V3",]$Y)

# CONT CHAR PATH
# 

library(ggpubr)
R <- getResultsInSubNetworks(DF,"PearsonCorrel","strength","threshold",0,"Limbic")
# 
g <- ggviolin(R, x = "Group", y = "Y", fill = "Group",
              palette = "npg",#color_palette,
              alpha = .6,
              width = .5,
              linetype = 'blank',
              add= "jitter",color = "Group")
print(g)
# 
llm  <- lme4::lmer(Y ~ Group + (1|ids),data = R)#,family = Gamma(link = "log"))
car::Anova(llm,type=2)
# #   
# library(sjPlot)
tab_model(llm,show.p=T)
# broom.mixed::tidy(llm)
# 
# shapiro.test(R$Y)
# shapiro.test(R[R$Group=="V1",]$Y)
# shapiro.test(R[R$Group=="V2",]$Y)
# shapiro.test(R[R$Group=="V3",]$Y)


# 
# 
g <- makeGraph(DF,"1","V1","PearsonCorrel",0)
# #
M <- as_adjacency_matrix(g,type="both",attr = "weight",sparse = F)
heatmap.2(M,trace = "none",dendrogram = "none",density.info="none",key=F,labRow = TRUE,labCol = TRUE)

V(g)$name[degree(g) > mean(degree(g)) + sd(degree(g))]



# LrsNet <- c("Default","SomMot","SalVentAttn","Vis","Limbic","Cont","DorsAttn")
# rsNet <- LrsNet[7]
# 
# i_DMN <- grep(rsNet,V(g)$name)
# #
# nodes_DMN <- V(g)$name[i_DMN]
# subg <- induced_subgraph(
#   g,
#   i_DMN
# )
# Msub <- as_adjacency_matrix(subg,type="both",attr = "weight",sparse = F)
# heatmap.2(Msub,trace = "none",notecex = 1,dendrogram = "none",density.info="none",key =F,main = rsNet,margins = c(15,15))
# 



# 
# weight <- "PearsonCorrel"
# m <- "global_eff"
# net <- "Default"
# G <- getResultsInSubNetworks(DF,weight,m,"threshold",0,net)


