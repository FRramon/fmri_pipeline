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
library(sjPlot)
library(tidyverse)
library(DescTools)
library(car)
library(gplots)
library(ggpubr)
library(lme4)
library(broom.mixed)
library(gridExtra)
library(readxl)
library(viridis)
library(ggseg)
library(ggsegSchaefer)
library(ggseg3d)
path <- "/Volumes/LaCie/derivatives/grouped_results"
data_path = "/Volumes/LaCie/derivatives/grouped_results/FC_schaefer.csv"
DF <- read.csv(data_path)
colnames(DF) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
DF <- DF[c("subject_id","visit_id","from","to","weight")]


DF$subject_id <- as.character(DF$subject_id)
DF$visit_id <- as.character(DF$visit_id)

DF$visit_id[DF$visit_id == "1"] <- "V1"
DF$visit_id[DF$visit_id == "2"] <- "V2"
DF$visit_id[DF$visit_id == "3"] <- "V3"

weight <- "PearsonCorrel"
real_weight <- "PearsonCorrel"

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

get_rich_club <- function(g){
  
  plot(degree(g))
  
  sorted_degrees <- sort(degree(g))
  
  top <- as.integer(0.2*length(sorted_degrees))
  print(top)
  
  ## Rich nodes
  richnodes <- tail(sorted_degrees,top)
  print(richnodes)
  # Local nodes
  nonrichnodes <- head(sorted_degrees,length(sorted_degrees)-top)
  
  ## Rich rich connections
  richsubg <- subgraph(g,richnodes)
  
  ## Local local connection
  localsubg <- subgraph(g,nonrichnodes)
  
  ## feeder
  gworich <- difference(g,richsubg)
  gfeeder <- difference(gworich,localsubg)
  
  ## All subnetworks
  #richfeeder
  richfeeder <- difference(g,localsubg)
  
  # feederlocal
  
  feederlocal <- difference(g,richsubg)
  
  list("rich" = richsubg,"local" = localsubg,"feeder" = gfeeder, "richfeeder" = richfeeder,"feederlocal" =  feederlocal)
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
      value <- mean_distance(g,weights = 1/(E(g)$weight))
      #value <- normalized_shortest_path(g)
    } else if (eval == 'global_eff'){
      # value <- global_efficiency(g)
      value <- global_efficiency(g,weights = 1/(E(g)$weight))
    } else if (eval == 'local_eff'){
      #value <- average_local_efficiency(g)
      value <- average_local_efficiency(g,weights = 1/(E(g)$weight))
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

# data_path = "/Volumes/LaCie/derivatives/grouped_results/FC_schaefer.csv"
# DF <- read.csv(data_path)


compute_Y_on_thresholds <- function(DF,weight,m,set_threshold){
  
  list_of_G <- list()
  mean_Y <- c()
  i <- 1
  for (t in set_threshold){
    print(paste0("Computing on threshold : ",t))
    
    G <- getResultsInSubNetworks(DF,weight,m,"density",t,"Limbic")
    list_of_G[[i]] <- G
    i <- i +1
    mean_Y <- c(mean_Y,mean(G$Y))
  }
  # for (i in 1:length(set_threshold)){
  # }
  
  multi_thresh_DF <- data.frame(thresholds = set_threshold, mean_y = mean_Y)
  
  
  
  result_list <- vector("list", nrow(list_of_G[[1]]))
  
  for (i in 1:nrow(list_of_G[[1]])) {
    row_info <- vector("list", length(list_of_G))
    for (j in 1:length(list_of_G)) {
      subject_id <- list_of_G[[j]]$ids[i]
      group <- list_of_G[[j]]$Group[i]
      y_value <- list_of_G[[j]]$Y[i]
      row_info[[j]] <- list(subject_id = subject_id, group = group, Y = y_value)
    }
    result_list[[i]] <- row_info
  }
  
  df <- data.frame(ids = character(), Group = character(),Y = numeric())
  df_only_values <- data.frame(X = set_threshold)
  for(i in 1:length(result_list)){
    Y_over_t <- c()
    for (j in 1:length(result_list[[i]])){
      Y_over_t <- c(Y_over_t,result_list[[i]][[j]]$Y)
    }
    
    df_only_values <- cbind(df_only_values,Y_over_t)
    
    Y_AUC <- AUC(set_threshold,Y_over_t)
    group <- result_list[[i]][[1]]$group
    subject_id <- result_list[[i]][[1]]$subject_id
    
    df <- df %>% 
      add_row(ids = subject_id, Group = group, Y = Y_AUC)
    
  }
  library(tidyverse)
  names(df_only_values)[2:ncol(df_only_values)] <- paste0("Y_over_t_", 1:(ncol(df_only_values)-1))
  df_long <- df_only_values %>%
    gather(key = "Y_variable", value = "Y_value", -X)
  
  
  
  df_means <- df_only_values %>%
    # Select the columns for each range and calculate the mean
    mutate(mean_V1 = rowMeans(select(., Y_over_t_1:Y_over_t_32)),
           mean_V2 = rowMeans(select(., Y_over_t_33:Y_over_t_63)),
           mean_V3 = rowMeans(select(., Y_over_t_64:Y_over_t_89)))  %>%
    # Select only the necessary columns
    select(X, mean_V1, mean_V2, mean_V3)

  df_means_long <- df_means %>%
    gather(key = "mean_variable", value = "mean_value", -X)
  
  # Plot the mean lines against X
  
  
  for (i in 1:length(list_of_G)){
    res.glmer <- lme4::lmer(Y ~ Group + (1|ids),data = list_of_G[[i]])
    print(Anova(res.glmer)$`Pr(>Chisq)`)
  }
  
  list("listG" = list_of_G,"dfallgroups" = multi_thresh_DF,"dfauc" = df,"dfonly" = df_only_values,"dfallplots" = df_long,"dfgroupplot" = df_means_long)
}

threshold_a <- 0.1
threshold_b <- 0.3

set_threshold = seq(threshold_a,threshold_b,0.1)

R <- compute_Y_on_thresholds(DF,"PearsonCorrel","global_eff",set_threshold)
g <- ggviolin(R$dfauc, x = "Group", y = "Y", fill = "Group",
              palette = "npg",#color_palette,
              alpha = .6,
              width = .5,
              linetype = 'blank',
              add= "jitter",color = "Group")
print(g)
# 
llm  <- lme4::lmer(Y ~ Group + (1|ids),data = R$dfauc)#,family = Gamma(link = "log"))
car::Anova(llm,type=2)
# #   
# library(sjPlot)
tab_model(llm,show.p=T)

######### CLINICAL CORRELATION #######

clinique <- read_excel("/Volumes/LaCie/Sync/Graph_project/results/Clinique_pourfrancois.xlsx")
clinique <- clinique[,c('subject_id',"MADRS_V1","HDRS_21_V1","remission","reponse")]
equivalence_table <- read.csv("/Users/francoisramon/Desktop/These/CONHECT/dMRI_pipeline/equivalence_table_Patients.csv")

dfauc <- R$dfauc
colnames(dfauc) <- c("functional_label","Group","Y")
DFV1 <- dfauc[dfauc$Group=="V1",]

df_with_demo <- merge(dfauc,equivalence_table,by = "functional_label")
df_with_demo_clinic <- merge(df_with_demo,clinique,by ="subject_id")

df_with_demo_clinic<- drop_na(df_with_demo_clinic)
DFV1 <- df_with_demo_clinic[df_with_demo_clinic$Group=="V1",]
DFV2 <- df_with_demo_clinic[df_with_demo_clinic$Group=="V2",]
DFV3 <- df_with_demo_clinic[df_with_demo_clinic$Group=="V3",]

ggbarplot(DFV1, x = "remission", y = "Y", add = c("mean_se", "jitter"),fill = "remission", color = "black",
          palette = c("#E7B800", "#FC4E07"))



res.lm <- lm(Y~MADRS_V1,data=DFV1)
Anova(res.lm)


ggplot(DFV1, aes(x = MADRS_V1, y = Y)) +
  geom_point() +
  stat_smooth(method = "lm", col = "skyblue")





# #
g <- makeGraph(DF,"1","V1","PearsonCorrel",0)


# # #
# M <- as_adjacency_matrix(g,type="both",attr = "weight",sparse = F)
# heatmap.2(M,trace = "none",dendrogram = "none",density.info="none",key=F,labRow = TRUE,labCol = TRUE)
#
#V(g)$name[degree(g) > mean(degree(g)) + sd(degree(g))]

LrsNet <- c("Default","SomMot","SalVentAttn","Vis","Limbic","Cont","DorsAttn")
rsNet <- LrsNet[1]
i_DMN <- grep(rsNet,V(g)$name)
# #
nodes_DMN <- V(g)$name[i_DMN]
subg <- induced_subgraph(
   g,
   i_DMN )


Msub <- as_adjacency_matrix(subg,type="both",attr = "weight",sparse = F)
heatmap.2(Msub,trace = "none",notecex = 1,dendrogram = "none",density.info="none",col = viridis(100),key =F,main = rsNet,margins = c(15,15))


G <- get_rich_club(g)

labels_rich <- V(G$rich)$name
# destr<- ggseg (atlas = ggsegSchaefer::schaefer7_400,
#                mapping = aes(fill = labels_rich)) +
#   theme(legend.justification = c(1, 0),
#         legend.position = "bottom" ,
#         legend.text = element_text(size = 5)) +
#   guides(fill = guide_legend(ncol = 3))
# 
subg <-g
dfdegree = data.frame(region = V(subg)$name,ispresent = rep(1,length(V(subg)$name)),degree = unname(degree(subg,V(subg)$name)),strength = unname(strength(subg,V(subg)$name)),betw = betweenness(subg,V(subg)$name))

plot(hoCort) +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 6)) +
  guides(fill = guide_legend(ncol = 3))

pl <- dfdegree %>% ggseg(atlas = ggsegSchaefer::schaefer7_400,position = "stacked",
  mapping = aes(fill = degree)) + 
  scale_fill_gradient2(low = "blue",mid = "white", high = "red", limits = range(dfdegree$degree)) 
pl

pl3d <- ggseg3d(.data = dfdegree,
                 atlas = schaefer7_400_3d,
                 surface = 'inflated',

                 colour = "degree", text = "degree",palette = heat.colors(100),na.colour = "lightgrey") %>%
  pan_camera("right lateral") %>%
  remove_axes()
pl3d
# 
# weight <- "PearsonCorrel"
# m <- "global_eff"
# net <- "Default"
# G <- getResultsInSubNetworks(DF,weight,m,"threshold",0,net)


