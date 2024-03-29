################################
# Author : Francois Ramon
# Date of creation : 01/15/2024

# Graph theory analysis of fMRI derived networks
################################

library(NetworkConhect)
library(igraph)
library(gplots)
library(ggpubr)
library(lme4)
library(broom.mixed)

######### 


#### Load functional dataset

path <- "/Volumes/LaCie/derivatives/grouped_results"
data_path = "/Volumes/LaCie/derivatives/grouped_results/FC_destrieux.csv"
DF <- read.csv(data_path)
colnames(DF) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
DF <- DF[c("subject_id","visit_id","from","to","weight")]

DF$subject_id <- as.character(DF$subject_id)
DF$visit_id <- as.character(DF$visit_id)

DF$visit_id[DF$visit_id == "1"] <- "V1"
DF$visit_id[DF$visit_id == "2"] <- "V2"
DF$visit_id[DF$visit_id == "3"] <- "V3"
DF_FC <- DF
### Load Structural dataset

structural_path <- "/Volumes/LaCie/Sync/Graph_project/data/stats_diffusion_metric_in_fullWM_FBC.xlsx"
DF_SC <- read_and_normalize_data(structural_path,"FBC")
#### Define model parameters

weight <- "PearsonCorrel"
weight_S <- "FBC"

threshold <- 0.20

### Run global analysis


equivalence_names <- read.csv("/Volumes/LaCie/equivalence_table_Patients.csv")
# equivalence_names$functional_label <- as.character(seq(1,44,1))
# write.csv(equivalence_names,"/Volumes/LaCie/equivalence_table_Patients.csv")

subject_list <- get_subject_ids(DF_FC,"V1")
subject_list <- subject_list[!grepl("12",subject_list) & !grepl("11",subject_list) ]
LUT <- getLUT("/Volumes/LaCie/Sync/code/FreeSurferColorLUT.txt")
spearman_matrix <- LUT$functional_labels
mean_r <- c()

for(sub in subject_list){
  

  print(sub)
  print(equivalence_names$conhect_label[match(sub,equivalence_names$functional_label)])
  g_F <- makeGraphFunc(DF_FC,sub,"V1",weight,0)

  g_S <- makeGraph(DF_SC,equivalence_names$conhect_label[match(sub,equivalence_names$functional_label)],"V1",weight_S,0)
  
  ## Truncate to keep only cortical nodes
  
  cortical_nodes <- V(g_S)$name[as.numeric((V(g_S))$name)>100]
  g_cort_S <- induced_subgraph(g_S,cortical_nodes)
  # 
  # ## Define a conversion look up table
  # 
  LUT <- LUT[LUT$No %in% cortical_nodes,]
  functional_nodes <- V(g_F)$name
  
  #print(functional_nodes)
# 
#   ## Rename structural nodes after functional nodes
  LUT$functional_labels <- functional_nodes
  #V(g_cort_S)
  V(g_cort_S)$name <- LUT$functional_labels[match(V(g_cort_S)$name, LUT$No)]

  ### Define matrixes


  MS <- as_adjacency_matrix(g_cort_S,type="both",attr = "weight",sparse = F)
  MF <- as_adjacency_matrix(g_F,type="both",attr = "weight",sparse = F)

  spearman_vector <- c()

  for(label in LUT$functional_labels){
    functional_vector <- MF[label,]
    structural_vector <- MS[label,][LUT$functional_labels]

    functional_vector[functional_vector == 0] <- NA
    structural_vector[structural_vector == 0] <- NA

    spearman_corr <- cor(functional_vector, structural_vector, method = "spearman",use = "pairwise.complete.obs")
    spearman_vector <- c(spearman_vector,spearman_corr)
    

  }
  
  spearman_matrix <- cbind(spearman_matrix,spearman_vector)
  mean_r <- c(mean_r,mean(spearman_vector))

  colors <- ifelse(spearman_vector <= 0, "skyblue", "lightcoral")
  barplot(spearman_vector,names.arg = LUT$functional_labels,las=2,col = colors)

}

coupling <- data.frame(subject_id = subject_list, SFC = mean_r)
colnames(equivalence_names) <- c("X","conhect_label","sub_number","subject_id")
coupling <- merge(coupling,equivalence_names, by = "subject_id")
#heatmap.2(spearman_matrix,trace="none",dendrogram = "none",key=T,col = viridis(100))

clinique <- read_excel("/Volumes/LaCie/Sync/Graph_project/results/Clinique_pourfrancois.xlsx")
clinique <- clinique[,c('subject_id',"MADRS_V1","HDRS_21_V1")]
colnames(clinique) <- c("conhect_label","MADRS_V1","HDRS_21_V1")
coupling_clinical <- merge(coupling,clinique,by = "conhect_label")


plot(coupling_clinical$MADRS_V1,coupling_clinical$SFC)
res.lm <- lm(SFC ~ HDRS_21_V1, data = coupling_clinical)

ggplot(coupling_clinical, aes(x = MADRS_V1, y = SFC)) + 
  geom_point() +
  stat_smooth(method = "lm", col = "blue")



## reorganiser selon les labels.
## 0 <- NA


