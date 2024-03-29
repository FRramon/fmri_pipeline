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

### Utils functions  ####

readLUT <- function(lut_path){
  doc <- read.table(lut_path)
  colnames(doc) <- c("No","labelname","r","g","b","A")
  df <- as.data.frame(doc)
  df[,c("No","labelname")]
}

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
# 
# G <- NetworkConhect::getData(DF, weight,"global_eff","density",threshold)
# 
# g <- ggviolin(G, x = "Group", y = "Y", fill = "Group",
#               palette = "npg",#color_palette,
#               alpha = .6,
#               width = .5,
#               linetype = 'blank',
#               add= "jitter",color = "Group")
# g
# 
# res.glmer <- lme4::lmer(Y ~ Group + (1|ids),data = G)
# Anova(res.glmer)
#res.lmer <- lme4::lmer(Y ~ Group + (1|ids),data = G)
#print(tidy(res.glmer,"fixed"),caption = paste("Analyzing ",m, " ~ Group + (1|id)"))
#print(tidy(res.lmer,"fixed"),caption = paste("Analyzing ",m, " ~ Group + (1|id)"))
g_F <- makeGraphFunc(DF_FC,"12","V1",weight,0)
g_S <- makeGraph(DF_SC,"1-E-101-MG","V1",weight_S,0)

cortical_nodes <- V(g_S)$name[as.numeric((V(g_S))$name)>100]
g_cort_S <- induced_subgraph(g_S,cortical_nodes)

## Define a conversion look up table

LUT <- getLUT("/Volumes/LaCie/Sync/code/FreeSurferColorLUT.txt")
LUT <- LUT[LUT$No %in% cortical_nodes,]
functional_nodes <- V(g_F)$name

## Rename structural nodes after functional nodes
LUT$functional_labels <- functional_nodes
V(g_cort_S)
V(g_cort_S)$name <- LUT$functional_labels[match(V(g_cort_S)$name, LUT$No)]


MS <- as_adjacency_matrix(g_cort_S,type="both",attr = "weight",sparse = F)
MF <- as_adjacency_matrix(g_F,type="both",attr = "weight",sparse = F)


cor(MS[1,],MF[1,],method = "spearman",  use = "pairwise.complete.obs")




vertex_list <- V(g)$name
#vertex_list <- vertex_list[!grepl("medial wall", vertex_list) & !grepl("background", vertex_list)]

gT <- sparseThresh(g,threshold)
M <- as_adjacency_matrix(g,type="both",attr = "weight",sparse = F)
heatmap.2(M,trace = "none",dendrogram = "none",key=T,col = viridis(100))

MT <- as_adjacency_matrix(gT,type="both",attr = "weight",sparse = F)
heatmap.2(MT,trace = "none",dendrogram = "none",key=T,col = viridis(100))


# ve <- unname(M[,390])
# matrixa<- matrix(ve,ncol = 400)
# image(matrixa, axes = FALSE, col = viridis(100))
# axis(1, at = 1:length(ve), labels = ve)

# nilearn_nodes <- c('L G_and_S_frontomargin', 'L G_and_S_occipital_inf', 'L G_and_S_paracentral', 'L G_and_S_subcentral', 'L G_and_S_transv_frontopol', 'L G_and_S_cingul-Ant', 'L G_and_S_cingul-Mid-Ant', 'L G_and_S_cingul-Mid-Post', 'L G_cingul-Post-dorsal', 'L G_cingul-Post-ventral', 'L G_cuneus', 'L G_front_inf-Opercular', 'L G_front_inf-Orbital', 'L G_front_inf-Triangul', 'L G_front_middle', 'L G_front_sup', 'L G_Ins_lg_and_S_cent_ins', 'L G_insular_short', 'L G_occipital_middle', 'L G_occipital_sup', 'L G_oc-temp_lat-fusifor', 'L G_oc-temp_med-Lingual', 'L G_oc-temp_med-Parahip', 'L G_orbital', 'L G_pariet_inf-Angular', 'L G_pariet_inf-Supramar', 'L G_parietal_sup', 'L G_postcentral', 'L G_precentral', 'L G_precuneus', 'L G_rectus', 'L G_subcallosal', 'L G_temp_sup-G_T_transv', 'L G_temp_sup-Lateral', 'L G_temp_sup-Plan_polar', 'L G_temp_sup-Plan_tempo', 'L G_temporal_inf', 'L G_temporal_middle', 'L Lat_Fis-ant-Horizont', 'L Lat_Fis-ant-Vertical', 'L Lat_Fis-post', 'L Pole_occipital', 'L Pole_temporal', 'L S_calcarine', 'L S_central', 'L S_cingul-Marginalis', 'L S_circular_insula_ant', 'L S_circular_insula_inf', 'L S_circular_insula_sup', 'L S_collat_transv_ant', 'L S_collat_transv_post', 'L S_front_inf', 'L S_front_middle', 'L S_front_sup', 'L S_interm_prim-Jensen', 'L S_intrapariet_and_P_trans', 'L S_oc_middle_and_Lunatus', 'L S_oc_sup_and_transversal', 'L S_occipital_ant', 'L S_oc-temp_lat', 'L S_oc-temp_med_and_Lingual', 'L S_orbital_lateral', 'L S_orbital_med-olfact', 'L S_orbital-H_Shaped', 'L S_parieto_occipital', 'L S_pericallosal', 'L S_postcentral', 'L S_precentral-inf-part', 'L S_precentral-sup-part', 'L S_suborbital', 'L S_subparietal', 'L S_temporal_inf', 'L S_temporal_sup', 'L S_temporal_transverse', 'R G_and_S_frontomargin', 'R G_and_S_occipital_inf', 'R G_and_S_paracentral', 'R G_and_S_subcentral', 'R G_and_S_transv_frontopol', 'R G_and_S_cingul-Ant', 'R G_and_S_cingul-Mid-Ant', 'R G_and_S_cingul-Mid-Post', 'R G_cingul-Post-dorsal', 'R G_cingul-Post-ventral', 'R G_cuneus', 'R G_front_inf-Opercular', 'R G_front_inf-Orbital', 'R G_front_inf-Triangul', 'R G_front_middle', 'R G_front_sup', 'R G_Ins_lg_and_S_cent_ins', 'R G_insular_short', 'R G_occipital_middle', 'R G_occipital_sup', 'R G_oc-temp_lat-fusifor', 'R G_oc-temp_med-Lingual', 'R G_oc-temp_med-Parahip', 'R G_orbital', 'R G_pariet_inf-Angular', 'R G_pariet_inf-Supramar', 'R G_parietal_sup', 'R G_postcentral', 'R G_precentral', 'R G_precuneus', 'R G_rectus', 'R G_subcallosal', 'R G_temp_sup-G_T_transv', 'R G_temp_sup-Lateral', 'R G_temp_sup-Plan_polar', 'R G_temp_sup-Plan_tempo', 'R G_temporal_inf', 'R G_temporal_middle', 'R Lat_Fis-ant-Horizont', 'R Lat_Fis-ant-Vertical', 'R Lat_Fis-post', 'R Pole_occipital', 'R Pole_temporal', 'R S_calcarine', 'R S_central', 'R S_cingul-Marginalis', 'R S_circular_insula_ant', 'R S_circular_insula_inf', 'R S_circular_insula_sup', 'R S_collat_transv_ant', 'R S_collat_transv_post', 'R S_front_inf', 'R S_front_middle', 'R S_front_sup', 'R S_interm_prim-Jensen', 'R S_intrapariet_and_P_trans', 'R S_oc_middle_and_Lunatus', 'R S_oc_sup_and_transversal', 'R S_occipital_ant', 'R S_oc-temp_lat', 'R S_oc-temp_med_and_Lingual', 'R S_orbital_lateral', 'R S_orbital_med-olfact', 'R S_orbital-H_Shaped', 'R S_parieto_occipital', 'R S_pericallosal', 'R S_postcentral', 'R S_precentral-inf-part', 'R S_precentral-sup-part', 'R S_suborbital', 'R S_subparietal', 'R S_temporal_inf', 'R S_temporal_sup', 'R S_temporal_transverse')
# out_nodes <- unique(DF$from)
# setdiff(nilearn_nodes,out_nodes)
# 


M <- as_adjacency_matrix(g,type="both",attr = "weight",sparse = F)
heatmap.2(M,trace = "none",dendrogram = "none",key=T,col = viridis(100))

# ### Nodes fs
# not_lateralized_LUT <- function(){
#   Node_list <- as.numeric(V(g)$name)
#   LUT <- readLUT("Users/francoisramon/Desktop/These/dMRI_pipeline/code/FreeSurferColorLUT.txt")
#   LUT <- LUT[LUT$No %in% sublist,]
#   #LUT$labelname <- substr(LUT$labelname,start=8,stop = nchar(LUT$labelname))
#   LUT
# }
# 
# LUT <- not_lateralized_LUT()
# length(grep("ctx_rh",LUT$labelname))
# ### nodes nilearn
# 
# nodes_nilearn <- data.frame(labelname = V(g)$name)
# #nodes_nilearn$labelname <- substr(nodes_nilearn$labelname,start = 3,stop =  nchar(nodes_nilearn$labelname))
# 
# left <- nodes_nilearn$labelname[grep("L ",nodes_nilearn$labelname)]
# right <-  nodes_nilearn$labelname[grep("R ",nodes_nilearn$labelname)]
# 
# sub_left <- substr(left,start = 3,stop =  nchar(left))
# sub_right <- substr(right,start = 3,stop =  nchar(right))                   
# 
# 
# 
# I = intersect(nodes_nilearn$labelname,LUT$labelname)
# S = setdiff(nodes_nilearn$labelname,LUT$labelname)

M <- read.csv("Volumes/LaCie/derivatives/Patients/sub-01/ses-002/functional_connectivity/sub-01_ses-002_destrieux_correlation_mat.csv",header=F)
M <- as.matrix(M,sparse=F)
M[M<0]<- 0
# gM <- graph_from_adjacency_matrix(M)