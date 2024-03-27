library(NetworkConhect)
library(igraph)
library(gplots)
# data_path = "/Volumes/LaCie/derivatives/grouped_results/SC_ROIs_all.csv"
# DF <- read.csv(data_path)
# colnames(DF) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
# DF <- DF[c("subject_id","visit_id","from","to","weight")]

DF$subject_id <- as.character(DF$subject_id)
DF$visit_id <- as.character(DF$visit_id)

weight <- "PearsonCorr"

g <- makeGraph(DF,"43","V3","FBC",0)

M <- as_adjacency_matrix(g,type="both",attr = "weight",sparse = F)
heatmap.2(M,trace = "none",dendrogram = "none",density.info="none")



#corr_mat <- read.csv("/Volumes/LaCie/derivatives/Patients/sub-01/ses-001/functional_connectivity/sub-01_ses-001_aal_correlation_mat.csv",header = F)
#labels <- c('L Aud', 'R Aud', 'Striate', 'L DMN', 'Med DMN', 'Front DMN', 'R DMN', 'Occ post', 'Motor', 'R DLPFC', 'R Front pol', 'R Par', 'R Post Temp', 'Basal', 'L Par', 'L DLPFC', 'L Front pol', 'L IPS', 'R IPS', 'L LOC', 'Vis', 'R LOC', 'D ACC', 'V ACC', 'R A Ins', 'L STS', 'R STS', 'L TPJ', 'Broca', 'Sup Front S', 'R TPJ', 'R Pars Op', 'Cereb', 'Dors PCC', 'L Ins', 'Cing', 'R Ins', 'L Ant IPS', 'R Ant IPS')

# labels_ho <- c('Frontal Pole', 'Insular Cortex', 'Superior Frontal Gyrus', 'Middle Frontal Gyrus', 'Inferior Frontal Gyrus, pars triangularis', 'Inferior Frontal Gyrus, pars opercularis', 'Precentral Gyrus', 'Temporal Pole', 'Superior Temporal Gyrus, anterior division', 'Superior Temporal Gyrus, posterior division', 'Middle Temporal Gyrus, anterior division', 'Middle Temporal Gyrus, posterior division', 'Middle Temporal Gyrus, temporooccipital part', 'Inferior Temporal Gyrus, anterior division', 'Inferior Temporal Gyrus, posterior division', 'Inferior Temporal Gyrus, temporooccipital part', 'Postcentral Gyrus', 'Superior Parietal Lobule', 'Supramarginal Gyrus, anterior division', 'Supramarginal Gyrus, posterior division', 'Angular Gyrus', 'Lateral Occipital Cortex, superior division', 'Lateral Occipital Cortex, inferior division', 'Intracalcarine Cortex', 'Frontal Medial Cortex', 'Juxtapositional Lobule Cortex (formerly Supplementary Motor Cortex)', 'Subcallosal Cortex', 'Paracingulate Gyrus', 'Cingulate Gyrus, anterior division', 'Cingulate Gyrus, posterior division', 'Precuneous Cortex', 'Cuneal Cortex', 'Frontal Orbital Cortex', 'Parahippocampal Gyrus, anterior division', 'Parahippocampal Gyrus, posterior division', 'Lingual Gyrus', 'Temporal Fusiform Cortex, anterior division', 'Temporal Fusiform Cortex, posterior division', 'Temporal Occipital Fusiform Cortex', 'Occipital Fusiform Gyrus', 'Frontal Opercular Cortex', 'Central Opercular Cortex', 'Parietal Opercular Cortex', 'Planum Polare', "Heschl's Gyrus (includes H1 and H2)", 'Planum Temporale', 'Supracalcarine Cortex', 'Occipital Pole')
# labels_aal <- c('Precentral_L', 'Precentral_R', 'Frontal_Sup_L', 'Frontal_Sup_R', 'Frontal_Sup_Orb_L', 'Frontal_Sup_Orb_R', 'Frontal_Mid_L', 'Frontal_Mid_R', 'Frontal_Mid_Orb_L', 'Frontal_Mid_Orb_R', 'Frontal_Inf_Oper_L', 'Frontal_Inf_Oper_R', 'Frontal_Inf_Tri_L', 'Frontal_Inf_Tri_R', 'Frontal_Inf_Orb_L', 'Frontal_Inf_Orb_R', 'Rolandic_Oper_L', 'Rolandic_Oper_R', 'Supp_Motor_Area_L', 'Supp_Motor_Area_R', 'Olfactory_L', 'Olfactory_R', 'Frontal_Sup_Medial_L', 'Frontal_Sup_Medial_R', 'Frontal_Med_Orb_L', 'Frontal_Med_Orb_R', 'Rectus_L', 'Rectus_R', 'Insula_L', 'Insula_R', 'Cingulum_Ant_L', 'Cingulum_Ant_R', 'Cingulum_Mid_L', 'Cingulum_Mid_R', 'Cingulum_Post_L', 'Cingulum_Post_R', 'Hippocampus_L', 'Hippocampus_R', 'ParaHippocampal_L', 'ParaHippocampal_R', 'Amygdala_L', 'Amygdala_R', 'Calcarine_L', 'Calcarine_R', 'Cuneus_L', 'Cuneus_R', 'Lingual_L', 'Lingual_R', 'Occipital_Sup_L', 'Occipital_Sup_R', 'Occipital_Mid_L', 'Occipital_Mid_R', 'Occipital_Inf_L', 'Occipital_Inf_R', 'Fusiform_L', 'Fusiform_R', 'Postcentral_L', 'Postcentral_R', 'Parietal_Sup_L', 'Parietal_Sup_R', 'Parietal_Inf_L', 'Parietal_Inf_R', 'SupraMarginal_L', 'SupraMarginal_R', 'Angular_L', 'Angular_R', 'Precuneus_L', 'Precuneus_R', 'Paracentral_Lobule_L', 'Paracentral_Lobule_R', 'Caudate_L', 'Caudate_R', 'Putamen_L', 'Putamen_R', 'Pallidum_L', 'Pallidum_R', 'Thalamus_L', 'Thalamus_R', 'Heschl_L', 'Heschl_R', 'Temporal_Sup_L', 'Temporal_Sup_R', 'Temporal_Pole_Sup_L', 'Temporal_Pole_Sup_R', 'Temporal_Mid_L', 'Temporal_Mid_R', 'Temporal_Pole_Mid_L', 'Temporal_Pole_Mid_R', 'Temporal_Inf_L', 'Temporal_Inf_R', 'Cerebelum_Crus1_L', 'Cerebelum_Crus1_R', 'Cerebelum_Crus2_L', 'Cerebelum_Crus2_R', 'Cerebelum_3_L', 'Cerebelum_3_R', 'Cerebelum_4_5_L', 'Cerebelum_4_5_R', 'Cerebelum_6_L', 'Cerebelum_6_R', 'Cerebelum_7b_L', 'Cerebelum_7b_R', 'Cerebelum_8_L', 'Cerebelum_8_R', 'Cerebelum_9_L', 'Cerebelum_9_R', 'Cerebelum_10_L', 'Cerebelum_10_R', 'Vermis_1_2', 'Vermis_3', 'Vermis_4_5', 'Vermis_6', 'Vermis_7', 'Vermis_8', 'Vermis_9', 'Vermis_10')
# corr_mat <- as.matrix(corr_mat)
# # 
# heatmap.2(corr_mat,trace = "none",dendrogram = "none",density.info="none")

# 
# pos_mat <- corr_mat
# pos_mat[pos_mat < 0] <- 0
# pos_mat <- as.matrix(pos_mat)
# 
# rownames(pos_mat) <- labels
# colnames(pos_mat) <- labels
# 
# heatmap.2(pos_mat,trace = "none",dendrogram = "none",density.info="none")
# 
# 
# g <- graph_from_adjacency_matrix(pos_mat,weighted = TRUE,mode = "undirected")
# 
# #gT <- sparseThresh(g,0.15)
# 
# #MT <- as_adjacency_matrix(gT,type = "both", attr = "weight", sparse = F)
# #heatmap.2(MT,trace = "none",dendrogram = "none",density.info="none")
# 
# L <- mean_distance(g, weights = E(g)$weight, directed = F)
# print(L)
# C <- transitivity(g,type = "globalundirected",weights = E(g)$weight )
# print(C)
# GE <- global_efficiency(g, weights = E(g)$weight)
# print(GE)
# SW <- C/L
# 



