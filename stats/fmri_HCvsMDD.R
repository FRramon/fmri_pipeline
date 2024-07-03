################################
# Author : Francois Ramon
# Date of creation : 01/15/2024

# Graph theory analysis of fMRI derived networks
# MDD vs HC
################################

library(NetworkConhect)
library(igraph)
library(gplots)
library(ggpubr)
library(lme4)
library(broom.mixed)
library(gridExtra)
library(tidyverse)
library(DescTools)
library(ggsignif)
######### 




#### Load functional dataset MDD

path <- "/Volumes/LaCie/derivatives/grouped_results"
data_path = "/Volumes/LaCie/derivatives/grouped_results/FC_destrieux.csv"
output_dir = "/Users/francoisramon/Desktop/These/CONHECT/fMRI_pipeline/results"
DF <- read.csv(data_path)
colnames(DF) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
DF <- DF[c("subject_id","visit_id","from","to","weight")]

DF$subject_id <- as.character(DF$subject_id)
DF$visit_id <- as.character(DF$visit_id)

DF$visit_id[DF$visit_id == "1"] <- "V1"
DF$visit_id[DF$visit_id == "2"] <- "V2"
DF$visit_id[DF$visit_id == "3"] <- "V3"
DF_FC <- DF

### Define model preference 

weight <- "PearsonCorrel"

#### Load functional dataset HC 

data_path_HC <- paste(path,"/HealthyVolunteers/FC_destrieux.csv",sep="")
output_dir = "/Users/francoisramon/Desktop/These/CONHECT/fMRI_pipeline/results"
DF_HC <- read.csv(data_path_HC)
colnames(DF_HC) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
DF_HC <- DF_HC[c("subject_id","visit_id","from","to","weight")]

DF_HC$subject_id <- as.character(DF_HC$subject_id)
DF_HC$visit_id <- as.character(DF_HC$visit_id)

DF_HC$visit_id[DF_HC$visit_id == "1"] <- "HC"
DF_HC$subject_id <- paste("HC_",DF_HC$subject_id,sep="")

### Define model preference 

weight <- "PearsonCorrel"

DF <- rbind(DF,DF_HC)

G <- NetworkConhect::getData(DF,"PearsonCorrel","global_eff","density",0.2)
ggplot(G, aes(Group, Y)) +
  geom_jitter( position = position_jitter(0.2),
               color = "blue") 


idsHC <- get_subject_ids(DF,"HC")
groupHC <- rep("HC",length(idsHC))
YHC <- computeSW(DF,"HC","PearsonCorrel","global_eff","density",0.2)
GHC <- data.frame(ids = idsHC,Y = YHC, Group = groupHC)

G <- rbind(G,GHC)

shapiro.test(G[G$Group == "V1",]$Y)
shapiro.test(G[G$Group == "HC",]$Y)


G_HCV1 <- G[G$Group=="V1" | G$Group == "HC",]
res.lm <- lm(Y~Group,data = G_HCV1)

ggplot(G_HCV1, aes(x = Group, y = Y, colour = Group)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  geom_signif(comparisons = list(c("V1","HC")),map_signif_level = TRUE)






