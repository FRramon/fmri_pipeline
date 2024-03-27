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
path <- "/Volumes/LaCie/derivatives/grouped_results"
data_path = "/Volumes/LaCie/derivatives/grouped_results/FC_msdl.csv"
DF <- read.csv(data_path)
colnames(DF) <- c('subject_id','visit_id','from_i','to_j','weight',"from","to")
DF <- DF[c("subject_id","visit_id","from","to","weight")]


DF$subject_id <- as.character(DF$subject_id)
DF$visit_id <- as.character(DF$visit_id)

DF$visit_id[DF$visit_id == "1"] <- "V1"
DF$visit_id[DF$visit_id == "2"] <- "V2"
DF$visit_id[DF$visit_id == "3"] <- "V3"

weight <- "FBC"
real_weight <- "PearsonCorrel"
threshold <- 0.1
G <- NetworkConhect::getData(DF, weight,"characteristic_path","density",threshold)

g <- ggviolin(G, x = "Group", y = "Y", fill = "Group",
              palette = "npg",#color_palette,
              alpha = .6,
              width = .5,
              linetype = 'blank',
              add= "jitter",color = "Group")
g

res.glmer <- lme4::glmer(Y ~ Group + (1|ids),data = G,family = Gamma(link = "log"))
#res.lmer <- lme4::lmer(Y ~ Group + (1|ids),data = G)
print(tidy(res.glmer,"fixed"),caption = paste("Analyzing ",m, " ~ Group + (1|id)"))
#print(tidy(res.lmer,"fixed"),caption = paste("Analyzing ",m, " ~ Group + (1|id)"))

g <- makeGraph(DF,"1","V3","FBC",0)
#gT <- sparseThresh(g,threshold)
M <- as_adjacency_matrix(g,type="both",attr = "weight",sparse = F)
heatmap.2(M,trace = "none",dendrogram = "none",density.info="none")
heatmap(M,col = terrain.colors(100))





