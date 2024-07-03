library(tidyr)
library(readxl)
library(dplyr)
library(readODS)


# Util function 
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Function to join clinical data
join_clinical_data <- function(clinique, D){
  D <- rename(D,subject_id = ids, session_id = Group)
  D <- inner_join(clinique,D,by = "subject_id")
  D
}

side_to_side_Y <- function(res_dataframe){

  df_V1 <- res_dataframe[res_dataframe$Group == "V1",]
  df_V2 <- res_dataframe[res_dataframe$Group== "V2",]
  df_V3 <- res_dataframe[res_dataframe$Group == "V3",]
 

  df_V12 <- left_join(df_V1,df_V2, by = "subject_id")
  df_V123 <- left_join(df_V12,df_V3, by = "subject_id")

  df_V123 <- rename(df_V123,Y_V1 = Y.x, Y_V2 = Y.y, Y_V3 = Y,i = X.x,subject_id = subject_id)

  #write.csv(df_V123,"/home/francoisramon/Desktop/side_to_side_Y.csv")
  df_V123[,c("i","subject_id","Y_V1","Y_V2","Y_V3")]

}

load_and_clean <- function(wd,filename){
  ## Choose the ODS file or the xlsx 
  ## Clean by taking only patients that went up to V3
  ## Remove oulier patient 147

  if (substrRight(filename,4) ==".ods"){
    clinical_data = read_ods(paste(wd,'/Graph_project/results/',filename,sep=""))
  }else{
    clinical_data = read_xlsx(paste(wd,'/Graph_project/results/',filename,sep=""))
  }
  clean_clinical_data = clinical_data[complete.cases(clinical_data$remission),]
  clean_clinical_data = clean_clinical_data[clean_clinical_data$subject_id != '1-E-147-FE',]
  
  clean_clinical_data
}

get_full_dataframe <- function(wd,metric,atlas){



  clinical_data <-  load_and_clean(wd,'full_clinical_data_conhect.ods')
  path_to_results <- paste0(wd,'/fMRI_pipeline/results/',atlas,'/auc_table_',metric,"_",atlas,'.csv',sep="")
  dresults <- read.csv(path_to_results)
 
  G <- read_excel(paste(wd,'/Graph_project/data/','Genres_conhect.xlsx',sep=""))
  A <- read_excel(paste(wd,'/Graph_project/data/','ages_conhect.xlsx',sep=""))
  
  dwA <- inner_join(clinical_data,A,by = "subject_id")
  dwAG <- inner_join(dwA,G,by = "subject_id")

  dwAG$MADRS_13 <- as.numeric(dwAG$MADRS_V1) - as.numeric(dwAG$MADRS_V3)
  dwAG$MADRS_12 <- as.numeric(dwAG$MADRS_V1) - as.numeric(dwAG$MADRS_V2)
  dwAG$MADRS_23 <- as.numeric(dwAG$MADRS_V2) - as.numeric(dwAG$MADRS_V3)
  

  dwAG$perc_MADRS_13 <- as.numeric(dwAG$MADRS_13) / as.numeric(dwAG$MADRS_V1)
  dwAG$perc_MADRS_23 <- as.numeric(dwAG$MADRS_23) / as.numeric(dwAG$MADRS_V2)


  
  #### Add Scores difference
  dD <- side_to_side_Y(dresults)
  dwAGD <- inner_join(dwAG,dD,by="subject_id")
  dwAGD$Y_V1 <- as.numeric(dwAGD$Y_V1)
  dwAGD$Y_V2 <- as.numeric(dwAGD$Y_V2)
  dwAGD$Y_V3 <- as.numeric(dwAGD$Y_V3)

  dwAGD$Y_13 <-  dwAGD$Y_V3 - dwAGD$Y_V1
  dwAGD$Y_12 <-  dwAGD$Y_V2 - dwAGD$Y_V1
  dwAGD$Y_23 <-  dwAGD$Y_V3 - dwAGD$Y_V2

  dwAGD$perc_Y_13 <- dwAGD$Y_13 / dwAGD$Y_V1
  dwAGD$perc_Y_12 <- dwAGD$Y_12 / dwAGD$Y_V1
  dwAGD$perc_Y_23 <- dwAGD$Y_23 / dwAGD$Y_V2


  Dclean <- filter_dataframe(dwAGD)

  Dclean <- Dclean[!is.na(Dclean$Y_13),]
  #write.csv(Dclean,paste("~/Desktop/These/Graph_project/full_tables/full_",metric,"_",weight,".csv",sep=""))


  Dclean
}
  
filter_dataframe <- function(full_dataframe){
#   "subject_id"         "Date_V1"            "HDRS_21_V1"        
#   "IDS_C_V1"           "MADRS_V1"           "Date_V2"           
#   "IRM2_N_ECT"         "HDRS_21_V2"         "perc_HDRS_V12"     
#   "remission_V2"       "reponse_V2"         "IDS_C_V2"          
#   "MADRS_V2"           "perc_MADRS_V12"     "remission_MADRS_V2"
#   "reponse_MADRS_V2"   "Date_V3"            "IRM3_N_ECT"        
#   "HDRS_21_V3"         "IDS_C_V3"           "MADRS_V3"          
#   "remission"          "reponse"            "X"                 
#   "Y"                  "session_id"         "age"               
#   "sex"                "MADRS_13"           "i"                 
#   "Y_V1"               "Y_V2"               "Y_V3"              
#   "Y_13"               "Y_12"               "Y_23"    

  D <- full_dataframe[,c("subject_id","age","sex","Date_V1","HDRS_21_V1","MADRS_V1","Y_V1","Date_V2","IRM2_N_ECT","HDRS_21_V2","perc_HDRS_V12","remission_V2","reponse_V2","MADRS_V2","MADRS_12","perc_MADRS_V12","remission_MADRS_V2","reponse_MADRS_V2","Y_V2","Y_12","perc_Y_12","Date_V3","IRM3_N_ECT","HDRS_21_V3","MADRS_V3","remission","reponse","MADRS_13","MADRS_23","perc_MADRS_13","perc_MADRS_23","Y_V3","Y_13","perc_Y_13","Y_23","perc_Y_23")]
  D
}



#D <- get_full_dataframe("/mnt/CONHECT_data","global_eff","FBC","V1")







