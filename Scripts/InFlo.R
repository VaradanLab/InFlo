###############################################INFLO-FINAL##################################################
#InFlo is a novel systems biology approach for characterizing the activities of complex signaling networks
#using a unique multidimensional framework integrating transcriptomic, genomic and/or epigenomic profiles 
#for any given biological sample.
############################################################################################################
#####Packages_required 
rm(list=ls())
source("/Projects/VaradanLab_Repositories/InFlo/Run_Conf_FIles/inFlo_RUN_CONF.txt")
source("/Projects/VaradanLab_Repositories/InFlo/Scripts/Engine.R")


Initial_chk <- File_chk()

Dir_create(Initial_chk)


run <- run_chk()
#############################################################################################################
run <- run_chk()

##############################################~~~~Data_Import~############################################## 
GE_Data <- Data_Read(GE_FILE)
CNV_Data <- Data_Read(CNV_FILE)
Samp_Info <-  read.table(SAMPLE_INFORMATION,header=T,check.names = F,stringsAsFactors = F)

tumor_samples <- Samp_Info[which(Samp_Info[,'Sample_Type']=="Tumor"),"Sample_Name"]
normal_samples <- Samp_Info[which(Samp_Info[,'Sample_Type']=="Normal"),"Sample_Name"]
############################################################################################ 

if(RNASeqV2){
  GE_WILCOX <<- DeSEQ_TEST(GE_Data)
  CNV_WILCOX <<- Wilcox_Test(CNV_Data)
}else{
  GE_WILCOX <<- Wilcox_Test(GE_Data)
  CNV_WILCOX <<- Wilcox_Test(CNV_Data)
}

if(GUASS){
  GE_WILCOX <- Guass_Fit(GE_Data)
  CNV_WILCOX <- Guass_Fit(CNV_Data)
}

PATHWAYS <- paste(anaTemp,"/pathways",sep="")

PRE_INFLO(GE_WILCOX,CNV_WILCOX,PATHWAYS)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   