###############################################INFLO-FINAL##################################################
#InFlo is a novel systems biology approach for characterizing the activities of complex signaling networks
#using a unique multidimensional framework integrating transcriptomic, genomic and/or epigenomic profiles 
#for any given biological sample.
############################################################################################################
#####Packages_required 
rm(list=ls())
require(doParallel)
args = commandArgs(trailingOnly=TRUE)
registerDoParallel(cores=(Sys.getenv("Num_Of_Cores")))
parallel_time <-system.time({
source(args[2])
source(paste(InFlo_Home,"/Scripts/Engine.R",sep=""))
Run <- run_chk()
if(Run){
Initial_chk <- File_chk()
Dir_create(Initial_chk)
RUN_PIDS <- PATH_PROCESS(PATHWAY_INFORMATION,PATHWAYS_DIR)
##############################################~~~~Data_Import~############################################## 
GE_Data <- Data_Read(GE_FILE)

CNV_Data <- try(Data_Read(CNV_FILE))
if(class(CNV_Data)=="try-error"){
CNV_Data <- NULL
}
Samp_Info <-  read.table(SAMPLE_INFORMATION,header=T,check.names = F,stringsAsFactors = F,sep="\t")
tumor_samples <- Samp_Info[which(Samp_Info[,'Sample_Type']=="Tumor"),"Sample_Name"]
normal_samples <- Samp_Info[which(Samp_Info[,'Sample_Type']=="Normal"),"Sample_Name"]
############################################################################################ 
if(RNASeqV2){
GE_WILCOX <<- DeSEQ_TEST(GE_Data)
if(is.null(CNV_Data)){
  CNV_WILCOX <<- matrix(1,nrow = length(rownames(GE_WILCOX)),ncol = length(colnames(GE_WILCOX)))
  rownames(CNV_WILCOX) <- rownames(GE_WILCOX)
  colnames(CNV_WILCOX) <- colnames(GE_WILCOX)
}else{
  CNV_WILCOX <<- Wilcox_Test(CNV_Data)
}
}else{
GE_WILCOX <<- Wilcox_Test(GE_Data)
if(is.null(CNV_Data)){
  CNV_WILCOX <<- matrix(1,nrow = length(rownames(GE_WILCOX)),ncol = length(colnames(GE_WILCOX)))
  rownames(CNV_WILCOX) <- rownames(GE_WILCOX)
  colnames(CNV_WILCOX) <- colnames(GE_WILCOX)
}else{
  CNV_WILCOX <<- Wilcox_Test(CNV_Data)
}
}

# if(Providing_Pvals){
#   GE_WILCOX <- GE_Data
#   CNV_WILCOX <<- matrix(1,nrow = length(rownames(GE_WILCOX)),ncol = length(colnames(GE_WILCOX)))
#   rownames(CNV_WILCOX) <- rownames(GE_WILCOX)
#   colnames(CNV_WILCOX) <- colnames(GE_WILCOX)
# }



# if(GUASS){
#   GE_WILCOX <- Guass_Fit(GE_Data)
#   CNV_WILCOX <- Guass_Fit(CNV_Data)
# }
PATHWAYS <- paste(anaPath,"/pathways/",sep="")
PRE_INFLO(GE_WILCOX,CNV_WILCOX,PATHWAYS)
InFlo(PATHWAYS,anaTemp)
Post_InFlo(anaTemp)
}
})[3]
getDoParWorkers()
parallel_time














                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
