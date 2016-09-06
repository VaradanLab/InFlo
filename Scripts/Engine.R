PATHWAYS_INTERACTION_NETWORK <- read.table("/Projects/InFlo_Analysis/Inflo-code/Pathways/PATHWAYS_INTERACTION_NETWORK.txt",sep="\t",header=T,check.names = F,stringsAsFactors = F)
PATHWAYS_COMPONENT_NETWORK <- read.table("/Projects/InFlo_Analysis/Inflo-code/Pathways/PATHWAYS_COMPONENT_NETWORK.txt",sep="\t",header=T,check.names = F,stringsAsFactors = F)

Interaction_Dir <- "/Projects/InFlo_Analysis/Inflo-code/ovarian_cancer_Inflo_data/Interactions"
Interaction_names <- read.delim("/Projects/InFlo_Analysis/Inflo-code/Pathways/names.tab",header = T,sep="\t")




###########~~~~~~~~~Download Required Package~~~~~~~~~~~~~~~~~########################
dwnPack <- function(x)
{
  if(!require(x,character.only = TRUE))
  {
    install.packages(x,dep=TRUE)
    if(!require(x,character.only = TRUE))
    {
      source("https://bioconductor.org/biocLite.R")
      biocLite(x)
      if(is.element(x, installed.packages()[,1])) stop("Package not found")
    }
  }
}

####################################################################################################
###########################################~~~~Run_Check~~~~~##############################################

#!/usr/bin/env Rscript

##Note : Check the directories again and again for any further changes. 

######################################################################################################################################################################

run_chk <- function()
{
  args <- commandArgs(trailingOnly = TRUE)
  print(args)
  
  #args <- c("-R","/Projects/InFlo/scripts/InFlo_CONF.txt")
  
  len_args <- length(args)
  if(len_args == 0)
  {
    writeLines(" Usage : Rscript Inflo_Run.R [-R] [InFlo_RUN_CONF(follow the template for creating CONF FILE : <Home>/inFlo/Run_Conf_FIles/inFlo_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
  }else{
    if(args[1]=="-H")
    {
      InFloHome <<- getwd()
      read_me <<- unlist(strsplit(InFloHome,split='/'))
      read_me <<- paste(read_me[-length(read_me)],collapse='/')
      writeLines(" Usage : Rscript Inflo_Run.R [-R] [InFlo_RUN_CONF(follow the template for creating CONF FILE : <Home>/inFlo/Run_Conf_FIles/inFlo_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
      system(paste("cat ",read_me,"/README.txt",sep=""))
      return(FALSE)
    }else if(args[1]=="-R")
    {
      if(len_args!=2)
      {
        writeLines(" Usage : Rscript Inflo_Run.R [-R] [InFlo_RUN_CONF(follow the template for creating CONF FILE : <Home>/inFlo/Run_Conf_FIles/inFlo_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
        return(FALSE)
      } else
      {
        if(file.exists(args[2]))
        {
          source(args[2])
          return(TRUE)
        }else{
          print("InFlo RUN CONFIG FILE MISSING")
        }
      }
    }else{
      writeLines(" Usage : Rscript Inflo_Run.R [-R] [InFlo_RUN_CONF(follow the template for creating CONF FILE : <Home>/inFlo/Run_Conf_FIles/inFlo_RUN_CONF.txt)]\n\t\t\t[-H : for help]")
      return(FALSE)
    }
  }
}




Dir_Create <- function()
################################################################################################################################

dwnPack("marray")
dwnPack("stringr")
dwnPack("plyr")
dwnPack("multtest")
dwnPack("permute")
dwnPack("IRanges")
dwnPack("GenomicRanges")
dwnPack("TCGAbiolinks")
dwnPack("tidyr")
dwnPack("dplyr")
dwnPack("survival")
dwnPack("DESeq2")
dwnPack("stats")
dwnPack("modeest")
dwnPack("mixtools")
################################################################################################################################

Data_Read <- function(X){
  Data2 <- read.table(X,sep="\t",header=T,check.names = F, stringsAsFactors = F)
  if(length(Data2[,'Gene_Name'])==length(unique(Data2[,'Gene_Name']))){
    rownames(Data2) <- Data2[,'Gene_Name']
    Data2[,"Gene_Name"] <- NULL
    return(Data2)
  }else{
      print("Duplicate_Gene_Names")
  }
}
  
Wilcox_Test <- function(X){
  Dat <- X
  tumor_samples <- intersect(tumor_samples,colnames(Dat))
  normal_samples <- intersect(normal_samples,colnames(Dat))
  if(length(normal_samples)>=3){
    tumor_matrix <- Dat[,tumor_samples]
    normal_matrix <- Dat[,normal_samples]
    #normal_matrix <- normal_matrix[rownames(tumor_matrix),]
    
    processed_tumor_matrix <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
    
    Dat <- cbind(tumor_matrix,normal_matrix)
    processed_tumor_matrix <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
    
    
    for(l in 1:length(tumor_samples)){
      print(l)
      tryCatch({processed_tumor_matrix[,l] <- apply(Dat,1,function(x){if(is.na(x[tumor_samples[l]])){NA}else{normal_vector <- x[normal_samples]; result <- wilcox.test(as.numeric(normal_vector),x[tumor_samples[l]])$p.value; if(x[tumor_samples[l]] <= median(as.numeric(normal_vector),na.rm=T)){-1*result}else{result}}})}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    colnames(processed_tumor_matrix) = colnames(tumor_matrix)
    rownames(processed_tumor_matrix) = rownames(tumor_matrix)
    
    return(processed_tumor_matrix)
  }
}




DeSEQ_TEST <- function(X){
  
  Dat <-X
  tumor_matrix <- Dat[,tumor_samples]
  normal_matrix <- Dat[,normal_samples]
  
  
  DEFSEQ_RES <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
  colnames(DEFSEQ_RES) = colnames(tumor_matrix)
  rownames(DEFSEQ_RES) = rownames(tumor_matrix)
  
  DEFSEQ_RES2 <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
  colnames(DEFSEQ_RES2) = colnames(tumor_matrix)
  rownames(DEFSEQ_RES2) = rownames(tumor_matrix)
  
  for(i in 1:length(tumor_matrix[1,])){
    print(paste(i,colnames(tumor_matrix)[i],as.character(Sys.time()),sep="_"))
    reqd_mat <- cbind.data.frame(tumor_matrix[,i], normal_matrix)
    reqd_mat <- data.matrix(reqd_mat)
    reqd_mat <- round(reqd_mat)
    colnames(reqd_mat)[1] <- colnames(tumor_matrix)[i]
    condition_mad <- matrix('Nor',ncol= length(reqd_mat[1,]),nrow=1)
    condition_mad[1] <- "Tum"
    #condition_mad[normal_indices] <- "Nor"
    samples <- data.frame(row.names = colnames(reqd_mat), condition = t(condition_mad))
    bckCDS <- DESeqDataSetFromMatrix(countData = reqd_mat,colData = samples, design = ~condition)
    bckCDS_1 <- DESeq(bckCDS)
    bck_Res <- results(bckCDS_1)
    bck_Res_df <- as.data.frame(bck_Res)
    bck_Res_df2 <- cbind(bck_Res_df,sign(bck_Res_df[,'log2FoldChange'])*bck_Res_df[,'pvalue'])
    colnames(bck_Res_df2)[7] <- "new_Pval" 
    DEFSEQ_RES[,colnames(tumor_matrix)[i]] <- bck_Res_df2[,'new_Pval']
    DEFSEQ_RES2[,colnames(tumor_matrix)[i]] <- bck_Res_df2[,'log2FoldChange']
  }
  return(DEFSEQ_RES)
}


PRE_INFLO <- function(X,Y,Z){
  exp <- X
  fac <- Y
  
  exp <- X
  exp.matrix <- as.matrix(exp)
  
  # Read in Fac data
  Fac <- Y
  Fac.matrix <- as.matrix(Fac)
  
  
  common_gene <- intersect(rownames(exp.matrix),rownames(Fac.matrix))
  common_Samples <- intersect(colnames(exp.matrix),colnames(Fac.matrix))
  
  exp.matrix2 <- exp.matrix[common_gene,common_Samples]
  exp.matrix2 <- exp.matrix2[order(rownames(exp.matrix2)),order(colnames(exp.matrix2))]
  exp.matrix <- exp.matrix2
  
  Fac.matrix2 <- Fac.matrix[common_gene,common_Samples]
  Fac.matrix2 <- Fac.matrix2[order(rownames(Fac.matrix2)),order(colnames(Fac.matrix2))]
  Fac.matrix <- Fac.matrix2
  # Read in pathway genes
  # the arguement should be a readable file (pid_1_pathway.tab)
  
  
  
  pathway_dir <- Y
  pathway_files <- list.files(pathway_dir, pattern = "_pathway.tab$")
  
  
  for(i in 1:length(pathway_files)){
    print(i)
    pathway_name<-pathway_files[i]
    pathway_gene_file <- paste(pathway_name,"_protein",sep="")
    
    # Uncomment if protein files don't already exist
    #command<- paste("./process_pathway.sh",pathway_name,">",pathway_gene_file)
    #system(command)
    
    pathway_genes_nonunique <- read.table(paste(pathway_dir,pathway_gene_file,sep=""))
    pathway_genes<-matrix(unique(pathway_genes_nonunique[,1]),length(unique(pathway_genes_nonunique[,1])))
    
    exp.matrix.pathway <-sapply(pathway_genes,function(x) matrix(exp.matrix[match(unlist(x),rownames(exp.matrix)),],ncol=ncol(exp.matrix)))
    exp.matrix.pathway<-t(exp.matrix.pathway)
    Fac.matrix.pathway <-sapply(pathway_genes,function(x) matrix(Fac.matrix[match(unlist(x),rownames(Fac.matrix)),],ncol=ncol(Fac.matrix)))
    Fac.matrix.pathway<-t(Fac.matrix.pathway)
    
    # set row and column names for the pathway data
    colnames(exp.matrix.pathway) <- colnames(exp.matrix)
    rownames(exp.matrix.pathway) <- pathway_genes
    
    colnames(Fac.matrix.pathway) <- colnames(Fac.matrix)
    rownames(Fac.matrix.pathway) <- pathway_genes
    
    
    # file names
    prefix<-strsplit(pathway_name,".tab")
    prefix <- paste(pathway_dir,prefix,sep="")
    mrna_file <-paste(prefix,"_mRNA.tab",sep="")
    copy_file<-paste(prefix,"_genome.tab",sep="") 
    
    # file headers
    names <- append("id",pathway_genes[,1])
    write.table(t(names),mrna_file,quote=F,sep='\t',row.names=F,col.names=F)
    write.table(t(names),copy_file,quote=F,sep='\t',row.names=F,col.names=F)
    
    # print to file
    write.table(t(exp.matrix.pathway),mrna_file,quote=F,sep='\t',col.names=F,append=T)
    write.table(t(Fac.matrix.pathway),copy_file,quote=F,sep='\t',col.names=F,append=T)
  }
}

InFlo <- function(X,Y){
  PATHWAY_DIR <- X
  RESULT_DIR <- Y
  setwd(INFLO_DIR)
  #Cmd <- 
}

Post_Info <- function(X){
  Inflo_Res_Dir <- X
  Inflo_means_Dir <- paste(Inflo_Res_Dir,"_MEAN",sep="")
  Inflo_Res_Files <- list.files(path = Inflo_Res_Dir, pattern = "-inconsistentOK-blackballing.txt")
  
  
  for(i in 1:length(Inflo_Res_Files)){
    print(i)
    pathway_file <- read.csv(paste(Inflo_Res_Dir,Inflo_Res_Files[i],sep="/"), sep = "\t", header = T, check.names = T)
    colnames(pathway_file)[1] <- "Patient"
    pathway_file[,length(pathway_file[1,])] <- NULL
    
    Without_Norm_av <- aggregate(.~pathway_file$Patient, data = pathway_file, FUN=mean, simplify = T)
    Without_Norm_av[,'Patient'] <- NULL
    Without_Norm_av[,'Sample.'] <- NULL
    colnames(Without_Norm_av)[1] <- "Patient" 
    Without_Norm_av_colnames <- colnames(Without_Norm_av)
    Without_Norm_av_colnames <- gsub("X","",Without_Norm_av_colnames)
    Without_Norm_av_colnames_2 <- Without_Norm_av_colnames
    Without_Norm_av_colnames_2[c(2:length(Without_Norm_av_colnames_2))] <- floor(as.numeric(Without_Norm_av_colnames_2[c(2:length(Without_Norm_av_colnames_2))]))
    colnames(Without_Norm_av) <- Without_Norm_av_colnames_2
    Without_Norm_av[,'Patient'] <- as.character(Without_Norm_av[,'Patient'])
    for(j in 2:length(Without_Norm_av[1,])){
      Without_Norm_av[,j] <- as.numeric(Without_Norm_av[,j])
    }
    av <- Without_Norm_av
    write.table(av,file = paste(Inflo_means_Dir,paste(Inflo_Res_Files[i],"mean",sep="_"),sep="/"),sep="\t",col.names = T, row.names = F, eol = "\n")
  }
  
  
  dir_mean = Inflo_means_Dir
  mean_filnames = list.files(dir_mean, pattern = "blackballing.txt_mean")
  
  Interaction_files <- list.files(Interaction_Dir, pattern = "-interaction-names.txt")
  Interaction_names_2 <- gsub("-interaction-names.txt","",Interaction_files)
  
  
  Inflo_mean_files <- list.files(Inflo_means_Dir,pattern="-inconsistentOK-blackballing.txt_mean")
  Inflo_mean_files_2 <- gsub("-inconsistentOK-blackballing.txt_mean","",Inflo_mean_files)
  
  
  req_files <- intersect(Inflo_mean_files_2,Interaction_names_2)
  res <- NULL
  for(i in 1:length(req_files)){
    
    print(i)
    int_file <- paste(req_files[i],"-interaction-names.txt",sep="")
    Int_names <- read.delim(paste(Interaction_Dir,int_file,sep="/"), sep="\t", header = F)
    Int_names[length(Int_names)] <- NULL
    Int_names[1] <- "Patient"
    Int_names[2] <- NULL
    #Int_names <- as.list(Int_names)
    #Int_names <- as.character(Int_names)
    #binned_file <- read.csv(paste(dir_mean,gsub("-interaction-names.txt","-inconsistentOK-blackballing_mean.txt",int_file),sep="/"),sep="\t", header = T)
    binned_file <- read.csv(paste(dir_mean,paste(req_files[i],"-inconsistentOK-blackballing.txt_mean",sep=""),sep="/"),sep="\t", header = F)
    #binned_file[,1] <- Samp_names[,2]
    Int_names <- t(Int_names)
    binned_file <- t(binned_file)
    binned_file2 <- cbind(Int_names,binned_file)
    pathway_num <- gsub("-interaction-names.txt","",int_file)
    pathway_num <- gsub("samples_","",pathway_num)
    pathway_name <- as.character(Interaction_names[which(Interaction_names[,"id"]==pathway_num),'name'])
    colnames(binned_file2) <- binned_file2[1,]
    colnames(binned_file2)[1] <- "Interaction"
    #binned_file2[,'Patient'] <- NULL
    binned_file3 <- binned_file2[c(3:length(binned_file2[,1])),]
    binned_file3 <- cbind(pathway_name,pathway_num,int_file,binned_file3)
    #colnames(binned_file) <- as.character(Int_names[1,])
    res <- rbind(res,binned_file3)
    
  }
  
  INFLO_DATA <- res
  INFLO_DATA <- as.data.frame(INFLO_DATA)
  INFLO_DATA[,'Patient'] <- NULL
  rownames(INFLO_DATA) <- NULL
  INFLO_DATA <- INFLO_DATA[complete.cases(INFLO_DATA),]
  INFLO_DATA$ID <- paste(INFLO_DATA[,'pathway_name'],INFLO_DATA[,'Interaction'],sep="*")
  INFLO_DATA$ID <- gsub(" ","_",INFLO_DATA$ID)
  rownames(INFLO_DATA) <- INFLO_DATA[,'ID']
  INFLO_DATA[,'ID'] <- NULL
  R1 <- NULL
  R1 <- data.frame(do.call('rbind', strsplit(as.character(INFLO_DATA$Interaction),':',fixed=TRUE)))
  colnames(R1) <- c("Target","Parents")
  R1$Parents <- gsub(" ","_",R1$Parents)
  R1$Target <- gsub(" ","_",R1$Target)
  R1$Parents <- gsub(",_","*",R1$Parents)
  R1$Parents <- gsub("\\(","",R1$Parents)
  R1$Parents <- gsub("\\)","",R1$Parents)
  
  INFLO_DATA <- cbind(R1,INFLO_DATA)
  
  INFLO_DATA <- INFLO_DATA[,c("pathway_name","pathway_num","int_file","Interaction","Target","Parents",colnames(INFLO_DATA)[c(7:length(INFLO_DATA[1,]))])]
  
  return(INFLO_DATA)
}

Guass_Fit <- function(X){
  Dat <- X
  tumor_samples <- intersect(tumor_samples,colnames(Dat))
  normal_samples <- intersect(normal_samples,colnames(Dat))
  if(length(normal_samples)>=3){
    tumor_matrix <- Dat[,tumor_samples]
    normal_matrix <- Dat[,normal_samples]
    #normal_matrix <- normal_matrix[rownames(tumor_matrix),]
    
    processed_tumor_matrix <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
    
    Dat <- cbind(tumor_matrix,normal_matrix)
    processed_tumor_matrix <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
    for(l in 1:length(tumor_samples)){
      print(l)
      
      fit_data <- cbind(tumor_matrix[,l],normal_matrix)
      tryCatch({processed_tumor_matrix[,l] <- apply(Dat,1,function(x){normalmixEM(na.omit(unlist(fit_data)),k=2,maxit=50,epsilon=0.01)$loglik})}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    colnames(processed_tumor_matrix) = colnames(tumor_matrix)
    rownames(processed_tumor_matrix) = rownames(tumor_matrix)
    
    return(processed_tumor_matrix)
    }
}
