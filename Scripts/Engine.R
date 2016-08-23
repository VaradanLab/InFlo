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
  Cmd <- 
  
  
  
  
  
  
}





