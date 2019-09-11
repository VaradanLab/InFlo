#!/usr/bin/env Rscript

######################!/usr/bin/env Rscript

################################
######################################################################################################################################################################
##The Following fucntion checks for all the required files. 
File_chk <- function()
{
  Run_Test <<- ifelse(dir.exists(InFlo_Home),TRUE,FALSE)
  InFlo_core_file <<- paste(InFlo_Home,"/InFlo/InFlo",sep="")
  Run_Test <<- ifelse(file.exists(InFlo_core_file),TRUE,FALSE)
  #Interaction_Dir <<- paste(InFlo_Home,"/Support_Files/Interactions/",sep="")
  #Interaction_Files <<- list.files(Interaction_Dir)
  #Interaction_names <<- read.delim(paste(InFlo_Home,"/Support_Files/names.tab",sep=""),header = T,sep="\t")
  #Run_Test <<- ifelse(length(Interaction_Files)>0,TRUE,FALSE)
  Run_Test <<- ifelse(file.exists(PATHWAY_INFORMATION),TRUE,FALSE)
  Run_Test <<- ifelse(dir.exists(PATHWAYS_DIR),TRUE,FALSE)
  Genes_Info <<- read.delim(paste(InFlo_Home,"/Support_Files/GENE_NAMES_CHECK.txt",sep=""),sep="\t",header = T, check.names = F, stringsAsFactors = F)
  Genes_Info <- Genes_Info[!duplicated(Genes_Info),]
  
  
  #Pathway_Comps <<- paste(InFlo_Home,"/Support_Files/PATHWAYS_COMPONENT_NETWORK.txt",sep="")
  #PATHWAYS_COMPONENT_NETWORK <<- read.table(Pathway_Comps,sep="\t",header=T,check.names = F,stringsAsFactors = F)
  #Pathway_Interactions <<- paste(InFlo_Home,"/Support_Files/PATHWAYS_INTERACTION_NETWORK.txt",sep="")
  #PATHWAYS_INTERACTION_NETWORK <<- read.table(Pathway_Interactions,header=T,check.names = F,stringsAsFactors = F)
  #Pathways <<- paste(InFlo_Home,"/pathways.zip",sep="")
  #Run_Test <<- ifelse(file.exists(Pathway_Comps),TRUE,FALSE)
  #Run_Test <<- ifelse(file.exists(Pathway_Interactions),TRUE,FALSE)
  ana <<- paste(InFlo_Home,"Analysis",sep='/')
  Run_Test <<- ifelse(dir.exists(ana),TRUE,FALSE)
  Run_Test <<- ifelse(file.exists(CNV_FILE),TRUE,FALSE)
  Run_Test <<- ifelse(file.exists(GE_FILE),TRUE,FALSE)
  #Run_Test <<- ifelse(file.exists(METH_FILE),TRUE,FALSE)
  return(Run_Test)
}

####################################################################################################
##The following function creates the analysis directory, copies the input_files and pathways to the new analysis directory with a unique time stamp.  
Dir_create <- function(Run_Test)
{
  if(Run_Test){
    if(is.null(PROJ_NAME))
    {
    anaPath <<- paste(ana,(paste("Analysis",toString(strftime(Sys.time(),format="%d_%m_%Y_%H_%M")),sep='_')),sep ="/")
    dir.create(anaPath,showWarnings = T)
    }else{
      anaPath <<- paste(ana,PROJ_NAME,sep ="/")
      dir.create(anaPath,showWarnings = T)
    }
    anaInput <<- paste(anaPath,"Input",sep='/')
    dir.create(anaInput,showWarnings = T)
    file.copy(c(CNV_FILE,GE_FILE),anaInput)
    anaTemp <<-paste(anaPath,"temp",sep='/')
    dir.create(anaTemp,showWarnings = T)
    #file.copy(Pathways,anaTemp)
    PathDir <<- paste(anaPath,"/pathways",sep="")
    dir.create(PathDir,showWarnings = T)
    #dir.create(anaTemp,showWarnings = T)
    #unzip(zipfile = Pathways,exdir = anaTemp)
    ResPath <<- paste(anaPath,"/Results",sep="")
    ResPathPath <<- paste(ResPath,"/PATHWAYS",sep="")
    dir.create(ResPath,showWarnings = T)
    dir.create(ResPathPath,showWarnings = T)
    GE_FILE <<- paste(anaInput,as.character(unlist(strsplit(GE_FILE,split = "/"))[length(unlist(strsplit(GE_FILE,split = "/")))]),sep="/")
    CNV_FILE <<- paste(anaInput,as.character(unlist(strsplit(CNV_FILE,split = "/"))[length(unlist(strsplit(CNV_FILE,split = "/")))]),sep="/")
  }
}

####################################################################################################
# The following function downloads all the required R packages. 
dwnPack <- function(x)
{
  if(!require(x,character.only = TRUE))
  {
    install.packages(x,dependencies =TRUE)
    if(!require(x,character.only = TRUE))
    {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
        BiocManager::install(x,dependencies = T)
      if(is.element(x, installed.packages()[,1])) stop("Package not found")
    }
  }
}
####################################################################################################
# The following function checks the user input for the config file. and returns possible options. 
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
####################################################################################################
###################################################################################################
for(pack in c("marray","stringr","plyr","multtest","permute","IRanges","GenomicRanges","TCGAbiolinks","tidyr","dplyr","survival","DESeq2","stats","modeest","mixtools","doParallel"))
{
  tryCatch({
    dwnPack(pack)}, error=function(e){cat("Error :",conditionMessage(e)," ")})
}




# dwnPack("marray")
# dwnPack("stringr")
# dwnPack("plyr")
# dwnPack("multtest")
# dwnPack("permute")
# dwnPack("IRanges")
# dwnPack("GenomicRanges")
# dwnPack("TCGAbiolinks")
# dwnPack("tidyr")
# dwnPack("dplyr")
# dwnPack("survival")
# dwnPack("DESeq2")
# dwnPack("stats")
# dwnPack("modeest")
# dwnPack("mixtools")
# dwnPack("doParallel")
###################################################################################################
# Function to read the Data files. and check for If the file consists of duplicate sample names or gene names. 
Data_Read <- function(X){
  Data2 <- read.table(X,sep="\t",header=T,check.names = F, stringsAsFactors = F)
  # if(class(Data2)=="try-error"){
  #   Data2 <- NULL
  # }
  #Data2 <- read.table(X,sep="\t",header=T,check.names = F, stringsAsFactors = F)
  # if(length(Data2[,'Gene_Name'])==length(unique(Data2[,'Gene_Name']))){
  #   rownames(Data2) <- Data2[,'Gene_Name']
  #   Data2[,"Gene_Name"] <- NULL
  #   return(Data2)
  # }else{
  #   print("Duplicate_Gene_Names")
  # }
}

####################################################################################################
# Function to perform Wilcox test between each tumor sample with the provided Normal samples. 
Wilcox_Test <- function(X){
  Dat <- X
  tumor_samples <- intersect(tumor_samples,colnames(Dat))
  normal_samples <- intersect(normal_samples,colnames(Dat))
  if(length(normal_samples)>=3){
    tumor_matrix <- Dat[,tumor_samples]
    normal_matrix <- Dat[,normal_samples]
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
####################################################################################################
PATH_PROCESS <- function(X,Y){
  
  Run_PathInfo <- X
  Run_PathDir <- Y
  Run_ReqdPaths <- read.delim(Run_PathInfo,sep="\t",header = T,check.names = F,stringsAsFactors = F)
  Run_ReqdPaths$Old_Location <- paste(Run_PathDir,Run_ReqdPaths[,"FileName"],sep="/")
  Run_ReqdPaths$New_Location <- paste(PathDir,"/","pid_",Run_ReqdPaths[,"PID"],"_pathway.tab",sep="")
  RUN_COMPONENTS <<- NULL
  RUN_INTERACTIONS <<- NULL
  RUN_PIDS <<- NULL
  AVAIL_PATHS <<- NULL
  for(i in 1:length(Run_ReqdPaths[,1])){
    tryCatch({
      if(file.copy(from = Run_ReqdPaths[i,"Old_Location"],to = Run_ReqdPaths[i,"New_Location"],overwrite = F,recursive = F)){
        print(i)
        RUN_PIDS <<- c(RUN_PIDS,Run_ReqdPaths[i,"PID"])
        path <- Run_ReqdPaths[i,"New_Location"]
        print(path)
        max_count_field <- max(count.fields(path,sep="\t"),na.rm = T)
        Pathway_data <- read.delim(path,sep = "\t",header = F,fill = T,strip.white = T,stringsAsFactors = F,col.names =c(1:3))
        colnames(Pathway_data) <- c("c1","c2","c3")
        Pathway_data$c2 <- gsub(" ","_",Pathway_data$c2)
        Pathway_data$c1 <- gsub(" ","_",Pathway_data$c1)
        #colnames(Pathway)
        
        Pathway_comp <- Pathway_data[which(Pathway_data[,3]==""),c(2,1)]
        Pathway_comp <- Pathway_comp[!duplicated(Pathway_comp),]
        Pathway_comp <- cbind(Run_ReqdPaths[i,c(1:3)],rownames(Pathway_comp),Pathway_comp)
        colnames(Pathway_comp)[c(4,5,6)] <- c("Parent_Index","Parent","Parent_Type")
        write.table(Pathway_comp[which(Pathway_comp[,"Parent_Type"]=="protein"),"Parent"],paste(path,"_protein",sep=""),row.names = F,col.names = F,quote = F)
        RUN_COMPONENTS <- rbind(RUN_COMPONENTS,Pathway_comp)
        
        
        
        Pathway_interactions <- Pathway_data[which(Pathway_data[,3]!=""),c(1,3,2)]
        Pathway_interactions <- cbind(Run_ReqdPaths[i,c(1:3)],Pathway_interactions)
        colnames(Pathway_interactions)[c(4:6)] <- c("Parent","Interaction","Target")
        COMP <- Pathway_comp[,c("Parent","Parent_Type","Parent_Index")]
        Pathway_interactions <- plyr:::join(COMP,Pathway_interactions,by="Parent",type="right",match="all")
        colnames(COMP) <- c("Target","Target_Type","Target_Index")
        Pathway_interactions <- plyr:::join(COMP,Pathway_interactions,by="Target",type="right",match="all")
        Pathway_interactions <- Pathway_interactions[,c("PID","FileName","Pathway_Name","Parent","Parent_Type","Parent_Index","Interaction","Target","Target_Type","Target_Index")]
        RUN_INTERACTIONS <- rbind(RUN_INTERACTIONS,Pathway_interactions)
        
        AVAIL_PATHS <- rbind(AVAIL_PATHS,Run_ReqdPaths[i,])
      }
    }, error=function(e){cat("ERROR :",conditionMessage(e),"\n")})
  }
  write.table(RUN_INTERACTIONS,paste(anaInput,"/INTERACTIONS.txt",sep=""),sep="\t",row.names = F,col.names = T)
  write.table(RUN_COMPONENTS,paste(anaInput,"/COMPONENTS.txt",sep=""),sep="\t",row.names = F,col.names = T)
  write.table(AVAIL_PATHS,paste(anaInput,"/AVAILABLE_PATHWAYS.txt",sep=""),sep="\t",row.names = F,col.names = T)
  return(RUN_PIDS)
}





####################################################################################################
## If user selects RNASeqV2 as true, the following function does a DeSeq test on the expression data. 
DeSEQ_TEST <- function(X){
  Dat <-X
  tumor_matrix <- Dat[,intersect(colnames(Dat),tumor_samples)]
  normal_matrix <- Dat[,intersect(colnames(Dat),normal_samples)]
  
  DEFSEQ_RES <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
  colnames(DEFSEQ_RES) = colnames(tumor_matrix)
  rownames(DEFSEQ_RES) = rownames(tumor_matrix)
  
  DEFSEQ_RES2 <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
  colnames(DEFSEQ_RES2) = colnames(tumor_matrix)
  rownames(DEFSEQ_RES2) = rownames(tumor_matrix)
  
  for(i in 1:length(tumor_matrix[1,])){
    tryCatch({
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
    }, error=function(e){cat("Error :",conditionMessage(e)," ")})
  }
  return(DEFSEQ_RES)
}

####################################################################################################
## The following funtion is a pre Inflo check that arranges the expression and CNV data in the same 
## order of sample names and gene names and creates the Pathways specific *_mRNA.tab files and
## *_Genome.tab files.
PRE_INFLO <- function(X,Y,Z){
  exp <- X
  exp.matrix <- as.matrix(exp)
  # Read in Fac data
  Reqd_Genes <- as.data.frame(rownames(exp))
  colnames(Reqd_Genes) <- "Gene"
  
  Reqd_Samples <- as.data.frame(colnames(exp))
  colnames(Reqd_Samples) <- "Sample"
  
  
  Fac <- Y
  Fac <- cbind(Gene=rownames(Fac),Fac)
  Fac <- plyr:::join(Reqd_Genes,as.data.frame(Fac),by="Gene",type="left",match="all")
  rownames(Fac) <- Fac[,"Gene"]
  Fac[,"Gene"] <- NULL
  Fac <- as.data.frame(t(Fac))
  Fac <- cbind(Sample=rownames(Fac),Fac)
  Fac <- plyr:::join(Reqd_Samples,as.data.frame(Fac),by="Sample",type="left",match="all")
  rownames(Fac) <- Fac[,"Sample"]
  Fac[,"Sample"] <- NULL
  Fac <- as.data.frame(t(Fac))
  k <- sapply(Fac, is.factor)
  Fac[k] <- lapply(Fac[k], function(x) as.numeric(as.character(x)))
  k <- NULL
  Fac[is.na(Fac)] <- 1
  Fac.matrix <- as.matrix(Fac)
  
  common_gene <- intersect(rownames(exp.matrix),rownames(Fac.matrix))
  common_Samples <- intersect(colnames(exp.matrix),colnames(Fac.matrix))
  
  exp.matrix2 <- exp.matrix[common_gene,common_Samples]
  exp.matrix2 <- exp.matrix2[order(rownames(exp.matrix2)),order(colnames(exp.matrix2))]
  exp.matrix <- as.data.frame(exp.matrix2)
  exp.matrix2 <- NULL
  
  Fac.matrix2 <- Fac.matrix[common_gene,common_Samples]
  Fac.matrix2 <- Fac.matrix2[order(rownames(Fac.matrix2)),order(colnames(Fac.matrix2))]
  Fac.matrix <- as.data.frame(Fac.matrix2)
  Fac.matrix2 <- NULL
  # Read in pathway genes
  # the arguement should be a readable file (pid_1_pathway.tab)
  
  pathway_dir <- Z
  pathway_files <- list.files(pathway_dir, pattern = "_pathway.tab$")
  
  for(i in 1:length(pathway_files)){
    print(i)
    pathway_name<-pathway_files[i]
    pathway_gene_file <- paste(pathway_name,"_protein",sep="")
    
    # Uncomment if protein files don't already exist
    #command<- paste("./process_pathway.sh",pathway_name,">",pathway_gene_file)
    #system(command)
    
    pathway_genes_nonunique <- read.delim(paste(pathway_dir,pathway_gene_file,sep=""),sep="\t",header = F,check.names = F,stringsAsFactors = F)
    colnames(pathway_genes_nonunique) <- "GENE"
    pathway_genes_nonunique$NEW_GENE <- sapply(pathway_genes_nonunique[,"GENE"],function(x) unlist(strsplit(as.character(x),split = c("-","_")))[1])
    colnames(pathway_genes_nonunique) <- c("Orig_Gene","GENE")
    pathway_genes_nonunique <- plyr:::join(pathway_genes_nonunique,Genes_Info,by="GENE",type="left",match="all")
    
    #pathway_genes<-matrix(unique(pathway_genes_nonunique[,1]),length(unique(pathway_genes_nonunique[,1])))
    ##############################################################################################
    exp.matrix.pathway <- exp.matrix[intersect(pathway_genes_nonunique[,"HGNC_SYMBOL"],rownames(exp.matrix)),]
    exp.matrix.pathway <- cbind(HGNC_SYMBOL=rownames(exp.matrix.pathway),exp.matrix.pathway)
    exp.matrix.pathway <- plyr:::join(as.data.frame(pathway_genes_nonunique),as.data.frame(exp.matrix.pathway),type="full",match="all",by="HGNC_SYMBOL")
    exp.matrix.pathway[,c("GENE","HGNC_SYMBOL")] <- NULL
    k <- sapply(exp.matrix.pathway, is.factor)
    exp.matrix.pathway[k] <- lapply(exp.matrix.pathway[k], function(x) as.numeric(as.character(x)))
    exp.matrix.pathway <- aggregate(exp.matrix.pathway[,c(2:length(exp.matrix.pathway[1,]))], by=list(exp.matrix.pathway$Orig_Gene),FUN=mean, na.rm=TRUE, na.action=NULL)
    rownames(exp.matrix.pathway) <- exp.matrix.pathway[,1]
    exp.matrix.pathway[,1] <- NULL
    exp.matrix.pathway <- exp.matrix.pathway[order(rownames(exp.matrix.pathway)),]
    exp.matrix.pathway <- as.data.frame(t(exp.matrix.pathway))
    #exp.matrix.pathway <- exp.matrix.pathway[,order(colnames(exp.matrix.pathway))]
    exp.matrix.pathway[,grep("-active",colnames(exp.matrix.pathway))] <- NA
    exp.matrix.pathway <- cbind(id=rownames(exp.matrix.pathway),exp.matrix.pathway)
    exp.matrix.pathway[exp.matrix.pathway =="NaN"] <- NA
    ###############################################################################################
    Fac.matrix.pathway <- Fac.matrix[intersect(pathway_genes_nonunique[,"HGNC_SYMBOL"],rownames(Fac.matrix)),]
    Fac.matrix.pathway <- cbind(HGNC_SYMBOL=rownames(Fac.matrix.pathway),Fac.matrix.pathway)
    Fac.matrix.pathway <- plyr:::join(as.data.frame(pathway_genes_nonunique),as.data.frame(Fac.matrix.pathway),type="full",match="all",by="HGNC_SYMBOL")
    Fac.matrix.pathway[,c("GENE","HGNC_SYMBOL")] <- NULL
    k <- sapply(Fac.matrix.pathway, is.factor)
    Fac.matrix.pathway[k] <- lapply(Fac.matrix.pathway[k], function(x) as.numeric(as.character(x)))
    Fac.matrix.pathway <- aggregate(Fac.matrix.pathway[,c(2:length(Fac.matrix.pathway[1,]))], by=list(Fac.matrix.pathway$Orig_Gene),FUN=mean, na.rm=TRUE, na.action=NULL)
    rownames(Fac.matrix.pathway) <- Fac.matrix.pathway[,1]
    Fac.matrix.pathway[,1] <- NULL
    Fac.matrix.pathway <- Fac.matrix.pathway[order(rownames(Fac.matrix.pathway)),]
    Fac.matrix.pathway <- as.data.frame(t(Fac.matrix.pathway))
    #Fac.matrix.pathway <- Fac.matrix.pathway[,order(colnames(Fac.matrix.pathway))]
    Fac.matrix.pathway[,grep("-active",colnames(Fac.matrix.pathway))] <- NA
    Fac.matrix.pathway <- cbind(id=rownames(Fac.matrix.pathway),Fac.matrix.pathway)
    Fac.matrix.pathway[Fac.matrix.pathway =="NaN"] <- NA
    #############################################################################################
    # # set row and column names for the pathway data
    # colnames(exp.matrix.pathway) <- colnames(exp.matrix)
    # rownames(exp.matrix.pathway) <- pathway_genes
    # 
    # colnames(Fac.matrix.pathway) <- colnames(Fac.matrix)
    # rownames(Fac.matrix.pathway) <- pathway_genes
    #############################################################################################
    
    # file names
    prefix<-strsplit(pathway_name,".tab")
    prefix <- paste(pathway_dir,prefix,sep="")
    mrna_file <-paste(prefix,"_mRNA.tab",sep="")
    copy_file<-paste(prefix,"_genome.tab",sep="") 

    # print to file
    write.table(exp.matrix.pathway,mrna_file,quote=F,sep='\t',col.names=T,row.names = F)
    write.table(Fac.matrix.pathway,copy_file,quote=F,sep='\t',col.names=T,row.names = F)
  }
}

####################################################################################################
## The following function runs the pathways specific InFlo core over the provided number of cores. 

InFlo <- function(X,Y){
  InFlo_core_file
  PATHWAY_DIR <- X
  RESULT_DIR <- Y
  for(i in RUN_PIDS){
    print(paste("Executive Pathway ",i))
    
    sys_cmd <- paste(InFlo_core_file," -p ",PATHWAY_DIR,"pid_",i,"_pathway.tab"," -c ",InFlo_Home,"/InFlo/em_simple.cfg -b ",PATHWAY_DIR,"pid_",i,"_pathway -r 1 > ",RESULT_DIR,"/samples_",i,"-inconsistentOK-test3.txt",sep="")
    print(sys_cmd)
    system(sys_cmd)
  }
}


####################################################################################################
## The following function comboines the pathway level Inflo results and combine them to give a
## final report that can be used for further downstream analysis. 
Post_InFlo <- function(X){
  Inflo_Res_Dir <- X
  tryCatch({
    Inflo_Res_Files <- list.files(path = Inflo_Res_Dir, pattern = "-inconsistentOK-test3.txt")
    Inflo_Res_ids <- sapply(Inflo_Res_Files,function(x) gsub("samples_","",gsub("-inconsistentOK-test3.txt","",as.character(x))))
    Inflo_Res_Files <- cbind(Inflo_Res_ids,Inflo_Res_Files)
    colnames(Inflo_Res_Files) <- c("id","File_name")
    rownames(Inflo_Res_Files) <- NULL
    #Int_names <- Interaction_names
    #Inflo_Res_Files <- plyr:::join(as.data.frame(Inflo_Res_Files),Int_names,by="id",type="inner",match="all")
    # l <- sapply(Inflo_Res_Files, is.factor)
    # Inflo_Res_Files[l] <- lapply(Inflo_Res_Files[l], as.character)
    res <- NULL
    for(i in 1:length(Inflo_Res_Files[,1])){
      print(i)
      tryCatch({
        pathway_file <- read.delim(file = paste(Inflo_Res_Dir,"/",Inflo_Res_Files[i,2],sep=""),sep ="\t", header =T, check.names=F,stringsAsFactors=F)
        colnames(pathway_file)[1] <- "Patient"
        pathway_file[,length(pathway_file[1,])] <- NULL
        
        Without_Norm_av <- aggregate(pathway_file,by=list(pathway_file$Patient),FUN=mean,na.rm = T)
        Without_Norm_av[,'Patient'] <- NULL
        Without_Norm_av[,'Sample.'] <- NULL
        colnames(Without_Norm_av)[1] <- "Patient" 
        # Without_Norm_av_colnames <- colnames(Without_Norm_av)
        # #Without_Norm_av_colnames <- gsub("X","",Without_Norm_av_colnames)
        # Without_Norm_av_colnames_2 <- Without_Norm_av_colnames
        # Without_Norm_av_colnames_2[c(2:length(Without_Norm_av_colnames_2))] <- floor(as.numeric(Without_Norm_av_colnames_2[c(2:length(Without_Norm_av_colnames_2))]))
        # colnames(Without_Norm_av) <- Without_Norm_av_colnames_2
        # Without_Norm_av[,'Patient'] <- as.character(Without_Norm_av[,'Patient'])
        for(j in 2:length(Without_Norm_av[1,])){
          Without_Norm_av[,j] <- as.numeric(as.character(Without_Norm_av[,j]))
        }
        av <- Without_Norm_av
        av[,"Sample#"] <- NULL
        #write.table(av,file = paste(Inflo_means_Dir,paste(Inflo_Res_Files[i,2],"mean",sep="_"),sep="/"),sep="\t",col.names = T, row.names = F, eol = "\n")
        rownames(av) <- av[,'Patient']
        av[,'Patient'] <- NULL
        binned_file <- as.data.frame(t(av))
        binned_file <- cbind(Inflo_Res_Files[i,1],rownames(binned_file),binned_file)
        colnames(binned_file)[c(1,2)] <- c("PID","Interaction")
        res <- rbind.data.frame(res,binned_file)
      }, error=function(e){cat("ERROR :",as.character(Inflo_Res_Files[i,1]),conditionMessage(e), "\n")})
    }
    INFLO_DATA <- NULL
    INFLO_DATA <- res
    INFLO_DATA <- as.data.frame(INFLO_DATA)
    #INFLO_DATA[,'Patient'] <- NULL
    rownames(INFLO_DATA) <- NULL
    INFLO_DATA <- INFLO_DATA[complete.cases(INFLO_DATA),]
    INFLO_DATA[,'Interaction'] <- as.character(INFLO_DATA[,'Interaction'])
    # INFLO_DATA$ID2 <- paste(INFLO_DATA[,'name'],INFLO_DATA[,'Interaction'],sep="*")
    # INFLO_DATA$ID2 <- gsub(" ","_",INFLO_DATA$ID2)
    # rownames(INFLO_DATA) <- INFLO_DATA[,'ID2']
    # INFLO_DATA[,'ID'] <- NULL
    R1 <- NULL
    R1 <- data.frame(do.call('rbind', strsplit(as.character(INFLO_DATA$Interaction),':',fixed=TRUE)))[c(1:2)]
    colnames(R1) <- c("Target","Parents")
    R1$Parents <- gsub(" ","_",R1$Parents)
    R1$Target <- gsub(" ","_",R1$Target)
    R1$Parents <- gsub(",_","*",R1$Parents)
    R1$Parents <- gsub("\\(","",R1$Parents)
    R1$Parents <- gsub("\\)","",R1$Parents)
    
    INFLO_DATA <- cbind(R1,INFLO_DATA)
    PATH_INFO <- read.delim(paste(anaInput,"/AVAILABLE_PATHWAYS.txt",sep=""),sep="\t",header = T,check.names = F,stringsAsFactors = F)
    PATH_INFO <- PATH_INFO[,c("PID","Pathway_Name")]
    INFLO_DATA <- plyr:::join(PATH_INFO,INFLO_DATA,by="PID",type="right",match="all")
    INFLO_DATA <- INFLO_DATA[,c("PID","Pathway_Name","Interaction","Target","Parents",colnames(INFLO_DATA)[c(6:length(INFLO_DATA[1,]))])]
    write.table(INFLO_DATA,paste(ResPath,"/InFLo_Results.txt",sep=""),sep="\t",row.names = F,col.names = T)
    #################################################################################################
    #INFLO_DATA <- read.delim("/Projects/BRC_Brain_Mets/BRC_BRAIN_CNS_METS_INFLO_RESULTS/Results/InFLo_Results.txt",sep="\t",header = T,check.names = F, stringsAsFactors = F)
    INFLO_DATA$Interaction <- sapply(INFLO_DATA$Interaction,function(x) gsub(" : ",":",as.character(x)))
    INFLO_DATA$Interaction <- sapply(INFLO_DATA$Interaction,function(x) gsub(": ",":",as.character(x)))
    INFLO_DATA$Interaction <- sapply(INFLO_DATA$Interaction,function(x) gsub(" :",":",as.character(x)))
    INFLO_DATA$Interaction <- sapply(INFLO_DATA$Interaction,function(x) gsub("::",":",as.character(x)))
    
    INFLO_DATA <- INFLO_DATA %>% mutate(Parents = strsplit(as.character(Parents), ",")) %>% unnest(Parents)
    INFLO_DATA <- INFLO_DATA[,c("PID","Pathway_Name","Interaction","Target","Parents",colnames(INFLO_DATA)[c(5:(length(INFLO_DATA[1,])-1))])]
    
    INFLO_DATA_PARENTS <- data.frame(do.call('rbind', str_split(as.character(INFLO_DATA$Parents),pattern = '\\*',n = 2)))
    colnames(INFLO_DATA_PARENTS) <- c("Parent","Parent_STATUS")
    
    
    #INFLO_DATA_PARENTS[] <- lapply(INFLO_DATA_PARENTS,as.character)
    
    
    INFLO_DATA_PARENTS$Parent_STATUS <- ifelse(INFLO_DATA_PARENTS[,"Parent"]==INFLO_DATA_PARENTS[,"Parent_STATUS"],NA,INFLO_DATA_PARENTS[,"Parent_STATUS"])
    INFLO_DATA_PARENTS$Parent_STATUS <- sapply(INFLO_DATA_PARENTS$Parent_STATUS,function(x) gsub("\\*","",as.character(x)))
    
    
    
    
    
    INFLO_DATA <- cbind(INFLO_DATA_PARENTS,INFLO_DATA)
    INFLO_DATA <- INFLO_DATA[,c("PID","Pathway_Name","Interaction","Target","Parent","Parent_STATUS",colnames(INFLO_DATA)[c(8:(length(INFLO_DATA[1,])))])]
    INFLO_DATA <- INFLO_DATA[!duplicated(INFLO_DATA),]
    
    
    write.table(INFLO_DATA,paste(ResPath,"/InFLo_Interaction_Results.txt",sep=""),sep="\t",row.names = F,col.names = T)
    
    
    #INFLO_DATA[,c("Parent_Type","Interaction","Pathway_Name")] <- NULL
    #INFLO_DATA <- cbind(ID=paste(INFLO_DATA[,"PID"],INFLO_DATA[,"Parent"],INFLO_DATA[,"Target"],sep="*"),INFLO_DATA)
    
    
    SIG_TAR <- INFLO_DATA[,c("Target",colnames(INFLO_DATA)[c(7:length(INFLO_DATA))])]
    SIG_TAR <- aggregate(.~Target, data=SIG_TAR, mean,na.action = na.omit)
    colnames(SIG_TAR)[1] <- "Node"
    
    SIG_PAR <- INFLO_DATA[,c("Parent",colnames(INFLO_DATA)[c(7:length(INFLO_DATA))])]
    SIG_PAR <- aggregate(.~Parent, data=SIG_PAR, mean,na.action = na.omit)
    colnames(SIG_PAR)[1] <-"Node"
    
    
    INFLO_DATA <- as.data.frame(rbind(SIG_TAR,SIG_PAR))
    INFLO_DATA <- INFLO_DATA[!duplicated(INFLO_DATA[,"Node"]),]
    #INFLO_DATA <- INFLO_DATA[,c("Node",BE_SAMPS,EAC_SAMPS)]
    write.table(INFLO_DATA,paste(ResPath,"/InFLo_Nodes_Results.txt",sep=""),sep="\t",row.names = F,col.names = T)
    #################################################################################################
    # 
    # INFLO_DATA_PLOT <- INFLO_DATA %>% mutate(Parents = strsplit(as.character(Parents), ",")) %>% unnest(Parents)
    # INFLO_DATA_PLOT <- INFLO_DATA_PLOT[,c("PID","Pathway_Name","Interaction","Target","Parents",colnames(INFLO_DATA_PLOT)[c(5:(length(INFLO_DATA_PLOT[1,])-1))])]
    # INFLO_DATA_PLOT_PARENTS <- data.frame(do.call('rbind', strsplit(as.character(INFLO_DATA_PLOT$Parents),'*',fixed=TRUE)))
    # colnames(INFLO_DATA_PLOT_PARENTS) <- c("Parent","Parent_STATUS")
    # INFLO_DATA_PLOT <- cbind(INFLO_DATA_PLOT_PARENTS,INFLO_DATA_PLOT)
    # INFLO_DATA_PLOT <- INFLO_DATA_PLOT[,c("PID","Pathway_Name","Interaction","Target","Parent","Parent_STATUS",colnames(INFLO_DATA_PLOT)[c(8:(length(INFLO_DATA_PLOT[1,])))])]
    # #RUN_COMPONENTS <- read.delim()
    # RUN_INTERACTIONS <- read.delim(paste(anaInput,"/INTERACTIONS.txt",sep=""),sep="\t",header = T,check.names = F,stringsAsFactors = F)
    # RUN_PIDS <- as.numeric(PATH_INFO[,1])
    # for(i in 1:length(RUN_PIDS)){
    #   print(i)
    #   Inflo_Path_Res <- INFLO_DATA_PLOT[which(INFLO_DATA_PLOT[,"PID"]==RUN_PIDS[i]),]
    #   Inflo_Path_Res <- Inflo_Path_Res[!duplicated(Inflo_Path_Res),]
    #   Inflo_Path_Res <- cbind(paste(Inflo_Path_Res[,'Parent'],Inflo_Path_Res[,'Target'],sep="*"),Inflo_Path_Res)
    #   colnames(Inflo_Path_Res)[1] <- "ID"
    #   
    #   Inflo_Path_Inter <- RUN_INTERACTIONS[which(RUN_INTERACTIONS[,"PID"]==RUN_PIDS[i]),]
    #   Inflo_Path_Inter <- Inflo_Path_Inter[!duplicated(Inflo_Path_Inter),c("Parent","Parent_Index","Parent_Type","Target","Target_Index","Target_Type","Interaction")]
    #   Inflo_Path_Inter <- cbind(paste(Inflo_Path_Inter[,'Parent'],Inflo_Path_Inter[,'Target'],sep="*"),Inflo_Path_Inter)
    #   colnames(Inflo_Path_Inter)[1] <- "ID"
    #   
    #   
    #   Inflo_Path_Res2 <- plyr:::join(Inflo_Path_Inter,Inflo_Path_Res,by="ID",type="right",match="all")
    #   Inflo_Path_Res2 <- Inflo_Path_Res2[complete.cases(Inflo_Path_Res2),]
    #   PATH_NET <- Inflo_Path_Res2[,c("Parent","Parent_Index","Parent_Type","Target","Target_Index","Target_Type","Interaction")]
    #   write.table(PATH_NET,paste(ResPathPath,"/",RUN_PIDS[i],"_NETWORK.txt",sep=""),sep="\t",row.names=F,col.names = T,quote = F)
    #   
    #   PATH_COMP_PARENTS <- Inflo_Path_Res2[,c("Parent","Parent_Type","Parent_Index",colnames(Inflo_Path_Res2)[c(15:(length(Inflo_Path_Res2[1,])))])]
    #   colnames(PATH_COMP_PARENTS)[c(1:3)] <- c("Component","Component_Type","Component_Index")
    #   PATH_COMP_TARGETS <- Inflo_Path_Res2[,c("Target","Target_Type","Target_Index",colnames(Inflo_Path_Res2)[c(15:(length(Inflo_Path_Res2[1,])))])]
    #   colnames(PATH_COMP_TARGETS)[c(1:3)] <- c("Component","Component_Type","Component_Index")
    #   
    #   PATH_COMP <- rbind(PATH_COMP_TARGETS,PATH_COMP_PARENTS)
    #   PATH_COMP <- PATH_COMP[!duplicated(PATH_COMP),]
    #   PATH_COMP2 <- aggregate(.~Component+Component_Type+Component_Index,data = PATH_COMP,FUN=mean, na.rm=TRUE)
    #   write.table(PATH_COMP2,paste(ResPathPath,"/",RUN_PIDS[i],"_COMPONENTS.txt",sep=""),sep="\t",row.names=F,col.names = T,quote = F)
    # }
    # 
    #return(INFLO_DATA)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
####################################################################################################


#####################Guassian-fit Modelling#########################################################
## We also tried fitting data over guassian models. but not using in the main pipepline. 
# Guass_Fit <- function(X){
#   Dat <- X
#   tumor_samples <- intersect(tumor_samples,colnames(Dat))
#   normal_samples <- intersect(normal_samples,colnames(Dat))
#   if(length(normal_samples)>=3){
#     tumor_matrix <- Dat[,tumor_samples]
#     normal_matrix <- Dat[,normal_samples]
#     #normal_matrix <- normal_matrix[rownames(tumor_matrix),]
#     
#     processed_tumor_matrix <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
#     
#     Dat <- cbind(tumor_matrix,normal_matrix)
#     processed_tumor_matrix <- matrix(NA,ncol=ncol(tumor_matrix),nrow=nrow(tumor_matrix))
#     for(l in 1:length(tumor_samples)){
#       print(l)
#       
#       fit_data <- cbind(tumor_matrix[,l],normal_matrix)
#       tryCatch({processed_tumor_matrix[,l] <- apply(Dat,1,function(x){2^as.numeric(normalmixEM(na.omit(unlist(fit_data)),k=2,maxit=50,epsilon=0.01)$loglik)})}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
#     }
#     colnames(processed_tumor_matrix) = colnames(tumor_matrix)
#     rownames(processed_tumor_matrix) = rownames(tumor_matrix)
#     
#     return(processed_tumor_matrix)
#     }
# }