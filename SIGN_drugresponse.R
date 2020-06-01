rm(list = ls())
library(SIGN)
library(parallel)
###############################################
## Identifying TSC between a target sample
## and set of sensitive and resistant samples
## in each arm
###############################################
BootSim <- function(Ref_ExpMat, TargetSam_ExpMat,BootSamNum, BootNum){
  
  TSCVec <- c()
  for(BootIter in 1:BootNum){
    SampleInd <- sample.int(ncol(Ref_ExpMat), size = BootSamNum, replace = FALSE)
    TSCVec <- c(TSCVec, TSC(TargetSam_ExpMat, Ref_ExpMat[,SampleInd]))
  }
  
  return(TSCVec)
}

TSC_ResSens <- function(ExpMat, TargetSamInd, Labels, BootSamNum_Sens, BootSamNum_Res,BootNum){
  
  TSCVec_Resistant <- BootSim(Ref_ExpMat=ExpMat[,which(
    Labels == "Resistant" & c(1:ncol(ExpMat)) != TargetSamInd)], 
    TargetSam_ExpMat=matrix(ExpMat[,TargetSamInd], ncol =1),BootSamNum=BootSamNum_Res, BootNum)
  
  TSCVec_Sensitive <- BootSim(Ref_ExpMat=ExpMat[,which(
    Labels == "Sensitive" & c(1:ncol(ExpMat)) != TargetSamInd)],
    TargetSam_ExpMat=matrix(ExpMat[,TargetSamInd], ncol =1),BootSamNum=BootSamNum_Sens, BootNum)
  ###############
  TSCVec <- c(median(TSCVec_Resistant), median(TSCVec_Sensitive))
  names(TSCVec) <- c("Resistant", "Sensitive")
  return(TSCVec)
}

###########################################
## Importing pathways
###########################################
Pathway_GeneList <- readRDS("data/c5_all_entrez_ProcessPathways.rds")
Pathway_GeneList <- Pathway_GeneList[which(unlist(lapply(Pathway_GeneList, length)) >= 10 & 
                                             unlist(lapply(Pathway_GeneList, length)) <= 30)]

###########################################
# generating random data for 
# 100 samples and 10980 genes
# (number of genes in the pathways)
###########################################
NeoExpMat <- matrix(runif(10980*100, min = 0, max = 1), ncol = 100)
rownames(NeoExpMat) <- unique(unlist(Pathway_GeneList))
RespVec <- c("Resistant", "Sensitive")[sample.int(2, size = 100, replace = TRUE)]
###########################################
# subsetting the pathways for the analysis
# for faster analysis 
# (for random data used in this code)
###########################################
Pathway_GeneList <- Pathway_GeneList[1:5]
###########################################
## Calculating TSC of pathways between 
## samples of NeoALTTO
###########################################
BootSamNum_Sens <- 5
BootSamNum_Res <- 5
BootNum <- 100
##
CindexVec <- c()
PredMat <- c()
for(PathwayIter in 1:length(Pathway_GeneList)){ #
  print(PathwayIter)
  #######
  MatchedGeneInd <- which(rownames(NeoExpMat) %in% Pathway_GeneList[[PathwayIter]])
  TSCMat <- c()
  for(SamIter in 1:ncol(NeoExpMat)){ #
    TSCMat <- rbind(TSCMat, TSC_ResSens(ExpMat=NeoExpMat[MatchedGeneInd,], TargetSamInd=SamIter,
                                        Labels=RespVec, BootSamNum_Sens=BootSamNum_Sens,
                                        BootSamNum_Res=BootSamNum_Res,
                                        BootNum=BootNum))
  }
  rownames(TSCMat) <- RespVec
  #############
  PredVec <- as.factor(as.character(apply(TSCMat, 1, function(X){colnames(TSCMat)[which(X == max(X))[1]]})))
  CindexVec <- c(CindexVec, Hmisc::rcorr.cens((TSCMat[,2]-TSCMat[,1]), as.factor(RespVec))["C Index"])
}

names(CindexVec) <- names(Pathway_GeneList)

