rm(list = ls())
######################
UnivarBio <- function(ExpMat, RepVec, PermutNum=100){
  
  CorVec <- c()
  FDRVec <- c()
  for(FeatIter in 1:nrow(ExpMat)){
    FeatValVec <- as.numeric(ExpMat[FeatIter,])
    RespVec_Pred <- as.numeric(RepVec)
    #####
    RemInd <- which(is.na(FeatValVec))
    if(length(RemInd) > 0){
      FeatValVec <- FeatValVec[-RemInd]
      RespVec_Pred <- RepVec[-RemInd]
    }
    Cor <-  as.numeric(Hmisc::rcorr.cens(RespVec_Pred, FeatValVec)["C Index"])
    CorVec <- c(CorVec,Cor)
    ####
    set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )
    
    if(Cor > 0.5){
      FDR <- max(1/(PermutNum+1), length(which(unlist(lapply(1:PermutNum, function(Iter){
        as.numeric(Hmisc::rcorr.cens(RespVec_Pred[
          sample.int(n = length(RespVec_Pred),size = length(RespVec_Pred),
                     replace = FALSE)], FeatValVec)["C Index"])})) > Cor))/PermutNum)
    }else{
      FDR <- max(1/(PermutNum+1), length(which(unlist(lapply(1:PermutNum, function(Iter){
        as.numeric(Hmisc::rcorr.cens(RespVec_Pred[
          sample.int(n = length(RespVec_Pred),size = length(RespVec_Pred),
                     replace = FALSE)], FeatValVec)["C Index"])})) < Cor))/PermutNum)
    }
    FDRVec <- c(FDRVec, FDR)
  }
  
  names(CorVec) <- rownames(ExpMat)
  names(FDRVec) <- rownames(ExpMat)
  #####
  PredList <- list(Cindex = CorVec,
                   FDR = FDRVec)
  
  return(PredList)
}
############################################
# generating random data for 
# 100 samples and 1000 genes
############################################
NeoExpMat <- matrix(runif(1e5, min = 0, max = 1), ncol = 100)

ERVec <- c("NEGATIVE", "POSITIVE")[sample.int(2, size = 100, replace = TRUE)]

RespVec <- c("0", "1")[sample.int(2, size = 100, replace = TRUE)]
####################
Gene_PredList <- UnivarBio(ExpMat=NeoExpMat, RepVec=RespVec, PermutNum=100)
#####
Gene_PredList_ERNEGATIVE <- UnivarBio(ExpMat=NeoExpMat[,which(ERVec == "NEGATIVE")],
                                      RepVec=RespVec[which(ERVec == "NEGATIVE")], PermutNum=100)
#####
Gene_PredList_ERPOSITIVE <- UnivarBio(ExpMat=NeoExpMat[,which(ERVec == "POSITIVE")],
                                      RepVec=RespVec[which(ERVec == "POSITIVE")], PermutNum=100)
#####
CindexList <- list(Gene_Pred_All=Gene_PredList,
                   Gene_Pred_ERNEGATIVE=Gene_PredList_ERNEGATIVE,
                   Gene_PRed_ERPOSITIVE=Gene_PredList_ERPOSITIVE)

