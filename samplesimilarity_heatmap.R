rm(list = ls())
######################
CorType <- "spearman"
TargetArm <-  "LAPATINIB ALONE"  ##"LAPATINIB IN COMBINATION WITH TRASTUZUMAB" ##"TRASTUZUMAB ALONE"  ##
############################################
# generating random data for 
# 100 samples and 1000 genes
############################################
NeoExpMat <- matrix(runif(1e5, min = 0, max = 1), ncol = 100)
colnames(NeoExpMat) <- paste("sample", seq(1,100), sep = "_")

ERVec <- c("NEGATIVE", "POSITIVE")[sample.int(2, size = 100, replace = TRUE)]

RespVec <- c("Resistant", "Sensitive")[sample.int(2, size = 100, replace = TRUE)]


Categories <- data.frame(Response = factor(RespVec),
                         ERstatus = factor(ERVec))
rownames(Categories) <- colnames(NeoExpMat)

CategoriesColor <- list(Response = c(Sensitive = 'forestgreen',Resistant = 'darkred'),
                        ERstatus = c(POSITIVE="red",NEGATIVE="blue"))
##########################
CorMat <- Hmisc::rcorr(NeoExpMat, type = CorType)$r
##########################


pheatmap::pheatmap(CorMat,
                   annotation_names_col = FALSE,
                   annotation_col = Categories,
                   annotation_colors = CategoriesColor ,
                   fontsize = 10,
                   clustering_method="ward.D2",
                   clustering_distance="euclidean",
                   clustering_distance_cols = "euclidean",
                   breaks = seq((min(CorMat)-0.01), 
                                (1+0.01),by=0.01), 
                   color = colorRampPalette(c(rep("darkblue",3),
                                              rep(rgb(120/255, 100/255,0/255),1),
                                              rep("gold",3.5),
                                              rep("orange",2.5),
                                              rep("red",1.5),
                                              rep("darkred",1)))(
                                                length(seq((min(CorMat)-0.01), 
                                                           (1+0.01),
                                                           by=0.01))-1))

