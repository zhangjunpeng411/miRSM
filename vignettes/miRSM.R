## ----style, echo=FALSE, results="asis", message=FALSE----------------------
BiocStyle::markdown()
knitr::opts_chunk$set(tidy = FALSE,
    warning = FALSE,
    message = FALSE)

## ----echo=FALSE, results='hide', message=FALSE-----------------------------
library(miRSM)

## ---- eval=TRUE, include=TRUE----------------------------------------------
ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
ceRExp <- data_tidy(ceRExp)
mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
mRExp <- data_tidy(mRExp)
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
modulegenes_WGCNA

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
#  ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
#  ceRExp <- data_tidy(ceRExp)
#  mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
#  mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
#  mRExp <- data_tidy(mRExp)
#  modulegenes_GFA <- module_GFA(ceRExp, mRExp)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
#  ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
#  ceRExp <- data_tidy(ceRExp)
#  mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
#  mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
#  mRExp <- data_tidy(mRExp)
#  modulegenes_igraph <- module_igraph(ceRExp, mRExp)

## ---- eval=TRUE, include=TRUE----------------------------------------------
ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
ceRExp <- data_tidy(ceRExp)
mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
mRExp <- data_tidy(mRExp)
modulegenes_ProNet <- module_ProNet(ceRExp, mRExp)
modulegenes_ProNet

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
#  ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
#  ceRExp <- data_tidy(ceRExp)
#  mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
#  mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
#  mRExp <- data_tidy(mRExp)
#  # Reimport NMF package to avoid conflicts with DelayedArray package
#  library(NMF)
#  modulegenes_NMF <- module_NMF(ceRExp, mRExp)

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
#  ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
#  ceRExp <- data_tidy(ceRExp)
#  mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
#  mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
#  mRExp <- data_tidy(mRExp)
#  modulegenes_biclust <- module_biclust(ceRExp, mRExp)
#  modulegenes_biclust

## ---- eval=TRUE, include=TRUE----------------------------------------------
miRExpcsv <- system.file("extdata","miRExp.csv",package="miRSM")
miRExp <- read.csv(miRExpcsv, header=FALSE, sep=",")
miRExp <- data_tidy(miRExp)
ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
ceRExp <- data_tidy(ceRExp)
mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
mRExp <- data_tidy(mRExp)
miRTargetcsv <- system.file("extdata","miRTarget.csv",package="miRSM")
miRTarget <- read.csv(miRTargetcsv, header=FALSE, sep=",")
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
# Identify miRNA sponge modules using cannonical correlation (CC)
miRSM_WGCNA_CC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                        modulegenes_WGCNA, nperms = 10, num_shared_miRNAs = 3, 
                        pvalue.cutoff = 0.05, method = "CC", CC.cutoff = 0.8)
# Identify miRNA sponge modules using sensitivity cannonical correlation (SCC) method
miRSM_WGCNA_SCC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                         modulegenes_WGCNA, nperms = 10, num_shared_miRNAs = 3, 
                         pvalue.cutoff = 0.05, method = "SCC", SCC.cutoff = 0.1)
miRSM_WGCNA_CC
miRSM_WGCNA_SCC

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  miRExpcsv <- system.file("extdata","miRExp.csv",package="miRSM")
#  miRExp <- read.csv(miRExpcsv, header=FALSE, sep=",")
#  miRExp <- data_tidy(miRExp)
#  ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
#  ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
#  ceRExp <- data_tidy(ceRExp)
#  mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
#  mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
#  mRExp <- data_tidy(mRExp)
#  miRTargetcsv <- system.file("extdata","miRTarget.csv",package="miRSM")
#  miRTarget <- read.csv(miRTargetcsv, header=FALSE, sep=",")
#  modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#  # Identify miRNA sponge modules using cannonical correlation (CC)
#  miRSM_WGCNA_CC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#                          modulegenes_WGCNA, nperms = 10, method = "CC")
#  miRSM_WGCNA_CC_genes <- miRSM_WGCNA_CC[[2]]
#  miRSM_WGCNA_CC_FEA <- module_FA(miRSM_WGCNA_CC_genes, Analysis.type = "FEA")
#  miRSM_WGCNA_CC_DEA <- module_FA(miRSM_WGCNA_CC_genes, Analysis.type = "DEA")

## --------------------------------------------------------------------------
sessionInfo()

