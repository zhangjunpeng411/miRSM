## ----style, echo=FALSE, results="asis", message=FALSE----------------------
BiocStyle::markdown()
knitr::opts_chunk$set(tidy = FALSE,
    warning = FALSE,
    message = FALSE)

## ----echo=FALSE, results='hide', message=FALSE-----------------------------
library(miRSM)

## ---- eval=TRUE, include=TRUE----------------------------------------------
data(ceRExp)
data(mRExp)
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
modulegenes_WGCNA

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  data(ceRExp)
#  data(mRExp)
#  modulegenes_GFA <- module_GFA(ceRExp, mRExp)

## ---- eval=TRUE, include=TRUE----------------------------------------------
data(ceRExp)
data(mRExp)
modulegenes_igraph <- module_igraph(ceRExp[, seq_len(100)],
                                    mRExp[, seq_len(100)])
modulegenes_igraph

## ---- eval=TRUE, include=TRUE----------------------------------------------
data(ceRExp)
data(mRExp)
modulegenes_ProNet <- module_ProNet(ceRExp[, seq_len(100)],
                                    mRExp[, seq_len(100)])
modulegenes_ProNet

## ---- eval=TRUE, include=TRUE----------------------------------------------
data(ceRExp)
data(mRExp)
# Reimport NMF package to avoid conflicts with DelayedArray package
library(NMF)
modulegenes_NMF <- module_NMF(ceRExp[, seq_len(100)],
                              mRExp[, seq_len(100)])
modulegenes_NMF

## ---- eval=TRUE, include=TRUE----------------------------------------------
data(ceRExp)
data(mRExp)
modulegenes_biclust <- module_biclust(ceRExp[, seq_len(100)],
                                      mRExp[, seq_len(100)])
modulegenes_biclust

## ---- eval=TRUE, include=TRUE----------------------------------------------
data(ceRExp)
data(mRExp)
data(miRExp)
data(miRTarget)
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
# Identify miRNA sponge modules using cannonical correlation (CC)
miRSM_WGCNA_CC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                        modulegenes_WGCNA, nperms = 10, num_shared_miRNAs                         = 3, pvalue.cutoff = 0.05, 
                        method = "CC", CC.cutoff = 0.8)
# Identify miRNA sponge modules using sensitivity cannonical correlation (SCC) method
miRSM_WGCNA_SCC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                         modulegenes_WGCNA, nperms = 10, 
                         num_shared_miRNAs = 3, 
                         pvalue.cutoff = 0.05, 
                         method = "SCC", SCC.cutoff = 0.1)
miRSM_WGCNA_CC
miRSM_WGCNA_SCC

## ---- eval=FALSE, include=TRUE---------------------------------------------
#  data(ceRExp)
#  data(mRExp)
#  data(miRExp)
#  data(miRTarget)
#  modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)
#  # Identify miRNA sponge modules using cannonical correlation (CC)
#  miRSM_WGCNA_CC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#                          modulegenes_WGCNA, nperms = 10,
#                          method = "CC")
#  miRSM_WGCNA_CC_genes <- miRSM_WGCNA_CC[[2]]
#  miRSM_WGCNA_CC_FEA <- module_FA(miRSM_WGCNA_CC_genes,                                                    Analysis.type ="FEA")
#  miRSM_WGCNA_CC_DEA <- module_FA(miRSM_WGCNA_CC_genes,
#                                  Analysis.type = "DEA")

## --------------------------------------------------------------------------
sessionInfo()

