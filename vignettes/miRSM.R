## ----style, echo=FALSE, results="asis", message=FALSE-------------------------
BiocStyle::markdown()
knitr::opts_chunk$set(tidy = FALSE,
    warning = FALSE,
    message = FALSE)

## ----echo=FALSE, results='hide', message=FALSE--------------------------------
suppressPackageStartupMessages(library(miRSM))

## ----eval=TRUE, include=TRUE--------------------------------------------------
data(BRCASampleData)

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(80)], 
                                  mRExp[, seq_len(80)])
modulegenes_WGCNA

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  modulegenes_GFA <- module_GFA(ceRExp[seq_len(20), seq_len(15)],
#                                mRExp[seq_len(20), seq_len(15)],
#                                iter.max = 2600)
#  modulegenes_GFA

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)],
                                    mRExp[, seq_len(10)])
modulegenes_igraph

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_ProNet <- module_ProNet(ceRExp[, seq_len(10)],
                                    mRExp[, seq_len(10)])
modulegenes_ProNet

## ----eval=TRUE, include=TRUE--------------------------------------------------
# Reimport NMF package to avoid conflicts with DelayedArray package
library(NMF)
modulegenes_NMF <- module_NMF(ceRExp[, seq_len(10)],
                              mRExp[, seq_len(10)])
modulegenes_NMF

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_clust <- module_clust(ceRExp[, seq_len(30)],
                                  mRExp[, seq_len(30)])
modulegenes_clust

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_biclust <- module_biclust(ceRExp[, seq_len(30)],
                                      mRExp[, seq_len(30)])
modulegenes_biclust

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)], 
                                  mRExp[, seq_len(10)])
# Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
miRSM_igraph_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                        modulegenes_igraph,
                        num_shared_miRNAs = 3, pvalue.cutoff = 0.05, 
                        method = "SRVC", MC.cutoff = 0.8,
                        SMC.cutoff = 0.01, RV_method = "RV")
miRSM_igraph_SRVC

## ----eval=TRUE, include=TRUE--------------------------------------------------
nsamples <- 3
modulegenes_igraph_all <- module_igraph(ceRExp[, 151:300], mRExp[, 151:300])
modulegenes_WGCNA_exceptk <- lapply(seq(nsamples), function(i) 
                                  module_WGCNA(ceRExp[-i, seq(150)], mRExp[-i, seq(150)]))
miRSM_igraph_SRVC_all <- miRSM(miRExp, ceRExp[, 151:300], mRExp[, 151:300], miRTarget,
                               modulegenes_igraph_all, method = "SRVC",
                               SMC.cutoff = 0.01, RV_method = "RV")
miRSM_WGCNA_SRVC_exceptk <- lapply(seq(nsamples), function(i) miRSM(miRExp[-i, ],      
                                   ceRExp[-i, seq(150)], mRExp[-i,  seq(150)], miRTarget,
                                   modulegenes_WGCNA_exceptk[[i]], method = "SRVC",
                                   SMC.cutoff = 0.01, RV_method = "RV"))
Modulegenes_all <- miRSM_igraph_SRVC_all[[2]]
Modulegenes_exceptk <- lapply(seq(nsamples), function(i) miRSM_WGCNA_SRVC_exceptk[[i]][[2]])
Modules_SS <- miRSM_SS(Modulegenes_all, Modulegenes_exceptk)
Modules_SS

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(150)],
#                                    mRExp[, seq_len(150)])
#  # Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
#  miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
#                           modulegenes_WGCNA, method = "SRVC",
#                           SMC.cutoff = 0.01, RV_method = "RV")
#  miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
#  miRSM_WGCNA_SRVC_FEA <- module_FA(miRSM_WGCNA_SRVC_genes, Analysis.type = 'FEA')
#  miRSM_WGCNA_SRVC_DEA <- module_FA(miRSM_WGCNA_SRVC_genes, Analysis.type = 'DEA')

## ----eval=TRUE, include=TRUE--------------------------------------------------
modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(150)], 
                                  mRExp[, seq_len(150)])
# Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget,
                         modulegenes_WGCNA, method = "SRVC",
                         SMC.cutoff = 0.01, RV_method = "RV")
miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
miRSM.CEA.pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_WGCNA_SRVC_genes)
miRSM.CEA.pvalue

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Using the built-in groundtruth from the miRSM package
#  Groundtruthcsv <- system.file("extdata", "Groundtruth_high.csv", package="miRSM")
#  Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",")
#  # Using the identified miRNA sponge modules based on WGCNA and sensitivity RV coefficient (SRVC)
#  miRSM.Validate <- module_Validate(miRSM_WGCNA_SRVC_genes, Groundtruth)

## ----eval=TRUE, include=TRUE--------------------------------------------------
# Using the identified miRNA sponge modules based on WGCNA and sensitivity RV coefficient (SRVC)
miRSM_WGCNA_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_WGCNA_SRVC_genes, resample = 10, method = "mean", test.method = "t.test")
miRSM_WGCNA_Coexpress

## ----eval=TRUE, include=TRUE--------------------------------------------------
# Using the identified miRNA sponge modules based on WGCNA and sensitivity RV coefficient (SRVC)
miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, miRTarget, miRSM_WGCNA_SRVC_genes)
miRSM_WGCNA_miRdistribute <- module_miRdistribute(miRSM_WGCNA_share_miRs)
head(miRSM_WGCNA_miRdistribute)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Using the identified miRNA sponge modules based on WGCNA and sensitivity RV coefficient (SRVC)
#  miRSM_WGCNA_miRtarget <- module_miRtarget(miRSM_WGCNA_share_miRs, miRSM_WGCNA_SRVC_genes)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Using the identified miRNA sponge modules based on WGCNA and sensitivity RV coefficient (SRVC)
#  miRSM_WGCNA_miRsponge <- module_miRsponge(miRSM_WGCNA_SRVC_genes)

## -----------------------------------------------------------------------------
sessionInfo()

