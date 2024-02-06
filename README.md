# miRSM R package

# Introduction
This package provides several utility functions to study miRNA sponge or ceRNA modules at single-sample and multi-sample levels, including popular methods for inferring gene modules (candidate miRNA sponge or ceRNA modules), and two functions to identify miRNA sponge modules at single-sample and multi-sample levels, as well as several functions to conduct modular analysis of miRNA sponge modules.

# Installation
```{r echo=FALSE, results='hide', message=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("miRSM")
```

# A quick example to use miRSM package
```{r echo=FALSE, results='hide', message=FALSE}
# Load miRSM package 
suppressPackageStartupMessages(library(miRSM))

# Load BRCA sample data
data(BRCASampleData)

# Identifying gene co-expression modules using WGCNA
modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(150)], 
                                  mRExp[, seq_len(150)])

# Identifying miRNA sponge modules using sensitivity RV coefficient (SRVC)
miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                        modulegenes_WGCNA, method = "SRVC",
                        SMC.cutoff = 0.01, RV_method = "RV")
                        
# Identifying sample-specific miRNA sponge modules
nsamples <- 3
modulegenes_igraph_all <- module_igraph(ceRExp[, 151:300], mRExp[, 151:300])
modulegenes_WGCNA_exceptk <- lapply(seq(nsamples), function(i) module_WGCNA(ceRExp[-i, seq(150)], mRExp[-i, seq(150)]))
miRSM_igraph_SRVC_all <- miRSM(miRExp, ceRExp[, 151:300], mRExp[, 151:300], miRTarget,
                               modulegenes_igraph_all, method = "SRVC",
                               SMC.cutoff = 0.01, RV_method = "RV")
miRSM_WGCNA_SRVC_exceptk <- lapply(seq(nsamples), function(i) miRSM(miRExp[-i, ], ceRExp[-i,                                    seq(150)], mRExp[-i,  seq(150)], miRTarget,
                                   modulegenes_WGCNA_exceptk[[i]], method = "SRVC",
                                   SMC.cutoff = 0.01, RV_method = "RV"))
Modulegenes_all <- miRSM_igraph_SRVC_all[[2]]
Modulegenes_exceptk <- lapply(seq(nsamples), function(i) miRSM_WGCNA_SRVC_exceptk[[i]][[2]])
Modules_SS <- miRSM_SS(Modulegenes_all, Modulegenes_exceptk)
                        
# Functional analysis of miRNA sponge modules
miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
miRSM_WGCNA_SRVC_FEA <- module_FA(miRSM_WGCNA_SRVC_genes,
                                Analysis.type ="FEA")
miRSM_WGCNA_SRVC_DEA <- module_FA(miRSM_WGCNA_SRVC_genes, 
                                Analysis.type = "DEA")
                                
# Cancer enrichment analysis of miRNA sponge modules
miRSM.CEA.pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_WGCNA_SRVC_genes)

# Validation of miRNA sponge interactions in miRNA sponge modules
Groundtruthcsv <- system.file("extdata", "Groundtruth_high.csv", package="miRSM")
Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",")
miRSM.Validate <- module_Validate(miRSM_WGCNA_SRVC_genes, Groundtruth)

# Co-expression analysis of miRNA sponge modules
miRSM_WGCNA_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_WGCNA_SRVC_genes, resample = 10, method = "mean", test.method = "t.test")

# Distribution analysis of sharing miRNAs
miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, miRTarget, miRSM_WGCNA_SRVC_genes)
miRSM_WGCNA_miRdistribute <- module_miRdistribute(miRSM_WGCNA_share_miRs)

# Predicting miRNA-target interactions
miRSM_WGCNA_miRtarget <- module_miRtarget(miRSM_WGCNA_share_miRs, miRSM_WGCNA_SRVC_genes)

# Identifying miRNA sponge interactions
miRSM_WGCNA_miRsponge <- module_miRsponge(miRSM_WGCNA_SRVC_genes)

```

# License
GPL-3
