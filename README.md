# miRSM R package

# Introduction
This package provides several functions to study miRNA sponge modules.

# Installation
```{r echo=FALSE, results='hide', message=FALSE}
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("miRSM")
```

# A quick example to use miRSM package
```{r echo=FALSE, results='hide', message=FALSE}
# Load miRSM package 
suppressMessages(library(miRSM))

# Load BRCA sample data
data(BRCASampleData)

# Identify ceRNA-mRNA co-expression modules using WGCNA
modulegenes_WGCNA <- module_WGCNA(ceRExp[, seq_len(150)], 
                                  mRExp[, seq_len(150)])

# Identify miRNA sponge modules using sensitivity RV coefficient (SRVC)
miRSM_WGCNA_SRVC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                        modulegenes_WGCNA, method = "SRVC",
                        SMC.cutoff = 0.01, RV_method = "RV")
                        
# Functional analysis of miRNA sponge modules
miRSM_WGCNA_SRVC_genes <- miRSM_WGCNA_SRVC[[2]]
miRSM_WGCNA_SRVC_FEA <- module_FA(miRSM_WGCNA_SRVC_genes,
                                Analysis.type ="FEA")
miRSM_WGCNA_SRVC_DEA <- module_FA(miRSM_WGCNA_SRVC_genes, 
                                Analysis.type = "DEA")
                                
# Cancer enrichment analysis of miRNA sponge modules
miRSM.CEA.pvalue <- module_CEA(ceRExp, mRExp, BRCA_genes, miRSM_WGCNA_SRVC_genes)

# Validation of miRNA sponge interactions in miRNA sponge modules
library(miRspongeR)
Groundtruthcsv <- system.file("extdata", "Groundtruth.csv", package="miRspongeR")
Groundtruth <- read.csv(Groundtruthcsv, header=TRUE, sep=",")
miRSM.Validate <- module_Validate(miRSM_WGCNA_SRVC_genes, Groundtruth)

# Co-expression analysis of miRNA sponge modules
miRSM_WGCNA_Coexpress <-  module_Coexpress(ceRExp, mRExp, miRSM_WGCNA_SRVC_genes, resample = 10, method = "mean")

# miRNA distribution analysis of sharing miRNAs
miRSM_WGCNA_share_miRs <-  share_miRs(miRExp, ceRExp, mRExp, miRTarget, miRSM_WGCNA_SRVC_genes)
miRSM_WGCNA_miRdistribute <- module_miRdistribute(miRSM_WGCNA_share_miRs)

# Predict miRNA-target interactions
miRSM_WGCNA_miRtarget <- module_miRtarget(miRSM_WGCNA_share_miRs, miRSM_WGCNA_SRVC_genes)

# Identify miRNA sponge interactions
miRSM_WGCNA_miRsponge <- module_miRsponge(ceRExp, mRExp, miRSM_WGCNA_SRVC_genes)

```

# License
GPL-3
