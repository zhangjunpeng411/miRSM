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

# Identify miRNA sponge modules using cannonical correlation (CC)
miRSM_WGCNA_CC <- miRSM(miRExp, ceRExp, mRExp, miRTarget, 
                        modulegenes_WGCNA, nperms = 100, 
                        method = "CC")
                        
# Functional analysis of miRNA sponge modules
miRSM_WGCNA_CC_genes <- miRSM_WGCNA_CC[[2]]
miRSM_WGCNA_CC_FEA <- module_FA(miRSM_WGCNA_CC_genes,
                                Analysis.type ="FEA")
miRSM_WGCNA_CC_DEA <- module_FA(miRSM_WGCNA_CC_genes, 
                                Analysis.type = "DEA")
```

# License
GPL-3
