#' ceRNA expression data 
#' 
#' @docType data 
#' @name ceRExp
#' @aliases ceRExp
#' @format ceRExp: A SummarizedExperiment object with 72 BRCA 
#' and 72 normal samples (rows) and 305 lncRNAs (columns). 
#' @details The matched breast invasive carcinoma (BRCA) miRNA, lncRNA 
#' and mRNA expression data is obtained from TCGA 
#' (http://cancergenome.nih.gov/). lncRNA expression data 
#' is regarded as ceRNA expression data. The data focuses on 
#' 72 individuals for which the complete sets of 
#' tumor and matched normal (i.e., normal tissue taken from the 
#' same patient) profiles are available. 
#' A lncRNA which has missing values in more than 
#' 10% of the samples is removed. The remaining missing values
#' are imputed using the k-nearest neighbours (KNN) algorithm
#' from the impute R package. We use the limma R package
#' to infer differentially expressed lncRNAs
#' between tumour and normal samples. After the analysis, 
#' we select top 305 lncRNAs which are differentially expressed
#' at a significant level (adjusted p-value < 1E-02,  
#' adjusted by Benjamini & Hochberg method).
NULL

#' mRNA expression data 
#' 
#' @docType data 
#' @name mRExp
#' @aliases mRExp
#' @format mRExp: A SummarizedExperiment 
#' object with 72 BRCA and 72 normal samples (rows) and 226 miRNAs 
#' (columns). 
#' @details The matched breast invasive carcinoma (BRCA) miRNA, lncRNA 
#' and mRNA expression data is obtained from TCGA 
#' (http://cancergenome.nih.gov/). The data focuses on 
#' 72 individuals for which the complete sets of 
#' tumor and matched normal (i.e., normal tissue taken from the 
#' same patient) profiles are available. 
#' A mRNA which has missing values in more than 
#' 10% of the samples is removed. The remaining missing values
#' are imputed using the k-nearest neighbours (KNN) algorithm
#' from the impute R package. We use the limma R package
#' to infer differentially expressed mRNAs
#' between tumour and normal samples. After the analysis, 
#' we select top 500 mRNAs which are differentially expressed
#' at a significant level (adjusted p-value < 1E-02,  
#' adjusted by Benjamini & Hochberg method).
NULL

#' miRNA expression data
#' 
#' @docType data 
#' @name miRExp
#' @aliases miRExp
#' @format miRExp: A SummarizedExperiment object with 72 BRCA and 72 normal
#' samples (rows) and 226 miRNAs (columns). 
#' @details The matched breast invasive carcinoma (BRCA) miRNA, lncRNA 
#' and mRNA expression data is obtained from TCGA 
#' (http://cancergenome.nih.gov/). The data focuses on 
#' 72 individuals for which the complete sets of 
#' tumor and matched normal (i.e., normal tissue taken from the 
#' same patient) profiles are available. 
#' A miRNA which has missing values in more than 
#' 10% of the samples is removed. The remaining missing values
#' are imputed using the k-nearest neighbours (KNN) algorithm
#' from the impute R package. We use the limma R package
#' to infer differentially expressed miRNAs, ceRNAs and mRNAs
#' between tumour and normal samples. After the analysis, 
#' we select top 226 miRNAs which are differentially expressed 
#' at a significant level (adjusted p-value < 1E-02,  
#' adjusted by Benjamini & Hochberg method).
NULL

#' miRNA-target ineractions 
#' 
#' @docType data 
#' @name miRTarget
#' @aliases miRTarget
#' @format miRTarget: A SummarizedExperiment object with
#' 29901 miRNA-target interactions.
#' @details The miRNA-target binding information is
#' from miRTarBase v7.0 (http://mirtarbase.mbc.nctu.edu.tw/php/index.php),
#' and LncBase v2.0
#' (http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=lncbasev2/index).
#' Among 226 miRNAs, 305 lncRNAs and 500 mRNAs 
#' which are differentially expressed, we obtain 29901 miRNA-target
#' interactions (including miRNA-lncRNA and miRNA-mRNA interactions).
#' @references Hastie T, Tibshirani R, Narasimhan B, Chu G. 
#' impute: Imputation for microarray data. 
#' R package version 1.54.0. doi:  10.18129/B9.bioc.impute.
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, 
#' Shi W, et al. limma powers differential expression 
#' analyses for RNA-sequencing and microarray studies. 
#' Nucleic Acids Res. 2015; 43(7):e47.
NULL