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
#' from miRTarBase v8.0 (http://mirtarbase.mbc.nctu.edu.tw/php/index.php),
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

#' BRCA genes 
#' 
#' @docType data 
#' @name BRCA_genes
#' @aliases BRCA_genes
#' @format BRCA_genes: A SummarizedExperiment object with
#' 4819 BRCA related genes (including lncRNAs and mRNAs).
#' @details The BRCA related lncRNAs are from LncRNADisease v2.0, Lnc2Cancer v2.0 
#' and MNDR v2.0. The BRCA related mRNAs are from DisGeNET v5.0 and COSMIC v86.
#' @references Bao Z, Yang Z, Huang Z, Zhou Y, Cui Q, Dong D. (2019) 
#' "LncRNADisease 2.0: an updated database of long non-coding RNA-associated diseases".
#' Nucleic Acids Res., 47(D1):D1034-D1037.
#' @references Cui T, Zhang L, Huang Y, Yi Y, Tan P, Zhao Y, Hu Y, Xu L, Li E, Wang D. 
#' (2018) "MNDR v2.0: an updated resource of ncRNA-disease 
#' associa-tions in mammals". Nucleic Acids Res., 46, D371-D374.
#' @references Gao Y, Wang P, Wang Y, Ma X, Zhi H, Zhou D, Li X, Fang Y, Shen W, Xu Y, 
#' Shang S, Wang L, Wang L, Ning S, Li X. (2019) "Lnc2Cancer v2.0: updated database of experimentally 
#' supported long non-coding RNAs in human cancers". Nucleic Acids Res., 47, D1028-D1033.
#' @references Forbes SA, Beare D, Boutselakis H, Bamford S, Bindal N, Tate J, Cole CG, 
#' Ward S, Dawson E, Ponting L, Stefancsik R, Harsha B, Kok CY, Jia M, Jubb H, Sondka Z, 
#' Thompson S, De T, Campbell PJ. (2017) "COSMIC: somatic cancer genetics at 
#' high-resolution". Nucleic Acids Res., 45, D777-D783
#' @references Pinero J, Bravo A, Queralt-Rosinach N, Gutierrez-Sacristan A, Deu-Pons J, 
#' Centeno E, Garcia-Garcia J, Sanz F, Furlong LI. 
#' (2017) "DisGeNET: a comprehensive platform integrating 
#' infor-mation on human disease-associated genes and variants". 
#' Nucleic Acids Res., 45, D833-D839.
NULL