#' miRNA expression data.
#'
#' A dataset containing the miRNA expressopn levels for 144 TCGA samples (72 cancer samples
#' vs 72 normal samples, https://cancergenome.nih.gov/)
#'
#' @format A SummarizedExperiment object with 144 samples (rows) and 226 miRNAs (columns).
"miRExp"

#' mRNA expression data.
#'
#' A dataset containing the mRNA expressopn levels for 144 TCGA samples (72 cancer samples
#' vs 72 normal samples, https://cancergenome.nih.gov/)
#'
#' @format A SummarizedExperiment object with 144 samples (rows) and 500 mRNAs (columns).
"mRExp"

#' ceRNA expression data.
#'
#' A dataset containing the lncRNA expressopn levels for 144 TCGA samples (72 cancer samples
#' vs 72 normal samples, https://cancergenome.nih.gov/)
#'
#' @format A SummarizedExperiment object with 144 samples (rows) and 305 lncRNAs (columns).
"ceRExp"

#' miRNA-target (including miRNA-mRNA and miRNA-lncRNA) binding information.
#'
#' A dataset containing the miRNA-target binding information from miRTarBase
#' (http://mirtarbase.mbc.nctu.edu.tw/php/index.php), and LncBase
#' (http://carolina.imis.athena-innovation.gr/diana_tools/web/index.php?r=lncbasev2/index)
#'
#' @format A SummarizedExperiment object with 29901 miRNA-target interactions.
"miRTarget"
