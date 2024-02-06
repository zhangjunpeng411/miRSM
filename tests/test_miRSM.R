suppressPackageStartupMessages(library(testthat))
suppressPackageStartupMessages(library(GSEABase))
suppressPackageStartupMessages(library(miRSM))

# Load datasets
data(BRCASampleData)

# Identify gene co-expression modules using igraph method
modulegenes_igraph <- module_igraph(ceRExp[, seq_len(10)], 
    mRExp[, seq_len(10)])


test_that("Test miRSM", {
    expect_equal(geneIds(module_igraph(ceRExp[, seq_len(10)], 
        mRExp[, seq_len(10)])), geneIds(modulegenes_igraph))
})
