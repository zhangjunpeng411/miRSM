library(testthat)
library(miRSM)

# Load datasets
data(ceRExp)
data(mRExp)
data(miRExp)
data(miRTarget)

# Identify gene co-expression modules using WGCNA method
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)


test_that("Test miRSM", {
    expect_equal(module_WGCNA(ceRExp, mRExp), modulegenes_WGCNA)
})
