library(testthat)
library(miRSM)

# Load datasets
miRExpcsv <- system.file("extdata","miRExp.csv",package="miRSM")
miRExp <- read.csv(miRExpcsv, header=FALSE, sep=",")
miRExp <- data_tidy(miRExp)
ceRExpcsv <- system.file("extdata","ceRExp.csv",package="miRSM")
ceRExp <- read.csv(ceRExpcsv, header=FALSE, sep=",")
ceRExp <- data_tidy(ceRExp)
mRExpcsv <- system.file("extdata","mRExp.csv",package="miRSM")
mRExp <- read.csv(mRExpcsv, header=FALSE, sep=",")
mRExp <- data_tidy(mRExp)
miRTargetcsv <- system.file("extdata","miRTarget.csv",package="miRSM")
miRTarget <- read.csv(miRTargetcsv, header=FALSE, sep=",")

# Identify gene co-expression modules using WGCNA method
modulegenes_WGCNA <- module_WGCNA(ceRExp, mRExp)


test_that("Test miRSM", {
    expect_equal(module_WGCNA(ceRExp, mRExp), modulegenes_WGCNA)    
})
