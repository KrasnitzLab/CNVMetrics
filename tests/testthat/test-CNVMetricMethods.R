### Unit tests for CNVMetricsMethods.R functions

library(CNVMetrics)
library(rtracklayer)

chrInfo <- Seqinfo(seqnames=c("chr1", "chr2", "chr3"),
                   seqlengths=c(10000, 20000, 1500), 
                   isCircular=c(FALSE, FALSE, FALSE),
                   genome="Alien")

### Tests prepareInformation() results

context("prepareInformation() results")

test_that("prepareInformation() must return error when segDirectory is not a Seqinfo", {
    
    error_message <- "chrInfo must be a Seqinfo object."
    
    expect_error(prepareInformation(segDirectory = "test", chrInfo = 33), 
                 error_message)
    expect_error(prepareInformation(segDirectory = "test", chrInfo = "allo"), 
                 error_message)
})

test_that("prepareInformation() must return error when segmentWithHeader is not a logical", {
    
    error_message02 <- "segmentWithHeader must be a logical."
    
    expect_error(prepareInformation(segDirectory = "test", chrInfo = chrInfo, 
                                    segmentWithHeader = "hi"), error_message02)
    expect_error(prepareInformation(segDirectory = "test", chrInfo = chrInfo, 
                                    segmentWithHeader = 444), error_message02)
})