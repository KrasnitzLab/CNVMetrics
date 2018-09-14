### Unit tests for CNVMetricsMethods.R functions

library(CNVMetrics)
library(rtracklayer)
library(GenomicRanges)

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



### Tests calculateWeigthedEucledianDistance() results

context("calculateWeigthedEucledianDistance() results")

test_that("calculateWeigthedEucledianDistance() must return good results", {
    
    dataGR <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
                                         end=c(29, 49, 55, 100, 110, 199, 220, 230)))
    
    values(dataGR) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
                                  File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
                                  File2 = c(NA, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
    
    expMatrix <- matrix(c(0, 0.387652570, 0.387652570, 0), byrow = T, 
                            ncol = 2, nrow = 2)
    colnames(expMatrix) <- c("File1", "File2")
    rownames(expMatrix) <- c("File1", "File2")
    
    expect_equivalent(calculateWeigthedEucledianDistance(dataGR), expMatrix, 
                        tolerance=0.000001)
    
})



