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



### Tests calculateWeightedEuclideanDistance() results

context("calculateWeightedEuclideanDistance() results")

test_that("calculateWeightedEuclideanDistance() must return an error when segmentData is not of good class", {
    
    error_message <- "segmentData must be a list marked as preMetricSegments class."
    
    dataGR <- list()
    
    dataGR$segments <- GRanges(seqnames = "chr1", 
                               ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
                                                end=c(29, 49, 55, 100, 110, 199, 220, 230)))
    
    values(dataGR$segments) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
                                File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
                                File2 = c(NA, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
    
    expect_error(calculateWeightedEuclideanDistance(dataGR), error_message) 
})


test_that("calculateWeightedEuclideanDistance() must return good results 01", {
    
    dataGR <- list()
    class(dataGR) <- "preMetricSegments"
    
    dataGR$segments <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
                                         end=c(29, 49, 55, 100, 110, 199, 220, 230)))
    
    values(dataGR$segments) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
                                  File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
                                  File2 = c(NA, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
    
    expMatrix <- matrix(c(0, 0.387652570, 0.387652570, 0), byrow = T, 
                            ncol = 2, nrow = 2)
    colnames(expMatrix) <- c("File1", "File2")
    rownames(expMatrix) <- c("File1", "File2")
    
    expect_equivalent(calculateWeightedEuclideanDistance(dataGR, minThreshold = 0.04), expMatrix, 
                        tolerance=0.000001)
})


test_that("calculateWeightedEuclideanDistance() must return good results 02", {
    
    dataGR <- list()
    class(dataGR) <- "preMetricSegments"
    
    dataGR$segments <- GRanges(seqnames = "chr1", 
                               ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
                                                end=c(29, 49, 55, 100, 110, 199, 220, 230)))
    
    values(dataGR$segments) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
                                         File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
                                         File2 = c(NA, 0.15, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
    
    expMatrix <- matrix(c(0, 0.174485599340559, 0.174485599340559, 0), byrow = T, 
                        ncol = 2, nrow = 2)
    colnames(expMatrix) <- c("File1", "File2")
    rownames(expMatrix) <- c("File1", "File2")
    
    expect_equivalent(calculateWeightedEuclideanDistance(dataGR, minThreshold = 0.07), expMatrix, 
                      tolerance=0.000001)
})


### Tests calculateOverlapRegionsMetric() results

context("calculateOverlapRegionsMetric() results")

test_that("calculateOverlapRegionsMetric() must return an error when segmentData has only one sample", {
    
    error_message <- "at least 2 samples must be present in the segmentData"
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1", 
        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
        end = c(2004603, 4577608, 31695808)), strand =  "*",
        state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    expect_error(calculateOverlapRegionsMetric(demo), error_message) 
})

test_that("calculateOverlapRegionsMetric() must return an error when segmentData has metadata status instead of state", {
    
    error_message <- paste0("at least one sample doesn't have a metadata column ", "
             called \'state\'")
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  status = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    demo[["sample02"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  status= c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    expect_error(calculateOverlapRegionsMetric(demo), error_message) 
})

test_that("calculateOverlapRegionsMetric() must return an error when segmentData doesn't have metadata state", {
    
    error_message <- paste0("at least one sample doesn't have a metadata column ", "
             called \'state\'")
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), 
                                  strand =  "*")
    
    demo[["sample02"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), 
                                  strand =  "*")
    
    expect_error(calculateOverlapRegionsMetric(demo), error_message) 
})
