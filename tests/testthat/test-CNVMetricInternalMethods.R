### Unit tests for CNVMetricsInternalMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)



### Tests prepareInformation() results

context("createSegments() results")

test_that("createSegments() must return expected results", {
    
    segFiles <- list()
    segFiles[[1]] <-  GRanges(seqnames = "chr1", 
                              ranges = IRanges(start = c(1,200), 
                                               end=c(100, 300)))
    values(segFiles[[1]]) <- DataFrame(score = c(0.1, 0.5), 
                                       source = c("File1", "File1"))
    segFiles[[3]] <-  GRanges(seqnames = "chr1", 
                              ranges = IRanges(start = c(50,150), 
                                               end=c(110, 250)))
    values(segFiles[[3]]) <- DataFrame(score = c(0.3, 0.4), 
                                       source = c("File2", "File2"))
    
    sourceFiles <- list()
    sourceFiles[[1]] <- "File1"
    sourceFiles[[3]] <- "File2"
    
    expected <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                            end=c(49, 100, 110, 199, 250, 300)))
    
    values(expected) <- DataFrame(included = c(rep(TRUE, 6)), 
                                 File1 = c(0.1, 0.1, NA, NA, 0.5, 0.5), 
                                 File2 = c(NA, 0.3, 0.3, 0.4, 0.4, NA))
                     
    results <- CNVMetrics:::createSegments(fileList = segFiles, 
                                           sourceList = sourceFiles, 
                                           bedExclusion = NULL)
    
    expect_identical(results, expected)
})

test_that("createSegments() with BED file must return expected results", {
    
    segFiles <- list()
    segFiles[[1]] <-  GRanges(seqnames = "chr1", 
                              ranges = IRanges(start = c(1,200), 
                                               end=c(55, 220)))
    values(segFiles[[1]]) <- DataFrame(score = c(0.1, 0.5), 
                                       source = c("File1", "File1"))
    segFiles[[3]] <-  GRanges(seqnames = "chr1", 
                              ranges = IRanges(start = c(30,180), 
                                               end=c(110, 230)))
    values(segFiles[[3]]) <- DataFrame(score = c(0.3, 0.4), 
                                       source = c("File2", "File2"))
    
    bedInfo <- GRanges(seqnames = "chr1", 
            ranges = IRanges(start = c(50), end=c(100)))
    
    sourceFiles <- list()
    sourceFiles[[1]] <- "File1"
    sourceFiles[[3]] <- "File2"
    
    expected <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
                                         end=c(29, 49, 55, 100, 110, 199, 220, 230)))
    
    values(expected) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
                                  File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
                                  File2 = c(NA, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
    
    results <- CNVMetrics:::createSegments(fileList = segFiles, 
                                           sourceList = sourceFiles, 
                                           bedExclusion = bedInfo)
    
    expect_identical(results, expected)
})

### Tests calculateRegressedValues() results

context("calculateRegressedValues() results")

test_that("calculateRegressedValues() must return expected results", {
    
    segment  <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                         end=c(49, 100, 110, 199, 250, 300)))
    
    y <- c(0.1, 0.1, NA, NA, 0.15, 0.2)
    x <- c(NA, 0.3, 0.3, 0.4, 0.4, 0.8)
    
    elementMetadata(segment) <- DataFrame(included = c(rep(TRUE, 6)), 
                                  File1 = y,  
                                  File2 = x)
    
    segmentData <- list()
    segmentData$segments <- segment
    segmentData$regression <- list()
    segmentData$regression[[1]] <- list()
    segmentData$regression[[1]][["y_used"]] <- "File1"
    segmentData$regression[[1]][["x_used"]] <- "File2"
    segmentData$regression[[1]][["lm"]] <- lm("y ~ x", data.frame(x=x, y=y))
    
    results <- CNVMetrics:::calculateRegressedValues(segmentData)
    
    
    
    segmentExp  <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                         end=c(49, 100, 110, 199, 250, 300)))
    
    y <- c(0.1, 0.1, NA, NA, 0.15, 0.2)
    x <- c(NA, 0.11428571428571400459, 0.11428571428571400459, 
           0.13214285714285700646, 0.13214285714285700646, 
           0.20357142857142898618)
    
    values(segmentExp) <- DataFrame(included = c(rep(TRUE, 6)), 
                                 File1 = y,  
                                 File2 = x)
    
    expected <- segmentData
    expected$regressedData <- segmentData$segments 
    
    elementMetadata(expected$regressedData) <- DataFrame(included = c(rep(TRUE, 6)), 
                                                         File1 = y,  
                                                         File2 = x)
    
    expect_equal(results, expected)
})


### Tests doRegression() results

context("doRegression() results")

test_that("doRegression() must return expected results", {
    
    segmentData <- list()
    
    segment  <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                         end=c(49, 100, 110, 199, 250, 300)))
    
    y <- c(0.1, 0.1, NA, NA, 0.15, 0.2)
    x <- c(NA, 0.3, 0.3, 0.4, 0.4, 0.8)
    
    elementMetadata(segment) <- DataFrame(included = c(rep(TRUE, 6)), 
                                          File1 = y,  
                                          File2 = x)
    segmentData$segments <- segment
    
    results <- CNVMetrics:::doRegression(segmentData)
    
    
    
    segmentExp  <- GRanges(seqnames = "chr1", 
                           ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                            end=c(49, 100, 110, 199, 250, 300)))

    elementMetadata(segmentExp) <- DataFrame(included = c(rep(TRUE, 6)), 
                                          File1 = y,  
                                          File2 = x)
    
    expected <- list()
    expected$segments <- segmentExp
    
    subData <- data.frame(x=x, y=y)
    
    
    segmentData <- list()
    segmentData$segments <- segment
    expected$regression <- list()
    expected$regression[[1]] <- list()
    expected$regression[[1]][["lm"]] <- lm("y ~ x", subData)
    expected$regression[[1]][["y_used"]] <- "File1"
    expected$regression[[1]][["x_used"]] <- "File2"
    
    expect_equal(results, expected)
})


