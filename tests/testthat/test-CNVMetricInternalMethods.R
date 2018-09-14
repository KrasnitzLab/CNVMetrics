### Unit tests for CNVMetricsInternalMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)



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
