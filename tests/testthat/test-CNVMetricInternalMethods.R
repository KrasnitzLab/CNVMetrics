### Unit tests for CNVMetricsInternalMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)


### Tests createDisjoinSegmentsForTwoSamples() results

context("createDisjoinSegmentsForTwoSamples() results")

test_that("createDisjoinSegmentsForTwoSamples() must return expected results", {

    
    segment1 <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                         end=c(49, 100, 110, 199, 250, 300)),
                        log2ratio=c(0.1222, 1.3211, 0.1212, -1.1111, -2.2222, -0.4444))
    
    
    segment2  <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 251),
                                         end=c(300, 100,  199,  400)),
                        log2ratio=c(0.002, 2.3111, 1.2222, -0.4444))
    
    results <- CNVMetrics:::createDisjoinSegmentsForTwoSamples(segmentDataSample1=segment1, 
                        segmentDataSample2=segment2, bedExclusion = NULL)
    
   
    expected <- GRanges(seqnames="chr1", 
                    ranges=IRanges(start=c(1, 50, 101, 111, 150, 200, 251, 301),
                                   end=c(49, 100, 110, 149, 199, 250, 300, 400)),
                    included=rep(TRUE, 8),
                    sample_1=c(0.1222, 1.3211, 0.1212, NA, -1.1111, -2.2222, -0.4444, NA),
                    sample_2=c(0.0020, 2.3111, 1.2222, 1.2222, 1.2222, 0.0020, -0.4444, -0.4444))
    
    
    expect_equal(results, expected)
})


test_that("createDisjoinSegmentsForTwoSamples() must return expected results when exclusion file used", {
    
    
    segment1 <- GRanges(seqnames = "chr1", 
                        ranges = IRanges(start=c(1, 50, 101, 150, 200, 251),
                                         end=c(49, 100, 110, 199, 250, 300)),
                        log2ratio=c(0.1222, 1.3211, 0.1212, -1.1111, -2.2222, -0.4444))
    
    
    segment2  <- GRanges(seqnames = "chr1", 
                         ranges = IRanges(start=c(1, 250, 301, 551),
                                          end=c(200, 300,  399,  800)),
                         log2ratio=c(0.002, 2.3111, 1.2222, -0.4444))
    
    exclusion <- GRanges(seqnames="chr1", ranges=IRanges(start=c(60, 444), end=c(102, 900)))
    
    results <- CNVMetrics:::createDisjoinSegmentsForTwoSamples(segmentDataSample1=segment1, 
                        segmentDataSample2=segment2, bedExclusion = exclusion)
    
    
    expected <- GRanges(seqnames="chr1", 
                        ranges=IRanges(start=c(1, 50, 101, 111, 150, 200, 201, 250, 251, 301, 551),
                                       end=c(49, 100, 110, 149, 199, 200, 249, 250, 300, 399, 800)),
                        included=c(TRUE, FALSE, FALSE, rep(TRUE, 7), FALSE),
                        sample_1=c(0.1222, 1.3211, 0.1212, NA, -1.1111, -2.2222, -2.2222, -2.2222, -0.4444, NA, NA),
                        sample_2=c(0.0020, 0.0020, 0.0020, 0.0020, 0.0020, 0.0020, NA, 2.3111,  2.3111, 1.2222, -0.4444))
    
    
    expect_equal(results, expected)
})


### Tests validateCalculateOverlapMetricParameters() results

context("validateCalculateOverlapMetricParameters() results")

test_that("validateCalculateOverlapMetricParameters() must return expected zero", {
    
    results <- CNVMetrics:::validateCalculateOverlapMetricParameters(states="LOH", 
                                        nJobs=1)
    
    expect_equal(results, 0L)
})

test_that("validateCalculateOverlapMetricParameters() must return error when nJobs is negative", {
    
    error_message <- "nJobs must be a positive integer"
    
    expect_error(CNVMetrics:::validateCalculateOverlapMetricParameters(states="LOH", 
            nJobs=-3), error_message)
})

test_that("validateCalculateOverlapMetricParameters() must return error when states is a list", {
    
    error_message <- paste0("the \'states\' argument must be a vector ", 
                                "of strings with at least one value")
    
    expect_error(CNVMetrics:::validateCalculateOverlapMetricParameters(states=list(A=1, B=2), 
            nJobs=1), error_message)
})


test_that("validateCalculateOverlapMetricParameters() must return error when states is a vector of numerics", {
    
    error_message <- paste0("the \'states\' argument must be a vector ", 
                                "of strings with at least one value")
    
    expect_error(CNVMetrics:::validateCalculateOverlapMetricParameters(states=c(1,3,4),
            nJobs=1), error_message)
})

