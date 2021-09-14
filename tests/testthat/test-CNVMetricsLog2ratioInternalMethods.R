### Unit tests for CNVMetricsLog2ratioInternalMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)

### Tests calculateWeightedEuclideanDistanceFor2Samples() results

context("calculateWeightedEuclideanDistanceFor2Samples() results")


test_that("calculateWeightedEuclideanDistanceFor2Samples() must return expected result for 2 identical samples", {
    
    ## Granges with the information for the 2 samples
    gr <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                    end = c(2004603, 4577608, 31695808)), 
                  strand = rep("+", 3),
                  included = rep(TRUE, 3),
                  sample_1 = c(1.1212, 2.2323, -0.9999),
                  sample_2 = c(1.1212, 2.2323, -0.9999))
    
    
    results <- CNVMetrics:::calculateWeightedEuclideanDistanceFor2Samples(segmentData = gr,
                                            minThreshold = 0.2)
    
    expected <- 1.0
    
    expect_equal(results, expected)
})

test_that("calculateWeightedEuclideanDistanceFor2Samples() must return expected result for 2 identical samples when using minThreshold", {
    
    ## Granges with the information for the 2 samples
    gr <- GRanges(seqnames = c("chr1", "chr2", "chr4", "chr5"), 
                  ranges =  IRanges(start = c(1905048, 4554832, 31686841, 90), 
                                    end = c(2004603, 4577608, 31695808, 1001)), 
                  strand = rep("+", 4),
                  included = rep(TRUE, 4),
                  sample_1 = c(1.1212, 2.2323, -0.8999, -1.1212),
                  sample_2 = c(0.9512, 2.2443, -0.9999, -1.3211))
    
    
    results <- CNVMetrics:::calculateWeightedEuclideanDistanceFor2Samples(segmentData = gr,
                                        minThreshold = 0.2)
    
    expected <- 1.0
    
    expect_equal(results, expected)
})

test_that("calculateWeightedEuclideanDistanceFor2Samples() must return expected result for 2 different samples", {
    
    ## Granges with the information for the 2 samples
    gr <- GRanges(seqnames = c("chr1", "chr2", "chr4", "chr5"), 
                  ranges =  IRanges(start = c(1905048, 4554832, 31686841, 90), 
                                    end = c(2004603, 4577608, 31695808, 1001)), 
                  strand = rep("+", 4),
                  included = rep(TRUE, 4),
                  sample_1 = c(2.1212, 2.2323, -0.8999, -1.1212),
                  sample_2 = c(0.0512, 4.2443, -1.9999, -2.3211))
    
    
    results <- CNVMetrics:::calculateWeightedEuclideanDistanceFor2Samples(segmentData = gr,
                                minThreshold = 0.2)
    
    expected <- 0.0733102557130918
    
    expect_equal(results, expected)
})



test_that("calculateWeightedEuclideanDistanceFor2Samples() must return expected result for 2 different samples and some regions not included", {
    
    ## Granges with the information for the 2 samples
    gr <- GRanges(seqnames = c("chr1", "chr2", "chr4", "chr5"), 
                  ranges =  IRanges(start = c(1905048, 4554832, 31686841, 90), 
                                    end = c(2004603, 4577608, 31695808, 1001)), 
                  strand = rep("+", 4),
                  included = c(rep(TRUE, 3), FALSE),
                  sample_1 = c(2.1212, 2.2323, -0.8999, -1.1212),
                  sample_2 = c(0.0512, 4.2443, -1.9999, -2.3211))
    
    
    results <- CNVMetrics:::calculateWeightedEuclideanDistanceFor2Samples(segmentData = gr,
                                        minThreshold = 0.2)
    
    expected <- 0.0765246076087433
    
    expect_equal(results, expected)
})


test_that("calculateWeightedEuclideanDistanceFor2Samples() must return expected result when not calculation possible", {
    
    ## Granges with the information for the 2 samples
    gr <- GRanges(seqnames = c("chr1", "chr2", "chr4", "chr5"), 
                  ranges =  IRanges(start = c(1905048, 4554832, 31686841, 90), 
                                    end = c(2004603, 4577608, 31695808, 1001)), 
                  strand = rep("+", 4),
                  included = c(rep(TRUE, 3), FALSE),
                  sample_1 = c(NA, 2.2323, NA, -1.1212),
                  sample_2 = c(0.0512, NA, -1.9999, NA))
    
    
    results <- CNVMetrics:::calculateWeightedEuclideanDistanceFor2Samples(segmentData = gr,
                                                    minThreshold = 0.2)
    
    expected <- NA
    
    expect_equal(results, expected)
})

