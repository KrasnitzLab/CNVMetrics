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


### Tests createDisjoinSegmentsForTwoSamples() results

context("createDisjoinSegmentsForTwoSamples() results")


test_that("createDisjoinSegmentsForTwoSamples() must return expected result when bedExclusion=NULL", {
    
    ## Granges with the information for the 2 samples
    gr1 <- GRanges(seqnames=c("chr1", "chr2", "chr4", "chr5"), 
                ranges=IRanges(start=c(1905048, 4554832, 31686841, 90), 
                                    end=c(2004603, 4577608, 31695808, 1001)), 
                strand=rep("+", 4),
                log2ratio=c(NA, 2.2323, NA, -1.1212))
    
    gr2 <- GRanges(seqnames=c("chr1", "chr2", "chr4", "chr5"), 
                   ranges=IRanges(start=c(2005048, 4564832, 31686841, 190), 
                                    end=c(23114603, 4567608, 31695808, 10001)), 
                strand=rep("+", 4),
                log2ratio=c(1.22, 2.2323, 3.33, -1.1212))    
    
    results <- CNVMetrics:::createDisjoinSegmentsForTwoSamples(segmentDataSample1=gr1,
                        segmentDataSample2=gr2, bedExclusion=NULL)
    
    expected <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr2", "chr2",
                            "chr4", "chr5", "chr5", "chr5"), 
                    ranges=IRanges(start=c(1905048, 2005048, 4554832, 4564832,
                                    4567609, 31686841, 90, 190, 1002), 
                                end=c(2004603, 23114603, 4564831, 4567608,
                                    4577608, 31695808, 189, 1001, 10001)), 
                    strand=rep("+", 9), included=rep(TRUE, 9),
                    sample_1=c(NA, NA, 2.2323, 2.2323, 2.2323, NA, -1.1212, 
                                -1.1212, NA),
                    sample_2=c(NA, 1.2200, NA, 2.2323, NA, 3.3300, NA, -1.1212,
                                -1.1212))
    
    expect_equal(results, expected)
})



test_that("createDisjoinSegmentsForTwoSamples() must return expected result for specific bedExclusion data", {
    
    ## Granges with the information for the 2 samples
    gr1 <- GRanges(seqnames=c("chr1", "chr2", "chr4", "chr5"), 
                    ranges=IRanges(start=c(190508, 4554832, 31686841, 900), 
                                end=c(2004603, 4577608, 31695808, 1001)), 
                    strand=rep("+", 4),
                    log2ratio=c(NA, 2.2323, NA, -1.1212))
    
    gr2 <- GRanges(seqnames=c("chr1", "chr2", "chr4", "chr5"), 
                    ranges=IRanges(start=c(2005048, 4564832, 31686841, 190), 
                                end=c(23114603, 4567608, 31695808, 10001)), 
                    strand=rep("+", 4),
                    log2ratio=c(1.22, 2.2323, 3.33, -1.1212))    
    
    bedExclusion <-  GRanges(seqnames=c("chr2",  "chr5"), 
                             ranges=IRanges(start=c(4566832, 90), 
                                            end=c(4568608, 801)), 
                             strand=rep("+", 2))  
    
    results <- CNVMetrics:::createDisjoinSegmentsForTwoSamples(segmentDataSample1=gr1,
                    segmentDataSample2=gr2, bedExclusion=bedExclusion)
    
    expected <- GRanges(seqnames=c("chr1", "chr1", "chr2", "chr2", "chr2",
                                        "chr4", "chr5", "chr5", "chr5"), 
                    ranges=IRanges(start=c(190508, 2005048, 4554832, 4564832,
                                        4567609, 31686841, 190, 900, 1002), 
                                    end=c(2004603, 23114603, 4564831, 4567608,
                                        4577608, 31695808, 899, 1001, 10001)), 
                                    strand=rep("+", 9), 
                    included=c(rep(TRUE, 3), rep(FALSE, 2), TRUE, FALSE, 
                                    rep(TRUE, 2)),
                    sample_1=c(NA, NA, 2.2323, 2.2323, 2.2323, NA, NA, 
                                    -1.1212, NA),
                    sample_2=c(NA, 1.2200, NA, 2.2323, NA, 3.3300, -1.1212, 
                                    -1.1212, -1.1212))

    expect_equal(results, expected)
})








