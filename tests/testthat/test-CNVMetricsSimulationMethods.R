### Unit tests for CNVMetricsMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)



#############################################################################
### Tests processSim() results
#############################################################################

context("processSim() results")

test_that("processSim() must return an error when nbSim is a vector of integers", {

    sample <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                log2ratio=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim = c(1, 3, 4)), error_message)
})

test_that("processSim() must return an error when nbSim is zero", {

    sample <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                log2ratio=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim = 0), error_message)
})

test_that("processSim() must return an error when nbSim is negative integer", {

    sample <- GRanges(seqnames = "chr1",
                ranges=IRanges(start = c(1995066, 31611222, 31690000),
                end=c(2204505, 31689898, 31895666)), strand=c("-", "+", "+"),
                state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                log2ratio=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim = -21), error_message)
})

test_that("processSim() must return an error when nbSim is a character string", {

    sample <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                end = c(2204505, 31689898, 31895666)),
                strand =  c("-", "+", "+"),
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                log2ratio=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim ="hello"), error_message)
})


test_that("processSim() must return an error when curSample doesn't have a log2ratio column", {

    sample <- GRanges(seqnames="chr1",
                      ranges=IRanges(start = c(1995066, 31611222, 31690000),
                                        end = c(2204505, 31689898, 31895666)),
                      strand=c("-", "+", "+"),
                      state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                      CN=c(0.5849625, 0.4444333, -1))

    error_message <- paste0("the sample must have a metadata column ",
                            "called \'log2ratio\'")

    expect_error(processSim(curSample=sample, nbSim=3), error_message)
})


test_that("processSim() must return an error when curSample is a character string", {

    error_message <- "the \'curSample\' argument must be a \'GRanges\' object"

    expect_error(processSim(curSample = "sample", nbSim = 10), error_message)
})


test_that("processSim() must return the expected result", {

    sample <- GRanges(seqnames = "chr1",
                      ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                                        end = c(2204505, 31689898, 31895666)),
                      strand =  c("-", "+", "+"),
                      state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                      log2ratio=c(0.5849625, 0.4444333, -1))

    expected <- data.frame(ID=c(rep("S1", 5), rep("S2", 5), rep("S3", 5)),
                    chr=c(rep("chr1", 15)), start=c(1995066, 31611222, 31612582,
                    31690000,  31698335, 1995066, 2204346, 31611222, 31690000,
                    31816977, 1995066, 2187998, 31611222, 31690000, 31814749),
                    end=c(2204505, 31612581, 31689898, 31698334, 31895669,
                          2204345,  2204505, 31689898, 31816976, 31895669,
                          2187997,  2204505, 31689898, 31814748, 31895669),
                    log2ratio=c(-1.000000,-1.000000, 0.4444333, 0.4444333, 0.5849625,
                                -1.000000, 0.5849625  , 0.5849625  , 0.5849625  , 0.4444333,
                                0.5849625  , -1.000000, -1.000000, -1.000000, 0.4444333),
                    state=c("DELETION", "DELETION", "AMPLIFICATION", "AMPLIFICATION",
                            "AMPLIFICATION", "DELETION", "AMPLIFICATION", "AMPLIFICATION",
                            "AMPLIFICATION", "AMPLIFICATION", "AMPLIFICATION", "DELETION",
                            "DELETION", "DELETION", "AMPLIFICATION"))

    set.seed(1212)

    result <- processSim(curSample = sample, nbSim = 3)

    expect_equal(result, expected)
})
