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
                CN=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim = c(1, 3, 4)), error_message)
})

test_that("processSim() must return an error when nbSim is zero", {

    sample <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                CN=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim = 0), error_message)
})

test_that("processSim() must return an error when nbSim is negative integer", {

    sample <- GRanges(seqnames = "chr1",
                ranges=IRanges(start = c(1995066, 31611222, 31690000),
                end=c(2204505, 31689898, 31895666)), strand=c("-", "+", "+"),
                state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                CN=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim = -21), error_message)
})

test_that("processSim() must return an error when nbSim is a character string", {

    sample <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                end = c(2204505, 31689898, 31895666)),
                strand =  c("-", "+", "+"),
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"),
                CN=c(0.5849625, 0.4444333, -1))

    error_message <- "nbSim must be a positive integer"

    expect_error(processSim(curSample = sample, nbSim ="hello"), error_message)
})



test_that("processSim() must return an error when curSample is a character string", {

    error_message <- "the \'curSample\' argument must be a \'GRanges\' object"

    expect_error(processSim(curSample = "sample", nbSim = 10), error_message)
})
