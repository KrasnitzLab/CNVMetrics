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





#############################################################################
### Tests simChr() results
#############################################################################

context("simChr() results")


test_that("simChr() must return the expected result when only 1 segment", {

    sample <- GRanges(seqnames = "chr1",
                      ranges =  IRanges(start = c(1995066),
                                        end = c(2204505)),
                      strand =  c("-"),
                      state = c("AMPLIFICATION"),
                      log2ratio=c(0.5849625))

    expected <- list()

    expected[[1]] <- data.frame(ID=c("S1"),
                           chr=c("chr1"), start=c(0),
                           end=c(1),
                           log2ratio=c(0.5849625),
                           state=c("AMPLIFICATION"))

    expected[[2]] <- data.frame(ID=c("S2"),
                                chr=c("chr1"), start=c(0),
                                end=c(1),
                                log2ratio=c(0.5849625),
                                state=c("AMPLIFICATION"))

    expected[[3]] <- data.frame(ID=c("S3"),
                                chr=c("chr1"), start=c(0),
                                end=c(1),
                                log2ratio=c(0.5849625),
                                state=c("AMPLIFICATION"))

    set.seed(1212)

    result <- CNVMetrics:::simChr(curSample = sample, chrCur="chr1", nbSim = 3)

    expect_equal(result, expected)
})

test_that("simChr() must return the expected result when only 2 segments", {

    sample <- GRanges(seqnames = rep("chr1", 2),
                      ranges =  IRanges(start = c(1995066, 4333332),
                                        end = c(2204505, 14333332)),
                      strand =  c("-", "*"),
                      state = c("AMPLIFICATION", "NEUTRAL"),
                      log2ratio=c(0.5849625, 0.00011))

    expected <- list()

    expected[[1]] <- data.frame(ID=rep("S1", 2),
                                chr=rep("chr1", 2), start=c(0, 0.019211531757713),
                                end=c(0.019211531757713, 1),
                                log2ratio=c(0.5849625, 0.00011),
                                state=c("AMPLIFICATION", "NEUTRAL"))

    expected[[2]] <- data.frame(ID=rep("S2", 2),
                                chr=rep("chr1", 2), start=c(0, 0.971207630270844),
                                end=c(0.971207630270844, 1),
                                log2ratio=c(0.00011, 0.5849625),
                                state=c("NEUTRAL", "AMPLIFICATION"))

    expected[[3]] <- data.frame(ID=rep("S3", 2),
                                chr=rep("chr1", 2), start=c(0, 0.988944448574609),
                                end=c(0.988944448574609, 1),
                                log2ratio=c(0.00011, 0.5849625),
                                state=c("NEUTRAL", "AMPLIFICATION"))

    set.seed(112)

    result <- CNVMetrics:::simChr(curSample = sample, chrCur="chr1", nbSim = 3)

    expect_equal(result, expected)
})


