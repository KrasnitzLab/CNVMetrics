### Unit tests for CNVMetricsOverlapInternalMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)
### Tests calculateSzymkiewicz() results

context("calculateSzymkiewicz() results")

test_that("calculateSzymkiewicz() must return expected result of 1 for identical GRanges", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), strand = rep("+", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), strand = rep("-", 3))
    
    results <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- 1.0
    
    expect_equal(results, expected)
})

test_that("calculateSzymkiewicz() must return expected result of 1 for one GRanges included in the other", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(905048, 1554832, 31686841), 
                                          end = c(2204603, 4577608,41695808)), 
                        strand = rep("+", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    
    results <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- 1.0
    
    expect_equal(results, expected)
})

test_that("calculateSzymkiewicz() must return expected result of NA when first GRanges empty", {
    
    
    sample01 <- GRanges()
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    
    results <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- NA
    
    expect_equal(results, expected)
})


test_that("calculateSzymkiewicz() must return expected result of NA when second GRanges empty", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    sample02 <- GRanges()
    
    results <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- NA
    
    expect_equal(results, expected)
})

test_that("calculateSzymkiewicz() must return expected result of 0.5", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(100, 10, 1000), 
                                          end = c(199, 19, 1999)), 
                        strand = rep("-", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(150, 15, 1500), 
                                          end = c(249, 24, 2499)), 
                        strand = rep("+", 3))
    
    results <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- 0.5
    
    expect_equal(results, expected)
})


test_that("calculateSzymkiewicz() must return expected result of 0", {
    
    ### Create a Seqinfo Object
    chrInfo <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chr4", "chr5"),
                       seqlengths=rep(2500, 5), 
                       isCircular=rep(FALSE, 5),
                       genome="Alien")
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(100, 10, 1000), 
                                          end = c(199, 19, 1999)), 
                        strand = rep("-", 3), seqinfo = chrInfo)
    sample02 <- GRanges(seqnames = c("chr2", "chr3", "chr5"), 
                        ranges =  IRanges(start = c(150, 15, 1500), 
                                          end = c(249, 24, 2499)), 
                        strand = rep("+", 3), seqinfo = chrInfo)
    
    results <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- 0
    
    expect_equal(results, expected)
})

test_that("calculateSzymkiewicz() must return expected result of NA when two GRanges empty", {
    
    sample01 <- GRanges()
    sample02 <- GRanges()
    
    result <- CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
    
    expected <- NA
    
    expect_equal(result, expected)
})


### Tests calculateSorensen() results

context("calculateSorensen() results")

test_that("calculateSorensen() must return expected result of 1 for identical GRanges", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), strand = rep("+", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), strand = rep("-", 3))
    
    results <- CNVMetrics:::calculateSorensen(sample01, sample02)
    
    expected <- 1.0
    
    expect_equal(results, expected)
})

test_that("calculateSorensen() must return expected result for one GRanges included in the other", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(905048, 1554832, 31686841), 
                                          end = c(2204603, 4577608,41695808)), 
                        strand = rep("+", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    
    result <- CNVMetrics:::calculateSorensen(sample01, sample02)
    
    expected <- 0.018157313600969
    
    expect_equal(result, expected)
})

test_that("calculateSorensen() must return expected result when first GRanges empty", {
    
    sample01 <- GRanges()
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    
    result <- CNVMetrics:::calculateSorensen(sample01, sample02)
    
    expected <- 0.000
    
    expect_equal(result, expected)
})


test_that("calculateSorensen() must return expected result of NA when two GRanges empty", {
    
    ## Create two empty GRanges
    sample01 <- GRanges()
    sample02 <- GRanges()
    
    result <- CNVMetrics:::calculateSorensen(sample01, sample02)
    
    expected <- NA
    
    expect_equal(result, expected)
})

test_that("calculateSorensen() must return expected result of 0.5", {
    
    ## Create two GRanges
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(100, 10, 1000), 
                                          end = c(199, 19, 1999)), 
                        strand = rep("-", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(150, 15, 1500), 
                                          end = c(249, 24, 2499)), 
                        strand = rep("+", 3))
    
    result <- CNVMetrics:::calculateSorensen(sample01, sample02)
    
    expected <- 0.5
    
    expect_equal(result, expected)
})


test_that("calculateSorensen() must return expected result of 0", {
    
    ### Create a Seqinfo Object
    chrInfo <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chr4", "chr5"),
                       seqlengths=rep(2500, 5), 
                       isCircular=rep(FALSE, 5),
                       genome="Alien")
    
    ## Create two GRanges with no overlap
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(100, 10, 1000), 
                                          end = c(199, 19, 1999)), 
                        strand = rep("-", 3), seqinfo = chrInfo)
    sample02 <- GRanges(seqnames = c("chr2", "chr3", "chr5"), 
                        ranges =  IRanges(start = c(150, 15, 1500), 
                                          end = c(249, 24, 2499)), 
                        strand = rep("+", 3), seqinfo = chrInfo)
    
    result <- CNVMetrics:::calculateSorensen(sample01, sample02)
    
    expected <- 0
    
    expect_equal(result, expected)
})


### Tests calculateJaccard() results

context("calculateJaccard() results")

test_that("calculateJaccard() must return expected result of 1 for identical GRanges", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), strand = rep("+", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), strand = rep("-", 3))
    
    results <- CNVMetrics:::calculateJaccard(sample01, sample02)
    
    expected <- 1.0
    
    expect_equal(results, expected)
})

test_that("calculateJaccard() must return expected result for one GRanges included in the other", {
    
    
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(905048, 1554832, 31686841), 
                                          end = c(2204603, 4577608,41695808)), 
                        strand = rep("+", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    
    result <- CNVMetrics:::calculateJaccard(sample01, sample02)
    
    expected <- 0.009161833946548
    
    expect_equal(result, expected)
})

test_that("calculateJaccard() must return expected result when first GRanges empty", {
    
    sample01 <- GRanges()
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                          end = c(2004603, 4577608, 31695808)), 
                        strand = rep("-", 3))
    
    result <- CNVMetrics:::calculateJaccard(sample01, sample02)
    
    expected <- 0.000
    
    expect_equal(result, expected)
})


test_that("calculateJaccard() must return expected result of NA when two GRanges empty", {
    
    ## Create two empty GRanges
    sample01 <- GRanges()
    sample02 <- GRanges()
    
    result <- CNVMetrics:::calculateJaccard(sample01, sample02)
    
    expected <- NA
    
    expect_equal(result, expected)
})

test_that("calculateJaccard() must return expected result", {
    
    ## Create two GRanges
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(100, 10, 1000), 
                                          end = c(199, 19, 1999)), 
                        strand = rep("-", 3))
    sample02 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(150, 15, 1500), 
                                          end = c(249, 24, 2499)), 
                        strand = rep("+", 3))
    
    result <- CNVMetrics:::calculateJaccard(sample01, sample02)
    
    expected <- 0.333333333333333
    
    expect_equal(result, expected)
})


test_that("calculateJaccard() must return expected result of 0", {
    
    ### Create a Seqinfo Object
    chrInfo <- Seqinfo(seqnames=c("chr1", "chr2", "chr3", "chr4", "chr5"),
                       seqlengths=rep(2500, 5), 
                       isCircular=rep(FALSE, 5),
                       genome="Alien")
    
    ## Create two GRanges with no overlap
    sample01 <- GRanges(seqnames = c("chr1", "chr2", "chr4"), 
                        ranges =  IRanges(start = c(100, 10, 1000), 
                                          end = c(199, 19, 1999)), 
                        strand = rep("-", 3), seqinfo = chrInfo)
    sample02 <- GRanges(seqnames = c("chr2", "chr3", "chr5"), 
                        ranges =  IRanges(start = c(150, 15, 1500), 
                                          end = c(249, 24, 2499)), 
                        strand = rep("+", 3), seqinfo = chrInfo)
    
    result <- CNVMetrics:::calculateJaccard(sample01, sample02)
    
    expected <- 0
    
    expect_equal(result, expected)
})
