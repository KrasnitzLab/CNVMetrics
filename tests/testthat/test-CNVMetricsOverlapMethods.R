### Unit tests for CNVMetricsOverlapMethods.R functions

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)


### Tests calculateOverlapMetric() results

context("calculateOverlapMetric() results")

test_that("calculateOverlapMetric() must return error when only one sample present", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1905048, 4554832, 31686841),
        end = c(2004603, 4577608, 31695808)), strand =  "*",
        state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))


    error_message <- "at least 2 samples must be present in the segmentData"
    
    expect_error(calculateOverlapMetric(segmentData = demo, 
                                                method = "sorensen"), 
                 error_message)
})

test_that("calculateOverlapMetric() must return error when method is available", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1905048, 4554832, 31686841),
                end = c(2004603, 4577608, 31695808)), strand =  "*",
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(1905048, 4554832, 31686841),
                end = c(2004603, 4577608, 31695808)), strand =  "*",
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    error_message <- gettextf("'arg' should be one of %s", 
             paste(dQuote(c("sorensen", "szymkiewicz")), 
                                                collapse = ", "))

    expect_error(calculateOverlapMetric(segmentData = demo, 
                                               method = "typo"), 
                 error_message)
})


test_that("calculateOverlapMetric() must return 1 when two samples identical with sorensen", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                    ranges =  IRanges(start = c(1905048, 4554832, 31686841),
                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                    state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                    ranges =  IRanges(start = c(1905048, 4554832, 31686841),
                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                    state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    result <- calculateOverlapMetric(segmentData=demo, 
                                                method="sorensen")
    
    expected <- list()
    expected$AMPLIFICATION <- matrix(c(NA, NA, 1, NA), nrow=2, ncol=2, 
                                    byrow=TRUE)
    colnames(expected$AMPLIFICATION) <- c("sample01", "sample02")
    rownames(expected$AMPLIFICATION) <- c("sample01", "sample02")
    
    expected$DELETION <- matrix(c(NA, NA, 1, NA), nrow=2, ncol=2, 
                                   byrow=TRUE)
    colnames(expected$DELETION) <- c("sample01", "sample02")
    rownames(expected$DELETION) <- c("sample01", "sample02")
    
    class(expected) <- "CNVMetric"
    attr(expected, 'metric') <- "sorensen"
    
    expect_equal(result, expected)
})


test_that("calculateOverlapMetric() must return expected results with sorensen", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(100, 300, 800),
                end = c(200, 500, 900)), strand =  "*",
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(150, 600, 1000),
                end = c(250, 700, 1500)), strand =  "*",
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample03"]] <- GRanges(seqnames = "chr1",
                ranges =  IRanges(start = c(50, 600, 1000),
                end = c(250, 700, 2000)), strand =  "*",
                state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    
    result <- calculateOverlapMetric(segmentData=demo, 
                                            method="sorensen")
    
    expected <- list()
    expected$AMPLIFICATION <- matrix(c(NA, NA, NA, 0.202380952380952, NA, NA, 
                0.334437086092715, 0.801587301587302, NA), nrow=3, ncol=3, 
                                    byrow=TRUE)
    colnames(expected$AMPLIFICATION) <- c("sample01", "sample02", "sample03")
    rownames(expected$AMPLIFICATION) <- c("sample01", "sample02", "sample03")
    
    expected$DELETION <- matrix(c(NA, NA, NA, 0, NA, NA, 
                                  0, 0.667110519307590, NA), 
                                nrow=3, ncol=3, 
                                byrow=TRUE)
    colnames(expected$DELETION) <- c("sample01", "sample02", "sample03")
    rownames(expected$DELETION) <- c("sample01", "sample02", "sample03")
    
    class(expected) <- "CNVMetric"
    attr(expected, 'metric') <- "sorensen"
    
    expect_equal(result, expected)
})


test_that("calculateOverlapMetric() must return expected results with szymkiewicz", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                    ranges =  IRanges(start = c(100, 300, 800),
                    end = c(200, 500, 900)), strand =  "*",
                    state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                    ranges =  IRanges(start = c(150, 600, 1000),
                    end = c(250, 700, 1500)), strand =  "*",
                    state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample03"]] <- GRanges(seqnames = "chr1",
                    ranges =  IRanges(start = c(50, 600, 1000),
                    end = c(250, 700, 2000)), strand =  "*",
                    state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    result <- calculateOverlapMetric(segmentData=demo, 
                                            metho="szymkiewicz")
    
    expected <- list()
    expected$AMPLIFICATION <- matrix(c(NA, NA, NA, 0.252475247524752, NA, NA, 
                                0.334437086092715, 1.000000000000000, NA), 
                                nrow = 3, ncol = 3, byrow = TRUE)
    colnames(expected$AMPLIFICATION) <- c("sample01", "sample02", "sample03")
    rownames(expected$AMPLIFICATION) <- c("sample01", "sample02", "sample03")
    
    expected$DELETION <- matrix(c(NA, NA, NA, 0.000000000000000, NA, NA, 
                                0.000000000000000, 1.000000000000000, NA), 
                                nrow = 3, ncol = 3, 
                                byrow = TRUE)
    colnames(expected$DELETION) <- c("sample01", "sample02", "sample03")
    rownames(expected$DELETION) <- c("sample01", "sample02", "sample03")
    
    class(expected) <- "CNVMetric"
    attr(expected, 'metric') <- "szymkiewicz"
    
    expect_equal(result, expected)
})


test_that("calculateOverlapMetric() must return expected results with jaccard", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(100, 300, 800),
                                                    end = c(200, 500, 900)), strand =  "*",
                                  state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(150, 600, 1000),
                                                    end = c(250, 700, 1500)), strand =  "*",
                                  state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample03"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(50, 600, 1000),
                                                    end = c(250, 700, 2000)), strand =  "*",
                                  state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    result <- calculateOverlapMetric(segmentData=demo, 
                                     metho="jaccard")
    
    expected <- list()
    expected$AMPLIFICATION <- matrix(c(NA, NA, NA, 0.112582781456954, NA, NA, 
                                       0.200795228628231, 0.668874172185430, NA), 
                                     nrow = 3, ncol = 3, byrow = TRUE)
    colnames(expected$AMPLIFICATION) <- c("sample01", "sample02", "sample03")
    rownames(expected$AMPLIFICATION) <- c("sample01", "sample02", "sample03")
    
    expected$DELETION <- matrix(c(NA, NA, NA, 0.000000000000000, NA, NA, 
                                  0.000000000000000, 0.500499500499501, NA), 
                                nrow = 3, ncol = 3, 
                                byrow = TRUE)
    colnames(expected$DELETION) <- c("sample01", "sample02", "sample03")
    rownames(expected$DELETION) <- c("sample01", "sample02", "sample03")
    
    class(expected) <- "CNVMetric"
    attr(expected, 'metric') <- "jaccard"
    
    expect_equal(result, expected)
})


test_that("calculateOverlapMetric() must return an error when segmentData has only one sample", {
    
    error_message <- "at least 2 samples must be present in the segmentData"
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    expect_error(calculateOverlapMetric(demo), error_message) 
})

test_that("calculateOverlapMetric() must return an error when segmentData has metadata status instead of state", {
    
    error_message <- paste0("at least one sample doesn't have a metadata column ", 
                            "called \'state\'")
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  status = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    demo[["sample02"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  status= c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    expect_error(calculateOverlapMetric(demo), error_message) 
})

test_that("calculateOverlapMetric() must return an error when segmentData doesn't have metadata state", {
    
    error_message <- paste0("at least one sample doesn't have a metadata column ", 
                            "called \'state\'")
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), 
                                  strand =  "*")
    
    demo[["sample02"]] <- GRanges(seqnames = "chr1", 
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
                                                    end = c(2004603, 4577608, 31695808)), 
                                  strand =  "*")
    
    expect_error(calculateOverlapMetric(demo), error_message) 
})

