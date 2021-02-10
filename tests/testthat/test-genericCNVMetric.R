### Unit tests for functions in genericCNVMetric.R file

library(CNVMetrics)
library(GenomicRanges)
library(S4Vectors)
library(IRanges)
library(GenomeInfoDb)

### Tests print() for CNVMetric class

context("print() CNVMetric class")

test_that("print() for CNVMetric object must return identical object", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1905048, 4554832, 31686841),
        end = c(2004603, 4577608, 31695808)), strand =  "*",
        state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))

    demo[["sample02"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1995066, 31611222, 31690000),
        end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
        state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))

    demo[["sample03"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1906069, 4558838),
        end = c(1909505, 4570601)), strand =  "*",
        state = c("AMPLIFICATION", "DELETION"))

    ## Calculating Sorensen metric
    expected <- calculateOverlapMetric(demo, method="sorensen")
    
    result <- print(expected)
    
    expect_equal(result, expected)
})


### Tests is() for CNVMetric class

context("is() CNVMetric class")

test_that("is() for CNVMetric object must return identical object", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841),
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                                                    end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
                                  state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    demo[["sample03"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(1906069, 4558838),
                                                    end = c(1909505, 4570601)), strand =  "*",
                                  state = c("AMPLIFICATION", "DELETION"))
    
    ## Calculating Sorensen metric
    metric <- calculateOverlapMetric(demo, method="sorensen")
    
    result <- is.CNVMetric(metric)
    
    expect_true(result)
})


