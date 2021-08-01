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


test_that("print() for CNVMetric object must return expected result when 6 samples", {
    
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
    
    demo[["sample04"]] <- GRanges(seqnames = "chr1",
            ranges =  IRanges(start = c(1906039, 4558538),
                    end = c(1909705, 4570801)), strand =  "*",
            state = c("AMPLIFICATION", "DELETION"))
    
    demo[["sample05"]] <- GRanges(seqnames = "chr1",
            ranges =  IRanges(start = c(1926039, 4557538),
                    end = c(2909705, 4570901)), strand =  "*",
            state = c("AMPLIFICATION", "DELETION"))
    
    
    demo[["sample06"]] <- GRanges(seqnames = "chr1",
            ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                    end = c(2204525, 31689998, 31895766)), strand =  c("-", "+", "+"),
            state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    ## Calculating Sorensen metric
    expected <- calculateOverlapMetric(demo, method="sorensen")
    
    result <- print(expected)
    
    expect_equal(result, expected)
})


test_that("print() for CNVMetric object must return expected result when 7 samples", {
    
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
    
    demo[["sample04"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1906039, 4558538),
            end = c(1909705, 4570801)), strand =  "*",
        state = c("AMPLIFICATION", "DELETION"))
    
    demo[["sample05"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1926039, 4557538),
            end = c(2909705, 4570901)), strand =  "*",
        state = c("AMPLIFICATION", "DELETION"))
    
    demo[["sample06"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1995066, 31611222, 31690000),
            end = c(2204525, 31689998, 31895766)), strand =  c("-", "+", "+"),
        state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    demo[["sample07"]] <- GRanges(seqnames = "chr1",
        ranges =  IRanges(start = c(1995096, 3161142, 31680000),
            end = c(2204510, 31689993, 31865766)), strand =  c("-", "+", "+"),
        state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    
    ## Calculating Szymkiewicz metric
    expected <- calculateOverlapMetric(demo, method="szymkiewicz")
    
    result <- print(expected)
    
    expect_equal(result, expected)
})


test_that("print() for CNVMetric object must return expected result with weighted Euclidean Distance metric", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(1905048, 4554832, 31686841),
                                                    end = c(2004603, 4577608, 31695808)), strand =  "*",
                                  log2ratio = c(0.3211, 0.4322, -0.9292))
    
    demo[["sample02"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(1995066, 31611222, 31690000),
                                                    end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
                                  log2ratio = c(2.2323, 1.4343, -1.1313))
    
    demo[["sample03"]] <- GRanges(seqnames = "chr1",
                                  ranges =  IRanges(start = c(1906069, 4558838),
                                                    end = c(1909505, 4570601)), strand =  "*",
                                  log2ratio = c(2.1212, -0.9898))
    
    ## Calculating Weighted Euclidean Distance metric
    expected <- calculateLog2ratioMetric(demo, method="weightedEuclideanDistance", 
                                            minThreshold=0.2, excludedRegions=NULL)
    
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


