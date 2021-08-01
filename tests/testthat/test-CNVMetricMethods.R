### Unit tests for CNVMetricsMethods.R functions

library(CNVMetrics)
library(rtracklayer)
library(GenomicRanges)

chrInfo <- Seqinfo(seqnames=c("chr1", "chr2", "chr3"),
                   seqlengths=c(10000, 20000, 1500), 
                   isCircular=c(FALSE, FALSE, FALSE),
                   genome="Alien")

# ### Tests prepareInformation() results
# 
# context("prepareInformation() results")
# 
# test_that("prepareInformation() must return error when segDirectory is not a Seqinfo", {
#     
#     error_message <- "chrInfo must be a Seqinfo object."
#     
#     expect_error(prepareInformation(segDirectory = "test", chrInfo = 33), 
#                  error_message)
#     expect_error(prepareInformation(segDirectory = "test", chrInfo = "allo"), 
#                  error_message)
# })
# 
# test_that("prepareInformation() must return error when segmentWithHeader is not a logical", {
#     
#     error_message02 <- "segmentWithHeader must be a logical."
#     
#     expect_error(prepareInformation(segDirectory = "test", chrInfo = chrInfo, 
#                                     segmentWithHeader = "hi"), error_message02)
#     expect_error(prepareInformation(segDirectory = "test", chrInfo = chrInfo, 
#                                     segmentWithHeader = 444), error_message02)
# })


# 
# ### Tests calculateWeightedEuclideanDistance() results
# 
# context("calculateWeightedEuclideanDistance() results")
# 
# test_that("calculateWeightedEuclideanDistance() must return an error when segmentData is not of good class", {
#     
#     error_message <- "segmentData must be a list marked as preMetricSegments class."
#     
#     dataGR <- list()
#     
#     dataGR$segments <- GRanges(seqnames = "chr1", 
#                                ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
#                                                 end=c(29, 49, 55, 100, 110, 199, 220, 230)))
#     
#     values(dataGR$segments) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
#                                 File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
#                                 File2 = c(NA, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
#     
#     expect_error(calculateWeightedEuclideanDistance(dataGR), error_message) 
# })
# 
# 
# test_that("calculateWeightedEuclideanDistance() must return good results 01", {
#     
#     dataGR <- list()
#     class(dataGR) <- "preMetricSegments"
#     
#     dataGR$segments <- GRanges(seqnames = "chr1", 
#                         ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
#                                          end=c(29, 49, 55, 100, 110, 199, 220, 230)))
#     
#     values(dataGR$segments) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
#                                   File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
#                                   File2 = c(NA, 0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
#     
#     expMatrix <- matrix(c(0, 0.387652570, 0.387652570, 0), byrow = TRUE, 
#                             ncol = 2, nrow = 2)
#     colnames(expMatrix) <- c("File1", "File2")
#     rownames(expMatrix) <- c("File1", "File2")
#     
#     expect_equivalent(calculateWeightedEuclideanDistance(dataGR, minThreshold = 0.04), expMatrix, 
#                         tolerance=0.000001)
# })
# 
# 
# test_that("calculateWeightedEuclideanDistance() must return good results 02", {
#     
#     dataGR <- list()
#     class(dataGR) <- "preMetricSegments"
#     
#     dataGR$segments <- GRanges(seqnames = "chr1", 
#                                ranges = IRanges(start=c(1, 30, 50, 56, 101, 180, 200, 221),
#                                                 end=c(29, 49, 55, 100, 110, 199, 220, 230)))
#     
#     values(dataGR$segments) <- DataFrame(included = c(rep(TRUE, 2), rep(FALSE, 2), rep(TRUE, 4)), 
#                                          File1 = c(0.1, 0.1, 0.1, NA, NA, NA, 0.5, NA), 
#                                          File2 = c(NA, 0.15, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4))
#     
#     expMatrix <- matrix(c(0, 0.174485599340559, 0.174485599340559, 0), 
#                         byrow = TRUE, ncol = 2, nrow = 2)
#     colnames(expMatrix) <- c("File1", "File2")
#     rownames(expMatrix) <- c("File1", "File2")
#     
#     expect_equivalent(calculateWeightedEuclideanDistance(dataGR, minThreshold = 0.07), expMatrix, 
#                       tolerance=0.000001)
# })
# 
# 
# 
# ### Tests plotMetric() results
# 
# context("plotMetric() results")
# 
# test_that("plotMetric() must return error when type wrong", {
#     
#     demo <- GRangesList()
#     demo[["sample01"]] <- GRanges(seqnames="chr1",
#                                   ranges =  IRanges(start=c(100, 300, 800), end=c(200, 500, 900)), 
#                                   strand =  "*", state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#     demo[["sample02"]] <- GRanges(seqnames="chr1",
#                                   ranges=IRanges(start=c(150, 600, 1000), end=c(250, 700, 1500)), 
#                                   strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#     
#     
#     metric <- calculateOverlapMetric(segmentData=demo, 
#                                      method="szymkiewicz")
#     
#     
#     error_message <- gettextf("'arg' should be one of %s", 
#                               paste(dQuote(c("ALL", "AMPLIFICATION", 
#                                              "DELETION")), 
#                                     collapse=", "))
#     
#     expect_error(plotMetric(metric=metric,  type="SAVE"), 
#                  error_message)
# })
# 
# 
# test_that("plotMetric() must return error when metric is not CNVMetric object", {
#     
#     error_message <- "\'metric\' must be a CNVMetric object."
#     
#     expect_error(plotMetric(metric="TEST01",  type="SAVE"), 
#                  error_message)
# })
# 
# 
# test_that("plotMetric() must return error when colorRange is vector of single letters", {
#     
#     
#     demo <- GRangesList()
#     demo[["sample01"]] <- GRanges(seqnames = "chr1",
#                                   ranges =  IRanges(start = c(100, 300, 800), end = c(200, 500, 900)), 
#                                   strand =  "*", state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#     demo[["sample02"]] <- GRanges(seqnames = "chr1",
#                                   ranges =  IRanges(start = c(150, 600, 1000), end = c(250, 700, 1500)), 
#                                   strand =  "*", state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#     
#     metric <- calculateOverlapMetric(segmentData = demo, 
#                                      method = "szymkiewicz")
#     
#     error_message <- "\'colorRange\' must be be a vector of 2 valid color names."
#     
#     expect_error(plotMetric(metric=metric,  type="AMPLIFICATION",
#                                    colorRange=c("a", "b")), 
#                  error_message)
# })
# 
# test_that("plotMetric() must return error when colorRange is vector of one color", {
#     
#     demo <- GRangesList()
#     demo[["sample01"]] <- GRanges(seqnames="chr1",
#                                   ranges=IRanges(start=c(100, 300, 800), end = c(200, 500, 900)), 
#                                   strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#     demo[["sample02"]] <- GRanges(seqnames="chr1",
#                                   ranges=IRanges(start=c(150, 600, 1000), end=c(250, 700, 1500)), 
#                                   strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#     
#     metric <- calculateOverlapMetric(segmentData = demo, 
#                                      method = "szymkiewicz")
#     
#     error_message <- "\'colorRange\' must be a vector of 2 color names."
#     
#     expect_error(plotMetric(metric=metric,  type="AMPLIFICATION",
#                                    colorRange=c("red")), 
#                  error_message)
# })


test_that("plotMetric() must return error when type==ALL and filename given", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames="chr1",
                                  ranges=IRanges(start=c(100, 300, 800), end = c(200, 500, 900)), 
                                  strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames="chr1",
                                  ranges=IRanges(start=c(150, 600, 1000), end=c(250, 700, 1500)), 
                                  strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    metric <- calculateOverlapMetric(segmentData = demo, 
                                     method = "szymkiewicz")
    
    error_message <- "\'type\' cannot be \'ALL\' when filename argument is used."
    
    expect_error(plotMetric(metric=metric,  type="ALL",
                                   filename="test.pdf"), 
                 error_message)
})

test_that("plotMetric() must return a gtable when graph for amplification", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames="chr1",
                                  ranges=IRanges(start=c(100, 300, 800), end = c(200, 500, 900)), 
                                  strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames="chr1",
                                  ranges=IRanges(start=c(150, 600, 1000), end=c(250, 700, 1500)), 
                                  strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    metric <- calculateOverlapMetric(segmentData = demo, 
                                     method = "szymkiewicz")
    
    result <- plotMetric(metric=metric,  type="AMPLIFICATION")
    
    expect_is(object=result, class="gtable")
})

test_that("plotMetric() must return a gtable when graph for deletion", {
    
    demo <- GRangesList()
    demo[["sample01"]] <- GRanges(seqnames="chr1",
                                  ranges=IRanges(start=c(100, 300, 800), end = c(200, 500, 900)), 
                                  strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    demo[["sample02"]] <- GRanges(seqnames="chr1",
                                  ranges=IRanges(start=c(150, 600, 1000), end=c(250, 700, 1500)), 
                                  strand="*", state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
    
    metric <- calculateOverlapMetric(segmentData = demo, 
                                     method = "szymkiewicz")
    
    result <- plotMetric(metric=metric,  type="DELETION")
    
    expect_is(object=result, class="gtable")
})



