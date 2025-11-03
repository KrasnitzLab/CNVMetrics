
#' @title Calculate metric using the log2ratio values between two samples.
#'
#' @description Calculate a specific metric using the level of
#' amplification/deletion, in log2 ratio, between two samples.
#'
#' @param entry a \code{list} which contains the row and column indexes
#' (always in this order) of
#' the metric in the final matrix. Those values correspond to the positions
#' of the two samples used
#' to calculate the metric in the \code{GRangesList} (\code{segmentData}).
#'
#' @param segmentData a \code{GRangesList} that contains a collection of
#' genomic ranges representing copy number events, including amplified/deleted
#' status, from at least 2 samples. All samples must have a metadata column
#' called '\code{log2ratio}' with the log2ratio values.
#' @param method a \code{character} string representing the metric to be
#' used ('\code{weightedEuclideanDistance}').
#'
#' @param minThreshold a single \code{numeric} setting the minimum value
#' to consider two segments as different during the metric calculation. If the
#' absolute difference is below or equal to threshold, the difference will be
#' replaced by zero.
#'
#' @param bedExclusion an optional \code{GRanges} containing the regions
#' that have to be excluded for the metric calculation or \code{NULL}.
#'
#' @details
#'
#' The method calculates a specified metric using overlapping
#' regions between the samples. Only regions corresponding to the type
#' specified by user are used in the calculation of the metric. The strand of
#' the regions is not taken into account while
#' calculating the metric.
#'
#' The Sorensen metric is calculated by dividing twice the size of
#' the intersection by the sum of the size of the two sets. If the sum of
#' the size of the two sets is zero; the value \code{NA} is
#' returned instead.
#'
#'
#' @return a \code{list} containing 1 entry:
#' \itemize{
#'  \item{\code{metric} a \code{data.frame}, which contains 3 columns. The 2
#' first columns, called \code{row} and \code{column} correspond to the
#' indexes of the metric in the final matrix. Those
#' 2 first columns match to the \code{entry} parameter. The third column,
#' called \code{metric},
#' contains the values of the specified metric for each combination.
#' If the metric cannot be calculated, \code{NA} is present. }
#' }
#'
#' @examples
#'
#' ## Load required package to generate the two samples
#' require(GenomicRanges)
#'
#' ## Create a GRangesList object with 3 samples
#' ## The stand of the regions doesn't affect the calculation of the metric
#' demo <- GRangesList()
#'
#' ## Generate two samples with log2value information as a metadata column
#' demo[["sample01"]]  <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(100, 201, 400),
#'     end=c(200, 350, 500)), strand="*",
#'     log2ratio=c(1.1111, 2.2222, -0.9999))
#' demo[["sample02"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(150, 200, 450),
#'     end=c(250, 350, 500)), strand="*",
#'     log2ratio=c(2.2121, 1.1212, -1.3939))
#'
#' ## The 2 samples used to calculate the metric
#' entries <- data.frame(row=c(2), col=c(1))
#'
#' ## Calculate weighted Euclidean distance
#' CNVMetrics:::calculateOneLog2valueMetricT(entry=entries,
#'     segmentData=demo, method="weightedEuclideanDistance",
#'     minThreshold=0.2, bedExclusion=NULL)
#'
#'
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
calculateOneLog2valueMetricT <- function(entry, segmentData, method,
                                            minThreshold, bedExclusion) {
    entries <- split(entry, entry$col)

    results <- list()

    for (data in entries) {

        result <- data
        result$metric <- rep(NA, nrow(result))

        sample02 <- segmentData[[data$col[1]]]

        for(i in seq_len(nrow(data))) {
            sample01 <- segmentData[[data$row[i]]]

            # Obtain the disjoint segments with log2ratio values in
            # metadata columns
            disjoinR <- createDisjoinSegmentsForTwoSamples(
                            segmentDataSample1=sample01,
                            segmentDataSample2=sample02,
                            bedExclusion=bedExclusion)

            if (length(sample01) > 0 && length(sample02) > 0) {
                result$metric[i] <- switch(method,
                            weightedEuclideanDistance=
                                calculateWeightedEuclideanDistanceFor2Samples(
                                        segmentData=disjoinR,
                                        minThreshold=minThreshold))
            }
        }

        results[[length(results) + 1]] <- result
    }

    results <- do.call(rbind, results)

    return(result <- list(metric=results))
}


#' @title Generate common segments to enable calculation of metrics on
#' two segmented samples.
#'
#' @description The two segments are gathered together, including excluded
#' regions when specified, and a disjoint operation is done to create a
#' collection of non-overlapping ranges. The ranges overlapping the excluded
#' regions are marked as so to be removed from future analysis. The log2value
#' of each samples are assigned to the new disjointed segments for each sample
#' in the metadata columns.
#'
#' @param segmentDataSample1 a \code{GRanges}, the segments from the first
#' sample.
#'
#' @param segmentDataSample2 a \code{GRanges}, the segments from the second
#' sample.
#'
#' @param bedExclusion a \code{GRanges}, the regions that must be
#' excluded from the analysis. Default: \code{NULL}.
#'
#' @return a \code{GRanges} containing the common segment information for the
#' two samples. The log2ration value are present, for the two samples, in
#' the metadata columns. When there is not log2ratio value for one sample,
#' NA is the assigned value. A metadata column also specifies if the segments
#' should be included in the analysis.
#'
#' @examples
#'
#' ## Load required package to generate the two samples
#' require(GenomicRanges)
#'
#' # Create first Granges representing first sample
#' sample01 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(100, 201, 400), end=c(200, 350, 500)),
#'     strand="*", log2ratio=c(0.3091175, 0.4582058, -0.3798390))
#'
#' # Create second Granges representing second sample
#' sample02 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(150, 200, 450), end=c(250, 350, 500)),
#'     strand="*", log2ratio=c(0.222174, 0.3282156, -0.2728292))
#'
#' # Create disjoint segment using the 2 samples and without any region
#' # excluded from the analysis (parameter bedExclusion set to null)
#' CNVMetrics:::createDisjoinSegmentsForTwoSamples(segmentDataSample1=sample01,
#'     segmentDataSample2=sample02, bedExclusion=NULL)
#'
#' @author Astrid Deschênes
#' @importFrom GenomicRanges disjoin findOverlaps elementMetadata
#' @importFrom S4Vectors queryHits subjectHits values<-
#' @importFrom magrittr %>%
#' @keywords internal
createDisjoinSegmentsForTwoSamples <- function(segmentDataSample1,
                                                segmentDataSample2,
                                                bedExclusion=NULL) {

    results <- disjoin(c(segmentDataSample1, segmentDataSample2))
    results$included <- TRUE

    ## Add information about excluded regions
    ## When a segment overlaps with an excluded region, it is marked as
    ## excluded
    if (!is.null(bedExclusion) && (length(bedExclusion) > 0)) {
        olaps <- findOverlaps(results, bedExclusion)

        if (length(olaps) > 0) {
            results[queryHits(olaps)]$included <- FALSE
        }
    }

    ## Assign the log2value of the 2 samples for each new segment as
    ## metadata columns
    segList <- list(segmentDataSample1, segmentDataSample2)
    for (i in seq_len(length(segList))) {
        olaps <- findOverlaps(results, segList[[i]])
        temp <- elementMetadata(results)
        sampleName <- paste0("sample_", i)
        temp[, sampleName] <- NA
        temp[queryHits(olaps), sampleName] <-
            segList[[i]]$log2ratio[subjectHits(olaps)]
        values(results) <- temp
    }

    return(results)
}



#' @title Calculate Weighted Euclidean distance-based metric between samples.
#'
#' @description The weighted Euclidean distance-based metric corresponds to the
#' euclidean distance between 2 samples multiplied by the natural logarithm
#' of the number of bases of the analyzed segment. The final metric is 1 over
#' 1 added to the
#' squared sum of the values obtained for all segments that are not
#' excluded of the analysis.
#'
#' @param segmentData a \code{list} marked as a \code{preMetricSegments}
#' \code{class} that contains the disjoint segment information from 2
#' samples and the log2ratio values of the samples in the metadata columns.
#'
#' @param minThreshold a single \code{numeric} setting the minimum value
#' to consider two segments as different for the metric calculation. If the
#' absolute difference is below or equal to threshold, the value will be
#' replaced by zero.
#'
#' @return a \code{numeric} representing the weighted euclidean distance
#' between the two samples. If the distance cannot be calculated as the two
#' samples don't share any segments with log2ratio value, the value NA is
#' assigned.
#'
#' @details
#'
#' The weighted euclidean distance is
#' \eqn{1/(1 + (\sum((x_i - y_i)^2 * log2(nbrBases_i))^0.5)}
#' where \code{x} and \code{y} are the
#' values of 2 samples for a specific segment \code{i} and \code{nbrBases} the
#' number of bases of the segment \code{i}.
#'
#'
#' @examples
#'
#' ## Load required package to generate the two samples
#' require(GenomicRanges)
#'
#' # Create first Granges representing first sample
#' sample01 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(100, 201, 400), end=c(200, 350, 500)),
#'     strand="*", log2ratio=c(0.3091175, 0.4582058, -0.3798390))
#'
#' # Create second Granges representing second sample
#' sample02 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(150, 200, 450), end=c(250, 350, 500)),
#'     strand="*", log2ratio=c(0.222174, 0.3282156, -0.2728292))
#'
#' # Create disjoint segment using the 2 samples and without any region
#' # excluded from the analysis (parameter bedExclusion set to null)
#' disjoinGRange <- CNVMetrics:::createDisjoinSegmentsForTwoSamples(
#'     segmentDataSample1=sample01, segmentDataSample2=sample02,
#'     bedExclusion=NULL)
#'
#' ## Calculate the weighted ecucidean distance between the two samples
#' CNVMetrics:::calculateWeightedEuclideanDistanceFor2Samples(
#'     segmentData=disjoinGRange, minThreshold=0.2)
#'
#' @author Astrid Deschênes
#' @importFrom GenomicRanges elementMetadata
#' @importFrom IRanges ranges width
#' @keywords internal
calculateWeightedEuclideanDistanceFor2Samples <- function(segmentData,
                                                            minThreshold) {

    names <- colnames(elementMetadata(segmentData))
    names <- names[names != "included"]

    incResults <- elementMetadata(segmentData[segmentData$included, ])
    temp01 <- incResults[, c(names[1])] - incResults[, c(names[2])]

    final <- NA

    ## Only calculate when at least one value is present
    if (!all(is.na(temp01))) {
        incWidth <- width(ranges(segmentData[segmentData$included, ]))

        ## Set values to zero when lower than threshold
        tempPos <- which(abs(temp01) <= minThreshold)
        if (length(tempPos) > 0) {
            temp01[tempPos] <- 0.0
        }

        ## Calculate metric
        temp01 <- temp01 * temp01 * log2(incWidth)
        final <- 1/(1 + sum(temp01, na.rm=TRUE) ^ (1/2))
    }

    return(final)
}
