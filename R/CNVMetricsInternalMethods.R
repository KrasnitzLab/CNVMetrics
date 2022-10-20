
#' @title Parameters validation for the \code{\link{calculateOverlapMetric}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{calculateOverlapMetric}} function.
#'
#' @param states a \code{vector} of \code{character} string with at least one
#' entry. The strings are representing the states that will be analyzed.
#'
#' @param nJobs a single positive \code{integer} specifying the number of
#' worker jobs to create in case of distributed computation.
#'
#' @return \code{0}.
#'
#' @examples
#'
#'
#' ## Return zero as all parameters are valid
#' CNVMetrics:::validateCalculateOverlapMetricParameters(
#'     states="GAIN", nJobs=1)
#'
#' @author Astrid Deschênes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validateCalculateOverlapMetricParameters <- function(states, nJobs) {

    ## Validate that nJobs is an positive integer
    if (!(isSingleInteger(nJobs) || isSingleNumber(nJobs)) ||
            as.integer(nJobs) < 1) {
            stop("nJobs must be a positive integer")
    }

    ## Validate that nJobs is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nJobs) != 1) {
            stop("nJobs must be 1 on a Windows system")
    }

    ## At least one state must be present
    if (!is.vector(states) | ! is.character(states) | length(states) < 1){
            stop("the \'states\' argument must be a vector of strings ",
                    "with at least one value")
    }

    return(0L)
}


#' @title Parameters validation for the \code{\link{calculateLog2ratioMetric}}
#' function
#'
#' @description Validation of all parameters needed by the public
#' \code{\link{calculateLog2ratioMetric}} function.
#'
#' @param minThreshold a single positive \code{numeric} setting the minimum
#' value to consider two segments as different during the metric calculation.
#' If the absolute difference is below or equal to threshold, the difference
#' will be replaced by zero.
#'
#' @param excludedRegions an optional \code{GRanges} containing the regions
#' that have to be excluded for the metric calculation or \code{NULL}.
#'
#' @param nJobs a single positive \code{integer} specifying the number of
#' worker jobs to create in case of distributed computation.
#'
#' @return \code{0}.
#'
#' @examples
#'
#'
#' ## Return zero as all parameters are valid
#' CNVMetrics:::validatecalculateLog2ratioMetricParameters(
#'     minThreshold=0.9, excludedRegions=NULL, nJobs=1)
#'
#' @author Astrid Deschênes
#' @importFrom S4Vectors isSingleInteger isSingleNumber
#' @keywords internal
validatecalculateLog2ratioMetricParameters <- function(minThreshold,
                                        excludedRegions, nJobs) {

    ## Validate that nJobs is an positive integer
    if (!(isSingleInteger(nJobs) || isSingleNumber(nJobs)) ||
        as.integer(nJobs) < 1) {
        stop("nJobs must be a positive integer")
    }

    ## Validate that nJobs is set to 1 on Windows system
    if (Sys.info()["sysname"] == "Windows" && as.integer(nJobs) != 1) {
        stop("nJobs must be 1 on a Windows system")
    }

    ## The minThreshold must be a positive numeric value
    if (!is.numeric(minThreshold) | minThreshold < 0.0) {
        stop("the \'minThreshold\' argument must be a positive numeric value")
    }

    ## The minThreshold must be a positive numeric value
    if (!is.null(excludedRegions) & !is(excludedRegions, "GRanges")) {
        stop("the \'excludedRegions\' argument must ",
                "a \'Granges\' object or NULL")
    }

    return(0L)
}


#' @title Plot one graph related to one set of metrics.
#'
#' @description Plot one heatmap of one set of metrics present in a
#' a \code{CNVMetric} object.
#'
#' @param metric a \code{CNVMetric} object containing the metrics calculated
#' by \code{calculateOverlapMetric}.
#'
#' @param type a \code{character} string indicating which graph to generate.
#' This should be (an unambiguous abbreviation of) one of
#' "\code{AMPLIFICATION}" or "\code{DELETION}" or "\code{LOG2RATIO}".
#'
#' @param show_colnames a \code{boolean} specifying if column names are
#' be shown.
#'
#' @param silent a \code{boolean} specifying if the plot should not be drawn.
#'
#' @param \ldots further arguments passed to
#' \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}} method.
#'
#' @return a \code{gtable} object containing the heatmap for the specified
#' metric.
#'
#' @seealso
#'
#' The default method  \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}}.
#'
#' @examples
#'
#' ## Load required package to generate the samples
#' require(GenomicRanges)
#'
#' ## Create a GRangesList object with 3 samples
#' ## The stand of the regions doesn't affect the calculation of the metric
#' demo <- GRangesList()
#' demo[["sample01"]] <- GRanges(seqnames = "chr1",
#'     ranges =  IRanges(start = c(1905048, 4554832, 31686841),
#'     end = c(2004603, 4577608, 31695808)), strand =  "*",
#'     state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#'
#' demo[["sample02"]] <- GRanges(seqnames = "chr1",
#'     ranges =  IRanges(start = c(1995066, 31611222, 31690000),
#'     end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
#'     state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#'
#' ## The amplified region in sample03 is a subset of the amplified regions
#' ## in sample01
#' demo[["sample03"]] <- GRanges(seqnames = "chr1",
#'     ranges =  IRanges(start = c(1906069, 4558838),
#'     end = c(1909505, 4570601)), strand =  "*",
#'     state = c("AMPLIFICATION", "DELETION"))
#'
#' ## Calculating Sorensen metric
#' metric <- calculateOverlapMetric(demo, method="sorensen")
#'
#' ## Plot amplification metrics using darkorange color
#' CNVMetrics:::plotOneMetric(metric, type="AMPLIFICATION",
#'     colorRange=c("white", "darkorange"), show_colnames=FALSE, silent=TRUE)
#'
#' @author Astrid Deschênes
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @importFrom methods hasArg
#' @importFrom stats as.dist
#' @import GenomicRanges
#' @encoding UTF-8
#' @keywords internal
plotOneMetric <- function(metric, type, colorRange, show_colnames,
                                            silent, ...)
{
    ## Extract matrix with metric values
    metricMat <- metric[[type]]

    ## Extract extra arguments
    dots <- list(...)

    ## Prepare matrix by filling upper triangle
    diag(metricMat) <- 1.0

    ## If clustering distances are not present in the arguments,
    ## the distance used is based on the samples distance
    if ((!("clustering_distance_cols" %in% names(dots))) &&
        (!("clustering_distance_rows" %in% names(dots)))) {
        ## Prepare matrix to be able to calculate distance
        metricMat[lower.tri(metricMat) & is.na(metricMat)] <- 0.0
        metricDist <- as.dist(1-metricMat)

        dots[["clustering_distance_cols"]] <- metricDist
        dots[["clustering_distance_rows"]] <- metricDist
    }

    ## Prepare matrix by filling upper triangle
    metricMat[upper.tri(metricMat)] <- t(metricMat)[upper.tri(metricMat)]
    metricMat[is.na(metricMat)] <- 0.0

    ## Prepare main title (might not be used if main argument given by user)
    if (!hasArg("main")) {
        metricInfo <- switch(attributes(metric)$metric,
                    "szymkiewicz"="Szymkiewicz-Simpson",
                    "sorensen"="Sorensen",
                    "jaccard"="Jaccard",
                    "weightedEuclideanDistance"="Weighted Euclidean Distance")
        dots[["main"]] <- paste0(type, " - ", metricInfo, " metric")
    }

    ## Create heatmap
    ## If color information given, that information is used to create graph
    ## If main title given, that information is used to create graph
    if (!hasArg("breaks") && !hasArg("color")) {
        ## Create color palette using colorRange parameter
        colors <- colorRampPalette(colorRange)(255)
        breaks <-  seq(0, 1, length.out=255)

        dots[["color"]] <- colors
        dots[["breaks"]] <- breaks
    }

    ## Add arguments
    dots[["mat"]] <- metricMat
    dots[["show_colnames"]] <- show_colnames
    dots[["silent"]] <- silent

    ## Create heatmap
    do.call(what="pheatmap", args=dots)[[4]]
}

