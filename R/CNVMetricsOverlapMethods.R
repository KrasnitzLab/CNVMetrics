
#' @title Calculate metric using overlapping amplified/deleted regions
#' 
#' @description Calculate a specific metric using overlapping 
#' amplified/deleted regions between to samples. The metric is calculated for
#' the amplified and deleted regions separately. When more than 2 samples are 
#' present, the metric is calculated for each sample pair.
#' 
#' @param segmentData a \code{GRangesList} that contains a collection of 
#' genomic ranges representing copy number events, including amplified/deleted 
#' status, from at least 2 samples. All samples must have a metadata column 
#' called '\code{state}' with amplified regions identified as 
#' '\code{AMPLIFICATION}' and deleted regions identified as '\code{DELETION}'; 
#' regions with different identifications will not be used in the
#' calculation of the metric. 
#' 
#' @param method a \code{character} string representing the metric to be used. 
#' This should be (an unambiguous abbreviation of) one of "sorensen" or 
#' "szymkiewicz". Default: "sorensen".
#' 
#' @details 
#' 
#' The two methods each estimate the overlap between paired samples. They use 
#' different metrics, all in the range [0, 1] with 0 indicating no overlap.
#' 
#' The available metrics are (written for two GRanges):
#' 
#' \code{sorensen}:
#' 
#' This metric is calculated by dividing twice the size of the intersection 
#' by the sum of the size of the two sets. 
#' With this metric, an overlap metric value of 1 is only obtained when the
#' two samples are identical. 
#' 
#' \code{szymkiewicz}:
#' 
#' This metric is calculated by dividing the size of the intersection 
#' by the size of the smallest set. With this metric, if one set is a 
#' subset of the other set, the overlap metric value is 1.
#' 
#' @return a \code{list} of class "\code{CNVMetric}". This list has
#' the following components:
#' \itemize{
#' \item{\code{AMPLIFICATION}}{ a lower-triangular \code{matrix} with the 
#'     results of the selected metric on the amplified regions for each paired
#'     samples. The value \code{NA} is present when the metric cannot be 
#'     calculated. The value \code{NA} is also present in the top-triangular 
#'     section of the matrix.
#'  }
#'  \item{\code{DELETION}}{ a lower-triangular \code{matrix} with the 
#'     results of the selected metric on the deleted regions for each paired
#'     samples. The value \code{NA} is present when the metric cannot be 
#'     calculated. The value \code{NA} is also present in the top-triangular 
#'     section of the matrix.
#' }}
#' 
#' The object has the following attributes (besides "class" equal 
#' to "CNVMetric"):
#' \itemize{
#' \item{\code{metric}}{ the metric used for the calculation.
#'  } 
#' \item{\code{names}}{ the names of the two matrix containing the metrics for
#' the amplified and deleted regions.
#' }}         
#' 
#' 
#' @references 
#' 
#' Sørensen, Thorvald. n.d. “A Method of Establishing Groups of Equal 
#' Amplitude in Plant Sociology Based on Similarity of Species and Its 
#' Application to Analyses of the Vegetation on Danish Commons.” 
#' Biologiske Skrifter, no. 5: 1–34.
#' 
#' Vijaymeena, M. K, and Kavitha K. 2016. “A Survey on Similarity Measures in 
#' Text Mining.” Machine Learning and Applications: An International 
#' Journal 3 (1): 19–28. doi: \url{https://doi.org/10.5121/mlaij.2016.3103}
#' 
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
#' calculateOverlapRegionsMetric(demo, method="sorensen")
#' 
#' ## Calculating Szymkiewicz-Simpson metric
#' calculateOverlapRegionsMetric(demo, method="szymkiewicz")
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges 
#' @export
calculateOverlapRegionsMetric <- function(segmentData, 
                                    method=c("sorensen", "szymkiewicz")) {
    
    method <- match.arg(method)
    
    names <- names(segmentData)
    nb <- length(names)
    
    ## At least 2 samples must be present
    if (nb < 2) {
        stop("at least 2 samples must be present in the segmentData")
    }
    
    ## All samples must have a metadata column called 'state' with
    ## AMPLIFICATION/DELETION status 
    if (!all(vapply(segmentData, 
                    FUN = function(x) {"state" %in% colnames(mcols(x))},
                    FUN.VALUE = logical(1)))) {
        stop(paste0("at least one sample doesn't have a metadata column ", 
                    "called \'state\'"))
    }
    
    results <- list()
    
    for(type in c("AMPLIFICATION", "DELETION")) {
        
        dataTMP <- matrix(rep(NA, nb^2), nrow=nb)
        rownames(dataTMP) <- names
        colnames(dataTMP) <- names
        
        for(i in seq_len(nb)[-1]) {
            for(j in seq_len(i-1)) {
                dataTMP[i, j] <- calculateOverlapMetric(
                    sample01=segmentData[[names[i]]], 
                    sample02=segmentData[[names[j]]],
                    method=method, type=type)
            }
        }
        results[[type]] <- dataTMP
    }
    
    # Return a list marked as an CNVMetric class containing:
    # 1- the metric results for the amplified regions
    # 2- the metric results for the deleted regions
    class(results) <- "CNVMetric"
    attr(results, 'metric') <- method
    
    return(results)
}


#' @title Plot metrics based on overlapping amplified/deleted regions
#' 
#' @description Plot a heatmap of the metrics based on overlapping 
#' amplified/deleted regions.
#' 
#' @param metric a \code{CNVMetric} object containing the metrics calculated
#' by \code{calculateOverlapRegionsMetric}.
#' 
#' @param type a \code{character} string indicating which graph to generate. 
#' This should be (an unambiguous abbreviation of) one of "\code{BOTH}", 
#' "\code{AMPLIFICATION}" or "\code{DELETION}". Default: "\code{BOTH}".
#' 
#' 
#' @return a \code{gtable} object containing the heatmap of the specified 
#' metric. TODO
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
#' metric <- calculateOverlapRegionsMetric(demo, method="sorensen")
#' 
#' ## Plot both amplification and deletion metrics
#' plotOverlapMetric(metric, type="BOTH")
#' 
#' @seealso 
#' 
#' The default method  \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}}.
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @importFrom grDevices colorRampPalette
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @import GenomicRanges
#' @export
plotOverlapMetric <- function(metric, 
                              type=c("BOTH", "AMPLIFICATION", "DELETION")) {
    
    ## Validate that the metric parameter is a CNVMetric object
    if (!is.CNVMetric(metric)) {
        stop("\'metric\' must be a CNVMetric object.")
    }
    
    ## Assign type parameter
    type <- match.arg(type)
    
    ## Extract the type of metric
    metricInfo <- attributes(metric)$metric
    
    plot_list <- list()
    
    colors <- colorRampPalette(brewer.pal(9, "Blues"))(255)
    breaks <-  seq(0, 1, length.out = 255)
    
    ## Amplification
    if (type %in% c("AMPLIFICATION", "BOTH")) {
        ampMatrix <- metric$AMPLIFICATION
        diag(ampMatrix) <- 1.0
        ampMatrix[upper.tri(ampMatrix)] <- t(ampMatrix)[upper.tri(ampMatrix)]
        
        rownames(ampMatrix) <-  rownames(metric$AMPLIFICATION)
        colnames(ampMatrix) <- NULL
        plot_list[["AMPLIFICATION"]] <- pheatmap(ampMatrix, cluster_rows=TRUE, 
                                                 cluster_cols=TRUE,
                                                 main="Amplification",
                                                 color=colors, 
                                                 breaks=breaks)[[4]]    
    }
    
    ## Deletion
    if (type %in% c("DELETION", "BOTH")) {
        delMatrix <- metric$DELETION
        diag(delMatrix) <- 1.0
        delMatrix[upper.tri(delMatrix)] <- t(delMatrix)[upper.tri(delMatrix)]
        
        rownames(delMatrix) <-  rownames(metric$DELETION)
        colnames(delMatrix) <- NULL
        
        plot_list[["DELETION"]] <- pheatmap(delMatrix, cluster_rows=TRUE, 
                                            cluster_cols=TRUE,
                                            main="Deletion",
                                            color=colors, breaks=breaks)[[4]]   
    }
    
    n_col <- ifelse(type == "BOTH", 2, 1)
    
    grid.arrange(arrangeGrob(grobs= plot_list, ncol=n_col))
}