#' 
#' #' @title Calculate metric using overlapping amplified/deleted regions
#' #' 
#' #' @description Calculate a specific metric using overlapping 
#' #' amplified/deleted regions between to samples. The metric is calculated for
#' #' the amplified and deleted regions separately. When more than 2 samples are 
#' #' present, the metric is calculated for each sample pair.
#' #' 
#' #' @param segmentData a \code{GRangesList} that contains a collection of 
#' #' genomic ranges representing copy number events, including amplified/deleted 
#' #' status, from at least 2 samples. All samples must have a metadata column 
#' #' called '\code{state}' with amplified regions identified as 
#' #' '\code{AMPLIFICATION}' and deleted regions identified as '\code{DELETION}'; 
#' #' regions with different identifications will not be used in the
#' #' calculation of the metric. 
#' #' 
#' #' @param method a \code{character} string representing the metric to be used. 
#' #' This should be (an unambiguous abbreviation of) one of "sorensen", 
#' #' "szymkiewicz" or "jaccard". Default: "sorensen".
#' #' 
#' #' @details 
#' #' 
#' #' The two methods each estimate the overlap between paired samples. They use 
#' #' different metrics, all in the range [0, 1] with 0 indicating no overlap.
#' #' 
#' #' The available metrics are (written for two GRanges):
#' #' 
#' #' \code{sorensen}:
#' #' 
#' #' This metric is calculated by dividing twice the size of the intersection 
#' #' by the sum of the size of the two sets. 
#' #' With this metric, an overlap metric value of 1 is only obtained when the
#' #' two samples are identical. 
#' #' 
#' #' \code{szymkiewicz}:
#' #' 
#' #' This metric is calculated by dividing the size of the intersection 
#' #' by the size of the smallest set. With this metric, if one set is a 
#' #' subset of the other set, the overlap metric value is 1.
#' #' 
#' #' \code{jaccard}:
#' #' 
#' #' This metric is calculated by dividing the size of the intersection 
#' #' by the size of the union of the two sets. With this metric, an overlap 
#' #' metric value of 1 is only obtained when the two samples are identical. 
#' #' 
#' #' @return an object of class "\code{CNVMetric}" which contains the calculated
#' #' metric. This object is a list with the following components:
#' #' \itemize{
#' #' \item{\code{AMPLIFICATION}}{ a lower-triangular \code{matrix} with the 
#' #'     results of the selected metric on the amplified regions for each paired
#' #'     samples. The value \code{NA} is present when the metric cannot be 
#' #'     calculated. The value \code{NA} is also present in the top-triangular 
#' #'     section, as well as the diagonal, of the matrix.
#' #' }
#' #' \item{\code{DELETION}}{ a lower-triangular \code{matrix} with the 
#' #'     results of the selected metric on the deleted regions for each paired
#' #'     samples. The value \code{NA} is present when the metric cannot be 
#' #'     calculated. The value \code{NA} is also present in the top-triangular 
#' #'     section, as well as the diagonal, of the matrix.
#' #' }}
#' #' 
#' #' The object has the following attributes (besides "class" equal 
#' #' to "CNVMetric"):
#' #' \itemize{
#' #' \item{\code{metric}}{ the metric used for the calculation.
#' #' } 
#' #' \item{\code{names}}{ the names of the two matrix containing the metrics for
#' #' the amplified and deleted regions.
#' #' }}         
#' #' 
#' #' 
#' #' @references 
#' #' 
#' #' Sørensen, Thorvald. n.d. “A Method of Establishing Groups of Equal 
#' #' Amplitude in Plant Sociology Based on Similarity of Species and Its 
#' #' Application to Analyses of the Vegetation on Danish Commons.” 
#' #' Biologiske Skrifter, no. 5: 1–34.
#' #' 
#' #' Vijaymeena, M. K, and Kavitha K. 2016. “A Survey on Similarity Measures in 
#' #' Text Mining.” Machine Learning and Applications: An International 
#' #' Journal 3 (1): 19–28. doi: \url{https://doi.org/10.5121/mlaij.2016.3103}
#' #' 
#' #' Jaccard, P. (1912), The Distribution of the Flora in the Alpine Zone.  
#' #' New Phytologist, 11: 37-50. 
#' #' doi: \url{https://doi.org/10.1111/j.1469-8137.1912.tb05611.x}
#' #' 
#' #' @examples
#' #'
#' #' ## Load required package to generate the samples
#' #' require(GenomicRanges)
#' #' 
#' #' ## Create a GRangesList object with 3 samples
#' #' ## The stand of the regions doesn't affect the calculation of the metric
#' #' demo <- GRangesList()
#' #' demo[["sample01"]] <- GRanges(seqnames = "chr1", 
#' #'     ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
#' #'     end = c(2004603, 4577608, 31695808)), strand =  "*",
#' #'     state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#' #' 
#' #' demo[["sample02"]] <- GRanges(seqnames = "chr1", 
#' #'     ranges =  IRanges(start = c(1995066, 31611222, 31690000), 
#' #'     end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
#' #'     state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#' #' 
#' #' ## The amplified region in sample03 is a subset of the amplified regions 
#' #' ## in sample01
#' #' demo[["sample03"]] <- GRanges(seqnames = "chr1", 
#' #'     ranges =  IRanges(start = c(1906069, 4558838), 
#' #'     end = c(1909505, 4570601)), strand =  "*",
#' #'     state = c("AMPLIFICATION", "DELETION"))
#' #' 
#' #' ## Calculating Sorensen metric
#' #' calculateOverlapMetric(demo, method="sorensen")
#' #' 
#' #' ## Calculating Szymkiewicz-Simpson metric
#' #' calculateOverlapMetric(demo, method="szymkiewicz")
#' #' 
#' #' @author Astrid Deschênes, Pascal Belleau
#' #' @import GenomicRanges 
#' #' @encoding UTF-8
#' #' @export
#' calculateOverlapMetric <- function(segmentData, 
#'                                     method=c("sorensen", "szymkiewicz", 
#'                                                 "jaccard")) {
#'     
#'     method <- match.arg(method)
#'     
#'     names <- names(segmentData)
#'     nb <- length(names)
#'     
#'     ## At least 2 samples must be present
#'     if (nb < 2) {
#'         stop("at least 2 samples must be present in the segmentData")
#'     }
#'     
#'     ## All samples must have a metadata column called 'state' with
#'     ## AMPLIFICATION/DELETION status 
#'     if (!all(vapply(segmentData, 
#'                     FUN = function(x) {"state" %in% colnames(mcols(x))},
#'                     FUN.VALUE = logical(1)))) {
#'         stop("at least one sample doesn't have a metadata column ", 
#'                     "called \'state\'")
#'     }
#'     
#'     results <- list()
#'     
#'     for(type in c("AMPLIFICATION", "DELETION")) {
#'         
#'         dataTMP <- matrix(rep(NA, nb^2), nrow=nb)
#'         rownames(dataTMP) <- names
#'         colnames(dataTMP) <- names
#'         
#'         for(i in seq_len(nb)[-1]) {
#'             for(j in seq_len(i-1)) {
#'                 dataTMP[i, j] <- calculateOneOverlapMetric(
#'                     sample01=segmentData[[names[i]]], 
#'                     sample02=segmentData[[names[j]]],
#'                     method=method, type=type)
#'             }
#'         }
#'         results[[type]] <- dataTMP
#'     }
#'     
#'     # Return a list marked as an CNVMetric class containing:
#'     # 1- the metric results for the amplified regions
#'     # 2- the metric results for the deleted regions
#'     class(results) <- "CNVMetric"
#'     attr(results, 'metric') <- method
#'     
#'     return(results)
#' }


#' #' @title Plot metrics based on overlapping amplified/deleted regions
#' #' 
#' #' @description Plot one heatmap (or two heatmaps) of the metrics based on 
#' #' overlapping amplified/deleted regions. The user can select to print the
#' #' heatmap related to amplified, deleted regions or both.
#' #' 
#' #' @param metric a \code{CNVMetric} object containing the metrics calculated
#' #' by \code{calculateOverlapMetric}.
#' #' 
#' #' @param type a \code{character} string indicating which graph to generate. 
#' #' This should be (an unambiguous abbreviation of) one of "\code{BOTH}", 
#' #' "\code{AMPLIFICATION}" or "\code{DELETION}". Default: "\code{BOTH}".
#' #' 
#' #' @param colorRange a \code{vector} of 2 \code{character} string 
#' #' representing the 2 colors that will be
#' #' assigned to the lowest (0) and highest value (1) in the heatmap. 
#' #' Default: \code{c("white", "darkblue")}.
#' #' 
#' #' @param show_colnames a \code{boolean} specifying if column names are 
#' #' be shown. Default: \code{FALSE}.
#' #' 
#' #' @param silent a \code{boolean} specifying if the plot should not be drawn. 
#' #' Default: \code{TRUE}.
#' #' 
#' #' @param \ldots further arguments passed to 
#' #' \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}} method. Beware that
#' #' the \code{filename} argument cannot be used when \code{type} is 
#' #'  "\code{BOTH}".
#' #' 
#' #' @return a \code{gtable} object containing the heatmap(s) of the specified 
#' #' metric(s).
#' #' 
#' #' @examples
#' #' 
#' #' ## Load required package to generate the samples
#' #' require(GenomicRanges)
#' #' 
#' #' ## Create a GRangesList object with 3 samples
#' #' ## The stand of the regions doesn't affect the calculation of the metric
#' #' demo <- GRangesList()
#' #' demo[["sample01"]] <- GRanges(seqnames = "chr1", 
#' #'     ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
#' #'     end = c(2004603, 4577608, 31695808)), strand =  "*",
#' #'     state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#' #' 
#' #' demo[["sample02"]] <- GRanges(seqnames = "chr1", 
#' #'     ranges =  IRanges(start = c(1995066, 31611222, 31690000), 
#' #'     end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
#' #'     state = c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#' #' 
#' #' ## The amplified region in sample03 is a subset of the amplified regions 
#' #' ## in sample01
#' #' demo[["sample03"]] <- GRanges(seqnames = "chr1", 
#' #'     ranges =  IRanges(start = c(1906069, 4558838), 
#' #'     end = c(1909505, 4570601)), strand =  "*",
#' #'     state = c("AMPLIFICATION", "DELETION"))
#' #' 
#' #' ## Calculating Sorensen metric
#' #' metric <- calculateOverlapMetric(demo, method="sorensen")
#' #' 
#' #' ## Plot both amplification and deletion metrics
#' #' plotOverlapMetric(metric, type="BOTH")
#' #' 
#' #' 
#' #' ## Extra parameters, used by pheatmap(), can also be passed to the function
#' #' ## Here, we have the metric values print to the cell while the 
#' #' ## row names and column names are removed
#' #' plotOverlapMetric(metric, type="DELETION", show_rownames=FALSE,
#' #'     show_colnames=FALSE, main="deletion", display_numbers=TRUE,
#' #'     number_format="%.2f")
#' #' 
#' #' @seealso 
#' #' 
#' #' The default method  \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}}.
#' #' 
#' #' @author Astrid Deschênes
#' #' @importFrom grDevices colorRampPalette col2rgb 
#' #' @importFrom pheatmap pheatmap
#' #' @importFrom gridExtra grid.arrange arrangeGrob
#' #' @importFrom methods hasArg
#' #' @import GenomicRanges
#' #' @encoding UTF-8
#' #' @export
#' plotOverlapMetric <- function(metric, 
#'                                 type=c("BOTH", "AMPLIFICATION", "DELETION"),
#'                                 colorRange=c(c("white", "darkblue")), 
#'                                 show_colnames=FALSE, silent=TRUE, ...) {
#'     
#'     ## Validate that the metric parameter is a CNVMetric object
#'     if (!is.CNVMetric(metric)) {
#'         stop("\'metric\' must be a CNVMetric object.")
#'     }
#'     
#'     ## Assign type parameter
#'     type <- match.arg(type)
#'     
#'     ## Validate that the filename argument is not used when
#'     ## type "BOTH" is selected
#'     if (type == "BOTH" &&  hasArg("filename")) {
#'         stop("\'type\' cannot be \'BOTH\' when filename argument is used.")
#'     }
#'     
#'     ## Validate that the color name has only one value
#'     if (!is.character(colorRange) || length(colorRange) < 2) {
#'         stop("\'colorRange\' must be a vector of 2 color names.")
#'     }
#'     
#'     ## Validate that the color name is valid
#'     tryCatch(col2rgb(colorRange), error = function(e) {
#'         stop("\'colorRange\' must be be a vector of 2 valid color names.")
#'     })
#'     
#'     ## Extract the type of metric
#'     metricInfo <- attributes(metric)$metric
#'     
#'     plot_list <- list()
#'     
#'     ## Amplification
#'     if (type %in% c("AMPLIFICATION", "BOTH")) {
#'         plot_list[["AMPLIFICATION"]] <- plotOneOverlapMetric(metric=metric,
#'                                             type="AMPLIFICATION", 
#'                                             colorRange=colorRange, 
#'                                             show_colnames=show_colnames, 
#'                                             silent=silent, ...)  
#'     }
#'     
#'     ## Deletion
#'     if (type %in% c("DELETION", "BOTH")) {
#'         plot_list[["DELETION"]] <- plotOneOverlapMetric(metric=metric,
#'                                             type="DELETION", 
#'                                             colorRange=colorRange, 
#'                                             show_colnames=show_colnames, 
#'                                             silent=silent, ...)  
#'     }
#'     
#'     n_col <- ifelse(type == "BOTH", 2, 1)
#'     
#'     grid.arrange(arrangeGrob(grobs=plot_list, ncol=n_col))
#' }