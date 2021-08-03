#' #' @title Read segment files and their values and exclusion file to generate 
#' #' a dataset that can be used to calculate the metrics.
#' #' 
#' #' @description The segments and their values are extracted from each file, 
#' #' with a "SEG" or "seg" extension, present in the 
#' #' specified directory. Each segment is associated to its original file by
#' #' having a source metadata field assigned. When specified, the segments
#' #' for a BED file are loaded and used to exclude regions from further analysis.
#' #' All segments are gathered together
#' #' and a disjoin operation is done to create a collection of non-overlapping 
#' #' ranges. Those ranges are going to be used to calculate the metrics.
#' #' 
#' #' @param segDirectory a \code{character} string, the path to the directory
#' #' containing the segment files to compare. Only files with the ".seg" (or 
#' #' ".SEG") extension will be used. At least 2 segment files are needed to
#' #' be able to used the metrics.
#' #' 
#' #' @param chrInfo a \code{Seqinfo} containing the name and the length of the
#' #' chromosomes to analyze. Only the chromosomes contained in this
#' #' \code{Seqinfo} will be analyzed.
#' #' 
#' #' @param bedExclusionFile a \code{character} string or \code{NULL}, the path 
#' #' to the BED file that contains the regions that must be excluded from the
#' #' analysis. If \code{NULL}, it is assumed that there is no region of
#' #' exclusion for the calculation of the metrics. Default: \code{NULL}.
#' #' 
#' #' @param segmentWithHeader a \code{logical}, when \code{TRUE}, the
#' #' segment files have all a header that should not be imported. 
#' #' Default: \code{FALSE}.
#' #' 
#' #' @return a \code{list} marked as a \code{preMetricSegments} \code{class}  
#' #' containing a \code{GRanges} with the segment information from all 
#' #' segment files.
#' #'
#' #' @details 
#' #' 
#' #' The supported BED format requires that the first three columns are: 
#' #' \itemize{
#' #' \item { 1. chrom - The name of the chromosome.}
#' #' \item { 2. chromStart - The starting position of the feature in the 
#' #' chromosome or scaffold. The first base in a chromosome is numbered 0.}
#' #' \item { 3. chromEnd - The ending position of the feature in the chromosome 
#' #' or scaffold. The chromEnd base is not included in the display of the 
#' #' feature. }}
#' #' 
#' #' @examples
#' #'
#' #' # Path to the exclusion BED file
#' #' bedFile <- system.file("extdata", "mm10.blacklist.bed", 
#' #'     package = "CNVMetrics")
#' #' 
#' #' ## Prepare information about the chromosomes
#' #' ## The files are from mouse samples, chr1 to chr4
#' #' require(rtracklayer)
#' #' 
#' #' ## Get the information for mouse genome mm10
#' #' ## Limit the information to chromosomes 1 to 4
#' #' mouseInfo <- SeqinfoForUCSCGenome("mm10")  
#' #' mouseInfo <- mouseInfo[c("chr1", "chr2", "chr3", "chr4"),]
#' #'             
#' #' ## Path to the directory containing the segmentation files
#' #' segDir <- system.file("extdata", package="CNVMetrics")
#' #' 
#' #' ## TODO
#' #' 
#' #' prepareInformation(segDirectory = segDir, chrInfo = mouseInfo,
#' #'     bedExclusionFile = bedFile, segmentWithHeader = TRUE)
#' #' 
#' #' @author Astrid Deschenes, Pascal Belleau
#' #' @importFrom rtracklayer import
#' #' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' #' @importFrom methods is
#' #' @export
#' prepareInformation <- function(segDirectory, chrInfo, bedExclusionFile = NULL,
#'                                 segmentWithHeader = FALSE) {
#'     
#'     ## Validate that the chrInfo is a Seqinfo object
#'     if (!is(chrInfo, "Seqinfo")) {
#'         stop("chrInfo must be a Seqinfo object.")
#'     }
#'     
#'     ## Validate that the segmentWithHeader is a logical 
#'     if (!is.logical(segmentWithHeader)) {
#'         stop("segmentWithHeader must be a logical.")
#'     }
#'     
#'     ## Remove ending slash when present
#'     if (is.character(segDirectory) && length(segDirectory) > 1 &&
#'             endsWith(segDirectory, "/")) {
#'         segDirectory <- segDirectory[seq_len(length(segDirectory)-1)]
#'     }
#'     
#'     filesList <- list.files(path = segDirectory, pattern = ".seg", 
#'                             all.files = FALSE,
#'                             full.names = FALSE, recursive = FALSE,
#'                     ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE)
#'     
#'     ## Validate that the directory contains at least one segment file
#'     if (length(filesList) == 0) {
#'         stop("There is not segment file (seg or SEG extension) in ", 
#'                 "the segDirectory.")
#'     }
#'     
#'     ## Validate that the directory contains at least one segment file
#'     if (length(filesList) < 2) {
#'         stop("At least 2 segment files (seg or SEG extension) are ", 
#'                     "needed in the segDirectory.")
#'     }
#'     
#'     ## Read BED file and keep only segments in selected chromosomes
#'     excludedRegions <- NULL
#'     if (!is.null(bedExclusionFile)) {
#'         excludedRegions <- import(bedExclusionFile, format="BED")
#'         
#'         # Extract list of chromosomes to keep
#'         seqToKeep <- seqlevels(excludedRegions)[seqlevels(excludedRegions) %in% 
#'                                                 seqlevels(chrInfo)]
#'         
#'         if (length(seqToKeep) > 0) {
#'             excludedRegions <- keepSeqlevels(excludedRegions, seqToKeep, 
#'                                 pruning.mode = "coarse")
#'         } else {
#'             excludedRegions <- NULL
#'         }
#'     }
#'     
#'     ## Read segment files
#'     segFiles <- list()
#'     sources <- list()
#'     for (position in seq(1, length(filesList))) {
#'         # Get file name
#'         segFile <- filesList[position]
#'         # Remove extension from file name
#'         segFileShort <- substr(segFile, 1, nchar(segFile) - 4)
#'     
#'         # Get file path
#'         segPath <- paste0(segDirectory, "/", segFile)
#'         
#'         # Read segments
#'         tempRanges <- readSEGFile(segPath, uniqueTag = segFileShort, 
#'                                         header = segmentWithHeader)
#'         
#'         # Extract list of chromosomes to keep
#'         seqToKeep <- seqlevels(tempRanges)[seqlevels(tempRanges) %in% 
#'                                                     seqlevels(chrInfo)]
#'         
#'         if (length(seqToKeep) > 0) {
#'             # Keep only segments in selected chromosomes
#'             tempRanges <- keepSeqlevels(tempRanges, seqToKeep, 
#'                                     pruning.mode = "coarse")
#'             
#'             segFiles[[position]] <- tempRanges
#'             sources[[position]] <- segFileShort
#'         }  
#'     }
#'     
#'     result <- list()
#'     
#'     # Create segments using disjoin
#'     result$segments <- createSegments(segFiles, sources, excludedRegions)
#'     
#'     # Do regression using the first file as the dependent value
#'     result <- doRegression(result)
#'     
#'     result <- calculateRegressedValues(result)
#'     
#'     class(result) <- "preMetricSegments"
#'     
#'     return(result)
#' }



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
#' This should be (an unambiguous abbreviation of) one of "sorensen", 
#' "szymkiewicz" or "jaccard". Default: "sorensen".
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
#' \code{jaccard}:
#' 
#' This metric is calculated by dividing the size of the intersection 
#' by the size of the union of the two sets. With this metric, an overlap 
#' metric value of 1 is only obtained when the two samples are identical. 
#' 
#' @return an object of class "\code{CNVMetric}" which contains the calculated
#' metric. This object is a list with the following components:
#' \itemize{
#' \item{\code{AMPLIFICATION}}{ a lower-triangular \code{matrix} with the 
#'     results of the selected metric on the amplified regions for each paired
#'     samples. The value \code{NA} is present when the metric cannot be 
#'     calculated. The value \code{NA} is also present in the top-triangular 
#'     section, as well as the diagonal, of the matrix.
#' }
#' \item{\code{DELETION}}{ a lower-triangular \code{matrix} with the 
#'     results of the selected metric on the deleted regions for each paired
#'     samples. The value \code{NA} is present when the metric cannot be 
#'     calculated. The value \code{NA} is also present in the top-triangular 
#'     section, as well as the diagonal, of the matrix.
#' }}
#' 
#' The object has the following attributes (besides "class" equal 
#' to "CNVMetric"):
#' \itemize{
#' \item{\code{metric}}{ the metric used for the calculation.
#' } 
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
#' Jaccard, P. (1912), The Distribution of the Flora in the Alpine Zone.  
#' New Phytologist, 11: 37-50. 
#' doi: \url{https://doi.org/10.1111/j.1469-8137.1912.tb05611.x}
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
#' calculateOverlapMetric(demo, method="sorensen")
#' 
#' ## Calculating Szymkiewicz-Simpson metric
#' calculateOverlapMetric(demo, method="szymkiewicz")
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges 
#' @encoding UTF-8
#' @export
calculateOverlapMetric <- function(segmentData, 
                                   method=c("sorensen", "szymkiewicz", 
                                            "jaccard")) {
    
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
        stop("at least one sample doesn't have a metadata column ", 
             "called \'state\'")
    }
    
    results <- list()
    
    for(type in c("AMPLIFICATION", "DELETION")) {
        
        dataTMP <- matrix(rep(NA, nb^2), nrow=nb)
        rownames(dataTMP) <- names
        colnames(dataTMP) <- names
        
        for(i in seq_len(nb)[-1]) {
            for(j in seq_len(i-1)) {
                dataTMP[i, j] <- calculateOneOverlapMetric(
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
#' called '\code{log2ratio}' with the log2ratio values.
#' 
#' @param method a \code{character} string representing the metric to be used. 
#' This should be (an unambiguous abbreviation of) one of 
#' "weightedEuclideanDistance". Default: "weightedEuclideanDistance".
#' 
#' @param minThreshold a single \code{numeric} setting the minimum value 
#' to consider two segments as different during the metric calculation. If the 
#' absolute difference is below or equal to threshold, the difference will be 
#' replaced by zero. Default: 0.2.
#'  
#' 
#' @param excludedRegions an optional \code{GRanges} containing the regions 
#' that have to be excluded for the metric calculation. Default: \code{NULL}.
#' 
#' @details 
#' 
#' The weighted euclidean distance is 
#' \eqn{(\sum((x_i - y_i)^2 * log(nbrBases_i))^0.5} 
#' where \code{x} and \code{y} are the
#' values of 2 samples for a specific segment \code{i} and \code{nbrBases} the 
#' number of bases of the segment \code{i}.
#' 
#' @return an object of class "\code{CNVMetric}" which contains the calculated
#' metric. This object is a list with the following components:
#' \itemize{
#' \item{\code{LOG2RATIO}}{ a lower-triangular \code{matrix} with the 
#'     results of the selected metric on the log2ratio values for each paired
#'     samples. The value \code{NA} is present when the metric cannot be 
#'     calculated. The value \code{NA} is also present in the top-triangular 
#'     section, as well as the diagonal, of the matrix.
#' }}
#' 
#' The object has the following attributes (besides "class" equal 
#' to "CNVMetric"):
#' \itemize{
#' \item{\code{metric}}{ the metric used for the calculation.
#' } 
#' \item{\code{names}}{ the names of the two matrix containing the metrics for
#' the amplified and deleted regions.
#' }}         
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
#'     log2ratio = c(2.5555, 1.9932, -0.9999))
#' 
#' demo[["sample02"]] <- GRanges(seqnames = "chr1", 
#'     ranges =  IRanges(start = c(1995066, 31611222, 31690000), 
#'     end = c(2204505, 31689898, 31895666)), strand =  c("-", "+", "+"),
#'     log2ratio = c(0.3422, 0.5454, -1.4444))
#' 
#' ## The amplified region in sample03 is a subset of the amplified regions 
#' ## in sample01
#' demo[["sample03"]] <- GRanges(seqnames = "chr1", 
#'     ranges =  IRanges(start = c(1906069, 4558838), 
#'     end = c(1909505, 4570601)), strand =  "*",
#'     log2ratio = c(3.2222, -1.3232))
#' 
#' ## Calculating Sorensen metric
#' calculateLog2ratioMetric(demo, method="weightedEuclideanDistance")
#' 
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges 
#' @encoding UTF-8
#' @export
calculateLog2ratioMetric <- function(segmentData, 
                                   method=c("weightedEuclideanDistance"),
                                   minThreshold=0.2, excludedRegions=NULL) {
    
    method <- match.arg(method)
    
    names <- names(segmentData)
    nb <- length(names)
    
    ## At least 2 samples must be present
    if (nb < 2) {
        stop("At least 2 samples must be present in segmentData")
    }
    
    ## All samples must have a metadata column called 'log2ratio' with
    ## log2ratio values
    if (!all(vapply(segmentData, 
                    FUN = function(x) {"log2ratio" %in% colnames(mcols(x))},
                    FUN.VALUE = logical(1)))) {
        stop("at least one sample doesn't have a metadata column ", 
             "called \'log2ratio\'")
    }
    
    results <- list()
    
    for(type in c("LOG2RATIO")) {
        
        dataTMP <- matrix(rep(NA, nb^2), nrow=nb)
        rownames(dataTMP) <- names
        colnames(dataTMP) <- names
        
        for(i in seq_len(nb)[-1]) {
            for(j in seq_len(i-1)) {
                
                dataTMP[i, j] <- calculateOneLog2valueMetric(
                    sample01=segmentData[[names[i]]], 
                    sample02=segmentData[[names[j]]],
                    method=method, minThreshold=minThreshold, 
                    bedExclusion=excludedRegions)
            }
        }
        results[[type]] <- dataTMP
    }
    
    # Return a list marked as an CNVMetric class containing:
    # 1- the metric results for the log2ratio 
    class(results) <- "CNVMetric"
    attr(results, 'metric') <- method
    
    return(results)
}



#' @title Plot metrics present in a \code{CNVMetric} object.
#' 
#' @description Plot one heatmap (or two heatmaps) of the metrics present in
#' a \code{CNVMetric} object. For the overlapping metrics, the user can select 
#' to print the heatmap related to amplified or deleted regions or both. The
#' \code{NA} values present in the metric matrix are transformed into zero for
#' the creation of the heatmap.
#' 
#' @param metric a \code{CNVMetric} object containing the metrics calculated
#' by \code{calculateOverlapMetric} or by \code{calculateLog2ratioMetric}.
#' 
#' @param type a \code{character} string indicating which graph to generate. 
#' This should be (an unambiguous abbreviation of) one of "\code{ALL}", 
#' "\code{AMPLIFICATION}" or "\code{DELETION}". This is useful for the 
#' overlapping metrics that have both the amplified and deleted metrics in the
#' \code{CNVMetric} object. Default: "\code{ALL}".
#' 
#' @param colorRange a \code{vector} of 2 \code{character} string 
#' representing the 2 colors that will be
#' assigned to the lowest (0) and highest value (1) in the heatmap. 
#' Default: \code{c("white", "darkblue")}.
#' 
#' @param show_colnames a \code{boolean} specifying if column names are 
#' be shown. Default: \code{FALSE}.
#' 
#' @param silent a \code{boolean} specifying if the plot should not be drawn. 
#' Default: \code{TRUE}.
#' 
#' @param \ldots further arguments passed to 
#' \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}} method. Beware that
#' the \code{filename} argument cannot be used when \code{type} is 
#'  "\code{BOTH}".
#' 
#' @return a \code{gtable} object containing the heatmap(s) of the specified 
#' metric(s).
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
#' ## Plot both amplification and deletion metrics
#' plotMetric(metric, type="ALL")
#' 
#' ## Extra parameters, used by pheatmap(), can also be passed to the function
#' ## Here, we have the metric values print to the cell while the 
#' ## row names and column names are removed
#' plotMetric(metric, type="DELETION", show_rownames=FALSE,
#'     show_colnames=FALSE, main="deletion", display_numbers=TRUE,
#'     number_format="%.2f")
#' 
#' @seealso 
#' 
#' The default method  \code{\link[pheatmap:pheatmap]{pheatmap::pheatmap()}}.
#' 
#' @author Astrid Deschênes
#' @importFrom grDevices colorRampPalette col2rgb 
#' @importFrom pheatmap pheatmap
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom methods hasArg
#' @import GenomicRanges
#' @encoding UTF-8
#' @export
plotMetric <- function(metric, 
                              type=c("ALL", "AMPLIFICATION", "DELETION"),
                              colorRange=c(c("white", "darkblue")), 
                              show_colnames=FALSE, silent=TRUE, ...) {
    
    ## Validate that the metric parameter is a CNVMetric object
    if (!is.CNVMetric(metric)) {
        stop("\'metric\' must be a CNVMetric object.")
    }
    
    ## Assign type parameter
    type <- match.arg(type)
    
    ## Validate that the filename argument is not used when
    ## type "ALL" is selected
    if (type == "ALL" &&  hasArg("filename")) {
        stop("\'type\' cannot be \'ALL\' when filename argument is used.")
    }
    
    ## Validate that the type argument is not "AMPLIFICATION"
    ## when the metric does not contain this type of metric
    if (type == "AMPLIFICATION" &&  ! "AMPLIFICATION" %in% names(metric)) {
        stop("\'type\' cannot be \'AMPLIFICATION\' for this metric object.")
    }
    
    ## Validate that the type argument is not "DELETION"
    ## when the metric does not contain this type of metric
    if (type == "DELETION" &&  ! "DELETION" %in% names(metric)) {
        stop("\'type\' cannot be \'DELETION\' for this metric object.")
    }
    
    ## Validate that the color name has only one value
    if (!is.character(colorRange) || length(colorRange) < 2) {
        stop("\'colorRange\' must be a vector of 2 color names.")
    }
    
    ## Validate that the color name is valid
    tryCatch(col2rgb(colorRange), error = function(e) {
        stop("\'colorRange\' must be be a vector of 2 valid color names.")
    })
    
    ## Extract the type of metric
    metricInfo <- attributes(metric)$metric
    
    plot_list <- list()
    
    nameList <- ifelse(type %in% c("AMPLIFICATION", "DELETION"), type, 
                                                            names(metric))
    
    for (name in nameList) {
        plot_list[[name]] <- plotOneMetric(metric=metric,
                                                type=name, 
                                                colorRange=colorRange, 
                                                show_colnames=show_colnames, 
                                                silent=silent, ...)  
    }
    
    
    n_col <- length(plot_list)
    
    grid.arrange(arrangeGrob(grobs=plot_list, ncol=n_col))
}



