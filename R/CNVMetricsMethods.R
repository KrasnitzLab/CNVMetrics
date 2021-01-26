#' @title Read segment files and their values and exclusion file to generate 
#' a dataset that can be used to calculate the metrics.
#' 
#' @description The segments and their values are extracted from each file, 
#' with a "SEG" or "seg" extension, present in the 
#' specified directory. Each segment is associated to its original file by
#' having a source metadata field assigned. When specified, the segments
#' for a BED file are loaded and used to exclude regions from further analysis.
#' All segments are gathered together
#' and a disjoin operation is done to create a collection of non-overlapping 
#' ranges. Those ranges are going to be used to calculate the metrics.
#' 
#' @param segDirectory a \code{character} string, the path to the directory
#' containing the segment files to compare. Only files with the ".seg" (or 
#' ".SEG") extension will be used. At least 2 segment files are needed to
#' be able to used the metrics.
#' 
#' @param chrInfo a \code{Seqinfo} containing the name and the length of the
#' chromosomes to analyze. Only the chromosomes contained in this
#' \code{Seqinfo} will be analyzed.
#' 
#' @param bedExclusionFile a \code{character} string or \code{NULL}, the path 
#' to the BED file that contains the regions that must be excluded from the
#' analysis. If \code{NULL}, it is assumed that there is no region of
#' exclusion for the calculation of the metrics. Default: \code{NULL}.
#' 
#' @param segmentWithHeader a \code{logical}, when \code{TRUE}, the
#' segment files have all a header that should not be imported. 
#' Default: \code{FALSE}.
#' 
#' @return a \code{list} marked as a \code{preMetricSegments} \code{class}  
#' containing a \code{GRanges} with the segment information from all 
#' segment files.
#'
#' @details 
#' 
#' The supported BED format requires that the first three columns are: 
#' \itemize{
#' \item { 1. chrom - The name of the chromosome.}
#' \item { 2. chromStart - The starting position of the feature in the 
#' chromosome or scaffold. The first base in a chromosome is numbered 0.}
#' \item { 3. chromEnd - The ending position of the feature in the chromosome 
#' or scaffold. The chromEnd base is not included in the display of the 
#' feature. }}
#' 
#' @examples
#'
#' # Path to the exclusion BED file
#' bedFile <- system.file("inst/extdata/exclusion_hg38_demo_chr1_and_2.bed", 
#'     package = "CNVMetrics")
#'     
#' # Path to the directory containing the segmentation files
#' segDir <- system.file("inst/extdata", package="CNVMetrics")
#' 
#' # TODO
#' 
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb keepSeqlevels seqlevels
#' @importFrom methods is
#' @export
prepareInformation <- function(segDirectory, chrInfo, bedExclusionFile = NULL,
                               segmentWithHeader=FALSE) {
    
    ## Validate that the chrInfo is a Seqinfo object
    if (!is(chrInfo, "Seqinfo")) {
        stop("chrInfo must be a Seqinfo object.")
    }
    
    ## Validate that the segmentWithHeader is a logical 
    if (!is.logical(segmentWithHeader)) {
        stop("segmentWithHeader must be a logical.")
    }
    
    ## Remove ending slash when present
    if (is.character(segDirectory) && length(segDirectory) > 1 &&
            endsWith(segDirectory, "/")) {
        segDirectory <- segDirectory[1:length(segDirectory)-1]
    }
    
    filesList <- list.files(path = segDirectory, pattern = ".seg", 
                            all.files = FALSE,
                            full.names = FALSE, recursive = FALSE,
               ignore.case = TRUE, include.dirs = FALSE, no.. = FALSE)
    
    ## Validate that the directory contains at least one segment file
    if (length(filesList) == 0) {
        stop(paste0("There is not segment file (seg or SEG extension) in ", 
             "the segDirectory."))
    }
    
    ## Validate that the directory contains at least one segment file
    if (length(filesList) < 2) {
        stop(paste0("At least 2 segment files (seg or SEG extension) are ", 
                    "needed in the segDirectory."))
    }
    
    ## Read BED file and keep only segments in selected chromosomes
    excludedRegions <- NULL
    if (!is.null(bedExclusionFile)) {
        excludedRegions <- import(bedExclusionFile, format="BED")
        
        # Extract list of chromosomes to keep
        seqToKeep <- seqlevels(excludedRegions)[seqlevels(excludedRegions) %in% 
                                               seqlevels(chrInfo)]
        
        if (length(seqToKeep) > 0) {
            excludedRegions <- keepSeqlevels(excludedRegions, seqToKeep, 
                                pruning.mode = "coarse")
        } else {
            excludedRegions <- NULL
        }
    }
    
    ## Read segment files
    segFiles <- list()
    sources <- list()
    for (position in seq(1, length(filesList))) {
        # Get file name
        segFile <- filesList[position]
        # Remove extension from file name
        segFileShort <- substr(segFile, 1, nchar(segFile) - 4)
    
        # Get file path
        segPath <- paste0(segDirectory, "/", segFile)
        
        # Read segments
        tempRanges <- readSEGFile(segPath, uniqueTag = segFileShort, 
                                     header = segmentWithHeader)
        
        # Extract list of chromosomes to keep
        seqToKeep <- seqlevels(tempRanges)[seqlevels(tempRanges) %in% 
                                                seqlevels(chrInfo)]
        
        if (length(seqToKeep) > 0) {
            # Keep only segments in selected chromosomes
            tempRanges <- keepSeqlevels(tempRanges, seqToKeep, 
                                    pruning.mode = "coarse")
            
            segFiles[[position]] <- tempRanges
            sources[[position]] <- segFileShort
        }  
    }
    
    result <- list()
    
    # Create segments using disjoin
    result$segments <- createSegments(segFiles, sources, excludedRegions)
    
    # Do regression using the first file as the dependent value
    result <- doRegression(result)
    
    result <- calculateRegressedValues(result)
    
    class(result) <- "preMetricSegments"
    
    return(result)
}


#' @title Calculate Weighted Euclidean distance metric between samples.
#' 
#' @description The weighted Euclidean distance metric corresponds to the
#' euclidean distance between 2 samples multiplied by the natural logarithm 
#' of the number of bases of the analyzed segment. The final metric is the 
#' squared sum of the values obtained for all segments that are not 
#' excluded of the analysis.
#' 
#' @param segmentData a \code{list} marked as a \code{preMetricSegments} 
#' \code{class} that contains the segment information from at least 2 segment 
#' files.
#' 
#' @param minThreshold a single \code{numeric} setting the minimum value 
#' to consider two segments as different durant the metric calculation. If the 
#' absolute difference is bellow or equal to threshold, the value will be 
#' replaced by zero. Default: 0.2.
#' 
#' @return a \code{matrix} TODO
#' 
#' @details 
#' 
#' The weighted euclidean distance is 
#' \eqn{(\sum((x_i - y_i)^2 * log(nbrBases_i))^0.5} 
#' where \code{x} and \code{y} are the
#' values of 2 samples for a specific segment \code{i} and \code{nbrBases} the 
#' number of bases of the segment \code{i}.
#' 
#' 
#' @examples
#'
#' # TODO
#' 
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @importFrom GenomicRanges elementMetadata
#' @importFrom IRanges ranges width
#' @export
calculateWeightedEuclideanDistance <- function(segmentData, minThreshold=0.2) {
    
    if (! is(segmentData, "preMetricSegments")) {
        stop("segmentData must be a list marked as preMetricSegments class.")
    }
    
    if (!is.numeric(minThreshold)) {
        stop("minThreshold must be a numerical value.")
    }
    
    segments <- segmentData$segments
    
    names <- colnames(elementMetadata(segments))
    names <- names[names != "included"]
    
    nbNames <- length(names)
    
    metric <- matrix(nrow = nbNames, ncol = nbNames, 
                        dimnames = rep(list(ID = names), 2))
    
    incWidth <- width(ranges(segments[segments$included, ]))
    incResults <- elementMetadata(segments[segments$included, ])
    
    for (i in 1:(nbNames-1)) {
        for (j in (i+1):nbNames) {
            temp01 <- incResults[, c(names[i])] - incResults[, c(names[j])]
            
            ## Set values to zero when lower than threshold
            tempPos <- which(abs(temp01) <= minThreshold) 
            if (length(tempPos) > 0) {
                temp01[tempPos] <- 0.0
            }
            
            temp01 <- temp01 * temp01 * log(incWidth)
            final <- sum(temp01, na.rm = TRUE) ^ (1/2)
            metric[names[i], names[j]] <- final
            metric[names[j], names[i]] <- final
        }
    }
 
    for (i in 1:nbNames) {
        metric[names[i], names[i]] <- 0.0
    }
    
    return(metric)   
}




#' @title Calculate metric using overlapping amplified/deleted regions
#' 
#' @description Calculate a specific metric using overlapping 
#' amplified/deleted regions between to samples. The metric is calculated for
#' the amplified and deleted regions separately. When more than 2 samples are 
#' present, the metric is calculated for each sample pair.
#' 
#' @param segmentData a \code{GRangesList} that contains the segment 
#' informations, including amplified/deleted status, from at least 2 samples. 
#' All samples must have a metadata column called 'state' with amplified 
#' regions identified as 'AMPLIFICATION' and deleted regions identified as
#' 'DELETION'; regions with different identifications will not be used in the
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
#' If \code{method} is "\code{sorensen}", the metric is calculated by dividing
#' twice the size of the intersection by the sum of the size of the two sets. 
#' With this metric, an overlap metric value of 1 is only obtained when the
#' two samples are identical. 
#' 
#' If \code{method} is "\code{szymkiewicz}", the metric is calculated by 
#' dividing the size of the intersection by the size of the smallest set. With
#' this metric, if one set is a subset of the other set, the overlap 
#' metric value is 1.
#' 
#' @return a \code{list} TODO
#' 
#' @references 
#' 
#' Sørensen, Thorvald. n.d. “A Method of Establishing Groups of Equal Amplitude 
#' in Plant Sociology Based on Similarity of Species and Its Application to 
#' Analyses of the Vegetation on Danish Commons.” Biologiske Skrifter, 
#' no. 5: 1–34.
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
    if (!all(sapply(segmentData, 
            FUN = function(x) {"state" %in% colnames(mcols(x))}))) {
        stop(paste0("at least one sample doesn't have a metadata column ", "
             called \'state\'"))
    }
    
    results <- list()
    
    for(type in c("AMPLIFICATION", "DELETION")) {
        
        dataTMP <- matrix(rep(NA, nb^2), nrow = nb)
        rownames(dataTMP) <- names
        colnames(dataTMP) <- names
        
        for(i in seq_len(nb)[-1]) {
            for(j in seq_len(i-1)) {
                dataTMP[i, j] <- calculateOverlapMetric(
                    sample01 = segmentData[[names[i]]], 
                    sample02 = segmentData[[names[j]]],
                    method=method, type=type)
            }
        }
        results[[type]] <- dataTMP
    }
        
    return(results)
}
