#' @title Read segment files and exclusion file to generate a dataset
#' that can be used to calculate the metrics.
#' 
#' @description TODO
#' 
#' @param segDirectory a \code{character} string, the path to the directory
#' containing the segment files to compared. Only files with the ".seg" (or 
#' ".SEG") extension will be used. At least 2 segment files are needed to
#' be able to used the metrics.
#' 
#' @param chrInfo a \code{Seqinfo} containing the name and the length of the
#' chromosomes to analyze. Only the chomosomes contained in this
#' \code{Seqinfo} will be analyzed.
#' 
#' @param segmentWithHeader a \code{logical}, when \code{TRUE}, the
#' segment files have all a header that should not be imported. 
#' Default: \code{FALSE}.
#'
#' @param bedExclusionFile a \code{character} string or \code{NULL}, the path 
#' to the BED file that contains the regions that must be excluded from the
#' analysis. If \code{NULL}, it is assume that there is no region of
#' exclusion for the calculation of the metrics. Default: \code{NULL}.
#' 
#' @return a \code{GRanges} containing the segment information from the file.
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
#' # TODO
#' 
#' @author Astrid Deschenes, Pascal Belleau
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
    
    ## Read BED file
    excludedRegions <- NULL
    if (!is.null(bedExclusionFile)) {
        # Read excluded segments
        excludedRegions <- import(bedExclusionFile, format="BED")
        
        # Keep only segments in selected chromosomes
        excludedRegions <- keepSeqlevels(excludedRegions, seqlevels(chrInfo), 
                                pruning.mode = "coarse")
    }
    
    ## Read segment files
    segFiles <- list()
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
        
        # Keep only segments in selected chromosomes
        tempRanges <- keepSeqlevels(tempRanges, seqlevels(chrInfo), 
                                    pruning.mode = "coarse")
        
        segFiles[[position]] <- tempRanges
    }
    
    result <- createSegments(segFiles, excludedRegions)
    
    return(result)
}