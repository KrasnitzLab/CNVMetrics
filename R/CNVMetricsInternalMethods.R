
#' @title TODO
#' 
#' @description TODO
#' 
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#'
#' @param bedExclusion a \code{GRanges}, the regions that must be
#' excluded from the analysis. Default: \code{NULL}.
#' 
#' @return a \code{GRanges} containing the segment information from the file.
#'
#' @examples
#'
#' # TODO
#' 
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom GenomicRanges disjoin findOverlaps
#' @importFrom S4Vectors queryHits
#' @importFrom magrittr %>%
#' @keywords internal
createSegments <- function(fileList, bedExclusion) {
    
    ## Add same columns to GRanges related to bed exclusion than the
    ## other GRanges
    if (!is.null(bedExclusion) && (length(bedExclusion) > 0)) {
        bedExclusion$score <- NA
        bedExclusion$source <- "exclusion"
        fileList[[length(fileList)+1]] <- bedExclusion   
    }
    
    results <- do.call("c", fileList) %>% disjoin()
    results$included <- TRUE
    
    ## Add information about excluded regions
    if (!is.null(bedExclusion) && (length(bedExclusion) > 0)) {
        olaps <- findOverlaps(results, bedExclusion)
        
        if (length(olaps) > 0) {
            results[queryHits(olaps)]$included <- FALSE
        }
    }
    
    return(results)
}