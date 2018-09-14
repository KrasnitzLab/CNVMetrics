
#' @title TODO
#' 
#' @description TODO
#' 
#' @param fileList a \code{list} of \code{GRanges}, the segments from multiple
#' files.
#' 
#' @param sourceList a \code{list}
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
#' @importFrom GenomicRanges disjoin findOverlaps elementMetadata
#' @importFrom S4Vectors queryHits subjectHits values<-
#' @importFrom magrittr %>%
#' @keywords internal
createSegments <- function(fileList, sourceList, bedExclusion) {
    
    fileData <- fileList
    
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
    
    for (i in 1:length(fileData)) {
        if (!is.null(fileData[[i]])) {
            olaps <- findOverlaps(results, fileData[[i]])
            temp <- elementMetadata(results)
            temp[, sourceList[[i]]] <- NA
            temp[queryHits(olaps), sourceList[[i]]] <- 
                fileData[[i]]$score[subjectHits(olaps)]
            values(results) <- temp
        }
    }
    
    return(results)
}