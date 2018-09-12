
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
#' @importFrom GenomicRanges disjoin 
#' @importFrom magrittr %>%
#' @internal
createSegments <- function(fileList, bedExclusion) {
    
    if (!is.null(bedExclusion)) {
        bedExclusion$source <- "exclusion"
        fileList[[length(fileList)+1]] <- bedExclusion   
    }
    
    results <- do.call("c", fileList) %>% disjoin()
    
    return(results)
}