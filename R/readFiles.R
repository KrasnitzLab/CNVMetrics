#' @title Read segment file and return a \code{GRanges}
#' 
#' @description Read a segment file and extract the relevant
#' information to generate a \code{GRanges}.
#' 
#' @param filepath a \code{character} string, the path to the file to 
#' read. Alternatively can be a connection.
#'
#' @param uniqueTag a \code{character} string, the unique name assigned to
#' the segment data associated to this file. Each imported file should
#' have a unique tag.
#'
#' @param header a \code{logical}, when \code{TRUE}, the
#' file has a header that should not be imported. Default: \code{FALSE}.
#'
#' @return a \code{GRanges} containing the segment information from the file.
#'
#' @examples
#'
#' # TODO
#' segFile <- system.file("inst/extdata/dem01.seg", 
#'     package = "CNVMetrics")
#' 
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom utils read.table
#' @export
readSEGFile <- function(filepath, uniqueTag, header=FALSE) {
    
    if (!is.character(uniqueTag)) {
        stop("uniqueTag must be a character string.")
    }
    
    if (!is.logical(header)) {
        stop("header must be a logical value.")
    }
    
    data <- read.table(filepath, header=header)
    
    if (nrow(data) < 6) {
        stop("The segment file must have at least 6 columns.")
    }
    
    dataRanges <- GRanges(seqnames = data[, 2], 
                            ranges=IRanges(start=data[, 3], end=data[, 4]), 
                            score=data[, 6], source=uniqueTag)
    
    return(dataRanges)
}