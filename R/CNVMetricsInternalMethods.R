
#' @title Generate common segments to enable calculation of metrics on 
#' multiple segment files.
#' 
#' @description All segments are gathered together, including exclusion 
#' segments when specified, and a disjoin operation is done to create a 
#' collection of non-overlapping ranges. The ranges included in the exclusion
#' segments are marked as so to be removed from futur analysis.
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


#' @title Linear regression analysis between each paired samples using
#' segment values.
#' 
#' @description Do an linear regression between each paired samples using
#' segment values only from the included segments. The segment file used as
#' reference is the same for all analysis.
#' 
#' @param segmentData a \code{list} of that 
#' contains the segments from multiple files. The \code{list} is composed of 
#' those entries:
#' \itemize{
#' \item a \code{segment} that contains the \code{GRanges} with the segment
#' information.
#' }
#' 
#' @return a \code{list} of that 
#' contains the segments from multiple files. The \code{list} is composed of 
#' those entries:
#' \itemize{
#' \item a \code{segment} that contains the \code{GRanges} with the segment
#' information.
#' \item a \code{regression} that contains the result of the paired 
#' regressions.
#' }
#' 
#' @examples
#'
#' # TODO
#' 
#' @author Astrid Deschenes, Pascal Belleau
#' @importFrom GenomicRanges elementMetadata
#' @importFrom stats lm
#' @keywords internal
doRegression <- function(segmentData) {
    
    segments <- segmentData$segments
    
    names <- colnames(elementMetadata(segments))
    names <- names[names != "included"]
    
    nbNames <- length(names)
    
    metric <- matrix(nrow = nbNames, ncol = nbNames, 
                     dimnames = rep(list(ID = names), 2))
    
    incResults <- elementMetadata(segments[segments$included, ])
    
    segmentData$regression <- list()
    
    for (i in 2:nbNames) {
        subData <- incResults[, c(names[1], names[i])]
        colnames(subData) <- c("y", "x")
        reg <- lm("y ~ x", data=subData)
        segmentData$regression[[i - 1]] <- list()
        segmentData$regression[[i - 1]][["lm"]] <- reg
        segmentData$regression[[i - 1]][["y_used"]] <- names[1]
        segmentData$regression[[i - 1]][["x_used"]] <- names[i]
    }
    
    return(segmentData)   
}


#' @title calculate the regressed values for all segment file.
#' 
#' @description Use the linear regression model obtained for each paired of
#' segment files (current file versus reference) to calculate the regressed 
#' values for all segment files except for the reference file.
#' 
#' @param segmentData a \code{list} of that 
#' contains the segments from multiple files. The \code{list} is composed of 
#' those entries:
#' \itemize{
#' \item a \code{segment} that contains the \code{GRanges} with the segment
#' information.
#' \item a \code{regression} that contains the result of the paired 
#' regressions. The number of entries corresponds to the number of paired
#' segment files considering that the reference file is always the same. Each 
#' entry is a \code{list} composed of those entries:
#' \itemize{
#' \item a \code{x_used} with the name of the file that segments values
#' are used as x values in the regression model.
#' \item a \code{y_used} with the name of the file that segments values
#' are used as y values in the regression model. This is the reference file. 
#' It is the same for all paired regressions.
#' \item a \code{lm} that contains the regression model.
#' }
#' }
#' 
#' @return a \code{list} of that 
#' contains the segments from multiple files. The \code{list} is composed of 
#' those entries:
#' \itemize{
#' \item a \code{segments} that contains the \code{GRanges} with the segment
#' information and values from each segment file.
#' \item a \code{regression} that contains the result of the paired 
#' regressions.
#' \item a \code{regressedData} that contains the \code{GRanges} with the 
#' segment information and the regressed values for each segment file.
#' }
#' 
#' @examples
#'
#' # TODO
#' 
#' @author Astrid Deschênes, Pascal Belleau
#' @importFrom GenomicRanges elementMetadata elementMetadata<-
#' @importFrom stats predict
#' @keywords internal
calculateRegressedValues <- function(segmentData) {
    
    segments <- elementMetadata(segmentData$segments)

    #segmentData$regressed <- segments
    
    for (i in 1:length(segmentData$regression)) {
        lmData <- segmentData$regression[[i]][["lm"]]
        xName <- segmentData$regression[[i]][["x_used"]]
        tempVal <- data.frame(x=segments[, xName])
        segments[, xName] <- as.vector(predict(lmData, newdata = tempVal))
    }
    
    segmentData$regressedData <- segmentData$segments
    elementMetadata(segmentData$regressedData) <- segments
    
    return(segmentData)
}




#' @title Calculate metric using overlapping amplified/deleted regions between 
#' two samples.
#' 
#' @description Calculate a specific metric using overlapping regions between
#' two samples. 
#' 
#' @param sample01 a \code{GRanges} TODO 
#' @param sample02 a \code{GRanges} TODO
#' @param method a string of \code{character} representing the metric to be
#' used.
#' @param type a string of \code{character} TODO
#' 
#' @return a \code{list} of TODO
#' 
#' @examples
#'
#' # TODO
#' 
#' @author Astrid Deschênes
#' @keywords internal
calculateOverlapMetric <- function(sample01, sample02, method, type) {
    
    sample01 <- sample01[sample01$state == type,]
    sample02 <- sample02[sample02$state == type,]
    
    result <- 0.0
    
    if (length(sample01) > 0 && length(sample02) > 0) { 
        result <- switch(method,
           sorensen = calculateSorensen(sample01, sample02),
           szymkiewicz = calculateSzymkiewicz(sample01, sample02))
    }
    
    return(result)
}
    

#' @title Calculate Sorensen metric using overlapping amplified/deleted 
#' regions between two samples.
#' 
#' @description Calculate Sorensen metric using overlapping regions between
#' two samples. 
#' 
#' @param sample01 a \code{GRanges} TODO 
#' @param sample02 a \code{GRanges} TODO
#' 
#' @return a \code{list} of TODO
#' 
#' @examples
#'
#' # TODO
#' 
#' require(GenomicRanges)
#'
#' sample01 <- GRanges(seqnames = "chr1", 
#'     ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
#'     end = c(2004603, 4577608, 31695808)), strand =  "*")
#' sample02 <- GRanges(seqnames = "chr1", 
#'     ranges =  IRanges(start = c(1995066, 31611222), 
#'     end = c(2204505, 31689898)), strand =  "*")
#'     
#' CNVMetrics:::calculateSorensen(sample01, sample02)
#'     
#' @author Astrid Deschênes
#' @importFrom GenomicRanges intersect width
#' @keywords internal
calculateSorensen <- function(sample01, sample02) {
    
    inter <- sum(width(intersect(sample01, sample02, ignore.strand=TRUE)))
    widthSample01 <- sum(width(sample01))
    widthSample02 <- sum(width(sample02))
    
    result <- ifelse((widthSample01 + widthSample02) > 0, 
                     (2.0 * inter)/(widthSample01 + widthSample02),
                     NA)
    return(result)
}
    
#' @title Calculate Szymkiewicz-Simpson metric using overlapping 
#' amplified/deleted regions between two samples.
#' 
#' @description Calculate Szymkiewicz-Simpson metric using overlapping 
#' regions between two samples. 
#' 
#' @param sample01 a \code{GRanges} TODO 
#' @param sample02 a \code{GRanges} TODO
#' 
#' @return a \code{list} of TODO
#' 
#' @examples
#'
#' # TODO
#' 
#' require(GenomicRanges)
#'
#' sample01 <- GRanges(seqnames = "chr1", 
#'     ranges =  IRanges(start = c(1905048, 4554832, 31686841), 
#'     end = c(2004603, 4577608, 31695808)), strand =  "*")
#' sample02 <- GRanges(seqnames = "chr1", 
#'     ranges =  IRanges(start = c(1995066, 31611222), 
#'     end = c(2204505, 31689898)), strand =  c("+", "-"))
#'     
#' CNVMetrics:::calculateSzymkiewicz(sample01, sample02)
#'     
#' @author Astrid Deschênes
#' @importFrom GenomicRanges intersect width
#' @keywords internal
calculateSzymkiewicz <- function(sample01, sample02) {
    
    inter <- sum(width(intersect(sample01, sample02, ignore.strand=TRUE)))
    widthSample01 <- sum(width(sample01))
    widthSample02 <- sum(width(sample02))
    
    
    result <- ifelse(min(widthSample01,widthSample02) > 0, 
                     inter/min(widthSample01,widthSample02),
                     NA)
    return(result)
}


    
    
    
    

