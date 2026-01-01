
#' @title Calculate metric using overlapping amplified/deleted regions between
#' two samples.
#'
#' @description Calculate a specific metric using overlapping
#' amplified/deleted regions between two samples.
#'
#' @param entry a \code{list} which contains the row and column indexes
#' (always in this order) of
#' the metric in the final matrix. Those values correspond to the positions
#' of the two samples used
#' to calculate the metric in the \code{GRangesList} (\code{segmentData}).
#'
#' @param segmentData a \code{GRangesList} that contains a collection of
#' genomic ranges representing copy number events, including amplified/deleted
#' status, from at least 2 samples. All samples must have a metadata column
#' called '\code{state}' with a state, in an character string format,
#' specified for each region (ex: DELETION, LOH, AMPLIFICATION, NEUTRAL, etc.).
#'
#' @param method a \code{character} string representing the metric to be
#' used ('\code{sorensen}' or '\code{szymkiewicz}'.
#'
#' @param type a \code{character} string representing the type of
#' copy number events to be used ('\code{AMPLIFICATION}' or '\code{DELETION}').
#'
#' @return a \code{list} containing 1 entry:
#' \itemize{
#' \item{\code{metric} a \code{data.frame}, which contains 3 columns. The 2
#' first columns, called \code{row} and \code{column} correspond to the
#' indexes of the metric in the final matrix. Those
#' 2 first columns match to the \code{entry} parameter. The third column,
#' called \code{metric},
#' contains the values of the specified metric for each combination.
#' If the metric cannot be calculated, \code{NA} is present.
#' }
#' }
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
#' demo[["sample01"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841, 32686222),
#'     end=c(2004603, 4577608, 31695808, 32689222)), strand="*",
#'     state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION", "LOH"))
#'
#' demo[["sample02"]] <- GRanges(seqnames="chr1",
#'     ranges= IRanges(start=c(1995066, 31611222, 31690000, 32006222),
#'     end=c(2204505, 31689898, 31895666, 32789233)),
#'     strand=c("-", "+", "+", "+"),
#'     state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION", "LOH"))
#'
#' ## The amplified region in sample03 is a subset of the amplified regions
#' ## in sample01
#' demo[["sample03"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1906069, 4558838),
#'     end=c(1909505, 4570601)), strand="*",
#'     state=c("AMPLIFICATION", "DELETION"))
#'
#' ## The 2 samples used to calculate the metric
#' entries <- data.frame(row=c(2, 3), col=c(1, 1))
#'
#' ## Calculate Sorensen metric for the amplified regions on samples 2 and 3
#' CNVMetrics:::calculateOneOverlapMetricT(entry=entries, segmentData=demo,
#'     method="sorensen", type="AMPLIFICATION")
#'
#' ## Calculate Szymkiewicz-Simpson metric for the amplified regions
#' ## in samples 1 and 2
#' ## Amplified regions of sample02 are a subset of the amplified
#' ## regions in sample01
#' CNVMetrics:::calculateOneOverlapMetricT(entry=entries, segmentData=demo,
#'     method="szymkiewicz", type="AMPLIFICATION")
#'
#' ## Calculate Sorensen metric for the deleted regions in samples 1 and 2
#' CNVMetrics:::calculateOneOverlapMetricT(entry=entries, segmentData=demo,
#'     method="sorensen", type="DELETION")
#'
#' @author Astrid Deschênes
#' @encoding UTF-8
#' @keywords internal
calculateOneOverlapMetricT <- function(entry, segmentData, method, type) {

    entries <- split(entry, entry$col)

    results <- list()

    for (data in entries) {

        result <- data
        result$metric <- rep(NA, nrow(result))

        sample02 <- segmentData[[data$col[1]]]
        sample02 <- sample02[sample02$state == type,]

        for(i in seq_len(nrow(data))) {
            sample01 <- segmentData[[data$row[i]]]
            sample01 <- sample01[sample01$state == type,]

            if (length(sample01) > 0 && length(sample02) > 0) {
                result$metric[i] <- switch(method,
                sorensen=calculateSorensen(sample01, sample02),
                    szymkiewicz=calculateSzymkiewicz(sample01, sample02),
                    jaccard=calculateJaccard(sample01, sample02))
            }
        }

        results[[length(results) + 1]] <- result
    }

    results <- do.call(rbind, results)

    return(result <- list(metric=results))
}

#' @title Calculate Sorensen metric
#'
#' @description Calculate Sorensen metric using overlapping regions between
#' two samples.
#'
#' @param sample01 a \code{GRanges} which contains a collection of
#' genomic ranges representing copy number events for the first sample.
#' @param sample02 a \code{GRanges} which contains a collection of
#' genomic ranges representing copy number events for the second sample.
#'
#' @details
#'
#' The method calculates the Sorensen metric using overlapping
#' regions between the samples. All regions present in both samples are used
#' for the calculation of the metric.
#'
#' The Sorensen metric is calculated by dividing twice the size of
#' the intersection by the sum of the size of the two sets. If the sum of
#' the size of the two sets is zero; the value \code{NA} is
#' returned instead. The strand of the regions is not taken into account while
#' calculating the intersection.
#'
#'
#' @return a \code{numeric}, the value of the Sorensen metric. If
#' the metric cannot be calculated, \code{NA} is returned.
#'
#' @references
#'
#' Sørensen, Thorvald. n.d. “A Method of Establishing Groups of Equal
#' Amplitude in Plant Sociology Based on Similarity of Species and Its
#' Application to Analyses of the Vegetation on Danish Commons.”
#' Biologiske Skrifter, no. 5: 1–34.
#'
#' @examples
#'
#' ## Load required package to generate the two samples
#' require(GenomicRanges)
#'
#' ## Generate two samples with identical sequence levels
#' sample01 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841),
#'     end=c(2004603, 4577608, 31695808)), strand="*")
#' sample02 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1995066, 31611222),
#'     end=c(2204505, 31689898)), strand="*")
#'
#' ## Calculate Sorensen metric
#' CNVMetrics:::calculateSorensen(sample01=sample01, sample02=sample02)
#'
#' @author Astrid Deschênes
#' @importFrom GenomicRanges intersect width
#' @encoding UTF-8
#' @keywords internal
calculateSorensen <- function(sample01, sample02) {

    ## Calculate intersection between the two sets as well as the
    ## total size of each set
    inter <- sum(as.numeric(width(intersect(sample01, sample02,
                                                ignore.strand=TRUE))))
    widthSample01 <- sum(as.numeric(width(sample01)))
    widthSample02 <- sum(as.numeric(width(sample02)))

    ## Calculate Sorensen metric if possible; otherwise NA
    result <- ifelse((widthSample01 + widthSample02) > 0,
                        (2.0 * inter)/(widthSample01 + widthSample02),
                        NA)
    return(result)
}

#' @title Calculate Szymkiewicz-Simpson metric
#'
#' @description Calculate Szymkiewicz-Simpson metric using overlapping
#' regions between two samples.
#'
#' @param sample01 a \code{GRanges} which contains a collection of
#' genomic ranges representing copy number events for the first sample.
#' @param sample02 a \code{GRanges} which contains a collection of
#' genomic ranges representing copy number events for the second sample.
#'
#' @details
#'
#' The method calculates the Szymkiewicz-Simpson metric using overlapping
#' regions between the samples. All regions present in both samples all used
#' for the calculation of the metric.
#'
#' The Szymkiewicz-Simpson metric is calculated by dividing the size of
#' the intersection by the smaller of the size of the two sets. If one sample
#' has a size of zero, the metric is not calculated; the value \code{NA} is
#' returned instead. The strand of the regions is not taken into account while
#' calculating the intersection.
#'
#' @return a \code{numeric}, the value of the Szymkiewicz-Simpson metric. If
#' the metric cannot be calculated, \code{NA} is returned.
#'
#' @references
#'
#' Vijaymeena, M. K, and Kavitha K. 2016. “A Survey on Similarity Measures in
#' Text Mining.” Machine Learning and Applications: An International
#' Journal 3 (1): 19–28. doi: \url{https://doi.org/10.5121/mlaij.2016.3103}
#'
#' @examples
#'
#' ## Load required package to generate the two samples
#' require(GenomicRanges)
#'
#' ## Generate two samples with identical sequence levels
#' sample01 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841),
#'     end=c(2004603, 4577608, 31695808)), strand="*")
#' sample02 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1995066, 31611222),
#'     end=c(2204505, 31689898)), strand=c("+", "-"))
#'
#' ## Calculate Szymkiewicz-Simpson metric
#' CNVMetrics:::calculateSzymkiewicz(sample01=sample01, sample02=sample02)
#'
#' @author Astrid Deschênes
#' @importFrom GenomicRanges intersect width
#' @encoding UTF-8
#' @keywords internal
calculateSzymkiewicz <- function(sample01, sample02) {

    ## Calculate intersection between the two sets as well as the
    ## total size of each set
    inter <- sum(as.numeric(width(intersect(sample01, sample02,
                                                ignore.strand=TRUE))))
    widthSample01 <- sum(as.numeric(width(sample01)))
    widthSample02 <- sum(as.numeric(width(sample02)))

    ## Calculate Szymkiewicz-Simpson metric if possible; otherwise NA
    result <- ifelse(min(widthSample01,widthSample02) > 0,
                        inter/min(widthSample01,widthSample02),
                        NA)
    return(result)
}



#' @title Calculate Jaccard metric
#'
#' @description Calculate Jaccard metric using overlapping regions between
#' two samples.
#'
#' @param sample01 a \code{GRanges} which contains a collection of
#' genomic ranges representing copy number events for the first sample.
#' @param sample02 a \code{GRanges} which contains a collection of
#' genomic ranges representing copy number events for the second sample.
#'
#' @details
#'
#' The method calculates the Jaccard metric using overlapping
#' regions between the samples. All regions present in both samples are used
#' for the calculation of the metric.
#'
#' The Jaccard metric is calculated by dividing the size of
#' the intersection by the size of the union of the two sets. If the
#' the size of the union of the two sets is zero; the value \code{NA} is
#' returned instead. The strand of the regions is not taken into account while
#' calculating the intersection.
#'
#'
#' @return a \code{numeric}, the value of the Jaccard metric. If
#' the metric cannot be calculated, \code{NA} is returned.
#'
#' @references
#'
#' Jaccard, P. (1912), The Distribution of the Flora in the Alpine Zone.
#' New Phytologist, 11: 37-50.
#' DOI: \url{https://doi.org/10.1111/j.1469-8137.1912.tb05611.x}
#'
#' @examples
#'
#' ## Load required package to generate the two samples
#' require(GenomicRanges)
#'
#' ## Generate two samples with identical sequence levels
#' sample01 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841),
#'     end=c(2004603, 4577608, 31695808)), strand="*")
#' sample02 <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1995066, 31611222),
#'     end=c(2204505, 31689898)), strand="*")
#'
#' ## Calculate Sorensen metric
#' CNVMetrics:::calculateJaccard(sample01=sample01, sample02=sample02)
#'
#' @author Astrid Deschênes
#' @importFrom GenomicRanges intersect width
#' @encoding UTF-8
#' @keywords internal
calculateJaccard <- function(sample01, sample02) {

    ## Calculate intersection between the two sets as well as the
    ## total size of each set
    inter <- sum(as.numeric(width(intersect(sample01, sample02,
                                                ignore.strand=TRUE))))
    widthSample01 <- sum(as.numeric(width(sample01)))
    widthSample02 <- sum(as.numeric(width(sample02)))

    ## Calculate Jaccard metric if possible; otherwise NA
    result <- ifelse((widthSample01 + widthSample02 - inter) > 0,
                        (inter)/(widthSample01 + widthSample02 - inter),
                        NA)
    return(result)
}


