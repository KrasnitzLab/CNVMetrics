
#' @title Calculate metric using overlapping amplified/deleted regions
#'
#' @description This function calculates a specific metric, as specified by
#' the user, using overlapping
#' regions of specific state between to samples. The metric is calculated for
#' each state separately. When more than 2 samples are
#' present, the metric is calculated for each sample pair. By default, the
#' function is calculating metrics for the AMPLIFICATION and DELETION states.
#' However, the user can specify the list of states to be analyzed.
#'
#' @param segmentData a \code{GRangesList} that contains a collection of
#' genomic ranges representing copy number events, including amplified/deleted
#' status, from at least 2 samples. All samples must have a metadata column
#' called '\code{state}' with a state, in an character string format,
#' specified for each region (ex: DELETION, LOH, AMPLIFICATION, NEUTRAL, etc.).
#'
#' @param states a \code{vector} of \code{character} string with at least one
#' entry. The strings are representing the states that will be analyzed.
#' Default: c('\code{AMPLIFICATION}', '\code{DELETION}').
#'
#' @param method a \code{character} string representing the metric to be used.
#' This should be (an unambiguous abbreviation of) one of "sorensen",
#' "szymkiewicz" or "jaccard". Default: "sorensen".
#'
#' @param nJobs a single positive \code{integer} specifying the number of
#' worker jobs to create in case of distributed computation.
#' Default: \code{1} and always \code{1} for Windows.
#'
#' @details
#'
#' The two methods each estimate the overlap between paired samples. They use
#' different metrics, all in the range [0, 1] with 0 indicating no overlap.
#' The \code{NA} is used when the metric cannot be calculated.
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
#' metric. This object is a list where each entry corresponds to one state
#' specified in the '\code{states}' parameter. Each entry is a \code{matrix}:
#' \itemize{
#' \item{\code{state}}{ a lower-triangular \code{matrix} with the
#'     results of the selected metric on the amplified regions for each paired
#'     samples. The value \code{NA} is present when the metric cannot be
#'     calculated. The value \code{NA} is also present in the top-triangular
#'     section, as well as the diagonal, of the matrix.
#' }
#' }
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
#' demo[["sample01"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841, 32686222),
#'     end=c(2004603, 4577608, 31695808, 32689222)), strand="*",
#'     state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION", "LOH"))
#'
#' demo[["sample02"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1995066, 31611222, 31690000, 32006222),
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
#' ## Calculating Sorensen metric for both AMPLIFICATION and DELETION
#' calculateOverlapMetric(demo, method="sorensen", nJobs=1)
#'
#' ## Calculating Szymkiewicz-Simpson metric on AMPLIFICATION only
#' calculateOverlapMetric(demo, states="AMPLIFICATION", method="szymkiewicz",
#'     nJobs=1)
#'
#' ## Calculating Jaccard metric on LOH only
#' calculateOverlapMetric(demo, states="LOH", method="jaccard", nJobs=1)
#'
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges
#' @importFrom BiocParallel multicoreWorkers SnowParam SerialParam bplapply bptry bpok
#' @encoding UTF-8
#' @export
calculateOverlapMetric <- function(segmentData,
                                    states=c("AMPLIFICATION", "DELETION"),
                                    method=c("sorensen", "szymkiewicz",
                                                "jaccard"),
                                    nJobs=1) {

    ## Select metric method to be used
    method <- match.arg(method)

    ## Validate some parameters
    validateCalculateOverlapMetricParameters(states=states, nJobs=nJobs)

    ## The cnv data must be in a GRangesList format
    if (!is(segmentData, "GRangesList")) {
        stop("the \'segmentData\' argument must be a \'GRangesList\' object")
    }

    names <- names(segmentData)
    nb <- length(names)

    ## At least 2 samples must be present
    if (nb < 2) {
        stop("at least 2 samples must be present in the segmentData")
    }

    ## All samples must have a metadata column called 'state'
    if (!all(vapply(segmentData,
                    FUN=function(x) {"state" %in% colnames(mcols(x))},
                    FUN.VALUE=logical(1)))) {
        stop("at least one sample doesn't have a metadata column ",
                "called \'state\'")
    }

    ## Select the type of parallel environment used for parallel processing
    nbrThreads <- as.integer(nJobs)
    if (nbrThreads == 1 || multicoreWorkers() == 1) {
        coreParam <- SerialParam()
    } else {
        seed <- sample(x=seq_len(999999), size=1)
        coreParam <- SnowParam(workers=nbrThreads, RNGseed=seed)
    }

    ## List that will contain the results
    results <- list()

    ## Loop for each state
    for(type in states) {
        ## Matrix to be filled with the metrics
        dataTMP <- matrix(rep(NA, nb^2), nrow=nb)
        rownames(dataTMP) <- names
        colnames(dataTMP) <- names

        ## Each combinaison that needs to be calculated
        ind <- which(lower.tri(dataTMP, diag=FALSE), arr.ind=TRUE)
        jobSplit <- rep(seq_len(nJobs),
                        each=ceiling(nrow(ind)/nJobs))[seq_len(nrow(ind))]
        entries <- split(as.data.frame(ind), jobSplit)

        ## Running each profile id on a separate thread
        processed <- bptry(bplapply(X=entries, FUN=calculateOneOverlapMetricT,
                                    segmentData=segmentData,
                                    method=method, type=type,
                                    BPPARAM=coreParam))
        ## Check for errors
        if (!all(bpok(processed))) {
            stop("At least one parallel task has failed.")
        }

        ## Assigned the metrics to the final matrix
        for (oneP in processed) {
            for(i in seq_len(nrow(oneP$metric))) {
                dataTMP[oneP$metric$row[i], oneP$metric$col[i]] <-
                                                    oneP$metric$metric[i]
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
#' @description This function calculates a specific metric, as specified
#' by the user, using overlapping amplified/deleted regions between to
#' samples. The metric is calculated for
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
#' @param minThreshold a single positive \code{numeric} setting the minimum
#' value to consider two segments as different during the metric calculation.
#' If the absolute difference is below or equal to threshold, the difference
#' will be replaced by zero. Default: 0.2.
#'
#' @param excludedRegions an optional \code{GRanges} containing the regions
#' that have to be excluded for the metric calculation. Default: \code{NULL}.
#'
#' @param nJobs a single positive \code{integer} specifying the number of
#' worker jobs to create in case of distributed computation.
#' Default: \code{1} and always \code{1} for Windows.
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
#' demo[["sample01"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841),
#'     end=c(2004603, 4577608, 31695808)), strand="*",
#'     log2ratio=c(2.5555, 1.9932, -0.9999))
#'
#' demo[["sample02"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1995066, 31611222, 31690000),
#'     end=c(2204505, 31689898, 31895666)), strand=c("-", "+", "+"),
#'     log2ratio=c(0.3422, 0.5454, -1.4444))
#'
#' ## The amplified region in sample03 is a subset of the amplified regions
#' ## in sample01
#' demo[["sample03"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1906069, 4558838),
#'     end=c(1909505, 4570601)), strand="*",
#'     log2ratio=c(3.2222, -1.3232))
#'
#' ## Calculating Sorensen metric
#' calculateLog2ratioMetric(demo, method="weightedEuclideanDistance", nJobs=1)
#'
#'
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges
#' @importFrom BiocParallel multicoreWorkers SnowParam SerialParam bplapply bptry bpok
#' @encoding UTF-8
#' @export
calculateLog2ratioMetric <- function(segmentData,
                                    method=c("weightedEuclideanDistance"),
                                    minThreshold=0.2, excludedRegions=NULL,
                                    nJobs=1) {

    ## Select metric method to be used
    method <- match.arg(method)

    ## Validate some parameters
    validatecalculateLog2ratioMetricParameters(minThreshold=minThreshold,
                                excludedRegions=excludedRegions, nJobs=nJobs)

    ## The CNV data must be in a GRangesList format
    if (!is(segmentData, "GRangesList")) {
        stop("the \'segmentData\' argument must be a \'GRangesList\' object")
    }

    names <- names(segmentData)
    nb <- length(names)

    ## At least 2 samples must be present
    if (nb < 2) {
        stop("at least 2 samples must be present in segmentData")
    }

    ## All samples must have a metadata column called 'log2ratio' with
    ## log2ratio values
    if (!all(vapply(segmentData,
                    FUN=function(x) {"log2ratio" %in% colnames(mcols(x))},
                    FUN.VALUE=logical(1)))) {
        stop("at least one sample doesn't have a metadata column ",
                "called \'log2ratio\'")
    }

    ## Select the type of parallel environment used for parallel processing
    nbrThreads <- as.integer(nJobs)
    if (nbrThreads == 1 || multicoreWorkers() == 1) {
        coreParam <- SerialParam()
    } else {
        seed <- sample(x=seq_len(999999), size=1)
        coreParam <- SnowParam(workers=nbrThreads, RNGseed=seed)
    }

    ## List that will contain the results
    results <- list()

    ## Loop for each state
    for(type in c("LOG2RATIO")) {
        ## Matrix to be filled with the metrics
        dataTMP <- matrix(rep(NA, nb^2), nrow=nb)
        rownames(dataTMP) <- names
        colnames(dataTMP) <- names

        ## Each combination that needs to be calculated
        ind <- which(lower.tri(dataTMP, diag=FALSE), arr.ind=TRUE)
        jobSplit <- rep(seq_len(nJobs),
                        each=ceiling(nrow(ind)/nJobs))[seq_len(nrow(ind))]
        entries <- split(as.data.frame(ind), jobSplit)

        ## Running each profile id on a separate thread
        processed <- bptry(bplapply(X=entries,
                                        FUN=calculateOneLog2valueMetricT,
                                        segmentData=segmentData,
                                        method=method,
                                        minThreshold=minThreshold,
                                        bedExclusion=excludedRegions,
                                        BPPARAM=coreParam))
        ## Check for errors
        if (!all(bpok(processed))) {
            stop("At least one parallel task has failed.")
        }

        ## Assigned the metrics to the final matrix
        for (oneP in processed) {
            for(i in seq_len(nrow(oneP$metric))) {
                dataTMP[oneP$metric$row[i], oneP$metric$col[i]] <-
                    oneP$metric$metric[i]
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


#' @title Plot metrics present in a \code{CNVMetric} object
#'
#' @description This function plots one heatmap (or two heatmaps) of the
#' metrics present in
#' a \code{CNVMetric} object. For the overlapping metrics, the user can select
#' to print the heatmap related to amplified or deleted regions or both. The
#' \code{NA} values present in the metric matrix are transformed into zero for
#' the creation of the heatmap.
#'
#' @param metric a \code{CNVMetric} object containing the metrics calculated
#' by \code{calculateOverlapMetric} or by \code{calculateLog2ratioMetric}.
#'
#' @param type a single \code{character} string indicating which graph
#' to generate. This should be a type present in the \code{CNVMetric} object or
#' "\code{ALL}". This is useful for the
#' overlapping metrics that have multiple types specified by the user.
#' Default: "\code{ALL}".
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
#'  "\code{ALL}".
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
#' demo[["sample01"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841),
#'     end=c(2004603, 4577608, 31695808)), strand="*",
#'     state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#'
#' demo[["sample02"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1995066, 31611222, 31690000),
#'     end=c(2204505, 31689898, 31895666)), strand=c("-", "+", "+"),
#'     state=c("AMPLIFICATION", "AMPLIFICATION", "DELETION"))
#'
#' ## The amplified region in sample03 is a subset of the amplified regions
#' ## in sample01
#' demo[["sample03"]] <- GRanges(seqnames="chr1",
#'     ranges=IRanges(start=c(1906069, 4558838),
#'     end=c(1909505, 4570601)), strand="*",
#'     state=c("AMPLIFICATION", "DELETION"))
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
plotMetric <- function(metric, type="ALL",
                        colorRange=c(c("white", "darkblue")),
                        show_colnames=FALSE, silent=TRUE, ...) {

    ## Validate that the metric parameter is a CNVMetric object
    if (!is.CNVMetric(metric)) {
        stop("\'metric\' must be a CNVMetric object.")
    }

    ## Validate that the filename argument is not used when
    ## type "ALL" is selected
    if (type == "ALL" &&  length(names(metric)) > 1 & hasArg("filename")) {
        stop("\'type\' cannot be \'ALL\' when filename argument is used.")
    }

    ## Validate that the type argument is a single character string
    if (!is.character(type) | length(type) > 1) {
        stop("the \'type\' must be a single character string")
    }

    ## Validate that the type argument is present in the object
    if (type != "ALL" & ! type %in% names(metric)) {
        stop("the specified \'type\' is not present in this metric object.")
    }

    ## Validate that the color name has only one value
    if (!is.character(colorRange) || length(colorRange) < 2) {
        stop("\'colorRange\' must be a vector of 2 color names.")
    }

    ## Validate that the color name is valid
    tryCatch(col2rgb(colorRange), error=function(e) {
        stop("\'colorRange\' must be be a vector of 2 valid color names.")
    })

    ## Extract the type of metric
    metricInfo <- attributes(metric)$metric

    plot_list <- list()

    nameList <- ifelse(type != "ALL", type, names(metric))

    for (name in nameList) {
        plot_list[[name]] <- plotOneMetric(metric=metric,
                                type=name, colorRange=colorRange,
                                show_colnames=show_colnames,
                                silent=silent, ...)
    }

    n_col <- length(plot_list)

    grid.arrange(arrangeGrob(grobs=plot_list, ncol=n_col))
}



