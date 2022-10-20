#' CNVMetrics: Copy number variant metrics
#'
#' The CNVMetrics package calculates similarity metrics to
#' facilitate copy number variant comparison among samples and/or methods.
#' Similarity metrics can be employed to compare CNV profiles of genetically
#' unrelated samples as well as those with a common genetic background.
#' Some metrics are based on the shared amplified/deleted regions while other
#' metrics rely on the level of amplification/deletion.
#' The data type used as input is a plain text file containing the genomic
#' position of the copy number variations, as well as the status and/or
#' the log2 ratio values.
#' Finally, a visualization tool is provided to explore resulting metrics.
#'
#' @docType package
#'
#' @name CNVMetrics-package
#'
#' @aliases CNVMetrics-package CNVMetrics
#'
#' @author Astrid Deschênes, Pascal Belleau, David A. Tuveson and
#' Alexander Krasnitz
#'
#' Maintainer:
#' Astrid Deschênes <adeschen@hotmail.com>
#'
#' @seealso
#' \itemize{
#'     \item \code{\link{calculateOverlapMetric}} {for calculating metric
#'     using overlapping amplified/deleted regions}
#'     \item \code{\link{calculateLog2ratioMetric}} {for calculating metric
#'     using log2ratio values}
#'     \item \code{\link{processSim}} {for generating simulations}
#'     \item \code{\link{plotMetric}} {for plotting metrics}
#' }
#'
#' @encoding UTF-8
#' @keywords package
NULL
