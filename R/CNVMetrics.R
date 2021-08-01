#' CNVMetrics: Copy number variant metrics
#'
#' The CNVMetrics package calculates similarity metrics to facilitate copy 
#' number variant comparison among samples and/or methods. Some metrics are 
#' based on the shared amplified/deleted regions while other metrics rely on 
#' the level of amplification/deletion (log2 ratio). Finally, a 
#' visualization tool is provided to explore resulting metrics.
#'
#' @docType package
#'
#' @name CNVMetrics-package
#'
#' @aliases CNVMetrics-package CNVMetrics
#'
#' @author Astrid Deschênes and
#' Pascal Belleau
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
#'     \item \code{\link{plotMetric}} {for plotting metrics}
#' }
#' 
#' @encoding UTF-8
#' @keywords package
NULL
