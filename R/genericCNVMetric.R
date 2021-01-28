#' @title Print CNVMetric object
#' 
#' @description Print a \code{CNVMetric} object and returns it 
#' \code{invisibly}.
#' 
#' @method print CNVMetric
#' 
#' @param x the output object from \code{calculateOverlapRegionsMetric} 
#'     function to be printed.
#' @param \ldots further arguments passed to or from other methods.
#' 
#' @return the argument \code{x}.
#' 
#' @seealso 
#' 
#' The default method \code{\link[base]{print.default}}.
#' 
#' @export
print.CNVMetric <- function(x, ...) {
    # Print title before printing the content of the regression object
    cat("CNV Metric done on overlapping regions\n")
    cat("Metric:\n")
    print(attr(x, 'metric'))
    cat("\nAMPLIFICATION:\n")
    print(x$AMPLIFICATION, ...)
    cat("\nDELETION:\n")
    print(x$DELETION, ...)
    invisible(x)
}



#' @method is CNVMetric
#'
#' @export
is.CNVMetric <- function(x, ...) {
    inherits(x, "CNVMetric")
}

