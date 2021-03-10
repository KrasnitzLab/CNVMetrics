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
    cat(attr(x, 'metric'))
    cat("\n")
    
    cat("\nAMPLIFICATION:\n")
    if (nrow(x$AMPLIFICATION) > 5) {
        print(x$AMPLIFICATION[seq_len(5),seq_len(5)], ...)
        if (nrow(x$AMPLIFICATION) == 6) {
            cat("[ -- omitted 1 row/column ]")
        } else {
            cat(paste0("[ -- omitted ", (nrow(x$AMPLIFICATION) - 5), 
                            " rows/columns ]"))
        }
    } else {
        print(x$AMPLIFICATION, ...) 
    }
    cat("\n")
    cat("\nDELETION:\n")
    if (nrow(x$DELETION) > 5) {
        print(x$DELETION[seq_len(5),seq_len(5)], ...)
        if (nrow(x$DELETION) == 6) {
            cat("[ -- omitted 1 row/column ]")
        } else {
            cat(paste0("[ -- omitted ", (nrow(x$DELETION) - 5), 
                        " rows/columns ]"))
        }
    } else {
        print(x$DELETION, ...)
    }
    cat("\n")
    invisible(x)
}



#' @method is CNVMetric
#'
#' @export
is.CNVMetric <- function(x, ...) {
    inherits(x, "CNVMetric")
}



