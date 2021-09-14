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
    cat("CNV Metrics\n")
    cat("Metric:\n")
    cat(attr(x, 'metric'))
    cat("\n")
    
    for(i in seq_len(length(x))) {
        cat(paste0("\n", names(x)[1], ":\n"))
        dataM <- x[[i]]
        if (nrow(dataM) > 5) {
            print(dataM[seq_len(5), seq_len(5)], ...)
            if (nrow(dataM) == 6) {
                cat("[ -- omitted 1 row/column ]")
            } else {
                cat(paste0("[ -- omitted ", (nrow(dataM) - 5), 
                            " rows/columns ]"))
            }
        } else {
            print(dataM, ...) 
        }
        cat("\n")
    }
    invisible(x)
}

#' @title Is an object of class \code{CNVMetric}
#' 
#' @description Functions to test inheritance relationships between an 
#' object and class \code{CNVMetric}.
#' 
#' @method is CNVMetric
#'
#' @param x an object.
#' 
#' @param \ldots further arguments passed to or from other methods.
#' 
#' @return a \code{logical}.
#' 
#' @importFrom methods is
#' @export
is.CNVMetric <- function(x, ...) {
    inherits(x, "CNVMetric")
}



