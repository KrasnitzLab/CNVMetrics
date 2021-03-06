% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CNVMetricsInternalMethods.R
\name{doRegression}
\alias{doRegression}
\title{Linear regression analysis between each paired samples using
segment values.}
\usage{
doRegression(segmentData)
}
\arguments{
\item{segmentData}{a \code{list} of that 
contains the segments from multiple files. The \code{list} is composed of 
those entries:
\itemize{
\item a \code{segment} that contains the \code{GRanges} with the segment
information.
}}
}
\value{
a \code{list} of that 
contains the segments from multiple files. The \code{list} is composed of 
those entries:
\itemize{
\item a \code{segment} that contains the \code{GRanges} with the segment
information.
\item a \code{regression} that contains the result of the paired 
regressions.
}
}
\description{
Do an linear regression between each paired samples using
segment values only from the included segments. The segment file used as
reference is the same for all analysis.
}
\examples{

# TODO

}
\author{
Astrid Deschenes, Pascal Belleau
}
\keyword{internal}
