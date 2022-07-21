
#' @title Generate a simulated chromosome based on a reference sample
#'
#' @description The function generates a list of simulated segments
#' that represent a simulated chromosome based on a reference sample
#' specified by the user. The function only accounts for the positions where
#' a segment is assigned. In addition, the total number of segments
#' is preserved. A Dirichlet distribution is used to assigned new sizes
#' to the segments with respect to the
#' relative initial size of the segment. Then, those new segments are shuffled
#' without replacement. The positions are replaced by values between zero and
#' one that represent the relative position in a chromosome where positions
#' without segment have been removed.
#' To ensure valuable results, the reference sample
#' should have segments covering
#' a good proportion of the chromosome; those should include NEUTRAL segments.
#'
#' @param curSample a \code{GRanges} that contains a collection of
#' genomic ranges representing copy number events, including amplified/deleted
#' status, from exactly one sample. The sample must have a metadata column
#' called '\code{state}' with a state, in an character string format,
#' specified for each region (ex: DELETION, LOH, AMPLIFICATION, NEUTRAL, etc.)
#' and a metadata column called '\code{CN}' that contains the log2 copy
#' number ratios.
#'
#' @param chrCur a \code{character} string representing the name of the
#' chromosome that is used as reference for the simulation.
#'
#' @param nbSim a single positive \code{integer} which is corresponding to
#' the number of simulations that will be generated.
#'
#' @details TODO
#'
#'
#' @return a code{list} containing one entry per simulation. Each entry is
#' a \code{data.frame} containing shuffled segments with 6 columns:
#' \itemize{
#' \item{\code{ID}}{ The name of the simulation. }
#' \item{\code{chr}}{ The name fo the chromosome. }
#' \item{\code{start}}{ The starting position of the segment; the positions
#' are between zero and one. The segment width is representing the
#' proportional size of the segment relative to the global segment size.}
#' \item{\code{end}}{ The ending position of the segment; the positions
#' are between zero and one. The segment width is representing the
#' proportional size of the segment relative to the global segment size. }
#' \item{\code{log2ratio}} { The log2 copy number ratio assigned to
#' the segment. }
#' \item{\code{state}} { The state of the region (ex: DELETION, LOH,
#' AMPLIFICATION, NEUTRAL, etc.). }
#' }
#'
#' @examples
#'
#' ## Load required package to generate the samples
#' require(GenomicRanges)
#'
#' ## Create one 'demo' genome with 2 chromosomes
#' ## in a GRanges object
#' ## The stand of the regions doesn't affect the calculation of the metric
#' sample01 <- GRanges(seqnames=c(rep("chr1", 4), rep("chr2", 3)),
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841, 32686222,
#'         1, 120331, 725531),
#'     end=c(2004603, 4577608, 31695808, 32689222, 117121,
#'         325555, 1225582)),
#'     strand="*",
#'     state=c("AMPLIFICATION", "NEUTRAL", "DELETION", "LOH",
#'         "DELETION", "NEUTRAL", "NEUTRAL"),
#'     log2ratio=(c(0.5849625, 0, -1, -1, -0.87777, 0, 0)))
#'
#'
#' ## Generates 10 simulated chromosomes (one chromosome per simulated sample)
#' ## based on chromosome 2 from the input sample.
#' ## The shuffled chromosomes have a start and an end between 0 an 1
#' CNVMetrics:::simChr(curSample=sample01, chrCur="chr2", nbSim=10)
#'
#' ## Generates 4 simulated chromosomes (one chromosome per simulated sample)
#' ## based on chromosome 1 from the input sample.
#' ## The shuffled chromosomes have a start and an end between 0 an 1
#' CNVMetrics:::simChr(curSample=sample01, chrCur="chr1", nbSim=4)
#'
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges
#' @importFrom rBeta2009 rdirichlet
#' @encoding UTF-8
#' @keywords internal
simChr <- function(curSample, chrCur, nbSim) {

    ## Extract genomic regions on the selected chromosome
    curSampleGR <- curSample[seqnames(curSample) == chrCur]
    ## Extract minininum and maximum positions
    minStart <- min(start(curSampleGR))
    maxEnd <- max(end(curSampleGR))

    listW <- width(curSampleGR)
    sizeT <- sum(listW)
    listWP <- listW / (sizeT * 0.001)
    nbSeg <- length(listW)
    listEvents <- curSampleGR$state
    listCN <- curSampleGR$log2ratio

    ## Create partitions representing the new segments using dirichlet with
    ## the proportion of the segment relative to the total size of the
    ## chromosome (covered with segment) as the shape parameters
    if(length(listWP) == 2) {
        tmp <- round(suppressWarnings(rdirichlet(nbSim, listWP)) * sizeT)
        partDir <- unname(cbind(tmp, sizeT - tmp))
    }else if(length(listWP) == 1) {
        ## The simulated segment has the same size that
        ## the only existing segment
        partDir <- matrix(rep(sizeT, nbSim), ncol = 1)
    }else{
        ## The simulated segments have a shape proportional to the
        ## relative size of the existing segments relative to the total size
        ## of the chromosome
        partDir <- round(rdirichlet(n=nbSim, shape=listWP) * sizeT)
    }

    # Adjust the length
    if(dim(partDir)[2] > 1) {
        partDir <- t(apply(partDir, 1, FUN = function(x,sizeT) {
            v <- sizeT - sum(x)
            x[which.max(x)] <- x[which.max(x)] + v
            return(x)
        }, sizeT=sizeT))

        ## Shuffle the segments (without replacement)
        orderPart <- t(replicate(nbSim, sample(seq_len(nbSeg), size=nbSeg,
                                                replace=FALSE)))
        partDir <- t(apply(cbind(partDir, orderPart), 1, FUN=function(x) {
                                x[x[(length(x)/2+1):(length(x))]]}))
        partDir <-  partDir/sizeT

        partCN <- t(apply(orderPart, 1, FUN=function(x, listCN) {
                                    listCN[x]}, listCN=listCN))
        partEvents <- t(apply(orderPart, 1, FUN=function(x, listEvents) {
                                    listEvents[x]}, listEvents=listEvents))
        partDir <- t(apply(partDir, 1, cumsum))
    }else{
        partDir <-  partDir / sizeT

        partCN <- matrix(rep(listCN, nbSim), ncol=1)
        partEvents <- matrix(rep(listEvents, nbSim), ncol=1)
    }


    ## Final returned list with all simulated samples
    res <- list()

    ## Create one entry per simulated samples
    for(i in seq_len(nbSim)) {
        ## A data.frame containing the information about the simulated
        ## segments
        ## The start and end positions are between zero and one.
        res[[i]] <- data.frame(ID=rep(paste0("S", i), ncol(partDir)),
                               chr=rep(chrCur, ncol(partDir)),
                               start=c(0, partDir[i, -1 * ncol(partDir)]),
                               end=partDir[i, ],
                               log2ratio=partCN[i, ],
                               state=partEvents[i, ],
                               stringsAsFactors=FALSE)
    }

    return(res)
}


#' @title TODO
#'
#' @description TODO
#'
#' @param curSample a \code{GRangesList} that contains a collection of
#' genomic ranges representing copy number events, including amplified/deleted
#' status, from at least 1 samples. The sample must have a metadata column
#' called '\code{state}' with a state, in an character string format,
#' specified for each region (ex: DELETION, LOH, AMPLIFICATION, NEUTRAL, etc.).
#'
#' @param dfChr data.frame TODO
#'
#' @param chrCur a \code{integer} TODO
#'
#' @details TODO
#'
#'
#' @return df TODO
#'
#'
#'
#' @examples
#' ## TODO
#'
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges
#' @encoding UTF-8
#' @keywords internal
processChr <- function(curSample, dfChr, chrCur) {
    curSampleGR <- curSample[seqnames(curSample) == chrCur]
    listStart <- start(curSampleGR)
    listOrd <- order(listStart)
    listStart <- listStart[listOrd]

    listEnd <- end(curSampleGR)
    listEnd <- listEnd[listOrd]

    minStart <- min(listStart)
    maxEnd <- max(listEnd)

    listW <- width(curSampleGR)[listOrd]
    sizeT <- sum(listW)

    listHole <- data.frame(start = listEnd[-1 * length(listEnd)],
                           end = listStart[-1])
    listHoles <- listHole[listHole$end - listHole$start > 100,]
    newChr <- list()


    curW <- as.integer(round((dfChr$end - dfChr$start) * sizeT))

    df <- data.frame(ID = dfChr[, "ID"],
                     chr = rep(chrCur ,nrow(dfChr)),
                     start = minStart + c(0, cumsum(curW[-1*length(curW)]) + 1),
                     end = minStart + cumsum(curW),
                     log2ratio = dfChr[, "log2ratio"],
                     state = dfChr[, "state"],
                     stringsAsFactors = FALSE)

    if(nrow(listHole) > 0){
        for(j in seq_len(nrow(listHole))){
            # z <- cbind(c(df$start+1, df$end+1, listHole$start[j]),
            #            c(seq_len(nrow(df)),
            #              -1 * seq_len(nrow(df)),
            #              0),
            #            c(rep(0, nrow(df)),
            #              rep(0, nrow(df)),
            #              1))
            # z <- z[order(z[,1]),]
            # pos <- cumsum(z[,2])[z[,3]==1]
            pos <- which(df$start < listHole$start[j] & df$end >= listHole$start[j])
            if(listHole$start[j] < df$end[pos]){
                tmp <- df[pos, ,drop = FALSE]
                df[pos,"end"] <- listHole$start[j]

                tmp$start <- listHole$end[j]
                tmp$end <- tmp$end + listHole$end[j] - listHole$start[j]
                if(pos < nrow(df)){
                    df$start[(pos+1):nrow(df)] <- df$start[(pos+1):nrow(df)] + listHole$end[j] - listHole$start[j]
                    df$end[(pos+1):nrow(df)] <- df$end[(pos+1):nrow(df)] + listHole$end[j] - listHole$start[j]
                }
                df <- rbind(df, tmp)
                df <- df[order(df$start),]
            }
        }
    }
    #newChr[[i]] <- df

    return(df)
}


#' @title Generate simulated samples with copy number profiles derived from
#' a specific sample
#'
#' @description The function uses the input sample to
#' simulate new samples. The simulated samples will possess similar
#' sizes of events, proportional to the original chromosome.
#' To generate realistic simulations, the specified sample
#' must contain segments covering the majority of the genome. Most
#' importantly, the NEUTRAL segments should be present.
#'
#' @param curSample a \code{GRanges} that contains a collection of
#' genomic ranges representing copy number events, including amplified/deleted
#' status, from one sample. The sample must have a metadata column
#' called '\code{state}' with a state, in an character string format,
#' specified for each region (ex: DELETION, LOH, AMPLIFICATION, NEUTRAL, etc.)
#' and a metadata column called '\code{CN}' that contains the log2 copy
#' number ratios.
#'
#' @param nbSim a single positive \code{integer} which is corresponding to the
#'  number of simulations that will be generated.
#'
#' @details TODO
#'
#'
#' @return  a \code{data.frame} containing the segments for each
#' simulated sample. The \code{data.frame} has 6 columns:
#' \itemize{
#' \item{\code{ID}}{ a \code{character} string, the name of the simulated
#' sample }
#' \item{\code{chr}}{ a \code{character} string, the name fo the chromosome }
#' \item{\code{start}}{ a \code{integer}, the starting position of the
#' segment }
#' \item{\code{end}}{ a \code{integer}, the ending position of the segment }
#' \item{\code{log2ratio}} { a \code{numerical}, the log2 copy number
#' ratio assigned to the segment }
#' \item{\code{state}} { a \code{character} string, the state of the segment
#' (ex: DELETION, AMPLIFICATION, NEUTRAL, etc.) }
#' }
#'
#'
#'
#' @examples
#'
#' ## Load required package to generate the samples
#' require(GenomicRanges)
#'
#' ## Create one 'demo' genome with 2 chromosomes and few segments
#' ## The stand of the regions doesn't affect the calculation of the metric
#' sample01 <- GRanges(seqnames=c(rep("chr1", 4), rep("chr2", 3)),
#'     ranges=IRanges(start=c(1905048, 4554832, 31686841, 32686222,
#'         1, 120331, 725531),
#'     end=c(2004603, 4577608, 31695808, 32689222, 117121,
#'         325555, 1225582)),
#'     strand="*",
#'     state=c("AMPLIFICATION", "NEUTRAL", "DELETION", "LOH",
#'         "DELETION", "NEUTRAL", "NEUTRAL"),
#'     log2ratio=c(0.5849625, 0, -1, -1, -0.87777, 0, 0))
#'
#' ## Generates 10 simulated genomes based on the 'demo' genome
#' simRes <- processSim(curSample=sample01, nbSim=10)
#'
#' @author Astrid Deschênes, Pascal Belleau
#' @import GenomicRanges
#' @importFrom S4Vectors isSingleInteger
#' @encoding UTF-8
#' @export
processSim <- function(curSample, nbSim) {

    ## Validate that nbSim is an positive integer
    if (!(isSingleInteger(nbSim) || isSingleNumber(nbSim)) ||
        as.integer(nbSim) < 1) {
        stop("nbSim must be a positive integer")
    }

    ## Validate that curSample is a GRanges object
    if (!is(curSample, "GRanges")) {
        stop("the \'curSample\' argument must be a \'GRanges\' object")
    }

    ## The sample must have a metadata column called 'log2ratio' with
    ## log2ratio values
    if (! "log2ratio" %in% colnames(mcols(curSample))) {
        stop("the sample must have a metadata column ",
             "called \'log2ratio\'")
    }

    ## The list of unique chromosomes in the sample
    listChr <- as.character(unique(seqnames(curSample)))

    ## Generated list of shuffled chromosomes
    resChr <- list()

    ## Create shuffled chromosomes, one chromosome at the time
    ## There is as many shuffled chromosomes created as the number
    ## of simulations specified
    for(chr in listChr) {
        resChr[[chr]] <- simChr(curSample=curSample, chrCur=chr, nbSim=nbSim)
    }

    ## Generated list of simulated samples
    df <- list()

    ##
    for(chr in listChr) {
        newChr <- list()
        for(i in seq_len(nbSim)){
            chrSel <- sample(x=seq_len(length(listChr)), 1)
            newChr[[i]] <- processChr(curSample=curSample,
                                        dfChr=resChr[[listChr[chrSel]]][[i]],
                                        chrCur=chr)
        }
        df[[chr]] <- do.call(rbind, newChr)
    }

    df <- do.call(rbind, df)
    rownames(df) <- seq_len(nrow(df))
    return(df)
}

