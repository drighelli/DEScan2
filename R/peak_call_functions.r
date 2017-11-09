#' This function calls peaks from bed or bam inputs using a variable window scan
#' with a poisson model using the surrounding 10kb as background.
#' files, filetype="bam", chr=1:19, fragmentLength=200,
#' @param files Character vector containing paths of files to be analyzed.
#' @param filetype Character, either "bam" or "bed" indicating format of input
#'                 file.
#' @param fragmentLength Integer indicating the average DNA fragment length bp.
#' @param minWin Integer indicating the minimum window size in units of 50 bp,
#'               i.e., min_win=2 resulting in a 100 bp window.
#' @param binSize Integer size in base pairs of the minimum window for scanning,
#'                50 is the default.
#' @param maxWin Integer indicating the maximum allowed window size in units of
#'               50 bp.
#' @param minCompWinWidth minimum bases width of a comparing window for Z-score
#' @param maxCompWinWidth maximum bases width of a comparing window for Z-score
#' @param zthresh Cuttoff value for z-scores. Only windows with greater z-scores
#'                will be retained.
#' @param minCount A small constant (usually no larger than one) to be added to
#'                 the counts prior to the log transformation to avoid
#'                 problems with log(0).
#' @param outputName A string, Name of the folder to save the Peaks (optional),
#'                   if the directory doesn't exist, it will be created.
#'                   default is "Peaks"
#' @param save Boolean, if TRUE files will be saved in a "./Peaks/chr*"
#'             directory created (if not already present) in the current
#'             working directory.
#' @return a matrix with peaks as rows and 4 columns describing the genomic
#'         coordinates (chr, start, end) as well as the associated z-score.
#' @export
#' @importFrom utils read.table write.table
#' @examples
#'
findPeaks <- function(files, filetype=c("bam", "bed"),
                      genomeName=NULL,
                      binSize=50, minWin=1, maxWin=20,
                      zthresh=5, minCount=0.1,
                      minCompWinWidth=5000,
                      maxCompWinWidth=10000,
                      outputName="Peaks", save=TRUE, verbose=FALSE,
                      fragmentLength=200
                      )
{

    if(length(files) == 0) {
        stop("You have to provide one or more input files!\nExiting.")
    }

    filetype <- match.arg(filetype)

    winVector <- c(minWin:maxWin)

    for (file in files) {

        bedGRanges <- constructBedRanges(filename=file, filetype=filetype,
                                         genomeName=genomeName)

        bedGrangesChrsList <- cutGRangesPerChromosome(bedGRanges)
        # chrGRanges <- bedGrangesChrsList[[1]]
        # bedGrangesChrsList<-bedGrangesChrsList[c(1,2)]
        # chrGRanges=bedGrangesChrsList[[2]]
        chrZRangesList <- GenomicRanges::GRangesList(
            lapply(bedGrangesChrsList, function(chrGRanges) { ## to parallelize

            runWinRleList <- computeCoverageMovingWindowOnChr(
                                                       chrBedGRanges=chrGRanges,
                                                       minWinWidth=minWin,
                                                       maxWinWidth=maxWin,
                                                       binWidth=binSize
                                                       )
            minCompRunWinRleList <- computeCoverageMovingWindowOnChr(
                                                    chrBedGRanges=chrGRanges,
                                                    minWinWidth=minCompWinWidth,
                                                    maxWinWidth=minCompWinWidth,
                                                    binWidth=binSize
                                                    )
            maxCompRunWinRleList <- computeCoverageMovingWindowOnChr(
                                                    chrBedGRanges=chrGRanges,
                                                    minWinWidth=maxCompWinWidth,
                                                    maxWinWidth=maxCompWinWidth,
                                                    binWidth=binSize
                                                    )
            lambdaChrRleList <- computeLambdaOnChr(
                                       chrGRanges=chrGRanges,
                                       winVector=winVector,
                                       minChrRleWComp=minCompRunWinRleList[[1]],
                                       minCompWinWidth=minCompWinWidth,
                                       maxChrRleWComp=maxCompRunWinRleList[[1]],
                                       maxCompWinWidth=maxCompWinWidth
                                       )
            z <- computeZ(lambdaChrRleList=lambdaChrRleList,
                          runWinRleList=runWinRleList,
                          chrLength=chrGRanges@seqinfo@seqlengths,
                          minCount=minCount, binSize=50
                          )
            new_s <- get_disjoint_max_win(z0=z[1:5000,], ###############################
                                          sigwin=fragmentLength/binSize,
                                          nmax=Inf, zthresh=zthresh,
                                          two_sided=FALSE, verbose=FALSE
                                          )
            chrZRanges <- createGranges(chrSeqInfo=chrGRanges@seqinfo,
                                      starts=as.numeric(rownames(z)[new_s[,1]]),
                                      widths=new_s[,2]*binSize,
                                      mcolname="z-score",
                                      mcolvalues=new_s[,3]
                                      )

            return(chrZRanges)
            })
        )

        zGRanges <- unlist(chrZRangesList)
        saveGRangesAsBed(GRanges=zGRanges, filepath=outputName, filename=file)
        fileGRangesList <- c(fileGRangesList, zGRanges)
    }
    names(fileGRangesList) <- files

    return(fileGRangesList)
}


#' Title
#'
#' @param chrSeqInfo
#' @param starts
#' @param widths
#' @param mcolname
#' @param mcolvalues
#'
#' @return
#' @export
#'
#' @examples
createGranges <- function(chrSeqInfo, starts, widths,
                          mcolname=NULL, mcolvalues=NULL) {

    gr <- GRanges(seqnames=as.character(chrSeqInfo@seqnames),
            ranges=IRanges(start=starts, width=widths),
            seqinfo=chrSeqInfo
          )

    if(!is.null(mcolname) )
    {
        if(!is.null(mcolvalues)
           &&
           (length(gr@ranges@start) == length(mcolvalues))
           )
        {
            mcols(gr)[[mcolname]] <- mcolvalues
        }
        else
        {
            warning("Cannot set mcolvalues! Vector length not matching Ranges")
        }
    }
    return(gr)
}

#' Computes Z-Scores returning the z matrix
#'
#' @param lambdaChrRleList an RleList of lambda values computed
#'                         by computeLambdaOnChr function
#'                         each element of the list is an Rle representing the
#'                         lambda for the moving window in the list position
#' @param runWinRleList an RleList of coverage values computed
#'                         by computeCoverageMovingWindowOnChr function
#'                         each element of the list is an Rle representing the
#'                         coverage for the moving window in the list position
#' @param chrLength the length of the chr in analysis
#' @param minCount A small constant (usually no larger than one) to be added to
#'                 the counts prior to the log transformation to avoid problems
#'                 with log(0).
#' @param binSize the size of the bin
#'
#' @return z a matrix of z scores for each window (column) and bin (row).
#'         where the rownames represent the starting base of each bin
#' @keywords internal
computeZ <- function(lambdaChrRleList, runWinRleList, chrLength,
                     minCount=0.1, binSize=50)
{
    # runWinRleM <- RleListToRleMatrix(runWinRleList)
    # lambdaChrRleM <- RleListToRleMatrix(lambdaChrRleList)

    lambdaChrRleMm <- matrix(unlist(lambdaChrRleList),
                             ncol=20, byrow=TRUE)

    runWinRleMm <- matrix(unlist(runWinRleList), ncol=20, byrow=TRUE)

    message("Computing Z-Score")
    z <- sqrt(2) * sign(runWinRleMm - lambdaChrRleMm) *
        sqrt(runWinRleMm *
                 log(pmax(runWinRleMm, minCount) / lambdaChrRleMm) -
                 (runWinRleMm - lambdaChrRleMm)
        )

    z <- binToChrCoordMatRowNames(binMatrix=z,
                                  chrLength=chrLength,
                                  binWidth=binSize)
    return(z)
}


#' find significant z score windows keeping the max value without intersections
#'
#' @param z0 Matrix containing z scores with bins as rows and windows size as
#'     columns
#' @param sigwin Integer indicating how many bins per fragment
#' @param nmax Integer indicating the maximum number of windows to return
#' @param zthresh Integer indicating the minimum z-score considered significant
#' @param two_sided #######
#' @param verbose
#' @return s. a matrix of integer containing founded peaks
#' @keywords internal
get_disjoint_max_win <- function(z0, sigwin=20, nmax=Inf,
                                 zthresh=-Inf, two_sided=FALSE, verbose=FALSE)
{
    message("Computing peaks...")

    s <- matrix(ncol=3, nrow=0)

    maxwin <- ncol(z0)

    if (two_sided) {
        z0 <- abs(z0)
    }

    while (TRUE) {
        inds <- which.max(z0) # find max z
        if (length(inds) == 0) break
        w <- ceiling(inds / nrow(z0)) # determine row

        t <- inds %% nrow(z0) # determine column ################

        if (t == 0) t <- nrow(z0)
        ## break loop once as max z below thresh
        if (z0[t, w] < zthresh) break

        s <- rbind(s, c(t, w, z0[t, w]))

        if (verbose) {
            cat("Maximizing window: ", t, ",", w, " Score=", z0[t, w], "\n")
        }
        st <- max(1, t - sigwin - maxwin + 1)
        ed <- min(t + w + sigwin - 1, nrow(z0))
        z0[st:ed, ] <- -Inf

        if (nrow(s) >= nmax) break
    }
    colnames(s) <- c("bin", "window", "z")

    message("...done!")
    return(s)
}



#' Title
#'
#' @param RleList
#' @param dimnames
#'
#' @return
#' @export
#'
#' @examples
RleListToRleMatrix <- function(RleList, dimnames=NULL)
{
    if(!is.null(dimnames)) {
        rlem <- DelayedArray::RleArray(rle=unlist(RleList),
                              dim=c(length(RleList[[1]]), length(RleList)),
                              dimnames=dimnames
                              )
    } else {
        rlem <- DelayedArray::RleArray(rle=unlist(RleList),
                              dim=c(length(RleList[[1]]), length(RleList)))
    }
    return(rlem)

}

#' Title
#'
#' @param chrGRanges
#' @param winVector
#' @param minChrRleWComp
#' @param minCompWinWidth
#' @param maxChrRleWComp
#' @param maxCompWinWidth
#'
#' @return
#' @export
#'
#' @examples
computeLambdaOnChr <- function(chrGRanges,
                                winVector=c(1:20),
                                minChrRleWComp,
                                minCompWinWidth=5000,
                                maxChrRleWComp,
                                maxCompWinWidth=10000)
{
    message("Computing lambdas")
    chrTotRds <- length(chrGRanges@ranges@start)
    chrTotBases <- (chrGRanges@ranges@start+chrGRanges@ranges@width)[chrTotRds]
                    - chrGRanges@ranges@start[1]

    chrLamb <- chrTotRds %*% t(winVector) / chrTotBases ## shouldn't it be for the size of bin? (bases size of a window)

    lamblMat <- rbind(chrLamb)[rep(1,length(minChrRleWComp)), ]

    lamblRleList <- IRanges::RleList(apply(X=lamblMat, MARGIN=2,
                                           function(x){as(x,"Rle")}))

    chrLamMinWComp <- as.vector(minChrRleWComp) %*% t(winVector)/minCompWinWidth
    chrLamMinWCompRleList <- apply(chrLamMinWComp, 2, S4Vectors::Rle)


    maxChrLamWComp <- as.vector(maxChrRleWComp) %*% t(winVector)/maxCompWinWidth
    maxChrLamWCompRleList <- apply(maxChrLamWComp, 2, S4Vectors::Rle)


    lamlocRleList <-  IRanges::RleList(lapply(winVector, function(idx) {
                                              pmax(maxChrLamWCompRleList[[idx]],
                                                   chrLamMinWCompRleList[[idx]],
                                                   lamblRleList[[idx]])
                                              })
                                       )
    return(lamlocRleList)

}


#' Title
#'
#' @param chrBedGRanges
#' @param minWinWidth
#' @param maxWinWidth
#' @param binWidth
#'
#' @return
#' @export
#'
#' @examples
computeCoverageMovingWindowOnChr <- function(chrBedGRanges, minWinWidth=1,
                                             maxWinWidth=20, binWidth=50)
{

    ## dividing chromosome in bins of binWidth dimention each
    binnedChromosome<-GenomicRanges::tileGenome(seqlengths=chrBedGRanges@seqinfo
                                                , tilewidth=binWidth
                                                , cut.last.tile.in.chrom=TRUE)

    ## computing coverage per single base on bed
    chrCoverage <- GenomicRanges::coverage(x=chrBedGRanges)

    ## computing coverage per each bin on chromosome
    message("Computing coverage on Chromosome ", chrBedGRanges@seqnames@values,
            " binned by ", binWidth, " bin dimension")

    chrCovRle <- binnedSumOnly(bins=binnedChromosome,
                                    numvar=chrCoverage,
                                    mcolname="bin_cov")


    wins <- minWinWidth:maxWinWidth
    runWinRleList <- IRanges::RleList(lapply(wins, function(win) {
        message("Running window ", win, " of ", maxWinWidth)
        oddRunSum(chrCovRle, k=win, endrule="constant")
    }))

    return(runWinRleList)
}

#' Title
#'
#' @param binMatrix
#' @param chrLength
#' @param binWidth
#'
#' @return
#' @export
#'
#' @examples
binToChrCoordMatRowNames <- function(binMatrix, chrLength, binWidth=50)
{

    ## computing bin in base ranges to add as rownames
    endBinRanges <- seq(from=binWidth-1, to=chrLength, by=binWidth)

    if(chrLength != endBinRanges[length(endBinRanges)])
    {
        ## it can only be lesser than the length of the chromosome
        ## adding the last missing bin
        endBinRanges <- c(endBinRanges, (chrLength-1))
    }
    ## computing the start of regions
    startBinRanges <- endBinRanges-(binWidth-1)
    if(dim(binMatrix)[1] != length(startBinRanges))
    {
        stop("something went wrong! matrix row dimension
             different than expected")
    }
    rownames(binMatrix) <- startBinRanges
    return(binMatrix)
}


#' Title
#'
#' @param bedGRanges
#'
#' @return
#' @export
#'
#' @examples
cutGRangesPerChromosome <- function(bedGRanges)
{
    interestedChrs <- bedGRanges@seqinfo@seqnames

    ##     bed19 <- keepSeqlevels(bed, "chr19")
    bedGRList <- lapply(interestedChrs, function(x)
    {
        bgr <- bedGRanges[which(bedGRanges@seqnames %in% x),]
        bgr@seqinfo <- bedGRanges@seqinfo[x]
        seqnames(bgr) <- droplevels(seqnames(bgr))
        bgr
    })
    names(bedGRList) <- interestedChrs
    ## intentionally left commented, GRangesList reconstruct the entire seqinfo,
    ## while we want it cutted per chromosomes
    # bedGRList <- GRangesList(bedGRList)

    return(bedGRList)
}




#' Title
#' http://crazyhottommy.blogspot.com/2016/02/compute-averagessums-on-granges-or.html
#' @param bins
#' @param numvar
#' @param mcolname
#'
#' @return
#' @export
#'
#' @examples
binnedSum <- function(bins, numvar, mcolname)
{
    stopifnot(is(bins, "GRanges"))
    stopifnot(is(numvar, "RleList"))
    stopifnot(identical(seqlevels(bins), names(numvar)))
    bins_per_chrom <- split(ranges(bins), seqnames(bins))
    sums_list <- lapply(names(numvar),
                        function(seqname) {
                            views <- Views(numvar[[seqname]],
                                           bins_per_chrom[[seqname]])
                            viewSums(views)
                        })
    new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
    mcols(bins)[[mcolname]] <- new_mcol

    return(bins)
}

#' Title
#'
#' @param bins
#' @param numvar
#' @param mcolname
#'
#' @return
#' @export
#'
#' @examples
binnedSumOnly <- function(bins, numvar, mcolname)
{
    binsGRanges <- binnedSum(bins=bins, numvar=numvar, mcolname=mcolname)
    ## coercing just the binned coverage to Rle
    chrCovRle <- as(mcols(binsGRanges)[[mcolname]],"Rle")
    return(chrCovRle)
}

#' Title
#'
#' @param x
#' @param k
#' @param endrule
#' @param na.rm
#'
#' @return
#' @export
#'
#' @examples
oddRunSum<-function(x, k, endrule = c("drop", "constant"), na.rm = FALSE)
{
    # endrule <- match.arg(endrule)
    n <- length(x)
    # k <- normargRunK(k = k, n = n, endrule = endrule)
    ans <- S4Vectors::.Call2("Rle_runsum", x, as.integer(k), as.logical(na.rm),
                             PACKAGE="S4Vectors")
    if (endrule == "constant") {
        j <- (k + 1L) %/% 2L

        S4Vectors::runLength(ans)[1L] <- S4Vectors::runLength(ans)[1L]+(j - 1L)
        if( (k %% 2L) == 0 ) j=j+1
        S4Vectors::runLength(ans)[S4Vectors::nrun(ans)] <-
                        S4Vectors::runLength(ans)[S4Vectors::nrun(ans)]+(j - 1L)
    }
    return(ans)
}



