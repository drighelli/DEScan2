#' This function calls peaks from bed or bam inputs using a variable window scan
#' with a poisson model using the surrounding maxCompWinWidth (10kb) as background.
#'
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
#' @param genomeName the code of the genome to use as reference for the input
#'                   files. (cfr. constructBedRanges function parameters)
#' @param onlyStdChrs a flag to work only with standard chromosomes
#'                    (cfr. constructBedRanges function parameters)
#' @param chr if not NULL, a character like "chr#" indicating the
#'            chromosomes to use
#'
#' @return A GRangesList where each element is a sample.
#'         each GRanges represents the founded peaks and attached the z-score
#'         of the peak as mcols.
#' @export
#' @examples
#'
findPeaks <- function(files, filetype=c("bam", "bed"),
                      genomeName=NULL,
                      binSize=50, minWin=1, maxWin=20,
                      zthresh=5, minCount=0.1,
                      minCompWinWidth=5000,
                      maxCompWinWidth=10000,
                      outputName="Peaks", save=TRUE, verbose=FALSE,
                      fragmentLength=200,
                      onlyStdChrs=TRUE,
                      chr=NULL
                      )
{

    if(!is.null(chr) && length(grep(pattern="chr", chr))!=length(chr))
        stop("Insert valid chr(s), use the \"chr#\" form!")
    if(!save) warning("Save is false, not saving results!")
    if(length(files) == 0)
        stop("You have to provide one or more input files!\nExiting.")
    filetype <- match.arg(filetype)

    fileGRangesList <- NULL
    winVector <- c(minWin:maxWin)


    for (fileIdx in 1:length(files))
    {
        file <- files[[fileIdx]]
        bedGRanges <- constructBedRanges(filename=file, filetype=filetype,
                                         genomeName=genomeName,
                                         onlyStdChrs=onlyStdChrs)

        bedGrangesChrsList <- cutGRangesPerChromosome(bedGRanges)
        if(!is.null(chr))
            bedGrangesChrsList <- keepRelevantChrs(bedGrangesChrsList, chr)

        chrZRangesList <- GenomicRanges::GRangesList(
            lapply(bedGrangesChrsList, function(chrGRanges) { ########## to parallelize on chromosomes

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
            new_s <- get_disjoint_max_win(z0=z,#[1:5000,], ###############################
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

            return(chrZRanges) ## one for each chromosome
            })
        )
        names(chrZRangesList) <- names(bedGrangesChrsList)
        if(save)
        {
            zGRangesToSave <- unlist(chrZRangesList)
            filename <- paste0(file, "_zt", zthresh, "_mnw", minWin,
                               "_mxw", maxWin, "_bin", binSize) #######################
            saveGRangesAsBed(GRanges=zGRangesToSave, filepath=outputName,
                             filename=filename) ########################################
        }

        fileGRangesList <- c(fileGRangesList, chrZRangesList)
    }
    names(fileGRangesList) <- files

    return(fileGRangesList)
}

#' computeZ: Computes Z-Scores returning the z matrix
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


#' get_disjoint_max_win find significant z score windows keeping the max value
#' without intersections
#'
#' @param z0 Matrix containing z scores with bins as rows and windows size as
#'     columns
#' @param sigwin Integer indicating how many bins per fragment
#' @param nmax Integer indicating the maximum number of windows to return
#' @param zthresh Integer indicating the minimum z-score considered significant
#' @param two_sided not used argument
#' @param verbose verbose flag
#' @return a matrix of integer containing founded peaks
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

        t <- inds %% nrow(z0) # determine column

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

#' computeLambdaOnChr: computes the lambdas on a chromosome for the
#' winVector windows and other two windows (min/maxCompWinWidth) to compare with
#'
#' @param chrGRanges the GRanges representing the reads of the chromosome
#' @param winVector the of width of the windows used to compute the coverage
#' @param minChrRleWComp and Rle object within coverage of window of width
#'                       minCompWinWidth
#' @param minCompWinWidth the width of the window used for the coverage of
#'                        minChrRleWComp
#' @param maxChrRleWComp and Rle object within coverage of window of width
#'                       minCompWinWidth
#' @param maxCompWinWidth the width of the window used for the coverage of
#'                        maxChrRleWComp
#'
#' @return an RleList where each element is a window of winVector, within an Rle
#'         representing the lambda computed for that window
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


#' computeCoverageMovingWindowOnChr: computes the coverage on a chromosomewith a
#' set of moving windows of dimensions minWinWidth:maxWinWidth
#'
#' @param chrBedGRanges a GRanges to compute the coverage
#' @param minWinWidth the minimum width of the window to use for the coverage
#' @param maxWinWidth the maximum width of the window to use for the coverage
#' @param binWidth the dimension of the bin in base number
#'
#' @return RleList where each element is a window within the Rle of its coverage
#' @export
#'
#' @examples
#'
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
    message("Computing coverage on Chromosome ",
            chrBedGRanges@seqnames@values,
            " binned by ", binWidth, " bin dimension")
    chrCovRle <- binnedSumOnly(bins=binnedChromosome,
                               numvar=chrCoverage,
                               mcolname="bin_cov")
    wins <- minWinWidth:maxWinWidth
    runWinRleList <- IRanges::RleList(
                        lapply(wins, function(win) {
                            message("Running window ", win, " of ", maxWinWidth)
                            evenRunSum(x=chrCovRle, k=win, endrule="constant")
                        })
                    )
    return(runWinRleList)
}

#' binToChrCoordMatRowNames: computes the starting range of the bins for the
#' binMatrix, taking in input the length of the chromosome of the matrix
#'
#' @param binMatrix a matrix where each row represents a bin
#' @param chrLength the length of the chromosome of the binMatrix
#' @param binWidth the width of the bin
#'
#' @return the binMatrix with start range as rownames
#'
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

#' binnedSum: this function computes the coverage over a binned chromosome,
#' starting from a per base computed coverage
#' credits: http://crazyhottommy.blogspot.com/2016/02/compute-averagessums-on-granges-or.html
#'
#' @param bins a GRanges object representing a chromosome binned
#' @param numvar an RleList representing the per base coverage over the chr
#' @param mcolname the name of column where the sum have to be stored
#'
#' @return the bins GRanges with the mcolname attached
#' @export
#'
#' @examples
#' ## dividing one chromosome in bins of 50 bp each
#' seqinfo <- GenomeInfoDb::Seqinfo("mm9")
#' bins <- GenomicRanges::tileGenome(seqlengths=seqinfo@chr1,
#'                                              tilewidth=50,
#'                                              cut.last.tile.in.chrom=TRUE)
#' gr <- GRanges(seqnames = Rle("Chr1", 100),
#' ranges =IRanges(start = seq(from=10, to=1000, by=10),
#' end=seq(from=20, to=1010, by = 10)))
#' ## computing coverage per single base on granges
#' cov <- GenomicRanges::coverage(x=gr)
#'
#' binnedCovGR <- binnedSum(bins, cov, "binned_cov")
#'
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

#' binnedSumOnly: it's useful just to coerce the bin coverage to an Rle object
#'
#' @param bins a GRanges object representing a chromosome binned
#' @param numvar an RleList representing the per base coverage over the chr
#' @param mcolname the name of column where the sum have to be stored
#'
#' @return an Rle within the per bin computed coverage
binnedSumOnly <- function(bins, numvar, mcolname)
{
    binsGRanges <- binnedSum(bins=bins, numvar=numvar, mcolname=mcolname)
    ## coercing just the binned coverage to Rle
    chrCovRle <- as(mcols(binsGRanges)[[mcolname]],"Rle")
    return(chrCovRle)
}

#' evenRunSum: this function computes a running sum over x with a window width k
#' (modified from S4Vectors package to work on even on even k, in such a case
#' it adds a length at the end of the output Rle)
#'
#' @param x an Rle object, typically a coverage object
#' @param k window dimension for the running sum over x
#' @param endrule refer to S4Vectors::runSum
#' @param na.rm refer to S4Vectors::runSum
#'
#' @return an Rle within the running sum over x with a win o length k
evenRunSum <- function(x, k, endrule = c("drop", "constant"), na.rm = FALSE)
{
    stopifnot(is(x, "Rle"))
    endrule <- match.arg(endrule)
    n <- length(x)
    # k <- normArgRunK(k = k, n = n, endrule = endrule)
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
