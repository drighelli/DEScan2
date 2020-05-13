#' findPeaks
#' @description This function calls peaks from bed or bam inputs using a
#' variable window scan with a poisson model using the surrounding
#' maxCompWinWidth (10kb) as background.
#'
#' @param files Character vector containing paths of files to be analyzed.
#' @param filetype Character, either "bam" or "bed" indicating format of input
#' file.
#' @param minWin Integer indicating the minimum window size in bases notation.
#' @param maxWin Integer indicating the maximum window size in bases notation.
#' @param binSize Integer size in bases of the minimum window for scanning,
#' 50 is the default.
#' @param minCompWinWidth minimum bases width of a comparing window for Z-score.
#' @param maxCompWinWidth maximum bases width of a comparing window for Z-score.
#' @param zthresh Cuttoff value for z-scores. Only windows with greater z-scores
#' will be kept, default is 10.
#' @param minCount A small constant (usually no larger than one) to be added to
#' the counts prior to the log transformation to avoid problems with log(0).
#' @param outputFolder A string, Name of the folder to save the Peaks (optional)
#' if the directory doesn't exist, it will be created. (Default is "Peaks")
#' @param save Boolean, if TRUE files will be saved in a "./Peaks/chr*"
#' directory created (if not already present) in the current working directory.
#' @param force a boolean flag indicating if to force output overwriting.
#' @param genomeName the code of the genome to use as reference for the input
#' files. (cfr. constructBedRanges function parameters)
#' @param onlyStdChrs a flag to work only with standard chromosomes.
#' (cfr. constructBedRanges function parameters).
#' @param chr if not NULL, a character like "chr#" indicating the chromosomes
#' to use.
#' @param sigwin an integer value used to compute the length of the signal
#' of a peak (default value is 10).
#' @param verbose if to show additional messages
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#'
#' @return A GRangesList where each element is a sample.
#' Each GRanges represents the founded peaks and attached the z-score of
#' the peak as mcols.
#' @export
#'
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomeInfoDb seqinfo seqlengths
#' @importFrom BiocParallel bpparam bplapply
#' @examples
#' bam.files <- list.files(system.file("extdata/bam", package = "DEScan2"),
#'                         full.names = TRUE)
#'
#' peaks <- findPeaks(files=bam.files[1], filetype="bam",
#'                         genomeName="mm9",
#'                         binSize=50, minWin=50, maxWin=1000,
#'                         zthresh=5, minCount=0.1, sigwin=10,
#'                         minCompWinWidth=5000, maxCompWinWidth=10000,
#'                         save=FALSE,
#'                         onlyStdChrs=TRUE,
#'                         chr=NULL,
#'                         verbose=FALSE)
#' head(peaks)
findPeaks <- function(files, filetype=c("bam", "bed"),
                        genomeName=NULL,
                        binSize=50, minWin=50, maxWin=1000,
                        zthresh=10, minCount=0.1,
                        minCompWinWidth=5000,
                        maxCompWinWidth=10000,
                        outputFolder="Peaks", save=TRUE, force=TRUE,
                        verbose=FALSE,
                        sigwin=10,
                        onlyStdChrs=TRUE,
                        chr=NULL,
                        BPPARAM=BiocParallel::bpparam())
{

    #if(!is.null(chr) && length(grep(pattern="chr", chr))!=length(chr))
    #    stop("Insert valid chr(s), use the \"chr#\" form!")
    if(!save) message("Save is false, not saving results!\n")
    if(length(files) == 0)
        stop("You have to provide one or more input files!\nExiting.")
    stopifnot((minWin %% binSize) == 0)
    stopifnot((maxWin %% binSize) == 0)
    stopifnot((minCompWinWidth %% binSize) == 0)
    stopifnot((maxCompWinWidth %% binSize) == 0)
    filetype <- match.arg(filetype)

    fileGRangesList <- NULL
    winVector <- c((minWin/binSize):(maxWin/binSize))

    for (fileIdx in seq_len(length(files)))
    {
        file <- files[[fileIdx]]
        bedGRanges <- constructBedRanges(filename=file, filetype=filetype,
                                            genomeName=genomeName,
                                            onlyStdChrs=onlyStdChrs,
                                            verbose=verbose)

        bedGrangesChrsList <- cutGRangesPerChromosome(bedGRanges)
        if(!is.null(chr))
            bedGrangesChrsList <- keepRelevantChrs(bedGrangesChrsList, chr)

        if(verbose) message("Calling Peaks on chromosomes...")

        chrZRangesList <- GenomicRanges::GRangesList(
            BiocParallel::bplapply(bedGrangesChrsList, function(chrGRanges)
            {

            runWinRleList <- computeCoverageMovingWindowOnChr(
                                                chrBedGRanges=chrGRanges,
                                                minWinWidth=minWin,
                                                maxWinWidth=maxWin,
                                                binWidth=binSize,
                                                verbose=verbose)
            minCompRunWinRleList <- computeCoverageMovingWindowOnChr(
                                                chrBedGRanges=chrGRanges,
                                                minWinWidth=minCompWinWidth,
                                                maxWinWidth=minCompWinWidth,
                                                binWidth=binSize,
                                                verbose=verbose)
            maxCompRunWinRleList <- computeCoverageMovingWindowOnChr(
                                                chrBedGRanges=chrGRanges,
                                                minWinWidth=maxCompWinWidth,
                                                maxWinWidth=maxCompWinWidth,
                                                binWidth=binSize,
                                                verbose=verbose)
            lambdaChrRleList <- computeLambdaOnChr(
                                    chrGRanges=chrGRanges,
                                    winVector=winVector,
                                    minChrRleWComp=minCompRunWinRleList[[1]],
                                    minCompWinWidth=minCompWinWidth,
                                    maxChrRleWComp=maxCompRunWinRleList[[1]],
                                    maxCompWinWidth=maxCompWinWidth,
                                    verbose=verbose)
            seqlength <- GenomeInfoDb::seqlengths(
                                        GenomeInfoDb::seqinfo(chrGRanges))[[1]]
            Z <- computeZ(lambdaChrRleList=lambdaChrRleList,
                            runWinRleList=runWinRleList,
                            chrLength=seqlength,
                            minCount=minCount, binSize=binSize,
                            verbose=verbose)
            newS <- c_get_disjoint_max_win(z0=Z,
                                    sigwin=sigwin,
                                    zthresh=zthresh,
                                    verbose=verbose)
            chrZRanges <- createGranges(
                                chrSeqInfo=GenomeInfoDb::seqinfo(chrGRanges),
                                starts=as.numeric(rownames(Z)[newS[,1]]),
                                widths=newS[,2]*binSize,
                                mcolname="z-score",
                                mcolvalues=newS[,3])

            return(chrZRanges) ## one for each chromosome
            }, BPPARAM=BPPARAM)
        )
        # names(chrZRangesList) <- names(bedGrangesChrsList)
        ZRanges <- unlist(chrZRangesList)
        # ZRanges <- sort(ZRanges)
        if(save)
        {
            # zGRangesToSave <- unlist(chrZRangesList)
            filename <- paste0(basename(file), "_zt", zthresh, "_mnw", minWin,
                            "_mxw", maxWin, "_sw", sigwin, "_bin", binSize)
            saveGRangesAsBed(GRanges=ZRanges, filepath=outputFolder,
                            filename=filename, verbose=verbose,
                            force=force)
        }
        fileGRangesList <- c(fileGRangesList, ZRanges)
    }
    names(fileGRangesList) <- files

    return(GenomicRanges::GRangesList(fileGRangesList))
}

#' computeZ
#' @description Computes Z-Scores returning the z matrix.
#'
#' @param lambdaChrRleList an RleList of lambda values computed by
#' computeLambdaOnChr function each element of the list is an Rle representing
#' the lambda for the moving window in the list position.
#' @param runWinRleList an RleList of coverage values computed.
#' by computeCoverageMovingWindowOnChr function each element of the list is an
#' Rle representing the coverage for the moving window in the list position.
#' @param chrLength the length of the chr in analysis.
#' @param minCount A small constant (usually no larger than one) to be added to
#' the counts prior to the log transformation to avoid problems with log(0).
#' @param binSize the size of the bin.
#' @param verbose verbose output.
#'
#' @return z a matrix of z scores for each window (column) and bin (row).
#' where the rownames represent the starting base of each bin.
#' @keywords internal.
computeZ <- function(lambdaChrRleList, runWinRleList, chrLength,
                        minCount=0.1, binSize=50, verbose=FALSE)
{
    # runWinRleM <- RleListToRleMatrix(runWinRleList)
    # lambdaChrRleM <- RleListToRleMatrix(lambdaChrRleList)

    # lambdaChrRleMm <- matrix(unlist(lambdaChrRleList),
    #                          ncol=20, byrow=TRUE)

    lambdaChrRleMm <- matrix(unlist(lambdaChrRleList),
                            ncol=length(lambdaChrRleList), byrow=FALSE)

    # runWinRleMm <- matrix(unlist(runWinRleList), ncol=20, byrow=TRUE)
    runWinRleMm <- matrix(unlist(runWinRleList),
                        ncol=length(runWinRleList), byrow=FALSE)

    if(verbose) message("Computing Z-Score")
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


#' c_get_disjoint_max_win
#' @description just a wrapper for the C function. Useful to modify indexes and
#' colnames.
#'
#' @param z0 the z matrix.
#' @param sigwin the sigwin.
#' @param nmax the nmax.
#' @param zthresh peaks lower than this value will not be kept.
#' @param verbose verbose flag.
#'
#' @return a matrix
#' @keywords internal
c_get_disjoint_max_win <- function(z0, sigwin=10, nmax=9999999,
                                    zthresh=10, verbose=FALSE)
{
    m <- rcpparma_get_disjoint_max_win(z0=z0,
                                        sigwin=sigwin,
                                        zthresh = zthresh,
                                        nmax=nmax,
                                        verbose=verbose)
    ## changing indexes to R style
    m[,1] <- m[,1]+1
    m[,2] <- m[,2]+1
    colnames(m) <- c("bin", "window", "z")
    return(m)
}

#' get_disjoint_max_win
#' @description find significant z score windows keeping the max value
#'                without intersections
#'
#' @param z0 Matrix containing z scores with bins as rows and windows size as
#'            columns.
#' @param sigwin Integer indicating how many bins per fragment.
#' @param nmax Integer indicating the maximum number of windows to return.
#' @param zthresh Integer indicating the minimum z-score considered significant.
#' @param two_sided not used argument.
#' @param verbose verbose flag.
#' @return a matrix of integer containing founded peaks
#' @keywords internal
get_disjoint_max_win <- function(z0, sigwin=20, nmax=Inf,
                                zthresh=-Inf, two_sided=FALSE, verbose=FALSE)
{
    if(verbose) message("Computing peaks...")

    s <- matrix(ncol=3, nrow=0)

    maxwin <- ncol(z0)

    if (two_sided) {
        z0 <- abs(z0)
    }
    i=1
    while (TRUE) {
        inds <- which.max(z0) # find max z
        if (length(inds) == 0) break
        w <- ceiling(inds / nrow(z0)) # determine row

        t <- inds %% nrow(z0) # determine column

        if (t == 0) t <- nrow(z0)
        ## break loop once as max z below thresh
        if (z0[t, w] < zthresh) break

        s <- rbind(s, c(t, w, z0[t, w]))
        if(verbose)
        {
            if((i %% 100) == 0)
            {
                message("Maximizing window: ", t, ",", w,
                        " Score=", z0[t, w], "\n")
            }
        }

        i=i+1


        st <- max(1, t - sigwin - maxwin + 1)
        ed <- min(t + w + sigwin - 1, nrow(z0))
        # message("st: ", st, " ed: ", ed, "\n")
        z0[st:ed, ] <- -Inf

        if (nrow(s) >= nmax) break
    }
    colnames(s) <- c("bin", "window", "z")

    if(verbose) message("...done!")
    return(s)
}

#' computeLambdaOnChr
#' @description computes the lambdas on a chromosome for the
#' winVector windows and other two windows (min/maxCompWinWidth) to compare with
#'
#' @param chrGRanges the GRanges representing the reads of the chromosome.
#' @param winVector the of width of the windows used to compute the coverage.
#' @param minChrRleWComp and Rle object within coverage of window of width
#'                        minCompWinWidth.
#' @param minCompWinWidth the width of the window used for the coverage of
#'                      minChrRleWComp in bases.
#' @param maxChrRleWComp and Rle object within coverage of window of width
#'                      minCompWinWidth.
#' @param maxCompWinWidth the width of the window used for the coverage of
#'                        maxChrRleWComp in bases.
#' @param binSize the size of the bin in bases.
#' @param verbose verbose flag.
#' @return an RleList where each element is a window of winVector, within an Rle
#'         representing the lambda computed for that window.
#'
#' @importFrom IRanges RleList
#' @importFrom S4Vectors Rle
#' @importFrom GenomicRanges end start
#'
#' @keywords internal
computeLambdaOnChr <- function(chrGRanges,
                                winVector=seq_len(20),
                                minChrRleWComp,
                                minCompWinWidth=5000,
                                maxChrRleWComp,
                                maxCompWinWidth=10000,
                                verbose=TRUE)
{
    if(verbose) message("Computing lambdas")
    chrTotRds <- length(chrGRanges)

    chrTotBases <- GenomicRanges::end(chrGRanges)[chrTotRds] -
                                            GenomicRanges::start(chrGRanges)[1]

    chrLamb <- chrTotRds %*% t(winVector) / chrTotBases

    lamblMat <- rbind(chrLamb)[rep(1,length(minChrRleWComp)), ]

    # lamblRleList <- IRanges::RleList(apply(X=lamblMat, MARGIN=2,
    #                                        function(x){as(x,"Rle")}))

    minchrLamWComp <- as.vector(minChrRleWComp) %*%
                                                t(winVector)/(minCompWinWidth)
    # minchrLamWCompRleList <- apply(minchrLamWComp, 2, S4Vectors::Rle)


    maxChrLamWComp <- as.vector(maxChrRleWComp) %*%
                                                t(winVector)/(maxCompWinWidth)
    # maxChrLamWCompRleList <- apply(maxChrLamWComp, 2, S4Vectors::Rle)

    # lamlocRleList <-  IRanges::RleList(lapply(winVector, function(win) {
    #                                       pmax(maxChrLamWCompRleList[[win]],
    #                                       minchrLamWCompRleList[[win]],
    #                                       chrLamb[win])
    #                                   }))

    lamlocRleList <-  IRanges::RleList(lapply(winVector, function(win) {
                                                    pmax(maxChrLamWComp[,win],
                                                    minchrLamWComp[,win],
                                                    chrLamb[win])
                                                    }))
    return(lamlocRleList)
}


#' computeCoverageMovingWindowOnChr
#' @description computes the coverage on a chromosomewith a
#' set of moving windows of dimensions minWinWidth:maxWinWidth
#'
#' @param chrBedGRanges a GRanges to compute the coverage
#' @param minWinWidth the minimum width of the window to use for the coverage
#' @param maxWinWidth the maximum width of the window to use for the coverage
#' @param binWidth the dimension of the bin in base number
#'
#' @return RleList where each element is a window within the Rle of its coverage
#'
#' @importFrom GenomicRanges tileGenome coverage
#' @importFrom IRanges RleList
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom GenomeInfoDb seqinfo seqnames
#' @importFrom S4Vectors runValue
#'
#' @keywords internal
computeCoverageMovingWindowOnChr <- function(chrBedGRanges, minWinWidth=50,
                                                maxWinWidth=1000, binWidth=50,
                                                verbose=TRUE)
{
    ## converting base notation to bin notation
    minWinWidth <- minWinWidth/binWidth
    maxWinWidth <- maxWinWidth/binWidth

    lengths <- GenomeInfoDb::seqinfo(chrBedGRanges)
    ## dividing chromosome in bins of binWidth dimention each
    binnedChromosome <- GenomicRanges::tileGenome(seqlengths=lengths,
                                                    tilewidth=binWidth,
                                                    cut.last.tile.in.chrom=TRUE)

    # olap <- GenomicAlignments::summarizeOverlaps(features=binnedChromosome,
    #                                               reads=chrBedGRanges,
    #                                               ignore.strand=TRUE,
    #                                               inter.feature=FALSE,
    #                                               mode="IntersectionNotEmpty")
    # #
    # # test the coverages with the old coverages using a small amount of reads
    # chrCovRle1 <- as(SummarizedExperiment::assays(olap)$counts, "Rle")

    ## computing coverage per single base on bed
    chrCoverage <- GenomicRanges::coverage(x=chrBedGRanges)
    ## computing coverage per each bin on chromosome
    if(verbose) message("Computing coverage on Chromosome ",
                S4Vectors::runValue(GenomeInfoDb::seqnames(chrBedGRanges)),
                " binned by ", binWidth, " bin dimension")

    chrCovRle <- binnedCovOnly(bins=binnedChromosome,
                                numvar=chrCoverage,
                                mcolname="bin_cov")

    wins <- minWinWidth:maxWinWidth
    runWinRleList <- IRanges::RleList(
                        lapply(wins, function(win) {
                            if(verbose) {
                                message("Running window ", (win*binWidth))
                            }
                            floor(evenRunMean(x=chrCovRle, k=win,
                                                endrule="constant"))
                        })
                    )
    return(runWinRleList)
}


#' binnedCoverage
#' @description  this function computes the coverage over a binned
#' chromosome, starting from a per base computed coverage.
#'
#' @param bins a GRanges object representing a chromosome binned.
#' @param numvar an RleList representing the per base coverage over the chr.
#' @param mcolname the name of column where the sum have to be stored.
#' @param covMethod a method to apply for the computing of the coverate
#' it can be one of "max", "mean", "sum", "min". ("max" is default)
#' @param roundingMethod a method to apply to round the computations
#' it can be one of "none", "floor", "ceiling", "round".
#' It's useful only when using covMethod="mean". ("none" is default)
#'
#' @return the bins GRanges with the mcolname attached
#' @export
#'
#' @importFrom GenomeInfoDb seqlevels seqnames
#' @importFrom S4Vectors split mcols
#' @importFrom IRanges ranges Views viewSums viewMaxs viewMeans viewMins
#' @examples
#' ## dividing one chromosome in bins of 50 bp each
#' seqinfo <- GenomeInfoDb::Seqinfo(genome="mm9")
#' bins <- GenomicRanges::tileGenome(
#'             seqlengths=GenomeInfoDb::seqlengths(seqinfo)[1],
#'             tilewidth=50,
#'             cut.last.tile.in.chrom=TRUE)
#' gr <- GenomicRanges::GRanges(seqnames = S4Vectors::Rle("chr1", 100),
#'             ranges=IRanges::IRanges(start = seq(from=10, to=1000, by=10),
#'             end=seq(from=20, to=1010, by = 10)))
## computing coverage per single base on granges
#' cov <- GenomicRanges::coverage(x=gr)
#' (binnedMaxCovGR <- binnedCoverage(bins, cov, "binned_cov"))
#' (binnedMeanCovGR <- binnedCoverage(bins, cov, "binned_cov",
#'                                 covMethod="mean", roundingMethod="floor"))
#' (binnedSumCovGR <- binnedCoverage(bins, cov, "binned_cov", covMethod="sum"))
binnedCoverage <- function(bins, numvar, mcolname,
                        covMethod=c("max", "mean", "sum", "min"),
                        roundingMethod=c("none", "floor", "ceiling", "round"))
{
    covMethod <- match.arg(covMethod)
    roundingMethod <- match.arg(roundingMethod)
    stopifnot(is(bins, "GRanges"))
    stopifnot(is(numvar, "RleList"))
    stopifnot(identical(GenomeInfoDb::seqlevels(bins), names(numvar)))

    binsChr <- S4Vectors::split(IRanges::ranges(bins),
                                GenomeInfoDb::seqnames(bins))
    binCovsR <- lapply(names(numvar), function(seqname)
    {
        views <- IRanges::Views(numvar[[seqname]], binsChr[[seqname]])
        binCovs <- switch(covMethod,
                            max=IRanges::viewMaxs(views),
                            mean=IRanges::viewMeans(views),
                            sum=IRanges::viewSums(views),
                            min=IRanges::viewMins(views)
        )
        binCovsR <- switch(roundingMethod,
                            none=binCovs,
                            floor=floor(binCovs),
                            ceiling=ceiling(binCovs),
                            round=round(binCovs)
        )
        return(binCovsR)
    })

    new_mcol <- unsplit(binCovsR, as.factor(GenomeInfoDb::seqnames(bins)))
    S4Vectors::mcols(bins)[[mcolname]] <- new_mcol

    return(bins)
}

#' binnedCovOnly
#' @description it's useful just to coerce the bin coverage to an Rle object
#'
#' @param bins a GRanges object representing a chromosome binned
#' @param numvar an RleList representing the per base coverage over the chr
#' @param mcolname the name of column where the sum have to be stored
#'
#' @importFrom S4Vectors mcols
#'
#' @return an Rle within the per bin computed coverage
#' @keywords internal
binnedCovOnly <- function(bins, numvar, mcolname)
{
    binsGRanges <- binnedCoverage(bins=bins, numvar=numvar, mcolname=mcolname,
                                    covMethod="max", roundingMethod="none")
    ## coercing just the binned coverage to Rle
    chrCovRle <- as(S4Vectors::mcols(binsGRanges)[[mcolname]], "Rle")
    return(chrCovRle)
}

#' evenRunSum
#' @description this function computes a running sum over x with a window
#' width k (modified from S4Vectors package to work on even k, in such a case
#' it adds a length at the end of the output Rle).
#'
#' @param x an Rle object, typically a coverage object.
#' @param k window dimension for the running sum over x.
#' @param endrule refer to S4Vectors::runSum.
#' @param na.rm refer to S4Vectors::runSum.
#'
#' @importFrom S4Vectors .Call2 runLength nrun
#'
#' @return an Rle within the running sum over x with a win of length k.
#' @keywords internal
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
        if( (k %% 2L) == 0 ) j=j+1 ##
        S4Vectors::runLength(ans)[S4Vectors::nrun(ans)] <-
                        S4Vectors::runLength(ans)[S4Vectors::nrun(ans)]+(j - 1L)
    }
    return(ans)
}


#' evenRunMean
#' @description this function computes a running mean over x with a window
#' width k (modified from S4Vectors package to work on even k, see evenRunSum).
#'
#' @param x an Rle object, typically a coverage object.
#' @param k window dimension for the running sum over x.
#' @param endrule refer to S4Vectors::runMean.
#' @param na.rm refer to S4Vectors::runMean.
#'
#' @importFrom S4Vectors Rle
#'
#' @return an Rle within the running mean over x with a win of length k.
#' @keywords internal
evenRunMean <- function(x, k, endrule = c("drop", "constant"), na.rm = FALSE)
{
    sums <- evenRunSum(x, k, endrule, na.rm)
    if (na.rm)
    {
        d <- S4Vectors::Rle(rep(1L, length(x)))
        d[is.na(x)] <- 0L
        means <- (sums / evenRunSum(d, k, endrule, na.rm))
    } else {
        means <- (sums / k)
    }
    return(means)
}



#' binToChrCoordMatRowNames
#' @description computes the starting range of the bins for the
#' binMatrix, taking in input the length of the chromosome of the matrix.
#'
#' @param binMatrix a matrix where each row represents a bin.
#' @param chrLength the length of the chromosome of the binMatrix.
#' @param binWidth the width of the bin.
#'
#' @return the binMatrix with start range as rownames.
#' @keywords internal
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
        stop("something went wrong! matrix row dimension",
                "different than expected")
    }
    rownames(binMatrix) <- startBinRanges
    return(binMatrix)
}

