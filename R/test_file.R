# binnedAverage <- function (bins, numvar, varname, na.rm = FALSE)
# {
#     if (!is(bins, "GRanges"))
#         stop("'x' must be a GRanges object")
#     if (!is(numvar, "RleList"))
#         stop("'numvar' must be an RleList object")
#     if (!identical(seqlevels(bins), names(numvar)))
#         stop("'seqlevels(bin)' and 'names(numvar)' must be identical")
#     viewMeans2 <- function(v, na.rm = FALSE) {
#         if (!isTRUEorFALSE(na.rm))
#             stop("'na.rm' must be TRUE or FALSE")
#         means <- viewMeans(v, na.rm = na.rm)
#         w0 <- width(v)
#         v1 <- trim(v)
#         w1 <- width(v1)
#         if (na.rm) {
#             na_count <- sum(is.na(v1))
#             w0 <- w0 - na_count
#             w1 <- w1 - na_count
#         }
#         means <- means * w1/w0
#         means[w0 != 0L & w1 == 0L] <- 0
#         means
#     }
#     bins_per_chrom <- split(ranges(bins), seqnames(bins))
#     means_list <- lapply(names(numvar), function(seqname) {
#         v <- Views(numvar[[seqname]], bins_per_chrom[[seqname]])
#         viewMeans2(v, na.rm = na.rm)
#     })
#     new_mcol <- unsplit(means_list, as.factor(seqnames(bins)))
#     mcols(bins)[[varname]] <- new_mcol
#     bins
# }

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
    bins
}

oddRunSum<-function(x, k, endrule = c("drop", "constant"), na.rm = FALSE)
{
    # endrule <- match.arg(endrule)
    n <- length(x)
    # k <- normargRunK(k = k, n = n, endrule = endrule)
    ans <- S4Vectors::.Call2("Rle_runsum", x, as.integer(k), as.logical(na.rm),
                  PACKAGE="S4Vectors")
    ans
    if (endrule == "constant") {
        j <- (k + 1L) %/% 2L

        runLength(ans)[1L] <- runLength(ans)[1L] + (j - 1L)
        if( (k %% 2L) == 0 ) j=j+1
        runLength(ans)[S4Vectors::nrun(ans)] <-
            runLength(ans)[S4Vectors::nrun(ans)] + (j - 1L)
    }
    return(ans)
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
        bgr
    })
    names(bedGRList) <- interestedChrs
    bedGRList <- GRangesList(bedGRList)

    return(bedGRList)
}


# chr19Coverage <- RleList(Coverage$chr19)
# names(chr19Coverage) = 'chr19'

# binnAvgChr19 <- GenomicRanges::binnedAverage(bins=mm9Chr19Tile50, numvar=RleList(chr19Coverage), varname="mean")
# bedGRanges <- constructBedRanges(filename=filename, filetype="bed", genomeName="mm9")
####################################################


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
    message("Computing coverage on Chromosome ", chrBedGRanges@seqinfo@seqnames,
            " binned by ", binWidth, " bin dimension")
    binnedCovChr <- binnedSum(bins=binnedChromosome,
                                numvar=chrCoverage,
                                mcolname="bin_cov")

    ## computing bin in base ranges to add as rownames
    endBinRanges <- seq(from=binWidth-1,
                        to=chrBedGRanges@seqinfo@seqlengths[1],
                        by=binWidth)

    if(chrBedGRanges@seqinfo@seqlengths != endBinRanges[length(endBinRanges)])
    {
        ## it can only be lesser than the length of the chromosome
        ## adding the last missing bin
        endBinRanges <- c(endBinRanges, (chrBedGRanges@seqinfo@seqlengths[1]-1))
    }
    ## computing the start of regions
    startBinRanges <- endBinRanges-(binWidth-1)



    ## coercing coverage to Rle
    chrCovRle <- as(binnedCovChr$bin_cov,"Rle")

    runWinRleList <- RleList()
    for(win in minWinWidth:maxWinWidth) {
        message("Running window ", win, " of ", maxWinWidth)
        runWinRleList[[win]] <- oddRunSum(chrCovRle, k=win, endrule="constant")

        # maxRuns <- max(maxRuns, sum(runLength(runwinRleList[[win]])) ) ##as control
        # message("maxRuns: ", maxRuns)
        ## maxRuns is equal to the number of bins ans is the number of
        ## chromosome bases divided by the binSize value as expected
        if(win == minWinWidth)
        {
            rleBinCov <- DelayedArray::RleArray(runWinRleList[[win]],
                                       dim=c(length(runWinRleList[[win]]), 1))
            rownames(rleBinCov) <- startBinRanges

        }
        else
        {
            rleBinCov <- cbind(rleBinCov,
                               DelayedArray::RleArray(runWinRleList[[win]],
                               dim=c(length(runWinRleList[[win]]), 1)))
        }

    }


    message("Running window 5000")
    win5k <- 5000
    runwinRle5K <- oddRunSum(chrCovRle, k=win5k, endrule="constant")
    rleBin5K <- DelayedArray::RleArray(runwinRle5K,
                                       dim=c(length(runwinRle5K), 1))
    message("Running window 10000")
    win10k <- 10000
    runwinRle10K <- oddRunSum(chrCovRle, k=win10k, endrule="constant")
    rleBin10K <- DelayedArray::RleArray(runwinRle10K,
                                       dim=c(length(runwinRle10K), 1))

    return(list(
        "rleBinCov"=rleBinCov,
        "rleBin5K"=rleBin5K,
        "rleBin10K"=rleBin10K)
    )
}


computeCoverageMovingWindow <- function(bedGRanges, minWinWidth=1,
                                        maxWinWidth=20, binWidth=50)
{
    bedGrangesChrsList <- cutGRangesPerChromosome(bedGRanges)

     chrRleMatricesList <- lapply(X=bedGrangesChrsList,
                                  computeCoverageMovingWindowOnChr,
                                  minWinWidth=minWinWidth,
                                  maxWinWidth=maxWinWidth,
                                  binWidth=binWidth) ##TO PARALLELIZE
    return(chrRleMatricesList)
}

# minWin=1
# maxWin=20
#
# genomeSeqinfo = bed@seqinfo
# str(genomeSeqinfo)
#
# genomeSeqinfoChr19mm9 = genomeSeqinfo["chr19",]
#
# binSize=50
#
# mm9Chr19Tile50 <- GenomicRanges::tileGenome(seqlengths=genomeSeqinfoChr19mm9, tilewidth=binSize, cut.last.tile.in.chrom=TRUE)
#
# bed19 <- keepSeqlevels(bed, "chr19")
#
# chr19Len=seqinfo(bed19)@seqlengths
# maxRuns <- 0
# runwinRleList <- RleList()
# Coverage <- GenomicRanges::coverage(x=bed19)#, width=(chr19Len+win-1))
# message("win: ", win)
# print(Coverage)
# binnedSumChr19 <- binnedSum(bins=mm9Chr19Tile50, numvar=Coverage, mcolname="coverage")
# coverageChr19Rle <- as(binnedSumChr19$coverage,"Rle")
#
#
# for(win in minWin:maxWin) {
#
#     runwinRleList[[win]] <- oddRunSum(coverageChr19Rle, k=win, endrule="constant")
#     print(win)
#
#     maxRuns <- max(maxRuns, sum(runLength(runwinRleList[[win]])) ) ##as control
#     message("maxRuns: ", maxRuns)
#     # if(win == minWin) {
#     #     rleBinCov <- DelayedArray::DelayedArray(runwinRleList[[win]], dim=c(length(runwinRleList[[win]] ), 1))
#     #
#     # } else {
#     #     rleBinCov <- cbind(rleBinCov, DelayedArray::RleArray(runwinRleList[[win]], dim=c(length(runwinRleList[[win]] ), 1)))
#     # }
#     ## maxRuns is equal to the number of bins ans is the number of chromosome
#     ## bases divided by the binSize value as expected
# }
#
# endBinRanges <- seq(from=binSize-1, to=genomeSeqinfoChr19mm9@seqlengths[1], by=binSize)
# if (genomeSeqinfoChr19mm9@seqlengths[1] != endBinRanges[length(endBinRanges)] ) {
#     ## it can only be lesser than the length of the chromosome
#     endBinRanges <- c(endBinRanges, (genomeSeqinfoChr19mm9@seqlengths[1] -1))
# }
# startBinRanges <- endBinRanges-49
#
# # rownames(rleBinCov) <- startBinRanges
#
# for(win in minWin:maxWin) {
#     if(win == minWin) {
#         rleBinCov <- DelayedArray::RleArray(runwinRleList[[win]], dim=c(length(runwinRleList[[win]] ), 1))
#         rownames(rleBinCov) <- startBinRanges
#     } else {
#         rleWinArr <- DelayedArray::RleArray(runwinRleList[[win]], dim=c(length(runwinRleList[[win]] ), 1))
#         rleBinCov <- cbind(rleBinCov, rleWinArr)
#     }
#     print(win)
# }
#
# rownames(rleBinCov)
#
# win5k <- 5000
# runwinRle5K <- oddRunSum(coverageChr19Rle, k=win5k, endrule="constant")
#
# win10k <- 10000
# runwinRle10K <- oddRunSum(coverageChr19Rle, k=win10k, endrule="constant")



# normargRunK <- function(k, n, endrule)
# {
#     if (!is.numeric(k))
#         stop("'k' must be a numeric vector")
#     if (k < 0)
#         stop("'k' must be positive")
#     if ((endrule != "drop") && (k %% 2 == 0)) {
#         k <- 1L + 2L * (k %/% 2L)
#         warning(paste("'k' must be odd when 'endrule != \"drop\"'!",
#                       "Changing 'k' to ", k))
#     }
#     if (k > n) {
#         k <- 1L + 2L * ((n - 1L) %/% 2L)
#         warning("'k' is bigger than 'n'! Changing 'k' to ", k)
#     }
#     as.integer(k)
# }

#
# str(runwinRleList[[1]])
#
# RleArray(rle=runwinRleList[[1]], dim=c(1226830, 1))
#
#  binnedSumChr19
#
# runsum(x=binnedSumChr19, k=1)
#
#
# seqlevels(mm9Chr19Tile50)
# names(chr19Coverage$chr19)
#
# # win = 50
# minStep = 1
# maxStep = 20

# for (step in minStep:maxStep) {
#     #compute seqinfo on single chr
#     tileOnChr <- GenomicRanges::tileGenome(seqlengths=seqInfoChr, tilewidth= step * win, cut.last.tile.in.chrom=TRUE)
# }
#
# averagePerBin <- function(x, binsize, mcolnames=NULL)
# {
#     if (!is(x, "GenomicRanges"))
#         stop("'x' must be a GenomicRanges object")
#     if (any(is.na(seqlengths(x))))
#         stop("'seqlengths(x)' contains NAs")
#     bins <- IRangesList(lapply(seqlengths(x),
#                                function(seqlen)
#                                    IRanges(breakInChunks(seqlen, binsize))))
#     ans <- as(bins, "GRanges")
#     seqinfo(ans) <- seqinfo(x)
#     if (is.null(mcolnames))
#         return(ans)
#     averageMCol <- function(colname)
#     {
#         cvg <- coverage(x, weight=colname)
#         views_list <- RleViewsList(
#             lapply(names(cvg),
#                    function(seqname)
#                        Views(cvg[[seqname]], bins[[seqname]])))
#         unlist(viewMeans(views_list), use.names=FALSE)
#     }
#     mcols(ans) <- DataFrame(lapply(mcols(x)[mcolnames], averageMCol))
#     ans
# }
# averagePerBin(x=bed, binsize=50, mcolnames="boh")
#
# bs <- binnedSum(bins=bins, numvar=numvar, mcolname="summed")
#
# bs[ which(bs$summed != 0), ]



# rm(mm9Chr19Tile1)

# mm9 <- Seqinfo(genome="mm9") <<<<< ----- seqinfo can return the lengths of several genomes
# mm9[mm9@seqnames[which( mm9@seqnames == "chr19")],]
# mm9["chr19",]

# class()
