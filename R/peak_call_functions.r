#' This function calls peaks from bed or bam inputs using a variable window scan
#' with a poisson model using the surrounding 10kb as background.
#' files, filetype="bam", chr=1:19, fragmentLength=200,
#' @param files Character vector containing paths of files to be analyzed.
#' @param filetype Character, either "bam" or "bed" indicating format of input
#' file.
#' @param chr An integer vector specifying which chromosomes to limit analysis
#' to.
#' @param fragmentLength Integer indicating the average DNA fragment length in bp.
#' @param readLength Integer indicating the read lengths of experiment in bp.
#' @param minWin Integer indicating the minimum window size in units of 50 bp,
#' i.e., min_win=2 resulting in a 100 bp window.
#' @param binSize Integer size in base pairs of the minimum window for scanning,
#' 50 is the default.
#' @param maxWin Integer indicating the maximum allowed window size in units of
#'    50 bp.
#' @param blocksize Integer indicating how much of the chromosome
#' will be analyzed at a time in order to avoid memory issues.
#' @param zthresh Cuttoff value for z-scores. Only windows with greater z-scores
#'    will be retained.
#' @param minCount A small constant (usually no larger than one) to be added to
#' the counts prior to the log transformation to avoid problems with log(0).
#' @param outputName A string, Name of the folder to save the Peaks (optional),
#' if the directory doesn't exist, it will be created. default is "Peaks"
#' @param save Boolean, if TRUE files will be saved in a "./Peaks/chr*"
#' directory created (if not already present) in the current working directory.
#' @return a matrix with peaks as rows and 4 columns describing the genomic
#' coordinates (chr, start, end) as well as the associated z-score.
#' @export
#' @importFrom utils read.table write.table
#' @examples
#' bam <- system.file("extdata", "test.bam", package="DEScan")
#' peaks <- findPeaks(bam, chr=1, filetype="bam")
#' head(peaks)

#### ONLY FOR SELF TEST
#bed.path <- 'testData/Bed/chr19'
#bed.files <- list.files(bed.path, full.names = TRUE)
#bam.path <- "testData/bams"
#bam.files <- list.files(bam.path, full.names = TRUE)
#
#filetype="bam";
# chr=19; fragmentLength=200;
# readLength=100;
# binSize=50;
# minWin=1; maxWin=20;
# blocksize=10000; zthresh=5; minCount=0.1;
# outputName="Peaks"; save=FALSE; verbose=TRUE
# file <- bam.files[1]
# genomeName=NULL
findPeaks <- function(files, filetype="bam", chr=1:19,
                      genomeName=NULL, fragmentLength=200,
                      readLength=100,  binSize=50, minWin=1, maxWin=20,
                      blocksize=10000, zthresh=5, minCount=0.1,
                      outputName="Peaks", save=TRUE, verbose=FALSE) {

        if(length(files) == 0) {
            stop("You have to provide one or more input files!\nExiting.")
        }

        for (file in files) {
            # if ( !is.character(chr) ) {
            #     chr <- paste0("chr", chr)
            # }
            # if (filetype == "bam") {
            #     bed <- readBamAsBed(file=file)#, chr=chr)
            # }
            # if (filetype == "bed") {
            #     bed <- readBedFile(filename=file)#, chr=chr)
            # }
            bedGRanges <- constructBedRanges(filename=file, filetype=filetype,
                                             genomeName=genomeName)


            chrRleListAllW <- computeCoverageMovingWindow(bedGRanges=bedGRanges,
                                                          minWinWidth=minWin,
                                                          maxWinWidth=maxWin,
                                                          binWidth=binSize)

            chrRleListW5K <- computeCoverageMovingWindow(bedGRanges=bedGRanges,
                                                          minWinWidth=5000,
                                                          maxWinWidth=5000,
                                                          binWidth=binSize)

            chrRleListW10K <- computeCoverageMovingWindow(bedGRanges=bedGRanges,
                                                         minWinWidth=10000,
                                                         maxWinWidth=10000,
                                                         binWidth=binSize)

            # Get lambdalocal.

            # # compute coverage using a 5kb window
            # headstart <- ceiling(5000 / gsize) * gsize
            # grid05k <- c(seq(grid0[1] - headstart,
            #                  grid0[1], by=gsize), grid0)
            # offset5k <- length(grid05k) - length(grid0)
            # c5k <- window_coverage(bed=bed, Fr=fr, Rr=rr, chr=chr, grid=grid05k,
            #                        fragmentLength=fragmentLength,
            #                        readLength=readLength,
            #                        maxWin=5000, minWin=5000, verbose=FALSE)
            #
            # # compute coverage using a 10kb window
            # headstart <- ceiling(10000 / gsize) * gsize
            # grid010k <- c(seq(grid0[1] - headstart, grid0[1], by=gsize),
            #               grid0)
            # offset10k <- length(grid010k) - length(grid0)
            # c10k <- window_coverage(bed=bed, Fr=fr, Rr=rr, grid=grid010k,
            #                         fragmentLength=fragmentLength, chr=chr,
            #                         readLength=readLength, maxWin=10000,
            #                         minWin=10000, verbose=FALSE)
            #
            # # determine lambda for 5k, 10k windows and baseline
            # lamloc <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
            # lam5k <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
            # lam10k <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
            # lambl <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
            # for (win in minWin:maxWin) {
            #     lam5k[, win - minWin + 1] <- c5k[offset5k + c(1:length(grid0)) - floor(win / 2)] * win / 5000
            #     lam10k[, win - minWin + 1] <- c10k[offset10k + c(1:length(grid0)) - floor(win / 2)] * win / 10000
            #     lambl[, win - minWin + 1] <- tot_rds * win / tot_base
            # }

            ## check!
            fr <- bedGRanges@ranges@start[which( as.vector(bedGRanges@strand) == "+")]
            rr <- bedGRanges@ranges@start[which( as.vector(bedGRanges@strand) == "-")]

            tot_rds <- length(fr) + length(rr)
            tot_base <- max(rr[length(rr)], fr[length(fr)]) - min(fr[1], rr[1])

            ## here, rewrite as a function for each chromosome
            lam5k <- lapply(chrRleListW5K, function(w5k) {
                    w5k %*% t(minWin:maxWin) / 5000
            })

            lam10k <- lapply(chrRleListW5K, function(w5k) {
                w5k %*% t(minWin:maxWin) / 10000
            })

            lambl_num <- Rle(tot_rds %*% t(minWin:maxWin) / tot_base)

            ## here only chr19
            lambl <- RleArray(rep(lambl_num, each=nrow(lam5k$chr19)),
                              dim=c(nrow(lam5k$chr19), ncol(lam5k$chr19)))

            lamloc <- pmax(lam5k$chr19, lam10k$chr19, lambl)
            ## calculate z score for each bin x window combination
            z <- sqrt(2) * sign(cmat - lamloc) *
                sqrt(cmat * log(pmax(cmat, minCount) / lamloc) -
                         (cmat - lamloc))

            ## find high z scores keeping one with no intersecting other
            ## bin/windows
            new_s <- get_disjoint_max_win(z0=z[1:blocksize_i, ],
                                          sigwin=fragmentLength / gsize, nmax=Inf,
                                          zthresh=zthresh, two_sided=FALSE,
                                          verbose=TRUE)
            ## convert new_s bins and width into genomic coordinates and
            ## append to s
            if (nrow(new_s) >= 1) {
                new_s[, 1] <- new_s[, 1] + block[1] - 1
                new_s <- cbind(grid[new_s[, 1, drop=FALSE]],
                               grid[new_s[, 1, drop=FALSE] +
                                        new_s[, 2, drop=FALSE] - 1],
                               new_s[, 3, drop=FALSE])
                s <- rbind(s, new_s)
            }

            elapsed <- proc.time() - ptm
            if (verbose) {
                cat("\tDone. That took ", format(elapsed[3] / 60, digits=1),
                    " minutes.\n")
            }
            if (finalblock) break

            # shift block by blocksize and repeat
            block <- c(block[1] + blocksize_i, block[2] + blocksize_i)
        }

    peaks <- cbind(rep(chr, dim(s)[1]), s)
    if (save == TRUE) {
        if (dir.exists(outputName) == FALSE) {
            dir.create(outputName)
        }
        if (dir.exists(paste0(outputName,"/", chr)) == FALSE) {
            dir.create(paste0(outputName,"/", chr))
        }
        fname <- strsplit(basename(file), split=".", fixed=TRUE)[[1]][1]
        fileprefix <- paste0(outputName,"/", chr, "/Peaks_", fname)


        save(peaks, fragmentLength, readLength, zthresh, minWin, maxWin,
             file=paste0(fileprefix, ".RData"))
    }

        return(peaks)
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
    message("Computing coverage on Chromosome ", chrBedGRanges@seqinfo@seqnames,
            " binned by ", binWidth, " bin dimension")
    # binnedCovChr <- binnedSum(bins=binnedChromosome,
    #                           numvar=chrCoverage,
    #                           mcolname="bin_cov")
    chrCovRle <- binnedSumOnly(bins=binnedChromosome,
                                    numvar=chrCoverage,
                                    mcolname="bin_cov")

    # ## coercing coverage to Rle
    # chrCovRle <- as(binnedCovChr$bin_cov,"Rle")


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

    idx <- minWinWidth:maxWinWidth
    runWinRleList <- RleList(lapply(idx, function(win) {
        oddRunSum(chrCovRle, k=win, endrule="constant")
    }))

    rleBinCov <- RleArray(unlist(runWinRleList),
                       dim=c(length(runWinRleList[[1]]), length(runWinRleList)))

    # # runWinRleList <- RleList()
    # for(win in minWinWidth:maxWinWidth) {
    #     message("Running window ", win, " of ", maxWinWidth)
    #     runWinRle <- oddRunSum(chrCovRle, k=win, endrule="constant")
    #
    #     # maxRuns <- max(maxRuns, sum(runLength(runwinRleList[[win]])) ) ##as control
    #     # message("maxRuns: ", maxRuns)
    #     ## maxRuns is equal to the number of bins ans is the number of
    #     ## chromosome bases divided by the binSize value as expected
    #     if(win == minWinWidth)
    #     {
    #         rleBinCov <- DelayedArray::RleArray(runWinRle,
    #                                             dim=c(length(runWinRle), 1))
    #         rownames(rleBinCov) <- startBinRanges
    #
    #     }
    #     else
    #     {
    #         rleBinCov <- cbind(rleBinCov,
    #                            DelayedArray::RleArray(runWinRle,
    #                                             dim=c(length(runWinRle), 1)))
    #     }
    #
    # }

    return(rleBinCov)
}


#' Title
#'
#' @param bedGRanges
#' @param minWinWidth
#' @param maxWinWidth
#' @param binWidth
#'
#' @return
#' @export
#'
#' @examples
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


#' compute window coverage.
#'
#' @param Fr forward reads.
#' @param Rr reverse reads.
#' @param grid
#' @param fraglen
#' @param rlen
#' @param min_win
#' @param max_win
#' @param verbose
#' @return cmat.
#' @keywords internal
#'
#'
# window_coverage <- function(Fr, Rr, grid, fragmentLength=200, readLength=100,
# minWin=1, maxWin=100, verbose=TRUE) {
window_coverage <- function(bed, Fr, Rr, grid, fragmentLength=200, readLength=100,
                            minWin=1, maxWin=100, chr, verbose=TRUE) {

    # cmat <- matrix(nrow=length(grid), ncol=maxWin - minWin + 1, data=0)
    # # dim(cmat)
    # for(i in minWin:maxWin) {
    #     bedCoverageWind <- GenomicRanges::coverage(x=bed, shift=(i-1), width=length(grid))# + 2*(i-1)))
    #
    #     bedCoverageWindChr <- bedCoverageWind[which(names(bedCoverageWind) %in% chr)]
    #
    #     print(bedCoverageWindChr$chr19)
    #
    #     bedCoverageSum <- S4Vectors::runsum(x=bedCoverageWindChr$chr19, k=i, endrule="drop")
    #
    #
    #     # cmat[,i] <- bedCoverageSum
    # }


    grid_st <- min(grid)
    ## Assume even grid. ! only if min_rds=0
    bin_size <- grid[2] - grid[1]
    ## the number of bins covered by each fragment
    fragbins <- fragmentLength / bin_size
    ngrid <- length(grid)

    ## initialize fragment count matrices for forward and reverse reads
    ffragctmat <- matrix(nrow=length(grid), ncol=maxWin - minWin + 1, data=0)
    rfragctmat <- matrix(nrow=length(grid), ncol=maxWin - minWin + 1, data=0)

    ## index forward reads within the current grid
    forward_index <- which(Fr > grid[1] - fragmentLength & Fr < grid[length(grid)])

    ## for each forward read within the given grid convert the base position
    ## to bin and add +1 count to the fragment count matrix (bins as rows,
    ## window size as columns) while moving window right along bins and
    ## extending left
    if (length(forward_index) > 0) {
        for (i in 1:length(forward_index)) {
            if (verbose && i %% 1000 == 0) {
                cat("Forward reads: ", i, "out of", length(forward_index), "\n") ## use message instead
            }
            fragst <- Fr[forward_index[i]]    # start of fragment
            # which bin contains start of fragment
            binst <- floor((fragst - grid_st) / bin_size) + 1
            for (j in minWin:maxWin) {

                st <- max(binst - j + 1, 1) # lower limit adj by current window size
                ed <- min(binst + fragbins - 1, ngrid) # upper limit
                ffragctmat[st:ed, j - minWin + 1] <- ffragctmat[st:ed, j -
                                                                minWin + 1] + 1
            }
        }
    }

    reverse_index <- which(Rr > grid[1] - readLength & (Rr + readLength - fragmentLength) <
                               grid[length(grid)])

    if (length(reverse_index) > 0) {
        for (i in 1:length(reverse_index)) {
            if (verbose && i %% 1000 == 0) {
                cat("Reverse reads: ", i, "out of", length(reverse_index), "\n") ## use message instead
            }
            ## start of fragment
            fragst <- Rr[reverse_index[i]] + readLength - fragmentLength
            ## bin containing start of fragment
            binst <- floor((fragst - grid_st) / bin_size) + 1
            for (j in minWin:maxWin) {
                st <- max(binst - j + 1, 1)
                ed <- min(binst + fragbins - 1, ngrid)
                rfragctmat[st:ed, j - minWin + 1] <- rfragctmat[st:ed, j -
                                                                minWin + 1] + 1
            }
        }
    }
    cmat <- ffragctmat + rfragctmat
    return(cmat)
}



#' make grid for enrichment scan.
#'
#' @param path Path to alignement files.
#' @param fraglen Integer indicating the average DNA fragment length in bp
#' @param rlen Integer indicating the read length in bp
#' @param min_bin Integer indicating the minimum window size to be used in
#'     constructing grid
#' @param min_rds This argument is currently not used.
#' @param gridStart Genomic Position to start the grid
#' @param grid_ed Genomic Position to end the grid
#' @return grid.
#' @keywords internal
makeGrid <- function(bed, fragmentLength, readLength,
                     binSize=10, minRds=0,
                     gridStart=NA, gridEnd=NA,
                     verbose=FALSE) {

    fr <- bed@ranges@start[which( as.vector(bed@strand) == "+")]
    rr <- bed@ranges@start[which( as.vector(bed@strand) == "-")]

    if (is.na(gridStart)) {
        gridStart <- min(fr, (rr + readLength - fragmentLength))
    }
    if (is.na(gridEnd)) {
        gridEnd <- max((fr + fragmentLength), (rr + readLength))
    }

    # Initial grid is even spaced:
    grid0 <- seq(gridStart, gridEnd, binSize)
    ## all ok till here!!! dr

    # Needs work
    if (minRds > 0) { # this portion of the code does not work!
        # Get the number of fragments that cover each bin.
        # fragct[i] is the number of fragments with positive overlap with
        # grid0[i]:grid0[i+1]
        ffragct <- matrix(nrow=length(grid0), ncol=length(fr), data=0)
        # Error: cannot allocate vector of size 2056.4 Gb #####################
        rfragct <- matrix(nrow=length(grid0), ncol=length(fr), data=0)
        for (k in seq_along(fr)) {
            for (i in 1:length(fr[[k]])) {
                if (i %% 10000 == 0) {
                    cat("Sample ", k, "Forward reads : ", i, "\n")
                }
                fragst <- fr[[k]][i]
                fraged <- fr[[k]][i] + fragmentLength - 1
                binst <- floor((fragst - gridStart) / binSize) + 1
                bined <- floor((fraged - gridStart) / binSize) + 1
                ffragct[binst:bined, k] <- ffragct[binst:bined, k] + 1
            }
        }
        for (k in 1:length(rr)) {
            for (i in 1:length(rr[[k]])) {
                if (i %% 10000 == 0) {
                    cat("Sample ", k, "Reverse reads : ", i, "\n")
                }
                fragst <- rr[[k]][i] + readLength - fragmentLength
                fraged <- rr[[k]][i] + readLength - 1
                binst <- floor((fragst - gridStart) / binSize) + 1
                bined <- floor((fraged - gridStart) / binSize) + 1
                rfragct[binst:bined, k] <- rfragct[binst:bined, k] + 1
            }
        }

        # Starting from the left most bin, join bins whose total fragment
        # coverage across all samples is less than minRds.
        totcov <- rowSums(ffragct) + rowSums(rfragct)
        remove_tick <- rep(FALSE, length(totcov))
        cursum <- 0
        for (i in 1:length(totcov)) {
            if (i %% 100000 == 0 && verbose) {
                cat(i, " out of ", length(totcov), "\n")
            }
            cursum <- cursum + totcov[i]
            if (cursum < minRds) {
                remove_tick[i] <- TRUE
            } else {
                cursum <- 0
            }
        }
        grid <- grid0[which(!remove_tick)]
    } else {
        grid <- grid0
    }
    return(grid)
}

#' find significant z score windows keeping the max value without intersections
#'
#' @param z0 Matrix containing z scores with bins as rows and windows size as
#'     columns
#' @param sigwin Integer indicating how many bins per fragment
#' @param nmax Integer indicating the maximum number of windows to return
#' @param zthresh Integer indicating the minimum z-score considered significant
#' @param two_sided
#' @param verbose
#' @return s.
#' @keywords internal
get_disjoint_max_win <- function(z0, sigwin=20, nmax=Inf,
                                zthresh=-Inf, two_sided=FALSE,
                                verbose=TRUE) {
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
            ## remove surrounding max window around max z from consideration
            st <- max(1, t - sigwin - maxwin + 1)
            ed <- min(t + w + sigwin - 1, nrow(z0))
            z0[st:ed, ] <- -Inf

            if (nrow(s) >= nmax) break
        }
        return(s)
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

        runLength(ans)[1L] <- runLength(ans)[1L] + (j - 1L)
        if( (k %% 2L) == 0 ) j=j+1
        runLength(ans)[S4Vectors::nrun(ans)] <-
            runLength(ans)[S4Vectors::nrun(ans)] + (j - 1L)
    }
    return(ans)
}



