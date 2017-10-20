#' This function calls peaks from bed or bam inputs using a variable window scan
#' with a poisson model using the surrounding 10kb as background.
#'
#' @param files Character vector containing paths of files to be analyzed.
#' @param filetype Character, either "bam" or "bed" indicating format of input
#' file.
#' @param chr An integer vector specifying which chromosomes to limit analysis
#' to.
#' @param fraglen Integer indicating the average DNA fragment length in bp.
#' @param rlen Integer indicating the read lengths of experiment in bp.
#' @param min_win Integer indicating the minimum window size in units of 50 bp,
#' i.e., min_win=2 resulting in a 100 bp window.
#'
#' @param min_bin Integer size in base pairs of the minimum window for scanning,
#' 50 is the default.
#' @param max_win Integer indicating the maximum allowed window size in units of
#'    50 bp.
#' @param blocksize Integer indicating how much of the chromosome
#' will be analyzed at a time in order to avoid memory issues.
#' @param zthresh Cuttoff value for z-scores. Only windows with greater z-scores
#'    will be retained.
#' @param min_count A small constant (usually no larger than one) to be added to
#' the counts prior to the log transformation to avoid problems with log(0).
#' @param output_name A string, Name of the folder to save the Peaks (optional),
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

##### ONLY FOR SELF TEST
# bed.path <- 'testData/Bed/chr19'
# bed.files <- list.files(bed.path, full.names = TRUE)
# bam.path <- "testData/bams"
# bam.files <- list.files(bam.path, full.names = TRUE)
#
# filetype="bam"; chr=1:19; fragmentLength=200;
# readLength=100;  minBin=50; minWin=1; maxWin=20;
# blocksize=10000; zthresh=5; minCount=0.1;
# outputName="Peaks"; save=FALSE; verbose=TRUE

findPeaks <- function(files, filetype="bam", chr=1:19, fragmentLength=200,
                      readLength=100,  minBin=50, minWin=1, maxWin=20,
                      blocksize=10000, zthresh=5, minCount=0.1,
                      outputName="Peaks", save=TRUE, verbose=FALSE) { ## NOT WORKING -> implementing GenomicRanges

        if(length(files) == 0) {
            stop("You have to provide one or more input files!\nExiting.")
        }

        for (file in files) {


            if ( !is.character(chr) ) {
                chr <- paste0("chr", chr)
            }

            if (filetype == "bam") {
                bed <- readBamAsBed(file=file, chr=chr)
            }
            if (filetype == "bed") {
                bed <- readBedFile(filename=file, chr=chr)
            }

            # grid=vector of integers describing the start
            # coordinates for each bin
            grid <- makeGrid(bed=bed, fragmentLength=fragmentLength,
                             readLength=readLength, minBin=minBin)

            ngrid <- length(grid)
            gsize <- grid[2] - grid[1]

            fr <- bed@ranges@start[which(bed@strand=="+")]
            rr <- bed@ranges@start[which(bed@strand=="-")]

            tot_rds <- length(fr) + length(rr)
            tot_base <- max(rr[length(rr)], fr[length(fr)]) - min(fr[1], rr[1])

            # do analysis blockwise due to memory issues
            # fill initial block
            block <- c(1, blocksize + maxWin)
            finalblock <- FALSE
            s <- matrix(ncol=3, nrow=0, data=0)
            #### ok till here dr
            while (TRUE) {
                blocksize_i <- blocksize
                ptm <- proc.time()

                # conditions to stop once as end of grid has been reached
                if (block[1] > ngrid) break
                if (block[2] >= ngrid) {
                    block[2] <- min(block[2], ngrid)
                    blocksize_i <- block[2] - block[1] + 1
                    finalblock <- TRUE
                }
                if (verbose) {
                    cat("Doing block ", block[1], " to ", block[2],
                        " out of ", ngrid, "\n")
                }
                grid0 <- grid[block[1]:block[2]]

                if (verbose) {
                    cat("\tComputing coverage stats...\n")
                }
                ## cmat=coverage matrix with bins as rows and window sizes as
                ## columns
                cmat <- window_coverage(Fr=fr, Rr=rr, grid=grid0,
                                        fragmentLength=fragmentLength,
                                        readLength=readLength, maxWin=maxWin,
                                        minWin=minWin, verbose=FALSE) ############### reinventing wheel?

                # Get lambdalocal.

                # compute coverage using a 5kb window
                headstart <- ceiling(5000 / gsize) * gsize
                grid05k <- c(seq(grid0[1] - headstart,
                                 grid0[1], by=gsize), grid0)
                offset5k <- length(grid05k) - length(grid0)
                c5k <- window_coverage(fr, rr, grid05k,
                                       fragmentLength=fragmentLength,
                                       readLength=readLength,
                                       maxWin=5000, minWin=5000, verbose=FALSE) ############### reinventing wheel?

                # compute coverage using a 10kb window
                headstart <- ceiling(10000 / gsize) * gsize
                grid010k <- c(seq(grid0[1] - headstart, grid0[1], by=gsize),
                              grid0)
                offset10k <- length(grid010k) - length(grid0)
                c10k <- window_coverage(fr, rr, grid010k,
                                        fragmentLength=fragmentLength,
                                        readLength=readLength, maxWin=10000,
                                        minWin=10000, verbose=FALSE) ############### reinventing wheel?

                # determine lambda for 5k, 10k windows and baseline
                lamloc <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
                lam5k <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
                lam10k <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
                lambl <- matrix(nrow=nrow(cmat), ncol=ncol(cmat), data=0)
                for (win in minWin:maxWin) {
                    lam5k[, win - minWin + 1] <- c5k[offset5k +
                                                c(1:length(grid0)) -
                                                floor(win / 2)] * win / 5000

                    lam10k[, win - minWin + 1] <- c10k[offset10k +
                                                 c(1:length(grid0)) -
                                                 floor(win / 2)] * win / 10000
                    lambl[, win - minWin + 1] <- tot_rds * win / tot_base
                }
                lamloc <- pmax(lam5k, lam10k, lambl)
                ## calculate z score for each bin x window combination
                z <- sqrt(2) * sign(cmat - lamloc) *
                    sqrt(cmat * log(pmax(cmat, minCount) / lamloc) -
                    (cmat - lamloc))

                ## find high z scores keeping one with no intersecting other
                ## bin/windows
                new_s <- get_disjoint_max_win(z0=z[1:blocksize_i, ],
                                              sigwin=fragmentLength / gsize, nmax=Inf,
                                              zthresh=zthresh, two_sided=FALSE,
                                              verbose=FALSE)
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
        }

        return(peaks)
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
                     minBin=10, minRds=0,
                     gridStart=NA, gridEnd=NA,
                     verbose=FALSE) {

    fr <- bed@ranges@start[which(bed@strand=="+")]
    rr <- bed@ranges@start[which(bed@strand=="-")]

    if (is.na(gridStart)) {
        gridStart <- min(fr, (rr + readLength - fragmentLength))
    }
    if (is.na(gridEnd)) {
        gridEnd <- max((fr + fragmentLength), (rr + readLength))
    }

    # Initial grid is even spaced:
    grid0 <- seq(gridStart, gridEnd, minBin)
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
                binst <- floor((fragst - gridStart) / minBin) + 1
                bined <- floor((fraged - gridStart) / minBin) + 1
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
                binst <- floor((fragst - gridStart) / minBin) + 1
                bined <- floor((fraged - gridStart) / minBin) + 1
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
window_coverage <- function(Fr, Rr, grid, fragmentLength=200, readLength=100,
                            minWin=1, maxWin=100, verbose=TRUE) {
        grid_st <- min(grid)
        ## Assume even grid. ! only if min_rds=0
        bin_size <- grid[2] - grid[1]
        ## the number of bins covered by each fragment
        fragbins <- fragmentLength / bin_size
        ngrid <- length(grid)

        ## initialize fragment count matrices for forward and reverse reads
        ffragctmat <- matrix(nrow=length(grid), ncol=maxWin - minWin + 1,
                                                 data=0)
        rfragctmat <- matrix(nrow=length(grid), ncol=maxWin - minWin + 1,
                                                 data=0)

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
                    # lower limit adj by current window size
                    st <- max(binst - j + 1, 1)
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

#' read a bam file into a bed like format.
#'
#' @param file Character indicating path to bam file.
#' @param chr Integer indicating which chromosome to read in.
#' @return GenomicRanges object.
#' @keywords internal
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import IRanges
readBamAsBed <- function(file, chr) {

    if ( !file.exists(paste0(file,".bai"))) {
        message("Indexing ", file, " bam file.")
        Rsamtools::indexBam(file)
    }

    bf <- Rsamtools::BamFile(file)
    sl <- GenomeInfoDb::seqlengths(Rsamtools::seqinfo(bf))

    message("Reading: ", file, "\nCromosome: ", chr, "\n")
    which <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, sl[chr]))
    param <- Rsamtools::ScanBamParam(what="mapq", which=which)
    ga <- GenomicAlignments::readGAlignments(file, index=file,
                                             param=param, use.names=FALSE)

    if(length(ga) == 0) {
        stop( file, " bam file doesn't contain ",
              chr, " chromosome(s).\nExiting." )
    }

    bed <- GenomicRanges::granges(x = ga)

    return(bed)
}


#' read a bed file into a GenomicRanges like format.
#'
#' @param file Character indicating path to bam file.
#' @param chr Integer indicating which chromosome(s) to read in.
#' @return GenomicRanges object
#' @keywords internal
#' @import tools
#' @import GenomicRanges
#' @import rtracklayer
readBedFile <- function(filename, chr) {

    if (tools::file_ext(filename) == "zip") {
        tmp <- utils::unzip(filename, list=T)$Name
        file <- unz(filename, tmp)
    } else {
        file <- filename
    }

    bed <- rtracklayer::import.bed(con = file)
    bed <- bed[which(bed@seqnames %in% chr), ]
    if(length(bed) == 0) {
        stop( filename, " bed file doesn't contain ",
              chr, " chromosome(s).\nExiting." )
    }

    bed <- GRanges(seqnames=bed@seqnames,
                   ranges=bed@ranges,
                   strand=bed@strand)
    return(bed)
}

