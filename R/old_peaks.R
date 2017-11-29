
find_Peaks_Old <- function(files, filetype = "bam", chr = 1:19, fraglen = 200,
                      rlen = 100, min_bin = 50, min_win = 1, max_win = 20, blocksize = 10000,
                      zthresh = 5, min_count = 0.1, output_name = "Peaks", save = TRUE, verbose = FALSE) {
    for (file in files) {
        if (is.character(chr) != TRUE) {
            chr <- paste0("chr", chr)
        }
        if (filetype == "bam") {
            bed <- read_bam(file, chr)
        }
        if (filetype == "bed") {
            if (tools::file_ext(file) == "zip") {
                tmp <- unzip(file, list = T)$Name
                file <- unz(file, tmp)
            }
            bed <- utils::read.table(file, sep = "\t", header = FALSE)[, c(1:3, 6)]
            colnames(bed) <- c("seqnames","start","end","strand")
            bed <- bed[bed[,1] == chr,]
        }

        # grid = vector of integers describing the start coordinates for each bin
        grid <- make_grid(bed, fraglen, rlen, min_bin = min_bin)
        ngrid <- length(grid)
        gsize <- grid[2] - grid[1]

        fr <- sort(bed[which(bed[, 4] == "+"), 2])
        rr <- sort(bed[which(bed[, 4] == "-"), 2])

        tot_rds <- length(fr) + length(rr)
        tot_base <- max(rr[length(rr)], fr[length(fr)]) - min(fr[1], rr[1])

        # do analysis blockwise due to memory issues
        # fill initial block
        block <- c(1, blocksize + max_win)
        finalblock <- FALSE
        s <- matrix(ncol = 3, nrow = 0, data = 0)
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
                cat("Doing block ", block[1], " to ", block[2], " out of ", ngrid, "\n")
            }
            grid0 <- grid[block[1]:block[2]]

            if (verbose) {
                cat("\tComputing coverage stats...\n")
            }
            # cmat = coverage matrix with bins as rows and window sizes as columns
            cmat <- window_coverage_old(Fr = fr, Rr = rr, grid = grid0, fraglen = fraglen,
                                    rlen = rlen, max_win = max_win, min_win = min_win,
                                    verbose = FALSE)

            # Get lambdalocal.

            # compute coverage using a 5kb window
            headstart <- ceiling(5000 / gsize) * gsize
            grid05k <- c(seq(grid0[1] - headstart, grid0[1], by = gsize), grid0)
            offset5k <- length(grid05k) - length(grid0)
            c5k <- window_coverage_old(fr, rr, grid05k, fraglen = fraglen, rlen = rlen,
                                   max_win = 5000, min_win = 5000, verbose = FALSE)
            # compute coverage using a 10kb window
            headstart <- ceiling(10000 / gsize) * gsize
            grid010k <- c(seq(grid0[1] - headstart, grid0[1], by = gsize), grid0)
            offset10k <- length(grid010k) - length(grid0)
            c10k <- window_coverage_old(fr, rr, grid010k, fraglen = fraglen,
                                    rlen = rlen, max_win = 10000,
                                    min_win = 10000, verbose = FALSE)
            # determine lambda for 5k, 10k windows and baseline
            lamloc <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
            lam5k <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
            lam10k <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
            lambl <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
            for (win in min_win:max_win) {
                lam5k[, win - min_win + 1] <- c5k[offset5k + c(1:length(grid0)) -
                                                      floor(win / 2)] * win / 5000
                lam10k[, win - min_win + 1] <- c10k[offset10k + c(1:length(grid0)) -
                                                        floor(win / 2)] * win / 10000
                lambl[, win - min_win + 1] <- tot_rds * win / tot_base
            }
            lamloc <- pmax(lam5k, lam10k, lambl)
            # calculate z score for each bin x window combination
            z <- sqrt(2) * sign(cmat - lamloc) *
                sqrt(cmat * log(pmax(cmat, min_count) / lamloc) - (cmat - lamloc))

            # find high z scores keeping one with no intersecting other bin/windows
            new_s <- get_disjoint_max_win_old(z0 = z[1:blocksize_i, ],
                                          sigwin = fraglen / gsize, nmax = Inf,
                                          zthresh = zthresh, two_sided = FALSE,
                                          verbose = FALSE)
            # convert new_s bins and width into genomic coordinates and append to s
            if (nrow(new_s) >= 1) {
                new_s[, 1] <- new_s[, 1] + block[1] - 1
                new_s <- cbind(grid[new_s[, 1, drop = FALSE]],
                               grid[new_s[, 1, drop = FALSE] +
                                        new_s[, 2, drop = FALSE] - 1],
                               new_s[, 3, drop = FALSE])
                s <- rbind(s, new_s)
            }

            elapsed <- proc.time() - ptm
            if (verbose) {
                cat("\tDone. That took ", format(elapsed[3] / 60, digits = 1),
                    " minutes.\n")
            }
            if (finalblock) break

            # shift block by blocksize and repeat
            block <- c(block[1] + blocksize_i, block[2] + blocksize_i)
        }

        peaks <- cbind(rep(chr, dim(s)[1]), s)
        if (save == TRUE) {
            if (dir.exists(output_name) == FALSE) {
                dir.create(output_name)
            }
            if (dir.exists(paste0(output_name,"/", chr)) == FALSE) {
                dir.create(paste0(output_name,"/", chr))
            }
            fname <- strsplit(basename(file), split = ".", fixed = TRUE)[[1]][1]
            fileprefix <- paste0(output_name,"/", chr, "/Peaks_", fname)


            save(peaks, fraglen, rlen, zthresh, min_win, max_win,
                 file = paste0(fileprefix, ".RData"))
        }
    }
    return(peaks)
}

make_grid <- function(bed, fraglen, rlen, min_bin = 10, min_rds = 0,
                      grid_st = NA, grid_ed = NA, verbose=FALSE) {
    fr <- sort(bed[which(bed[, 4] == "+"), 2])
    rr <- sort(bed[which(bed[, 4] == "-"), 2])

    if (is.na(grid_st)) {
        grid_st <- min(fr, (rr + rlen - fraglen))
    }
    if (is.na(grid_ed)) {
        grid_ed <- max((fr + fraglen), (rr + rlen))
    }

    # Initial grid is even spaced:
    grid0 <- seq(grid_st, grid_ed, min_bin)

    # Needs work
    if (min_rds > 0) { # this portion of the code does not work!
        # Get the number of fragments that cover each bin.
        # fragct[i] is the number of fragments with positive overlap with
        # grid0[i]:grid0[i+1]
        ffragct <- matrix(nrow = length(grid0), ncol = length(fr), data = 0)
        # Error: cannot allocate vector of size 2056.4 Gb #######################
        rfragct <- matrix(nrow = length(grid0), ncol = length(fr), data = 0)
        for (k in seq_along(fr)) {
            for (i in 1:length(fr[[k]])) {
                if (i %% 10000 == 0) cat("Sample ", k, "Forward reads : ", i, "\n")
                fragst <- fr[[k]][i]
                fraged <- fr[[k]][i] + fraglen - 1
                binst <- floor((fragst - grid_st) / min_bin) + 1
                bined <- floor((fraged - grid_st) / min_bin) + 1
                ffragct[binst:bined, k] <- ffragct[binst:bined, k] + 1
            }
        }
        for (k in 1:length(rr)) {
            for (i in 1:length(rr[[k]])) {
                if (i %% 10000 == 0) cat("Sample ", k, "Reverse reads : ", i, "\n")
                fragst <- rr[[k]][i] + rlen - fraglen
                fraged <- rr[[k]][i] + rlen - 1
                binst <- floor((fragst - grid_st) / min_bin) + 1
                bined <- floor((fraged - grid_st) / min_bin) + 1
                rfragct[binst:bined, k] <- rfragct[binst:bined, k] + 1
            }
        }

        # Starting from the left most bin, join bins whose total fragment coverage
        # across all samples is less than min_rds.
        totcov <- rowSums(ffragct) + rowSums(rfragct)
        remove_tick <- rep(FALSE, length(totcov))
        cursum <- 0
        for (i in 1:length(totcov)) {
            if (i %% 100000 == 0 && verbose) cat(i, " out of ", length(totcov), "\n")
            cursum <- cursum + totcov[i]
            if (cursum < min_rds) {
                remove_tick[i] <- TRUE
            } else {
                cursum <- 0
            }
        }
        grid <- grid0[which(!remove_tick)]
    } else {
        grid <- grid0
    }
    grid
}

window_coverage_old <- function(Fr, Rr, grid, fraglen = 200, rlen = 100,
                            min_win = 1, max_win = 100,
                            verbose = TRUE) {
    grid_st <- min(grid)
    bin_size <- grid[2] - grid[1] # Assume even grid. ! only if min_rds = 0
    fragbins <- fraglen / bin_size # the number of bins covered by each fragment
    ngrid <- length(grid)

    # initialize fragment count matrices for forward and reverse reads
    ffragctmat <- matrix(nrow = length(grid), ncol = max_win - min_win + 1,
                         data = 0)
    rfragctmat <- matrix(nrow = length(grid), ncol = max_win - min_win + 1,
                         data = 0)

    # index forward reads within the current grid
    forward_index <- which(Fr > grid[1] - fraglen & Fr < grid[length(grid)])

    # for each forward read within the given grid convert the base position to
    # bin and add +1 count to the fragment count matrix (bins as rows, window
    # size as columns) while moving window right along bins and extending left
    if (length(forward_index) > 0) {
        for (i in 1:length(forward_index)) {
            if (verbose && i %% 1000 == 0) {
                cat("Forward reads: ", i, "out of", length(forward_index), "\n")
            }
            fragst <- Fr[forward_index[i]]  # start of fragment
            # which bin contains start of fragment
            binst <- floor((fragst - grid_st) / bin_size) + 1
            for (j in min_win:max_win) {
                st <- max(binst - j + 1, 1) # lower limit adj by current window size
                ed <- min(binst + fragbins - 1, ngrid) # upper limit
                ffragctmat[st:ed, j - min_win + 1] <- ffragctmat[st:ed, j -
                                                                     min_win + 1] + 1
            }
        }
    }

    reverse_index <- which(Rr > grid[1] - rlen &
                               (Rr + rlen - fraglen) < grid[length(grid)])
    if (length(reverse_index) > 0) {
        for (i in 1:length(reverse_index)) {
            if (verbose && i %% 1000 == 0) {
                cat("Reverse reads: ", i, "out of", length(reverse_index), "\n")
            }
            fragst <- Rr[reverse_index[i]] + rlen - fraglen  # start of fragment
            # bin containing start of fragment
            binst <- floor((fragst - grid_st) / bin_size) + 1
            for (j in min_win:max_win) {
                st <- max(binst - j + 1, 1)
                ed <- min(binst + fragbins - 1, ngrid)
                rfragctmat[st:ed, j - min_win + 1] <- rfragctmat[st:ed, j -
                                                                     min_win + 1] + 1
            }
        }
    }
    cmat <- ffragctmat + rfragctmat
    cmat
}

get_disjoint_max_win_old <- function(z0, sigwin = 20, nmax = Inf,
                                 zthresh = -Inf, two_sided = FALSE,
                                 verbose = TRUE) {
    s <- matrix(ncol = 3, nrow = 0)
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
        if (z0[t, w] < zthresh) break # break loop once as max z below thresh
        s <- rbind(s, c(t, w, z0[t, w]))
        if (verbose) {
            cat("Maximizing window: ", t, ",", w, " Score=", z0[t, w], "\n")
        }
        # remove surrounding max window around max z from consideration
        st <- max(1, t - sigwin - maxwin + 1)
        ed <- min(t + w + sigwin - 1, nrow(z0))
        z0[st:ed, ] <- -Inf

        if (nrow(s) >= nmax) break
    }
    s
}

read_bam <- function(file, chr) {

    if (file.exists(paste0(file,".bai")) == FALSE) {
        Rsamtools::indexBam(file)
    }
    bf <- Rsamtools::BamFile(file)
    sl <- GenomeInfoDb::seqlengths(Rsamtools::seqinfo(bf))

    cat("Reading:", file, chr, "\n")
    which <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, sl[chr]))
    param <- Rsamtools::ScanBamParam(what = "mapq", which = which)
    ga <- GenomicAlignments::readGAlignments(file, index = file,
                                             param = param, use.names = FALSE)
    bed <- as.data.frame(ga)[, c("seqnames", "start", "end" ,"strand")]
    bed <- droplevels.data.frame(bed)
    # bed$start <- bed$start - 1 #1-based indexing used instead of 0-based
    bed
}
