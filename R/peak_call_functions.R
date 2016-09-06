findPeaks <- function(files, chr = 1:19, fraglen = 200, rlen = 100, min_win = 1,
                      max_win = 100, blocksize = 10000, zthresh = 5,
                      min_count = 0.1, grid = NA, filetype = "bam") {
  for (file in files) {

    if (is.character(chr) != T) {
      chr <- paste0("chr", chr)
    }
    if (filetype == "bam") {
      bed <- read_bam(file, chr)
    }
    if (filetype == "bed") {
      bed <- read.table(file, sep = "\t", header = F)[, c(1:3, 6)]
      colnames(bed) <- c("seqnames","start","end","strand")
      bed <- bed[bed[,1] == chr,]
    }

    grid <- makeGrid(bed, fraglen, rlen, min_bin = min_win * 50)
    ngrid <- length(grid)
    gsize <- grid[2] - grid[1]

    fr = sort(bed[which(bed[,4] == "+"),2])
    rr = sort(bed[which(bed[,4] == "-"),2])
    tot_rds <- length(fr) + length(rr)
    tot_base <- max(rr[length(rr)], fr[length(fr)]) - min(fr[1], rr[1])

    block <- c(1, blocksize + max_win)
    finalblock <- FALSE
    s <- matrix(ncol = 3, nrow = 0, data = 0)
    while (TRUE) {
      ptm <- proc.time()

      if (block[1] > ngrid) break
      if (block[2] >= ngrid) {
        block[2] <- min(block[2], ngrid)
        blocksize <- block[2] - block[1] + 1
        finalblock <- TRUE
      }
      cat("Doing block ", block[1], " to ", block[2], " out of ", ngrid, "\n")
      grid0 <- grid[block[1]:block[2]]

      cat("\tComputing coverage stats...\n")
      c <- window_coverage(fr, rr, grid0, fraglen = fraglen, rlen = rlen,
                           max_win = max_win, min_win = min_win,
                           verbose = FALSE)

      # Get lambdalocal.
      headstart <- ceiling(5000 / gsize) * gsize
      grid05k <- c(seq(grid0[1] - headstart, grid0[1], by = gsize), grid0)
      offset5k <- length(grid05k) - length(grid0)
      c5k <- window_coverage(fr, rr, grid05k, fraglen = fraglen,
                                    rlen = rlen, max_win = 5000, min_win = 5000,
                                    verbose = FALSE)
      headstart <- ceiling(10000 / gsize) * gsize
      grid010k <- c(seq(grid0[1] - headstart, grid0[1], by = gsize), grid0)
      offset10k <- length(grid010k) - length(grid0)
      c10k <- window_coverage(fr, rr, grid010k, fraglen = fraglen,
                                     rlen = rlen, max_win = 10000,
                                     min_win = 10000, verbose = FALSE)
      lamloc <- matrix(nrow = nrow(c), ncol = ncol(c), data = 0)
      lam5k <- matrix(nrow = nrow(c), ncol = ncol(c), data = 0)
      lam10k <- matrix(nrow = nrow(c), ncol = ncol(c), data = 0)
      lambl <- matrix(nrow = nrow(c), ncol = ncol(c), data = 0)
      for (win in min_win:max_win) {
        lam5k[, win - min_win + 1] <- c5k[offset5k + c(1:length(grid0)) -
                                          floor(win / 2)] * win / 5000
        lam10k[, win - min_win + 1] <- c10k[offset10k + c(1:length(grid0)) -
                                            floor(win / 2)] * win / 10000
        lambl[, win - min_win + 1] <- tot_rds * win / tot_base
      }
      lamloc <- pmax(lam5k, lam10k, lambl)
      z <- sqrt(2) * sign(c - lamloc) *
        sqrt(c * log(pmax(c, min_count) / lamloc) - (c - lamloc))

      new_s <- get_disjoint_max_win(z0 = z[1:blocksize, ],
                                          sigwin = fraglen / gsize, nmax = 200,
                                          zthresh = zthresh, two_sided = FALSE,
                                          verbose = FALSE)
      if (nrow(new_s) >= 1) {
        new_s[, 1] <- new_s[, 1] + block[1] - 1
        new_s <- cbind(grid[new_s[, 1, drop = FALSE]],
                     grid[new_s[, 1, drop = FALSE] +
                            new_s[, 2, drop = FALSE] - 1],
                     new_s[, 3, drop = FALSE])
        s <- rbind(s, new_s)
      }

      elapsed <- proc.time() - ptm
      cat("\tDone. That took ", format(elapsed[3] / 60, digits = 1),
          " minutes.\n")
      if (finalblock) break
      block <- c(block[1] + blocksize, block[2] + blocksize)
    }

    if (dir.exists("Peaks") == F) {
      dir.create("Peaks")
    }
    if (dir.exists(paste0("Peaks/", chr)) == F) {
      dir.create(paste0("Peaks/", chr))
    }
    fname <- strsplit(basename(file), split = ".", fixed = T)[[1]][1]
    fileprefix = paste0("Peaks/", chr, "/Peaks_", fname, "_", chr)

    s <- cbind(rep(chr, dim(s)[1]), s)
    save(s, fraglen, rlen, zthresh, min_win, max_win,
         file = paste0(fileprefix, ".RData"))
  }
}

read_bam <- function(file, chr) {

  if (file.exists(paste0(file,".bai")) == F) {
    Rsamtools::indexBam(file)
  }
  bf <- Rsamtools::BamFile(file)
  sl <- GenomeInfoDb::seqlengths(Rsamtools::seqinfo(bf))

  cat("Reading:", file, chr, "\n")
  which <- GenomicRanges::GRanges(chr, IRanges::IRanges(1, sl[chr]))
  param <- Rsamtools::ScanBamParam(what = "mapq", which = which)
  ga <- GenomicAlignments::readGAlignments(file, index = file,
                                           param = param, use.names = F)
  bed <- as.data.frame(ga)[, c("seqnames", "start", "end" ,"strand")]
  bed <- droplevels.data.frame(bed)
  bed$start <- bed$start - 1 #1-based indexing used instead of 0-based
  bed
}

window_coverage <- function(Fr, Rr, grid, fraglen = 200, rlen = 100,
                                    min_win = 1, max_win = 100,
                                    verbose = TRUE) {
  grid_st <- min(grid)
  bin_size <- grid[2] - grid[1] # Assume even grid.
  fragbins <- fraglen / bin_size # the number of bins covered by each fragment.
  ngrid <- length(grid)

  ffragctmat <- matrix(nrow = length(grid), ncol = max_win - min_win + 1,
                       data = 0)
  rfragctmat <- matrix(nrow = length(grid), ncol = max_win - min_win + 1,
                       data = 0)
  tot <- 0

  ix <- which(Fr > grid[1] - fraglen & Fr < grid[length(grid)])
  if (length(ix) > 0) {
    for (i in 1:length(ix)) {
      if (verbose && i %% 1000 == 0) {
        cat("Forward reads: ", i, "out of", length(ix), "\n")
      }
      fragst <- Fr[ix[i]]  # start of fragment
      # bin containing start of fragment
      binst <- floor((fragst - grid_st) / bin_size) + 1
      for (j in min_win:max_win) {
        st <- max(binst - j + 1, 1)
        ed <- min(binst + fragbins - 1, ngrid)
        ffragctmat[st:ed, j - min_win + 1] <- ffragctmat[st:ed, j -
                                                           min_win + 1] + 1
      }
    }
  }
  tot <- tot + length(ix)
  ix <- which(Rr > grid[1] - rlen & (Rr + rlen - fraglen) < grid[length(grid)])
  if (length(ix) > 0) {
    for (i in 1:length(ix)) {
      if (verbose && i %% 1000 == 0) {
        cat("Reverse reads: ", i, "out of", length(ix), "\n")
      }
      fragst <- Rr[ix[i]] + rlen - fraglen  # start of fragment
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
  tot <- tot + length(ix)
  c <- ffragctmat + rfragctmat
  c
}

get_disjoint_max_win <- function(z0, sigwin = 20, nmax = 20,
                                         zthresh = -Inf, two_sided = FALSE,
                                         verbose = TRUE) {
  s <- matrix(ncol = 3, nrow = 0)
  maxwin <- ncol(z0)
  if (two_sided) {
    z0 <- abs(z0)
  }
  while (TRUE) {
    inds <- which.max(z0)
    if (length(inds) == 0) break
    w <- ceiling(inds / nrow(z0))
    t <- inds %% nrow(z0)
    if (t == 0) t <- nrow(z0)
    if (z0[t, w] < zthresh) break
    s <- rbind(s, c(t, w, z0[t, w]))
    if (verbose) {
      cat("Maximizing window: ", t, ",", w, " Score=", z0[t, w], "\n")
    }
    for (i in 1:maxwin) {
      st <- max(1, t - sigwin - i + 1)
      ed <- min(t + w + sigwin - 1, nrow(z0))
      z0[st:ed, ] <- -Inf
    }
    if (nrow(s) >= nmax) break
  }
  s
}

makeGrid <- function(bed, fraglen, rlen, min_bin = 10, min_rds = 0,
                     grid_st = NA, grid_ed = NA) {

  fr = sort(bed[which(bed[,4] == "+"),2])
  rr = sort(bed[which(bed[,4] == "-"),2])

  if (is.na(grid_st)) {
    grid_st <- Inf
    for (i in 1:length(fr)) grid_st <- min(grid_st, fr[[i]][1])
    for (i in 1:length(rr)) grid_st <- min(grid_st, rr[[i]][1] + rlen - fraglen)
  }
  if (is.na(grid_ed)) {
    grid_ed <- -Inf
    for (i in 1:length(fr)) {
      grid_ed <- max(grid_ed, fr[[i]][length(fr[[i]])] + fraglen)
    }
    for (i in 1:length(rr)) {
      grid_ed <- max(grid_ed, rr[[i]][length(rr[[i]])] + rlen)
    }
  }

  # Initial grid is even spaced:
  grid0 <- seq(grid_st, grid_ed, min_bin)

  if (min_rds > 0) {
    # Get the number of fragments that cover each bin.
    # fragct[i] is the number of fragments with positive overlap with
    # grid0[i]:grid0[i+1]
    ffragct <- matrix(nrow = length(grid0), ncol = length(fr), data = 0)
    rfragct <- matrix(nrow = length(grid0), ncol = length(fr), data = 0)
    for (k in 1:length(fr)) {
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
      if (i %% 100000 == 0) cat(i, " out of ", length(totcov), "\n")
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
