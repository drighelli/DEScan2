#
# testPeaksManually <- function() {
#
#
# file="inst/extdata/Bed/P43615_Sample_FC1_Input_fwd_chr19_Smartfilter.bed.zip"
#
# chr="chr19"; filetype="bed"; fragmentLength=200;
# binSize=50; minWin=50; maxWin=1000; genomeName="mm9";
# minCompWinWidth=5000; maxCompWinWidth=10000;
# zthresh=5; minCount=0.1; verbose=TRUE; save=FALSE
# onlyStdChrs=TRUE
# fileGRangesList <- NULL
# winVector <- c( (minWin/50):(maxWin/50))
#
#
#     bedGRanges <- constructBedRanges(filename=file, filetype=filetype,
#                                      genomeName=genomeName,
#                                      onlyStdChrs=onlyStdChrs)
#
#     bedGrangesChrsList <- cutGRangesPerChromosome(bedGRanges)
#     if(!is.null(chr)) bedGrangesChrsList <- keepRelevantChrs(bedGrangesChrsList, chr)
#
#     if(verbose) message("Calling Peaks on chromosomes...")
#     chrGRanges=bedGrangesChrsList[[1]]
#
#     chrBedGRanges=chrGRanges
#     minWinWidth=minWin
#     maxWinWidth=maxWin
#     binWidth=binSize
#             runWinRleList <- computeCoverageMovingWindowOnChr(
#                 chrBedGRanges=chrGRanges,
#                 minWinWidth=minWin,
#                 maxWinWidth=maxWin,
#                 binWidth=binSize
#             )
#             minCompRunWinRleList <- computeCoverageMovingWindowOnChr(
#                 chrBedGRanges=chrGRanges,
#                 minWinWidth=minCompWinWidth,
#                 maxWinWidth=minCompWinWidth,
#                 binWidth=binSize
#             )
#             maxCompRunWinRleList <- computeCoverageMovingWindowOnChr(
#                 chrBedGRanges=chrGRanges,
#                 minWinWidth=maxCompWinWidth,
#                 maxWinWidth=maxCompWinWidth,
#                 binWidth=binSize
#             )
#             minChrRleWComp=minCompRunWinRleList[[1]]
#             maxChrRleWComp=maxCompRunWinRleList[[1]]
#             lambdaChrRleList <- computeLambdaOnChr(
#                 chrGRanges=chrGRanges,
#                 winVector=winVector,
#                 minChrRleWComp=minCompRunWinRleList[[1]],
#                 minCompWinWidth=minCompWinWidth,
#                 maxChrRleWComp=maxCompRunWinRleList[[1]],
#                 maxCompWinWidth=maxCompWinWidth
#             )
#             Z <- computeZ(lambdaChrRleList=lambdaChrRleList,
#                           runWinRleList=runWinRleList,
#                           chrLength=chrGRanges@seqinfo@seqlengths,
#                           minCount=minCount, binSize=binSize
#             )
#             newS <- get_disjoint_max_win(z0=Z[1:70000,],
#                                           sigwin=fragmentLength/binSize,
#                                           nmax=Inf, zthresh=zthresh,
#                                           two_sided=FALSE, verbose=verbose
#             )
#             chrZRanges <- createGranges(chrSeqInfo=chrGRanges@seqinfo,
#                                         starts=as.numeric(rownames(Z)[newS[,1]]),
#                                         widths=newS[,2]*binSize,
#                                         mcolname="z-score",
#                                         mcolvalues=newS[,3]
#             )
#
#
# }
#
#
# oldpeaks <- function()
# {
#     source("R/old_peaks.R")
#
#     file="testData/Bed/chr19/P43615_Sample_FC1_Input_fwd_chr19_Smartfilter.bed.zip"
#     chr = 19; filetype = "bed"; fraglen = 200;
#     min_win=1;
#     rlen = 100; min_bin = 50; max_win = 20; blocksize = 10000;
#     zthresh = 5; min_count = 0.1; verbose = FALSE; save = FALSE
#
#
#         chr <- paste0("chr", chr)
#         if (tools::file_ext(file) == "zip") {
#             tmp <- unzip(file, list = T)$Name
#             file <- unz(file, tmp)
#         }
#         bed <- utils::read.table(file, sep = "\t", header = FALSE)[, c(1:3, 6)]
#         colnames(bed) <- c("seqnames","start","end","strand")
#         bed <- bed[bed[,1] == chr,]
#
#     # grid = vector of integers describing the start coordinates for each bin
#     grid <- make_grid(bed, fraglen, rlen, min_bin = min_bin)
#     ngrid <- length(grid)
#     gsize <- grid[2] - grid[1]
#
#     fr <- sort(bed[which(bed[, 4] == "+"), 2]) ## only starts
#     rr <- sort(bed[which(bed[, 4] == "-"), 2]) ## only starts
#
#     tot_rds <- length(fr) + length(rr)
#     tot_base <- max(rr[length(rr)], fr[length(fr)]) - min(fr[1], rr[1]) ## max start - min start
#
#     # do analysis blockwise due to memory issues
#     # fill initial block
#     block <- c(1, blocksize + max_win)
#     finalblock <- FALSE
#     s <- matrix(ncol = 3, nrow = 0, data = 0)
#         blocksize_i <- blocksize
#
#         # conditions to stop once as end of grid has been reached
#             cat("Doing block ", block[1], " to ", block[2], " out of ", ngrid, "\n")
#
#         grid0 <- grid[block[1]:block[2]]
#
#
#         # cmat = coverage matrix with bins as rows and window sizes as columns
#         cmat <- window_coverage_old(Fr = fr, Rr = rr, grid = grid, fraglen = fraglen,
#                                     rlen = rlen, max_win = max_win, min_win = min_win,
#                                     verbose = FALSE)
#
#         # Get lambdalocal.
#
#         # compute coverage using a 5kb window
#         headstart <- ceiling(5000 / gsize) * gsize
#         # grid05k <- c(seq(grid0[1] - headstart, grid0[1], by = gsize), grid0)
#         grid05k <-  grid0
#         offset5k <- length(grid05k) - length(grid0)
#         c5k <- window_coverage_old(fr, rr, grid05k, fraglen = fraglen, rlen = rlen,
#                                    max_win = 5000, min_win = 5000, verbose = FALSE)
#         # compute coverage using a 10kb window
#         headstart <- ceiling(10000 / gsize) * gsize
#         # grid010k <- c(seq(grid0[1] - headstart, grid0[1], by = gsize), grid0)
#         grid010k <- grid0
#         offset10k <- length(grid010k) - length(grid0)
#         c10k <- window_coverage_old(fr, rr, grid010k, fraglen = fraglen,
#                                     rlen = rlen, max_win = 10000,
#                                     min_win = 10000, verbose = FALSE)
#         # determine lambda for 5k, 10k windows and baseline
#         lamloc <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
#         lam5k <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
#         lam10k <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
#         lambl <- matrix(nrow = nrow(cmat), ncol = ncol(cmat), data = 0)
#         for (win in min_win:max_win) {
#             # lam5k[, win - min_win + 1] <- c5k[offset5k + c(1:length(grid0)) - floor(win / 2)] * win / 5000
#             # lam10k[, win - min_win + 1] <- c10k[offset10k + c(1:length(grid0)) - floor(win / 2)] * win / 10000
#             lam5k[, win - min_win + 1] <- c5k[c(1:length(grid0))] * win / 5000
#             lam10k[, win - min_win + 1] <- c10k[c(1:length(grid0))] * win / 10000
#             lambl[, win - min_win + 1] <- tot_rds * win / tot_base
#         }
#         lamloc <- pmax(lam5k, lam10k, lambl)
#         # calculate z score for each bin x window combination
#         z <- sqrt(2) * sign(cmat - lamloc) *
#             sqrt(cmat * log(pmax(cmat, min_count) / lamloc) - (cmat - lamloc))
#         if(length(which(z>100)) >0 ) print(z[which(z>100, arr.ind=TRUE)])
#         # find high z scores keeping one with no intersecting other bin/windows
#         new_ss <- get_disjoint_max_win_old(z0 = z[1:blocksize_i, ],
#                                       sigwin = fraglen / gsize, nmax = Inf,
#                                       zthresh = zthresh, two_sided = FALSE,
#                                       verbose = FALSE)
#         # convert new_s bins and width into genomic coordinates and append to s
#         new_s=matrix(nrow=nrow(new_ss), ncol=ncol(new_ss))
#         if (nrow(new_ss) >= 1) {
#             new_s[, 1] <- new_ss[, 1] + block[1] - 1
#             new_s <- cbind(grid[new_ss[, 1, drop = FALSE]],
#                            ## when the window is equal to 1 then the end is equal to the start,
#                            ## but the window should be multiplicated by the binsize
#                            grid[new_ss[, 1, drop = FALSE] + new_ss[, 2, drop = FALSE] - 1],
#                            new_ss[, 3, drop = FALSE])
#             s <- rbind(s, new_s)
#         }
#
#
#
# }
#
#
