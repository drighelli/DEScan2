# library("DEScan")
context("Peak Finding")

test_that("Test if new findPeaks works with bam and bed files", {
    bam.path <- system.file(file.path("extdata","tests","inputdata","bam"),
                            package="DEScan2")
    bam.file <- list.files(bam.path, full.names=TRUE, pattern=".bam$")

    binSize=50
    minWin=50
    maxWin=100
    zthresh=10
    minCount=0.1
    minCompWinWidth=5000
    maxCompWinWidth=10000
    outputName="Peaks"
    savef=FALSE
    verb=FALSE
    onlyStdChrs=TRUE
    chr=NULL
    filetype="bam"
    sigwin=10
    osc=TRUE
    genName="mm9"
    if( length(bam.file) != 0 )
    {
        bampeaksGRL <- findPeaks(files=bam.file, filetype=filetype,
                            binSize=binSize,
                            minWin=minWin, maxWin=maxWin,
                            zthresh=zthresh, minCount=minCount,
                            minCompWinWidth=minCompWinWidth,
                            maxCompWinWidth=maxCompWinWidth,
                            save=savef, verbose=verb,genomeName=genName,
                            sigwin=sigwin, onlyStdChrs=osc)
        names(bampeaksGRL) <- "FC1"

        grl.path <- system.file(
                            file.path("extdata","tests","peaks","FC1_GRL.rds"),
                            package="DEScan2")
        bampeaksRef <- readRDS(grl.path)
        expect_identical(bampeaksGRL, bampeaksRef)
    }
    else
    {
        warning("bam file does not exist!")
    }

    bedpath <- system.file(file.path("extdata","tests","inputdata","bed"),
                                package="DEScan2")
    bedfile <- list.files(bedpath, full.names=TRUE, pattern=".bed.zip$")
    if( file.exists(bedfile) ) {
        bedpeaksGRL <- findPeaks(files=bedfile, filetype="bed",
                                 binSize=binSize,
                                 minWin=minWin, maxWin=maxWin,
                                 zthresh=zthresh, minCount=minCount,
                                 minCompWinWidth=minCompWinWidth,
                                 maxCompWinWidth=maxCompWinWidth,
                                 save=savef, verbose=verb,genomeName=genName,
                                 sigwin=sigwin, onlyStdChrs=osc)
        names(bedpeaksGRL) <- "FC1"
        expect_identical(bedpeaksGRL, bampeaksGRL)
    } else {
        warning("bed file does not exist!")
    }

})

# test_that("Test disjoint function R & C", {
#     zpath <- system.file(file.path("extdata","tests","z"), package="DEScan2")
#     zfile <- list.files(zpath, full.names=TRUE)
#     zzz <- readRDS(zfile)
#
#     sigw=10
#     zthr=10
#     nmax=1e5
#     verb=TRUE
#     cTime <- system.time({
#         smatC <- c_get_disjoint_max_win(z0=zzz, sigwin=sigw, zthresh=zthr,
#                                                nmax=nmax, verbose=verb)
#     })
#     rTime <- system.time({
#         smatR <- get_disjoint_max_win(z0=zzz, sigwin=sigw, zthresh=zthr,
#                                       nmax=nmax, verbose=verb)
#     })
#     message("cTime (secs): ", round(cTime[1], 3),
#             " rTime (secs): ", round(rTime[1], 3))
#     expect_identical(smatC, smatR)
# })
