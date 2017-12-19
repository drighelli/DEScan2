# library("DEScan")
context("Peak Finding")

test_that("Test if new findPeaks works with bam and bed files", {

    bamfile <- "testData/bams/head/H_P43615_Sample_FC4_Input_fwd_chr19_Smartfilter.bam"

    binSize=50
    minWin=50
    maxWin=1000
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
    if( length(bamfile) != 0 )
    {
        bampeaksGRL <- findPeaks(files=bamfile, filetype=filetype, binSize=binSize,
                            minWin=minWin, maxWin=maxWin,
                            zthresh=zthresh, minCount=minCount,
                            minCompWinWidth=minCompWinWidth,
                            maxCompWinWidth=maxCompWinWidth,
                            save=savef, verbose=verb,genomeName=genName,
                            sigwin=sigwin, onlyStdChrs=osc)
        bampeaksRef <- readRDS("testData/peaks/H_P43615_Sample_FC1_Input_fwd_chr19_SmartfilterGRL.RDS")
        expect_identical(bampeaksGRL, bampeaksRef)
    }
    else
    {
        warning("bam file does not exist!")
    }
    bedfile <- "testData/bed/head/P43615_Sample_FC1_Input_fwd_chr19_Smartfilter.bam.bed.bed.zip"
    if( file.exists(bedfile) ) {
        bedpeaksGRL <- findPeaks(files=bedfile, filetype="bed", binSize=binSize,
                                 minWin=minWin, maxWin=maxWin,
                                 zthresh=zthresh, minCount=minCount,
                                 minCompWinWidth=minCompWinWidth,
                                 maxCompWinWidth=maxCompWinWidth,
                                 save=savef, verbose=verb,genomeName=genName,
                                 sigwin=sigwin, onlyStdChrs=osc)
        expect_identical(bedpeaksGRL, bampeaksGRL)
    } else {
        warning("bed file does not exist!")
    }

})

test_that("Test disjoint function R & C", {
    zzz <- readRDS("testData/z/z_matrix.rds")
    sigw=10
    zthr=10
    nmax=Inf
    verb=TRUE
    cTime <- system.time({
        smatC <- c_get_disjoint_max_win(z0=zzz, sigwin=sigw, zthresh=zthr,
                                               nmax=nmax, verbose=verb)
    })
    rTime <- system.time({
        smatR <- get_disjoint_max_win(z0=zzz, sigwin=sigw, zthresh=zthr,
                                      nmax=nmax, verbose=verb)
    })
    message("cTime (mins): ", cTime/60, " rTime (mins): ", rTime/60)
    expect_identical(smatC, smatR)
})
