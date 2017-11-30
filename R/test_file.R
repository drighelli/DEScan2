
testPeaksManually <- function() {


file="inst/extdata/Bed/P43615_Sample_FC1_Input_fwd_chr19_Smartfilter.bed.zip"

chr="chr19"; filetype="bed"; fragmentLength=200;
binSize=50; minWin=1; maxWin=20; genomeName="mm9";
minCompWinWidth=5000; maxCompWinWidth=10000;
zthresh=5; minCount=0.1; verbose=TRUE; save=FALSE
onlyStdChrs=TRUE
fileGRangesList <- NULL
winVector <- c(minWin:maxWin)


    bedGRanges <- constructBedRanges(filename=file, filetype=filetype,
                                     genomeName=genomeName,
                                     onlyStdChrs=onlyStdChrs)

    bedGrangesChrsList <- cutGRangesPerChromosome(bedGRanges)
    if(!is.null(chr)) bedGrangesChrsList <- keepRelevantChrs(bedGrangesChrsList, chr)

    if(verbose) message("Calling Peaks on chromosomes...")
    chrGRanges=bedGrangesChrsList[[1]]


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
            Z <- computeZ(lambdaChrRleList=lambdaChrRleList,
                          runWinRleList=runWinRleList,
                          chrLength=chrGRanges@seqinfo@seqlengths,
                          minCount=minCount, binSize=binSize
            )
            newS <- get_disjoint_max_win(z0=Z[1:70000,],
                                          sigwin=fragmentLength/binSize,
                                          nmax=Inf, zthresh=zthresh,
                                          two_sided=FALSE, verbose=verbose
            )
            chrZRanges <- createGranges(chrSeqInfo=chrGRanges@seqinfo,
                                        starts=as.numeric(rownames(z)[newS[,1]]),
                                        widths=newS[,2]*binSize,
                                        mcolname="z-score",
                                        mcolvalues=newS[,3]
            )


}

# chr = 19; filetype = "bed"; fraglen = 200;
# rlen = 100; min_bin = 50; max_win = 20; blocksize = 10000;
# zthresh = 5; min_count = 0.1; verbose = FALSE; save = FALSE
