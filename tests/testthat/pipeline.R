test_that("finalRegions", {
    library("devtools")
    document()
    load_all()
    bam.files <- list.files(path="testData/bams", pattern="bam", full.names = T)
    bam.files <- bam.files[-grep(pattern=".bai", x=bam.files)]

    peaksGRL <- findPeaks(files=bam.files, chr="chr19", filetype="bam",
                        fragmentLength=200,
                        binSize=50, minWin=50, maxWin=1000, genomeName="mm9",
                        minCompWinWidth=5000, maxCompWinWidth=10000,
                        zthresh=5, minCount=0.1, verbose=TRUE, save=FALSE)
    saveRDS(object=peaksGRL, file="peaksGRL_12_13_17.rds")

    regionsGR <- finalRegions(peakSamplesGRangesList=peaksGRL,
                                zThreshold=20, minCarriers=4,
                                saveFlag=FALSE, verbose=TRUE)
    saveRDS(object=regionsGR, file="regionsGR_12_13_17.rds")
    # regionsGR <- load()
    bamFilePath <- system.file("testData/bams", package="DEScan2")
    finalRegions <- DEScan2::countFinalRegions(regionsGRanges=regionsGR,
                                                readsFilePath="testData/bams",
                                                fileType="bam",
                                                minCarriers=1,
                                                genomeName="mm9",
                                                onlyStdChrs=TRUE,
                                                ignStrandSO=TRUE)
    saveRDS(object=finalRegions, file="finalRegions_12_13_17.rds")
})
