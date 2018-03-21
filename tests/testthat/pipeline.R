# test_that("whole pipeline", {
#     library("devtools")
#
#     document()
#     load_all()
#     library("DEScan2")
#     bam.files <- list.files(path="inst/extdata/bam/", pattern="bam$",
#                             full.names = TRUE)
#
#     # files=bam.files; chr="chr19"; filetype="bam";
#     # fragmentLength=200;
#     # binSize=50; minWin=50; maxWin=1000; genomeName="mm9";
#     # minCompWinWidth=5000; maxCompWinWidth=10000;
#     # zthresh=5; minCount=0.1; verbose=TRUE; save=TRUE;
#     # outputFolder="testData/Peaks/new_peaks_head"
#     # onlyStdChrs=TRUE
#
#     peaksGRL <- findPeaks(files=bam.files, filetype="bam",
#                               genomeName="mm9",
#                               binSize=50, minWin=50, maxWin=1000,
#                               zthresh=5, minCount=0.1, sigwin=10,
#                               minCompWinWidth=5000, maxCompWinWidth=10000,
#                               save=FALSE,
#                               outputFolder="testData/Peaks/new_peaks_head",
#                               force=TRUE,
#                               onlyStdChrs=TRUE,
#                               chr=NULL,
#                               verbose=TRUE)
#
#     saveRDS(object=peaksGRL, file="inst/extdata/peaks/RData/peaksGRL_all_files.rds")
#     regionsGR <- finalRegions(peakSamplesGRangesList=peaksGRL,
#                                 zThreshold=10, minCarriers=3,
#                                 saveFlag=FALSE,
#                                 outputFolder="testData/regions/prova",
#                                 verbose=TRUE)
#     saveRDS(object=regionsGR, file=paste0("inst/extdata/regions/regions.rds"))
#
#     # regionsGR <- load()
#     # bamFilePath <- system.file("testData/bams", package="DEScan2")
#
#     finalRegions <- countFinalRegions(regionsGRanges=regionsGR,
#                                         readsFilePath="inst/extdata/bam/onereg/",
#                                         fileType="bam",
#                                         minCarriers=3,
#                                         genomeName="mm9",
#                                         onlyStdChrs=TRUE,
#                                         ignStrandSO=TRUE, saveFlag=FALSE,
#                                         savePath="testData/regions/prova/",
#                                         verbose=TRUE)
#     saveRDS(object=finalRegions, file="inst/extdata/tests/counts/finalRegionsSE.rds")
#     filtered <- NOISeq::filtered.data(dataset=finalRegions, factor=c(rep(1, 8)),
#                                     norm=FALSE, method=3, cpm=1500)
#     rownames(filtered)
# })
