# test_that("finalRegions", {
#     library("devtools")
#     document()
#     load_all()
#     bam.files <- list.files(path="testData/bams/head", pattern="bam",
#                             full.names = TRUE)
#     bam.files <- bam.files[-grep(pattern=".bai", x=bam.files)]
#
#     # files=bam.files; chr="chr19"; filetype="bam";
#     # fragmentLength=200;
#     # binSize=50; minWin=50; maxWin=1000; genomeName="mm9";
#     # minCompWinWidth=5000; maxCompWinWidth=10000;
#     # zthresh=5; minCount=0.1; verbose=TRUE; save=TRUE;
#     # outputFolder="testData/Peaks/new_peaks_head"
#     # onlyStdChrs=TRUE
#     bam.files <- list.files(system.file("extdata/bam", package="DEScan2"),
#                             full.names=TRUE, pattern="bam$")
#
#     peaksGRL <- findPeaks(files=bam.files, filetype="bam",
#                               genomeName="mm9",
#                               binSize=50, minWin=50, maxWin=1000,
#                               zthresh=5, minCount=0.1, sigwin=10,
#                               minCompWinWidth=5000, maxCompWinWidth=10000,
#                               save=TRUE,
#                               outputFolder="testData/Peaks/new_peaks_head",
#                               force=TRUE,
#                               onlyStdChrs=TRUE,
#                               chr=NULL,
#                               verbose=TRUE)
#
#     save(object=peaksGRL, file="testData/Peaks/new_peaks_head/peaksGRL_head_12_15_17.RData")
#     save(object=peaksGRL, file="testData/Peaks/new_peaks_head/peaksGRL_zt5_mnw50_mxw1000_bin50.RData")
#     regionsGR <- finalRegions(peakSamplesGRangesList=peaksGRL,
#                                 zThreshold=20, minCarriers=4,
#                                 saveFlag=TRUE,
#                                 outputFolder="testData/regions/prova",
#                                 verbose=TRUE)
#     # saveRDS(object=regionsGR, file=paste0("testData/regions/prova/regions_",datename ,"_zt20_minK4.RDS"))
#
#     # regionsGR <- load()
#     # bamFilePath <- system.file("testData/bams", package="DEScan2")
#
#     finalRegions <- countFinalRegions(regionsGRanges=regionsGR,
#                                         readsFilePath="inst/extdata/bam",
#                                         fileType="bam",
#                                         minCarriers=1,
#                                         genomeName="mm9",
#                                         onlyStdChrs=TRUE,
#                                         ignStrandSO=TRUE, saveFlag=TRUE,
#                                         savePath="testData/regions/prova/",
#                                         verbose=TRUE)
#     saveRDS(object=finalRegions, file="finalRegions_.rds")
#     filtered <- NOISeq::filtered.data(dataset=finalRegions, factor=c(rep(1, 8)),
#                                     norm=FALSE, method=3, cpm=1500)
#     rownames(filtered)
# })
