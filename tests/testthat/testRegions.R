# library(DEScan)
context("Region Alignment")

test_that("finalRegions is consistent", {
  peaks <- system.file("extdata/peaks", package = "DEScan2")

  # finalRegions <- finalRegions(peaks, zthresh = 30, min_carriers = 4, chromosome = 19, keep_files = F)

  # expect_equal_to_reference(finalRegions, file = "final_regions.rds")
})
#
#
# test_that("just load data for chippeakanno", {
#     # peak_path <- file.path(system.file("extdata", "", package = "DEScan"), "Peaks", "chr19")
#
#     # sall <- loadPeaks(peak_path, verbose=TRUE)
#
#     # grl <- convertSallToGrl(sall)
#     samplePeaksGRangelist <- readRDS("testData/new_files/samplePeaksGRangelist.rds")
#     foundedRegions <- findOverlapsOverSamples(samplePeaksGRangelist=samplePeaksGRangelist)
#
#     # print(str(examp_gr12))
#     # print(examp_gr12)
#     # finalRegions <- finalRegions(peaks, zthresh = 30, min_carriers = 4, chromosome = 19, keep_files = F)
#
#     # expect_equal_to_reference(finalRegions, file = "final_regions.rds")
# })
test_that("finalRegions", {
    ## need some new output of findPeaks
})

test_that("countMatrix is consistent",
{
    regions = readRDS("testData/regions/regionsGR_12_13_17.rds")
        files = list.files("testData/bams/", pattern=".bam$")
    fr = countFinalRegions(regionsGRanges=regions,
                            readsFilePath="testData/bams",
                            fileType="bam",
                            genomeName="mm9",
                            saveFlag=FALSE)
    # .loadPeaks(peakdirname=peak.path)
    # peakSamplesGRangesList=peaksGRL;
    # zThreshold=20; minCarriers=4;
    # saveFlag=FALSE; verbose=FALSE
    #
    # zedPeaksSamplesGRList <- GenomicRanges::GRangesList(
    #     lapply(peakSamplesGRangesList, function(sample)
    #     {
    #         return(sample[which(sample$`z-score` >= zThreshold),])
    #     }))
    #
    # zedPeaksChrsGRList <- fromSamplesToChrsGRangesList(zedPeaksSamplesGRList)
    # chrSampleGRList = zedPeaksChrsGRList[[1]]
    # findOverlapsOverSamples(chrSampleGRList)
    # #### to parallelize over chrs
    # overlappedPeaksGRList <- GenomicRanges::GRangesList(
    #     lapply(zedPeaksChrsGRList, function(chrSampleGRList) {
    #         return(findOverlapsOverSamples(chrSampleGRList))
    #     }))




})

