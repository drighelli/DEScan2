# library(DEScan)
context("Region Alignment")

test_that("finalRegions is consistent", {
  peaks <- system.file("extdata", "", package = "DEScan")

  # finalRegions <- finalRegions(peaks, zthresh = 30, min_carriers = 4, chromosome = 19, keep_files = F)

  # expect_equal_to_reference(finalRegions, file = "final_regions.rds")
})


test_that("just load data for chippeakanno", {
    peak_path <- file.path(system.file("extdata", "", package = "DEScan"), "Peaks", "chr19")

    sall <- loadPeaks(peak_path, verbose=TRUE)

    grl <- convertSallToGrl(sall)

    examp_gr12 <- findOverlapsOverSamples(samplePeaksGRangelist=grl)
    print(str(examp_gr12))
    print(examp_gr12)
    # finalRegions <- finalRegions(peaks, zthresh = 30, min_carriers = 4, chromosome = 19, keep_files = F)

    # expect_equal_to_reference(finalRegions, file = "final_regions.rds")
})
