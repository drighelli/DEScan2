# library(DEScan)
context("Region Alignment")

test_that("finalRegions is consistent", {
  peaks <- system.file("extdata", "", package = "DEScan")

  # finalRegions <- finalRegions(peaks, zthresh = 30, min_carriers = 4, chromosome = 19, keep_files = F)

  # expect_equal_to_reference(finalRegions, file = "final_regions.rds")
})
