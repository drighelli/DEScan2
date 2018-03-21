context("Region Alignment")

test_that("finalRegions is consistent", {
    peaks.path <- system.file(file.path("extdata","peaks","RData"),
                            package="DEScan2")
    peaks.file <- list.files(peaks.path, full.names=TRUE, pattern="rds$")
    peaksGRL <- readRDS(peaks.file)
    finalRegions <- finalRegions(peakSamplesGRangesList=peaksGRL,
                            zThreshold=10, minCarriers=3,
                            saveFlag=FALSE, verbose=TRUE)
    regionspath <- system.file(file.path("extdata","regions"),
                               package="DEScan2")
    regionsfile <- list.files(regionspath, full.names=TRUE, pattern="rds$")
    expect_equal_to_reference(finalRegions, file=regionsfile)
})
