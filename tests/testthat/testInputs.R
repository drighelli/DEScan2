# library("DEScan")
context("Input reading")
test_that("Check if bed and bam file produces the same object", {
    bed.path <- system.file("extdata/tests/inputdata/bed", package="DEScan2")
    bed.files <- list.files(path=bed.path, full.names=TRUE, pattern=".bed")
    bam.path <- system.file("extdata/tests/inputdata/bam", package="DEScan2")
    bam.files <- list.files(bam.path, full.names=TRUE, pattern=".bam$")
    bam.file <- bam.files[1]
    bed.file <- bed.files[1]
    bamRange <- constructBedRanges(filename=bam.file, filetype="bam",
                                   genomeName="mm9", onlyStdChrs=TRUE,
                                   verbose=FALSE)
    bedRange <- constructBedRanges(filename=bed.file, filetype="bed.zip",
                                   genomeName="mm9", onlyStdChrs=TRUE,
                                   verbose=FALSE)
    expect_identical(bamRange, bedRange)
})

