# library("DEScan")
context("Input reading")
test_that("Loading same data bed and bam file produces the same object", {
    bed.path <- system.file("extdata/Bed/", package = "DEScan")
    bed.files <- list.files(bed.path, full.names = TRUE)
    bam.path <- system.file("extdata/Bam/chr19", package = "DEScan")
    bam.files <- list.files(bam.path, full.names = TRUE)
    bam.file <- bam.files[1]
    bed.file <- bed.files[1]
    bamRange <- readBamAsBed(file=bam.file)
    bedRange <- readBedFile(filename=bed.file)
    expect_identical(bamRange@ranges, bedRange@ranges)
})

