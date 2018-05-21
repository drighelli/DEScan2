# library("DEScan")
context("Input reading")
test_that("Check if bed and bam file produces the same object", {
    bed.path <- system.file(file.path("extdata","tests","inputdata","bed"),
                            package="DEScan2")
    bed.files <- list.files(path=bed.path, full.names=TRUE, pattern=".bed")
    bam.path <- system.file(file.path("extdata","tests","inputdata","bam"),
                            package="DEScan2")
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


test_that("Check readFilesAsGRangesList with and without peaks ", {
    bed.path <- system.file(file.path("extdata","peaks", "bed"),
                            package="DEScan2")

    grl <- readFilesAsGRangesList(filePath=bed.path, fileType="bed.zip",
                        genomeName="mm10", onlyStdChrs=TRUE, arePeaks=TRUE,
                        verbose=FALSE)
    grl1 <- readRDS(system.file(file.path("extdata","tests","peaks",
                            "peaks_all_samples.rds"), package="DEScan2"))
    expect_identical(grl, grl1)

    bam.path <- system.file(file.path("extdata","tests","inputdata","bam"),
                            package="DEScan2")

    grl <- readFilesAsGRangesList(filePath=bam.path, fileType="bam",
                          genomeName="mm10", onlyStdChrs=TRUE, arePeaks=FALSE,
                          verbose=FALSE)

    grl1 <- readRDS(system.file(file.path("extdata","tests","inputdata",
                            "FCGRL.rds"), package="DEScan2"))

    expect_identical(grl1, grl)
})
