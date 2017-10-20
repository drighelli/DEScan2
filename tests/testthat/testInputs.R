test_that("Loading same data bed and bam file produces the same object", {
    bed.path <- system.file("extdata/Bed/chr19", package = "DEScan")
    bed.files <- list.files(bed.path, full.names = TRUE)
    bam.path <- system.file("extdata/Bam/chr19", package = "DEScan")
    bam.files <- list.files(bam.path, full.names = TRUE)

    bam.file <- bam.files[1]

    bed.file <- bed.files[1]

    expect_identical(readBamAsBed(file=bam.file, chr="chr9"), readBedFile(filename=bed.file, chr="chr9"))
})

