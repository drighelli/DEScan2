library("DEScan")
context("Peak Finding")
#
# test_that("findPeaks is consistent with both bed and bam inputs", {
#   test_bam <- system.file("extdata", "test.bam", package = "DEScan")
#
#   bam_peaks <- findPeaks(test_bam, chr = 1, fraglen = 200, rlen = 100, min_win = 1,
#                      max_win = 100, blocksize = 10000, zthresh = 5,
#                      min_count = 0.1, grid = NA, filetype = "bam", save = F)
#
#   expect_equal_to_reference(bam_peaks, file = "bam_peaks.rds")
#
#   test_bed <- system.file("extdata", "test.bed", package = "DEScan")
#
#   bed_peaks <- findPeaks(test_bed, chr = 1, fraglen = 200, rlen = 100, min_win = 1,
#                      max_win = 100, blocksize = 10000, zthresh = 5,
#                      min_count = 0.1, grid = NA, filetype = "bed", save = F)
#
#   expect_equal_to_reference(bed_peaks, file = "bed_peaks.rds")
#
#   expect_equal(bam_peaks, bed_peaks)
# })


test_that("Test if the new findPeaks is consistent with the older one", {
    bam.path <- system.file("extdata/Bam/chr19", package = "DEScan")
    bam.files <- list.files(bam.path, full.names = TRUE)
    testBam <- bam.files[1]
    # print(testBam)

    bamPeaksOld <- find_Peaks_Old(files=testBam, chr = 19, fraglen = 200,
                                rlen = 100, min_win = 1, max_win = 100,
                                blocksize = 10000, zthresh = 5,
                                min_count = 0.1,
                                filetype = "bam", save = FALSE)

    bamPeaksNew <- DEScan::findPeaks(files=testBam, chr = 19, fragmentLength = 200, readLength = 100, minWin = 1, maxWin = 100, blocksize = 10000, zthresh = 5, minCount = 0.1, filetype = "bam", save = FALSE)

    expect_equal(bamPeaksOld, bamPeaksNew)
})
