#' count reads falling within the final regions.
#'
#' @param regions regions of interest to count reads across.
#' @param bam_file_path path to bam files to be used for read counting.
#' @param min_carriers integer indicating minumum number of replicates a region must appear in to be utilized for count matrix construction.
#' @return Matrix containing read counts with regions as rows and samples as columns.
#' @export
#' @importFrom utils read.table write.table
#' @import Rsubread
countFinalRegions <- function(regions, bam_file_path, min_carriers = 2) {
  bam_files <- list.files(bam_file_path, full.names = TRUE)
  bam_files <- bam_files[-grep("bai", bam_files)]

  #final_regions <- read.table(region_file, sep = "\t", header = TRUE)
  regions <- regions[regions[[5]] >= min_carriers, ]
  region_anno <- cbind("GeneID" = rownames(regions), regions[,1:3],
                      "Strand" = rep("*", dim(regions)[1]))
  colnames(region_anno)[2:4] <- c("Chr", "Start", "End")
  rnames <- paste(region_anno$Chr,
                  paste(region_anno$Start, region_anno$End, sep = "-"), sep = ":")
  rownames(region_anno) <- rnames
  count <- Rsubread::featureCounts(bam_files, annot.ext = region_anno, nthreads = 8)
  countmat <- count$counts
  rownames(countmat) <- rnames
  countmat <- countmat[, sort(colnames(countmat))]
  utils::write.table(countmat, file = "countmatrix.txt", sep = "\t", row.names = TRUE,
              col.names = TRUE, quote = FALSE)
  countmat
}
