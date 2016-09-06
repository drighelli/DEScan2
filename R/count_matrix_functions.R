countFinalRegions <- function(region_file, bam_files, min_carriers = 2) {
  final_regions <- read.table(region_file, sep = "\t", header = T)
  final_regions <- final_regions[final_regions$NumCarriers >= min_carriers, ]
  region_anno <- cbind("GeneID" = rownames(final_regions), final_regions[,1:3],
                      "Strand" = rep("*", dim(final_regions)[1]))
  colnames(region_anno)[2:4] <- c("Chr", "Start", "End")
  rnames <- paste(region_anno$Chr,
                  paste(region_anno$Start, region_anno$End, sep = "-"), sep = ":")
  rownames(region_anno) <- rnames
  count <- Rsubread::featureCounts(bam_files, annot.ext = region_anno, nthreads = 8)
  countmat <- count$counts
  rownames(countmat) <- rnames
  countmat <- countmat[, sort(colnames(countmat))]
  write.table(countmat, file = "countmatrix.txt", sep = "\t", row.names = T,
              col.names = T, quote = F)
  countmat
}
