#' Align peaks to form common regions then filter regions for presence in
#' multiple replicates
#'
#' @param peak_path Path to peak files.
#' @param zthresh Integer indicating minimum z-score considered significant
#' @param min_carriers Integer indiciating the minimum number of replicates a
#'   region must be present in to be retained for testing
#' @param save_file Character indicating whether to save region files in the
#'   "bed" or "RData" format. Default is "bed".
#' @param keep_files Logical indicating whether to erase chromosome files after
#'   concatenating into a genome wide file.
#' @param chromosome Integer indicating which chromsomes to align. Defaults to
#'   all chromosomes.
#' @return Matrix containing the regions as rows with columns for genomic
#'   coordinates, z-score and number of carriers.
#' @export
#' @importFrom utils read.table write.table
finalRegions <- function(peak_path, zthresh = 30, min_carriers = 2,
                         chr = 1:19, save_file = "bed",
                         keep_files = TRUE, verbose = FALSE) {
    if (length(chr) == 1) {
      single_chromosome <- TRUE
    } else { single_chromosome <- FALSE }

    for (i in chr) {
      chromosome <- paste0("chr", i)
      #print(paste(peak_path, chromosome, sep = "/"))
      sall <- loadPeaks(paste(peak_path, chromosome, sep = "/"), verbose = verbose)
      sall_sorted <- matrix(nrow = 0, ncol = ncol(sall[[1]]) + 1)
      for (i in 1:length(sall)) {
        sel <- which(as.numeric(sall[[i]][, 4]) >= zthresh)
        sall_sorted <- rbind(sall_sorted,
                             cbind(sall[[i]][sel,], rep(i, length(sel))))
      }
      ord <- order(as.numeric(sall_sorted[,2]))
      sall_sorted <- sall_sorted[ord,]

      common_regions <- merge_overlapping_intervals(sall_sorted, verbose = verbose)
      names(common_regions) <- c("Start","End","AvgZ","NumCarriers")

      sel <- which(common_regions[, 4] >= min_carriers)
      final_regions <- common_regions[sel,]
      final_regions <- cbind(rep(chromosome, nrow(final_regions)),
                             final_regions)
      len <- final_regions[, 3] - final_regions[, 4]
      if (verbose) {
      cat("Regions scanned: ", nrow(final_regions), "\n")
      cat("Bases scanned (in MB): ", sum(len)/1000000, "\n")
      }

      if (save_file == "RData") {
        save(final_regions, chromosome,
             file = paste0("FinalRegions/FinalRegions_", chromosome, ".RData"))
      }
      if (save_file == "bed") {
        fname <- paste0("FinalRegions/FinalRegions_", chromosome, ".bed")
        utils::write.table(x = final_regions,
                           file = fname, sep = "\t", col.names = FALSE,
                           row.names = FALSE, quote = FALSE)
      }
    }
    if (single_chromosome == FALSE) {
      final_regions <- catRegions(path = "FinalRegions", type = save_file,
                                  keep_files = keep_files)
    }
    colnames(final_regions) <- c("Chr", "Start", "End", "AvgZ", "NumCarriers")
    return(final_regions)
}

#' merge overlapping peaks across replicates.
#'
#' @param s.
#' @return s2.
#' @keywords internal
merge_overlapping_intervals <- function(s, verbose) {
    avg_vec <- as.numeric(s[,4])
    count_vec <- as.numeric(s[,5])
    s <- cbind(as.numeric(s[, 2]) - 200, #200
                as.numeric(s[, 3]) + 200) #mess

    if (is.null(avg_vec)) {
      avg_vec <- rep(0, nrow(s))
    }
    if (is.null(count_vec)) {
      count_vec <- c(1:nrow(s))
    }

    s2 <- matrix(nrow = 0, ncol = 4)
    n <- nrow(s)
    i <- 1
    j <- 2
    while (i <= n) {
      #if (i %% 100 == 0 & verbose) cat(i, "\n")
      interval <- s[i, ]
      avg_over <- avg_vec[i]
      count_over <- count_vec[i]

      while (j <= n && s[j, 1] <= interval[2]) {
        interval <- c(interval[1], max(interval[2], s[j, 2]))
        avg_over <- c(avg_over, avg_vec[j])
        count_over <- c(count_over, count_vec[j])
        j <- j + 1
      }
      avg <- mean(avg_over, na.rm = TRUE)
      count <- unique(count_over)
      s2 <- rbind(s2, c(interval, avg, length(count)))
      i <- j
      j <- i + 1
    }
    return(as.data.frame(s2))
}

#' load peak files.
#'
#' @param peakdirname
#' @return sall.
#' @keywords internal
loadPeaks <- function(peakdirname, verbose = verbose) {
    all.files <- list.files(peakdirname, pattern = "Peaks")
    if (verbose) {
      cat("Found", length(all.files),"peak files.\n")
    }
    peaks_all <- vector("list", length(all.files))
    for (i in 1:length(all.files)) {
      load(paste0(peakdirname, "/", all.files[i]))
      if (ncol(peaks) == 3) {
        chr <- strsplit(peakdirname, split = "/")[[1]][2]
        peaks <- cbind(rep(chr, nrow(peaks)), peaks)
      }
      peaks_all[[i]] <- peaks
      if (verbose) {
        cat("File: ", all.files[i], " number of regions:", nrow(peaks), "\n")
      }
    }
    peaks_all
}

#' concatenate region files from various chromosomes into a genome
#' wide region file.
#'
#' @param path.
#' @param type.
#' @param keep_files
#' @return c.
#' @keywords internal
#' @importFrom utils read.table write.table
catRegions <- function(path = "FinalRegions", type = "RData",
                       keep_files = TRUE) {
    files <- list.files(path, recursive = TRUE, full.names = TRUE,
                        pattern = type)
    final_regions_allchr <- NULL
    for (f in files) {
      if (type == "RData") {
        load(f)
      }
      if (type == "bed") {
        final_regions <- utils::read.table(f, sep = "\t", header = FALSE)
        colnames(final_regions) <- c("chr","start","end","AvgZ","NumCarriers")
      }
      final_regions_allchr <- rbind(final_regions_allchr, final_regions)
    }
    utils::write.table(final_regions_allchr,
                file = paste0(path, "/FinalRegions_allChr.bed"),
                sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
    if (keep_files == FALSE) {
      file.remove(files)
      unlink(list.dirs(path)[-1], recursive = TRUE)
    }
    return(final_regions_allchr)
}
