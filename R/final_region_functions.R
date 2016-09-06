finalRegions <- function(sall, zthresh = 30, min_carriers = 2,
                         savefile = "RData") {
  sall_sorted <- matrix(nrow = 0, ncol = ncol(sall[[1]]) + 1)
  for (i in 1:length(sall)) {
    sel <- which(as.numeric(sall[[i]][, 4]) >= zthresh)
    sall_sorted <- rbind(sall_sorted, cbind(sall[[i]][sel,], rep(i, length(sel))))
  }
  chr <- unique(sall_sorted[,1])
  ord <- order(sall_sorted[, 1], sall_sorted[,2])
  sall_sorted <- sall_sorted[ord,]

  common_regions <- merge_overlapping_intervals(sall_sorted)
  names(common_regions) <- c("Start","End","AvgZ","NumCarriers")

  sel <- which(common_regions[, 4] >= min_carriers)
  final_regions <- common_regions[sel,]
  final_regions <- cbind(rep(chr, nrow(final_regions)), final_regions)
  len <- final_regions[, 3] - final_regions[, 4]
  cat("Regions scanned: ", nrow(final_regions), "\n")
  cat("Bases scanned (in MB): ", sum(len)/1000000, "\n")

  if (savefile == "RData") {
    save(final_regions, chr,
         file = paste0("FinalRegions/FinalRegions_", chr, ".RData"))
  }
  if (savefile == "bed") {
    write.table(x = final_regions,
                file = paste0("FinalRegions/FinalRegions_", chr, ".bed"),
                sep = "\t", col.names = F, row.names = F, quote = F)
  }
  final_regions
}

loadPeaks <- function(peakdirname) {
  all.files <- list.files(peakdirname, pattern = "Peaks")
  cat("Found", length(all.files),"peak files.\n")
  sall <- vector("list", length(all.files))
  for (i in 1:length(all.files)) {
    load(paste0(peakdirname, "/", all.files[i]))
    #s <- S
    if (ncol(s) == 3) {
      chr <- strsplit(peakdirname, split = "/")[[1]][2]
      s <- cbind(rep(chr, nrow(s)), s)
    }
    sall[[i]] <- s
    cat("File: ", all.files[i], " number of regions:", nrow(s), "\n")
  }
  sall
}

catPeaks <- function(path = "Peaks", keep_files = T) {

  files <- list.files(path, recursive = T, full.names = T)

  basenames <- unique(unlist(lapply(basename(files),
                       function(x) {strsplit(x, split = "_")[[1]][2]})))
  for (basename in basenames) {
    s_allchr <- NULL
    for (f in files[grep(basename, files)]) {
      load(f)
      s_allchr <- rbind(s_allchr, s)
    }
    s <- s_allchr
    save(s, fraglen, rlen, zthresh, min_win, max_win,
         file = paste0("Peaks/", basename, "_Peaks_allChr.RData"))
  }
  if (keep_files == F) {
    file.remove(files)
    unlink(list.dirs(path)[-1], recursive = T)
  }
}

catRegions <- function(path = "FinalRegions", type = "RData", keep_files = T) {
  files <- list.files(path, recursive = T, full.names = T, pattern = type)
  final_regions_allchr <- NULL
  for (f in files) {
    if (type == "RData") {
      load(f)
    }
    if (type == "bed") {
      final_regions <- read.table(f, sep = "\t", header = F)
      colnames(final_regions) <- c("chr","start","end","AvgZ","NumCarriers")
    }
    final_regions_allchr <- rbind(final_regions_allchr, final_regions)
  }

  write.table(final_regions_allchr,
              file = paste0(path, "/FinalRegions_allChr.bed"),
              sep = "\t", col.names = F, row.names = F, quote = F)
  if (keep_files == F) {
    file.remove(files)
    unlink(list.dirs(path)[-1], recursive = T)
  }
}

merge_overlapping_intervals <- function(s) {

  avg_vec <- as.numeric(s[,4])
  count_vec <- as.numeric(s[,5])
  s <- cbind(as.numeric(s[, 2]) - 200,
             as.numeric(s[, 3]) + 200)

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
    if (i %% 100 == 0) cat(i, "\n")
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
  as.data.frame(s2)
}

