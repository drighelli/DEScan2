#' Align peaks to form common regions then filter regions for presence in
#' multiple replicates
#'
#' @param peak_path Path to peak files.
#' @param zthresh Integer indicating minimum z-score considered significant
#' @param min_carriers Integer indiciating the minimum number of replicates a
#'     region must be present in to be retained for testing
#' @param save_file Character indicating whether to save region files in the
#'     "bed" or "RData" format. Default is "bed".
#' @param keep_files Logical indicating whether to erase chromosome files after
#'     concatenating into a genome wide file.
#' @param chromosome Integer indicating which chromsomes to align. Defaults to
#'     all chromosomes.
#' @return Matrix containing the regions as rows with columns for genomic
#'     coordinates, z-score and number of carriers.
#' @export
#' @importFrom utils read.table write.table
finalRegions <- function(peak_path, zthresh=20, min_carriers=2,
                       chr=1:19, save_file="bed",
                       keep_files=FALSE, verbose=FALSE) {
        if (length(chr) == 1) {
            single_chromosome <- TRUE
        } else { single_chromosome <- FALSE }

        if (dir.exists("FinalRegions") == FALSE) {
                dir.create("FinalRegions")
        }

        for (i in chr) {
            chromosome <- paste0("chr", i)
            #print(paste(peak_path, chromosome, sep="/"))
            sall <- loadPeaks(paste(peak_path, chromosome, sep="/"),
                            verbose=verbose) ## a list of each file peak (for the chr)
            sall_sorted <- matrix(nrow=0, ncol=ncol(sall[[1]]) + 1) ## empty matrix ## as GRange
            for (i in 1:length(sall)) { ## for each file
                sel <- which(as.numeric(sall[[i]][, 4]) >= zthresh) ## choose peaks over the threshold
                sall_sorted <- rbind(sall_sorted, cbind(sall[[i]][sel,], rep(i, length(sel))))
            }
            ord <- order(as.numeric(sall_sorted[,2])) ## based on start range
            sall_sorted <- sall_sorted[ord,]
            common_regions <- merge_overlapping_intervals(sall_sorted, verbose=verbose) ## intervals
            names(common_regions) <- c("Start","End","AvgZ","NumCarriers")
            sel <- which(common_regions[, 4] >= min_carriers) ## seleziona solo le regioni sotto la soglia voluta di numcarriers
            final_regions <- common_regions[sel,]
            final_regions <- cbind(rep(chromosome, nrow(final_regions)), final_regions) ## create final region dataframe
            len <- final_regions[, 3] - final_regions[, 4] ##
            if (verbose) {
            cat("Regions scanned: ", nrow(final_regions), "\n")
            cat("Bases scanned (in MB): ", sum(len)/1000000, "\n")
            }

            if (save_file == "RData") {
                save(final_regions, chromosome,
                    file=paste0("FinalRegions/FinalRegions_",
                                chromosome, ".RData"))
            }
            if (save_file == "bed") {
                fname <- paste0("FinalRegions/FinalRegions_", chromosome, ".bed")
                utils::write.table(x=final_regions,
                                   file=fname, sep="\t", col.names=FALSE,
                                   row.names=FALSE, quote=FALSE)
            }
        }
        if (single_chromosome == FALSE) {
            final_regions <- catRegions(path="FinalRegions", type=save_file,
                                      keep_files=keep_files)
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
        avg_vec <- as.numeric(s[,4]) ## z-score
        count_vec <- as.numeric(s[,5]) ## samples
        s <- cbind(as.numeric(s[, 2]) - 200, #200 ## start - 200 (extend start)
                                as.numeric(s[, 3]) + 200) #mess end + 200 (extend end)
        if (is.null(avg_vec)) {
            avg_vec <- rep(0, nrow(s))
        }
        if (is.null(count_vec)) {
            count_vec <- c(1:nrow(s))
        }
        s2 <- matrix(nrow=0, ncol=4)
        n <- nrow(s)
        i <- 1
        j <- 2
        while (i <= n) { #all over the dimension of s
            #if (i %% 100 == 0 & verbose) cat(i, "\n")
            interval <- s[i, ] ## peak-i
            avg_over <- avg_vec[i] ## zscore
            count_over <- count_vec[i] ## sample
            while (j <= n && s[j, 1] <= interval[2]) { ## fin quando lo start del picco-j è incluso nel picco-i
                interval <- c(interval[1], max(interval[2], s[j, 2]))  ## picco con i è lo start del picco-i e la end il massimo tra i due (i-j) end
                avg_over <- c(avg_over, avg_vec[j]) ## accoda zscore
                count_over <- c(count_over, count_vec[j]) ## accoda campione
                j <- j + 1
                # print(j)
            }
            avg <- mean(avg_over, na.rm=TRUE) ## calcola media degli zscore
            count <- unique(count_over) ## prendi ogni campione una sola volta
            s2 <- rbind(s2, c(interval, avg, length(count))) ## accoda alla nuova matrice compreso il numero di samples per trovare l'intervallo
            i <- j
            j <- i + 1
            # print(i)
            # print(j)
        }
        return(as.data.frame(s2))
}

#' load peak files.
#'
#' @param peakdirname
#' @return sall.
#' @keywords internal
loadPeaks <- function(peakdirname, verbose=verbose) {
        all.files <- list.files(peakdirname, pattern="Peaks")
        if (verbose) {
            cat("Found", length(all.files),"peak files.\n")
        }
        peaks_all <- vector("list", length(all.files))
        for (i in 1:length(all.files)) {
            load(paste0(peakdirname, "/", all.files[i]))
            if (ncol(peaks) == 3) {
                chr <- strsplit(peakdirname, split="/")[[1]][2]
                peaks <- cbind(rep(chr, nrow(peaks)), peaks)
            }
            peaks_all[[i]] <- peaks
            if (verbose) {
                cat("File: ", all.files[i], " number of regions:", nrow(peaks), "\n") ## use message instead
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
catRegions <- function(path="FinalRegions", type="RData",
                                             keep_files=TRUE) {
        files <- list.files(path, recursive=TRUE, full.names=TRUE,
                                                pattern=type)
        final_regions_allchr <- NULL
        for (f in files) {
            if (type == "RData") {
                load(f)
            }
            if (type == "bed") {
                final_regions <- utils::read.table(f, sep="\t", header=FALSE)
                colnames(final_regions) <- c("chr","start","end","AvgZ",
                                           "NumCarriers")
            }
            final_regions_allchr <- rbind(final_regions_allchr, final_regions)
        }
        utils::write.table(final_regions_allchr,
                           file=paste0(path, "/FinalRegions_allChr.bed"),
                           sep="\t", col.names=TRUE, row.names=FALSE,
                           quote=FALSE)
        if (keep_files == FALSE) {
            file.remove(files)
            unlink(list.dirs(path)[-1], recursive=TRUE)
        }
        return(final_regions_allchr)
}

giveUniqueNamesToPeaksOverSamples <- function(samplePeaksGRangelist)
{
    stopifnot(is(samplePeaksGRangelist, "GRangesList"))
    ## this is just for naming the peaks
    ## total number of samples
    ns <- length(samplePeaksGRangelist)
    ## total number of decimals for the samples
    ncs <- nchar(as.character(ns))
    ## total number of decimals for the peaks taking the max number of peaks
    ncp <- nchar(as.character(max(unlist(
                                         lapply(samplePeaksGRangelist, length)
                                         )
                                    )))
    sFormat <- paste0("s%0", ncs,"d")
    format <- paste0("s%0", ncs,"d_p%0", ncp, "d")
    listNames <- character()
    samplePeaksGRangelista <- lapply(
        seq_along(samplePeaksGRangelist),
        function(x, i)
        {
            peakNames <- sprintf(format, i, 1:length(x[[i]]))
            # print(peakNames)
            names(x[[i]]) <- peakNames
            return(x[[i]])
        },
        x=samplePeaksGRangelist
    )

    names(samplePeaksGRangelista) <- sprintf(sFormat,
                                             seq_along(samplePeaksGRangelist))
    ####
    return(samplePeaksGRangelista)
}
#
# initMergedPeaks <- function(mergedGRanges)
# {
#     stopifnot(is(mergedGRanges, "GRanges"))
#     ncp <- length(mergedGRanges@ranges)
#     format <- paste0("s%0", ncs,"d_p%0", ncp, "d")
# }

findOverlapsOverSamples <- function(samplePeaksGRangelist,
                                    minOverlap=0L,
                                    maxGap=200)
{
    stopifnot(is(samplePeaksGRangelist, "GRangesList"))

    namedSamplePeaksGRL <- giveUniqueNamesToPeaksOverSamples(samplePeaksGRangelist)

    namedSamplePeaksGRL <- lapply(namedSamplePeaksGRL, function(x)
    {
        mcols(x)[["n-peaks"]] <-  1
        mcols(x)[["k-carriers"]] <-  1
        return(x)
    })


    message("Computing overlapping reagions over samples...")
    startTime <- Sys.time()
    for(i in 2:length(namedSamplePeaksGRL))
    {
        if( i == 2 ) {
            gri <- namedSamplePeaksGRL[[1]]
            foundedPeaks <- NULL
        } else {
            gri <- foundedPeaks
        }

        grj <- namedSamplePeaksGRL[[i]]

        grij <- ChIPpeakAnno::findOverlapsOfPeaks(gri,
                                                  grj,
                                                  minoverlap=minOverlap,
                                                  maxgap=maxGap,
                                                  connectedPeaks="merge"
                                                  )

        mmpeaks <- grij$peaksInMergedPeaks
        ## cleaning peaks names
        mrgPks <- grij$mergedPeaks
        mrgPksNms <- as.list(mrgPks$peakNames)
        stTime <- Sys.time()
        newcols <- lapply(mrgPksNms, function(l)
        {
            idx <- which(names(mmpeaks) %in% l)
            scores <- as.numeric(mmpeaks$`z-score`[idx])
            nPeaks <- as.numeric(mmpeaks$`n-peaks`[idx])
            kCarr <- as.numeric(mmpeaks$`k-carriers`[idx])
            np = sum(nPeaks) ## total number of peaks found
            ## it's necessary to rescale the score on the basis of the peaks
            ## found from previous computations
            mmzp <- sum(scores*nPeaks)/np
            ## the carriers are just the number of samples
            k <- max(kCarr)+1
            as.data.frame(cbind(mmzp, np, k))
        })
        endTime <- Sys.time()
        print((endTime-stTime))
        newcols1 <- data.table::rbindlist(newcols)
        mcols(mrgPks) <-  S4Vectors::DataFrame(newcols1)
        colnames(mcols(mrgPks)) <- c("z-score", "n-peaks", "k-carriers")
        ## peaks uniques
        unqPks <- grij$uniquePeaks
        ## putting together all the peaks
        foundedPeaks <- unlist(GRangesList(unqPks, mrgPks))

        names(foundedPeaks@ranges) <- NULL
    }
    endingTime <- Sys.time()
    print((endingTime - startTime))
    message("...done!")
    return(foundedPeaks)
}

convertSallToGrl <- function(sall)
{
    lgr <- list()
    for(sample in sall)
    {
        gr <- GRanges(seqnames=sample[,1],
                      ranges=IRanges(start=as.numeric(sample[,2]),
                                     end=as.numeric(sample[,3])
                      )
        )
        mcols(gr)[["z-score"]] <- sample[,4]
        lgr <- c(lgr, gr)
    }

    grl <- GRangesList(lgr)
    return(grl)
}
