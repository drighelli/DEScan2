#' finalRegions
#' @description Align peaks to form common regions then filter regions for
#' presence in multiple replicates taking in input a GRangesList where each
#' element is a sample of called peaks.
#'
#' @param peakSamplesGRangesList GRangesList where each element is a sample
#' of called peaks. A z-score mcols values is needed for each GRanges.
#' (tipically returned by findPeaks function).
#' @param zThreshold a threshold for the z score.
#' @param minCarriers a threshold of minimum samples (carriers) for overlapped
#'                    regions.
#' @param saveFlag a flag for saving results in a bed file.
#' @param outputName the directory name to store the bed file.
#'
#' @return a GRanges of overlapping peaks selected with z-score,
#' n-peaks, k-carriers as mcols.
#' @export
#' @importFrom GenomicRanges GRangesList
# @examples TBW
finalRegions <- function(peakSamplesGRangesList, zThreshold=20, minCarriers=2,
                                    saveFlag=TRUE, outputName="overlappedPeaks",
                                    verbose=FALSE)
{
    stopifnot(is(peakSamplesGRangesList, "GRangesList"))

    if(verbose) message("computing final regions on ",
                        length(peakSamplesGRangesList), " samples...")
    zedPeaksSamplesGRList <- GenomicRanges::GRangesList(
                    lapply(peakSamplesGRangesList, function(sample)
                    {
                        return(sample[which(sample$`z-score` >= zThreshold),])
                    }))

    zedPeaksChrsGRList <- fromSamplesToChrsGRangesList(zedPeaksSamplesGRList)

    #### to parallelize over chrs
    overlappedPeaksGRList <- GenomicRanges::GRangesList(
        lapply(zedPeaksChrsGRList, function(chrSampleGRList) {
            return(findOverlapsOverSamples(chrSampleGRList))
    }))
    overlappedPeaksGR <- unlist(overlappedPeaksGRList)
    idxK <- which(overlappedPeaksGR$`k-carriers` >= minCarriers)
    overlMinKPeaksGR <- overlappedPeaksGR[idxK,]

    if(saveFlag) {
        ##################### check
        filename <- paste0(file, "_zt", zThreshold, "_minK", minCarriers)
        saveGRangesAsBed(GRanges=overlMinKPeaksGR, filepath=outputName,
                        filename=filename, verbose=verbose)
    }
    if(verbose) message("done\n")
    return(overlMinKPeaksGR)
}

#' giveUniqueNamesToPeaksOverSamples
#' @description given a GRangesList of samples assigns unique names to the peaks
#' of each sample.
#' @param samplePeaksGRangelist a GRangeList of peaks, one GRanges for each
#' sample.
#'
#' @return a GRangesList of samples within renamed peaks for each element.
#' @keywords internal
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
    return(samplePeaksGRangelista)
}

#' initMergedPeaksNames
#' @description given a GRanges of merged peaks assigns them new names.
#' @param mergedGRanges A GRanges object.
#' (Tipically Generated in findOverlapsOverSamples function )
#'
#' @return a granges of renamed peaks.
#' @keywords internal
initMergedPeaksNames <- function(mergedGRanges)
{
    stopifnot(is(mergedGRanges, "GRanges"))
    ncp <- nchar(length(mergedGRanges@ranges))
    nk <- nchar(max(mergedGRanges$`k-carriers`))
    np <- nchar(max(mergedGRanges$`n-peaks`))
    format <-  paste0("p%0", ncp, "d_np%0", np, "d_k%0", nk,"d")

    peakNames <- sprintf(format, 1:length(mergedGRanges@ranges),
                        mergedGRanges$`n-peaks`, mergedGRanges$`k-carriers`)
    names(mergedGRanges) <- peakNames
    return(mergedGRanges)
}

#' findOverlapsOverSamples
#' @description given in input a GRangeList where each element is a sample
#' computes the coverage extending a both direction window of prefixed length.
#'
#' @param samplePeaksGRangelist given a granges list of samples finds
#' the overlapping regions between them.
#' @param extendRegions the number of bases to extend each region at its start
#' and end.
#' @param minOverlap the minimum overlap each peak needs to have.
#' (see ChipPeakAnno::findOverlapsOfPeaks)
#' @param maxGap the maximum gap admissible between the peaks.
#' (see ChipPeakAnno::findOverlapsOfPeaks)
#'
#' @return a GRanges of peaks overlapped and unique between samples.
#' @export
#'
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom ChIPpeakAnno findOverlapsOfPeaks
#' @importFrom data.table rbindlist
#' @importFrom GenomicRanges GRangesList
#'
# @examples TBW
findOverlapsOverSamples <- function(samplePeaksGRangelist,
                                    extendRegions=200,
                                    minOverlap=0L,
                                    maxGap=-1L
                                    )
{
    stopifnot(is(samplePeaksGRangelist, "GRangesList"))

    namedSamplePeaksGRL <- giveUniqueNamesToPeaksOverSamples(
                                                        samplePeaksGRangelist)

    namedSamplePeaksGRL <- lapply(namedSamplePeaksGRL, function(x)
    {
        x <- x[as.numeric(S4Vectors::mcols(x)[["z-score"]]) >= 20]
        S4Vectors::mcols(x)[["n-peaks"]] <-  1
        S4Vectors::mcols(x)[["k-carriers"]] <-  1
        BiocGenerics::start(x) <- BiocGenerics::start(x) - extendRegions
        idxNeg <- which( BiocGenerics::start(x) < 0 )
        if(length(idxNeg) > 0)
        {
            warning("extendRegions of ", extendRegions,
                    " is too high for region(s) start ", idxNeg,
                    " forcing these starts to 0")
            BiocGenerics::start(x)[idxNeg] <- 0
        }
        BiocGenerics::end(x) <- BiocGenerics::end(x) + extendRegions
        ## it must return just one chromosome
        idxHigh <- which( BiocGenerics::end(x) > GenomeInfoDb::seqlengths(x) )
        if(length(idxHigh) > 0) {
            warning("extendRegions of ", extendRegions,
                    " is too high for region(s) end ", idxHigh,
                    " forcing these ends to ", GenomeInfoDb::seqlengths(x))
            BiocGenerics::end(x)[idxHigh] <- GenomeInfoDb::seqlengths(x)
        }

        return(x)
    })


    message("Computing overlapping reagions over samples...")
    # startTime <- Sys.time()
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
                                                connectedPeaks="merge")

        mmpeaks <- grij$peaksInMergedPeaks
        if(length(grij$peaksInMergedPeaks) == 0) stop("No merged peaks found!")
        ## cleaning peaks names
        mrgPks <- grij$mergedPeaks
        mrgPksNms <- as.list(mrgPks$peakNames)
        # stTime <- Sys.time()
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
        # endTime <- Sys.time()
        # print((endTime-stTime))
        newcols1 <- data.table::rbindlist(newcols)
        S4Vectors::mcols(mrgPks) <-  S4Vectors::DataFrame(newcols1)
        colnames(S4Vectors::mcols(mrgPks)) <- c("z-score", "n-peaks",
                                                "k-carriers")
        ## peaks uniques
        unqPks <- grij$uniquePeaks
        ## putting together all the peaks
        foundedPeaks <- unlist(GenomicRanges::GRangesList(unqPks, mrgPks))

        foundedPeaks <- initMergedPeaksNames(foundedPeaks)
    }
    # endingTime <- Sys.time()
    # print((endingTime - startTime))
    message("...done!")
    # save(foundedPeaks, file="testData/new_files/foundedPeaks.RData")
    return(foundedPeaks)
}


.convertSallToGrl <- function(sall)
{
    lgr <- list()
    for(sample in sall)
    {
        gr <- GenomicRanges::GRanges(seqnames=sample[,1],
                        ranges=IRanges::IRanges(start=as.numeric(sample[,2]),
                                                end=as.numeric(sample[,3])
                                                )
        )
        S4Vectors::mcols(gr)[["z-score"]] <- as.numeric(sample[,4])
        lgr <- c(lgr, gr)
    }

    grl <- GenomicRanges::GRangesList(lgr)
    return(grl)
}