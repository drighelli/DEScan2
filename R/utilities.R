#' readBamAsBed
#' @description read a bam file into a bed like format.
#'               forcing UCSC format for chromosomes names.
#' @param file Character indicating path to bam file.
#' @return GRanges object.
#'
#' @export
#' @importFrom GenomicRanges granges
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqlevelsStyle
#' @examples
#' files <- list.files(system.file("extdata/bam", package="DEScan2"),
#'                     full.names=TRUE)
#' gr <- readBamAsBed(files[1])
readBamAsBed <- function(file)
{
    ga <- GenomicAlignments::readGAlignments(file, index=file)
    gr <- GenomicRanges::granges(x=ga)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    return(gr)
}


#' readBedFile
#' @description read a bed file into a GenomicRanges like format.
#'              forcing UCSC format for chromosomes names.
#'
#' @param filename the bed filename.
#' @param arePeaks a flag indicating if the the bed file represents peaks.
#'
#' @return GRanges object
#' @export
#' @importFrom tools file_ext
#' @importFrom utils unzip
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.bed ranges strand
#' @importFrom GenomeInfoDb seqlevelsStyle seqnames
#' @importFrom S4Vectors mcols
#' @examples
#' bedFile <- list.files(system.file("extdata/bed",package="DEScan2"),
#'                         full.names=TRUE)
#' gr <- readBedFile(bedFile)
readBedFile <- function(filename, arePeaks=FALSE)
{
    if (tools::file_ext(filename) == "zip") {
        tmp <- utils::unzip(filename, list=TRUE)$Name
        file <- base::unz(filename, tmp)
    } else {
        file <- filename
    }

    bed <- rtracklayer::import.bed(con=file)
    if(!arePeaks) {
        bed <- GenomicRanges::GRanges(seqnames=GenomeInfoDb::seqnames(bed),
                                        ranges=rtracklayer::ranges(bed),
                                        strand=rtracklayer::strand(bed))
    } else {
        cidx <- grep("name", colnames(S4Vectors::mcols(bed)))
        if(length(cidx) > 0 )
        {
            S4Vectors::mcols(bed) <- S4Vectors::mcols(bed)[,-cidx, drop=FALSE]
        }
    }

    GenomeInfoDb::seqlevelsStyle(bed) <- "UCSC"
    return(bed)
}

#' setGRGenomeInfo
#' given a genome code (i.e. "mm9","mm10","hg19","hg38") retrieve the SeqInfo of
#' that genome and assigns it to the input GRanges. Finally filters out those
#' Infos not necessary to the GRanges.
#'
#' @param GRanges a GRanges object.
#' @param genomeName a genome code (i.e. "mm9")
#' @param verbose verbose output
#'
#' @return a GRanges object with the seqinfo of the genome code
#' @export
#' @importFrom S4Vectors runValue
#' @importFrom GenomeInfoDb seqnames Seqinfo seqinfo
#' @importFrom glue glue_collapse
#' @examples
#' library("GenomicRanges")
#' gr <- GRanges(
#'         seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'         ranges=IRanges(1:10, end=10),
#'         strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'         seqlengths=c(chr1=11, chr2=12, chr3=13))
#' mm9gr <- setGRGenomeInfo(GRanges=gr, genomeName="mm9", verbose=TRUE)
setGRGenomeInfo <- function(GRanges, genomeName=NULL, verbose=FALSE)
{
    stopifnot(is(GRanges, "GRanges"))
    stopifnot(!is.null(genomeName))
    if(length(genomeName)>1) stop("Please provide just one genome code!\n")

    uniqueSeqnames <- droplevels(unique(S4Vectors::runValue(
                                    GenomeInfoDb::seqnames(GRanges))))

    if(verbose) message("Get seqlengths from genome ", genomeName)
    tryCatch({genomeInfo <- GenomeInfoDb::Seqinfo(genome=genomeName)},
        error=function(e)
        {
            stop("Unable to retrieve the genome ", genomeName, " returned: ", e)
        }
    )

    seqNamesIdx <- which(GenomeInfoDb::seqnames(genomeInfo) %in% uniqueSeqnames)
    if(length(seqNamesIdx) != 0)
    {
        sqi <- genomeInfo[GenomeInfoDb::seqnames(genomeInfo)[seqNamesIdx]]
        # sqi <- sqi[GenomeInfoDb::seqnames(sqi)[
        #                                 order(GenomeInfoDb::seqnames(sqi))],]
        GenomeInfoDb::seqlevels(GRanges) <-
                                        GenomeInfoDb::seqlevelsInUse(GRanges)
        GRanges <- GenomeInfoDb::sortSeqlevels(GRanges)
        tryCatch({GenomeInfoDb::seqinfo(GRanges) <- sqi},
                warning=function(w)
                {
                    warning(paste0("The genome ", genomeName,
                                " you chose maybe it's not ",
                                "the rightest one for these reads. ",
                                "\nPlease check the genome version and retry!",
                                "\nMoreover: ", w))
                },
                error=function(e)
                {
                    stop("The genome code ", genomeName, " returned: ", e)
                }
        )

    }
    else
    {
        stop("Cannot find the ", glue::glue_collapse(uniqueSeqnames, " "),
            " in genome ", genomeName,
            " Maybe a problem of chromosome labels")
    }
    tryCatch({GenomeInfoDb::seqnames(GRanges) <- droplevels(
                                    GenomeInfoDb::seqnames(GRanges))},
            warning=function(w)
            {
                warning(paste0("The genome ", genomeName,
                    " you chose maybe it's not ",
                    "the rightest one for these reads. ",
                    "\nPlease check the genome version and retry!",
                    "\nMoreover: ", w))
            },
            error=function(e)
            {
                stop("The genome code ", genomeName, " returned: ", e)
            }
    )
    return(GRanges)
}

#' constructBedRanges
#' @description Constructs a GRanges object from a bam/bed/bed.zip file in a
#' consistent way.
#' @param filename the complete file path of a bam?bed file.
#' @param filetype the file type bam/bed/bed.zip/narrow/broad.
#' @param genomeName the name of the genome used to map the reads (i.e. "mm9").
#' N.B. if NOT NULL the GRanges Seqinfo will be forced to genomeName Seqinfo
#' (needs Internet access, but strongly suggested!)
#' @param onlyStdChrs flag to keep only standard chromosome.
#' @param arePeaks flag indicating if the file contains peaks.
#' @param verbose flag to obtain verbose output.
#' @return a GRanges object.
#' @export
#' @importFrom GenomeInfoDb keepStandardChromosomes seqinfo Seqinfo seqnames
#' keepSeqlevels
#' @importFrom GenomicRanges sort
#' @examples
#' files <- list.files(system.file("extdata/bam/", package="DEScan2"),
#'                     pattern="bam$", full.names=TRUE)
#' bgr <- constructBedRanges(files[1], filetype="bam", genomeName="mm9",
#'                             onlyStdChrs=TRUE)
#' bgr
constructBedRanges <- function(filename,
                                filetype=c("bam", "bed", "bed.zip",
                                            "narrow", "broad"),
                                genomeName=NULL,
                                onlyStdChrs=FALSE,
                                arePeaks=FALSE,
                                verbose=FALSE)
{
    filetype <- match.arg(filetype)
    if(verbose) message("processing ", filename)

    if(filetype == "bam") {
        bedGRanges <- readBamAsBed(file=filename)
    } else {
        bedGRanges <- readBedFile(filename=filename, arePeaks=arePeaks)
    }
    if(onlyStdChrs)
    {
        bedGRanges <- GenomeInfoDb::keepStandardChromosomes(x=bedGRanges,
                                                        pruning.mode="coarse")
    }

    uniqueSeqnames <- droplevels(unique(GenomeInfoDb::seqnames(bedGRanges)))

    if( !is.null(genomeName) )
    {
        bedGRanges <- setGRGenomeInfo(GRanges=bedGRanges,
                                        genomeName=genomeName,
                                        verbose=verbose)
        return(bedGRanges)
    }

    # checking bed seqnames, useful in peak calling algorithm
    veclengths <- as.vector(
                    GenomeInfoDb::seqlengths(
                        GenomeInfoDb::seqinfo(bedGRanges)))
    vecnames <- GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(bedGRanges))

    if( (sum(is.na(veclengths)) > 0) || (length(vecnames) == 0 ))
    {
        if(!arePeaks)
        {
            stop("No seqlengths present in file ", filename,
                "\nDEScan2 needs seqlenghts to work properly.",
                "\nPlease provide a genomeName to setup the GRanges!")
        }
    } else if(length(uniqueSeqnames) <
            length(GenomeInfoDb::seqnames(GenomeInfoDb::seqinfo(bedGRanges))))
    {
        if(verbose) message("Keeping only necessary seqInfos")
        bedGRanges <- GenomeInfoDb::keepSeqlevels(bedGRanges, uniqueSeqnames)
    }

    bedGRanges <- GenomicRanges::sort(bedGRanges, ignore.strand=TRUE)
    return(bedGRanges)
}

#' readFilesAsGRangesList
#' @description Takes in input the path of bam/bed files to process and stores
#' them in a GRangesList object, named with filePath/filenames.
#' (for lazy people)
#' @param filePath the path of input files.
#' @param fileType the type of the files (bam/bed/bed.zip/narrow/broad).
#' @param genomeName the genome code to associate to the files. (recommended)
#' (i.e. "mm9", "hg17")
#' @param onlyStdChrs a flag to keep only standard chromosomes.
#' @param arePeaks a flag indicating if the files contain peaks.
#' @param verbose verbose output flag.
#'
#' @return a GRangesList object
#' @importFrom GenomicRanges GRangesList
#' @export
#'
#' @examples
#' files.path <- system.file("extdata/bam", package="DEScan2")
#' grl <- readFilesAsGRangesList(filePath=files.path, fileType="bam",
#'                                 genomeName="mm9", onlyStdChrs=TRUE,
#'                                 verbose=TRUE)
#' class(grl)
#' names(grl)
#' grl
readFilesAsGRangesList <- function(filePath, fileType=c("bam", "bed", "bed.zip",
                                                        "narrow", "broad"),
                            genomeName=NULL, onlyStdChrs=TRUE, arePeaks=TRUE,
                            verbose=TRUE)
{

    fileType <- match.arg(fileType)
    stopifnot(is.character(filePath))

    files <- list.files(filePath, pattern=fileType, full.names=TRUE)
    if(fileType == "bam") {
        idx <- grep(pattern="bai", x=files)
        if(length(idx) != 0 ) files <- files[-idx]
    }
    grl <- GenomicRanges::GRangesList(
                            lapply(files, constructBedRanges,
                                    filetype=fileType,
                                    genomeName=genomeName,
                                    onlyStdChrs=onlyStdChrs,
                                    arePeaks=arePeaks,
                                    verbose=verbose)
                            )
    names(grl) <- basename(files)
    return(grl)
}

#' saveGRangesAsBed
#' @description save a GRanges object as bed file.
#'
#' @param GRanges the GRanges object.
#' @param filepath the path to store the files.@
#' @param filename the name to give to the files.
#' @param force force overwriting.
#' @param verbose verbose output flag.
#' @importFrom rtracklayer export.bed
#' @importFrom GenomeInfoDb sortSeqlevels seqnames
#' @importFrom S4Vectors mcols
#' @importFrom BiocGenerics start end
#'
#' @return none
#' @export
#' @examples
#' library("GenomicRanges")
#' gr <- GRanges(
#'         seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'         ranges=IRanges(1:10, end=10),
#'         strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'         seqlengths=c(chr1=11, chr2=12, chr3=13))
#'
#' saveGRangesAsBed(GRanges=gr, filepath=tempdir(), filename=tempfile(),
#'                     verbose=TRUE)
saveGRangesAsBed <- function(GRanges, filepath=tempdir(), filename=tempfile(),
                            force=FALSE, verbose=FALSE)#, extraCols=NULL)
{
    stopifnot(is(GRanges, "GRanges"))

    if(!exists(filepath))
    {
        dir.create(path=filepath, showWarnings=FALSE, recursive=TRUE)
    }
    filePathName <- file.path(filepath, paste0(basename(filename), ".bed"))

    if(file.exists(filePathName))
    {
        if(!force)
        {
            stop(filePathName, " already exists!\nNot overwriting!")
        }
        else
        {
            message("overwriting ", filePathName)
        }
    }


    GRanges <- GenomeInfoDb::sortSeqlevels(GRanges)
    GRanges <- sort(GRanges)
    if(length(unique(names(GRanges))) < length(GRanges))
    {
        nn <- paste0(GenomeInfoDb::seqnames(GRanges), ":",
                        BiocGenerics::start(GRanges), "-",
                        BiocGenerics::end(GRanges))
        names(GRanges) <- nn
    }
    if(length(which(colnames(S4Vectors::mcols(GRanges)) %in% "z-score")) > 0)
        if(length(which(colnames(S4Vectors::mcols(GRanges)) %in% "score")) == 0)
        S4Vectors::mcols(GRanges)$score <- S4Vectors::mcols(GRanges)$`z-score`

    # if(!is.null(extraCols))
    # {
    #     idx <- which(extraCols %in% colnames(S4Vectors::mcols(GRanges)))
    #     if(length(idx) != 0)
    #     {
    #         extraCols <- extraCols[idx]
    #         rtracklayer::export.bed(object=GRanges, con=filePathName,
    #                                 extraCols=extraCols)
    #     }
    #     else
    #     {
    #         rtracklayer::export.bed(object=GRanges, con=filePathName)
    #     }
    #
    # }
    # else
    # {
        rtracklayer::export.bed(object=GRanges, con=filePathName)
    # }


    if(verbose) message("file ", filePathName, " written on disk!")
}


#' saveGRangesAsTsv
#' @description save a GRanges object as tsv file.
#'
#' @param GRanges the GRanges object.
#' @param filepath the path to store the files.
#' @param filename the name to give to the files.
#' @param verbose verbose output flag.
#' @param force force overwriting.
#'
#' @importFrom utils write.table
#'
#' @return none
#' @export
#' @examples
#' gr <- GRanges(
#'         seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'         ranges=IRanges(1:10, end=10),
#'         strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'         seqlengths=c(chr1=11, chr2=12, chr3=13))
#' saveGRangesAsTsv(gr, verbose=TRUE)
saveGRangesAsTsv <- function(GRanges, filepath=tempdir(), filename=tempfile(),
                            force=FALSE, verbose=FALSE)
{
    stopifnot(is(GRanges, "GRanges"))

    if(!exists(filepath))
    {
        dir.create(path=filepath, showWarnings=FALSE, recursive=TRUE)
    }
    filename <- basename(filename)
    filePathName <- file.path(filepath, paste0(filename, ".tsv"))
    if(file.exists(filePathName))
    {
        if(!force)
        {
            stop(filePathName, " already exists!\nNot overwriting!")
        }
        else
        {
            message("overwriting", filePathName)
        }
    }
    if(!is.null(names(GRanges)))
    {
        rownames <- names(GRanges)
    } else {
        rownames <- NULL
    }
    grdf <- as.data.frame(GRanges, row.names=rownames)
    utils::write.table(x=grdf, file=filePathName, quote=FALSE,
                sep="\t", row.names=TRUE, col.names=NA)
    if(verbose) message("file ", filePathName, " written on disk!")
}


#' RleListToRleMatrix
#' @description a wrapper to create a RleMatrix from a RleList object.
#'
#' @param RleList an RleList object with all elements of the same length.
#' @param dimnames the names for dimensions of RleMatrix (see DelayedArray pkg).
#'
#' @return a RleMatrix from DelayedArray package.
#'
#' @importFrom  DelayedArray RleArray
#' @export
#' @examples
#' library("DelayedArray")
#' lengths <-  c(3, 1, 2)
#' values <- c(15, 5, 20)
#' el1 <- S4Vectors::Rle(values=values, lengths=lengths)
#'
#' el2 <- S4Vectors::Rle(values=sort(values), lengths=lengths)
#'
#' rleList <- IRanges::RleList(el1, el2)
#' names(rleList) <- c("one", "two")
#' (rleMat <- RleListToRleMatrix(rleList))
RleListToRleMatrix <- function(RleList, dimnames=NULL)
{
    lengths <- unlist(lapply(RleList, length), use.names=FALSE)
    stopifnot(all.equal(lengths, rep(lengths[1], length(lengths))))
    if(!is.null(dimnames)) {
        rlem <- DelayedArray::RleArray(data=unlist(RleList, use.names=FALSE),
                                        dim=c(length(RleList[[1]]),
                                        length(RleList)),
                                        dimnames=dimnames
        )
    } else {
        rlem <- DelayedArray::RleArray(data=unlist(RleList, use.names=FALSE),
                                        dim=c(length(RleList[[1]]),
                                        length(RleList)))
    }
    return(rlem)
}


#' createGranges
#' @description a simplified wrapper function to create a GRanges object.
#'
#' @param chrSeqInfo a seqinfo object.
#' @param starts the start ranges.
#' @param widths the width of each range.
#' @param mcolname the name for the mcol attribute.
#' @param mcolvalues the values for the mcol attribute.
#'
#' @return a GRanges object.
#'
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom rtracklayer start
#' @importFrom S4Vectors mcols
#' @importFrom GenomeInfoDb seqlengths isCircular seqnames
#' @export
#' @examples
#' chrSeqInfo <- GenomeInfoDb::Seqinfo(genome="mm9")["chr1"]
#' starts=sample(seq_len(100), 10)
#' widths=starts+10;
#' mcolname <- "z-score";
#' mcolvalues <- sample(seq_len(100), 10)
#' chrGR <- createGranges(chrSeqInfo=chrSeqInfo, starts=starts, widths=widths,
#'               mcolname=mcolname, mcolvalues=mcolvalues)
createGranges <- function(chrSeqInfo, starts, widths,
                            mcolname=NULL, mcolvalues=NULL) {
    stopifnot(is(chrSeqInfo, "Seqinfo"))
    stopifnot(identical(length(starts), length(widths)))

    maxlengths <- starts+widths
    slen <- as.numeric(GenomeInfoDb::seqlengths(chrSeqInfo))
    iscirc <- as.logical(GenomeInfoDb::isCircular(chrSeqInfo))
    idxm <- which(maxlengths >= slen)
    if( (length(idxm) > 0) && (!iscirc) )
    {

        warning("GRanges object contains ", length(idxm), " out-of-bound range",
                " located on sequence ",
                GenomeInfoDb::seqnames(chrSeqInfo), ".",
                " A non-circular sequence!",
                "\nTrimming out-of-bound range to the admitted seqlength."
                )
        widths[idxm] <- slen - starts[idxm]
    }

    gr <- GenomicRanges::GRanges(seqnames=as.character(
                            GenomeInfoDb::seqnames(chrSeqInfo)),
                            ranges=IRanges::IRanges(start=starts, width=widths),
                            seqinfo=chrSeqInfo)

    if(!is.null(mcolname) )
    {
        if(!is.null(mcolvalues) &&
            (length(rtracklayer::start(gr)) == length(mcolvalues))
        )
        {
            S4Vectors::mcols(gr)[[mcolname]] <- mcolvalues
        } else {
            warning("Cannot set mcols values!",
                    " Vector length not matching Ranges")
        }
    }
    return(gr)
}


#' cutGRangesPerChromosome
#' @description  takes in input a GRanges object, producing a LIST of
#' GRanges, one for each chromosome.
#'
#' @param GRanges a GRanges object.
#'
#' @return a named list of GRanges, one for each chromosome.
#'
#' @importFrom GenomeInfoDb seqnames seqinfo seqlevels seqlevelsInUse
#' @importFrom GenomicRanges GRangesList
#' @importFrom GenomicAlignments levels
#' @importFrom S4Vectors runValue
#' @export
#'
#' @examples
#' library("GenomicRanges")
#' gr <- GRanges(
#'         seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'         ranges=IRanges(1:10, end=10),
#'         strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'         seqlengths=c(chr1=11, chr2=12, chr3=13))
#' (grchrlist <- cutGRangesPerChromosome(gr))
cutGRangesPerChromosome <- function(GRanges)
{
    stopifnot(is(GRanges, "GRanges"))

    interestedChrs <- GenomicAlignments::levels(S4Vectors::runValue(
                                            GenomeInfoDb::seqnames(GRanges)))

    GRList <- lapply(interestedChrs, function(x)
    {
        bgr <- GRanges[ GenomeInfoDb::seqnames(GRanges) == x ]
        if(length(bgr) > 0)
        {
            GenomeInfoDb::seqlevels(bgr) <- GenomeInfoDb::seqlevelsInUse(bgr)
            GenomeInfoDb::seqinfo(bgr) <- GenomeInfoDb::seqinfo(GRanges)[x]
            return(bgr)
        }
    })
    names(GRList) <- interestedChrs

    ## intentionally left commented, GRangesList reconstruct the entire seqinfo,
    ## while we want it cutted per chromosomes
    ## GRList <- GRangesList(GRList)
    return(GRList)
}

#' keepRelevantChrs
#' @description  subselect a list of GRanges created with
#' cutGRangesPerChromosome returning only the relevant chromosomes GRanges.
#' @param chrGRangesList where each element is a chromosome,
#' tipically created with cutGRangesPerChromosome.
#' @param chr a character vector of chromosomes names of the form "chr#".
#'
#' @return the input chrGRangesList with only the relevat chromosomes.
#'
#' @export
#' @examples
#' library("GenomicRanges")
#' gr1 <- GRanges(
#'             seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'             ranges=IRanges(1:10, end=10),
#'             strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'             seqlengths=c(chr1=11, chr2=12, chr3=13))
#' grlc <- cutGRangesPerChromosome(gr1)
#' (grlChr <- keepRelevantChrs(grlc, c("chr1", "chr3")))
keepRelevantChrs <- function(chrGRangesList, chr=NULL)
{
    # if(!is.null(chr) && length(grep(pattern="chr", chr)) != length(chr))
    #     stop("Insert valid chr(s), use the \"chr#\" form!")
    # stopifnot(is(chrGRangesList, "GRangesList"))

    idxs <- which(names(chrGRangesList) %in% chr)
    if(length(idxs) == 0)
        stop("Something went wrong in the chr subselection!",
                "\nPlease check the Chromosomes names!")

    chrGRangesList <- chrGRangesList[idxs]

    return(chrGRangesList)
}

#' fromSamplesToChrsGRangesList
#' @description converts a GRangesList orgnized per samples to a GRangesList
#' organized per Chromosomes where each element is a GRangesList of samples.
#'
#' @param samplesGRangesList a GRangesList of samples.
#' Tipically generaed by findPeaks function.
#'
#' @return A GRangesList of chromosomes where each element is a GRanges list
#' of samples.
#' @export
#' @importFrom GenomicRanges GRangesList
#'
#' @examples
#' library("GenomicRanges")
#' gr1 <- GRanges(
#'             seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'             ranges=IRanges(1:10, end=10),
#'             strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'             seqlengths=c(chr1=11, chr2=12, chr3=13))
#' gr2 <- GRanges(
#'             seqnames=Rle(c("chr1", "chr4", "chr1", "chr3"), c(1, 3, 2, 4)),
#'             ranges=IRanges(1:10, end=10),
#'             strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'             seqlengths=c(chr1=11, chr4=12, chr3=13))
#' sgrl <- GRangesList(gr1, gr2)
#' names(sgrl) <- c("samp1", "samp2")
#' (chrGrlSampGr <- fromSamplesToChrsGRangesList(sgrl))
fromSamplesToChrsGRangesList <- function(samplesGRangesList)
{
    stopifnot(is(samplesGRangesList, "GRangesList"))
    samplesChrList <- divideEachSampleByChromosomes(samplesGRangesList)
    sampChromsTab <- generateDFofSamplesPerChromosomes(samplesChrList)

    allChrs <- unique(unlist(strsplit(sampChromsTab$chromosomes, split=";")))

    chrsSamplesList <- lapply(allChrs, function(chr)
    {
        sampNames <- sampChromsTab$samples[grep(chr, sampChromsTab$chromosomes)]
        samplesList <- samplesChrList[sampNames]
        chrList <- GenomicRanges::GRangesList(
            lapply(samplesList, function(samp)
            {
                idx <- grep(paste0(chr,"$"), names(samp))
                return(samp[[idx]]) ## it can be only one chr
            }))
        return(chrList)
    })
    names(chrsSamplesList) <- allChrs
    return(chrsSamplesList)
}

#' divideEachSampleByChromosomes
#' @description taken in input a grangeslist of samples, generate a list of
#' samples where each element has a GRangesList each element of the GRangesList
#' represents a single chromosome.
#' @param samplesGRangesList a GRangesList of samples.
#'
#' @return list of samples where each element is a list of chromosomes and each
#' of these elements is a GRanges.
#' @export
#'
#' @examples
#' library("GenomicRanges")
#' gr1 <- GRanges(
#'             seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'             ranges=IRanges(1:10, end=10),
#'             strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'             seqlengths=c(chr1=11, chr2=12, chr3=13))
#' gr2 <- GRanges(
#'             seqnames=Rle(c("chr1", "chr4", "chr1", "chr3"), c(1, 3, 2, 4)),
#'             ranges=IRanges(1:10, end=10),
#'             strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'             seqlengths=c(chr1=11, chr4=12, chr3=13))
#' sgrl <- GRangesList(gr1, gr2)
#' names(sgrl) <- c("samp1", "samp2")
#' (sampChrGrl <- divideEachSampleByChromosomes(sgrl))
divideEachSampleByChromosomes <- function(samplesGRangesList)
{

    samplesChrList <- lapply(samplesGRangesList, function(gr)
    {
        sampleChrGRList <- cutGRangesPerChromosome(gr)
        idx <- unlist(lapply(sampleChrGRList, is.null))
        idx <- !idx
        sampleChrGRList <- sampleChrGRList[which(idx)]
        return(sampleChrGRList)
    })
    return(samplesChrList)
}

#' generateDFofSamplesPerChromosomes
#' @description generates a dataframe where each row is a sample (1st col) and
#' a string with its chromosomes separated by ";" (2nd col)
#' (useful to fromSamplesToChromosomesGRangesList function).
#' @param samplesChrGRList a GRangesList of samples each divided by chromosome.
#'
#' @return a dataframe  where each row is a sample (1st col) and
#' a string with its chromosomes separated by ";" (2nd col).
#'
#' @importFrom plyr ldply
#' @keywords internal
generateDFofSamplesPerChromosomes <- function(samplesChrGRList)
{
    ## generates a dataframe where each row is a sample (1st col)
    ## and a string with its chromosomes separated by ";" (2nd col)
    sampChromsTab <- plyr::ldply(samplesChrGRList, function(sgrl)
    {
        paste(names(sgrl), collapse=";")

    })
    colnames(sampChromsTab) <- c("samples", "chromosomes")
    return(sampChromsTab)
}
