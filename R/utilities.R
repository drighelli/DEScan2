#' readBamAsBed
#' @description read a bam file into a bed like format.
#'               forcing UCSC format for chromosomes names.
#' @param file Character indicating path to bam file.
#' @return GRanges object.
#'
#' @keywords internal
#' @importFrom GenomicRanges granges
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomeInfoDb seqlevelsStyle
readBamAsBed <- function(file) {
    message("processing ", file)
    ga <- GenomicAlignments::readGAlignments(file, index=file)
    gr <- GenomicRanges::granges(x = ga)
    GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
    return(gr)
}


#' readBedFile
#' @description read a bed file into a GenomicRanges like format.
#'              forcing UCSC format for chromosomes names.
#'
#' @param filename
#'
#' @return GRanges object
#' @keywords internal
#' @importFrom tools file_ext
#' @importFrom utils unzip
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqlevelsStyle
#'
readBedFile <- function(filename) {
    message("processing ", filename)
    if (tools::file_ext(filename) == "zip") {
        tmp <- utils::unzip(filename, list=T)$Name
        file <- base::unz(filename, tmp)
    } else {
        file <- filename
    }

    bed <- rtracklayer::import.bed(con = file)
    bed <- GenomicRanges::GRanges(seqnames=bed@seqnames,
                                  ranges=bed@ranges,
                                  strand=bed@strand)
    GenomeInfoDb::seqlevelsStyle(bed) <- "UCSC"
    return(bed)
}


#' constructBedRanges
#' @description Constructs a GRanges object from a bam/bed file in a consistent way
#'
#' @param filename the complete file path of a bam?bed file
#' @param filetype the file type bam/bed
#' @param genomeName the name of the genome used to map the reads (i.e. mm9)
#'                   N.B. if NOT NULL the GRanges Seqinfo will be forced to
#'                   genomeName Seqinfo
#' @param onlySrdChrs a flag to keep only standard chromosomes from the input
#'                    files.
#'
#' @return a GRanges object
#' @importFrom GenomeInfoDb keepStandardChromosomes seqinfo Seqinfo
#' @importFrom glue collapse
constructBedRanges <- function(filename,
                               filetype=c("bam", "bed"),
                               genomeName=NULL, onlyStdChrs=FALSE)
{
    filetype <- match.arg(filetype)

    if(filetype == "bam") {
        bedGRanges <- readBamAsBed(file=filename)
    } else {
        bedGRanges <- readBedFile(filename=filename)
    }
    if(onlyStdChrs)
    {
        bedGRanges <- GenomeInfoDb::keepStandardChromosomes(x=bedGRanges,
                                                         pruning.mode="coarse")
    }
    uniqueSeqnames <- droplevels(unique(bedGRanges@seqnames))

    if( !is.null(genomeName) )
    {

        message("Get seqlengths from genome ", genomeName)
        genomeInfo <- GenomeInfoDb::Seqinfo(genome=genomeName)
        seqNamesIdx <- which(genomeInfo@seqnames %in% uniqueSeqnames)
        if(length(seqNamesIdx) != 0)
        {
            bedGRanges@seqinfo <- genomeInfo[genomeInfo@seqnames[seqNamesIdx]]
        }
        else
        {
            stop("Cannot find the ", glue::collapse(uniqueSeqnames, " "),
                 " in genome ", genomeName, "
                 Maybe a problem of chromosome labels")
        }
        return(bedGRanges)
    }

    # checking bed seqnames, useful in peak calling algorithm
    if( (sum(is.na(GenomeInfoDb::seqinfo(bedGRanges)@seqlengths)) > 0) ||
        (length(GenomeInfoDb::seqinfo(bedGRanges)@seqnames) == 0 ))
    {
        stop("No seqlengths present in file ", filename,
             "\nPlease provide the correct genomeName to setup the GRanges!")
    }
    else if(length(uniqueSeqnames) < length(GenomeInfoDb::seqinfo(bedGRanges)@seqnames))
    {
        message("Keeping only necessary seqInfos")
        bedGRanges@seqinfo <- bedGRanges@seqinfo[as.character(uniqueSeqnames)]
        bedGRanges@seqnames <- droplevels(bedGRanges@seqnames)
    }

    return(bedGRanges)
}



#' saveGRangesAsBed
#' @description save a GRanges object as bed file and RData file
#'
#' @param GRanges the GRanges object
#' @param filepath the path to store the files
#' @param filename the name to give to the files
#'
#' @importFrom rtracklayer export.bed
#' @return none
#' @examples TBW
saveGRangesAsBed <- function(GRanges, filepath, filename)
{
    stopifnot(is(GRanges, "GRanges"))
    ## add some parameters
    dir.create(path=filepath, showWarnings=FALSE, recursive=TRUE)
    filePathName <- file.path(filepath, paste0(filename, "_peaks"))
    if(file.exists(filePathName)) {stop(filePathName, " already exists!\n"
                                        , "Not overwriting!")}
    GRanges <- sortSeqlevels(GRanges)
    GRanges <- sort(GRanges)
    rtracklayer::export.bed(object=GRanges, con=paste0(filePathName, ".bed"))
    # save(GRanges, file=paste0(filePathName, ".RData"))
}


#' RleListToRleMatrix
#' @description a wrapper to create a RleMatrix from a RleList object
#'
#' @param RleList an RleList object
#' @param dimnames the names for dimensions of RleMatrix (see DelayedArray pkg)
#'
#' @return a RleMatrix from DelayedArray package
#' @importFrom  DelayedArray RleArray
#' @export
#' @examples TBW
RleListToRleMatrix <- function(RleList, dimnames=NULL)
{
    if(!is.null(dimnames)) {
        rlem <- DelayedArray::RleArray(rle=unlist(RleList),
                                       dim=c(length(RleList[[1]]),
                                       length(RleList)),
                                       dimnames=dimnames
        )
    } else {
        rlem <- DelayedArray::RleArray(rle=unlist(RleList),
                                       dim=c(length(RleList[[1]]),
                                       length(RleList)))
    }
    return(rlem)

}


#' createGranges
#' @description a simplified wrapper function to create a GRanges object
#'
#' @param chrSeqInfo a seqinfo object
#' @param starts the start ranges
#' @param widths the width of each range
#' @param mcolname the name for the mcol attribute
#' @param mcolvalues the values for the mcol attribute
#'
#' @return a GRanges object
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom S4Vectors mcols
#' @export
#' @examples TBW
createGranges <- function(chrSeqInfo, starts, widths,
                          mcolname=NULL, mcolvalues=NULL) {
    stopifnot(is(chrSeqInfo, "Seqinfo"))
    stopifnot(identical(length(starts), length(widths)))

    gr <- GenomicRanges::GRanges(seqnames=as.character(chrSeqInfo@seqnames),
                                 ranges=IRanges::IRanges(start=starts, width=widths),
                                 seqinfo=chrSeqInfo)

    if(!is.null(mcolname) )
    {
        if(!is.null(mcolvalues)
           &&
           (length(gr@ranges@start) == length(mcolvalues))
        )
        {
            S4Vectors::mcols(gr)[[mcolname]] <- mcolvalues
        }
        else
        {
            warning("Cannot set mcols values! Vector length not matching Ranges")
        }
    }
    return(gr)
}


#' cutGRangesPerChromosome
#' @description  takes in input a GRanges object, producing a LIST of
#' GRanges, one for each chromosome
#'
#' @param GRanges a GRanges object
#'
#' @return a named list of GRanges, one for each chromosome
#' @importFrom GenomeInfoDb seqnames
#' @export
#'
#' @examples
#' gr <- GRanges(
#'       seqnames=Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'       ranges=IRanges(1:10, end=10),
#'       strand=Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)),
#'       seqlengths=c(chr1=11, chr2=12, chr3=13))
#' (grchrlist <- cutGRangesPerChromosome(gr))
#'
cutGRangesPerChromosome <- function(GRanges)
{
    stopifnot(is(GRanges, "GRanges"))
    interestedChrs <- GRanges@seqinfo@seqnames

    GRList <- lapply(interestedChrs, function(x)
    {
        bgr <- GRanges[ GenomeInfoDb::seqnames(GRanges) == x ]
        if(length(bgr) > 0)
        {
            bgr@seqinfo <- GRanges@seqinfo[x]
            GenomeInfoDb::seqnames(bgr) <- droplevels(GenomeInfoDb::seqnames(bgr))
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
#' cutGRangesPerChromosome returning only the relevant chromosomes GRanges
#'
#' @param GRangesList where each element is a chromosome,
#'                    tipically created with cutGRangesPerChromosome
#' @param chr a character vector of chromosomes names of the form "chr#"
#'
#' @return the input GRangesList with only the relevat chromosomes
#' @export
#' @examples TBW
keepRelevantChrs <- function(GRangesList, chr)
{
    if(!is.null(chr) && length(grep(pattern="chr", chr))!=length(chr))
        stop("Insert valid chr(s), use the \"chr#\" form!")
    stopifnot(is(GRangesList, "GRangesList"))

    idxs <- which(names(GRangesList) %in% chr)
    if(length(idxs) == 0)
        stop("Something went wrong in the chr subselection!",
             "\nPlease check the Chromosomes names!")
    GRangesList <- GRangesList[[idxs]]
    return(GRangesList)
}

#' fromSamplesToChromosomesGRangesList
#' @description converts a GRangesList orgnized per samples to a GRangesList
#'              organized per Chromosomes where each element
#'              is a GRangesList of samples
#' @param samplesGRangesList a GRangesList of samples.
#'                           Tipically generaed by findPeaks
#'
#' @return A GRangesList of chromosomes where each element is a GRanges list
#'         of samples
#' @export
#'
#' @examples TBW
fromSamplesToChromosomesGRangesList <- function(samplesGRangesList)
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
                idx <- grep(chr, names(samp))
                return(samp[[idx]]) ## it can be only one chr
            }))
        return(chrList)
    })
    names(chrsSamplesList) <- allChrs
    return(chrsSamplesList)
}

#' divideEachSampleByChromosomes
#' @description taken in input a grangeslist of samples, generate a list of samples
#'              where each element has a GRangesList each element of the
#'              GRangesList represents a single chromosome
#' @param samplesGRangesList a GRangesList of samples
#'
#' @return list of samples where each element is a list of chromosomes and each
#'          of these elements is a granges
#' @export
#'
#' @examples TBW
divideEachSampleByChromosomes <- function(samplesGRangesList)
{
    ## taken in input a grangeslist of samples, generate a list of samples
    ## where each element has a GRangesList
    ## each element of the GRangesList represents a single chromosome
    samplesChrList <- lapply(samplesGRangesList, function(gr)
    {   ## list of samples where each element is a list of
        ## chromosomes and each of these elements is a granges
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
#'              a string with its chromosomes separated by ";" (2nd col)
#'              (useful to fromSamplesToChromosomesGRangesList function)
#' @param samplesChrGRList a GRangesList of samples each divided by chromosome
#'
#' @return a dataframe  where each row is a sample (1st col) and
#'         a string with its chromosomes separated by ";" (2nd col)
#'
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
