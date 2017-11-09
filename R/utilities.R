#' readBamAsBed: read a bam file into a bed like format.
#'
#' @param file Character indicating path to bam file.
#' @return GenomicRanges object.
#'
#' @keywords internal
#' @import GenomicRanges
#' @import GenomicAlignments
#'
readBamAsBed <- function(file) {
    message("processing ", file)
    ga <- GenomicAlignments::readGAlignments(file, index=file)
    gr <- GenomicRanges::granges(x = ga)

    return(gr)
}


#' readBedFile: read a bed file into a GenomicRanges like format.
#'
#' @param file Character indicating path to bam file.
#' @return GenomicRanges object
#' @keywords internal
#' @import tools
#' @import GenomicRanges
#' @import rtracklayer
#'
readBedFile <- function(filename) {
    message("processing ", filename)
    if (tools::file_ext(filename) == "zip") {
        tmp <- utils::unzip(filename, list=T)$Name
        file <- unz(filename, tmp)
    } else {
        file <- filename
    }

    bed <- rtracklayer::import.bed(con = file)
    bed <- GenomicRanges::GRanges(seqnames=bed@seqnames,
                                  ranges=bed@ranges,
                                  strand=bed@strand)
    return(bed)
}



#' Constructs a GRanges object from a bam/bed file in a consistent way
#'
#' @param filename the complete file path of a bam?bed file
#' @param filetype the file type bam/bed
#' @param genomeName the name of the genome used to map the reads (i.e. mm9)
#'                   N.B. if NOT NULL the GRanges Seqinfo will be forced to
#'                   genomeName Seqinfo
#' @param onlySrdChrs a flag to keep only standard chromosomes from the input
#'                    files.
#'
#' @return bedRanges a GRanges object
#' @export
#'
#' @examples
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
        genomeInfo <- Seqinfo(genome=genomeName)
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
    if( (sum(is.na(seqinfo(bedGRanges)@seqlengths)) > 0) ||
        (length(seqinfo(bedGRanges)@seqnames) == 0 ))
    {
        stop("No seqlengths present in file ", filename,
             "\nPlease provide the correct genomeName to setup the GRanges!")
    }
    else if(length(uniqueSeqnames) < length(seqinfo(bedGRanges)@seqnames))
    {
        message("Keeping only necessary seqInfos")
        bedGRanges@seqinfo <- bedGRanges@seqinfo[as.character(uniqueSeqnames)]
        bedGRanges@seqnames <- droplevels(bedGRanges@seqnames)
    }

    return(bedGRanges)
}



#' Title
#'
#' @param GRanges
#' @param filepath
#' @param filename
#'
#' @return
#' @export
#'
#' @examples
saveGRangesAsBed <- function(GRanges, filepath, filename)
{
    stopifnot(is(GRanges, "GRanges"))
    ## add some parameters
    filePathName <- file.path(filepath, paste0(filename, "_peaks.bed"))
    if(file.exists(filePathName)) {stop(filePathName, " already exists!
                                        Not overwriting!")}
    rtracklayer::export.bed(object=GRanges, con=filePathName)
}


#' RleListToRleMatrix: a wrapper to create a RleMatrix from a RleList object
#'
#' @param RleList an RleList object
#' @param dimnames the names for dimensions of RleMatrix (see DelayedArray pkg)
#'
#' @return a RleMatrix from DelayedArray package
#' @export
#'
#' @examples
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


#' createGranges: a simplified wrapper function to create a GRanges object
#'
#' @param chrSeqInfo a seqinfo object
#' @param starts the start ranges
#' @param widths the width of each range
#' @param mcolname the name for the mcol attribute
#' @param mcolvalues the values for the mcol attribute
#'
#' @return a GRanges object
#' @export
#'
#' @examples
createGranges <- function(chrSeqInfo, starts, widths,
                          mcolname=NULL, mcolvalues=NULL) {
    stopifnot(is(chrSeqInfo, "seqinfo"))
    stopifnot(identical(length(starts), length(widths)))

    gr <- GRanges(seqnames=as.character(chrSeqInfo@seqnames),
                  ranges=IRanges(start=starts, width=widths),
                  seqinfo=chrSeqInfo
    )

    if(!is.null(mcolname) )
    {
        if(!is.null(mcolvalues)
           &&
           (length(gr@ranges@start) == length(mcolvalues))
        )
        {
            mcols(gr)[[mcolname]] <- mcolvalues
        }
        else
        {
            warning("Cannot set mcolvalues! Vector length not matching Ranges")
        }
    }
    return(gr)
}
