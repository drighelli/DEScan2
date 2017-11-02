#' read a bam file into a bed like format.
#'
#' @param file Character indicating path to bam file.
#' @return GenomicRanges object.
#' @keywords internal
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import IRanges
readBamAsBed <- function(file) {
    message("processing ", filename)
    ga <- GenomicAlignments::readGAlignments(file, index=file)
    gr <- GenomicRanges::granges(x = ga)

    return(gr)
}


#' read a bed file into a GenomicRanges like format.
#'
#' @param file Character indicating path to bam file.
#' @return GenomicRanges object
#' @keywords internal
#' @import tools
#' @import GenomicRanges
#' @import rtracklayer
readBedFile <- function(filename) {
    message("processing ", filename)
    if (tools::file_ext(filename) == "zip") {
        tmp <- utils::unzip(filename, list=T)$Name
        file <- unz(filename, tmp)
    } else {
        file <- filename
    }

    bed <- rtracklayer::import.bed(con = file)
    bed <- GRanges(seqnames=bed@seqnames,
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
#'
#' @return bedRanges a GRanges object
#' @export
#'
#' @examples
constructBedRanges <- function(filename
                               , filetype=c("bam", "bed")
                               , genomeName=NULL)
{
    filetype <- match.arg(filetype)

    if(filetype == "bam") {
        bedGRanges <- readBamAsBed(filename)
    } else {
        bedGRanges <- readBedFile(filename)
    }

    uniqueSeqnames <- droplevels(unique(bedGRanges@seqnames))

    if( !is.null(genomeName) )
    {

        message("Taking seqlengths from genome ", genomeName)
        genomeInfo <- Seqinfo(genome=genomeName)
        seqNamesIdx <- which(genomeInfo@seqnames %in% uniqueSeqnames)
        if(length(seqNamesIdx) != 0)
        {
            bedGRanges@seqinfo<-genomeInfo[genomeInfo@seqnames[seqNamesIdx]]
        }
        else
        {
            stop("Cannot find the ", uniqueSeqnames, " in genome ", genomeName)
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
        message("Keeping just necessary seqInfos")
        bedGRanges@seqinfo <- bedGRanges@seqinfo[as.character(uniqueSeqnames)]
    }

    return(bedGRanges)
}
