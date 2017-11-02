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

    return(bed)
}



#' Constructs a GRanges object from a bam/bed file
#'
#' @param filename the complete file path of a bam?bed file
#' @param filetype the file type bam/bed
#' @param genomeName the name of the genome used to map the reads (i.e. mm9)
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
        bedRanges <- readBamAsBed(filename)
    } else {
        bedRanges <- readBedFile(filename)
    }
    # checking bed seqnames, useful in peak calling algorithm
    if( sum(is.na(seqinfo(bedRanges)@seqlengths)) > 0 )
    {
        message("No seqlengths present in file ", filename)
        if( !is.null(genomeName) )
        {
            genomeInfo <- Seqinfo(genome=genomeName)
            ## insert also the possibility to download the seqinfo on request
            seqNamesIdx <- which(genomeInfo@seqnames %in%
                                     seqinfo(bedRanges)@seqnames)
            if(length(seqNamesIdx) != 0)
            {
                bedRanges@seqinfo <-genomeInfo[genomeInfo@seqnames[seqNamesIdx]]
            }
            else
            {
                bedRanges@seqinfo <- genomeInfo
            }

        }
        else
        {
            stop("Please provide a genomeName in order to setup the GRanges!")
        }
    }

    return(bedRanges)
}
