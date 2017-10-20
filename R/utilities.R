#' read a bam file into a bed like format.
#'
#' @param file Character indicating path to bam file.
#' @param chr Integer indicating which chromosome to read in.
#' @return GenomicRanges object.
#' @keywords internal
#' @import GenomicRanges
#' @import Rsamtools
#' @import GenomeInfoDb
#' @import GenomicAlignments
#' @import IRanges
readBamAsBed <- function(file, chr) {
    ga <- GenomicAlignments::readGAlignments(file, index=file)
    gr <- GenomicRanges::granges(x = ga)
    gr <- gr[which(gr@seqnames %in% chr),]

    if(length(gr) == 0) {
        stop( file, " bam file doesn't contain ",
              chr, " chromosome(s).\nExiting." )
    }
    return(gr)
}


#' read a bed file into a GenomicRanges like format.
#'
#' @param file Character indicating path to bam file.
#' @param chr Integer indicating which chromosome(s) to read in.
#' @return GenomicRanges object
#' @keywords internal
#' @import tools
#' @import GenomicRanges
#' @import rtracklayer
readBedFile <- function(filename, chr) {

    if (tools::file_ext(filename) == "zip") {
        tmp <- utils::unzip(filename, list=T)$Name
        file <- unz(filename, tmp)
    } else {
        file <- filename
    }

    bed <- rtracklayer::import.bed(con = file)
    bed <- bed[which(bed@seqnames %in% chr), ]
    if(length(bed) == 0) {
        stop( filename, " bed file doesn't contain ",
              chr, " chromosome(s).\nExiting." )
    }

    bed <- GRanges(seqnames=bed@seqnames,
                   ranges=bed@ranges,
                   strand=bed@strand)
    return(bed)
}
