
#' buildValleys
#' @description It builds a GRanges of regions not present (Valleys)
#' in the @regionsGR parameter.
#' @param regionsGR a GRanges object
#' @param genomeCode the code of the genome to compare the regions to
#'
#' @return a GRanges object of the specular (respect to the genome) regions of
#' @regionsGR object
#'
#' @export
#'
#' @examples
#' filename <- system.file("extdata/regions/regions.rds", package="DEScan2")
#' regionsGR <- readRDS(file=filename)
#' (valleysGR <- buildValleys(regionsGR, "mm9"))
#'
buildValleys <- function(regionsGR, genomeCode)
{
    stopifnot(is(regionsGR, "GRanges"))
    ref.genome <- as(rtracklayer::SeqinfoForUCSCGenome(genomeCode), "GRanges")
    # si <- GenomeInfoDb::Seqinfo(genome=genomeCode)
    ref.genome <- DEScan2::keepRelevantChrs(chrGRangesList=ref.genome,
                                        chr=unique(seqnames(regionsGR)))
    valleys <- setdiff(ref.genome, regionsGR)

    return(valleys)
}


#
# reads.path <- system.file("extdata/bam", package="DEScan2")
# finalRegionsSE <- countFinalRegions(regionsGRanges=valleysGR,
#         readsFilePath=reads.path, fileType="bam", minCarriers=1,
#         genomeName="mm9", onlyStdChrs=TRUE, ignStrandSO=TRUE,
#         saveFlag=FALSE, verbose=TRUE)
#
# rowRanges(finalRegionsSE)

