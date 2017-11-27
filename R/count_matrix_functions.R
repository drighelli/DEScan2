#' count reads falling within the final regions.
#'
#' @param regionsGRanges
#' @param readsFilePath
#' @param fileType
#' @param minCarriers
#' @param genomeName
#' @param onlyStdChrs
#' @param saveFlag
#' @param savePath
#'
#' @return Matrix containing read counts with regions as rows and samples as
#'     columns.
#' @export
#' @importFrom utils read.table write.table
#' @import Rsubread
countFinalRegions <- function(regionsGRanges, readsFilePath,
                              fileType=c("bam", "bed"),
                              minCarriers=2,
                              #minCoverage=0,
                              genomeName=NULL,
                              onlyStdChrs=FALSE,
                              saveFlag=FALSE,
                              savePath="finalRegions")
{
    match.arg(fileType)
    stopifnot(is(regionsGRanges, "GRanges"))

    idxK <- which(regionsGRanges$`k-carriers` >= minCarriers)
    regionsGRanges <- regionsGRanges[idxK,]

    readsFiles <- list.files(readsFilePath, full.names=TRUE, pattern=fileType)
    idxBai <- grep("bai", readsFiles)
    if(length(idxBai) > 0) readsFiles <- readsFiles[-idxBai]

    summRegDF <- #plyr::adply(readsFiles, 1, function(file)
        sapply(readsFiles, function(file)
    {
        fileReads <- constructBedRanges(filename=as.character(file),
                                        filetype=fileType,
                                        genomeName=genomeName,
                                        onlyStdChrs=onlyStdChrs)
        summReg <- GenomicAlignments::summarizeOverlaps(features=regionsGRanges,
                                                        reads=fileReads)
        ##the mode is to compare with the feature counts method outpuut
        return(SummarizedExperiment::assay(summReg))
        ##not working properly
    })

    # if(minCoverage != 0)
    # {
    #     which(summRegDF > minCoverage, arr.ind=TRUE)
    # }

    ### saving functionality

    return(summRegDF)
    ########  old version modified
    #  regions<-read.table(region_file, sep="\t", header=TRUE)
    #  regions<-regions[regions[[5]] >= min_carriers, ]
    #
    #  region_anno <- as.data.frame(cbind("GeneID"=names(regionsGRanges),
    #                     "Chr"=as.vector(GenomeInfoDb::seqnames(regionsGRanges)),
    #                     "Start"=BiocGenerics::start(regionsGRanges),
    #                     "End"=BiocGenerics::end(regionsGRanges),
    #                     "Strand"=rep("*", length(regionsGRanges))))
    #
    #  # region_anno<-cbind("GeneID"=rownames(regions), regions[,1:3],
    #                                       # "Strand"=rep("*", dim(regions)[1]))
    #  # colnames(region_anno)[2:4]<-c("Chr", "Start", "End")
    #  rnames<-paste(region_anno$Chr, paste(region_anno$Start, region_anno$End, sep="-"), sep=":")
    #  rownames(region_anno)<-rnames
    #  if (verbose == FALSE) {
    #          dummy<-capture.output(
    #                  count<-Rsubread::featureCounts(bam_files,
    #                                                 annot.ext=region_anno)
    #                  )
    #  }
    #  else
    # {
    #      count<-Rsubread::featureCounts(readsFiles, annot.ext=region_anno)
    #  }
    #  countmat<-count$counts
    #  rownames(countmat)<-rnames
    #  countmat<-countmat[, sort(colnames(countmat))]
    #  utils::write.table(countmat, file="countmatrix.txt", sep="\t",
    #                     row.names=TRUE, col.names=TRUE, quote=FALSE)
    #  countmat
}
