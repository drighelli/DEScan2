


atacDEScan2Peaks <- readFilesAsGRangesList(
    filePath="../../F1000_ATAC_analysis/results/atac_bwa/peaks/",
    fileType="bed", genomeName="mm10", onlyStdChrs=TRUE, arePeaks=TRUE)

names(atacDEScan2Peaks) <- gsub(".sorted.markdup.bam_zt10_mnw50_mxw1000_sw10_bin50.bed","", names(atacDEScan2Peaks))

valleysByFile <- GRangesList(lapply(atacDEScan2Peaks, function(peaksGR)
{
    buildValleys(regionsGR=peaksGR, genomeCode="mm10")
}))

valleyRegions <- finalRegions(peakSamplesGRangesList=valleysByFile,
            zThreshold=0,
            scorecolname="score",
            extendRegions=0,
            minCarriers=8,
            saveFlag=FALSE,
            verbose=TRUE)

bam.path <- "/Volumes/Elements/DEScan2/ATAC-seq_tim/bam_bwa/rmdup"
valleyCounts <- countFinalRegions(regionsGRanges=valleyRegions, readsFilePath=bam.path,
                  fileType="bam", minCarriers=1, genomeName="mm10",
                  onlyStdChrs=TRUE, saveFlag=FALSE,
                  verbose=TRUE)
aaa <- assay(valleyCounts)
colnames(aaa) <- basename(colnames(aaa))
colnames(aaa) <- gsub(pattern=".sorted.markdup.bam.rmdup.bam", replacement="", colnames(aaa))
plotPCA(aaa)

############## in alternativa a ipvalues usare i quantili della distribuzione
############## degli scores oppure "adattivo" calcolando i quantili per ogni sample
scores <- lapply(peaks, function(peak) {
    return(peak$score)
})

scores <- unlist(scores)

quantile(scores)
