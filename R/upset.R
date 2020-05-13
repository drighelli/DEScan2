

#' plotUpsetRanges
#'
#' @description It computes a matrix of 0/1 by using all the regions within
#' each GRanges of the @namedGRList.
#' The unified regions varies by the @methodRows selected, it could be one of
#' merge, intersect (as defined in GenomicRanges), union (as defined in
#' GenomicRanges).
#' It assigns 1 in the matrix when the findOverlaps (with default parameters)
#' method finds an overlap between the unified region on the row (query)
#' and the GRanges represented by the column (subject)..
#'
#'
#' @param namedGRList a named GRangesList to plot
#' @param labelsList optional alternative list of labels (default is NULL)
#' @param methodRows the method to use for computing the intersections.
#' @param showPlot boolean for showing (TRUE) or returning (FALSE) the upset.
#' @param title a title to add to the plot. Works only with showPlot=TRUE.
#' (Waiting for upset function to implement title by itself).
#' @param ... additional parameters to pass to the UpSetR::upset function

#'
#' @return an upset object if the @showPlot flag is FALSE, NULL otherwhise
#'
#' @importFrom GenomicRanges insersect union findOverlaps
#' @importFrom UpSetR upset
#' @importFrom S4Vectors from
#' @export
#'
#' @examples
#' tbd
plotUpsetRanges <- function(namedGRList, labelsList=NULL,
                            methodRows=c("merge", "intersect", "union"),
                            showPlot=TRUE, title=NULL, ...)
{
    stopifnot(is(namedGRList, "GRangesList"))
    stopifnot(!is.null(names(namedGRList)))
    if(!is.null(labelsList)) names(namedGRList) <- labelsList
    match.arg(methodRows)

    unitedGRanges <- unlist(namedGRList,"GRangesList")
    names(unitedGRanges) <- NULL
    switch(methodRows,
           merge={
               unitedGRanges <- unique(unitedGRanges)
           },
           intersect={
               unitedGRanges <- GenomicRanges::intersect(unitedGRanges,
                                                         unitedGRanges)
           },
           union={
               unitedGRanges <- GenomicRanges::union(unitedGRanges, unitedGRanges)
           },
           {
               stop("Weird methodRows value found!")
           }
    )

    naMat <- matrix(nrow=length(unitedGRanges), ncol=length(namedGRList))
    overlapsList <- lapply(seq_along(namedGRList), function(i)
    {
        overl <- GenomicRanges::findOverlaps(query=unitedGRanges,
                                            subject=namedGRList[[i]])
        naMat[S4Vectors::from(overl), i] <- 1
        return(naMat[,i])
    })

    overlapsMat <- matrix(unlist(overlapsList), ncol=length(namedGRList))
    colnames(overlapsMat) <- names(namedGRList)
    idx <- which(is.na(overlapsMat), arr.ind=TRUE)
    if(length(idx) != 0) overlapsMat[idx] <- 0

    ups <- UpSetR::upset(as.data.frame(overlapsMat), ... )

    if(showPlot)
    {
        print(ups)
        if(!is.null(title)) ## to remove in case of title update of upset funct
        {
            grid::grid.text(title, x=0.50, y=0.97, gp=grid::gpar(fontsize=20))
        }
    } else {
        return(ups)
    }
}


#' plotUpsetGenes
#'
#' @description It computes a matrix of 0/1 by using all the genes within
#' each element of the @namedList.
#' It computes a unique list of genes from the lists inside the @namedList
#' argument.
#'
#' @param namedList a named list of elements containing gene lists.
#' @param labelsList optional alternative list of labels (default is NULL)
#' @param showPlot boolean for showing (TRUE) or returning (FALSE) the upset.
#' @param title a title to add to the plot. Works only with showPlot=TRUE.
#' (Waiting for upset function to implement title by itself).
#' @param ... additional parameters to pass to the UpSetR::upset function
#'
#' @return an upset object if the @showPlot flag is FALSE, NULL otherwhise
#'
#' @importFrom UpSetR upset
#' @export
#'
#' @examples
#' tbd
plotUpsetGenes <- function(namedList, labelsList=NULL,
                            showPlot=TRUE, title=NULL, ...)
{
    stopifnot(!is.null(names(namedList)))
    if(!is.null(labelsList)) names(namedList) <- labelsList

    unitedGenes <- unique(unlist(namedList))


    geneMat <- matrix(nrow=length(unitedGenes), ncol=length(namedList), data=0)
    rownames(geneMat) <- unitedGenes
    for(i in 1:length(namedList))
    {
        geneMat[rownames(geneMat) %in% namedList[[i]], i] <- 1
    }

    colnames(geneMat) <- names(namedList)

    ups <- UpSetR::upset(as.data.frame(geneMat), ...)

    if(showPlot)
    {
        print(ups)
        if(!is.null(title)) ## to remove in case of title update of upset funct
        {
            grid.text(title, x=0.50, y=0.97, gp=gpar(fontsize=20))
        }
    } else {
        return(ups)
    }
}


