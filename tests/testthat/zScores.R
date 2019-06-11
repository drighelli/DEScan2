context("Z Scores")

#' computeZnbLRT
#'
#' @param mu1 a matrix m x n
#' @param mu2 a matrix m x n
#' @param theta a matrix m x n (vector of theta has to be replicated m times)
#' @param minCount a minimal constant to avoid zeros computations
#'
#' @return a matrix m x n of z values
#' @keyword internal
computeZnbLRT <- function(mu1, mu2, theta, minCount=0.1)
{
    zlrt <-  sqrt(2) * sign(mu1 - mu2) *
        sqrt(
            log( (theta + mu2) / (theta + mu1) ) *
                (theta + mu1) -
                mu1 * log(mu2 / pmax(mu1, minCount))
        )
    zlrt
}


computeZLRT <- function(mu1, mu2, minCount=0.1)
{
    z <- sqrt(2) * sign(mu1 - mu2) *
        sqrt(mu1 *
                 log(pmax(mu1, minCount) / mu2) -
                 (mu1 - mu2)
        )
    z
}


z.inputs <- load(system.file(
                file.path("extdata","tests","z", "Z_score_inputs.RData"),
                package="DEScan2"))

chrLength=seqlength


lambdaChrRleMm <- matrix(unlist(lambdaChrRleList),
                         ncol=length(lambdaChrRleList), byrow=FALSE)

# runWinRleMm <- matrix(unlist(runWinRleList), ncol=20, byrow=TRUE)
runWinRleMm <- matrix(unlist(runWinRleList),
                      ncol=length(runWinRleList), byrow=FALSE)


if(verbose) message("Computing Z-Score")
z <- sqrt(2) * sign(runWinRleMm - lambdaChrRleMm) *
    sqrt(runWinRleMm *
             log(pmax(runWinRleMm, minCount) / lambdaChrRleMm) -
             (runWinRleMm - lambdaChrRleMm)
    )
z <- binToChrCoordMatRowNames(binMatrix=z,
                                chrLength=chrLength,
                                binWidth=binSize)


if(verbose) message("Computing NB Z-Score")
if(verbose) message("estimating dispersions")
matwin <- as.matrix(runWinRleList)
phiWins1 <- edgeR::estimateDisp(matwin,
                               # lib.size=pmax(apply(matwin, 2, sum) != 0, minCount))
                               lib.size = rep(1, ncol(matwin)),
                               trend.method = "none")
thetaWins <- 1/phiWins1$tagwise.dispersion
if(verbose) message("computing scores")
zNB <-  sqrt(2) * sign(runWinRleMm - lambdaChrRleMm) *
    sqrt(
        log( (thetaWins + lambdaChrRleMm) / (thetaWins + runWinRleMm) ) *
            (thetaWins + runWinRleMm) -
            runWinRleMm * log(lambdaChrRleMm/ pmax(runWinRleMm, minCount))
    )

zNB <- binToChrCoordMatRowNames(binMatrix=zNB,
                              chrLength=chrLength,
                              binWidth=binSize)


mu2MeanWins <- apply(lambdaChrRleMm, 2, mean)
nobs <- nrow(runWinRleMm) # num of bins
nbinom <- lapply(c(1:20), function(i) {
    rnbinom(n=nobs, size=thetaWins[i], mu=mu2MeanWins[i])
})
nbinom <- matrix(unlist(nbinom), ncol=ncol(lambdaChrRleMm), byrow=TRUE)

mu1MeanWins <- apply(runWinRleMm, 2, mean)
nobs <- nrow(runWinRleMm) # num of bins
nbinom1 <- lapply(c(1:20), function(i) {
    rnbinom(n=nobs, size=thetaWins[i], mu=mu1MeanWins[i])
})

nbinom1 <- matrix(unlist(nbinom1), ncol=ncol(lambdaChrRleMm), byrow=TRUE)

zNBrand <-  sqrt(2) * sign(nbinom1 - nbinom) *
    sqrt(
        log( (thetaWins + nbinom) / (thetaWins + nbinom1) ) *
            (thetaWins + nbinom1) -
            nbinom1 * log(nbinom/ pmax(runWinRleMm, minCount))
    )
zNBrand

myTheta <- c(0.0001, 0.001, 0.01, 0.1, 1, 5, 10, 50, 100, 250, 500, 750, 1000,
            2500, 5000, 10000, 50000, 100000, 500000, 1000000)

lambda1ChrRleMm <- matrix(rep(lambdaChrRleMm[,1], 20), ncol=20, byrow=F)
runWin1RleMm <- matrix(rep(runWinRleMm[,1], 20), ncol=20, byrow=F)


############ THETA DEVE ESSERE MOLTIPLICATO NEL MODO GIUSTO. creare matrice di theta?


myThetaNew <- matrix(rep(myTheta, dim(lambdaChrRleMm)[1]), ncol=20, nrow=dim(lambdaChrRleMm)[1], byrow=TRUE)
mythetaznbfake <- computeZnbLRT(mu1=runWin1RleMm, mu2=lambda1ChrRleMm, theta=myThetaNew)
summary(mythetaznbfake)
mythetaznb <- computeZnbLRT(mu1=runWinRleMm, mu2=lambdaChrRleMm, theta=myThetaNew)
# mythetaznbaa <- computeZnbLRT(mu1=runWinRleMm, mu2=lambdaChrRleMm, theta=myTheta)
summary(mythetaznb)
mythetaznbmu1 <- computeZnbLRT(mu1=nbinom1, mu2=lambdaChrRleMm, theta=myTheta)
summary(mythetaznbmu1)
mythetaznbmu2 <- computeZnbLRT(mu1=nbinom, mu2=lambdaChrRleMm, theta=myTheta)
summary(mythetaznbmu2)
mythetaznbmu1mu2 <- computeZnbLRT(mu1=nbinom1, mu2=nbinom, theta=myTheta)
summary(mythetaznbmu1mu2)


oldZ <- computeZLRT(mu1=runWinRleMm, mu2=lambdaChrRleMm)
summary(oldZ)
oldZmu1 <- computeZLRT(mu1=nbinom1, mu2=lambdaChrRleMm)
summary(oldZmu1)
oldZmu2 <- computeZLRT(mu1=runWinRleMm, mu2=nbinom)
summary(oldZmu2)
oldZmu1mu2 <- computeZLRT(mu1=nbinom1, mu2=nbinom)
summary(oldZmu1mu2)
