## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  warning = FALSE
)

## ----setup1, eval = FALSE-----------------------------------------------------
#  devtools::install_github('sarahleavitt/nbTransmission')

## ----setup2-------------------------------------------------------------------
library('nbTransmission')

## ----seed---------------------------------------------------------------------
set.seed(0)

## ----indData------------------------------------------------------------------
data(indData)
str(indData)

## ----indToPair----------------------------------------------------------------
pairDataRaw <- indToPair(indData,
                         indIDVar = "individualID",
                         separator = "_",
                         dateVar = "infectionDate",
                         units = "days",
                         ordered = FALSE)
str(pairDataRaw)

## ----indCovariates------------------------------------------------------------
prop.table(table(indData$X1))
prop.table(table(indData$X2))
prop.table(table(indData$X3))
prop.table(table(indData$X4))

## ----pairCovariates-----------------------------------------------------------
prop.table(table(pairData$transmission, pairData$Z1), 1)
prop.table(table(pairData$transmission, pairData$Z2), 1)
prop.table(table(pairData$transmission, pairData$Z3), 1)
prop.table(table(pairData$transmission, pairData$Z4), 1)

## ----ordered------------------------------------------------------------------
pairDataOrdered <- pairData[pairData$infectionDate.2 >= pairData$infectionDate.1, ]

## ----snpClose-----------------------------------------------------------------
pairDataOrdered$snpClose <- ifelse(pairDataOrdered$snpDist < 3, TRUE,
                                   ifelse(pairDataOrdered$snpDist > 10, FALSE, NA))
table(pairDataOrdered$snpClose, useNA = "ifany")
prop.table(table(pairDataOrdered$snpClose, useNA = "ifany"))

## ----probabilities, results="hide"--------------------------------------------
resGen <- nbProbabilities(orderedPair = pairDataOrdered,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                            label = "SNPs", l = 1,
                            n = 10, m = 1, nReps = 1)

## ----probabilities2-----------------------------------------------------------
str(resGen)

## ----combine------------------------------------------------------------------
nbResults <- merge(resGen$probabilities, pairDataOrdered, by = "pairID", all = TRUE)
tapply(nbResults$pScaled, nbResults$transmission, summary)

## ----results, fig.width=7, fig.height=3---------------------------------------
library(ggplot2)
ggplot(data = nbResults, aes(x = pScaled)) +
  geom_histogram(bins = 20) +
  facet_wrap(~transmission, scales = "free_y")

## ----covarTab-----------------------------------------------------------------
library(knitr)
#Exponentiating the log odds ratios and creating a table of odds ratios
orTab <- resGen$estimates
orTab$orMean <- round(exp(orTab$logorMean), 2)
orTab$orCILB <- round(exp(orTab$logorCILB), 2)
orTab$orCIUB <- round(exp(orTab$logorCIUB), 2)
orTab$Value <- gsub("[A-z0-9]+\\:", "", orTab$level)
orTab$Variable <- gsub("\\:[A-z0-9+-<=>]+", "", orTab$level)
orTab <- orTab[, c("Variable", "Value", "orMean", "orCILB", "orCIUB")]

## ----covarFig, fig.height=5, fig.width=7--------------------------------------
ggplot(data = orTab, aes(x = Value, y = orMean, ymin = orCILB,
                           ymax = orCIUB)) +
  geom_point(size = 2) +
  geom_errorbar(width = 0.3) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  facet_wrap(~Variable, scales = "free_y") +
  ylab("Odds ratio with 95% confidence interval") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 360),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
        legend.position = "bottom") +
  scale_y_log10(breaks = c(0.01, 0.1, 1, 10)) +
  coord_flip()


## ----clustHC------------------------------------------------------------------
#Clustering the probabilities
clustHC <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "hc_absolute", cutoff = 0.05)
table(clustHC$cluster)
#Subsetting to just the top cluster
topClustHC <- clustHC[clustHC$cluster == 1, ]

## ----clustExample, fig.height=3, fig.width=7----------------------------------
ggplot(data = clustHC[clustHC$individualID.2 %in% c(32, 40), ],
        aes(x = pRank, y = pScaled, color = cluster, shape = transmission)) +
   geom_point() +
   facet_wrap(~individualID.2, scales = "free") +
   theme(legend.position = "none")

## ----clustKD------------------------------------------------------------------
#Clustering the probabilities
clustKD <- clusterInfectors(nbResults, indIDVar = "individualID", pVar = "pScaled",
                           clustMethod = "kd", cutoff = 0.01)
table(clustKD$cluster)
#Subsetting to just the top cluster
topClustKD <- clustKD[clustKD$cluster == 1, ]

## ----savepar, include=FALSE---------------------------------------------------
oldpar <- par(no.readonly=TRUE)

## ----networkFull, fig.height=5, fig.width=5-----------------------------------
par(mar = c(0, 0, 0.2, 0))
nbNetwork(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
pVar = "pScaled", clustMethod = "none")

## ----networkTop, fig.height=5, fig.width=5------------------------------------
par(mar = c(0, 0, 0.2, 0))
nbNetwork(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05)

## ----heatmap, fig.height=5, fig.width=5---------------------------------------
par(mar = c(0, 0, 1, 0))
nbHeatmap(nbResults, indIDVar = "individualID", dateVar = "infectionDate",
          pVar = "pScaled", clustMethod = "hc_absolute", cutoff = 0.05)


## ----rInitial-----------------------------------------------------------------
rInitial <- estimateR(nbResults,
                      dateVar = "infectionDate",
                      indIDVar = "individualID",
                      pVar = "pScaled",
                      timeFrame = "months")
str(rInitial)

## ----cutting, fig.height = 3.5, fig.width = 6---------------------------------
plotRt(rInitial)

## ----rFinal, results = "hide"-------------------------------------------------
rFinal <- estimateR(nbResults, dateVar = "infectionDate",
             indIDVar = "individualID", pVar = "pScaled",
             timeFrame = "months", rangeForAvg = c(25, 125),
             bootSamples = 10, alpha = 0.05)

## ----rFinal2------------------------------------------------------------------
rFinal$RtAvgDf

## ----fig.height = 3.5, fig.width = 6------------------------------------------
plotRt(rFinal, includeRtAvg = TRUE, includeRtCI = TRUE, includeRtAvgCI = TRUE)

## ----probNoT, results = "hide"------------------------------------------------
resGenNoT <- nbProbabilities(orderedPair = pairDataOrdered,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z1", "Z2", "Z3", "Z4"),
                            label = "SNPs", l = 1,
                            n = 10, m = 1, nReps = 1)
nbResultsNoT <- merge(resGenNoT[[1]], pairDataOrdered, by = "pairID", all = TRUE)

## ----siEstimate, results = "hide"---------------------------------------------
siPars <- estimateSI(nbResultsNoT,
                     indIDVar = "individualID",
                     timeDiffVar = "infectionDiffY",
                     pVar = "pScaled",
                     clustMethod = "hc_absolute",
                     cutoff = 0.05,
                     initialPars = c(2, 2))

## ----siEstimateRes------------------------------------------------------------
siPars

## ----siEstimateCI, results = "hide"-------------------------------------------
siParsCI <- estimateSI(nbResultsNoT,
                       indIDVar = "individualID",
                       timeDiffVar = "infectionDiffY",
                       pVar = "pScaled",
                       clustMethod = "hc_absolute",
                       cutoff = 0.05,
                       initialPars = c(2, 2),
                       bootSamples = 5)

## ----siEstimateCIRes----------------------------------------------------------
siParsCI

## ----siPlot, fig.height = 4, fig.width = 6, warning=FALSE---------------------
truePairs <- nbResultsNoT[nbResultsNoT$transmission == TRUE, ]
ggplot(data = nbResultsNoT, aes(x = infectionDiffY)) +
  geom_histogram(data = truePairs, aes(y = ..density..), bins = 40) +
  scale_x_continuous(name = "Serial Interval (years)", limits = c(0, 20)) +
  geom_line(aes(y = dgamma(infectionDiffY, shape = siPars$shape, scale = siPars$scale)))

## ----revertpar, include=FALSE-------------------------------------------------
par(oldpar)

