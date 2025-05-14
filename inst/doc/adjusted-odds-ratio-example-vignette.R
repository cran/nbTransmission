## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup1, eval = FALSE-----------------------------------------------------
#  devtools::install_github('sarahleavitt/nbTransmission')

## ----setup--------------------------------------------------------------------
library(nbTransmission)

## ----seed---------------------------------------------------------------------
set.seed(0)

## ----makeData-----------------------------------------------------------------
#load individual level data
data(pairData)

#create ordered pair level dataset
pairDataOrdered <- pairData[pairData$infectionDate.2 >= pairData$infectionDate.1, ]

#create SNP distances thresholds 
#SNP distances less than 3 are considered transmission links
#SNP distances greater than 10 are non-transmission links
# SNP distances 4-10 are considered indeterminate and only used in the prediction dataset
pairDataOrdered$snpClose <- ifelse(pairDataOrdered$snpDist < 3, TRUE,
                                   ifelse(pairDataOrdered$snpDist > 10, FALSE, NA))

## ----probabilities, results="hide"--------------------------------------------
unadjRes <- nbProbabilities(orderedPair = pairDataOrdered,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                            label = "SNPs", l = 1,
                            n = 5, m = 1, nReps = 5, 
                            orType = "univariate")

## ----probabilities2-----------------------------------------------------------
str(unadjRes)

## ----covarTab-----------------------------------------------------------------
library(knitr)
#Exponentiating the log odds ratios and creating a table of odds ratios
orTab <- unadjRes$estimates
orTab$orMean <- round(exp(orTab$logorMean), 2)
orTab$orCILB <- round(exp(orTab$logorCILB), 2)
orTab$orCIUB <- round(exp(orTab$logorCIUB), 2)
orTab$Value <- gsub("[A-z0-9]+\\:", "", orTab$level)
orTab$Variable <- gsub("\\:[A-z0-9+-<=>]+", "", orTab$level)
orTabPrint1 <- orTab[, c("Variable", "Value", "orMean", "orCILB", "orCIUB")]

## ----covarFig, fig.height=5, fig.width=7--------------------------------------
library(ggplot2)
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


## ----CI90pct------------------------------------------------------------------
orTab$logor90CILB <- orTab$logorMean - 1.645*orTab$logorSE
orTab$logor90CIUB <- orTab$logorMean + 1.645*orTab$logorSE
orTab$or90CILB <- exp(orTab$logor90CILB)
orTab$or90CIUB <- exp(orTab$logor90CIUB)

## ----sig90--------------------------------------------------------------------
orTab$sig90pctCI <- ifelse(orTab$or90CILB < 1 & orTab$or90CIUB > 1, "Not significant", "Significant")

orTabPrint2 <- orTab[, c("Variable", "Value", "orMean", "or90CILB", "or90CIUB", "sig90pctCI")]


## ----adj_probabilities, results="hide"----------------------------------------
adjRes <- nbProbabilities(orderedPair = pairDataOrdered,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z2", "Z4", "timeCat"),
                            label = "SNPs", l = 1,
                            n = 5, m = 1, nReps = 5, 
                            orType = "adjusted", nBS = 10)

## ----adj_covarTab-------------------------------------------------------------
#Exponentiating the log odds ratios and creating a table of odds ratios
orTabAdj <- adjRes$estimates
orTabAdj <- orTabAdj[-c(orTabAdj$level == "(Intercept)"), ]
orTabAdj$orMean <- round(exp(orTabAdj$logorMean), 2)
orTabAdj$orCILB <- round(exp(orTabAdj$logorCILB), 2)
orTabAdj$orCIUB <- round(exp(orTabAdj$logorCIUB), 2)
orTabAdj$Variable <- substr(orTabAdj$level, 1, nchar(orTabAdj$level)-1)
orTabAdj$Value <- substr(orTabAdj$level, nchar(orTabAdj$level), nchar(orTabAdj$level))
orTabAdjPrint <- orTabAdj[, c("Variable", "Value", "orMean", "orCILB", "orCIUB")]

## ----timeCatTab---------------------------------------------------------------
table(pairDataOrdered$timeCat, pairDataOrdered$snpClose)


## ----timeCatfix---------------------------------------------------------------
pairDataOrdered$timeCat2 <- as.factor(ifelse(pairDataOrdered$timeCat == 6, 5, pairDataOrdered$timeCat ))

## ----adj_probabilities_rerun, results="hide"----------------------------------
adjRes2 <- nbProbabilities(orderedPair = pairDataOrdered,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z2", "Z4", "timeCat2"),
                            label = "SNPs", l = 1,
                            n = 5, m = 1, nReps = 5, 
                            orType = "adjusted", nBS = 10)

orTabAdj2 <- adjRes2$estimates
orTabAdj2 <- orTabAdj2[-c(orTabAdj2$level == "(Intercept)"), ]
orTabAdj2$orMean <- round(exp(orTabAdj2$logorMean), 2)
orTabAdj2$orCILB <- round(exp(orTabAdj2$logorCILB), 2)
orTabAdj2$orCIUB <- round(exp(orTabAdj2$logorCIUB), 2)
orTabAdj2$Variable <- substr(orTabAdj2$level, 1, nchar(orTabAdj2$level)-1)
orTabAdj2$Value <- substr(orTabAdj2$level, nchar(orTabAdj2$level), nchar(orTabAdj2$level))
orTabAdjPrint2 <- orTabAdj2[, c("Variable", "Value", "orMean", "orCILB", "orCIUB")]

