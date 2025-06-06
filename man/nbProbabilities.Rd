% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcProbabilities.R
\name{nbProbabilities}
\alias{nbProbabilities}
\title{Estimates relative transmission probabilities}
\usage{
nbProbabilities(
  orderedPair,
  indIDVar,
  pairIDVar,
  goldStdVar,
  covariates,
  label = "",
  l = 1,
  n = 10,
  m = 1,
  nReps = 10,
  orType = "univariate",
  nBS = 100,
  pSampled = 1,
  progressBar = TRUE
)
}
\arguments{
\item{orderedPair}{The name of the ordered pair-level dataset with the covariates.}

\item{indIDVar}{The name (in quotes) of the column with the individual ID.
(data frame \code{orderedPair} must have columns called \code{<indIDVar>.1} and \code{<indIDVar>.2}).}

\item{pairIDVar}{The name (in quotes) of the column with the unique pair ID variable.}

\item{goldStdVar}{The name (in quotes) of the column with a logical vector defining
training links/non-links}

\item{covariates}{A character vector containing the covariate column names (in quotes).
All covariates need to be categorical factor variables.}

\item{label}{An optional label string for the run.}

\item{l}{Laplace smoothing parameter that is added to each cell.}

\item{n}{The number of folds for nxm cross validation (should be at least 10).}

\item{m}{The number of times to create n folds in nxm cross validation.}

\item{nReps}{The number of times to randomly select the "true" infector (should be at least 10).}

\item{orType}{Takes value \code{"univariate"} or \code{"adjusted"}. \code{"univariate"} produces contingency table
odds ratios and \code{"adjusted"} produces adjusted odds ratios from a bootstrapped multivariable logistic regression.}

\item{nBS}{Number of bootstrap samples to run in each cross-validation fold/iteration (default is 100). Only
relevant when \code{orType = "adjusted"}.}

\item{pSampled}{Proportion of unlinked cases to include in bootstrap sample (default is 1, i.e.a true
bootstrap). Only relevant when \code{orType = "adjusted"}.}

\item{progressBar}{A logical indicating if a progress bar should be printed (default is TRUE).}
}
\value{
List containing two data frames:
\enumerate{
  \item \code{probabilities} - a data frame of transmission probabilities. Column names:
     \itemize{
       \item \code{label} - the optional label of the run.
       \item \code{<pairIDVar>} - the pair ID with the name specified.
       \item \code{pAvg} - the mean transmission probability for the pair over all iterations.
       \item \code{pSD} - the standard deviation of the transmission probability for the pair
        over all iterations.
       \item \code{pScaled} - the mean relative transmission probability for the pair over.
        all iterations: pAvg scaled so that the probabilities for all infectors per infectee add to 1.
       \item \code{pRank} - the rank of the probability of the the pair out of all pairs for that
       infectee (in case of ties all values have the minimum rank of the group).
       \item \code{nEstimates} - the number of probability estimates that contributed to pAvg. This
       represents the number of prediction datasets this pair was included in over the \code{nxm}
       cross prediction repeated \code{nReps} times.
     }
  \item \code{estimates} - a data frame with the contribution of covariates. Column names:
     \itemize{
       \item \code{label} - the optional label of the run
       \item \code{level} - the covariate name and level
       \item \code{nIter} - the number of iterations included in the estimates: \code{n*m*nReps}
       \item \code{logorMean} - the mean value of the log odds ratio across iterations
       \item \code{logorSE} - the standard error of the log odds ratio across iterations
       \item \code{logorCILB} - the lower bound of the 95% confidence interval of the log odds ratio
        across iterations
       \item \code{logorCIUB} - the upper bound of the 95% confidence interval of the log odds ratio
        across iterations
     }
}
}
\description{
The function \code{nbProbabilities} uses naive Bayes and an interative estimation
procedure to estimate relative transmission probabilities
}
\details{
This algorithm takes a dataset of ordered possible infector-infectee pairs in an
infectious disease outbreak or cluster and estimates the relative probability the cases are
linked by direct transmission using a classification technique called naive Bayes (NB).
NB is a simple machine learning algorithm that uses Bayes rule to estimate the
probability of an outcome in a prediction dataset given a set of covariates from
the observed frequencies in a training dataset.

The input dataset - \code{orderedPair} - should represent ordered pairs of cases
(where the potential infector was observed before the infectee) and have
a unique identifier for each pair (\code{pairIDVar}) as well as the individual ids that are
included in the pair (\code{<indIDVar>.1} and \code{<indIDVar>.2}). If cases are concurrent
(meaning the order cannot determined) both orders can be included.

A subset of pairs should also have pathogen WGS, contact investigation,  or some other
'gold standard' defined by \code{goldStdVar} which should be a logical vector with
\code{TRUE} indicating links, \code{FALSE} nonlinks, and \code{NA} if
the pair cannot be used to train (does not have the information or is indeterminate).
These pairs will be used to a training dataset of probable links and non/links.
The covariates can be any categorical variables and could represent
spatial, clinical, demographic, and temporal characteristics of the case pair.

Because the outcomes in the training set represent probable and not certain
transmission events and a given case could have mulitple probable infectors,
the algorithm uses an iterative estimation procedure. This procedure randomly chooses one
link of all of the possible links to include in the training dataset \code{nReps}
times, and then uses \code{mxn} cross prediction to give all pairs a turn
in the prediction dataset.

The output of this function is a list of two dataframes: one with the estimates of the
transmission probabilities (\code{probabilities}) and the other with the contribution of
the covariates to the probabilities in the form of odds ratios (\code{estimates}). The
95\% confidence intervals reported for these odds ratios use Rubin's Rules, a technique developed
for multiple imputation, to pool the error across all iterations.

This function generates odds ratios describing the associations between covariates in the training data
and outcome defined in the gold standard variable (\code{goldStdVar}) argument. Unadjusted odds ratios are the default.
These odds ratios are produced using contingency table methods. Adjusted odds ratios are calculated via bootstrapped
logistic regression to produce non-parametric standard errors. The bootstrap is controlled by parameters \code{nBS},
the number of bootstrap samples to run, and \code{pSampled}, the proportion of unlinked cases to include in the bootstrap
sample. \code{pSampled} is recommended only for large datasets in which it is computationally unfeasible to run a full
bootstrap. Sensitivity analyses should be run to determine an adequate value for \code{pSampled}.
}
\examples{
## Use the pairData dataset which represents a TB-like outbreak
# First create a dataset of ordered pairs
orderedPair <- pairData[pairData$infectionDiffY >= 0, ]

## Create a variable called snpClose that will define probable links
# (<3 SNPs) and nonlinks (>12 SNPs) all pairs with between 2-12 SNPs
# will not be used to train.
orderedPair$snpClose <- ifelse(orderedPair$snpDist < 3, TRUE,
                        ifelse(orderedPair$snpDist > 12, FALSE, NA))
table(orderedPair$snpClose)

## Running the algorithm
#NOTE should run with nReps > 1.
resGen <- nbProbabilities(orderedPair = orderedPair,
                            indIDVar = "individualID",
                            pairIDVar = "pairID",
                            goldStdVar = "snpClose",
                            covariates = c("Z1", "Z2", "Z3", "Z4", "timeCat"),
                            label = "SNPs", l = 1,
                            n = 10, m = 1, nReps = 1)

## Merging the probabilities back with the pair-level data
nbResults <- merge(resGen[[1]], orderedPair, by = "pairID", all = TRUE)

}
\references{
Barnard J. and Rubin D. Small-Sample Degrees of Freedom with Multiple Imputation
\emph{Biometrika}. 1999 Dec;86(4):948-55.
}
