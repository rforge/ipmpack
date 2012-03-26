\name{makeSurvObjManyCov}
\alias{makeSurvObjManyCov}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Functions to build survival objects}
\description{
  A function to fit logistic regressions defining
  survival following user defined formulas (e.g., size+size^2, etc)
  to build survival objects for which methods to build an IPM  object are defined.
}
\usage{
makeSurvObjManyCov(dataf,explanatoryVariables="size+size2+covariate1")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{dataf}{a dataframe with columns `size' and
    `surv' (`size' is size at t, `surv' is survival to t+1 coded as 0
    or 1); dataf may include a number of covariates for which
    column names must follow the sequence `covariate1', `covariate2', ...}
    \item{explanatoryVariables}{a formula describing the desired
      explanatory variables (interactions, etc) in classical R
      style, i.e. separated by `+', `*', `:'. Possible covariates
      include `size', 'size2' (size^2), `size3'
      (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
      `covariate1' `covariate2', etc.}
   }

\value{
  An object of class survObjManyCov
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
  See manual for details on building case-specific survival objects. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{makeSurvObj}}
}
\examples{
#generate data
dff <- generateData()
#make simple logistic regression survival object relating survival to size at t
sv1 <- makeSurvObjManyCov(dff, explanatoryVariables="size")
#assess fit for model with discrete environmental classes fitted
dff$covariate1 <- sample(1:3,size=nrow(dff),replace=TRUE)
sv1 <- makeSurvObjManyCov(dff, explanatoryVariables="size+covariate1")

}
