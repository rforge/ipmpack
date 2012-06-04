\name{trackPopStructManyCov}
\alias{trackPopStructManyCov}
\alias{plotResultsStochStruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Tracks population structure in a stochastic environment.
}
\description{
  Provided a set of time-varying covariates, projects the population forward, and stores state of the population. 
}
\usage{
trackPopStructManyCov(covariate, nRunIn, tMax, growthObj, survObj, fecObj, 
nBigMatrix, minSize, maxSize, nMicrosites, integrateType, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covariate}{
    matrix with tMax rows, and as many columns as there are relevant covariates.
}
  \item{nRunIn}{
    numeric, number of initial samples to discard.
}
  \item{tMax}{
    numeric, total number of time-steps to run (same as ncol(covariate)).
}
  \item{growthObj}{
    a growth object, defined to correspond to covariate definition
    (indexing used to make the growth object must match up).
}
  \item{survObj}{
    a survival object, defined to correspond to covariate definition in covariate.
}
  \item{fecObj}{
    a fecundity object, defined to correspond to covariate definition in covariate.
}
  \item{nBigMatrix}{
    numeric, number of size bins in the IPM.
}
  \item{minSize}{
    numeric, minimum size in the IPM
}
  \item{maxSize}{
    numeric, maximum size in the IPM.
}
  \item{nMicrosites}{
    vector, if sum(nMicrosites) > 0 then density dependence is assumed to operate on seedling establishment, and if length(nMicrosites) > 1, then the number of microsites available for establishment at time t is nMicrosites[min(t,length(nMicrosites))].
}
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density
    function); other option is "cumul" (which uses the cumulative density function).
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }

}
\value{
\item{rc }{matrix of the numbers of individuals in each size and seed
  class (row) over time (columns).}
\item{IPM.here}{IPM constructed corresponding to pop structure and
  covariates at tMax.}
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{stochGrowthRateManyCov}}
}
\examples{

### NOT RUN ####

# Generate relevant data, build objects
#dff <- generateDataStoch()
#print(head(dff))
#gr1 <- makeGrowthObjManyCov(dff, explanatoryVariables = "size+size2+covariate1")
#sv1 <- makeSurvObjManyCov(dff, explanatoryVariables = "size+size2+covariate2")
#fv1 <- makeFecObj(dff, explanatoryVariables = "size+size2",Transform="log")

# Generate time series of covariates for which growth rate is required
# here set to be seasonal environment. NOTE that the order of covariates
# is the same as in dff
#tVals <- seq(1, 20, by = 1/12)
#covTest <- (1 + 0.5*sin(2*pi*tVals))
#covMatTest <- data.frame(covariate1 = rnorm(length(covTest), covTest, 0.5) - 1, 
#covariate2 = rnorm(length(covTest),covTest,0.5) - 1, 
#covariate3 = rnorm(length(covTest), covTest,0.5) - 1)

#st <- trackPopStructManyCov(covariate = covMatTest, nRunIn = 2*10, 
#tMax = length(tVals), growthObj = gr1, survObj = sv1, fecObj = fv1, 
#nBigMatrix = 100, minSize = 2*min(dff$size, na.rm = TRUE), 
#maxSize = 1.5*max(dff$size, na.rm = TRUE), nMicrosites = 10)

# Function to plot output
#plotResultsStochStruct(tVals = tVals, st = st, covTest = covTest)

}
