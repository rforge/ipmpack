\name{StochGrowthRateManyCov}
\alias{StochGrowthRateManyCov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for estimating lambdas or invasion rate with many varying covariates
}
\description{
Iterates a population through a time-series of covariates according to
growth, survival and fertility objects, and
calculates the stochastic population rate of increase if no
density-dependence is specified, or the rate of invasion if density
dependence is specified.
}
\usage{
StochGrowthRateManyCov(covariate, n.runin, Tmax,
growthObj, survObj, fecObj, n.big.matrix, minsize, maxsize, 
n.microsites, integrate.type,correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{covariate}{
    matrix with Tmax rows, and as many columns as there are relevant covariates
}
  \item{n.runin}{
    numeric, number of initial samples to discard
}
  \item{Tmax}{
    numeric, total number of time-steps to run (same as ncol(covariate))
}
  \item{growthObj}{
    a growth object, defined to correspond to covariate definition
    (indexing used to make the growth object must match up)
}
  \item{survObj}{
    a survival object, defined to correspond to covariate definition in covariate
}
  \item{fecObj}{
    a fecundity object, defined to correspond to covariate definition in covariate
}
  \item{n.big.matrix}{
    numeric, number of size bins in the IPM
}
  \item{minsize}{
    numeric, minimum size in the IPM
}
  \item{maxsize}{
    numeric, maximum size in the IPM
}
  \item{n.microsites}{
    vector, if sum(n.microsites)> 0 then density dependence is assumed to operate on
    seedling establishment, and if length(n.microsites)>1, then the
    number of microsites available for establishment at time t is n.microsites[min(t,length(n.microsites))] 
}
  \item{integrate.type}{
    integration type, defaults to "midpoint" (which uses prob. density
    function); other option is "cumul" (which uses the cumul. density function)
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }
}
\details{
Forms of density dependence beyond density dependence in seedling
establishment not yet defined. 
}
\value{
numeric, converging on log lambda_s (no density dependence) or invasion
rate (density dependence) for large enough Tmax, and if
covariate distribution is stationary. 
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{TrackPopStructManyCov}},\code{\link{TrackPopStructManyCov}}
}
\examples{

#generate relevant data, build objects
dff <- generateDataStoch()
print(head(dff))
gr1 <- makeGrowthObjGeneralManyCov(dff,explanatoryVariables="size+size2+covariate1")
sv1 <- makeSurvObjGeneralManyCov(dff,explanatoryVariables="size+size2+covariate2")
fv1 <- makeFecObjGeneralManyCov(dff,explanatoryVariables="size+size2")

#generate time-series of covariates for which population growth rate is required
# here set to be seasonal environment. NAMES of covariates must be same as in dff
tvals <- seq(1,20,by=1/12)
covtest <- (1+0.5*sin(2*pi*tvals))
covmattest <- data.frame(covariate1=rnorm(length(covtest),covtest,0.5)-1,
                    covariate2=rnorm(length(covtest),covtest,0.5)-1,
                    covariate3=rnorm(length(covtest),covtest,0.5)-1)

#calculate
r<-StochGrowthRateManyCov(covariate=covmattest,n.runin=5*10,Tmax=length(tvals),
                   growthObj=gr1,survObj=sv1,fecObj=fv1,
                   n.big.matrix=100,minsize=2*min(dff$size,na.rm=TRUE),
                   maxsize=1.5*max(dff$size, na.rm = TRUE),n.microsites=10)
r

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
