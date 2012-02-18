\name{TrackPopStructManyCov}
\alias{TrackPopStructManyCov}
\alias{plotResultsStochStruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to track population structure in a stochastic environment
}
\description{
  Provided a set of time-varying covariates, projects the population
  forward, and stores state of the population. 
}
\usage{
TrackPopStructManyCov(covariate, n.runin, Tmax, growthObj, survObj, fecObj, 
n.big.matrix, minSize, maxSize, n.microsites, integrate.type, correction)
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
  \item{minSize}{
    numeric, minimum size in the IPM
}
  \item{maxSize}{
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

}
\value{
\item{rc }{matrix of the numbers of individuals in each size and seed
  class (row) over time (columns)}
\item{IPM.here}{IPM constructed corresponding to pop structure and
  covariates at Tmax}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{StochGrowthRateManyCov}}
}
\examples{

#generate relevant data, build objects
dff <- generateDataStoch()
print(head(dff))
gr1 <- makeGrowthObjManyCov(dff,explanatoryVariables="size+size2+covariate1")
sv1 <- makeSurvObjManyCov(dff,explanatoryVariables="size+size2+covariate2")
fv1 <- makeFecObj(dff,explanatoryVariables="size+size2")

#generate time-series of covariates for which growth rate is required
# here set to be seasonal environment. NOTE that the order of covariates
# is the same as in dff. 
tvals <- seq(1,20,by=1/12)
covtest <- (1+0.5*sin(2*pi*tvals))
covmattest <- data.frame(covariate1=rnorm(length(covtest),covtest,0.5)-1,
                    covariate2=rnorm(length(covtest),covtest,0.5)-1,
                    covariate3=rnorm(length(covtest),covtest,0.5)-1)


st<-TrackPopStructManyCov(covariate=covmattest,n.runin=2*10,Tmax=length(tvals),
                   growthObj=gr1,survObj=sv1,fecObj=fv1,
                   n.big.matrix=100,minSize=2*min(dff$size,na.rm=TRUE),
                   maxSize=1.5*max(dff$size, na.rm = TRUE),n.microsites=10)


#function to plot output
plotResultsStochStruct(tvals=tvals,st=st,covtest=covtest)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
