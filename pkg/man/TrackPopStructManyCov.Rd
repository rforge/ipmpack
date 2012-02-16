\name{TrackPopStructManyCov}
\alias{TrackPopStructManyCov}
\alias{TrackPopStructManyCovSeedBank}
\alias{plotResultsStochStruct}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to track population structure in a stochastic environment
}
\description{
  Provided a set of time-varying covariates, projects the population
  forward, and stores state of the population, for populations structured
  by size alone, or size and a seedbank. 
}
\usage{
TrackPopStructManyCov(covariate, n.runin, Tmax, growthObj, survObj, fecObj, n.big.matrix, minsize, maxsize, n.microsites)
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
}
\details{
If TrackPopStructManyCovSeedBank is used, the IPM matrix constructed
includes a seedbank; otherwise the population is just structured by
size. Forms of density dependence beyond density dependence in seedling
establishment not yet defined. 
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
%%  ~~who you are~~
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
dff <- generateDataArabi()
print(head(dff))
gr1 <- makeGrowthObjSubsetCovariates(dff,1:3)
sv1 <- makeSurvObjSubsetCovariates(dff,1:3)
fv1 <- makeFecObjBasicSubsetCovariates(dff,1:3,1:3)
fv1@seed.mortality <- 0.5
fv1@germ.params <- c(-2,8)

#generate time-series of covariates for which growth rate is required
# here set to be seasonal environment. NOTE that the order of covariates
# is the same as in dff. 
tvals <- seq(1,20,by=1/12)
covtest <- (1+0.5*sin(2*pi*tvals))
covmattest <- cbind(rnorm(length(covtest),covtest,0.5)-1,
                    rnorm(length(covtest),covtest,0.5)-1,
                    rnorm(length(covtest),covtest,0.5)-1)

st<-TrackPopStructManyCovSeedBank(covariate=covmattest,n.runin=12*10,Tmax=length(tvals),
                   growthObj=gr1,survObj=sv1,fecObj=fv1,
                   n.big.matrix=100,minsize=2*min(dff$size),maxsize=1.5*minsize=max(dff$size, na.rm = TRUE),n.microsites=10)


#function to plot output
plotResultsStochStruct(tvals=tvals,st=st,covtest=covtest)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
