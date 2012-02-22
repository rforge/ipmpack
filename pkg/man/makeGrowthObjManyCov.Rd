\name{makeGrowthObjManyCov}
\alias{makeGrowthObjManyCov}

\title{Functions to build growth objects for varying environments
  (discretely, continuously, with potential for several covariates)}
\description{
  A function that fits regressions defining growth (following next size,
  size increment, or log size increment) with covariates and from these build growth objects for which methods to build an IPM
  object are defined.
}
\usage{
makeGrowthObjManyCov(dataf, explanatoryVariables="size+size2+covariate1", responseType="sizeNext",regType="constantVar")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizeNext' (`size' is size at t, `sizeNext' is size at t+1);
    facultatively, dataf may include a number of covariates for which
    column names must follow the sequence `covariate1', `covariate2', ...
    For models fitting growth increment, 'incr' or 'logincr' may be directly provided as
    a column in the dataframe, otherwise it is calculated as
    dataf$sizeNext-dataf$size; or log(dataf$sizeNext-dataf$size). }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate1', `covariate2' through to whatever the desired number of covariates is.}
  \item{responseType}{possible values include `sizeNext' (size in one
    discrete time step) and `incr' and `logincr'}
  \item{regType}{possible values inclue `constantVar' (defines wether a
  normal regression is fitted) or `declineVar' (declining variance with
  growth e.g., fitted using a gls)}
  
  }
\details{
 
}
\value{
An object of class growthObjManyCov, or growthObjManyCov.incr, or growthObjManyCov.logincr;
or growthObjManyCov.declinevar, or growthObjManyCov.incr.declinevar, or growthObjManyCov.logincr.declinevar.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
See manual for details on building case-specific growth objects. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeGrowthObjSubsetCovariates}}, \code{\link{makeSurvObj}}, \code{\link{makeFecObj}}
}
\examples{
#generate data
dff<-generateDataStoch()
#make simple linear regression growth object relating size to size at t+1
gr1 <- makeGrowthObjManyCov(dataf=dff,
explanatoryVariables="size+covariate1+covariate3", responseType="sizeNext")
#same but relating size to incr with slightly different covariate set
gr1 <- makeGrowthObjManyCov(dataf=dff, explanatoryVariables="size+covariate1", responseType="incr")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
