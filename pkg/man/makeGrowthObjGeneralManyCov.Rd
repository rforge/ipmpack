\name{makeGrowthObjGeneralManyCov}
\alias{makeGrowthObjGeneralManyCov}

\title{Functions to build growth objects for varying environments
  (discretely, continuously, with potential for several covariates)}
\description{
  A function that fits regressions defining growth (following next size,
  size increment, or log size increment) with covariates and from these build growth objects for which methods to build an IPM
  object are defined.
}
\usage{
makeGrowthObjGeneralManyCov(dataf, explanatoryVariables="size+covariate1+covariate2", responseType="sizenext")
makeGrowthObjGeneralManyCov(dataf, explanatoryVariables="size+size2+covariate3", responseType="incr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizenext' (`size' is size at t, `sizenext' is size at t+1);
    facultatively, dataf may include a number of covariates for which
    column names must follow the sequence `covariate1', `covariate2', ...
    For models fitting growth increment, 'incr' or 'logincr' may be directly provided as
    a column in the dataframe, otherwise it is calculated as
    dataf$sizenext-dataf$size; or log(dataf$sizenext-dataf$size). }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) in classical R
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate1', `covariate2' through to whatever the desired number of covariates is.}
  \item{responseType}{possible values include `sizenext' (size in one
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
%%  ~~who you are~~
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
dff<-generateDataArabi()
#make simple linear regression growth object relating size to size at t+1
gr1 <- makeGrowthObjGeneralManyCov(dataf=dff,
explanatoryVariables="size+covariate1+covariate3", responseType="sizenext")
#same but relating size to incr
gr1 <- makeGrowthObjGeneralManyCov(dataf=dff, explanatoryVariables="size+covariate3", responseType="incr")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
