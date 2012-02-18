\name{makeGrowthObjGeneral}
\alias{makeGrowthObjGeneral}
\alias{makeGrowthObjHossfeld}

\alias{picGrow}

\title{Functions to build growth objects for constant environments, or
  discretely varying environments}
\description{
  A function that fits regressions defining growth (following next size,
  size increment, or log size increment) and from these build growth objects for which methods to build an IPM
  object are defined.
}
\usage{
makeGrowthObjGeneral(dataf, explanatoryVariables="size", responseType="sizenext",regType="constantVar")
makeGrowthObjGeneral(dataf, explanatoryVariables="size+size2", responseType="incr",regType="constantVar")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizenext'(`size' is size at t, `sizenext' is size at t+1);
    facultatively, dataf may include `covariate'
    and `covariatel' for a single discrete covariate, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. For models fitting growth increment, 'incr' or 'logincr' may be directly provided as
    a column in the dataframe, otherwise they are calculated as
    dataf$sizenext-dataf$size or log(dataf$sizenext-dataf$size),
  respectively. }
  \item{explanatoryVariables}{a formula describing the desired
    explanatory variables (interactions, etc) according to the R notation for \code{formula}.
    style, i.e. separated by `+', `*', `:'. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate'.}
  \item{responseType}{possible values include `sizenext' (size in one
    discrete time step) and `incr' and `logincr'}
  \item{regType}{possible values include `constantVar' or `declineVar'}
  
  }
\details{
 
}
\value{
An object of class growthObj, or growthObj.incr, or growthObj.logincr;
or growthObj.declinevar, or growthObj.incr.declinevar, or growthObj.logincr.declinevar.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
See manual for details on building case-specific growth objects. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeGrowthObjGeneralManyCov}}, \code{\link{makeSurvObjGeneral}}, \code{\link{makeFecObj}}
}
\examples{
#generate data
dff<-generateData()
#make simple linear regression growth object relating size to size at t+1
(gr1 <- makeGrowthObjGeneral(dataf=dff, explanatoryVariables="size", responseType="sizenext"))
#same but relating size to incr
(gr1 <- makeGrowthObjGeneral(dataf=dff, explanatoryVariables="size", responseType="incr"))
#same but relating size to incr and discrete covariate
(gr1 <- makeGrowthObjGeneral(dataf=dff, explanatoryVariables="size+covariate", responseType="incr"))


#assess fit 
picGrow(dff,gr1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
