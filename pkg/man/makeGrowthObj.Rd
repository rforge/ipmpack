\name{makeGrowthObj}
\alias{makeGrowthObj}

\title{Functions to build growth objects for constant environments, or
  discretely varying environments}
\description{
  A function that fits regressions that define growth (following next size,
  size increment, or log size increment) and from these build growth objects for which methods to build an IPM
  object are defined.
}
\usage{
makeGrowthObj(dataf, Formula = sizeNext ~ size, 
	regType = "constantVar", Family = "gaussian", link = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataf}{a dataframe with columns `size' and
    `sizeNext'(`size' is size at t, `sizeNext' is size at t+1);
    facultatively, dataf may include `covariate'
    and `covariatel' for a single discrete covariate, indicating values at t,
    and at t+1, respectively; these must take values of sequential integers, starting at
    `1'. For models fitting growth increment, 'incr' or 'logincr' may be directly provided as
    a column in the dataframe, otherwise they are calculated as
    dataf$sizeNext - dataf$size or log(dataf$sizeNext - dataf$size),
  respectively. }
  \item{Formula}{a formula describing the desired
    explanatory variables (interactions, etc) according to the R notation for \code{formula}.
    style, i.e. separated by `+', `*', `:' and response variable. Possible covariates include `size', 'size2' (size^2), `size3'
    (size^3),`logsize' (log(size)), `logsize2' (log(size)^2), and
    `covariate' (if this name is used, the assumption is made that this is a discrete covariate from which compound matrices may be constructed); or any other covariate available in dataf.}
  \item{regType}{possible values include `constantVar' or `declineVar'}
  \item{Family}{possible values include `gaussian', `poisson', `negbin'}
   \item{link}{defaults to NULL, currently only relevant with Family="negbin", and only "log" and "identity" are permitted}
  
  }
\value{
An object of class growthObj, or growthObjPois, growthObjIncr, or growthObjLogIncr;
or growthObjDeclineVar, or growthObjIncrDeclineVar, or growthObjLogIncrDeclineVar. 
These are S4 objects which contain the slots: 
   \item{fit}{an object of class \code{lm} or \code{glm} or \code{gls} that can be used with predict in the survival methods}
Slots can be listed by using slotNames(growthObj)   

}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
See manual for details on building case-specific growth objects; which may contain different slots, etc.  
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 \code{\link{makeSurvObj}}, \code{\link{makeFecObj}}
}
\examples{
#generate data
dff <- generateData()
#make simple linear regression growth object relating size to size at t+1
(gr1 <- makeGrowthObj(dataf = dff, Formula = sizeNext~size))
#same but relating size to incr
(gr1 <- makeGrowthObj(dataf = dff, Formula = incr ~ size))

#assess fit 
picGrow(dff,gr1)

#same but relating size to incr and discrete covariate
(gr1 <- makeGrowthObj(dataf = dff, Formula = incr ~ size + covariate, 
	regType = "declineVar"))
#with declining increment
(gr1 <- makeGrowthObj(dataf = dff, Formula = incr ~ size + covariate, 
	regType = "declineVar"))


}
