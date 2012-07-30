\name{coerceGrowthObj}
\alias{coerceGrowthObj}
\alias{coerceSurvivalObj}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to coerce growth or survival objects, i.e., impose user-defined parameters 
}
\description{
supplied with a growth and survival object, over-writes coefficients, and for growth, the sd of growth 
}
\usage{

coerceGrowthObj(growthObj,coeff,sd)
coerceSurvObj(growthObj,coeff,sd)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{growthObj}{an object of class growthObj}
    \item{survObj}{an object of class survObj}
    \item{coeff}{a numeric vector}
    \item{sd}{a numeric vector of length 1}
}
  
\details{
  These functions can be used to impose coefficients and sd on growth and survival objects where direct fitting is not desired
}

\value{
  an object of class growthObj / survObj
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeSurvObj}}, \code{\link{makeGrowthObj}}

}
\examples{
dff<-generateData()

gr1 <- makeGrowthObj(dataf=dff,
Formula=sizeNext~size,regType="constantVar")

#halve the slope
gr2 <- coerceGrowthObj(gr1,coeff=c(gr1@fit$coefficients[1],gr1@fit$coefficients[2]*0.5),sd=gr1@sd)

par(mfrow=c(1,2),pty="s")
picGrow(dff,gr1)
picGrow(dff,gr2)


}
