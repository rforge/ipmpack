\name{createGrowthObj}
\alias{createGrowthObj}
\alias{createSurvObj}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to create growth or survival objects with user-defined parameters 
}
\description{
supplied with a Formula, coefficients, and for growth a standard
deviation, creates a growth or survival object 
}
\usage{

createGrowthObj(Formula=sizeNext~size, coeff=c(1,1), sd=1)
createSurvObj(Formula=surv~size, coeff=c(1,1))

}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Formula}{an object of class formula}
    \item{coeff}{a numeric vector of length the number of coefficients
      supplied by Formula}
    \item{sd}{a numeric vector of length 1}
}
  
\details{
  These functions create growth and survival objects where direct
  fitting is not desired; note that currently they are not capable of
  managing Formula that encompass factors. 
}

\value{
  an object of class growthObj / survObj
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeSurvObj}}, \code{\link{makeGrowthObj}}, \code{\link{coerceGrowthObj}},\code{\link{coerceSurvObj}}

}
\examples{
dff<-generateData()

#for growth
gr1 <- createGrowthObj(Formula=sizeNext~size,coeff=c(0.1,0.9),sd=1)

#halve the slope
gr2 <- createGrowthObj(Formula=sizeNext~size,coeff=c(0.1,0.45),sd=1)


#for survival
sv1 <- createSurvObj(Formula=surv~size,coeff=c(-10,3))



}
