\name{createCompoundFmatrix}
\alias{createCompoundFmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds a compound F matrix.
}
\description{
  Uses fecundity object, and environmental transition objects to construct a matrix defining probabilities for transitions between sizes due to fecundity given both a continuous state (e.g. size) and environmental state, as well as a discrete stage if necessary (e.g. seedbank). Currently only pre-census fertility can be handled.  
}
\usage{
createCompoundFmatrix(nEnvClass = 2, nBigMatrix = 50,
minSize = -1, maxSize = 50, envMatrix, fecObj, integrateType, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nEnvClass}{
    numeric, number of environmental classes, defaults to 2.
  }
  \item{nBigMatrix}{
    numeric, number of size bins in the T matrix, defaults to 50.
  }
  \item{minSize}{
    numeric, minimum size of the T matrix, defaults to -1.
  }
  \item{maxSize}{
    numeric, maximum size of the T matrix, defaults to 50.
  }
  \item{envMatrix}{
    envMatrix object defining transitions between environmental states
    for each size.
  }
  \item{fecObj}{
    fecundity object.
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  an object of class IPMmatrix with dimensions
  nBigMatrix*nEnvClass, or if discrete transitions exist
  (nBigMatrix+nDisc)*nEnvClass

}
\references{
For information on F matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.

For midpoint: Zuidema, Jongejans, Chien, During & Schieving. Integral projection models for trees: a new parameterization method and a validation of model output. Journal of Ecology 98, p345-355.

For multiple-vital rate integration on fecundity: Yang, Jongejans, Yang & Bishop. 2011. The effect of consumers and mutualists of Vaccinum membranaceum at Mount St. Helens: dependence on successional context. PLoS One 10, p1-11.

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createCompoundFmatrix}},\code{\link{createIPMTmatrix}}
}
\examples{
# Data with only continuous stage and two habitats
dff <- generateData()
Fmatrix <-createCompoundFmatrix(nBigMatrix = 20, minSize = min(dff$size, na.rm = TRUE), maxSize = max(dff$size, na.rm=TRUE), envMatrix = makeEnvObj(dff), fecObj = makeFecObj(dff, explanatoryVariables = "size+size2+covariate", Transform = "log"))

image(1:nrow(Fmatrix), 1:ncol(Fmatrix), t(log(Fmatrix)), xlab="Continuous state (e.g. size) at t", ylab="Continuous state (e.g. size) at t+1", axes = FALSE)
axis(1, at = 1:nrow(Fmatrix), lab = round(rep(Fmatrix@meshpoints, Fmatrix@nEnvClass), 2))
axis(2,at = 1:nrow(Fmatrix), lab = round(rep(Fmatrix@meshpoints,Fmatrix@nEnvClass), 2))
abline(h = length(Fmatrix@meshpoints)*(1:Fmatrix@nEnvClass))
abline(v = length(Fmatrix@meshpoints)*(1:Fmatrix@nEnvClass))

# Data with continuous and discrete stages
dff <- generateDataDiscrete()
dff$covariate <- sample(1:3, size = nrow(dff), replace = TRUE)
dff$covariateNext <- sample(1:3, size = nrow(dff), replace = TRUE)
fv1 <- makeFecObj(dff, explanatoryVariables = "size", Transform = "log", offspringSplitter=data.frame(continuous = 0.9, dormant = 0.1))
Fmatrix <- createCompoundFmatrix(minSize = min(dff$size, na.rm=TRUE),
maxSize = max(dff$size, na.rm = TRUE), envMatrix = makeEnvObj(dff), fecObj = fv1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
