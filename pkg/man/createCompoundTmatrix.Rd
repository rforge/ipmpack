\name{createCompoundTmatrix}
\alias{createCompoundTmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds a compound T matrix.
}
\description{
  Uses growth, survival, discreteTrans, and environmental transition objects to construct a matrix defining probabilities for transitions between continuous stages (e.g. size) due to growth and survival given both size and environmental state and discrete stages. 
}
\usage{
createCompoundTmatrix(nEnvClass = 2, nBigMatrix = 50,
minSize = -1, maxSize = 50, envMatrix, growObj, survObj,
discreteTrans = 1, integrateType = "midpoint",correction = "none")
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
  \item{growObj}{
    growth object.
  }
  \item{survObj}{
    survival object.
  }
  \item{discreteTrans}{
    object of class discreteTrans, or numeric.
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
This structure can also be used to define size x age IPMs, where the transition 
between ages is reflected by a similar matrix.
}
\value{
  an object of class IPMmatrix with dimensions
  nBigMatrix * nEnvClass, or if discrete transitions exist
  (nBigMatrix + nDisc) * nEnvClass
}
\references{
For information on T matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.

For habitat x stage modeling: Horvitz and Schemske. 1995. Spatiotemporal variation in demographic transitions of a tropical understory herb. Projection matrix analysis. Ecological Monographs 65, p155-192.

For age x size modeling: Garcia, Dahlgren and Ehrlen. 2011. No evidence of senescence in a 300-year-old mountain herb. Journal of Ecology 99, p1424Ð1430.

For general information:

Easterling, Ellner and Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81, p694-708.

Ellner and Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167, p410-428.
%% ~put references to the literature/web site here ~
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
Tmatrix <- createCompoundTmatrix(minSize = min(dff$size,na.rm = TRUE), maxSize = max(dff$size,na.rm = TRUE), envMatrix = makeEnvObj(dff), growObj = makeGrowthObj(dff, explanatoryVariables = "size+size2+covariate"), survObj = makeSurvObj(dff, explanatoryVariables = "size+size2+covariate"))

image(1:nrow(Tmatrix), 1:ncol(Tmatrix), t(log(Tmatrix)), xlab = "Continuous stage (e.g. size) at t", ylab = "Continuous stage (e.g. size) at t+1", axes = FALSE)
axis(1, at = 1:nrow(Tmatrix), lab = round(rep(Tmatrix@meshpoints,Tmatrix@nEnvClass), 2))
axis(2, at = 1:nrow(Tmatrix), lab = round(rep(Tmatrix@meshpoints,Tmatrix@nEnvClass), 2))
abline(h = length(Tmatrix@meshpoints) * (1:Tmatrix@nEnvClass))
abline(v = length(Tmatrix@meshpoints) * (1:Tmatrix@nEnvClass))

# Data with continuous and discrete stages
dff <- generateDataDiscrete()
dff$covariate <- sample(1:3, size = nrow(dff), replace = TRUE)
dff$covariateNext <- sample(1:3, size = nrow(dff), replace = TRUE)
discM <- makeDiscreteTrans(dff)
Tmatrix <- createCompoundTmatrix(minSize = min(dff$size, na.rm = TRUE), maxSize = max(dff$size, na.rm = TRUE), envMatrix = makeEnvObj(dff), growObj = makeGrowthObj(dff, explanatoryVariables = "size+size2+covariate"), survObj = makeSurvObj(dff, explanatoryVariables = "size+size2+covariate"), discreteTrans = discM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
