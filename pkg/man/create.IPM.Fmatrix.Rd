\name{create.IPM.Fmatrix}
\alias{create.IPM.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here. 
\title{
  Builds F matrices.
}
\description{
Uses fecundity objects to construct a matrix defining per-capita contribution to recruitment stages (e.g., propagules [seed, spore], seedlings, calves) by reproductive stages due to sexual reproduction.
}
\usage{
create.IPM.Fmatrix(fecObj, nEnvClass = 1, nBigMatrix = 50, minSize = -1, maxSize = 50, chosen.cov = 1, integrateType="midpoint", correction="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fecObj}{
    fecundity object.
  }
  \item{nEnvClass}{
    numeric, number of environmental classes, always = 1 for
    non-compound matrices.
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
  \item{chosen.cov}{
    numeric indicating level of the discrete covariate, or vector
    setting values where multiple covariates are modeled.
  }
  \item{integrateType}{
    integration type, defaults to "midpoint" (which uses probability density
    function); other option is "cumul" (which uses the cumulative density function)
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }
 
}
\details{
%%
}
\value{
  an object of class IPM.matrix of dimensions nBigMatrix or nBigMatrix+ndiscrete classes (defined by fecObj@OffspringSplitter-1).
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
\code{\link{create.IPM.Tmatrix}}
}
\examples{
# Data with only continuous stage and one habitat
dff <- generateData()
fv1 <- makeFecObj(dff, explanatoryVariables = "size", Transform = "log")
Fmatrix <- create.IPM.Fmatrix(fecObj = fv1, nBigMatrix = 20, minSize = min(dff$size, na.rm = TRUE), maxSize = max(dff$size, na.rm = TRUE))

slotNames(Fmatrix)

image(Fmatrix@meshpoints, Fmatrix@meshpoints, t(Fmatrix), xlab = "Continuous (e.g. size) stage at t", ylab = "Continous (e.g. size) stage at t+1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
