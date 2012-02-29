\name{growthModelComp}
\alias{growthModelComp}
\alias{survModelComp}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compares growth and survival objects built from different covariate sets.
}
\description{
Function compares model fits for growth and survival objects built with different linear combinations of covariates.  Growth can have multiple response forms. 
Returns a list containing a summary table of covariates and scores, and a nother list containing all of the growth (or survival) objects used in the comparison.}
\usage{
growthModelComp <- function(dataf, expVars = c("1", "size", "size + size2"), regressionType = "constantVar", respType = "sizeNext", testType = "AIC", makePlot = FALSE, mainTitle = "")
		
survModelComp <- function(dataf, expVars = c("1", "size", "size + size2"), regressionType = "constantVar", testType = "AIC", makePlot = FALSE, mainTitle = "")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{expVars}{
    vector, list of covariates.  Defaults to \code{c("1",  "size", "size + size2")}.
  }
  \item{regressionType}{
    character string identifying whether the type of regression run will have constant or changing variance.  Defaults to \code{constantVar}.
  }
  \item{respType}{
    character string identifying the response variable for the growthModelComp.  wither \code{"sizeNext"}, \code{"logincr"} or \code{"incr"}. Defaults to \code{sizeNext}. 
  }
  \item{testType}{
    character string identifying the metric used to compare models.  Can be any string that uses \code{loglike} from the lm or glm object.  For example \code{"AIC"} or \code{"BIC"}.  Defaults to \code{"AIC"}.
  }
  \item{envMatrix}{
    envMatrix object defining transitions between environmental states
    for each size.
  }
  \item{fecObj}{
    fecundity object
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
  an object of class IPM.matrix with dimensions
  nBigMatrix*nEnvClass, or if discrete transitions exist
  (nBigMatrix+n.disc)*nEnvClass

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
\code{\link{create.compound.Fmatrix}},\code{\link{create.IPM.Tmatrix}}
}
\examples{
# Data with only continuous stage and two habitats
dff <- generateData()
Fmatrix <-create.compound.Fmatrix(nBigMatrix = 20, minSize = min(dff$size, na.rm = TRUE), maxSize = max(dff$size, na.rm=TRUE), envMatrix = makeEnvObj(dff), fecObj = makeFecObj(dff, explanatoryVariables = "size+size2+covariate", Transform = "log"))

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
Fmatrix <- create.compound.Fmatrix(minSize = min(dff$size, na.rm=TRUE),
maxSize = max(dff$size, na.rm = TRUE), envMatrix = makeEnvObj(dff), fecObj = fv1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
