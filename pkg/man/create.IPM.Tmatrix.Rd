\name{create.IPM.Tmatrix}
\alias{create.IPM.Tmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds T matrices.
}
\description{
  Uses growth and survival objects to construct a matrix defining probabilities for transitions between sizes due to growth and survival. Extensions for transition to discrete classes are possible.
}
\usage{
create.IPM.Tmatrix(n.env.class = 1, nBigMatrix = 50,
minSize = -1, maxSize = 50, chosen.cov = 1, growObj, survObj, discreteTrans, 
integrateType, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.env.class}{
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
    setting values where multiple covarites are modeled.
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
  The number of bins (nBigMatrix) is combined with the minimum and
  maximum size to define the meshpoints of the IPM. Bins should be
  sufficient and the size range should encompass the size range of the
  data. If a "discreteTrans" exists, then discrete stages will be added
  to the IPM structure. If multiple discrete covariate levels are available,
  chosen.cov identifies the covariate value for which an IPM is
  required; if a series of covariates are being modeled, chosen.cov is a
  vector of these covariates, and growth will reflect these values.  
}
\value{
  an object of class IPM.matrix with dimensions nBigMatrix*nBigMatrix, or
  nbig.matrix+nrow(discreteTrans).
}
\references{
Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81, p694-708.

Ellner & Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167, p410-428.

For information on T matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.
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
\code{\link{create.IPM.Fmatrix}}, \code{\link{diagnosticsTmatrix}}, \code{\link{makediscreteTrans}}  
}
\examples{
dff <- generateData()
Tmatrix <- create.IPM.Tmatrix(minSize = min(dff$size, na.rm=TRUE), maxSize=max(dff$size, na.rm=TRUE), growObj = makeGrowthObj(dff), survObj = makeSurvObj(dff))

slotNames(Tmatrix)

image(Tmatrix@meshpoints, Tmatrix@meshpoints, t(Tmatrix), xlab="Continuous (e.g. size) stage at t", ylab = "Continuous (e.g. size) stage at t+1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
