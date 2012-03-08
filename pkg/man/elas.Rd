\name{elas}
\alias{elas}
\alias{sens}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates matrix element sensitivities and elasticities.
}
\description{
  Estimates sensitivities and elasticities of each element of a discretized IPM.
}
\usage{
elas(A)
sens(A)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{A}{
    a matrix defining all transitions (growth, survival, fecundity) between sizes/stages.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a matrix.
}
\references{
Caswell, 2001, Matrix Population Models: construction, analysis,
interpretation. 2nd ed. Sinauer. p206-256.

de Kroon, Plaisier, van Groenendael & Caswell. 1986. Elasticity: the relative contribution of demographic parameters to population growth rate. Ecology 67, p1427-1431.

de Kroon, van Groenendael & Ehrlen. 2000. Elasticities: a review of methods and model limitations. Ecology 81, p607-618.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
Modified following code developed by Mark Rees & Dylan Childs
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}
}
\examples{
dff <- generateData()
Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm=TRUE),
maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff),
survObj = makeSurvObj(dff))

Fmatrix <- createIPMFmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), fecObj = makeFecObj(dff))

IPM <- Tmatrix + Fmatrix

par(mfrow = c(1, 2))

senst <- sens(IPM)
image(Tmatrix@meshpoints, Tmatrix@meshpoints,t(senst), main = "Sensitivity", xlab = "Continuous (e.g. size) stage in t", ylab = "Continuous (e.g. size) stage in t+1")

elast <- elas(IPM)
image(Tmatrix@meshpoints, Tmatrix@meshpoints, t(elast), main = "Elasticity", xlab = "Continuous (e.g. size) stage in t", ylab = "Continuous (e.g. size) stage in t+1")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
