\name{meanLifeExpect}
\alias{meanLifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Calculates mean life expectancy.
}
\description{
  Provided a T matrix, which defines survival transitions across stages, this function outputs a vector defining life expectancy in units of the time step used (see convertIncrement()), for each of the size bins.
}
\usage{
meanLifeExpect(IPMmatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPMmatrix}{
    an IPMmatrix object defining survival transitions.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A vector of life expectancies each corresponding to a value of the size bins defined by Tmatrix@meshpoints.
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p118-120.

Cochran & Ellner. 1992. Simple methods for calculating age-based life history parameters for stage-structured populations. Ecological Monographs 62, p345-364.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LifeExpect}}, \code{\link{createIPMTmatrix}} 
}
\examples{
# With a single continuous state variable (e.g. size):
dff <- generateData()
Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), growObj=makeGrowthObj(dff),
survObj = makeSurvObj(dff))
meanLifeExpect(Tmatrix)

Tmatrix <- createIPMTmatrix(minSize = min(dff$size, na.rm = TRUE),
maxSize = max(dff$size, na.rm = TRUE), growObj = makeGrowthObj(dff),
survObj = makeSurvObj(dff))

plot(meanLifeExpect(Tmatrix), ylab = "Mean life expectancy", xlab = "Continuous (e.g. size) stage", type = "l", col="dark gray", ylim = c(0,max(meanLifeExpect(Tmatrix))))

# With continuous (e.g. size) and discrete (e.g. seedbank) stages:
dff <- generateDataDiscrete()
dff$covariate <- sample(1:3, size = nrow(dff), replace = TRUE)
dff$covariateNext <- sample(1:3, size = nrow(dff), replace = TRUE)
discM <- makeDiscreteTrans(dff)
Tmatrix <- createCompoundTmatrix(minSize = min(dff$size, na.rm = TRUE), maxSize = max(dff$size, na.rm = TRUE), envMatrix = makeEnvObj(dff),
growObj = makeGrowthObj(dff, explanatoryVariables = "size+size2+covariate"),
survObj = makeSurvObj(dff, explanatoryVariables = "size+size2+covariate"), discreteTrans = discM)
meanLifeExpect(Tmatrix)

plot(meanLifeExpect(Tmatrix), 
ylab = "Mean life expectancy", 
xlab = "Continuous (e.g. size) and discrete (e.g. seedbank) stages", type = "l", col = "dark gray", ylim = c(0,max(meanLifeExpect(Tmatrix))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
