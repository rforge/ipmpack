\name{lifeExpect}
\alias{lifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates life expectancy in a discretely stochastic environment.
}
\description{
  Provided a compound T matrix that describes survival transitions across continuous (e.g. size) stage and across discrete environmental states, this function provides a matrix defining life expectancy in units of the time step used in the census, for each of the size bins in each of the environments (rows). 
}
\usage{
lifeExpect(IPMmatrix, envMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPMmatrix}{
    an IPMmatrix object describing growth and survival transitions
    across size and environment.
}
  \item{envMatrix}{
    an envMatrix object describing transition between environments.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  matrix, each row corresponds to life expectancy in a different environments, columns correspond to sizes provided in Tmatrix@meshpoints.
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p110-132.

Cochran & Ellner. 1995. Simple methods for calculating age-based life history parameters for stage-structured populations. Ecological Monographs 62, p345-364.

Tuljapurkar & Horvitz, 2006. From stage to age in variable environments. Life expectancy and survivorship. Ecology 87, p1497-1509.
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{meanLifeExpect}} \code{\link{createCompoundTmatrix}}
}
\examples{

# With a compound matrix (discrete stages not yet defined)
dff <- generateData()
Tmatrix <- createCompoundTmatrix(minSize = min(dff$size, na.rm = TRUE), 
	maxSize = max(dff$size, na.rm = TRUE), 
	growObj = makeGrowthObj(dff, explanatoryVariables = "size+size2+covariate"), 
		survObj = makeSurvObj(dff), envMatrix = makeEnvObj(dff))
le <- lifeExpect(Tmatrix, makeEnvObj(dff))

plot(Tmatrix@meshpoints, le[1, ], ylab = "Life expectancy", 
xlab = "Continuous (e.g. size) stage in t", 
	ylim = c(0,max(le)), type = "l", col = "blue")
points(Tmatrix@meshpoints, le[2, ], type = "l", col = "dark green", lty = 2)
legend(min(dff$size, na.rm = TRUE), 
min(dff$sizeNext, na.rm = TRUE)+3,c("Habitat 1","Habitat 2"),
	col=c("blue","dark green"),lty = c(1,2), merge = TRUE)
	
	


	
	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
