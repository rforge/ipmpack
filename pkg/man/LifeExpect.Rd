\name{LifeExpect}
\alias{LifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate Life Expectancy in a discretely stochastic environment
}
\description{
  Provided a compound Tmatrix that describes survival transitions across size and
  across discrete environmental states, this function provides a matrix defining life-expectancy in
  units of the time-step used, for each of the size bins in each of the
  environments (rows). 
}
\usage{
LifeExpect(IPM.matrix, env.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object describing growth and survival transitions
    across size and environment
}
  \item{env.matrix}{
    an env.matrix object describing transition between environments
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  matrix, each row corresponds to life expectancy in a different environments, columns
  correspond to sizes provided in Tmatrix@meshpoints
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p110-132

Cochran & Ellner. 1995. Simple methods for calculating age-based life history parameters for stage-structured populations. Ecological Monographs 62 p345-364.

Tuljapurkar & Horvitz, 2006. From stage to age in variable environments. Life expectancy and survivorship. Ecology 87 p1497-1509
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MeanLifeExpect}} \code{\link{create.compound.Tmatrix}}
}
\examples{
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
PassageTime(5,Tmatrix)

Tmatrix <-create.compound.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),
growObj=makeGrowthObj(dff,explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObj(dff), envMatrix=makeEnvObj(dff))
le <- LifeExpect(Tmatrix, makeEnvObj(dff))
plot(Tmatrix@meshpoints,le[1,],type="l", ylim=range(le))
points(Tmatrix@meshpoints,le[2,],type="l",col=2)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
