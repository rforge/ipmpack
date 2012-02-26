\name{LifeExpect}
\alias{LifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate life expectancy in a discretely stochastic environment
}
\description{
  Provided a compound T matrix that describes survival transitions across size and
  across discrete environmental states, this function provides a matrix defining life expectancy in
  units of the time-step used in the census, for each of the size bins in each of the
  environments (rows). 
}
\usage{
LifeExpect(IPM.matrix, envMatrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object describing growth and survival transitions
    across size and environment
}
  \item{envMatrix}{
    an envMatrix object describing transition between environments
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
##with a single continuous state variable (e.g. size)
dff<-generateData()
Tmatrix<-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
plot(Tmatrix@meshpoints,le[1,],ylab="Life expectancy",xlab="Size",type="l", ylim=range(le))
points(Tmatrix@meshpoints,le[2,],type="l",col=2)


##with a compound matrix (discrete stages not yet defined)
Tmatrix<-create.compound.Tmatrix(minSize=min(dff$size,na.rm=TRUE),maxSize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff,explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObj(dff),envMatrix=makeEnvObj(dff))
le<-LifeExpect(Tmatrix,makeEnvObj(dff))
plot(Tmatrix@meshpoints,le[1,],ylab="Life expectancy",xlab="Size",type="l", ylim=range(le))
points(Tmatrix@meshpoints,le[2,],type="l",col=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
