\name{VarLifeExpect}
\alias{VarLifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to calculate variation in life expectancy in a discretely stochastic environment
}
\description{
  Provided a T matrix, defining survival transitions across size, this function provides a vector with variance in life expectancy in units of the time-step used, for each of the size bins
}
\usage{
VarLifeExpect(IPM.matrix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{IPM.matrix}{
    an IPM.matrix object defining survival transitions
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  a vector of variance in life-expectancies each corresponding to Tmatrix@meshpoints

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
\code{\link{MeanLifeExpect}}, \code{\link{create.IPM.Tmatrix}} 

}
\examples{
##with a single continuous state variable (e.g. size)
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
vle<-VarLifeExpect(Tmatrix)

plot(Tmatrix@meshpoints,vle,ylab="Variation life expectancy",xlab="Size",type="l",ylim=c(0,max(vle)))


##with a compound matrix (discrete stages not yet defined)
Tmatrix<-create.compound.Tmatrix(minSize=min(dff$size,na.rm=TRUE),maxSize=max(dff$size,na.rm=TRUE),growObj=makeGrowthObj(dff,explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObj(dff),envMatrix=makeEnvObj(dff))
le<-LifeExpect(Tmatrix,makeEnvObj(dff))

plot(Tmatrix@meshpoints,vle,ylab="Variation life expectancy",xlab="Size",type="l", ylim=range(vle))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
