\name{MeanLifeExpect}
\alias{MeanLifeExpect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to calculate mean life expectancy
}
\description{
  Provided a T matrix, which defines survival transitions across stages, this function outputs a vector defining life expectancy in units of the time-step used (see convertIncrement()), for each of the size bins
}
\usage{
MeanLifeExpect(IPM.matrix)
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
  A vector of life expectancies each corresponding to a value of the size bins defined by Tmatrix@meshpoints
}
\references{
Caswell, 2001. Matrix population models: analysis, construction and
interpretation. 2nd ed. Sinauer. p118-120 

Cochran & Ellner. 1992. Simple methods for calculating age-based life history parameters for stage-structured populations. Ecological Monographs 62, p345-364

}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{LifeExpect}}, \code{\link{create.IPM.Tmatrix}} 
}
\examples{
##with a single continuous state variable (e.g. size)
dff<-generateData()
Tmatrix<-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
MeanLifeExpect(Tmatrix)

Tmatrix<-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))

plot(MeanLifeExpect(Tmatrix),ylab="Mean life expectancy",xlab="Continuous stage (e.g. size)",type="l",col="dark gray",ylim=c(0,max(MeanLifeExpect(Tmatrix))))


##with discrete stages (e.g. seed bank)
dff<-generateDataDiscrete()
dff$covariate <- sample(1:3,size=nrow(dff),replace=TRUE)
dff$covariateNext <- sample(1:3,size=nrow(dff),replace=TRUE)
discM<-makeDiscreteTrans(dff)
Tmatrix<-create.compound.Tmatrix(minSize=min(dff$size,na.rm=TRUE),
maxSize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff),
growObj=makeGrowthObj(dff, explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObj(dff,explanatoryVariables="size+size2+covariate"),DiscreteTrans=discM)
MeanLifeExpect(Tmatrix)

plot(MeanLifeExpect(Tmatrix),ylab="Mean life expectancy",xlab="Continuous (e.g. size) and discrete (e.g. seedbank) stages",type="l",col="dark gray",ylim=c(0,max(MeanLifeExpect(Tmatrix))))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
