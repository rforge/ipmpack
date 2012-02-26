\name{predictFutureDistribution}
\alias{predictFutureDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict size distribution in the future giving current population's size distribution
}
\description{
  Function to project a population forwards using an IPM and a starting
  environment. The IPM may be structured by size alone, or by size and
  the environment. 
}
\usage{
predictFutureDistribution(startingSizes, IPM, n.time.steps, startingEnv = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{startingSizes}{
    vector containing the sizes of the desired starting population
}
  \item{IPM}{
    an IPM.matrix obejct (Tmatrix if only interested in survival
    projections, Tmatrix+Fmatrix otherwise)
}
  \item{n.time.steps}{
    a numeric defining the number of time-steps for which projection is required
}
  \item{startingEnv}{
    vector defining the desired starting environment, of length one or
    length startingSizes; ignored if no
    environmental states are provided; otherwise if the length is less
    than startingSizes assumes all individuals start in the same
    environment, given by startingEnv[1]. 
}
}
\details{
Currently this does not accept IPMs with discrete stages (e.g. seedbank)
}
\value{
 \item{n.new.dist0}{Starting frequency distribution along meshpoints in IPM.matrix}
 \item{n.new.dist}{Final frequency distribution}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gomez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

#define starting pop of interest
startPop <- rnorm(1000,4,1)

#build T and F matrix
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObj(dff),
survObj=makeSurvObj(dff))
Fmatrix <-create.IPM.Fmatrix(minSize=min(dff$size, na.rm = TRUE),
maxSize=max(dff$size, na.rm = TRUE), fecObj=makeFecObj(dff))

#make an IPM.matrix object containing Tmatrix+Fmatrix by replacing the Tmatrix
IPM <- Tmatrix
IPM@.Data <- Tmatrix+Fmatrix

#project five steps
a5 <- predictFutureDistribution(startingSizes=startPop,IPM=IPM,
n.time.steps=5,startingEnv=1)

#plot
 plot(IPM@meshpoints,a1$n.new.dist,type="l",
 ylim=range(c(a1$n.new.dist,a1$n.new.dist0)), xlab="Size", ylab="Number of individuals")
 points(IPM@meshpoints,a5$n.new.dist0,type="l",lty=2,col="red")
legend("topright", legend=c("Start distribution", "Distribution in 5 steps"), col=c(1,2),lty=c(1,2))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
