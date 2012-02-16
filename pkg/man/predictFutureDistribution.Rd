\name{predictFutureDistribution}
\alias{predictFutureDistribution}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Predict size distibution in the future giving current population
}
\description{
  Function to project a population forwards using an IPM and a starting
  environment. The IPM may be structured by size alone, or by size and
  the environment. 
}
\usage{
predictFutureDistribution(starting.sizes, IPM, n.time.steps, starting.env = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{starting.sizes}{
    vector containing the sizes of the desired starting population
}
  \item{IPM}{
    an IPM.matrix obejct (Tmatrix if only interested in survival
    projections, Tmatrix+Fmatrix otherwise)
}
  \item{n.time.steps}{
    a numeric defining the number of time-steps for which projection is required
}
  \item{starting.env}{
    vector defining the desired starting environment, of length one or
    length starting.sizes; ignored if no
    environmental states are provided; otherwise if the length is less
    than starting.sizes assumes all individuals start in the same
    environment, given by starting.env[1]. 
}
}
\details{
}
\value{
 \item{n.new.dist0}{Starting frequency distribution along meshpoints in IPM.matrix}
 \item{n.new.dist}{Final frequency distribution}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
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
startpop <- rnorm(1000,4,1)

#build T and F matrix
dff <- generateData()
Tmatrix <-create.IPM.Tmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE),growObj=makeGrowthObjGeneral(dff),
survObj=makeSurvObjGeneral(dff))
Fmatrix <-create.IPM.Fmatrix(minsize=min(dff$size, na.rm = TRUE),
maxsize=max(dff$size, na.rm = TRUE), fecObj=makeFecObjGeneral(dff))

#make an IPM.matrix object containing Tmatrix+Fmatrix by replacing the Tmatrix
IPM <- Tmatrix
IPM@.Data <- Tmatrix+Fmatrix

#project
a1 <- predictFutureDistribution(starting.sizes=startpop, IPM=IPM,
n.time.steps=4, starting.env = 1)

#plot
 plot(IPM@meshpoints,a1$n.new.dist,type="l",
 ylim=range(c(a1$n.new.dist,a1$n.new.dist0)), xlab="Size", ylab="")
 points(IPM@meshpoints,a1$n.new.dist0,type="l",col=2)
legend("topright", legend=c("start dist", "end dist"), col=c(2,1),lty=1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
