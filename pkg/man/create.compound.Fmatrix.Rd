\name{create.compound.Fmatrix}
\alias{create.compound.Fmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build compound F matrix
}
\description{
  Uses fecundity object, and environmental
  transition objects to construct a matrix
  defining probabilities for transitions between sizes due fertility
  given both size and environmental state and discrete stages. 
}
\usage{
create.compound.Fmatrix(n.env.class = 2, n.big.matrix = 50,
minsize = -1, maxsize = 50, envMatrix, fecObj, integrate.type, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.env.class}{
    numeric, number of environmental classes, defaults to 2
  }
  \item{n.big.matrix}{
    numeric, number of size bins in the T matrix, defaults to 50
  }
  \item{minsize}{
    numeric, minimum size of the T matrix, defaults to -1
  }
  \item{maxsize}{
    numeric, maximum size of the T matrix, defaults to 50
  }
  \item{envMatrix}{
    env.matrix object defining transitions between environmental states
    for each size
  }
  \item{fecObj}{
    fertility object
  }
  \item{integrate.type}{
    integration type, defaults to "midpoint" (which uses prob. density
    function); other option is "cumul" (which uses the cumul. density function)
  }
  \item{correction}{
    correction type, defaults to "none"; option is "constant" which will
    multiply every column of the IPM by a constant sufficient to adjust
    values to those predicted for survival at that size. 
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  an object of class IPM.matrix with dimensions
  n.big.matrix*n.env.class, or if discrete transitions exist
  (n.big.matrix+n.disc)*n.env.class

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean McMahon, Roberto Salguero-Gómez & Eelke Jongejans
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{create.compound.Fmatrix}},\code{\link{create.IPM.Tmatrix}}
}
\examples{
#without discrete stages
dff<-generateData()
Fmatrix <-create.compound.Fmatrix(n.big.matrix=20,minsize=min(dff$size, na.rm=TRUE),
maxsize=max(dff$size, na.rm=TRUE),envMatrix=makeEnvObj(dff),
fecObj=makeFecObjGeneral(dff, explanatoryVariables="size+size2+covariate"))

image(1:nrow(Fmatrix),1:ncol(Fmatrix),log(Fmatrix), xlab="Size at t+1",
ylab="Size at t", axes=FALSE)
axis(1,at=1:nrow(Fmatrix),lab=round(rep(Fmatrix@meshpoints,Fmatrix@n.env.class),2))
axis(2,at=1:nrow(Fmatrix),lab=round(rep(Fmatrix@meshpoints,Fmatrix@n.env.class),2))
abline(h=length(Fmatrix@meshpoints)*(1:Fmatrix@n.env.class))
abline(v=length(Fmatrix@meshpoints)*(1:Fmatrix@n.env.class))



#with discrete stages
dff<-generateDataDiscrete()
dff$covariate <- sample(1:3,size=nrow(dff),replace=TRUE)
dff$covariatenext <- sample(1:3,size=nrow(dff),replace=TRUE)
fv1 <- makeFecObjGeneral(dff,
explanatoryVariables="size",Transform="none",
offspring.splitter=data.frame(continuous=0.9,dormant=0.1))
Fmatrix <-create.compound.Fmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff),
fecObj=fv1)




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
