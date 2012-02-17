\name{create.compound.Tmatrix}
\alias{create.compound.Tmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to build compound T matrix
}
\description{
  Uses growth, survival, discreteTrans, and environmental
  transition objects to construct a matrix
  defining probabilities for transitions between sizes due to growth and survival
  given both size and environmental state and discrete stages. 
}
\usage{
create.compound.Tmatrix(n.env.class = 2, n.big.matrix = 50,
minsize = -1, maxsize = 50, envMatrix, growObj, survObj, integrate.type, correction)
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
  \item{growObj}{
    growth object
  }
  \item{survObj}{
    survival object
  }
  \item{DiscreteTrans}{
    object of class DiscreteTrans, or numeric
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
%%  ~~who you are~~
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
Tmatrix <-create.compound.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff),
growObj=makeGrowthObjGeneral(dff, explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObjGeneral(dff,explanatoryVariables="size+size2+covariate"))


#with discrete stages
dff<-generateDataDiscrete()
dff$covariate <- sample(1:3,size=nrow(dff),replace=TRUE)
dff$covariatenext <- sample(1:3,size=nrow(dff),replace=TRUE)
discM <- makeDiscreteTrans(dff)
Tmatrix <-create.compound.Tmatrix(minsize=min(dff$size,na.rm=TRUE),
maxsize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff),
growObj=makeGrowthObjGeneral(dff, explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObjGeneral(dff,explanatoryVariables="size+size2+covariate"),DiscreteTrans=discM)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
