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
minSize = -1, maxSize = 50, envMatrix, growObj, survObj, integrate.type, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n.env.class}{
    numeric, number of environmental classes, defaults to 2
  }
  \item{n.big.matrix}{
    numeric, number of size bins in the T matrix, defaults to 50
  }
  \item{minSize}{
    numeric, minimum size of the T matrix, defaults to -1
  }
  \item{maxSize}{
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
Note that this structure can also be used to define size x age IPMs, where the transition 
between ages is reflected by a similar matrix
}
\value{
  an object of class IPM.matrix with dimensions
  n.big.matrix*n.env.class, or if discrete transitions exist
  (n.big.matrix+n.disc)*n.env.class
}
\references{

For information on T matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112

For habitat x stage modeling: Horvitz & Schemske. 1995. Spatiotemporal variation in demographic transitions of a tropical understory herb. Projection matrix analysis. Ecological Monographs 65 p155-192

For age x size modeling: Garcia, Dahlgren & Ehrlen. 2011. No evidence of senescence in a 300-year-old mountain herb. Journal of Ecology 99 p1424Ð1430

For general information:

Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81 p694-708

Ellner & Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167 p410-428

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
\code{\link{create.compound.Fmatrix}},\code{\link{create.IPM.Tmatrix}}
}
\examples{
#without discrete stages
dff<-generateData()
Tmatrix <-create.compound.Tmatrix(minSize=min(dff$size,na.rm=TRUE), maxSize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff), growObj=makeGrowthObj(dff, explanatoryVariables="size+size2+covariate"), survObj=makeSurvObj(dff,explanatoryVariables="size+size2+covariate"))

image(1:nrow(Tmatrix),1:ncol(Tmatrix),t(log(Tmatrix)), xlab="Size at t",ylab="Size at t+1",axes=FALSE)
axis(1,at=1:nrow(Tmatrix),lab=round(rep(Tmatrix@meshpoints,Tmatrix@n.env.class),2))
axis(2,at=1:nrow(Tmatrix),lab=round(rep(Tmatrix@meshpoints,Tmatrix@n.env.class),2))
abline(h=length(Tmatrix@meshpoints)*(1:Tmatrix@n.env.class))
abline(v=length(Tmatrix@meshpoints)*(1:Tmatrix@n.env.class))

#with discrete stages
dff<-generateDataDiscrete()
dff$covariate <- sample(1:3,size=nrow(dff),replace=TRUE)
dff$covariateNext <- sample(1:3,size=nrow(dff),replace=TRUE)
discM <- makeDiscreteTrans(dff)
Tmatrix <-create.compound.Tmatrix(minSize=min(dff$size,na.rm=TRUE),
maxSize=max(dff$size,na.rm=TRUE),envMatrix=makeEnvObj(dff),
growObj=makeGrowthObj(dff, explanatoryVariables="size+size2+covariate"),
survObj=makeSurvObj(dff,explanatoryVariables="size+size2+covariate"),DiscreteTrans=discM)

image(1:nrow(Tmatrix),1:ncol(Tmatrix),t(log(Tmatrix)), xlab="Size at t",ylab="Size at t+1",axes=FALSE)
axis(1,at=1:nrow(Tmatrix),lab=round(rep(Tmatrix@meshpoints,Tmatrix@n.env.class),2))
axis(2,at=1:nrow(Tmatrix),lab=round(rep(Tmatrix@meshpoints,Tmatrix@n.env.class),2))
abline(h=length(Tmatrix@meshpoints)*(1:Tmatrix@n.env.class))
abline(v=length(Tmatrix@meshpoints)*(1:Tmatrix@n.env.class))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
