\name{createIntegerFmatrix}
\alias{createIntegerFmatrix}
\alias{createIntegerPmatrix}
%- Also NEED an '\alias' for EACH other topic documented here. 
\title{
Builds P and F matrices built off regressions fitted to discrete variables (with probability mass functions rather than pdfs).
}
\description{
Uses growth, survival and fecundity objects to construct a matrix defining per-capita contribution to recruitment stages (e.g., propagules [seed, spore], seedlings, calves) by reproductive stages due to sexual reproduction (for the F matrix); and contributions via and growth and survival to continuous stages. 
}
\usage{
createIntegerFmatrix(fecObj, nEnvClass = 1, meshpoints=1:20, 
	chosenCov = data.frame(covariate=1), 
	preCensus = TRUE, survObj = NULL, growObj = NULL)

createIntegerPmatrix(nEnvClass = 1, 
		meshpoints=1:20,
		chosenCov = data.frame(covariate = 1), 
		growObj, survObj, 
		discreteTrans = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fecObj}{
    fecundity object.
  }
  \item{nEnvClass}{
    numeric, number of environmental classes, always = 1 for
    non-compound matrices.
  }
  \item{meshpoints}{
    numeric, identifying meshpoints
  }
  \item{chosenCov}{
    data-frame indicating level of the discrete covariate, or 
    range of values where multiple covariates are modeled.
  }
  \item{preCensus}{
    logical (TRUE or FALSE), indicating whether the fecundity object should
    represent an interval between pre-breeding or a post-breeding censusses.
    defaults to TRUE (pre-breeding census), meaning that all reproduction 
    and offspring rates required for the F matrix are embedded in fecObj. 
    Alternatively, an F matrix based on post-breeding census (preCensus=FALSE)
    requires additional survObj and growObj, to cover the survival and growth
    of the parents until the reproduction event. 
  }
  \item{survObj}{
    survival object, describing the survival of parents from a census until 
    the reproduction event starts (at some point during the inter-census time step. 
  }
  \item{growObj}{
    growth object, describing the growth of parents that survive until the 
    reproduction event starts. Warning: this growth object is still ignored in 
    createIPMFmatrix in the current version of IPMpack. It will become functional 
    in coming versions.  
  }
   \item{discreteTrans}{
    object of class discreteTrans, or numeric.
  }
 
}
\value{
  an object of class IPMmatrix with dimensions length(meshpoints)*length(meshpoints), or
  length(meshpoints)+nrow(discreteTrans).
}
\details{
 do check whether the Pmatrix adequately reflects survival. 
}
\references{
For information on F and P matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
With thanks to Alden Griffith.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{makeDiscreteTransInteger}}, \code{\link{makeFecObjInteger}}
}
\examples{

# Generate data with only continuous stage and discretize sizes
dff <- generateData()
dff$fec[dff$fec==0] <- NA
dff$size <- pmax(floor(dff$size+10),0)
dff$sizeNext <- pmax(floor(dff$sizeNext+10),0)

# make objects
gr1<-makeGrowthObj(dff,Family="poisson")
sv1<-makeSurvObj(dff)
fv1 <- makeFecObjInteger(dff, Formula = log(fec)~size, Transform = "log")

#plot the new growth object
par(mfrow=c(1,3),bty="l")
plot(dff$size,dff$sizeNext, xlab="stage at t", ylab="stage at t+1")
points(1:20,predict(gr1@fit,newdata=data.frame(size=1:20), type="response"),type="l",col=2)

#make the P and T matrix
Fmatrix <- createIntegerFmatrix(fecObj = fv1, meshpoints=1:20)
Pmatrix <- createIntegerPmatrix(meshpoints=1:20,
			growObj=gr1, survObj=sv1, 
			discreteTrans = 1)


#plot them
image(Fmatrix@meshpoints, Fmatrix@meshpoints, t(Fmatrix), 
	xlab = "Stage at t", 
		ylab = "Stage at t+1")

image(Pmatrix@meshpoints, Pmatrix@meshpoints, t(Pmatrix), 
	xlab = "Stage at t", 
		ylab = "Stage at t+1")


#same but with negative binomial instead of poisson for stage transitions
gr1<-makeGrowthObj(dff,Family="negbin")

#plot
par(mfrow=c(1,3),bty="l")
plot(dff$size,dff$sizeNext, xlab="stage at t", ylab="stage at t+1")
points(1:20,predict(gr1@fit[[1]],newdata=data.frame(size=1:20),type="response"),type="l",col=2)

#make the P and T matrix
Fmatrix <- createIntegerFmatrix(fecObj = fv1, meshpoints=1:20)
Pmatrix <- createIntegerPmatrix(meshpoints=1:20,
			growObj=gr1, survObj=sv1, 
			discreteTrans = 1)


#plot them
image(Fmatrix@meshpoints, Fmatrix@meshpoints, t(Fmatrix), 
	xlab = "Stage at t", 
		ylab = "Stage at t+1")

image(Pmatrix@meshpoints, Pmatrix@meshpoints, t(Pmatrix), 
	xlab = "Stage at t", 
		ylab = "Stage at t+1")
}
