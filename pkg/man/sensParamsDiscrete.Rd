\name{sensParamsDiscrete}
\alias{sensParamsDiscrete}
\alias{sensParamsDiscreteR0}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates sensitivity and elasticity of lambda to parameters
  underlying the discrete survival and growth transition section of an IPM.
}
\description{
  Uses perturbation to estimate the sensitivity and elasticity of
  all the discrete transition parameters linked to survival and growth
  underlying an IPM. 
}
\usage{
sensParamsDiscrete(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, delta=1e-4)

sensParamsDiscreteR0(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
    discreteTrans = 1, integrateType = "midpoint", correction = "none", preCensus = TRUE, delta=1e-4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObj}{
    a growth object.
  }
  \item{survObj}{
    a survival object.
  }
  \item{fecObj}{
    a fecundity object (not necessary for life
    expectancy analysis).
  }
  \item{nBigMatrix}{
    numeric, number of bins of size used in the IPM matrix.
  }
  \item{minSize}{
    numeric, minimum size used for meshpoints of the IPM matrix.
  }
  \item{maxSize}{
    numeric, maximum size used for meshpoints of the IPM matrix.
  }
  \item{discreteTrans}{
	matrix of discrete transitions; or 1 if there is none
  }
  \item{integrateType}{
    integration type.
  }
  \item{correction}{
    correction.
  }
  \item{preCensus}{
    boolean saying if preCensus fertility or not (not necessary for life
    expectancy analysis)
  }
  \item{delta}{
    size of the perturbation desired
  }
  \item{chosenBin}{
    for analysis of life expectancy, which bin in the IPM Life
    expectancy should be compared for
  }
  }
\value{
  \item{slam }{a vector of sensitivities of lambda with names
    corresponding to parameters.}
  \item{elam }{a vector of elasticities to lambda with names
    corresponding to parameters.}
}
\references{
Rees & Rose. 2002. Evolution of flowering strategies in Oenothera
glazioviana: an integral projection model approach. Proceedings of the Royal Society London Seres B 269, p1509-1515.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}
\note{
Modified following code developed by Rees & Rose 2002 (above).  
}
\details{
  The values returned by sensParamDiscrete are calculated by first calculating
    lambda for the chosen IPM; then modifying the focal parameter c by a
    very small amount, c.new=c*(1+delta) (the default for delta =1e-4, but users may specify the value that they want). The function then rebuilds the T and F matrices, and re-calculates lambda. Sensitivity is calculated as: 
  
  sens = df(x)/dx = (lam.new-lam)/(c*delta)
  
  i.e., the function estimates the degree to which a small change in the parameter results in a small change in lambda; and elasticity is calculated as: 
  
  elas = sens*c/lam = (lam.new-lam)/(lam*delta)
  
  which corresponds to the proportional change in lambda as an outcome of the proportional change in the parameter.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}, \code{\link{elas}}
}
\examples{


dff <- generateData()
rc<-sensParamsDiscrete(growObj = makeGrowthObj(dff),
survObj = makeSurvObj(dff),
fecObj = makeFecObj(dff,Transform="log"),
nBigMatrix=50,minSize=-2,maxSize=10,
discreteTrans=makeDiscreteTrans(dff))

barplot(rc[[1]],main = expression("Parameter sensitivity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 

barplot(rc[[2]],main = expression("Parameter elasticity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 

}
