\name{sensParams}
\alias{sensParams}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates sensitivity and elasticity of lambda (or R0, or Life
  expectancy of a chosen bin) to parameters underlying an IPM.
}
\description{
  Uses perturbation to estimate the sensitivity and elasticity of
  all the parameters underlying an IPM. 
}
\usage{
sensParams(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize,  
	chosenCov = data.frame(covariate = 1),discreteTrans = 1, 
	integrateType = "midpoint", correction = "none", preCensus = TRUE, 
	delta=1e-4, response="lambda", chosenBin=1)

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
  \item{chosenCov}{
level or value of the covariate(s) at which sensitivity estimation is desired
  }
  \item{response}{
	 whether lambda, R0 or life expectancy of a desired bin (lifeExpect with chosenBin) is required 
  }
  \item{chosenBin}{
which is required for lifeExpect
  }
  
  }
\value{
  \item{sens }{a vector of sensitivities of lambda or other variable with names
    corresponding to parameters.}
  \item{elas }{a vector of elasticities to lambda or other variable with names
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
  The values returned by sensParam are calculated by first calculating
  lambda for the chosen IPM; then modifying the focal parameter c by a
  very small amount, c.new=c*(1+delta) (the default for delta =1e-4, but users may specify the value that they want). The function then rebuilds the T and F matrices, and re-calculates lambda. Sensitivity is calculated as: 
  
  sens = df(x)/dx = (lam.new-lam)/(c*delta)
  
  i.e., the function estimates the degree to which a small change in the parameter results in a small change in lambda; and elasticity is calculated as: 
  
  elas = sens*c/lam = (lam.new-lam)/(lam*delta)
  
  which corresponds to the proportional change in lambda as an outcome
  of the proportional change in the parameter; analagous calculations
  are used for R0 and life expectancy. 

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}, \code{\link{elas}}
}
\examples{
dff <- generateData()

#lambda
res <- sensParams(growObj = makeGrowthObj(dff), 
survObj = makeSurvObj(dff), fecObj = makeFecObj(dff, Transform="log"), nBigMatrix = 50, 
minSize = min(dff$size, na.rm=TRUE), maxSize = max(dff$size, na.rm = TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$sens, 
main = expression("Parameter sensitivity of population growth rate "* lambda), 
las = 2, cex.names = 0.5) 
barplot(res$elas, 
main = expression("Parameter elasticity of population growth rate "* lambda), 
las = 2, cex.names = 0.5) 

#R0
resR0 <- sensParams(growObj = makeGrowthObj(dff), 
survObj = makeSurvObj(dff), fecObj = makeFecObj(dff, Transform="log"), nBigMatrix = 50, 
minSize = min(dff$size, na.rm=TRUE), maxSize = max(dff$size, na.rm = TRUE), response="R0")

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(resR0$sens, 
main = expression("Parameter sensitivity of net reproductive rate R"[0]), 
las = 2, cex.names = 0.5) 
barplot(resR0$elas, 
main = expression("Parameter elasticity of net reproductive rate R"[0]), 
las = 2, cex.names = 0.5) 

#life expectancy
resLE <- sensParams(growObj = makeGrowthObj(dff), 
survObj = makeSurvObj(dff),  nBigMatrix = 50, 
minSize = min(dff$size, na.rm=TRUE), maxSize = max(dff$size, na.rm =
TRUE), chosenBin=1, response="lifeExpect")

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(resLE$sens, 
main = expression("Parameter sensitivity of Life Expectancy"*eta[0]), 
las = 2, cex.names = 0.5) 
barplot(resLE$elas, 
main = expression("Parameter elasticity of Life expectancy"*eta[0]), 
las = 2, cex.names = 0.5) 

# Same as lambda above, but with two fecundity functions
dff$fec2 <- dff$fec>0 #create binomial describing e.g., prob of flowering
dff$fec[dff$fec==0] <- NA #take out zeros to avoid infinity when fit on a log scale
fv1 <- makeFecObj(dff, Formula = c(fec~size+size2,fec2~size), Transform=c("log","none"),
Family = c("gaussian","binomial"))

res <- sensParams(growObj=makeGrowthObj(dff), survObj = makeSurvObj(dff), 
fecObj = fv1, nBigMatrix = 50, minSize = min(dff$size, na.rm = TRUE), 
maxSize = max(dff$size, na.rm = TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$sens, 
main = expression("Parameter sensitivity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 
barplot(res$elas, 
main = expression("Parameter elasticity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 

# Same but with two fecundity functions and a constant
fv1@fecConstants[1] <-0.5
res <- sensParams(growObj = makeGrowthObj(dff), survObj = makeSurvObj(dff), 
fecObj = fv1, nBigMatrix = 50, minSize = min(dff$size, na.rm = TRUE), 
maxSize = max(dff$size, na.rm = TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$sens, 
main = expression("Parameter sensitivity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 
barplot(res$elas, 
main = expression("Parameter elasticity of population growth rate " *lambda), 
las = 2, cex.names = 0.5)
}
