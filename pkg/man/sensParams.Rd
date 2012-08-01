\name{sensParams}
\alias{sensParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimates sensitivity and elasticity of parameters underlying an IPM.
}
\description{
  Uses perturbation to estimate the sensitivity and elasticity of
  all the parameters underlying an IPM  (compound T and F matrices not possible yet). If discrete stages are present, the sensitivity analysis can be run, but parameters linked to movement into and out of discrete stages will not be addressed. 
}
\usage{
sensParams(growObj, survObj, fecObj, nBigMatrix, minSize, maxSize, 
discreteTrans,integrateType, correction)
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
    a fecundity object.
  }
  \item{nBigMatrix}{
    numeric, number of bins of size used in the T matrix.
  }
  \item{minSize}{
    numeric, minimum size used for meshpoints of the T matrix.
  }
  \item{maxSize}{
    numeric, maximum size used for meshpoints of the T matrix.
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}, \code{\link{elas}}
}
\examples{
dff <- generateData()
res <- sensParams(growObj = makeGrowthObj(dff), 
survObj = makeSurvObj(dff), fecObj = makeFecObj(dff, Transform="log"), nBigMatrix = 50, 
minSize = min(dff$size, na.rm=TRUE), maxSize = max(dff$size, na.rm = TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$slam, 
main = expression("Parameter sensitivity of population growth rate "* lambda), 
las = 2, cex.names = 0.5) 
barplot(res$elam, 
main = expression("Parameter elasticity of population growth rate "* lambda), 
las = 2, cex.names = 0.5) 

# Same but with two fecundity functions
dff$fec2 <- dff$fec>0 #create binomial describing e.g., prob of flowering
dff$fec[dff$fec==0] <- NA #take out zeros to avoid infinity when fit on a log scale
fv1 <- makeFecObj(dff, Formula = c(fec~size+size2,fec2~size), Transform=c("log","none"),
Family = c("gaussian","binomial"))

res <- sensParams(growObj=makeGrowthObj(dff), survObj = makeSurvObj(dff), 
fecObj = fv1, nBigMatrix = 50, minSize = min(dff$size, na.rm = TRUE), 
maxSize = max(dff$size, na.rm = TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$slam, 
main = expression("Parameter sensitivity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 
barplot(res$elam, 
main = expression("Parameter elasticity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 

# Same but with two fecundity functions and a constant
fv1@fecConstants[1] <-0.5
res <- sensParams(growObj = makeGrowthObj(dff), survObj = makeSurvObj(dff), 
fecObj = fv1, nBigMatrix = 50, minSize = min(dff$size, na.rm = TRUE), 
maxSize = max(dff$size, na.rm = TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$slam, 
main = expression("Parameter sensitivity of population growth rate " *lambda), 
las = 2, cex.names = 0.5) 
barplot(res$elam, 
main = expression("Parameter elasticity of population growth rate " *lambda), 
las = 2, cex.names = 0.5)
}
