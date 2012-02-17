\name{sensParams}
\alias{sensParams}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Estimating sensitivity and elasticity of parameters underlying an IPM
}
\description{
  Uses perturbation to estimate the sensitivity and elasticity of
  all the parameters underlying an IPM that is structured only by size
  (compound T and F matrixes not possible yet)
}
\usage{
sensParams(growObj, survObj, fecObj, n.big.matrix, minsize, maxsize, integrate.type, correction)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{growObj}{
    a growth object
  }
  \item{survObj}{
    a survival object
  }
  \item{fecObj}{
    a fecundity object
  }
  \item{n.big.matrix}{
    numeric, number of bins of size used in the T matrix
  }
  \item{minsize}{
    numeric, minimum size used for meshpoints of the T matrix
  }
  \item{maxsize}{
    numeric, maximum size used for meshpoints of the T matrix
  }
  \item{integrate.type}{
    integration type
  }
  \item{correction}{
    correction
  }

  
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{slam }{a vector of sensitivities of lambda with names
    corresponding to parameters}
  \item{elam }{a vector of elasticities to lambda with names
    corresponding to parameters}
}
\references{
Rees & Rose, 2002 Evolution of flowering strategies in Oenothera
glazioviana: an integral projection model
approach. Proc. Roy. Soc. Lond. Ser. B. 269 p1509-1515
}
\author{
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sens}}, \code{\link{elas}}
}
\examples{
dff <- generateData()
res <- sensParams(growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff),
fecObj=makeFecObjGeneral(dff), n.big.matrix=50, minsize=min(dff$size, na.rm=TRUE),
maxsize=max(dff$size, na.rm=TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$slam, main = expression("Parameter sensitivity of " * 
 lambda), las = 2, cex.names = 0.5) 
barplot(res$elam, main = expression("Parameter elasticity of " * 
 lambda), las = 2, cex.names = 0.5) 

#same but with two fertility functions WON'T WORK ALTHOUGH DOES DIRECTLY???
dff$fec2 <- dff$fec>0 #create binomial describing e.g., prob of flowering
fv1 <- makeFecObjGeneral(dff,
explanatoryVariables=c("size+size2","size"), Family=c("gaussian","binomial"))

res <- sensParams(growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff),
fecObj=fv1, n.big.matrix=50, minsize=min(dff$size, na.rm=TRUE),
maxsize=max(dff$size, na.rm=TRUE))

par(mfrow = c(2, 1), bty = "l", pty = "m") 
barplot(res$slam, main = expression("Parameter sensitivity of " * 
 lambda), las = 2, cex.names = 0.5) 
barplot(res$elam, main = expression("Parameter elasticity of " * 
 lambda), las = 2, cex.names = 0.5) 


#same but with two fertility functions and a constant
%fv1@fec.constants[1] <- 0.5
%res <- sensParams(growObj=makeGrowthObjGeneral(dff), survObj=makeSurvObjGeneral(dff),
%fecObj=fv1, n.big.matrix=50, minsize=min(dff$size, na.rm=TRUE),
%maxsize=max(dff$size, na.rm=TRUE))

%par(mfrow = c(2, 1), bty = "l", pty = "m") 
%barplot(res$slam, main = expression("Parameter sensitivity of " * 
% lambda), las = 2, cex.names = 0.5) 
%barplot(res$elam, main = expression("Parameter elasticity of " * 
% lambda), las = 2, cex.names = 0.5) 




}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
