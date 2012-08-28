\name{createIntegerPmatrix}
\alias{createIntegerPmatrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Builds P matrices where the continuous measurement is discrete.
}
\description{
  Uses growth and survival objects to construct a matrix defining probabilities for transitions between sizes due to growth and survival. Extensions for transition to discrete classes are possible.
}
\usage{
createIntegerPmatrix(nEnvClass = 1, 
		meshpoints=1:20,
		chosenCov = data.frame(covariate = 1), 
		growObj, survObj, 
		discreteTrans = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nEnvClass}{
    numeric, number of environmental classes, always = 1 for
    non-compound matrices.
  }
  \item{meshpoints}{
    numeric, identity of the discrete integers
  }
  \item{chosenCov}{
    data-frame indicating level of the discrete covariate, or 
    range of values where multiple covariates are modeled.
  }
  \item{growObj}{
    growth object.
  }
  \item{survObj}{
    survival object.
  }
  \item{discreteTrans}{
    object of class discreteTrans, or numeric.
  }
}
\details{
This is the discrete analogue of createIPMPatrix; using probability mass functions instead of pdfs
}
\value{
  an object of class IPMmatrix with dimensions length(meshpoints)*length(meshpoints), or
  length(meshpoints)+nrow(discreteTrans).
}
\references{
Easterling, Ellner & Dixon. 2000. Size-specific sensitivity: a new structured population model. Ecology 81, p694-708.

Ellner & Rees. 2006. Integral projection models for species with complex demography. The American Naturalist 167, p410-428.

For information on P matrix: Caswell. 2001. Matrix population models: construction, analysis, and interpretation. 2nd ed. Sinauer. p110-112.

%% ~put references to the literature/web site here ~
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{createIntegerFmatrix}}, \code{\link{createIPMPmatrix}}  
}
\examples{

}
