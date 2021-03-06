\name{coerceMatrixIPM}
\alias{coerceMatrixIPM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Coerces a matrix population model to a IPMmatrix object.
}
\description{
Forces a MPM (a matrix population model, or survival-growth transition matrix, or similar) into an IPMmatrix object for use with passage time,
life expectancy methods, etc. 
}
\usage{
coerceMatrixIPM(amat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{amat}{
    A matrix population model; or subset (P matrix, F matrix) with any number of bins.
}
}
\value{
  an object of class IPMmatrix.
}
\author{
C. Jessica E. Metcalf, Sean M. McMahon, Roberto Salguero-Gomez & Eelke Jongejans.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

# Build a simple matrix model with three classes, 
# no retrogression and stasis only in the third stage, 
# which is the only reproductive one in the life cycle.
aMat <- matrix(0, 3, 3)
aMat[cbind(2:3, 1:2)] <- 0.5
aMat[cbind(3, 3)] <- 0.8
aMat[1, 3] <-5

IPMamat <- coerceMatrixIPM(aMat)
}
