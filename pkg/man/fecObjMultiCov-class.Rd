\name{fecObjMultiCov-class}
\Rdversion{1.1}
\docType{class}
\alias{fecObjMultiCov-class}

\title{Class \code{"fecObjMultiCov"}}
\description{
A class object description for fertility objects with multiple covariates
%%  ~~ A concise (1-5 lines) description of what the class is. ~~
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("fecObjMultiCov", ...)}.
%%  ~~ describe objects here ~~ 
}
\section{Slots}{
  \describe{
    \item{\code{fitFec}:}{Object of class \code{"list"} ~~ }
    \item{\code{fecConstants}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{offspringSplitter}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{meanOffspringSize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{sdOffspringSize}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{offspringTypeRates}:}{Object of class \code{"data.frame"} ~~ }   
    \item{\code{fecByDiscrete}:}{Object of class \code{"data.frame"} ~~ }
    \item{\code{Transform}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
No methods defined with class "fecObjMultiCov" in the signature.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
showClass("fecObjMultiCov")
}
\keyword{classes}
